# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
# This file is part of code_aster.
#
# code_aster is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# code_aster is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
# --------------------------------------------------------------------

from math import exp, log, sqrt

import numpy
from libaster import ConvergenceError, IntegrationError, SolverError

from ...Cata.Syntax import _F
from ...Messages import MessageLog
from ...Utilities import cmp, force_list, logger, no_new_attributes, MPI
from ..Basics import EventId, Observer


class TimeStepper(Observer):
    """This object deals with the time steps.

    It gives the list of the time steps to be calculated. The initial time is
    not in the list since it is already known. The final time is.
    The initial time may be *None*, undeterminated. In this case, the first
    increment is also *None*.

    Arguments:
        times (list[float]): List of time steps to be calculated.
        epsilon (float, optional): Value used to check equality between two times
            (default: 1.e-16).
        initial (float, optional): Initial time (default: 0.0).
        final (float, optional): Final time (default: the last given).
    """

    _times = _eps = _current = _initial = _final = _last = None
    _actions = _state = None
    _split = _maxLevel = _minStep = _maxStep = None
    __setattr__ = no_new_attributes(object.__setattr__)

    default_increment = 1.0e-12
    maxNbSteps = 1.0e6

    def __init__(self, times, epsilon=default_increment, initial=0.0, final=None):
        super().__init__()
        times = list(times)
        if sorted(times) != times:
            raise ValueError("the time steps must be ordered")
        self._times = times
        self._eps = epsilon
        self._initial = initial
        self._final = final
        self._actions = []
        self._split = []
        self._maxLevel = -1
        self._minStep = 1.0e-99
        self._maxStep = 1.0e99
        self._resetState()
        logger.debug("TimeStepper.init: %s, %s, %s", initial, final, times)
        self._check_bounds()
        self._check_epsilon()

    @property
    def null_increment(self):
        """float: delta between two steps to be considered as null."""
        return self._eps

    @null_increment.setter
    def null_increment(self, value):
        """Setter for null_increment"""
        self._eps = value
        self._check_epsilon()

    def copy(self):
        """Return a copy of the object.

        Returns:
            TimeStepper: copy of the object.
        """
        new = TimeStepper(self._times, initial=self._initial, final=self._final, epsilon=self._eps)
        new._maxLevel = self._maxLevel
        new._minStep = self._minStep
        new._maxStep = self._maxStep
        for act in self._actions:
            new.register_event(act.copy())
        new.register_default_error_event()
        return new

    def _check_bounds(self):
        """Remove out of bounds values."""
        times = self._times
        if self._initial is not None:
            while times and self.cmp(times[0], self._initial) <= 0:
                times.pop(0)
        if self._final is not None:
            while times and self.cmp(times[-1], self._final) > 0:
                times.pop()
            if not times or self.cmp(self._final, times[-1]) > 0:
                times.append(self._final)
        self._current = 0
        if not times:
            # empty list
            self._last = -1
            self._final = self._initial
            return
        self._last = len(times) - 1
        self._final = times[-1]

    def _check_epsilon(self):
        """Remove values consistency with epsilon."""
        times = numpy.array(self._times)
        incr = times[1:] - times[:-1]
        if (incr < self._eps).any():
            raise ValueError(
                f"inconsistent values: increment is lesser than epsilon: {self._eps:f}"
            )

    def setInitial(self, time):
        """Define the initial time. Lesser values are removed.

        Calling `setInitial` resets the next calculated time at the first
        position.

        Arguments:
            time (float): First time to be used.
        """
        self._initial = time
        self._check_bounds()

    def setFinal(self, time=None, current=None):
        """Limit the sequence to the times lower than `time`.

        If `time` is not already in the sequence, it is appended.
        If `time` is not provided, the final time is set to the last one of
        the sequence.
        Calling `setFinal` resets the next calculated time at the first
        position. If `current` is provided, moves the current time at this step.

        Arguments:
            time (float, optional): Last time to be used.
        """
        self._final = time
        self._check_bounds()
        if current is not None:
            for idx, value in enumerate(self._times):
                if self.cmp(value, current) == 0:
                    self._current = idx

    def size(self):
        """Return the total number of steps in the list.

        Returns:
            int: Number of steps.
        """
        return len(self._times)

    def remaining(self):
        """Return the number of steps not yet calculated.

        Returns:
            int: Number of steps.
        """
        return self._last - self._current + 1

    def isFinished(self):
        """Tell if there are steps to be computed.

        Returns:
            bool: *True* if there is no step to be computed, *False* otherwise.
        """
        return self._current > self._last

    def _insert(self, index, time):
        """Inserts the step to given index. The caller must check for already
        existing time.

        Arguments:
            index (int): index to insert new time.
            time (float): time value to insert.
        """
        # print("\ninsert at", index, time, self._current, end=" ")
        if index < self._current:
            raise KeyError("can not insert a step before the current time.")
        if self.size() >= TimeStepper.maxNbSteps:
            logger.error(MessageLog.GetText("F", "ADAPTATION_13"))
        self._times.insert(index, time)
        self._last = len(self._times) - 1
        # print("\n->", self._current, self._last, self._times, flush=True)

    def getInitial(self):
        """Returns the initial time (not calculated).

        Returns:
            float: Initial time value.
        """
        return self._initial

    def getFinal(self):
        """Returns the last time to be calculated.

        Returns:
            float: Final time value.
        """
        return self._final

    def getPrevious(self):
        """Returns the previous calculated step.

        Returns:
            float: Previous time value.
        """
        if self.isFinished():
            raise IndexError("no more timesteps")
        if self._current == 0:
            return self._initial

        return self._times[self._current - 1]

    def getCurrent(self):
        """Returns the current step, this to be calculated.

        Returns:
            float: Current time value, this to be calculated.
        """
        if self.isFinished():
            raise IndexError("no more timesteps")
        return self._times[self._current]

    def getIncrement(self):
        """Returns the increment from the previous calculated step.

        Returns:
            float: increment from the previous time value.
        """
        prev = self.getInitial() if self._current == 0 else self.getPrevious()
        if prev is None:
            return None
        return self.getCurrent() - prev

    def getNextIncrement(self):
        """Returns the next increment to next step to be calculated.

        Returns:
            float: increment to the next time value.
        """
        if self._current >= self._last:
            return None
        return self._times[self._current + 1] - self.getCurrent()

    def completed(self):
        """Register the current step as completed successfully."""
        self._resetState()
        if self.isFinished():
            raise IndexError("no more timesteps")
        self._current += 1
        if not self.isFinished():
            last = self.getCurrent()
            logger.debug("check splitting level: %s %s", self._split, last)
            while self._split and self.cmp(last, self._split[-1]) >= 0:
                self._split.pop()

    def _resetState(self):
        self._state = {}
        self._state.setdefault("history", {})
        self._state.setdefault("criteria", {})
        self._state.setdefault("converged", {})

    def notify(self, event):
        """Receive notification from event, store the data from a future use
        by actions.

        Arguments:
            event (EventSource): Object that sends the notification.
        """
        state = self._state
        eid, data = event.get_state()
        if eid & EventId.IterationSolver:
            logger.debug("+ received from an iteration solver: %s", data)
            if "PRED" not in data.get("matrix", ""):
                hist = state["history"]
                crit = state["criteria"]
                for para in ("ITER_GLOB_MAXI", "RESI_GLOB_RELA", "RESI_GLOB_MAXI"):
                    hist.setdefault(para, []).append(data[para].value)
                    crit[para] = data[para].reference
        else:
            raise TypeError(f"unsupported event: eid={eid}")

    @property
    def splitting_level(self):
        """int: The current splitting level."""
        return len(self._split)

    def _append_split(self, value):
        """Store the next target for unsplitting.

        Arguments:
            value (int): Level.
        """
        self._split.append(value)

    def cmp(self, time1, time2):
        """Compare two times using epsilon.

        Arguments:
            time1 (float): first argument.
            time2 (float): second argument.

        Returns:
            int: -1 if time1 < time2, 0 if time1 == time2, +1 if time1 > time2
            using epsilon.
        """
        return cmp(time1, time2, abs_tol=self._eps)

    def __repr__(self):
        return f"<TimeStepper(from {self._initial} to {self._final}, size {self.size()}: {self._times})>"

    @classmethod
    def from_keywords(cls, **args):
        """Initialize a TimeStepper from user keywords as provided to
        the INCREMENT common set of keywords.

        Arguments:
            args (dict): keywords as for INCREMENT.

        Returns:
            TimeStepper: a new TimeStepper or a copy of the object from INCREMENT.
        """
        assert "LIST_INST" in args, "THER_NON_LINE not yet supported!"
        try:
            stp = args["LIST_INST"].stepper.copy()
            times = stp._times
        except AttributeError:
            # ListOfFloats
            times = args["LIST_INST"].getValues()
            stp = None
        initial = times[0]
        eps = args.get("PRECISION", cls.default_increment)
        if "INST_INIT" in args:  # because None has a special meaning
            initial = args["INST_INIT"]
        if args.get("NUME_INST_INIT"):
            initial = times[args["NUME_INST_INIT"]]
        final = args.get("INST_FIN")
        if args.get("NUME_INST_FIN"):
            final = times[args["NUME_INST_FIN"]]
        if stp is None:
            stp = TimeStepper(times, initial=initial, final=final, epsilon=eps)
        else:
            stp.setInitial(initial)
            stp.setFinal(final)
            stp.null_increment = eps
        stp.register_default_error_event()
        return stp

    @staticmethod
    def command_factory(args):
        """Create a TimeStepper from DEFI_LIST_INST keywords.

        *Transitional function during migration from legacy operators that need
        a TimesList object and the ones are used a TimeStepper.*

        Argumentss:
            args (dict): User keywords

        Returns:
            TimeStepper: a new TimeStepper.
        """
        args = _F(args)
        definition = args["DEFI_LIST"][0]
        if "VALE" in definition:
            times = definition["VALE"]
        elif "LIST_INST" in definition:
            times = definition["LIST_INST"].getValues()
        else:
            # this option does not exist with AUTO
            result = definition["RESULTAT"]
            div = definition["SUBD_PAS"]
            orig = [result.getTime(idx) for idx in result.getIndexes()]
            times = [orig.pop(0)]
            for step in orig:
                times.extend(numpy.linspace(times[-1], step, div + 1)[1:])
        stp = TimeStepper(times, initial=None)

        for fail in args["ECHEC"]:
            if fail["EVENEMENT"] == "ERREUR":
                event = TimeStepper.Error()
            elif fail["EVENEMENT"] == "DELTA_GRANDEUR":
                event = TimeStepper.MaximumIncrement(
                    fieldName=fail["NOM_CHAM"],
                    component=fail["NOM_CMP"],
                    maxValue=fail["VALE_REF"],
                    group=fail["GROUP_MA"] or fail["GROUP_NO"],
                )
            else:  # not yet supported, ignored
                continue
            if fail["ACTION"] == "ARRET":
                act = TimeStepper.Interrupt(event)
            elif fail["ACTION"] == "DECOUPE":
                if fail["SUBD_METHODE"] == "MANUEL":
                    stp._maxLevel = max(stp._maxLevel, fail["SUBD_NIVEAU"])
                    act = TimeStepper.Split(
                        event, nbSubSteps=fail["SUBD_PAS"], minStep=fail["SUBD_PAS_MINI"]
                    )
                else:
                    assert fail["SUBD_METHODE"] == "AUTO"
                    # TODO not supported yet
                    act = TimeStepper.AutoSplit(event, minStep=fail["SUBD_PAS_MINI"])
            else:  # not supported yet, ignored
                continue
            stp.register_event(act)

        if args["METHODE"] == "AUTO":
            if "PAS_MINI" in definition:
                stp._minStep = definition["PAS_MINI"]
            if "PAS_MAXI" in definition:
                stp._maxStep = definition["PAS_MAXI"]
            maxNbSteps = definition["NB_PAS_MAXI"]
            stp.register_event(TimeStepper.Finalize(TimeStepper.MaximumNbOfSteps(maxNbSteps)))
            for adapt in args["ADAPTATION"]:
                if adapt["EVENEMENT"] == "AUCUN":
                    continue
                elif adapt["EVENEMENT"] == "SEUIL":
                    # not supported yet, ignored
                    continue
                else:
                    event = TimeStepper.Always()
                if adapt["MODE_CALCUL_TPLUS"] == "FIXE":
                    act = TimeStepper.AdaptConst(event, factor=(1.0 + adapt["PCENT_AUGM"] / 100.0))
                elif adapt["MODE_CALCUL_TPLUS"] == "ITER_NEWTON":
                    act = TimeStepper.AdaptFromNbIter(event, nbRef=adapt["NB_ITER_NEWTON_REF"])
                elif adapt["MODE_CALCUL_TPLUS"] == "DELTA_GRANDEUR":
                    act = TimeStepper.AdaptIncrement(
                        event,
                        fieldName=adapt["NOM_CHAM"],
                        component=adapt["NOM_CMP"],
                        maxValue=adapt["VALE_REF"],
                        group=adapt["GROUP_MA"] or adapt["GROUP_NO"],
                    )
                else:
                    # IMPLEX not supported yet, ignored
                    pass
                stp.register_event(act)

        stp.register_default_error_event()
        return stp

    # event management
    def register_event(self, action):
        """Register an action to react to an event.

        Args:
            action (~TimeStepper.Action): Type of action.
        """
        self._actions.append(action)

    def register_default_error_event(self):
        """Register a default action for Error event if there is no one."""
        if not [act for act in self._actions if isinstance(act.event, TimeStepper.Error)]:
            self._maxLevel = 3
            self.register_event(TimeStepper.Split(TimeStepper.Error(), nbSubSteps=4))
        elif self._maxLevel < 0:
            # AutoSplit and no other Split defined: only use minStep
            self._maxLevel = int(1e6)

    def failed(self, exc):
        """React to a raised exception.

        Args:
            exc (AsterError): The exception just raised.
        """
        for act in self._actions:
            if not isinstance(act.event, TimeStepper.Error):
                continue
            if act.event.is_raised(exception=exc):
                act.call(
                    timeStepper=self,
                    exception=exc,
                    residuals=self._state["history"].get("RESI_GLOB_RELA", []),
                    criteria=self._state["criteria"],
                )
                return
        raise TypeError("should not pass here!")  # because of default error

    def check_event(self, phys_state):
        """Check if an event should be raised.

        Arguments:
            phys_state (PhysicalState): Current physical state.

        Returns:
            bool: *False* if something went wrong, *True* if the step is ok.
        """
        # compute increment
        delta = phys_state.getCurrentDelta()
        # check events & actions
        if not self._check_error_posteriori(delta):
            return False
        self._check_adapt(delta)
        return True

    def _check_error_posteriori(self, delta):
        """Check for ErrorPosteriori events."""
        for act in self._actions:
            if not isinstance(act.event, TimeStepper.ErrorPosteriori):
                continue
            if act.event.is_raised(timeStepper=self, delta=delta):
                return act.call(timeStepper=self)
        return True

    def _check_adapt(self, delta):
        """Check for AdaptAction actions."""
        # last step?
        if self.remaining() == 1:
            return
        delta_t = 2.1e6
        currIncr = self.getIncrement()
        for act in self._actions:
            if not isinstance(act, TimeStepper.AdaptAction):
                continue
            if delta_t > 2.0e6:
                logger.info(MessageLog.GetText("I", "ADAPTATION_1"))
            delta_t = min(delta_t, 1.1e6)
            if act.event.is_raised(delta=delta):
                try:
                    dt_i = act.call(timeStepper=self, delta=delta)
                    delta_t = min(delta_t, dt_i * currIncr)
                    logger.info(
                        MessageLog.GetText("I", "ADAPTATION_2", valk=act.name, valr=delta_t)
                    )
                except ValueError:
                    logger.info(MessageLog.GetText("I", "ADAPTATION_3", valk=act.name))
                    raise
            else:
                logger.info(MessageLog.GetText("I", "ADAPTATION_3", valk=act.name))
        if delta_t < 1.0e6:
            logger.info(MessageLog.GetText("I", "ADAPTATION_5", valr=delta_t))
            nextIncr = self.getNextIncrement()
            if self.cmp(delta_t, nextIncr) > 0:
                delta_t = nextIncr
            logger.info(MessageLog.GetText("I", "ADAPTATION_6", valr=delta_t))
            if delta_t > self._maxStep:
                dict_args = dict(valr=(delta_t, self._maxStep, self._maxStep))
                logger.info(MessageLog.GetText("I", "ADAPTATION_12", **dict_args))
                delta_t = self._maxStep
            if delta_t <= self._minStep:
                logger.error(MessageLog.GetText("F", "ADAPTATION_11", valr=delta_t))
            new = self.getCurrent() + delta_t
            index = self._current + 1
            if self.cmp(new, self._times[index]) < 0:
                self._insert(index, new)
        elif delta_t < 2.0e6:
            logger.info(MessageLog.GetText("I", "ADAPTATION_4", valr=currIncr))
        return True

    def checkMaxLevel(self, min):
        """Check that 'maxLevel' is greater than 'min'
        (for compatibility with STAT_NON_LINE).

        Arguments:
            min (int): Minimal value.
        """
        if -1 < self._maxLevel < min:
            raise ValueError(
                f"DEFI_LIST_INST/SUBD_NIVEAU is {self._maxLevel}, "
                f"it must be greater than or equal to {min}."
            )

    class Event:
        """Event raised in case of error."""

        # TODO
        # ResidualDivergence: raised if resi_t+2 > min(resi_t, resi_t+1) (DIVE_RESI)
        # MaximumResidual: raised if RESI_GLOB_MAXI > value before the end of iterations (RESI_MAXI)
        # COLLISION, INTERPENETRATION, INSTABILITE

        __setattr__ = no_new_attributes(object.__setattr__)

        def is_raised(self, **context):
            """Tell if the event is raised in this context.

            Returns:
                bool: *True* if the event is raised, *False* otherwise.
            """
            raise NotImplementedError("must be subclassed!")

    class Always(Event):
        """Event raised at each step."""

        @staticmethod
        def is_raised(**context):
            """Tell if the event is raised in this context.

            Returns:
                bool: Always *True* .
            """
            return True

    class Error(Event):
        """Event raised in case of generic error."""

        @staticmethod
        def is_raised(**context):
            """Tell if the event is raised in this context.

            Returns:
                bool: *True* if the event is raised, *False* otherwise.
            """
            return isinstance(
                context.get("exception"), (ConvergenceError, IntegrationError, SolverError)
            )

    class ErrorPosteriori(Event):
        """Event that may raise an error after the convergency."""

    class MaximumIncrement(ErrorPosteriori):
        """Event raised when the increment of a component exceeds a value
        (DELTA_GRANDEUR keyword).

        Arguments:
            fieldName (str): Name of the field to be checked.
            component (str): Component name.
            maxValue (float): Maximum increment.
            group (*misc*): Restrict the checking to a part of the model.
        """

        _fieldName = _cmp = _value = _group = None
        __setattr__ = no_new_attributes(object.__setattr__)

        def __init__(self, fieldName, component, maxValue, group=None):
            super().__init__()
            self._fieldName = fieldName
            self._cmp = component
            self._value = maxValue
            self._group = force_list(group or [])

        def is_raised(self, **context):
            """Tell if the event is raised in this context.

            Returns:
                bool: *True* if the event is raised, *False* otherwise.
            """
            delta = context.get("delta")
            if not delta:
                return False
            field = delta.get(self._fieldName)
            if not field:
                return False
            array = numpy.array(field.getValuesWithDescription(self._cmp, self._group)[0])
            maxIncr = max(abs(array))
            maxIncr = MPI.ASTER_COMM_WORLD.allreduce(maxIncr, MPI.MAX)
            raised = maxIncr > self._value
            logger.debug("check delta of %s: %s >? %s: %s", self._cmp, maxIncr, self._value, raised)
            if raised:
                logger.info(MessageLog.GetText("I", "MECANONLINE10_24"))
            return raised

    class MaximumNbOfSteps(ErrorPosteriori):
        """Event raised after a number of steps have been completed.

        Arguments:
            maxNbSteps (int): Maximum number of time steps.
        """

        _value = None
        __setattr__ = no_new_attributes(object.__setattr__)

        def __init__(self, maxNbSteps):
            super().__init__()
            self._value = maxNbSteps

        def is_raised(self, **context):
            """Tell if the event is raised in this context.

            Returns:
                bool: *True* if the event is raised, *False* otherwise.
            """
            stp = context.get("timeStepper")
            # + 1 because 'completed()' not yet called for the current
            done = stp.size() - stp.remaining() + 1
            return done >= self._value

    class Action:
        """This object allow to execute an action on a TimeStepper when an event
        occurs."""

        _event = None
        __setattr__ = no_new_attributes(object.__setattr__)

        def __init__(self, event) -> None:
            self._event = event

        def copy(self):
            """Return a copy of the object."""
            return self.__class__(self._event)

        @property
        def event(self):
            """TimeStepper.Event: event to react."""
            return self._event

        def call(self, **context):
            """Execute the action.

            Arguments:
                context (dict): Context of the event.

            Returns:
                bool: *False* in case of "error", *True* otherwise.
            """
            raise NotImplementedError("must be subclassed!")

    class AdaptAction(Action):
        """Action that provides a multiplicative factor for the next timestep."""

        def call(self, **context):
            """Execute the action.

            Arguments:
                context (dict): Context of the event.

            Returns:
                float: multiplicative factor.
            """
            raise NotImplementedError("must be subclassed!")

    class Interrupt(Action):
        """This action stops the calculation (keyword value: ARRET)."""

        def call(self, **context):
            """Execute the action. Raises an exception.

            Arguments:
                context (dict): Context of the event.
            """
            logger.info(MessageLog.GetText("I", "MECANONLINE10_30"))
            raise context["exception"]

    class Finalize(Action):
        """This action finalizes the calculation without error."""

        def call(self, **context):
            """Execute the action.

            Arguments:
                context (dict): Context of the event.

            Returns:
                bool: always *True*.
            """
            stp = context.get("timeStepper")
            step = stp.getCurrent()
            stp.setFinal(step, current=step)
            logger.info(MessageLog.GetText("I", "ADAPTATION_13"))
            return True

    class Split(Action):
        """This action adds intermediate timesteps (keyword value: DECOUPE).

        Arguments:
            nbSubSteps (int): Number of sub-steps.
            minStep (float): Minimal value of a step.
        """

        _nbSubSteps = _minStep = None
        __setattr__ = no_new_attributes(object.__setattr__)

        def __init__(self, event, nbSubSteps, minStep=1.0e-12):
            super().__init__(event)
            assert nbSubSteps > 1, nbSubSteps
            self._nbSubSteps = int(nbSubSteps)
            self._minStep = minStep

        def copy(self):
            """Return a copy of the object."""
            return self.__class__(self._event, self._nbSubSteps, self._minStep)

        def call(self, **context):
            """Execute the action.

            Arguments:
                context (dict): Context of the event.
            """
            logger.info(MessageLog.GetText("I", "SUBDIVISE_1"))
            stp = context["timeStepper"]
            last = stp.getCurrent()
            stp._append_split(last)
            if stp.splitting_level > stp._maxLevel:
                stop = TimeStepper.Interrupt(self._event)
                stop.call(exception=SolverError("SUBDIVISE_17", (), (stp._maxLevel,), ()))

            time_step = stp.getIncrement() / self._nbSubSteps
            logger.info(
                MessageLog.GetText("I", "SUBDIVISE_10", (), (self._nbSubSteps,), (last, time_step))
            )
            if time_step < self._minStep:
                stop = TimeStepper.Interrupt(self._event)
                stop.call(exception=SolverError("SUBDIVISE_53", (), (), (self._minStep,)))
            if time_step < stp.null_increment:
                logger.warning(MessageLog.GetText("A", "SUBDIVISE_54"))
            new = stp.getCurrent()
            for _ in range(max(0, self._nbSubSteps - 1)):
                new -= time_step
                stp._insert(stp._current, new)

    class AutoSplit(Action):
        """This action adds intermediate timesteps (keyword value: DECOUPE)
        using less parameters than the *Split* action.

        Arguments:
            minStep (float): Minimal value of a step.
        """

        _minStep = None
        _manual = None
        __setattr__ = no_new_attributes(object.__setattr__)

        def __init__(self, event, minStep):
            super().__init__(event)
            self._minStep = minStep
            self._manual = TimeStepper.Split(event, nbSubSteps=4, minStep=minStep)

        def copy(self):
            """Return a copy of the object."""
            return self.__class__(self._event, self._minStep)

        def call(self, **context):
            """Execute the action.

            Arguments:
                context (dict): Context of the event.
            """
            logger.info(MessageLog.GetText("I", "SUBDIVISE_2"))
            logger.info(MessageLog.GetText("I", "EXTRAPOLATION_1"))
            stp = context["timeStepper"]
            last = stp.getCurrent()
            residuals = context["residuals"]
            if len(residuals) < 6:
                logger.info(MessageLog.GetText("I", "EXTRAPOLATION_3"))
                split = TimeStepper.Split(self._event, nbSubSteps=4, minStep=self._minStep)
                ctxt = context.copy()
                ctxt["exception"] = SolverError("EXTRAPOLATION_3")
                split.call(**ctxt)
                return

            stp._append_split(last)
            if stp.splitting_level > stp._maxLevel:
                stop = TimeStepper.Interrupt(self._event)
                stop.call(exception=SolverError("SUBDIVISE_17", (), (stp._maxLevel,), ()))
            # ignore the first three points
            residuals = residuals[3:]
            nbSteps, ratio = self._splittingRatio(residuals, context["criteria"])

            dt0 = ratio * stp.getIncrement()
            dtn = (1.0 - ratio) * stp.getIncrement() / (nbSteps - 1)
            logger.info(MessageLog.GetText("I", "SUBDIVISE_11", (), (nbSteps,), (last, dt0, dtn)))
            if dtn < self._minStep:
                stop = TimeStepper.Interrupt(self._event)
                stop.call(exception=SolverError("SUBDIVISE_53", (), (), (self._minStep,)))
            if dtn < stp.null_increment:
                logger.warning(MessageLog.GetText("A", "SUBDIVISE_54"))
            new = stp.getCurrent()
            for _ in range(nbSteps - 1):
                new -= dtn
                stp._insert(stp._current, new)

        @staticmethod
        def _splittingRatio(residuals, criteria):
            """Compute the ratio of the first step."""
            cresi = criteria["RESI_GLOB_RELA"]
            minIter = criteria["ITER_GLOB_MAXI"]
            maxIter = criteria["ITER_GLOB_MAXI"]
            # nmdcae
            xdet, xa0, xa1 = TimeStepper.AutoSplit._extrapol(residuals)
            nbSteps = 4
            ratio0 = 24.0 / ((3.0 * nbSteps + 1) ** 2 - 1.0)
            if xdet < 1.0e-16:
                ratio = ratio0
            else:
                ciblen = (xa0 + xa1 * log(cresi)) / xdet
                if 1.2 * ciblen < minIter:
                    ratio = ratio0
                else:
                    if xa1 < 1.0e-16:
                        ratio = ratio0
                    else:
                        if ciblen - maxIter <= -10.0 * xa1 / xdet:
                            ratio = exp((ciblen - maxIter) * xdet / xa1)
                        else:
                            ratio = exp(-10.0)
                        ratio = 0.48485 * ratio
                        xxbb = (-1.0 + sqrt(1.0 + 24.0 / ratio)) / 3.0
                        if xxbb < 2.0:
                            nbSteps = 2
                            ratio = 0.5
                        else:
                            nbSteps = round(xxbb)
            logger.debug("splittingRatio: %s, %s", nbSteps, ratio)
            return nbSteps, ratio

        @staticmethod
        def _extrapol(residuals):
            """Extrapolation: (xa0 + iter*xa1) / xdet"""
            # nmdcrg
            xn = 0.0
            sx = 0.0
            sy = 0.0
            sxx = 0.0
            syx = 0.0
            for i, resi in enumerate(residuals):
                xx = log(resi)
                if i >= len(residuals) - 3:
                    weight = 2.0
                else:
                    weight = 1.0
                xn += weight
                sx += weight * xx
                sy += weight * (i + 3)
                sxx += weight * xx**2
                syx += weight * xx * (i + 3)
            xdet = -(sx**2) + sxx * xn
            xa0 = sxx * sy - sx * syx
            xa1 = -sx * sy + syx * xn
            return xdet, xa0, xa1

    class AdaptConst(AdaptAction):
        """This action returns a constant multiplicative factor for the next timestep.

        Arguments:
            factor (float): Multiplicative factor.
        """

        _factor = None
        __setattr__ = no_new_attributes(object.__setattr__)

        def __init__(self, event, factor):
            super().__init__(event)
            self._factor = factor

        def copy(self):
            """Return a copy of the object."""
            return self.__class__(self._event, self._factor)

        def call(self, **context):
            """Execute the action.

            Arguments:
                context (dict): Context of the event.
            """
            return self._factor

    class AdaptFromNbIter(AdaptAction):
        """This action returns a multiplicative factor for the next timestep
        calculated from the number of iterations of the step.

        Arguments:
            nbRef (float): Reference number of iterations.
        """

        _nbRef = None
        name = "ITER_NEWTON"
        __setattr__ = no_new_attributes(object.__setattr__)

        def __init__(self, event, nbRef):
            super().__init__(event)
            self._nbRef = nbRef

        def copy(self):
            """Return a copy of the object."""
            return self.__class__(self._event, self._nbRef)

        def call(self, **context):
            """Execute the action.

            Arguments:
                context (dict): Context of the event.
            """
            stp = context["timeStepper"]
            nbIter = stp._state["converged"].get("ITER_GLOB_MAXI", self._nbRef - 1)
            return sqrt(self._nbRef / (nbIter + 1))

    class AdaptIncrement(AdaptAction):
        """This action returns a multiplicative factor for the next timestep
        calculated from the increment of the field component.

        Arguments:
            fieldName (str): Name of the field to be checked.
            component (str): Component name.
            maxValue (float): Maximum increment.
            group (*misc*): Restrict the checking to a part of the model.
        """

        _fieldName = _cmp = _value = _group = None
        name = "DELTA_GRANDEUR"

        def __init__(self, event, fieldName, component, maxValue, group=None):
            super().__init__(event)
            self._fieldName = fieldName
            self._cmp = component
            self._value = maxValue
            self._group = force_list(group or [])

        def copy(self):
            """Return a copy of the object."""
            return self.__class__(self._event, self._fieldName, self._cmp, self._value, self._group)

        def call(self, **context):
            """Execute the action.

            Arguments:
                context (dict): Context of the event.
            """
            delta = context.get("delta")
            if not delta:
                raise ValueError
            field = delta.get(self._fieldName)
            if not field:
                raise ValueError
            array = numpy.array(field.getValuesWithDescription(self._cmp, self._group)[0])
            nonzero = array[numpy.flatnonzero(array)]
            factor = numpy.min(self._value / numpy.abs(nonzero))
            factor = MPI.ASTER_COMM_WORLD.allreduce(factor, MPI.MIN)
            logger.debug("check delta of %s / %s: %s", self._cmp, self._value, factor)
            return float(factor)

    # ITER_SUPPL
    # AUTRE_PILOTAGE
    # ADAPT_COEF_PENA
    # CONTINUE
