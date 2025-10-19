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

from libaster import deleteTemporaryObjects

from ...Utilities import no_new_attributes, profile
from ..Basics import ProblemType
from .base_step_solver import BaseStepSolver


class ThermalStepSolver(BaseStepSolver):
    """Solves a step, loops on iterations."""

    problem_type = ProblemType.Thermal

    _theta = None
    __setattr__ = no_new_attributes(object.__setattr__)

    @classmethod
    def create(cls, param=None):
        """Setup a solver for the given problem.

        Arguments:
            param (dict) : user keywords.

        Returns:
            *StepSolver*: A relevant *StepSolver* object.
        """
        calc_type, theta = param["TYPE_CALCUL"], None
        if calc_type == "TRAN":
            assert param["SCHEMA_TEMPS"]["SCHEMA"] == "THETA"
            theta = param["SCHEMA_TEMPS"]["THETA"]
        return cls(theta=theta)

    def __init__(self, theta=0.57):
        super(ThermalStepSolver, self).__init__()
        self._theta = theta

    def initialize(self):
        """Initialization."""
        super().initialize()
        self.state.primal_step = self.state.createPrimal(self.problem, 0.0)

    @profile
    def solve(self):
        """Solve a step.

        Raises:
            *ConvergenceError* exception in case of error.
        """
        logManager = self.createLoggingManager()
        logManager.printConvTableEntries()

        # Solve current iteration
        self._iterations_solv.setLoggingManager(logManager)
        self._iterations_solv.initialize()
        self.current_matrix = self._iterations_solv.solve(self.current_matrix)

        deleteTemporaryObjects()
        logManager.printConvTableEnd()
