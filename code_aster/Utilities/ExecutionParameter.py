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

# person_in_charge: mathieu.courtois at edf.fr

"""
:py:mod:`ExecutionParameter` --- Management of the execution parameters
***********************************************************************

A singleton object :py:class:`.ExecutionParameter` is created during the
initialization. Its main feature is to parse and store the arguments read from
the command line or passed to the :py:func:`~code_aster.CodeCommands.debut.init`
function.
It also stores Python objects that have to be available from :py:mod:`libaster`.
They will be available though properties of the :py:class:`.ExecutionParameter`
object.

"""

import json
import os
import os.path as osp
import platform
import re
import shutil
import sys
import warnings
from argparse import SUPPRESS, Action, ArgumentParser
from contextlib import contextmanager

try:
    import yaml
except ImportError:
    yaml = None

import libaster

from run_aster.export import Export

from .as_timer import Timer
from .base_utils import Singleton, no_new_attributes
from .compatibility import deprecate
from .logger import DEBUG, INFO, logger
from .options import Options
from .strfunc import convert
from .version import get_version_desc

# aster_pkginfo will only be available after installation
try:
    from .aster_pkginfo import version_info
except ImportError:
    version_info = ()


DEFAULT_MEMORY_LIMIT = 2047 if "32" in platform.architecture()[0] else 4096
DEFAULT_TIME_LIMIT = 86400
RCDIR = osp.abspath(
    osp.join(osp.dirname(__file__), os.pardir, os.pardir, os.pardir, os.pardir, "share", "aster")
)


class ExecutionParameter(metaclass=Singleton):
    """This class stores and provides the execution parameters.

    The execution parameters are defined by reading the command line or using
    the method `set_option()`.

    Attributes:
        _args (dict): Command line arguments and execution parameters.
        _bool (int): Bitwise combination of boolean parameters.
        _catalc (*CataLoiComportement*): Object that gives access to the
            catalog of behaviors.
        _unit (*LogicalUnitFile*): Class that manages the logical units.
        _syntax (*CommandSyntax*): Class that passes user keywords up to
            Fortran operator.
    """

    _singleton_id = "Utilities.ExecutionParameter"
    _argv = _args = _bool = None
    _timer = _export = _command_counter = None
    _global = {}
    global_objects = (
        "catalc",
        "logical_unit",
        "syntax",
        "print_header",
        "checksd",
        "testresu_print",
        "iniran",
        "getran",
    )

    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self):
        """Initialization of attributes"""
        self._argv = None
        self._args = {}
        # options also used in F90
        self._args["dbgjeveux"] = 0
        self._args["jxveri"] = 0
        self._args["sdveri"] = 0

        self._args["memory"] = 0.0
        self._args["tpmax"] = 0.0
        self._args["maxbase"] = 0.0
        self._args["numthreads"] = 0

        self._args["stage_number"] = 0
        self._args["max_check"] = 0
        self._args["max_print"] = 0

        # boolean (on/off) options
        self._bool = Options.Null
        # Export object and input/output files
        self._export = None
        self._args["link"] = []

        self._args["rcdir"] = RCDIR
        # TODO probably to be removed?
        self._args["repmat"] = "."
        self._args["repdex"] = "."

        self._args["hook_post_exec"] = None

        self._computed()

    def _computed(self):
        """Fill some "computed" values"""
        # hostname
        self._args["hostname"] = platform.node()
        # ex. i686/x86_64
        self._args["processor"] = platform.machine()
        # ex. Linux
        self._args["system"] = platform.system()
        # ex. 32bit/64bit
        self._args["architecture"] = platform.architecture()[0]
        # ex. 2.6.32...
        self._args["osrelease"] = platform.release()
        self._args["osname"] = platform.platform()
        version = version_info.version
        keys = ("parentid", "branch", "date", "from_branch", "changes", "uncommitted")
        self._args.update(list(zip(keys, version_info[1:])))
        self._args["version"] = ".".join(str(i) for i in version)
        self._args["versMAJ"] = version[0]
        self._args["versMIN"] = version[1]
        self._args["versSUB"] = version[2]
        self._args["exploit"] = version_info.branch.startswith("v")
        self._args["versionD0"] = "%2d.%02d.%02d" % version
        self._args["versLabel"] = get_version_desc()

        self._timer = None
        self._command_counter = 0

    def set_argv(self, argv):
        """Store command line arguments.

        Useful for interactive execution with IPython that "executes and exits"
        to remember the arguments passed to the code_aster script.

        Arguments:
            argv (list[str]): Arguments passed to code_aster script.
        """
        self._argv = (argv or [])[:]

    def set_option(self, option, value):
        """Set the value of an execution parameter.

        Arguments:
            option (str): Name of the parameter.
            value (misc): Parameter value.
        """
        # Options must at least be declared by __init__
        if option in self._args:
            self._args[option] = value
            if option == "tpmax":
                libaster.set_option(option, value)
            elif option == "numthreads":
                os.environ["OMP_NUM_THREADS"] = str(value)
        elif value is not None:
            if value:
                self.enable(Options.by_name(option))
            else:
                self.disable(Options.by_name(option))

    def get_option(self, option):
        """Return the value of an execution parameter. If the *option* is
        unknown, it returns *None*.

        Arguments:
            option (str): Name of the parameter.

        Returns:
            misc: Parameter value.
        """
        logger.debug("get_option %r", option)
        if option.startswith("prog:"):
            value = get_program_path(re.sub("^prog:", "", option))
        elif option in self._args:
            value = self._args[option]
        else:
            try:
                value = self.option & Options.by_name(option)
            except AttributeError:
                value = None
        if isinstance(value, str):
            value = convert(value)
        logger.debug("return for %r: %s %s", option, value, type(value))
        return value

    def enable(self, option):
        """Enable a boolean option.

        Arguments:
            option (int): An 'Options' value.
        """
        self._bool |= option

        # for options that required an action
        if option & Options.Debug:
            logger.setLevel(DEBUG)
        if option & Options.ShowDeprecated:
            # disabled by default in python>=2.7
            warnings.simplefilter("default")
            warnings.filterwarnings(
                "error", message=".*EXTR_COMP.*removed", category=DeprecationWarning
            )

    def disable(self, option):
        """Disable a boolean option.

        Arguments:
            option (int): An 'Options' value.
        """
        self._bool = (self._bool | option) ^ option

        # for options that required an action
        if option & Options.Debug:
            logger.setLevel(INFO)
        if option & Options.ShowDeprecated:
            warnings.resetwarnings()

    @property
    def option(self):
        """Attribute that holds the boolean options.

        Returns:
            int: Bitwise combination of boolean execution parameters.
        """
        return self._bool

    @property
    def timer(self):
        """Attribute that holds the timer object.

        Returns:
            Timer: Timer object.
        """
        return self._timer

    @property
    def export(self):
        """Export: Attribute that holds the 'export' property."""
        return self._export

    def parse_args(self, argv=None):
        """Parse the command line arguments to set the execution parameters.

        Arguments:
            argv (list[str]): List of arguments inserted before those set with
                `set_argv()` or `sys.argv`.
        """
        # command arguments parser
        parser = ArgumentParser(
            description="execute a Code_Aster study", prog="Code_Aster{called by Python}"
        )
        parser.add_argument(
            "-g",
            "--debug",
            dest="Debug",
            action="store_const",
            const=1,
            default=0,
            help="add debug informations",
        )

        parser.add_argument(
            "--abort",
            dest="Abort",
            action="store_const",
            const=1,
            default=0,
            help="abort execution in case of error (testcase mode, by default "
            "raise an exception)",
        )
        parser.add_argument(
            "--test",
            dest="TestMode",
            action="store_const",
            const=1,
            default=0,
            help="set execution in testcase mode",
        )
        parser.add_argument(
            "--slave",
            dest="SlaveMode",
            action="store_const",
            const=1,
            default=0,
            help="slave mode, try not to exit in case of error",
        )
        parser.add_argument(
            "--interactive_interpreter",
            dest="InteractiveInterpreter",
            action="store_const",
            const=1,
            default=0,
            help="interactive execution, do not automatically try to close the "
            "database in case of exception",
        )

        parser.add_argument(
            "--hook_post_exec",
            action="store",
            default=None,
            help="[for debugging only] path to a function called by 'post_exec'",
        )

        parser.add_argument(
            "--dbgjeveux",
            action="store_const",
            const=1,
            default=0,
            help="turn on some additional checkings in the memory management",
        )

        parser.add_argument("--jxveri", action="store_const", const=1, default=0, help="")
        parser.add_argument("--sdveri", action="store_const", const=1, default=0, help="")
        parser.add_argument(
            "--impr_macro",
            dest="ShowChildCmd",
            action="store_const",
            const=1,
            default=0,
            help="show syntax of commands called by macro-commands",
        )

        parser.add_argument(
            "--memory",
            action=MemoryAction,
            type=float,
            default=DEFAULT_MEMORY_LIMIT,
            help="memory limit in MB used for code_aster objects "
            "(default: {0} MB)".format(DEFAULT_MEMORY_LIMIT),
        )
        parser.add_argument(
            "--memjeveux", dest="memory", action=MemoryAction, type=float, help=SUPPRESS
        )

        parser.add_argument(
            "--tpmax",
            action="store",
            type=float,
            default=DEFAULT_TIME_LIMIT,
            help="time limit of the execution in seconds "
            "(default: {0} s)".format(DEFAULT_TIME_LIMIT),
        )
        parser.add_argument(
            "--maxbase",
            action="store",
            type=float,
            default=None,
            help="size limit in MB for code_aster out-of-core files (glob.*, " "default: 2 TB)",
        )
        parser.add_argument(
            "--max_base", dest="maxbase", action="store", type=float, default=None, help=SUPPRESS
        )
        parser.add_argument(
            "--numthreads", action="store", type=int, default=1, help="maximum number of threads"
        )

        parser.add_argument(
            "--rcdir",
            dest="rcdir",
            action="store",
            type=str,
            metavar="DIR",
            default=RCDIR,
            help="directory containing resources (material properties, "
            "additional data files...). Defaults to {0}".format(RCDIR),
        )
        parser.add_argument(
            "--rep_mat", dest="repmat", action="store", metavar="DIR", default=".", help=SUPPRESS
        )
        parser.add_argument(
            "--rep_dex", dest="repdex", action="store", metavar="DIR", default=".", help=SUPPRESS
        )

        parser.add_argument(
            "--stage_number",
            action="store",
            type=int,
            default=1,
            metavar="NUM",
            help="Stage number in the Study",
        )
        parser.add_argument(
            "--start",
            dest="ForceStart",
            action="store_const",
            const=1,
            default=0,
            help="turn on to start a new calculation (force initialization with DEBUT)",
        )
        parser.add_argument(
            "--continue",
            dest="Continue",
            action="store_const",
            const=1,
            default=0,
            help="turn on to continue a previous execution (force initialization with POURSUITE)",
        )
        parser.add_argument(
            "--last",
            dest="LastStep",
            action="store_const",
            const=1,
            default=0,
            help="to be used for the last step of a study",
        )
        parser.add_argument(
            "--save_db",
            dest="SaveBase",
            action="store_const",
            const=1,
            default=0,
            help="to be enabled if the database will be saved",
        )
        parser.add_argument(
            "--use_legacy_mode",
            dest="UseLegacyMode",
            action="store",
            default=1,
            help="use (=1) or not (=0) the legacy mode for macro-commands " "results. (default: 1)",
        )

        parser.add_argument(
            "--deprecated",
            dest="ShowDeprecated",
            action="store_const",
            const=1,
            default=1,
            help="turn on deprecation warnings",
        )
        parser.add_argument(
            "--no-deprecated",
            dest="ShowDeprecated",
            action="store_const",
            const=0,
            help="turn off deprecation warnings (default: on)",
        )

        parser.add_argument(
            "--max_check",
            action="store",
            type=int,
            default=500,
            help="maximum number of occurrences to be checked, " "next are ignored",
        )
        parser.add_argument(
            "--max_print",
            action="store",
            type=int,
            default=500,
            help="maximum number of keywords or values printed in " "commands echo",
        )

        parser.add_argument(
            "--link",
            action="append",
            default=[],
            help="define a new link to an input or output file",
        )

        if not self._argv:
            self.set_argv(sys.argv)
        args, ignored = parser.parse_known_args(list(argv or []) + self._argv)

        logger.debug("Ignored arguments: %r", ignored)
        logger.debug("Read options: %r", vars(args))
        if "-max_base" in " ".join(ignored):
            deprecate("-max_base", case=4, help="Use '--max_base' instead.")

        # at most one of...
        if args.ForceStart and args.Continue:
            logger.warning("'--start' and '--continue' are incompatible, '--continue' is removed.")
            args.Continue = 0

        # default value with DEBUT
        if not args.Continue:
            if not args.maxbase:
                # use default from envima.c
                args.maxbase = -1
        else:
            # --max_base must not be changed with POURSUITE
            if args.maxbase:
                logger.warning("'--max_base' is ignored with POURSUITE.")
                args.maxbase = None

        # assign parameter values
        for opt, value in list(vars(args).items()):
            self.set_option(opt, value)

        self._timer = Timer(format="aster", limit=self._args["max_print"])

        # store Export object
        # TODO may be passed directly as argument
        self._export = Export(from_string=" ", check=False)
        for link in self._args["link"]:
            self._export.import_file_argument(link)

        # For convenience DEBUG can be set from environment
        if int(os.getenv("DEBUG", 0)) >= 1:
            self.enable(Options.Debug)

    def sub_tpmax(self, tsub):
        """Reduce the cpu time limit of `tsub`."""
        self.set_option("tpmax", self.get_option("tpmax") - tsub)

    def incr_command_counter(self):
        """Increment the counter of Commands.

        Returns:
            int: Current value of the counter.
        """
        self._command_counter += 1
        return self._command_counter

    # register objects callable from libaster
    def register_global_object(self, key, obj):
        """Register an object to be callable from libaster.

        Arguments:
            key (str): Access key.
            obj (*misc*): Registered object, usually a class or a function.
        """
        self._global[key] = obj

    @property
    def catalc(self):
        """Attribute that holds the catalog of behavior."""
        return self._global["catalc"]

    @property
    def logical_unit(self):
        """Attribute that holds the logical units manager."""
        return self._global["logical_unit"]

    @property
    def syntax(self):
        """Attribute that holds the command syntax class."""
        return self._global["syntax"]

    @property
    def print_header(self):
        """func: Attribute that holds the 'print_header' function."""
        return self._global["print_header"]

    @property
    def checksd(self):
        """func: Attribute that holds the 'checksd' property."""
        return self._global["checksd"]

    @property
    def testresu_print(self):
        """func: Attribute that holds the 'testresu_print' property."""
        return self._global["testresu_print"]

    @property
    def iniran(self):
        """func: Attribute that holds the 'iniran' function."""
        return self._global["iniran"]

    @property
    def getran(self):
        """func: Attribute that holds the 'getran' function."""
        return self._global["getran"]


class MemoryAction(Action):
    """Specific action to store the memory limit argument."""

    def __call__(self, parser, namespace, value, optstr):
        """Check and store the memory limit.

        Arguments:
            See :py:func:`argparse.Action`.
        """
        factor = 8 if "64" in platform.architecture()[0] else 4
        if optstr == "--memjeveux":
            value = value * factor
        namespace.memory = value


def get_program_path(program):
    """Return the path to *program* as stored by 'waf configure'.

    Returns:
        str: Path stored during configuration or *program* itself otherwise.
    """
    if getattr(get_program_path, "_cache", None) is None:
        prog_cfg = {}
        fname = osp.join(os.environ["ASTER_DATADIR"], "external_programs.yaml")
        if not osp.isfile(fname) or not yaml:
            fname = osp.splitext(fname)[0] + ".json"
        if osp.isfile(fname):
            with open(fname, "rb") as fcfg:
                if osp.splitext(fname)[-1] == ".yaml":
                    assert (
                        yaml
                    ), f"yaml not available, can not use {fname}, please convert it into '.json' instead"
                    prog_cfg = yaml.load(fcfg, Loader=yaml.Loader)
                else:
                    prog_cfg = json.load(fcfg)
        get_program_path._cache = prog_cfg

    programs = get_program_path._cache
    value = programs.get(program, shutil.which(program))
    if not value or not osp.isfile(value):
        value = program
    return value


@contextmanager
def disable_fpe():
    """Context manager to locally disable FPE interception."""
    libaster.matfpe(-1)
    yield
    libaster.matfpe(1)
