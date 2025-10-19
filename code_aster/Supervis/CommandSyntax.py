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

"""
:py:mod:`CommandSyntax` --- Interface between C++/Fortran operators and user syntax
***********************************************************************************

The :py:class:`CommandSyntax` provides an interface between operators originally
defined in Fortran (:file:`opXXXX` subroutines) and the user syntax.
The C/Fortran interface is defined in :file:`aster_module.c`.

It is also used in some C++ objects to emulate a user command.
The C++ interface is available from C++ ``CommandSyntax`` class.

The method :py:meth:`~CommandSyntax.define` stores the user keywords.
The execution of an operator (by :py:func:`libaster.call_oper` for example) registers
the current :py:class:`CommandSyntax` instance that can be requested by the
operator using the functions:

    - :py:meth:`~CommandSyntax.getres`,

    - :py:meth:`~CommandSyntax.getvtx`,

    - :py:meth:`~CommandSyntax.getltx`,

    - :py:meth:`~CommandSyntax.getvid`,

    - :py:meth:`~CommandSyntax.getvis`,

    - :py:meth:`~CommandSyntax.getvr8`,

    - :py:meth:`~CommandSyntax.getvc8`,

    - :py:meth:`~CommandSyntax.getfac`,

    - :py:meth:`~CommandSyntax.getexm`,

    - :py:meth:`~CommandSyntax.getmjm`.

"""

from ..Cata import Commands
from ..Cata.SyntaxChecker import CheckerError, SyntaxCheckerVisitor
from ..Objects import DataStructure
from ..Utilities import (
    force_list,
    import_object,
    is_complex,
    is_float,
    is_int,
    is_str,
    logger,
    value_is_sequence,
)
from ..Utilities.outputs import command_text
from .typeaster import typeaster

# WARNING:
#   user keywords dict may be very big, uncomment 'logger.debug' lines
#   only during debugging


class CommandSyntax:
    """This class describes the syntax of command for compatibility
    with the fortran code that use the legacy supervisor.
    """

    _currentCommand = None

    @classmethod
    def setCurrentCommand(cls, command):
        """Register current command.

        Arguments:
            command (*CommandSyntax*): *CommandSyntax* to register.
        """
        # The object is registered with `register_sh_etape()` by `aster_oper`
        # if called from Python or by the constructor of `CommandSyntax` if
        # called from C++.
        cls._currentCommand = command

    @classmethod
    def getCurrentCommand(cls):
        """Return the current command.

        Returns:
            *CommandSyntax*: Current registered *CommandSyntax*.
        """
        return cls._currentCommand

    def __init__(self, name="unNamed", cata=None):
        """Create a new command or part of syntax.

        Arguments:
            name (str): Command name (ex: "DEBUT") or path to catalog (ex:
                "code_aster.Cata.Commands.DEBUT").
            cata (*PartOfSyntax*): Command description.
        """
        self._name = name.split(".")[-1]
        self._resultName = " "
        self._resultType = " "
        self._resultValue = None
        self._definition = None
        logger.debug("Syntax: new command is %r", name)
        currentCommand = self.getCurrentCommand()
        # only FIN is allowed to free the current "in failure" command
        if self._name == "FIN" and currentCommand is not None:
            currentCommand.free()
            currentCommand = self.getCurrentCommand()
        assert currentCommand is None, "CommandSyntax {} must be freed".format(
            currentCommand.getName()
        )
        self.setCurrentCommand(self)
        if not cata:
            if "." not in name:
                cata = getattr(Commands, name, None)
            else:
                cata = import_object(name)
            if not cata:
                logger.debug("CommandSyntax: catalog not found for %r", name)
        self._commandCata = cata

    def free(self):
        """Reset the current command pointer as soon as possible"""
        # `currentCommand` must be reset before the garbage collector will do it
        # logger.debug("Syntax: del command %r", name)
        self.setCurrentCommand(None)

    def __repr__(self):
        """Representation of the command"""
        return "Command {!r}, returns {!r} <{!r}>\n` syntax: {}".format(
            self._name, self._resultName, self._resultType, self._definition
        )

    def debugPrint(self):
        """Representation of the command"""
        logger.debug("%r", self)

    def setResult(self, sdName, sdType):
        """Register the result of the command: name and type.

        Arguments:
            sdName (str): Name of the result created by the Command.
            sdType (str): Type of the result created by the Command.
        """
        self._resultName = sdName
        self._resultType = sdType

    def define(self, dictSyntax, add_default=True, check_syntax=False):
        """Register the keywords values.

        Arguments:
            dictSyntax (dict): User keywords.
            add_default (bool, optional): Tell if default keywords have to be
                added or not.
            check_syntax (bool, optional): Check that the keywords validity.
                Only when it has not already been checked by ExecuteCommand object.
        """
        if self._commandCata is not None and add_default:
            # logger.debug("define0 %r: %r", self._name, dictSyntax)
            self._commandCata.addDefaultKeywords(dictSyntax)
        if check_syntax:
            checker = SyntaxCheckerVisitor()
            try:
                self._commandCata.accept(checker, dictSyntax)
            except (CheckerError, KeyError, TypeError, ValueError) as exc:
                text = command_text(self._name, dictSyntax)
                logger.error("Invalid keywords:\n  Error: %s\n  Keywords:\n%s", exc, text)
        self._definition = dictSyntax
        # logger.debug("define1 %r: %r", self._name, self._definition)

    def getName(self):
        """Return the command name.

        Returns:
            str: Command name.
        """
        return self._name

    def getResultName(self):
        """Return the name of the result of the Command.

        Returns:
            str: Name of the result of the Command.
        """
        return self._resultName

    def getResultType(self):
        """Return the type of the result of the command.

        Returns:
            str: Type name of the result of the Command.
        """
        return self._resultType

    def getResultValue(self):
        """Return the value of the result of the Command.

        Returns:
            str: Name of the result of the Command.
        """
        return self._resultValue

    def _getFactorKeyword(self, factName):
        """Return the occurrences of a factor keyword.

        Arguments:
            factName (str): Name of the factor keyword.

        Returns:
            list: List of occurences of a factor keyword, *None* if it is not
                provided by the user.
        """
        # a factor keyword may be empty: {} (None means 'does not exist')
        dictDef = self._definition.get(factName, None)
        # logger.debug("factor keyword %r: %r", factName, dictDef)
        if dictDef is None:
            return None
        if isinstance(dictDef, dict):
            dictDef = [dictDef]
        return dictDef

    def _getFactorKeywordOccurrence(self, factName, occurrence):
        """Return the definition of an occurrence of a factor keyword.

        Arguments:
            factName (str): Name of the factor keyword.
            occurrence (int): Index of the occurrence (start from 0).

        Returns:
            dict: Occurrence of a factor keyword, *None* if it is not
                provided by the user."""
        dictDef = self._getFactorKeyword(factName)
        if dictDef is None:
            return None
        try:
            return dictDef[occurrence]
        except Exception:
            return None

    def _getDefinition(self, factName, occurrence):
        """Return the definition of a factor keyword or of the top-level
        if `factName` is blank.

        Arguments:
            factName (str): Name of the factor keyword.
            occurrence (int): Index of the occurrence (start from 0).

        Returns:
            dict: User keywords under a factor keyword or the Command.
        """
        if not factName.strip():
            dictDef = self._definition
        else:
            dictDef = self._getFactorKeywordOccurrence(factName, occurrence)
            if not dictDef:
                dictDef = {}
        assert isinstance(dictDef, dict), "syntax not defined"
        return dictDef

    def _getCataDefinition(self, factName):
        """Return the definition of a factor keyword in the catalog or of the
        top-level if `factName` is blank.

        Arguments:
            factName (str): Name of the factor keyword.

        Returns:
            *CataDefinition*: Definition of the factor keyword or the Command.
        """
        if not self._commandCata:
            logger.debug("CommandSyntax: catalog is not available")
            return None
        factName = factName.strip()
        catadef = self._commandCata.definition
        if not factName:
            return catadef
        keywords = catadef.factor_keywords.get(factName)
        if not keywords:
            return None
        return keywords.definition

    def getFactorKeywordNbOcc(self, factName):
        """Return the number of occurrences of a factor keyword in the
        user's keywords.

        Arguments:
            factName (str): Name of the factor keyword.

        Returns:
            int: Number of occurrences of a factor keyword.
        """
        dictDef = self._getFactorKeyword(factName)
        # logger.debug("_getFactorKeyword %r: %r", factName, dictDef)
        if dictDef is None:
            return 0
        # logger.debug("getFactorKeywordNbOcc: len(dictDef) = %d", len(dictDef))
        return len(dictDef)

    getfac = getFactorKeywordNbOcc

    def existsFactorAndSimpleKeyword(self, factName, occurrence, simpName):
        """Tell if the couple ( factor keyword, simple keyword ) exists in the
        user keywords.

        Arguments:
            factName (str): Name of the factor keyword.
            occurrence (int): Index of the occurrence (start from 0).
            simpName (str): Name of the simple keyword.

        Returns:
            int: 1 if the keyword exists, else 0.
        """
        dictDef = self._getDefinition(factName, occurrence)
        value = dictDef.get(simpName, None)
        if value is None or isinstance(value, dict):
            return 0
        return 1

    def getValue(self, factName, occurrence, simpName):
        """Return the values of a (simple) keyword.

        Arguments:
            factName (str): Name of the factor keyword.
            occurrence (int): Index of the occurrence (start from 0).
            simpName (str): Name of the simple keyword.

        Returns:
            list[misc]: List of the values provided by the user.
        """
        if not self.existsFactorAndSimpleKeyword(factName, occurrence, simpName):
            return []
        value = self._getDefinition(factName, occurrence)[simpName]
        value = force_list(value)
        # logger.debug("getValue: %s %s %s", factName, simpName, value)
        return value

    def getltx(self, factName, simpName, occurrence, maxval, lenmax):
        """Wrapper function to return length of strings.

        Arguments:
            factName (str): Name of the factor keyword.
            simpName (str): Name of the simple keyword.
            occurrence (int): Index of the occurrence (start from 0).
            maxval (int): Maximum number of values read.
            lenmax (int): Maximum of lengths.

        Returns:
            list[misc]: List of the values provided by the user.
        """
        value = self.getValue(factName, occurrence, simpName)
        value = _check_strings(factName, simpName, value)
        size = len(value)
        if size > maxval:
            size = -size
        length = [min(len(i), lenmax) for i in value]
        return size, tuple(length[:maxval])

    def getvid(self, factName, simpName, occurrence, maxval):
        """Wrapper function to return a list of results.

        Arguments:
            factName (str): Name of the factor keyword.
            occurrence (int): Index of the occurrence (start from 0).
            simpName (str): Name of the simple keyword.
            maxval (int): Maximum number of values read.

        Returns:
            int, list: Returns two values ``(size, values)``.
            ``size`` is the number of the values provided by the user.
            If ``size > maxval``, ``-size`` is returned.
            ``values`` is a list of result names.
        """
        value = self.getValue(factName, occurrence, simpName)
        value = [i.getName() if hasattr(i, "getName") else i for i in value]
        size = len(value)
        if size > maxval:
            size = -size
        return size, tuple(value[:maxval])

    def getvtx(self, factName, simpName, occurrence, maxval):
        """Wrapper function to return a list of strings.

        Arguments:
            factName (str): Name of the factor keyword.
            occurrence (int): Index of the occurrence (start from 0).
            simpName (str): Name of the simple keyword.
            maxval (int): Maximum number of values read.

        Returns:
            int, list: Returns two values ``(size, values)``.
            ``size`` is the number of the values provided by the user.
            If ``size > maxval``, ``-size`` is returned.
            ``values`` is a list of strings.
        """
        value = self.getValue(factName, occurrence, simpName)
        value = _check_strings(factName, simpName, value)
        size = len(value)
        if size > maxval:
            size = -size
        return size, tuple(value[:maxval])

    def getvis(self, factName, simpName, occurrence, maxval):
        """Wrapper function to return a list of integers.

        Arguments:
            factName (str): Name of the factor keyword.
            occurrence (int): Index of the occurrence (start from 0).
            simpName (str): Name of the simple keyword.
            maxval (int): Maximum number of values read.

        Returns:
            int, list: Returns two values ``(size, values)``.
            ``size`` is the number of the values provided by the user.
            If ``size > maxval``, ``-size`` is returned.
            ``values`` is a list of integers.
        """
        value = self.getValue(factName, occurrence, simpName)
        if len(value) > 0 and not is_int(value[0], onvalue=True):
            raise TypeError("integer expected, got %s" % type(value[0]))
        size = len(value)
        if size > maxval:
            size = -size
        return size, tuple([round(i) for i in value[:maxval]])

    def getvr8(self, factName, simpName, occurrence, maxval):
        """Wrapper function to return a list of float numbers.

        Arguments:
            factName (str): Name of the factor keyword.
            occurrence (int): Index of the occurrence (start from 0).
            simpName (str): Name of the simple keyword.
            maxval (int): Maximum number of values read.

        Returns:
            int, list: Returns two values ``(size, values)``.
            ``size`` is the number of the values provided by the user.
            If ``size > maxval``, ``-size`` is returned.
            ``values`` is a list of floats.
        """
        value = self.getValue(factName, occurrence, simpName)
        if len(value) > 0:
            try:
                float(value[0])
            except TypeError:
                raise TypeError("float expected, got %s" % type(value[0]))
        size = len(value)
        if size > maxval:
            size = -size
        return size, tuple(value[:maxval])

    def getvc8(self, factName, simpName, occurrence, maxval):
        """Wrapper function to return a list of complex numbers.

        Arguments:
            factName (str): Name of the factor keyword.
            occurrence (int): Index of the occurrence (start from 0).
            simpName (str): Name of the simple keyword.
            maxval (int): Maximum number of values read.

        Returns:
            int, list: Returns two values ``(size, values)``.
            ``size`` is the number of the values provided by the user.
            If ``size > maxval``, ``-size`` is returned.
            ``values`` is a list of complex numbers.
        """
        values = self.getValue(factName, occurrence, simpName)
        if len(values) > 0:
            if isinstance(values[0], str):
                values = (values,)
            toReturn = []
            for value in values:
                if type(value) in (list, tuple):
                    assert value[0] in ("RI", "MP")
                    toReturn.append(tuple(value))
                else:
                    toReturn.append(value)

            size = len(toReturn)
            if size > maxval:
                size = -size
            return size, tuple(toReturn[:maxval])
        size = len(values)
        if size > maxval:
            size = -size
        return size, tuple(values[:maxval])

    def getvpy(self, factName, simpName, occurrence):
        """Wrapper function to return a list of Python objects.

        Arguments:
            factName (str): Name of the factor keyword.
            occurrence (int): Index of the occurrence (start from 0).
            simpName (str): Name of the simple keyword.

        Returns:
            list: Tuple containing the objects provided by the user.
        """
        try:
            values = tuple(self.getValue(factName, occurrence, simpName))
        except Exception as exc:
            print("DEBUG: exception:", str(exc))
            values = ()
        return len(values), values

    def getres(self):
        """Return the name and type of the result, and the command name.

        Returns:
            str: Name of the Command result (name of the jeveux object).
            str: Type name of the result.
            str: Command name.
        """
        jev, typ, cmd = self.getResultName(), self.getResultType(), self.getName()
        # logger.debug("Command %s: result name %r, type %r", cmd, jev, typ)
        return jev, typ, cmd

    def setres(self, value):
        """Define a value for special commands that returns a builtin type
        (*int* or *float*).

        Arguments:
            value (int|float): Value returned.
        """
        self._resultValue = value

    def getexm(self, factName, simpName):
        """Tell if the couple ( factor keyword, simple keyword ) exists in the
        Command catalog. It supposes that all blocs conditions are verified.

        Arguments:
            factName (str): Name of the factor keyword.
            simpName (str): Name of the simple keyword.

        Returns:
            int: 1 if the keyword exists, else 0.
        """
        catadef = self._getCataDefinition(factName)
        # logger.debug("getexm: %s / %s", factName, simpName)
        if not catadef:
            return 0
        if not simpName.strip():
            return 1
        keywords = catadef.simple_keywords
        # logger.debug("getexm: simple keywords: %s", list(keywords.keys()))
        return int(keywords.get(simpName) is not None)

    def getmjm(self, factName, occurrence, maxval):
        """Return the list of simple keywords provided by the user under a
        factor keyword.

        Arguments:
            factName (str): Name of the factor keyword.
            occurrence (int): Index of the occurrence (start from 0).
            maxval (int): Maximum number of values returned.

        Returns:
            list[str], list[str]: Two lists: one for the names of simple
            keywords, alphabetically sorted, and one for the type names
            of the keywords.
        """
        userkw = self._getDefinition(factName, occurrence)
        if not userkw:
            return (), ()
        # logger.debug("getmjm: user keywords: %s", userkw)
        catadef = self._getCataDefinition(factName).simple_keywords
        lkeywords = sorted(userkw.keys())
        kws, types = [], []
        for kw in lkeywords:
            obj = userkw[kw]
            if obj is None:
                continue
            # ignore factor keyword: legacy getmjm returned typ='MCList'
            if kw not in catadef:
                continue
            kws.append(kw)
            typ = typeaster(catadef[kw].definition["typ"])
            if value_is_sequence(obj) and not is_complex(obj):
                obj = obj[0]
            if is_complex(obj) and typ == "C8":
                pass
            elif is_float(obj) and typ in ("R8", "C8"):
                pass
            elif is_int(obj, onvalue=True) and typ in ("IS", "R8", "C8"):
                pass
            elif is_str(obj) and typ == "TX":
                pass
            elif is_str(obj) and "FORMULE" in typ:
                typ = typ[0]
            elif isinstance(obj, DataStructure):
                typ = obj.getType()
            else:
                raise TypeError("unsupported type: {0!r} {1}".format(obj, type(obj)))
            types.append(typ)
        return kws, types

    def getmat(self):
        raise NotImplementedError("'getmat' is not yet supported.")

    def getoper(self):
        """Return the operator number.

        Returns:
            int: Number of the fortran operator subroutine.
        """
        return self._commandCata.definition["op"]


def _check_strings(factName, simpName, value):
    """Check that keyword values are strings.

    Arguments:
        value (str): Values extracted from a keyword.

    Returns:
        list[str]: String values or names for DataStructure objects.
    """
    if len(value) > 0 and not is_str(value[0]):
        try:
            value2 = []
            for i in range(len(value)):
                value2.append(value[i].getName())
            value = value2
        except AttributeError:
            raise TypeError(
                "string expected for {0}/{1}, got {2}".format(factName, simpName, type(value[0]))
            )
    return value
