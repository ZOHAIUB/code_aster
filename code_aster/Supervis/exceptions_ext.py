# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

# person_in_charge: mathieu.courtois@edf.fr
"""
Definition of code_aster exceptions
***********************************

The exceptions classes are created from C++.
The module :py:mod:`code_aster.Supervis.exceptions_ext` defines convenient
methods for formatting and adds them to the base class :py:class:`code_aster.AsterError`.

The class hierarchy for code_aster exceptions is:

.. code-block:: none

    AsterError
     |
     +-- ConvergenceError
     |      Raised in case of convergence problem.
     |
     +-- IntegrationError
     |      Raised during local integration (of the behavior for example).
     |
     +-- SolverError
     |      Raised during solving phases.
     |
     +-- ContactError
     |      Raised during contact algorithms.
     |
     +-- TimeLimitError
            Raised when the time limit is reached.

The arguments of *AsterError* exceptions are the *message identifier (idmess)*
and, optionally and in this order, *list of strings (valk)*,
*list of ints (vali)* and *list of floats (valr)*.
"""

from libaster import (
    AsterError,
    ContactError,
    ConvergenceError,
    IntegrationError,
    SolverError,
    TimeLimitError,
)

from ..Messages import message_exception
from ..Utilities import convert


def format(exc, code):
    """Return the related message formatted according to `code`.

    Arguments:
        exc (AsterError): Exception raised.
        code (str): Level of the message (one of 'AIESF', see `Utmess`).

    Returns;
        str: Message of the exception.
    """
    try:
        txt = message_exception(code, exc)
    except:
        txt = str(exc.args)
    return convert(txt)


def format_exception(exc):
    """Return the text to show an exception.

    Arguments:
        exc (AsterError): Exception raised.

    Returns;
        str: String representation of the exception.
    """
    return format(exc, "Z")


def get_idmess(exc):
    """Return the message identifier of the exception.

    Arguments:
        exc (AsterError): Exception raised.

    Returns;
        str: Message identifier (example: ``CALCUL_12``).
    """
    return len(exc.args) >= 1 and exc.args[0] or ""


AsterError.__repr__ = format_exception
AsterError.__str__ = format_exception
AsterError.format = format
AsterError.id_message = property(get_idmess)
AsterError.message = property(format_exception)
