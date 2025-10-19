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

"""
:py:class:`FullTransientResult` --- results result
****************************************************
"""

from libaster import FullResult

from ..Utilities import injector
from ..ObjectsExt.result_ext import ResultStateBuilder


class FullResultStateBuilder(ResultStateBuilder):
    """Class that returns the internal state of a *FullResult* to be pickled."""

    def save(self, result):
        """Return the internal state of a *FullResult* to be pickled.

        Arguments:
            result (*FullResult*): The *FullResult* object to be pickled.

        Returns:
            *InternalStateBuilder*: The internal state itself.
        """
        # numbering
        self._st["numbering"] = result.getDOFNumbering()
        super().save(result)

        return self

    def restore(self, result):
        """Restore the *FullResult* content from the previously saved internal
        state.

        Arguments:
            result (*FullResult*): The *DataStructure* object to be restored.
        """
        result.setDOFNumbering(self._st["numbering"])
        super().restore(result)


@injector(FullResult)
class ExtentedFullResult:
    """Object for FullResult."""

    cata_sdj = "SD.sd_resultat.sd_resultat"
    internalStateBuilder = FullResultStateBuilder
