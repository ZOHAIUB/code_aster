# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
:py:class:`ConstantFieldOnCellsReal` --- Results container
**************************************************
"""

import aster
from libaster import ConstantFieldOnCellsReal, ConstantFieldOnCellsChar16

from ..Utilities import injector


@injector(ConstantFieldOnCellsReal)
class ExtendedConstantFieldOnCellsReal:
    cata_sdj = "SD.sd_carte.sd_carte"

    def __getinitargs__(self):
        """Returns the argument required to reinitialize a ConstantFieldOnCellsReal
        object during unpickling.
        """
        return (self.getName(), self.getMesh())


@injector(ConstantFieldOnCellsChar16)
class ExtendedConstantFieldOnCellsChar16:
    cata_sdj = "SD.sd_carte.sd_carte"

    def __getinitargs__(self):
        """Returns the argument required to reinitialize a ConstantFieldOnCellsChar16
        object during unpickling.
        """
        return (self.getName(), self.getMesh())
