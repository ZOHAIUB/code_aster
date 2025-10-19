# coding: utf-8

# Copyright (C) 1991 - 2023  EDF R&D                www.code-aster.org
#
# This file is part of Code_Aster.
#
# Code_Aster is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Code_Aster is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Code_Aster.  If not, see <http://www.gnu.org/licenses/>.

# person_in_charge: nicolas.tardieu@edf.fr

from ..Objects import XfemCrack
from ..Supervis import ExecuteCommand


class XFEMCrackDefinition(ExecuteCommand):
    """Command that defines :class:`~code_aster.Objects.XfemCrack`."""

    command_name = "DEFI_FISS_XFEM"

    def create_result(self, keywords):
        """Initialize the result.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """
        self._result = XfemCrack(keywords["MAILLAGE"])

    def post_exec(self, keywords):
        if "TYPE_DISCONTINUITE" in keywords:
            if keywords["TYPE_DISCONTINUITE"] == "INTERFACE":
                self._result.setDiscontinuityType("Interface")
            elif keywords["TYPE_DISCONTINUITE"] == "COHESIF":
                self._result.setDiscontinuityType("Cohesive")
        self._result.updateInternalState()


DEFI_FISS_XFEM = XFEMCrackDefinition.run
