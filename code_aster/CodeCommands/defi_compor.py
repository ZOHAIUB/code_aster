# coding: utf-8

# Copyright (C) 1991 - 2022  EDF R&D                www.code-aster.org
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

# person_in_charge: mathieu.courtois@edf.fr

from ..Objects import BehaviourDefinition
from ..Supervis import ExecuteCommand


class BehaviourDefinitionOperator(ExecuteCommand):
    """Execute legacy operator DEFI_COMPOR."""

    command_name = "DEFI_COMPOR"

    def create_result(self, keywords):
        """Create the result.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """
        self._result = BehaviourDefinition()


DEFI_COMPOR = BehaviourDefinitionOperator.run
