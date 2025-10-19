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

# person_in_charge: nicolas.sellenet@edf.fr

from ..Supervis import ExecuteCommand


class RestReduitComplet(ExecuteCommand):
    """Command REST_REDUIT_COMPLET"""

    command_name = "REST_REDUIT_COMPLET"

    def create_result(self, keywords):
        """Initialize the result.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """
        self._result = type(keywords["RESULTAT_REDUIT"])()

    def post_exec(self, keywords):
        """Execute the command.

        Arguments:
            keywords (dict): User's keywords.
        """
        self._result.setModel(keywords["MODELE"])
        try:
            materField = keywords["RESULTAT_REDUIT"].getMaterialField()
            if materField:
                self._result.setMaterialField(materField)
        except RuntimeError:
            pass
        self._result.build()


REST_REDUIT_COMPLET = RestReduitComplet.run
