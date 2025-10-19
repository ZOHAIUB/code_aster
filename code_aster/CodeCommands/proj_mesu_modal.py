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

# person_in_charge: natacha.bereux@edf.fr

from ..Objects import (
    GeneralizedModeResult,
    HarmoGeneralizedResult,
    TransientGeneralizedResult,
    GeneralizedDOFNumbering,
)
from ..Supervis import ExecuteCommand


class ProjMesuModal(ExecuteCommand):
    """Command PROJ_MESU_MODAL"""

    command_name = "PROJ_MESU_MODAL"

    def create_result(self, keywords):
        """Initialize the result.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """
        mesure = keywords["MODELE_MESURE"]["MESURE"]
        if mesure.getType() == "DYNA_TRANS":
            self._result = TransientGeneralizedResult()
        elif mesure.getType() == "DYNA_HARMO":
            self._result = HarmoGeneralizedResult()
        elif mesure.getType() == "MODE_MECA":
            self._result = GeneralizedModeResult()
        elif mesure.getType() == "MODE_MECA_C":
            self._result = GeneralizedModeResult()
        else:
            raise TypeError("Type not allowed")

    def post_exec(self, keywords):
        """Execute the command.

        Arguments:
            keywords (dict): User's keywords.
        """
        self._result.setGeneralizedDOFNumbering(
            GeneralizedDOFNumbering(self._result.getName() + ".NUGEN")
        )
        base = keywords["MODELE_CALCUL"]["BASE"]
        if isinstance(self._result, TransientGeneralizedResult) or isinstance(
            self._result, HarmoGeneralizedResult
        ):
            self._result.setDOFNumbering(base.getDOFNumbering())
        else:
            self._result.setDOFNumbering(base.getDOFNumbering())
            self._result.setMesh(base.getMesh())
            self._result.build()


PROJ_MESU_MODAL = ProjMesuModal.run
