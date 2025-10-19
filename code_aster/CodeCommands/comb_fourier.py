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

# person_in_charge: nicolas.sellenet@edf.fr

from ..Objects import CombinedFourierResult
from ..Supervis import ExecuteCommand


class CombFourier(ExecuteCommand):
    """Command that creates the :class:`~code_aster.Objects.CombinedFourierResult`"""

    command_name = "COMB_FOURIER"

    def create_result(self, keywords):
        """Initialize the result.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """
        self._result = CombinedFourierResult()

    def post_exec(self, keywords):
        """Execute the command.

        Arguments:
            keywords (dict): User's keywords.
        """
        model = keywords["RESULTAT"].getModel()
        if model is not None:
            self._result.setModel(model)
        else:
            self._result.setMesh(keywords["RESULTAT"].getMesh())

        self._result.build(
            keywords["RESULTAT"].getFiniteElementDescriptors(),
            keywords["RESULTAT"].getEquationNumberings(),
        )


COMB_FOURIER = CombFourier.run
