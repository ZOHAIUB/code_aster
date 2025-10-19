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

from ..Objects import GeneralizedAssemblyMatrixComplex, GeneralizedAssemblyMatrixReal
from ..Supervis import ExecuteCommand


class GeneralizedMatrixBuilder(ExecuteCommand):
    """Command that creates the :class:`~code_aster.Objects.GeneralizedAssemblyMatrixReal`
    and `~code_aster.Objects.GeneralizedAssemblyMatrixComplex`"""

    command_name = "ASSE_MATR_GENE"

    def create_result(self, keywords):
        """Initialize the result.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """

        if keywords["METHODE"] == "INITIAL":
            self._result = GeneralizedAssemblyMatrixReal()
        elif keywords["OPTION"] == "RIGI_GENE_C":
            self._result = GeneralizedAssemblyMatrixComplex()
        else:
            self._result = GeneralizedAssemblyMatrixReal()

    def post_exec(self, keywords):
        """Store references to ElementaryMatrix objects.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords, changed
                in place.
        """
        self._result.setGeneralizedDOFNumbering(keywords["NUME_DDL_GENE"])


ASSE_MATR_GENE = GeneralizedMatrixBuilder.run
