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

from ..Objects import (
    AssemblyMatrixDisplacementComplex,
    AssemblyMatrixDisplacementReal,
    AssemblyMatrixPressureComplex,
    AssemblyMatrixPressureReal,
    AssemblyMatrixTemperatureComplex,
    AssemblyMatrixTemperatureReal,
    GeneralizedAssemblyMatrixComplex,
    GeneralizedAssemblyMatrixReal,
)
from ..Supervis import ExecuteCommand


class MatrixCombination(ExecuteCommand):
    """Command that creates the
    :class:`~code_aster.Objects.AssemblyMatrixDisplacementReal`
    :class:`~code_aster.Objects.AssemblyMatrixDisplacementComplex`
    :class:`~code_aster.Objects.GeneralizedAssemblyMatrixReal`
    :class:`~code_aster.Objects.GeneralizedAssemblyMatrixComplex`
    """

    command_name = "COMB_MATR_ASSE"

    def create_result(self, keywords):
        """Initialize the result.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """
        self._result = None
        if "COMB_C" in keywords:
            if len(keywords["COMB_C"]) > 0:
                if "MATR_ASSE" in keywords["COMB_C"][0]:
                    matrix = keywords["COMB_C"][0]["MATR_ASSE"]
                    if isinstance(matrix, AssemblyMatrixDisplacementReal) or isinstance(
                        matrix, AssemblyMatrixDisplacementComplex
                    ):
                        self._result = AssemblyMatrixDisplacementComplex()
                    if isinstance(matrix, GeneralizedAssemblyMatrixReal) or isinstance(
                        matrix, GeneralizedAssemblyMatrixComplex
                    ):
                        self._result = GeneralizedAssemblyMatrixComplex()
                    if isinstance(matrix, AssemblyMatrixTemperatureReal) or isinstance(
                        matrix, AssemblyMatrixTemperatureComplex
                    ):
                        self._result = AssemblyMatrixTemperatureComplex()
                    if isinstance(matrix, AssemblyMatrixPressureReal) or isinstance(
                        matrix, AssemblyMatrixPressureComplex
                    ):
                        self._result = AssemblyMatrixPressureComplex()
        elif "COMB_R" in keywords:
            if len(keywords["COMB_R"]) > 0:
                if "MATR_ASSE" in keywords["COMB_R"][0]:
                    matrix = keywords["COMB_R"][0]["MATR_ASSE"]
                    if isinstance(matrix, AssemblyMatrixDisplacementReal) or isinstance(
                        matrix, AssemblyMatrixDisplacementComplex
                    ):
                        self._result = AssemblyMatrixDisplacementReal()
                    if isinstance(matrix, GeneralizedAssemblyMatrixReal) or isinstance(
                        matrix, GeneralizedAssemblyMatrixComplex
                    ):
                        self._result = GeneralizedAssemblyMatrixReal()
                    if isinstance(matrix, AssemblyMatrixTemperatureReal) or isinstance(
                        matrix, AssemblyMatrixTemperatureComplex
                    ):
                        self._result = AssemblyMatrixTemperatureReal()
                    if isinstance(matrix, AssemblyMatrixPressureReal) or isinstance(
                        matrix, AssemblyMatrixPressureComplex
                    ):
                        self._result = AssemblyMatrixPressureReal()
        elif "CALC_AMOR_GENE" in keywords:
            self._result = GeneralizedAssemblyMatrixReal()
        if self._result is None:
            raise NotImplementedError()

    def post_exec(self, keywords):
        """Execute the command.

        Arguments:
            keywords (dict): User's keywords.
        """
        comb = keywords.get("COMB_R")
        if comb is None:
            comb = keywords.get("COMB_C")
        if comb is not None:
            if isinstance(self._result, GeneralizedAssemblyMatrixReal) or isinstance(
                self._result, GeneralizedAssemblyMatrixComplex
            ):
                dofNum = comb[0]["MATR_ASSE"].getGeneralizedDOFNumbering()
                self._result.setGeneralizedDOFNumbering(dofNum)
            else:
                dofNum = comb[0]["MATR_ASSE"].getDOFNumbering()
                self._result.setDOFNumbering(dofNum)


COMB_MATR_ASSE = MatrixCombination.run
