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

from ..Objects import FullTransientResult
from ..Supervis import ExecuteCommand


class RestCondTran(ExecuteCommand):
    """Command REST_COND_TRAN"""

    command_name = "REST_COND_TRAN"

    def create_result(self, keywords):
        """Initialize the result.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """
        reuse = keywords.get("reuse")
        if reuse is not None:
            self._result = reuse
        else:
            if keywords["TYPE_RESU"] == "DYNA_TRANS":
                self._result = FullTransientResult()
            else:
                self._result = type(keywords["RESULTAT"])()

    def post_exec(self, keywords):
        """Execute the command.

        Arguments:
            keywords (dict): User's keywords.
        """
        if keywords.get("BASE_MODALE"):
            modeMeca = keywords["BASE_MODALE"]
        else:
            modeMeca = keywords["MACR_ELEM_DYNA"].getMechanicalMode()
        dofNum = modeMeca.getDOFNumbering()
        fnds = []
        if dofNum:
            if keywords["TYPE_RESU"] == "DYNA_TRANS":
                self._result.setDOFNumbering(dofNum)
            else:
                self._result.setModel(dofNum.getModel())
                for i in dofNum.getFiniteElementDescriptors():
                    self._result.addFiniteElementDescriptor(i)
                fnds.append(dofNum.getEquationNumbering())

        self._result.build(fnds=fnds)


REST_COND_TRAN = RestCondTran.run
