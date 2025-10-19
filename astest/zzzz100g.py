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

# person_in_charge: mathieu.courtois at edf.fr
import numpy as NP

from code_aster.Cata.Commons import *
from code_aster.Cata.DataStructure import *
from code_aster.Cata.Syntax import *
from code_aster.Supervis.ExecuteCommand import UserMacro


def user_func(X, vx, vy):
    return NP.interp(X, vx, vy)


def macro_test_ops(self, VALE, **kwargs):
    """Macro de test."""
    from code_aster.Commands import FORMULE

    vale_x = NP.arange(5.0)
    vale_y = NP.arange(5.0) * 2 * VALE

    result = FORMULE(
        NOM_PARA="X",
        VALE="user_func(X, vale_x, vale_y)",
        user_func=user_func,
        vale_x=vale_x,
        vale_y=vale_y,
    )

    return result


MACRO_TEST_CATA = FORM(
    nom="MACRO_TEST",
    op=macro_test_ops,
    sd_prod=formule,
    fr=tr("Macro de test"),
    regles=(UN_PARMI("VALE", "VALE_C"),),
    VALE=SIMP(statut="o", typ="R"),
)

MACRO_TEST = UserMacro("MACRO_TEST", MACRO_TEST_CATA, macro_test_ops)
