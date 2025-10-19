# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
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

# person_in_charge: francesco.bettonte at edf.fr

from ..Commons import *
from ..Language.DataStructure import *
from ..Language.Syntax import *


def calc_homo_sdprod(self, CORR_MECA, CORR_THER, **kwargs):
    if kwargs.get("__all__"):
        return ([table_sdaster],)

    if CORR_THER is not None:
        self.type_sdprod(CORR_THER, evol_ther_dict)
    if CORR_MECA is not None:
        self.type_sdprod(CORR_MECA, evol_elas_dict)

    return table_sdaster


CALC_MATE_HOMO = MACRO(
    nom="CALC_MATE_HOMO",
    op=OPS("code_aster.MacroCommands.MateHomo.mate_homo_ops.mate_homo_ops"),
    sd_prod=calc_homo_sdprod,
    reentrant="n",
    fr=tr("Calcul des paramètres elastiques équivalents par homogénéisation périodique"),
    MAILLAGE=SIMP(statut="o", typ=maillage_sdaster),
    AFFE=FACT(
        statut="o",
        max="**",
        regles=UN_PARMI("TOUT", "GROUP_MA"),
        TOUT=SIMP(statut="f", typ="TXM", into=("OUI",)),
        GROUP_MA=SIMP(statut="f", typ=grma, validators=NoRepeat(), max="**"),
        MATER=SIMP(statut="o", typ=mater_sdaster, max=1),
    ),
    VARC=FACT(
        statut="o",
        max=1,
        NOM_VARC=SIMP(statut="o", typ="TXM", into=("TEMP", "IRRA")),
        VALE=SIMP(statut="o", typ="R", min=1, max="**"),
    ),
    TYPE_HOMO=SIMP(
        statut="o", typ="TXM", into=("MASSIF", "PLAQUE", "PLAQUE_CT_MINDLIN", "PLAQUE_CT_TOURATIER")
    ),
    INFO=SIMP(statut="f", typ="I", defaut=1, into=(1, 2)),
    CORR_MECA=SIMP(statut="f", typ=CO, validators=NoRepeat()),
    b_corr_massif=BLOC(
        condition="""equal_to("TYPE_HOMO", "MASSIF")""",
        CORR_THER=SIMP(statut="f", typ=CO, validators=NoRepeat()),
    ),
)
