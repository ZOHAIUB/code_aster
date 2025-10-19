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

# person_in_charge: mickael.abbas at edf.fr

from ..Commons import *
from ..Language.DataStructure import *
from ..Language.Syntax import *

AFFE_CHAR_SECH_F = OPER(
    nom="AFFE_CHAR_SECH_F",
    op=34,
    sd_prod=char_ther,
    fr=tr(
        "Affectation de charges et conditions aux limites pour le séchage fonction d'un (ou plusieurs)"
        " paramètres (temps, ...)"
    ),
    reentrant="n",
    regles=(AU_MOINS_UN("SECH_IMPO", "FLUX_NL"),),
    MODELE=SIMP(statut="o", typ=(modele_sdaster)),
    DOUBLE_LAGRANGE=SIMP(statut="f", typ="TXM", into=("OUI", "NON"), defaut="OUI"),
    SECH_IMPO=FACT(
        statut="f",
        max="**",
        regles=(AU_MOINS_UN("TOUT", "GROUP_MA", "GROUP_NO"),),
        TOUT=SIMP(statut="f", typ="TXM", into=("OUI",)),
        GROUP_NO=SIMP(statut="f", typ=grno, validators=NoRepeat(), max="**"),
        GROUP_MA=SIMP(statut="f", typ=grma, validators=NoRepeat(), max="**"),
        SANS_GROUP_MA=SIMP(statut="f", typ=grma, validators=NoRepeat(), max="**"),
        SANS_GROUP_NO=SIMP(statut="f", typ=grno, validators=NoRepeat(), max="**"),
        SECH=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
    ),
    FLUX_NL=FACT(
        statut="f",
        max="**",
        regles=(UN_PARMI("TOUT", "GROUP_MA", TOUT="OUI"),),
        TOUT=SIMP(statut="f", typ="TXM", into=("OUI",)),
        GROUP_MA=SIMP(statut="f", typ=grma, validators=NoRepeat(), max="**"),
        FLUN=SIMP(statut="o", typ=(fonction_sdaster,)),
    ),
    INFO=SIMP(statut="f", typ="I", defaut=1, into=(1, 2)),
    translation={"AFFE_CHAR_SECH_F": "Assign variable drying load"},
)
