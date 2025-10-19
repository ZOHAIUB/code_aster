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


import aster
from ..Cata.Syntax import _F
from ..CodeCommands import (
    AFFE_MATERIAU,
    AFFE_MODELE,
    CREA_CHAMP,
    CREA_RESU,
    CREA_TABLE,
    DEFI_GROUP,
    DEFI_MATERIAU,
    IMPR_TABLE,
    POST_ELEM,
    POST_RELEVE_T,
)


def fctZeroUn(listIN):
    """Fonction qui renvoie une liste de 0/1 en fonction du signe des éléments
    de la liste listIN:"""
    listOUT = []
    for n in listIN:
        if n > 0:
            listOUT.append(1)
        else:
            listOUT.append(0)
    return listOUT


def post_decollement_ops(self, RESULTAT, NOM_CHAM, NOM_CMP, GROUP_MA, INFO, **args):
    """
    Corps de la macro POST_DECOLLEMENT
    """

    # On importe les definitions des commandes a utiliser dans la macro

    # on recupere le concept maillage
    MAILLAGE = RESULTAT.getModel().getMesh()

    # Creation du groupe de noeuds 'PDECOL'
    DEFI_GROUP(
        reuse=MAILLAGE,
        MAILLAGE=MAILLAGE,
        DETR_GROUP_NO=_F(NOM="PDECOL"),
        CREA_GROUP_NO=_F(GROUP_MA=GROUP_MA, NOM="PDECOL"),
        ALARME="NON",
    )

    # model restreint au GROUP_MA
    __model = AFFE_MODELE(
        MAILLAGE=MAILLAGE, AFFE=_F(GROUP_MA=GROUP_MA, PHENOMENE="MECANIQUE", MODELISATION="3D")
    )

    # Calcul de la surface du GROUP_MA : __surf
    __unit = CREA_CHAMP(
        OPERATION="AFFE",
        TYPE_CHAM="NOEU_NEUT_R",
        MODELE=__model,
        AFFE=_F(GROUP_NO="PDECOL", NOM_CMP="X1", VALE=1.0),
    )

    __chpg0 = CREA_CHAMP(
        PROL_ZERO="OUI", MODELE=__model, OPERATION="DISC", TYPE_CHAM="ELGA_NEUT_R", CHAM_GD=__unit
    )

    __mater0 = DEFI_MATERIAU(ELAS=_F(E=210000000.0, NU=0.3))

    __chmat0 = AFFE_MATERIAU(MODELE=__model, MAILLAGE=MAILLAGE, AFFE=_F(TOUT="OUI", MATER=__mater0))

    __resu0 = CREA_RESU(
        OPERATION="AFFE",
        TYPE_RESU="EVOL_ELAS",
        AFFE=_F(
            NOM_CHAM="VARI_ELGA", CHAM_MATER=__chmat0, MODELE=__model, CHAM_GD=__chpg0, INST=0.0
        ),
    )

    __tbSurf0 = POST_ELEM(
        RESULTAT=__resu0,
        INST=0.0,
        MODELE=__model,
        INTEGRALE=_F(NOM_CHAM="VARI_ELGA", NOM_CMP="X1", GROUP_MA=GROUP_MA, TYPE_MAILLE="2D"),
    )

    __surf = __tbSurf0.EXTR_TABLE().values()["INTE_X1"][0]

    __linst = aster.GetResu(RESULTAT.getName(), "VARI_ACCES")["INST"]

    # Calcul de la surface des noeuds décollés
    __pct = []

    for inst in __linst:
        __dep = CREA_CHAMP(
            OPERATION="EXTR",
            RESULTAT=RESULTAT,
            TYPE_CHAM="NOEU_" + NOM_CHAM[:4] + "_R",
            INST=inst,
            NOM_CHAM=NOM_CHAM,
        )

        __tb1 = POST_RELEVE_T(
            ACTION=_F(
                OPERATION="EXTRACTION",
                GROUP_NO="PDECOL",
                INTITULE=GROUP_MA,
                CHAM_GD=__dep,
                NOM_CMP=NOM_CMP,
            )
        )

        __col = fctZeroUn(__tb1.EXTR_TABLE().values()[NOM_CMP])

        __tb2 = CREA_TABLE(
            LISTE=(
                _F(LISTE_K=__tb1.EXTR_TABLE().values()["NOEUD"], PARA="NOEUD"),
                _F(LISTE_R=__col, PARA="X1"),
            )
        )

        __ch = CREA_CHAMP(OPERATION="EXTR", TYPE_CHAM="NOEU_NEUT_R", TABLE=__tb2, MAILLAGE=MAILLAGE)

        __chg = CREA_CHAMP(
            MODELE=__model, OPERATION="DISC", TYPE_CHAM="ELGA_NEUT_R", PROL_ZERO="OUI", CHAM_GD=__ch
        )

        __resu = CREA_RESU(
            OPERATION="AFFE",
            TYPE_RESU="EVOL_ELAS",
            AFFE=_F(
                NOM_CHAM="VARI_ELGA", CHAM_MATER=__chmat0, MODELE=__model, CHAM_GD=__chg, INST=0.0
            ),
        )

        __tb3 = POST_ELEM(
            RESULTAT=__resu,
            INST=0.0,
            MODELE=__model,
            INTEGRALE=_F(NOM_CHAM="VARI_ELGA", NOM_CMP="X1", GROUP_MA=GROUP_MA, TYPE_MAILLE="2D"),
        )

        __su2 = __tb3.EXTR_TABLE().values()["INTE_X1"][0]

        __pct.append(100.0 * __su2 / __surf)

    C_out = CREA_TABLE(LISTE=(_F(LISTE_R=__linst, PARA="INST"), _F(LISTE_R=__pct, PARA="%DECOL")))
    if INFO > 1:
        IMPR_TABLE(UNITE=6, TABLE=C_out)
    return C_out
