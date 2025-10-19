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
from ..CodeCommands import CALC_CHAMP, CALC_TABLE, CREA_TABLE, MACR_LIGN_COUPE
from ..Objects.table_py import Table
from ..Messages import UTMESS, MasquerAlarme, RetablirAlarme


def post_coque_ops(self, RESULTAT, COOR_POINT, CHAM, NUME_ORDRE=None, INST=None, **args):
    """
    macro post_coque
    """

    # On importe les definitions des commandes a utiliser dans la macro

    MasquerAlarme("MODELISA4_9")

    assert RESULTAT.getType() in ("EVOL_ELAS", "EVOL_NOLI")
    dico = aster.GetResu(RESULTAT.getName(), "CHAMPS")
    dico2 = aster.GetResu(RESULTAT.getName(), "VARI_ACCES")
    # si ni INST ni NUME_ORDRE ne sont presents, on prend le premier
    # instant calcule
    if not INST and not NUME_ORDRE:
        INST = dico2["INST"][0]

    if NUME_ORDRE:
        if NUME_ORDRE not in dico2["NUME_ORDRE"]:
            UTMESS("F", "POST0_25", vali=NUME_ORDRE)
    else:
        if INST not in dico2["INST"]:
            UTMESS("F", "POST0_26", valr=INST)
    #
    if NUME_ORDRE:
        if CHAM == "EFFORT":
            if NUME_ORDRE not in dico["EFGE_ELNO"]:
                if NUME_ORDRE in dico["DEPL"]:
                    CALC_CHAMP(
                        RESULTAT=RESULTAT,
                        reuse=RESULTAT,
                        CONTRAINTE="EFGE_ELNO",
                        NUME_ORDRE=NUME_ORDRE,
                    )
                else:
                    UTMESS("F", "POST0_19", vali=NUME_ORDRE)
        elif CHAM == "DEFORMATION":
            if NUME_ORDRE not in dico["DEGE_ELNO"]:
                if NUME_ORDRE in dico["DEPL"]:
                    CALC_CHAMP(
                        RESULTAT=RESULTAT,
                        reuse=RESULTAT,
                        DEFORMATION="DEGE_ELNO",
                        NUME_ORDRE=NUME_ORDRE,
                    )
                else:
                    UTMESS("F", "POST0_19", vali=NUME_ORDRE)

    dico = aster.GetResu(RESULTAT.getName(), "CHAMPS")

    # Appel MACR_LIGN_COUPE :
    motscles = {}
    if CHAM == "EFFORT":
        motscles["NOM_CHAM"] = "EFGE_ELNO"
    if CHAM == "DEFORMATION":
        motscles["NOM_CHAM"] = "DEGE_ELNO"

    if CHAM == "EFFORT":
        motscles["LIGN_COUPE"] = []
        iocc = 0
        for m in COOR_POINT:
            iocc = iocc + 1
            lst = m["COOR"]
            if len(lst) == 4 and lst[3] != 0.0:
                UTMESS("A", "POST0_21", vali=iocc, valr=lst[3])
            lst = lst[0:3]
            motscles["LIGN_COUPE"].append(
                _F(TYPE="SEGMENT", NB_POINTS=2, COOR_ORIG=lst, COOR_EXTR=lst, DISTANCE_MAX=10.0)
            )
        __tabl = MACR_LIGN_COUPE(RESULTAT=RESULTAT, **motscles)

    if CHAM == "DEFORMATION":
        motscles["LIGN_COUPE"] = []
        iocc = 0
        for m in COOR_POINT:
            iocc = iocc + 1
            lst = m["COOR"]
            if len(lst) != 4:
                UTMESS("F", "POST0_22", vali=iocc)
            else:
                lst = lst[0:3]
                motscles["LIGN_COUPE"].append(
                    _F(TYPE="SEGMENT", NB_POINTS=2, COOR_ORIG=lst, COOR_EXTR=lst, DISTANCE_MAX=10.0)
                )
        __tabl = MACR_LIGN_COUPE(RESULTAT=RESULTAT, **motscles)

    tab2 = __tabl.EXTR_TABLE()
    if NUME_ORDRE:
        tab3 = tab2.NUME_ORDRE == NUME_ORDRE
    else:
        tab3 = tab2.INST == INST
    tab2 = tab3

    tab4 = Table()
    ilig = 0
    for ligne in tab2:
        ilig = ilig + 1
        if (ilig % 2) == 0:
            tab4.append(ligne)
    tab4 = tab4[tab2.para]
    #
    #  on cree une table(dege) bidon qu'on va surcharger
    #
    if CHAM == "DEFORMATION":
        motscles["NOM_CHAM"] = "DEGE_ELNO"
        motscles["LIGN_COUPE"] = []
        tabz = []
        iocc = 0
        for m in COOR_POINT:
            iocc = iocc + 1
            lst = m["COOR"]
            z = lst[3]
            tabz.append(z)
            lst = lst[0:3]
            motscles["LIGN_COUPE"].append(
                _F(TYPE="SEGMENT", NB_POINTS=2, COOR_ORIG=lst, COOR_EXTR=lst, DISTANCE_MAX=10.0)
            )
        __tabeps = MACR_LIGN_COUPE(RESULTAT=RESULTAT, **motscles)
        __teps = CALC_TABLE(
            TABLE=__tabeps,
            ACTION=(
                _F(OPERATION="RENOMME", NOM_PARA=("EXX", "EPXX")),
                _F(OPERATION="RENOMME", NOM_PARA=("EYY", "EPYY")),
                _F(OPERATION="RENOMME", NOM_PARA=("EXY", "EPZZ")),
                _F(OPERATION="RENOMME", NOM_PARA=("KXX", "EPXY")),
                _F(OPERATION="RENOMME", NOM_PARA=("KYY", "EPXZ")),
                _F(OPERATION="RENOMME", NOM_PARA=("KXY", "EPYZ")),
                _F(
                    OPERATION="EXTR",
                    NOM_PARA=(
                        "INTITULE",
                        "NOM_CHAM",
                        "NUME_ORDRE",
                        "INST",
                        "ABSC_CURV",
                        "COOR_X",
                        "COOR_Y",
                        "COOR_Z",
                        "EPXX",
                        "EPYY",
                        "EPZZ",
                        "EPXY",
                        "EPXZ",
                        "EPYZ",
                    ),
                ),
            ),
        )

        tabep2 = __teps.EXTR_TABLE()
        if NUME_ORDRE:
            tabep3 = tabep2.NUME_ORDRE == NUME_ORDRE
        else:
            tabep3 = tabep2.INST == INST
        tabep2 = tabep3

        tabep4 = Table()
        ilig = 0
        for ligne in tabep2:
            ilig = ilig + 1
            if (ilig % 2) == 0:
                tabep4.append(ligne)
        tabep4 = tabep4[tabep2.para]

        iligout = 0
        for ligout in tabep4:
            iligout = iligout + 1
            iligin = 0
            for ligin in tab4:
                iligin = iligin + 1
                if iligout == iligin:
                    ligout["EPXX"] = ligin["EXX"] + ligin["KXX"] * tabz[iligout - 1]
                    ligout["EPYY"] = ligin["EYY"] + ligin["KYY"] * tabz[iligout - 1]
                    ligout["EPXY"] = ligin["EXY"] + ligin["KXY"] * tabz[iligout - 1]
                    ligout["EPZZ"] = 0.0
                    ligout["EPXZ"] = ligin["GAX"] * 0.5
                    ligout["EPYZ"] = ligin["GAY"] * 0.5

    if CHAM == "EFFORT":
        dprod = tab4.dict_CREA_TABLE()
    elif CHAM == "DEFORMATION":
        dprod = tabep4.dict_CREA_TABLE()

    tabout = CREA_TABLE(TYPE_TABLE="TABLE", **dprod)
    RetablirAlarme("MODELISA4_9")
    return tabout
