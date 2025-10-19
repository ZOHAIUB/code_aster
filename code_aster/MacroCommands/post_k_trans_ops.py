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

from ..Messages import UTMESS

from ..Cata.Syntax import _F
from ..CodeCommands import RECU_GENE, CREA_TABLE
from ..Utilities.misc import get_titre_concept


def post_k_trans_ops(self, **args):
    """
    Ecriture de la macro post_k_trans
    """
    EnumTypes = (list, tuple)

    RESU_TRANS = args.get("RESU_TRANS")
    K_MODAL = args.get("K_MODAL")
    NUME_ORDRE = args.get("NUME_ORDRE")
    LIST_ORDRE = args.get("LIST_ORDRE")
    INST = args.get("INST")
    LIST_INST = args.get("LIST_INST")

    # ------------------------------------------------------------------
    # On importe les definitions des commandes a utiliser dans la macro

    # ------------------------------------------------------------------
    TABK = K_MODAL["TABL_K_MODA"]

    __kgtheta = TABK

    tablin = __kgtheta.EXTR_TABLE()

    #  sif_arg = args['tablin']
    if "K3" in tablin.para:
        DIME = 3
    else:
        DIME = 2

    # -----------------------------------------
    #
    # Verification de cohérence sur le nombre de modes
    #
    # RESULTAT TRANSITOIRE
    nmodtr = RESU_TRANS.getNumberOfModes()
    # BASE MODALE
    if DIME == 2:
        n_mode = len((__kgtheta.EXTR_TABLE())["K1"])
        nbno = 1
    else:
        n_mode = max((__kgtheta.EXTR_TABLE())["NUME_MODE"].values()["NUME_MODE"])
        nbno = max((__kgtheta.EXTR_TABLE())["NUM_PT"].values()["NUM_PT"])
        labsc = (__kgtheta.EXTR_TABLE())["ABSC_CURV"].values()["ABSC_CURV"][0:nbno]
    if nmodtr != n_mode:
        n_mode = min(nmodtr, n_mode)
        UTMESS("A", "RUPTURE0_50", valk=RESU_TRANS.getName(), vali=n_mode)

    #
    # Traitement des mots clés ORDRE/INST/LIST_INST et LIST_ORDRE
    #
    l0_inst = RESU_TRANS.getTimes()
    l0_ord = RESU_TRANS.getIndexes()
    nbtrans = len(l0_ord)
    li = [[l0_ord[i], l0_inst[i]] for i in range(nbtrans)]
    ln = [[l0_ord[i], i] for i in range(nbtrans)]
    lo = [[l0_inst[i], l0_ord[i]] for i in range(nbtrans)]
    li = [(i[0], i[1:]) for i in li]
    ln = [(i[0], i[1:]) for i in ln]
    lo = [(i[0], i[1:]) for i in lo]
    d_ord = dict(lo)
    d_ins = dict(li)
    d_num = dict(ln)

    l_ord = []
    l_inst = []
    if LIST_ORDRE or NUME_ORDRE:
        if NUME_ORDRE:
            if type(NUME_ORDRE) not in EnumTypes:
                NUME_ORDRE = (NUME_ORDRE,)
            ltmp = list(NUME_ORDRE)
        elif LIST_ORDRE:
            ltmp = LIST_ORDRE.getValues()
        for ord in ltmp:
            if ord in l0_ord:
                l_ord.append(ord)
                l_inst.append(d_ins[ord][0])
            else:
                UTMESS("A", "RUPTURE0_51", vali=ord, valk=nomresu)
    elif LIST_INST or INST:
        CRITERE = args["CRITERE"]
        PRECISION = args["PRECISION"]
        if INST:
            if type(INST) not in EnumTypes:
                INST = (INST,)
            ltmp = list(INST)
        elif LIST_INST:
            ltmp = LIST_INST.getValues()
        for ins in ltmp:
            if CRITERE == "RELATIF" and ins != 0.0:
                match = [x for x in l0_inst if abs((ins - x) / ins) < PRECISION]
            else:
                match = [x for x in l0_inst if abs(ins - x) < PRECISION]
            if len(match) == 0:
                UTMESS("A", "RUPTURE0_38", valr=ins)
            elif len(match) >= 2:
                UTMESS("A", "RUPTURE0_39", valr=ins)
            else:
                l_inst.append(match[0])
                l_ord.append(d_ord[match[0]][0])
    else:
        l_ord = l0_ord
        l_inst = l0_inst
    nbarch = len(l_ord)
    if nbarch == 0:
        UTMESS("F", "RUPTURE0_54")

    #
    # Calcul des K(t)
    #
    K1mod = [None] * n_mode * nbno
    K2mod = [None] * n_mode * nbno
    K1t = [None] * nbarch * nbno
    K2t = [None] * nbarch * nbno
    if DIME == 3:
        K3mod = [None] * n_mode * nbno
        K3t = [None] * nbarch * nbno
        k1 = "K1"
        k2 = "K2"
        k3 = "K3"
    else:
        k1 = "K1"
        k2 = "K2"

    for x in range(0, nbno):
        for k in range(0, n_mode):
            K1mod[k * nbno + x] = __kgtheta[k1, k * nbno + x + 1]
            K2mod[k * nbno + x] = __kgtheta[k2, k * nbno + x + 1]
            if DIME == 3:
                K3mod[k * nbno + x] = __kgtheta[k3, k * nbno + x + 1]

        for num in range(0, nbarch):
            K1t[num * nbno + x] = 0.0
            K2t[num * nbno + x] = 0.0
            if DIME == 3:
                K3t[num * nbno + x] = 0.0
            vect_gen = RECU_GENE(RESU_GENE=RESU_TRANS, INST=l_inst[num], NOM_CHAM="DEPL")
            coef = vect_gen.getValues()
            for k in range(0, n_mode):
                num_ord = d_num[l_ord[num]][0]
                alpha = coef[k]
                K1t[num * nbno + x] = K1t[num * nbno + x] + alpha * K1mod[k * nbno + x]
                K2t[num * nbno + x] = K2t[num * nbno + x] + alpha * K2mod[k * nbno + x]
                if DIME == 3:
                    K3t[num * nbno + x] = K3t[num * nbno + x] + alpha * K3mod[k * nbno + x]

    titre = get_titre_concept()
    if DIME == 2:
        tabout = CREA_TABLE(
            LISTE=(
                _F(LISTE_I=l_ord, PARA="NUME_ORDRE"),
                _F(LISTE_R=l_inst, PARA="INST"),
                _F(LISTE_R=K1t, PARA=k1),
                _F(LISTE_R=K2t, PARA=k2),
            ),
            TITRE=titre,
        )
    else:
        lo = []
        li = []
        for i in range(nbarch):
            for j in range(nbno):
                lo.append(l_ord[i])
                li.append(l_inst[i])
        tabout = CREA_TABLE(
            LISTE=(
                _F(LISTE_I=lo, PARA="NUME_ORDRE"),
                _F(LISTE_R=li, PARA="INST"),
                _F(LISTE_I=list(range(nbno)) * nbarch, PARA="NUM_PT"),
                _F(LISTE_R=labsc * nbarch, PARA="ABSC_CURV"),
                _F(LISTE_R=K1t, PARA=k1),
                _F(LISTE_R=K2t, PARA=k2),
                _F(LISTE_R=K3t, PARA=k3),
            ),
            TITRE=titre,
        )

    # ------------------------------------------------------------------
    return tabout
