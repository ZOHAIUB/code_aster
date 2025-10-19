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

from ..Cata.Syntax import _F
from ..CodeCommands import (
    ASSE_MATRICE,
    ASSE_VECTEUR,
    CALC_CHAMP,
    CALC_MATR_ELEM,
    CALC_VECT_ELEM,
    CREA_RESU,
    FACTORISER,
)
from ..CodeCommands import NUME_DDL as NUME_DDL_CMD
from ..CodeCommands import RESOUDRE
from ..Messages import UTMESS
from ..Objects import PhysicalProblem, DiscreteComputation


def macro_elas_mult_ops(
    self,
    MODELE,
    CAS_CHARGE,
    CHAM_MATER=None,
    CARA_ELEM=None,
    NUME_DDL=None,
    CHAR_MECA_GLOBAL=None,
    LIAISON_DISCRET=None,
    SOLVEUR=None,
    **args
):
    """
    Ecriture de la macro MACRO_ELAS_MULT
    """
    args = _F(args)

    if args["reuse"] is not None:
        changed = args.get("RESULTAT")
        if changed is None:
            UTMESS("F", "SUPERVIS2_79", valk="RESULTAT")

    ielas = 0
    ifour = 0
    for m in CAS_CHARGE:
        if m["NOM_CAS"]:
            ielas = 1  # mot clé NOM_CAS      présent sous CAS_CHARGE
            tyresu = "MULT_ELAS"
        else:
            ifour = 1  # mot clé MODE_FOURIER présent sous CAS_CHARGE
            tyresu = "FOURIER_ELAS"
    if ielas == 1 and ifour == 1:
        UTMESS("F", "ELASMULT0_1")

    if NUME_DDL in self.sdprods or NUME_DDL is None:
        # Si le concept NUME_DDL est dans self.sdprods ou n est pas nommé
        # il doit etre  produit par la macro
        # il faudra donc appeler la commande NUME_DDL
        lnume = 1
    else:
        lnume = 0

    if ielas == 1:
        motscles = {}
        if CHAR_MECA_GLOBAL:
            motscles["CHARGE"] = CHAR_MECA_GLOBAL
        if CHAM_MATER:
            motscles["CHAM_MATER"] = CHAM_MATER
        if CARA_ELEM:
            motscles["CARA_ELEM"] = CARA_ELEM
        __nomrig = CALC_MATR_ELEM(OPTION="RIGI_MECA", MODELE=MODELE, **motscles)

        if lnume:
            # On peut passer des mots cles egaux a None. Ils sont ignores
            motscles = {}
            if NUME_DDL is not None:
                num = NUME_DDL_CMD(MATR_RIGI=__nomrig, **motscles)
                self.register_result(num, NUME_DDL)
            else:
                _num = NUME_DDL_CMD(MATR_RIGI=__nomrig, **motscles)
                num = _num
        else:
            num = NUME_DDL

        __nomras = ASSE_MATRICE(MATR_ELEM=__nomrig, NUME_DDL=num)

        if SOLVEUR["METHODE"] == "MUMPS":
            __nomras = FACTORISER(
                reuse=__nomras,
                MATR_ASSE=__nomras,
                NPREC=SOLVEUR["NPREC"],
                STOP_SINGULIER=SOLVEUR["STOP_SINGULIER"],
                METHODE="MUMPS",
                PCENT_PIVOT=SOLVEUR["PCENT_PIVOT"],
                RENUM=SOLVEUR["RENUM"],
                TYPE_RESOL=SOLVEUR["TYPE_RESOL"],
                PRETRAITEMENTS=SOLVEUR["PRETRAITEMENTS"],
                ELIM_LAGR=SOLVEUR["ELIM_LAGR"],
                GESTION_MEMOIRE=SOLVEUR["GESTION_MEMOIRE"],
                ACCELERATION=SOLVEUR["ACCELERATION"],
                LOW_RANK_SEUIL=SOLVEUR["LOW_RANK_SEUIL"],
            )
        else:
            __nomras = FACTORISER(
                reuse=__nomras,
                MATR_ASSE=__nomras,
                NPREC=SOLVEUR["NPREC"],
                STOP_SINGULIER=SOLVEUR["STOP_SINGULIER"],
                METHODE=SOLVEUR["METHODE"],
                RENUM=SOLVEUR["RENUM"],
            )
    #
    # boucle sur les items de CAS_CHARGE

    nomchn = []
    lcharg = []
    iocc = 0
    for m in CAS_CHARGE:
        iocc = iocc + 1

        if "SOUS_TITRE" not in m:
            m["SOUS_TITRE"] = None
        if "VECT_ASSE" not in m:
            m["VECT_ASSE"] = None

        # calcul de lcharg : liste des listes de char_meca (mots clé CHAR_MECA
        # et CHAR_MECA_GLOBAL)
        xx1 = m["CHAR_MECA"]
        if type(xx1) not in (tuple, list):
            xx1 = (xx1,)
        xx2 = CHAR_MECA_GLOBAL
        if type(xx2) not in (tuple, list):
            xx2 = (xx2,)
        lchar1 = []
        for chargt in tuple(xx1) + tuple(xx2):
            if chargt:
                lchar1.append(chargt)
        lcharg.append(lchar1)
        assert len(lchar1) > 0

        if ifour:
            motscles = {}
            if CHAR_MECA_GLOBAL:
                motscles["CHARGE"] = CHAR_MECA_GLOBAL
            if CHAM_MATER:
                motscles["CHAM_MATER"] = CHAM_MATER
            if CARA_ELEM:
                motscles["CARA_ELEM"] = CARA_ELEM
            motscles["MODE_FOURIER"] = m["MODE_FOURIER"]
            __nomrig = CALC_MATR_ELEM(OPTION="RIGI_MECA", MODELE=MODELE, **motscles)

            if lnume:
                _num = NUME_DDL_CMD(MATR_RIGI=__nomrig)
                num = _num
                lnume = 0

            __nomras = ASSE_MATRICE(MATR_ELEM=__nomrig, NUME_DDL=num)

            if SOLVEUR["METHODE"] == "MUMPS":
                __nomras = FACTORISER(
                    reuse=__nomras,
                    MATR_ASSE=__nomras,
                    NPREC=SOLVEUR["NPREC"],
                    STOP_SINGULIER=SOLVEUR["STOP_SINGULIER"],
                    METHODE="MUMPS",
                    PCENT_PIVOT=SOLVEUR["PCENT_PIVOT"],
                    RENUM=SOLVEUR["RENUM"],
                    TYPE_RESOL=SOLVEUR["TYPE_RESOL"],
                    PRETRAITEMENTS=SOLVEUR["PRETRAITEMENTS"],
                    ELIM_LAGR=SOLVEUR["ELIM_LAGR"],
                    GESTION_MEMOIRE=SOLVEUR["GESTION_MEMOIRE"],
                    ACCELERATION=SOLVEUR["ACCELERATION"],
                    LOW_RANK_SEUIL=SOLVEUR["LOW_RANK_SEUIL"],
                )
            else:
                __nomras = FACTORISER(
                    reuse=__nomras,
                    MATR_ASSE=__nomras,
                    NPREC=SOLVEUR["NPREC"],
                    STOP_SINGULIER=SOLVEUR["STOP_SINGULIER"],
                    METHODE=SOLVEUR["METHODE"],
                    RENUM=SOLVEUR["RENUM"],
                )

        if m["VECT_ASSE"] is None:
            motscles = {}
            l_calc_varc = False
            if CHAM_MATER:
                motscles["CHAM_MATER"] = CHAM_MATER
                if CHAM_MATER.hasExternalStateVariable():
                    l_calc_varc = True
            if CARA_ELEM:
                motscles["CARA_ELEM"] = CARA_ELEM
            if ifour:
                motscles["MODE_FOURIER"] = m["MODE_FOURIER"]
            if len(lchar1) > 0:
                motscles["CHARGE"] = lchar1
            __nomvel = CALC_VECT_ELEM(OPTION="CHAR_MECA", **motscles)

            # chargement du aux variables de commandes
            if l_calc_varc:
                __nomasv = ASSE_VECTEUR(VECT_ELEM=__nomvel, NUME_DDL=num)

                __phys_pb = PhysicalProblem(MODELE, CHAM_MATER)
                __phys_pb.setDOFNumbering(num)
                __phys_pb.computeReferenceExternalStateVariables()
                __varc = __phys_pb.getExternalStateVariables(1.0)

                nh = 0
                if "MODE_FOURIER" in motscles:
                    nh = motscles["MODE_FOURIER"]

                __nomasv += DiscreteComputation(__phys_pb).getExternalStateVariablesForces(
                    1.0, __varc, mode=nh
                )

            else:
                __nomasv = ASSE_VECTEUR(VECT_ELEM=(__nomvel,), NUME_DDL=num)
        else:
            __nomasv = m["VECT_ASSE"]

        if SOLVEUR["METHODE"] == "MUMPS":
            __nomchn = RESOUDRE(
                MATR=__nomras,
                CHAM_NO=__nomasv,
                TITRE=m["SOUS_TITRE"],
                RESI_RELA=SOLVEUR["RESI_RELA"],
                POSTTRAITEMENTS=SOLVEUR["POSTTRAITEMENTS"],
            )
        else:
            __nomchn = RESOUDRE(MATR=__nomras, CHAM_NO=__nomasv, TITRE=m["SOUS_TITRE"])

        nomchn.append(__nomchn)

    # fin de la boucle sur les items de CAS_CHARGE
    #

    motscles = {}
    iocc = 0
    motscle2 = {}
    if CHAM_MATER:
        motscle2["CHAM_MATER"] = CHAM_MATER
    if CARA_ELEM:
        motscle2["CARA_ELEM"] = CARA_ELEM
    if ielas:
        motscles["AFFE"] = []
        for m in CAS_CHARGE:
            if len(lcharg[iocc]) > 0:
                motscles["AFFE"].append(
                    _F(
                        MODELE=MODELE,
                        NOM_CHAM="DEPL",
                        CHAM_GD=nomchn[iocc],
                        NOM_CAS=m["NOM_CAS"],
                        CHARGE=lcharg[iocc],
                        **motscle2
                    )
                )
            else:
                motscles["AFFE"].append(
                    _F(
                        NOM_CHAM="DEPL",
                        MODELE=MODELE,
                        CHAM_GD=nomchn[iocc],
                        NOM_CAS=m["NOM_CAS"],
                        **motscle2
                    )
                )
            iocc = iocc + 1
    else:
        motscles["AFFE"] = []
        for m in CAS_CHARGE:
            if len(lcharg[iocc]) > 0:
                motscles["AFFE"].append(
                    _F(
                        MODELE=MODELE,
                        NOM_CHAM="DEPL",
                        CHAM_GD=nomchn[iocc],
                        NUME_MODE=m["MODE_FOURIER"],
                        TYPE_MODE=m["TYPE_MODE"],
                        CHARGE=lcharg[iocc],
                        **motscle2
                    )
                )
            else:
                motscles["AFFE"].append(
                    _F(
                        MODELE=MODELE,
                        NOM_CHAM="DEPL",
                        CHAM_GD=nomchn[iocc],
                        NUME_MODE=m["MODE_FOURIER"],
                        TYPE_MODE=m["TYPE_MODE"],
                        **motscle2
                    )
                )
            iocc = iocc + 1

    if args.get("reuse"):
        motscles["reuse"] = args.get("reuse")
        motscles["RESULTAT"] = args.get("reuse")
    nomres = CREA_RESU(OPERATION="AFFE", TYPE_RESU=tyresu, **motscles)

    #
    # boucle sur les items de CAS_CHARGE pour SIEF_ELGA

    iocc = 0
    for m in CAS_CHARGE:
        iocc = iocc + 1

        if m["OPTION"] == "SIEF_ELGA":
            motscles = {}
            if ielas:
                motscles["NOM_CAS"] = m["NOM_CAS"]
            else:
                motscles["NUME_MODE"] = m["MODE_FOURIER"]
            CALC_CHAMP(reuse=nomres, RESULTAT=nomres, CONTRAINTE="SIEF_ELGA", **motscles)

    # fin de la boucle sur les items de CAS_CHARGE
    #
    return nomres
