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

import math

from ..Messages import UTMESS

from ..Cata.Syntax import _F
from ..CodeCommands import (
    AFFE_CARA_ELEM,
    AFFE_CHAR_MECA,
    AFFE_MATERIAU,
    AFFE_MODELE,
    CALC_CHAMP,
    CALC_TABLE,
    CREA_CHAMP,
    CREA_RESU,
    DEFI_FONCTION,
    IMPR_RESU,
    LIRE_MAILLAGE,
    MODI_MAILLAGE,
    MODI_REPERE,
    POST_RELEVE_T,
    STAT_NON_LINE,
)
from ..Helpers import FileAccess, LogicalUnitFile
from ..Utilities import is_sequence
from .Utils.calc_point_mat import CALC_POINT_MAT


def simu_point_mat_ops(
    self,
    MATER,
    INCREMENT=None,
    SIGM_IMPOSE=None,
    EPSI_IMPOSE=None,
    SIGM_INIT=None,
    EPSI_INIT=None,
    VARI_INIT=None,
    NEWTON=None,
    CONVERGENCE=None,
    MASSIF=None,
    ANGLE=None,
    COMPORTEMENT=None,
    INFO=None,
    ARCHIVAGE=None,
    SUPPORT=None,
    **args
):
    """Simulation de la reponse d'un point materiel"""

    # On importe les definitions des commandes a utiliser dans la macro
    # Le nom de la variable doit etre obligatoirement le nom de la commande

    # -- Tests de cohérence
    __fonczero = DEFI_FONCTION(
        NOM_PARA="INST", VALE=(0, 0, 10, 0), PROL_DROITE="CONSTANT", PROL_GAUCHE="CONSTANT"
    )
    EPS = {}
    SIG = {}
    itetra = 0
    CMP_EPS = ["EPXX", "EPYY", "EPZZ", "EPXY", "EPXZ", "EPYZ"]
    CMP_SIG = ["SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ"]

    if COMPORTEMENT:
        lcomp = COMPORTEMENT[0]

    if SUPPORT is not None:
        if SUPPORT == "ELEMENT":
            itetra = 1
    if itetra == 0:
        if lcomp["DEFORMATION"] != "PETIT":
            if "GRAD_IMPOSE" in args:
                if args["GRAD_IMPOSE"] is not None:
                    if lcomp["DEFORMATION"] != "SIMO_MIEHE":
                        UTMESS("F", "COMPOR2_22", valk=lcomp["DEFORMATION"])
                    itetra = 0
                else:
                    itetra = 1
                    UTMESS("A", "COMPOR2_1", valk=lcomp["DEFORMATION"])

    # ===============================================================
    # cas ou il n'y a pas d'élement fini : appel à CALC_POINT_MAT
    # ===============================================================
    if itetra == 0:
        isig = 0
        ieps = 0
        igrd = 0
        ic1c2 = 0
        if SIGM_IMPOSE:
            SIG = dict(SIGM_IMPOSE[0])
            isig = 1
        if EPSI_IMPOSE:
            EPS = dict(EPSI_IMPOSE[0])
            ieps = 1
        if "GRAD_IMPOSE" in args:
            if args["GRAD_IMPOSE"] is not None:
                FIJ = dict(args["GRAD_IMPOSE"][0])
                igrd = 1
        if "MATR_C1" in args:
            if args["MATR_C1"] is not None:
                ic1c2 = 1
        if "MATR_C2" in args:
            if args["MATR_C2"] is not None:
                ic1c2 = 1

        motscles = {}
        if igrd:
            for i in list(FIJ.keys()):
                motscles[i] = FIJ[i]
        elif ic1c2:
            if "MATR_C1" in args:
                if args["MATR_C1"] is not None:
                    motscles["MATR_C1"] = args["MATR_C1"]
            if "MATR_C2" in args:
                if args["MATR_C2"] is not None:
                    motscles["MATR_C2"] = args["MATR_C2"]
            if "VECT_IMPO" in args:
                if args["VECT_IMPO"] is not None:
                    motscles["VECT_IMPO"] = args["VECT_IMPO"]
        else:
            nbsig = 6
            for index in range(nbsig):
                iks = CMP_SIG[index]
                ike = CMP_EPS[index]
                inds = 0
                inde = 0
                if ieps:
                    if EPS.get(ike) is not None:
                        inde = 1
                if isig:
                    if SIG.get(iks) is not None:
                        inds = 1
                if inde * inds != 0:
                    UTMESS("F", "COMPOR2_2", valk=iks)
                if inde == 1:
                    motscles[ike] = EPS[ike]
                elif inds == 1:
                    motscles[iks] = SIG[iks]
                else:
                    motscles[iks] = __fonczero
        #      Etat initial
        etatinit = 0
        if SIGM_INIT is not None:
            motscles["SIGM_INIT"] = SIGM_INIT
        if EPSI_INIT is not None:
            motscles["EPSI_INIT"] = EPSI_INIT
        if VARI_INIT is not None:
            motscles["VARI_INIT"] = VARI_INIT
        if NEWTON is not None:
            motscles["NEWTON"] = NEWTON
        if CONVERGENCE is not None:
            motscles["CONVERGENCE"] = CONVERGENCE
        if MASSIF is not None:
            motscles["MASSIF"] = MASSIF
        if COMPORTEMENT:
            motscles["COMPORTEMENT"] = COMPORTEMENT
        #      -- Deroulement du calcul
        motscles["INCREMENT"] = INCREMENT

        if "FORMAT_TABLE" in args:
            if args["FORMAT_TABLE"] is not None:
                motscles["FORMAT_TABLE"] = args["FORMAT_TABLE"]

        if "OPER_TANGENT" in args:
            if args["OPER_TANGENT"] is not None:
                motscles["OPER_TANGENT"] = args["OPER_TANGENT"]

        if "NB_VARI_TABLE" in args:
            if args["NB_VARI_TABLE"] is not None:
                motscles["NB_VARI_TABLE"] = args["NB_VARI_TABLE"]

        if ARCHIVAGE:
            motscles["ARCHIVAGE"] = ARCHIVAGE

            #     variables de commande
        if "AFFE_VARC" in args:
            if args["AFFE_VARC"] is not None:
                lvarc = args["AFFE_VARC"]
                nbvarc = len(lvarc)
                lmotcle = []
                for ivarc in range(nbvarc):
                    dico = {}
                    if str(lvarc[ivarc]["NOM_VARC"]) == "M_ZIRC":
                        dico["NOM_VARC"] = "ALPHPUR"
                        dico["VALE_FONC"] = lvarc[ivarc]["V1"]
                        lmotcle.append(dico)
                        dico = {}
                        dico["NOM_VARC"] = "ALPHBETA"
                        dico["VALE_FONC"] = lvarc[ivarc]["V2"]
                        lmotcle.append(dico)
                        dico = {}
                        dico["NOM_VARC"] = "BETA"
                        dico["VALE_FONC"] = lvarc[ivarc]["V3"]
                        lmotcle.append(dico)
                    elif str(lvarc[ivarc]["NOM_VARC"]) == "M_ACIER":
                        dico["NOM_VARC"] = "PFERRITE"
                        dico["VALE_FONC"] = lvarc[ivarc]["V1"]
                        lmotcle.append(dico)
                        dico = {}
                        dico["NOM_VARC"] = "PPERLITE"
                        dico["VALE_FONC"] = lvarc[ivarc]["V2"]
                        lmotcle.append(dico)
                        dico = {}
                        dico["NOM_VARC"] = "PBAINITE"
                        dico["VALE_FONC"] = lvarc[ivarc]["V3"]
                        lmotcle.append(dico)
                        dico = {}
                        dico["NOM_VARC"] = "PMARTENS"
                        dico["VALE_FONC"] = lvarc[ivarc]["V4"]
                        lmotcle.append(dico)
                        dico = {}
                        dico["NOM_VARC"] = "PAUSTENI"
                        dico["VALE_FONC"] = lvarc[ivarc]["V5"]
                        lmotcle.append(dico)
                        dico = {}
                        dico["NOM_VARC"] = "PCOLDSUM"
                        dico["VALE_FONC"] = lvarc[ivarc]["V6"]
                        lmotcle.append(dico)
                    else:
                        dico["NOM_VARC"] = lvarc[ivarc]["NOM_VARC"]
                        dico["VALE_FONC"] = lvarc[ivarc]["VALE_FONC"]

                        if (
                            str(lvarc[ivarc]["NOM_VARC"]) == "TEMP"
                            or str(lvarc[ivarc]["NOM_VARC"]) == "SECH"
                        ):
                            dico["VALE_REF"] = lvarc[ivarc]["VALE_REF"]
                        lmotcle.append(dico)
                motscles["AFFE_VARC"] = lmotcle
        if lcomp["RELATION"] == "META_LEMA_ANI":
            UTMESS("F", "COMPOR2_91", valk=lcomp["RELATION"])

        Titre = "CALC_POINT_MAT"
        if ARCHIVAGE is not None:
            #         on ne prend en compte que ARCHIVAGE / LIST_INST
            if ARCHIVAGE.get("LIST_INST") is not None:
                __REP1 = CALC_POINT_MAT(INFO=INFO, MATER=MATER, ANGLE=ANGLE, **motscles)
                lr8 = ARCHIVAGE["LIST_INST"]
                lr = lr8.getValues()
                REPONSE = CALC_TABLE(
                    TABLE=__REP1,
                    TITRE=Titre,
                    ACTION=_F(
                        OPERATION="FILTRE",
                        NOM_PARA="INST",
                        VALE=lr,
                        PRECISION=ARCHIVAGE["PRECISION"],
                    ),
                )
            else:
                REPONSE = CALC_POINT_MAT(INFO=INFO, MATER=MATER, ANGLE=ANGLE, **motscles)
        else:
            REPONSE = CALC_POINT_MAT(INFO=INFO, MATER=MATER, ANGLE=ANGLE, **motscles)

    # ===============================================================
    # cas ou on fait le calcul sur un TETRA4 A UN SEUL POINT DE GAUSS
    # ===============================================================
    elif itetra == 1:
        EPS = {}
        SIG = {}
        MODELISATION = "3D"
        if "MODELISATION" in args:
            if args["MODELISATION"] is not None:
                MODELISATION = args["MODELISATION"]

        if MODELISATION == "3D":
            nbsig = 6
            CMP_EPS = ["EPXX", "EPYY", "EPZZ", "EPXY", "EPXZ", "EPYZ"]
            CMP_SIG = ["SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ"]
        else:
            nbsig = 3
            CMP_EPS = ["EPXX", "EPYY", "EPXY"]
            CMP_SIG = ["SIXX", "SIYY", "SIXY"]

        if SIGM_IMPOSE:
            SIG = dict(SIGM_IMPOSE[0])
            for i in CMP_SIG:
                if SIG.get(i) is None:
                    SIG[i] = __fonczero
        else:
            for i in range(nbsig):
                SIG[CMP_SIG[i]] = __fonczero

        if EPSI_IMPOSE:
            EPS = dict(EPSI_IMPOSE[0])
        else:
            for i in range(nbsig):
                EPS[CMP_EPS[i]] = None

        for index in range(nbsig):
            iks = CMP_SIG[index]
            ike = CMP_EPS[index]
            if EPS.get(ike) is not None and SIG[iks] != __fonczero:
                UTMESS("F", "COMPOR2_3", valk=str(iks) + " " + str(ike))

        #     -- Definition du maillage
        if MODELISATION == "3D":
            texte_ma = """
           COOR_3D
             P0  0.0   0.0   0.0
             P1  1.0   0.0   0.0
             P2  0.0   1.0   0.0
             P3  0.0   0.0   1.0
           FINSF
           TRIA3
             F1   P0 P3 P2
             F2   P0 P1 P3
             F3   P0 P2 P1
             F4   P1 P2 P3
           FINSF
           TETRA4
             VOLUME = P0 P1 P2 P3
           FINSF
           GROUP_MA
           TOUT  VOLUME
           FINSF
           GROUP_MA
           VOLUME  VOLUME
           FINSF
           GROUP_MA
           F1  F1
           FINSF
           GROUP_MA
           F2  F2
           FINSF
           GROUP_MA
           F3  F3
           FINSF
           GROUP_MA
           F4  F4
           FINSF
           GROUP_NO
           TOUT  P1 P2 P0 P3
           FINSF
           GROUP_NO
           P0  P0
           FINSF
           GROUP_NO
           P1  P1
           FINSF
           GROUP_NO
           P2  P2
           FINSF
           GROUP_NO
           P3  P3
           FINSF
           FIN
         """

        else:
            texte_ma = """
           COOR_2D
             P0  0.0   0.0
             P1  1.0   0.0
             P2  0.0   1.0
           FINSF
           SEG2
             S1   P2 P0
             S2   P0 P1
             S3   P1 P2
           FINSF
           TRIA3
             VOLUME = P0 P1 P2
           FINSF
           GROUP_MA
           TOUT  VOLUME
           FINSF
           GROUP_MA
           VOLUME  VOLUME
           FINSF
           GROUP_MA
           S1  S1
           FINSF
           GROUP_MA
           S2  S2
           FINSF
           GROUP_MA
           S3  S3
           FINSF
           GROUP_NO
           TOUT  P1 P2 P0
           FINSF
           GROUP_NO
           P0  P0
           FINSF
           GROUP_NO
           P1  P1
           FINSF
           GROUP_NO
           P2  P2
           FINSF
           FIN
         """

        with open("simu.mail", "w") as fi_mail:
            fi_mail.write(texte_ma)

        logical_unit = LogicalUnitFile.open("simu.mail", access=FileAccess.Old)
        umail = logical_unit.unit

        __MA = LIRE_MAILLAGE(FORMAT="ASTER", UNITE=umail)
        logical_unit.release()

        if MODELISATION == "3D":
            __MO = AFFE_MODELE(
                MAILLAGE=__MA,
                AFFE=_F(
                    GROUP_MA=("VOLUME", "F1", "F2", "F3", "F4"),
                    PHENOMENE="MECANIQUE",
                    MODELISATION="3D",
                ),
                DISTRIBUTION=_F(METHODE="CENTRALISE"),
            )
            # ANGLE : rotation de ANGLE autour de Z uniquement, et seulement pour les déformations
            # imposées.
            if ANGLE != 0.0:
                __MA = MODI_MAILLAGE(
                    reuse=__MA, MAILLAGE=__MA, ROTATION=_F(POIN_1=(0.0, 0.0), ANGLE=ANGLE)
                )
                c = math.cos(ANGLE * math.pi / 180.0)
                s = math.sin(ANGLE * math.pi / 180.0)
                __C_RIGIDE = AFFE_CHAR_MECA(
                    MODELE=__MO,
                    DDL_IMPO=_F(GROUP_NO="P0", DX=0, DY=0.0, DZ=0.0),
                    LIAISON_DDL=(
                        _F(
                            GROUP_NO=("P1", "P1", "P2", "P2"),
                            DDL=("DX", "DY", "DX", "DY"),
                            COEF_MULT=(s, -c, c, s),
                            COEF_IMPO=0,
                        ),
                        _F(
                            GROUP_NO=("P1", "P3", "P3"),
                            DDL=("DZ", "DX", "DY"),
                            COEF_MULT=(-1.0, c, s),
                            COEF_IMPO=0,
                        ),
                        _F(
                            GROUP_NO=("P2", "P3", "P3"),
                            DDL=("DZ", "DX", "DY"),
                            COEF_MULT=(-1.0, -s, c),
                            COEF_IMPO=0,
                        ),
                    ),
                )
            else:
                #     -- Mouvement de corps rigide
                __C_RIGIDE = AFFE_CHAR_MECA(
                    MODELE=__MO,
                    DDL_IMPO=_F(GROUP_NO="P0", DX=0, DY=0, DZ=0),
                    LIAISON_DDL=(
                        _F(GROUP_NO=("P2", "P1"), DDL=("DX", "DY"), COEF_MULT=(1, -1), COEF_IMPO=0),
                        _F(GROUP_NO=("P3", "P1"), DDL=("DX", "DZ"), COEF_MULT=(1, -1), COEF_IMPO=0),
                        _F(GROUP_NO=("P3", "P2"), DDL=("DY", "DZ"), COEF_MULT=(1, -1), COEF_IMPO=0),
                    ),
                )
        else:
            # MODELISATION 2D
            __MO = AFFE_MODELE(
                MAILLAGE=__MA,
                AFFE=_F(
                    GROUP_MA=("VOLUME", "S1", "S2", "S3"),
                    PHENOMENE="MECANIQUE",
                    MODELISATION=MODELISATION,
                ),
                DISTRIBUTION=_F(METHODE="CENTRALISE"),
            )
            # ANGLE : rotation de ANGLE autour de Z uniquement, et seulement pour les déformations
            # imposées.
            if ANGLE != 0.0:
                __MA = MODI_MAILLAGE(
                    reuse=__MA, MAILLAGE=__MA, ROTATION=_F(POIN_1=(0.0, 0.0), ANGLE=ANGLE)
                )
                c = math.cos(ANGLE * math.pi / 180.0)
                s = math.sin(ANGLE * math.pi / 180.0)
                __C_RIGIDE = AFFE_CHAR_MECA(
                    MODELE=__MO,
                    DDL_IMPO=_F(GROUP_NO="P0", DX=0, DY=0.0),
                    LIAISON_DDL=(
                        _F(
                            GROUP_NO=("P1", "P1", "P2", "P2"),
                            DDL=("DX", "DY", "DX", "DY"),
                            COEF_MULT=(s, -c, c, s),
                            COEF_IMPO=0,
                        ),
                    ),
                )
            else:
                __C_RIGIDE = AFFE_CHAR_MECA(
                    MODELE=__MO,
                    DDL_IMPO=_F(GROUP_NO="P0", DX=0, DY=0),
                    LIAISON_DDL=(
                        _F(GROUP_NO=("P2", "P1"), DDL=("DX", "DY"), COEF_MULT=(1, -1), COEF_IMPO=0),
                    ),
                )

        #     --MASSIF : orientation du materiau (monocristal, orthotropie)
        if MASSIF:
            ANGMAS = dict(MASSIF[0])
            if ANGMAS["ANGL_REP"] is None:
                __CARA = AFFE_CARA_ELEM(
                    MODELE=__MO, MASSIF=_F(GROUP_MA="VOLUME", ANGL_EULER=ANGMAS["ANGL_EULER"])
                )
            else:
                __CARA = AFFE_CARA_ELEM(
                    MODELE=__MO, MASSIF=_F(GROUP_MA="VOLUME", ANGL_REP=ANGMAS["ANGL_REP"])
                )

        #     -- Chargement en deformation

        __E = [None] * nbsig

        __E[0] = AFFE_CHAR_MECA(
            MODELE=__MO, LIAISON_OBLIQUE=_F(GROUP_NO="P1", DX=1, ANGL_NAUT=ANGLE)
        )

        __E[1] = AFFE_CHAR_MECA(
            MODELE=__MO, LIAISON_OBLIQUE=_F(GROUP_NO="P2", DY=1, ANGL_NAUT=ANGLE)
        )

        if MODELISATION == "3D":
            __E[2] = AFFE_CHAR_MECA(
                MODELE=__MO, LIAISON_OBLIQUE=_F(GROUP_NO="P3", DZ=1, ANGL_NAUT=ANGLE)
            )

            __E[3] = AFFE_CHAR_MECA(
                MODELE=__MO, LIAISON_OBLIQUE=_F(GROUP_NO="P1", DY=1, ANGL_NAUT=ANGLE)
            )

            __E[4] = AFFE_CHAR_MECA(
                MODELE=__MO, LIAISON_OBLIQUE=_F(GROUP_NO="P1", DZ=1, ANGL_NAUT=ANGLE)
            )

            __E[5] = AFFE_CHAR_MECA(
                MODELE=__MO, LIAISON_OBLIQUE=_F(GROUP_NO="P2", DZ=1, ANGL_NAUT=ANGLE)
            )

        else:
            c = math.cos(ANGLE * math.pi / 180.0)
            s = math.sin(ANGLE * math.pi / 180.0)
            __E[2] = AFFE_CHAR_MECA(
                MODELE=__MO, LIAISON_OBLIQUE=_F(GROUP_NO="P1", DY=1, ANGL_NAUT=ANGLE)
            )

        #     -- Chargement en contrainte

        __S = [None] * nbsig

        if MODELISATION == "3D":
            r33 = 3**-0.5
            __S[0] = AFFE_CHAR_MECA(
                MODELE=__MO, FORCE_FACE=(_F(GROUP_MA="F1", FX=-1), _F(GROUP_MA="F4", FX=r33))
            )

            __S[1] = AFFE_CHAR_MECA(
                MODELE=__MO, FORCE_FACE=(_F(GROUP_MA="F2", FY=-1), _F(GROUP_MA="F4", FY=r33))
            )

            __S[2] = AFFE_CHAR_MECA(
                MODELE=__MO, FORCE_FACE=(_F(GROUP_MA="F3", FZ=-1), _F(GROUP_MA="F4", FZ=r33))
            )

            __S[3] = AFFE_CHAR_MECA(
                MODELE=__MO,
                FORCE_FACE=(
                    _F(GROUP_MA="F1", FY=-1),
                    _F(GROUP_MA="F2", FX=-1),
                    _F(GROUP_MA="F4", FX=r33, FY=r33),
                ),
            )

            __S[4] = AFFE_CHAR_MECA(
                MODELE=__MO,
                FORCE_FACE=(
                    _F(GROUP_MA="F1", FZ=-1),
                    _F(GROUP_MA="F3", FX=-1),
                    _F(GROUP_MA="F4", FX=r33, FZ=r33),
                ),
            )

            __S[5] = AFFE_CHAR_MECA(
                MODELE=__MO,
                FORCE_FACE=(
                    _F(GROUP_MA="F2", FZ=-1),
                    _F(GROUP_MA="F3", FY=-1),
                    _F(GROUP_MA="F4", FY=r33, FZ=r33),
                ),
            )

        else:
            r22 = 2**-0.5
            __S[0] = AFFE_CHAR_MECA(
                MODELE=__MO, FORCE_CONTOUR=(_F(GROUP_MA="S1", FX=-1), _F(GROUP_MA="S3", FX=r22))
            )

            __S[1] = AFFE_CHAR_MECA(
                MODELE=__MO, FORCE_CONTOUR=(_F(GROUP_MA="S2", FY=-1), _F(GROUP_MA="S3", FY=r22))
            )

            __S[2] = AFFE_CHAR_MECA(
                MODELE=__MO,
                FORCE_CONTOUR=(
                    _F(GROUP_MA="S1", FY=-1),
                    _F(GROUP_MA="S2", FX=-1),
                    _F(GROUP_MA="S3", FX=r22, FY=r22),
                ),
            )

        #     -- Construction de la charge

        l_char = [_F(CHARGE=__C_RIGIDE)]

        for i in range(nbsig):
            ike = CMP_EPS[i]
            if EPS.get(ike):
                l_char.append(_F(CHARGE=__E[i], FONC_MULT=EPS.get(ike)))

        for i in range(nbsig):
            iks = CMP_SIG[i]
            l_char.append(_F(CHARGE=__S[i], FONC_MULT=SIG.get(iks)))

        #     variables de commande
        mcvarc = []
        if "AFFE_VARC" in args:
            if args.get("AFFE_VARC") is not None:
                lvarc = args["AFFE_VARC"]
                nbvarc = len(lvarc)
                for ivarc in range(nbvarc):
                    dico = {}
                    if str(lvarc[ivarc]["NOM_VARC"]) == "SECH":
                        typech = "NOEU_TEMP_R"
                        labsc = lvarc[ivarc]["VALE_FONC"].Absc()
                        lordo = lvarc[ivarc]["VALE_FONC"].Ordo()
                        l_affe_cham = []
                        __CHV = [None] * len(labsc)
                        for it, time in enumerate(labsc):
                            __CHV[it] = (
                                CREA_CHAMP(
                                    TYPE_CHAM=typech,
                                    OPERATION="AFFE",
                                    MAILLAGE=__MA,
                                    AFFE=_F(GROUP_MA="VOLUME", NOM_CMP="TEMP", VALE=lordo[it]),
                                ),
                            )
                            dicoch = {}
                            dicoch["CHAM_GD"] = __CHV[it]
                            dicoch["INST"] = time
                            dicoch["NOM_CHAM"] = "TEMP"
                            l_affe_cham.append(dicoch)
                        __EVOV = CREA_RESU(
                            OPERATION="AFFE", TYPE_RESU="EVOL_VARC", AFFE=l_affe_cham
                        )
                    elif str(lvarc[ivarc]["NOM_VARC"]) == "M_ZIRC":
                        typech = "ELNO_NEUT_R"
                        labsc = lvarc[ivarc]["V1"].Absc()
                        lordo1 = lvarc[ivarc]["V1"].Ordo()
                        lordo2 = lvarc[ivarc]["V2"].Ordo()
                        lordo3 = lvarc[ivarc]["V3"].Ordo()
                        lordo4 = lvarc[ivarc]["V4"].Ordo()
                        lordo5 = lvarc[ivarc]["V5"].Ordo()
                        l_affe_cham = []
                        __CHV = [None] * len(labsc)
                        __CHN = [None] * len(labsc)
                        for it, time in enumerate(labsc):
                            __CHN[it] = (
                                CREA_CHAMP(
                                    TYPE_CHAM=typech,
                                    OPERATION="AFFE",
                                    PROL_ZERO="OUI",
                                    MODELE=__MO,
                                    AFFE=(
                                        _F(GROUP_MA="VOLUME", NOM_CMP="X1", VALE=lordo1[it]),
                                        _F(GROUP_MA="VOLUME", NOM_CMP="X2", VALE=lordo2[it]),
                                        _F(GROUP_MA="VOLUME", NOM_CMP="X3", VALE=lordo3[it]),
                                        _F(GROUP_MA="VOLUME", NOM_CMP="X4", VALE=lordo4[it]),
                                    ),
                                ),
                            )
                            dicoch = {}
                            __CHV[it] = CREA_CHAMP(
                                OPERATION="ASSE",
                                TYPE_CHAM="ELNO_VARI_R",
                                MODELE=__MO,
                                PROL_ZERO="OUI",
                                ASSE=(
                                    _F(
                                        CHAM_GD=__CHN[it],
                                        NOM_CMP="X1",
                                        NOM_CMP_RESU="V1",
                                        GROUP_MA="TOUT",
                                    ),
                                    _F(
                                        CHAM_GD=__CHN[it],
                                        NOM_CMP="X2",
                                        NOM_CMP_RESU="V2",
                                        GROUP_MA="TOUT",
                                    ),
                                    _F(
                                        CHAM_GD=__CHN[it],
                                        NOM_CMP="X3",
                                        NOM_CMP_RESU="V3",
                                        GROUP_MA="TOUT",
                                    ),
                                    _F(
                                        CHAM_GD=__CHN[it],
                                        NOM_CMP="X4",
                                        NOM_CMP_RESU="V4",
                                        GROUP_MA="TOUT",
                                    ),
                                    _F(
                                        CHAM_GD=__CHN[it],
                                        NOM_CMP="X5",
                                        NOM_CMP_RESU="V5",
                                        GROUP_MA="TOUT",
                                    ),
                                ),
                            )
                            dicoch["CHAM_GD"] = __CHV[it]
                            dicoch["INST"] = time
                            dicoch["NOM_CHAM"] = str("META_ELNO")
                            l_affe_cham.append(dicoch)
                        __EVOV = CREA_RESU(
                            OPERATION="AFFE", TYPE_RESU="EVOL_THER", AFFE=l_affe_cham
                        )
                        IMPR_RESU(FORMAT="RESULTAT", RESU=_F(RESULTAT=__EVOV))
                    elif str(lvarc[ivarc]["NOM_VARC"]) == "M_ACIER":
                        typech = "ELNO_NEUT_R"
                        labsc = lvarc[ivarc]["V1"].Absc()
                        lordo1 = lvarc[ivarc]["V1"].Ordo()
                        lordo2 = lvarc[ivarc]["V2"].Ordo()
                        lordo3 = lvarc[ivarc]["V3"].Ordo()
                        lordo4 = lvarc[ivarc]["V4"].Ordo()
                        lordo5 = lvarc[ivarc]["V5"].Ordo()
                        lordo6 = lvarc[ivarc]["V6"].Ordo()
                        lordo7 = lvarc[ivarc]["V7"].Ordo()
                        lordo8 = lvarc[ivarc]["V8"].Ordo()
                        lordo9 = lvarc[ivarc]["V9"].Ordo()
                        l_affe_cham = []
                        __CHV = [None] * len(labsc)
                        __CHN = [None] * len(labsc)
                        for it, time in enumerate(labsc):
                            __CHN[it] = (
                                CREA_CHAMP(
                                    TYPE_CHAM=typech,
                                    OPERATION="AFFE",
                                    PROL_ZERO="OUI",
                                    MODELE=__MO,
                                    AFFE=(
                                        _F(GROUP_MA="VOLUME", NOM_CMP="X1", VALE=lordo1[it]),
                                        _F(GROUP_MA="VOLUME", NOM_CMP="X2", VALE=lordo2[it]),
                                        _F(GROUP_MA="VOLUME", NOM_CMP="X3", VALE=lordo3[it]),
                                        _F(GROUP_MA="VOLUME", NOM_CMP="X4", VALE=lordo4[it]),
                                        _F(GROUP_MA="VOLUME", NOM_CMP="X5", VALE=lordo5[it]),
                                        _F(GROUP_MA="VOLUME", NOM_CMP="X6", VALE=lordo6[it]),
                                        _F(GROUP_MA="VOLUME", NOM_CMP="X7", VALE=lordo7[it]),
                                        _F(GROUP_MA="VOLUME", NOM_CMP="X8", VALE=lordo8[it]),
                                        _F(GROUP_MA="VOLUME", NOM_CMP="X9", VALE=lordo9[it]),
                                    ),
                                ),
                            )
                            dicoch = {}
                            __CHV[it] = CREA_CHAMP(
                                OPERATION="ASSE",
                                TYPE_CHAM="ELNO_VARI_R",
                                MODELE=__MO,
                                PROL_ZERO="OUI",
                                ASSE=(
                                    _F(
                                        CHAM_GD=__CHN[it],
                                        NOM_CMP="X1",
                                        NOM_CMP_RESU="V1",
                                        GROUP_MA="TOUT",
                                    ),
                                    _F(
                                        CHAM_GD=__CHN[it],
                                        NOM_CMP="X2",
                                        NOM_CMP_RESU="V2",
                                        GROUP_MA="TOUT",
                                    ),
                                    _F(
                                        CHAM_GD=__CHN[it],
                                        NOM_CMP="X3",
                                        NOM_CMP_RESU="V3",
                                        GROUP_MA="TOUT",
                                    ),
                                    _F(
                                        CHAM_GD=__CHN[it],
                                        NOM_CMP="X4",
                                        NOM_CMP_RESU="V4",
                                        GROUP_MA="TOUT",
                                    ),
                                    _F(
                                        CHAM_GD=__CHN[it],
                                        NOM_CMP="X5",
                                        NOM_CMP_RESU="V5",
                                        GROUP_MA="TOUT",
                                    ),
                                    _F(
                                        CHAM_GD=__CHN[it],
                                        NOM_CMP="X6",
                                        NOM_CMP_RESU="V6",
                                        GROUP_MA="TOUT",
                                    ),
                                    _F(
                                        CHAM_GD=__CHN[it],
                                        NOM_CMP="X7",
                                        NOM_CMP_RESU="V7",
                                        GROUP_MA="TOUT",
                                    ),
                                    _F(
                                        CHAM_GD=__CHN[it],
                                        NOM_CMP="X8",
                                        NOM_CMP_RESU="V8",
                                        GROUP_MA="TOUT",
                                    ),
                                    _F(
                                        CHAM_GD=__CHN[it],
                                        NOM_CMP="X9",
                                        NOM_CMP_RESU="V9",
                                        GROUP_MA="TOUT",
                                    ),
                                ),
                            )
                            dicoch["CHAM_GD"] = __CHV[it]
                            dicoch["INST"] = time
                            dicoch["NOM_CHAM"] = str("META_ELNO")
                            l_affe_cham.append(dicoch)
                        __EVOV = CREA_RESU(
                            OPERATION="AFFE", TYPE_RESU="EVOL_THER", AFFE=l_affe_cham
                        )
                        IMPR_RESU(FORMAT="RESULTAT", RESU=_F(RESULTAT=__EVOV))
                    else:
                        typech = "NOEU_" + str(lvarc[ivarc]["NOM_VARC"]) + "_R"
                        labsc = lvarc[ivarc]["VALE_FONC"].Absc()
                        lordo = lvarc[ivarc]["VALE_FONC"].Ordo()
                        l_affe_cham = []
                        __CHV = [None] * len(labsc)
                        for it, time in enumerate(labsc):
                            __CHV[it] = (
                                CREA_CHAMP(
                                    TYPE_CHAM=typech,
                                    OPERATION="AFFE",
                                    MAILLAGE=__MA,
                                    AFFE=_F(
                                        GROUP_MA="VOLUME",
                                        NOM_CMP=lvarc[ivarc]["NOM_VARC"],
                                        VALE=lordo[it],
                                    ),
                                ),
                            )
                            dicoch = {}
                            dicoch["CHAM_GD"] = __CHV[it]
                            dicoch["INST"] = time
                            dicoch["NOM_CHAM"] = str(lvarc[ivarc]["NOM_VARC"])
                            l_affe_cham.append(dicoch)
                        __EVOV = CREA_RESU(
                            OPERATION="AFFE", TYPE_RESU="EVOL_VARC", AFFE=l_affe_cham
                        )
                    dico["GROUP_MA"] = "VOLUME"
                    dico["EVOL"] = __EVOV
                    if str(lvarc[ivarc]["NOM_VARC"]) == "M_ZIRC":
                        dico["NOM_VARC"] = "M_ZIRC"
                    elif str(lvarc[ivarc]["NOM_VARC"]) == "M_ACIER":
                        dico["NOM_VARC"] = "M_ACIER"
                    else:
                        dico["NOM_VARC"] = lvarc[ivarc]["NOM_VARC"]
                        if lvarc[ivarc].get("VALE_REF") is not None:
                            dico["VALE_REF"] = lvarc[ivarc]["VALE_REF"]
                    mcvarc.append(dico)
        #      -- Materiau et modele
        if len(mcvarc) > 0:
            __CHMAT = AFFE_MATERIAU(
                MAILLAGE=__MA, AFFE=_F(GROUP_MA="VOLUME", MATER=MATER), AFFE_VARC=mcvarc
            )
        else:
            __CHMAT = AFFE_MATERIAU(MAILLAGE=__MA, AFFE=_F(GROUP_MA="VOLUME", MATER=MATER))

        #     Etat initial
        SIGINI = {}
        VARINI = {}
        LCSIG = []
        LVSIG = []
        init_dico = {}
        etatinit = 0

        #     --contraintes initiales
        if SIGM_INIT:
            etatinit = 1
            SIGINI = dict(SIGM_INIT[0])
            for i in list(SIGINI.keys()):
                if SIGINI.get(i) is not None:
                    LCSIG.append(i)
                    LVSIG.append(SIGINI[i])

            __SIG_INIT = CREA_CHAMP(
                MAILLAGE=__MA,
                OPERATION="AFFE",
                TYPE_CHAM="CART_SIEF_R",
                AFFE=_F(TOUT="OUI", NOM_CMP=LCSIG, VALE=LVSIG),
            )
            init_dico["SIGM"] = __SIG_INIT

        #     --variables internes initiales
        if VARI_INIT:
            etatinit = 1
            lnomneu = []
            lnomvar = []
            VARINI = dict(VARI_INIT[0])
            if not is_sequence(VARINI["VALE"]):
                VARINI["VALE"] = [VARINI["VALE"]]
            nbvari = len(VARINI["VALE"])
            for i in range(nbvari):
                lnomneu.append("X" + str(i + 1))
                lnomvar.append("V" + str(i + 1))

            __NEUT = CREA_CHAMP(
                OPERATION="AFFE",
                TYPE_CHAM="CART_N480_R",
                MAILLAGE=__MA,
                AFFE=_F(GROUP_MA="VOLUME", NOM_CMP=lnomneu, VALE=VARINI["VALE"]),
            )

            __VAR_INIT = CREA_CHAMP(
                MODELE=__MO,
                OPERATION="ASSE",
                TYPE_CHAM="ELGA_VARI_R",
                ASSE=_F(TOUT="OUI", CHAM_GD=__NEUT, NOM_CMP=lnomneu, NOM_CMP_RESU=lnomvar),
            )
            init_dico["VARI"] = __VAR_INIT

        # --deformations initiales
        if EPSI_INIT:
            etatinit = 1
            EPSINI = {}
            LIST_AFFE = []
            mon_dico = {}
            mon_dico["GROUP_NO"] = "P0"
            mon_dico["NOM_CMP"] = ("DX", "DY", "DZ")
            mon_dico["VALE"] = (0.0, 0.0, 0.0)
            LIST_AFFE.append(mon_dico)

            EPSINI = dict(EPSI_INIT[0])
            mon_dico = {}
            mon_dico["GROUP_NO"] = "P1"
            mon_dico["NOM_CMP"] = "DX"
            mon_dico["VALE"] = EPSINI["EPXX"]
            LIST_AFFE.append(mon_dico)
            mon_dico = {}
            mon_dico["GROUP_NO"] = "P2"
            mon_dico["NOM_CMP"] = "DY"
            mon_dico["VALE"] = EPSINI["EPYY"]
            LIST_AFFE.append(mon_dico)
            if MODELISATION == "3D":
                mon_dico = {}
                mon_dico["GROUP_NO"] = "P3"
                mon_dico["NOM_CMP"] = "DZ"
                mon_dico["VALE"] = EPSINI["EPZZ"]
                LIST_AFFE.append(mon_dico)
                mon_dico = {}
                mon_dico["GROUP_NO"] = "P1"
                mon_dico["NOM_CMP"] = "DY"
                mon_dico["VALE"] = EPSINI["EPXY"]
                LIST_AFFE.append(mon_dico)
                mon_dico = {}
                mon_dico["GROUP_NO"] = "P2"
                mon_dico["NOM_CMP"] = "DX"
                mon_dico["VALE"] = EPSINI["EPXY"]
                LIST_AFFE.append(mon_dico)
                mon_dico = {}
                mon_dico["GROUP_NO"] = "P1"
                mon_dico["NOM_CMP"] = "DZ"
                mon_dico["VALE"] = EPSINI["EPXZ"]
                LIST_AFFE.append(mon_dico)
                mon_dico = {}
                mon_dico["GROUP_NO"] = "P3"
                mon_dico["NOM_CMP"] = "DX"
                mon_dico["VALE"] = EPSINI["EPXZ"]
                LIST_AFFE.append(mon_dico)
                mon_dico = {}
                mon_dico["GROUP_NO"] = "P2"
                mon_dico["NOM_CMP"] = "DZ"
                mon_dico["VALE"] = EPSINI["EPYZ"]
                LIST_AFFE.append(mon_dico)
                mon_dico = {}
                mon_dico["GROUP_NO"] = "P3"
                mon_dico["NOM_CMP"] = "DY"
                mon_dico["VALE"] = EPSINI["EPYZ"]
                LIST_AFFE.append(mon_dico)
            else:
                mon_dico = {}
                mon_dico["GROUP_NO"] = ("P1",)
                mon_dico["NOM_CMP"] = "DY"
                mon_dico["VALE"] = EPSINI["EPXY"]
                LIST_AFFE.append(mon_dico)
                mon_dico = {}
                mon_dico["GROUP_NO"] = "P2"
                mon_dico["NOM_CMP"] = "DX"
                mon_dico["VALE"] = EPSINI["EPXY"]
                LIST_AFFE.append(mon_dico)
            __DEP_INI = CREA_CHAMP(
                MAILLAGE=__MA, OPERATION="AFFE", TYPE_CHAM="NOEU_DEPL_R", AFFE=LIST_AFFE
            )
            init_dico["DEPL"] = __DEP_INI

        #     -- Deroulement du calcul
        motscles = {}
        if COMPORTEMENT:
            motscles["COMPORTEMENT"] = COMPORTEMENT

        if lcomp["RELATION"] == "META_LEMA_ANI":
            UTMESS("A", "COMPOR2_92", valk=lcomp["RELATION"])

        motscles["CONVERGENCE"] = CONVERGENCE

        motscles["NEWTON"] = NEWTON

        if "RECH_LINEAIRE" in args:
            if args.get("RECH_LINEAIRE") is not None:
                motscles["RECH_LINEAIRE"] = args["RECH_LINEAIRE"]

        motscles["INCREMENT"] = INCREMENT

        if ARCHIVAGE:
            motscles["ARCHIVAGE"] = ARCHIVAGE

        if "SUIVI_DDL" in args:
            if args["SUIVI_DDL"] is not None:
                motscles["SUIVI_DDL"] = args["SUIVI_DDL"]

        if etatinit == 1:
            if MASSIF:
                __EVOL1 = STAT_NON_LINE(
                    INFO=INFO,
                    CARA_ELEM=__CARA,
                    MODELE=__MO,
                    CHAM_MATER=__CHMAT,
                    ETAT_INIT=init_dico,
                    EXCIT=l_char,
                    **motscles
                )
            else:
                __EVOL1 = STAT_NON_LINE(
                    INFO=INFO,
                    MODELE=__MO,
                    CHAM_MATER=__CHMAT,
                    ETAT_INIT=init_dico,
                    EXCIT=l_char,
                    **motscles
                )

        else:
            if MASSIF:
                __EVOL1 = STAT_NON_LINE(
                    INFO=INFO,
                    MODELE=__MO,
                    CARA_ELEM=__CARA,
                    CHAM_MATER=__CHMAT,
                    EXCIT=l_char,
                    **motscles
                )
            else:
                __EVOL1 = STAT_NON_LINE(
                    INFO=INFO, MODELE=__MO, CHAM_MATER=__CHMAT, EXCIT=l_char, **motscles
                )

        if lcomp["DEFORMATION"] != "PETIT":
            nomepsi = "EPSG_ELNO"
        else:
            nomepsi = "EPSI_ELNO"

        __EVOL1 = CALC_CHAMP(
            reuse=__EVOL1,
            RESULTAT=__EVOL1,
            CONTRAINTE="SIGM_ELNO",
            DEFORMATION=nomepsi,
            VARI_INTERNE="VARI_ELNO",
        )

        if MODELISATION == "3D":
            angles = (ANGLE, 0, 0)
            __EVOL = MODI_REPERE(
                RESULTAT=__EVOL1,
                MODI_CHAM=(
                    _F(NOM_CHAM="DEPL", TYPE_CHAM="VECT_3D"),
                    _F(NOM_CHAM="SIGM_ELNO", TYPE_CHAM="TENS_3D"),
                    _F(NOM_CHAM=nomepsi, TYPE_CHAM="TENS_3D"),
                ),
                REPERE="UTILISATEUR",
                AFFE=_F(ANGL_NAUT=angles, TOUT="OUI"),
            )
        else:
            angles = ANGLE
            __EVOL = MODI_REPERE(
                RESULTAT=__EVOL1,
                MODI_CHAM=(
                    _F(NOM_CHAM="DEPL", TYPE_CHAM="VECT_2D"),
                    _F(NOM_CHAM="SIGM_ELNO", TYPE_CHAM="TENS_2D"),
                    _F(NOM_CHAM=nomepsi, TYPE_CHAM="TENS_2D"),
                ),
                REPERE="UTILISATEUR",
                AFFE=_F(ANGL_NAUT=angles, TOUT="OUI"),
            )

        #     -- Recuperation des courbes

        __REP_VARI = POST_RELEVE_T(
            ACTION=(
                _F(
                    INTITULE="VARI_INT",
                    RESULTAT=__EVOL1,
                    NOM_CHAM="VARI_ELNO",
                    TOUT_CMP="OUI",
                    OPERATION="EXTRACTION",
                    GROUP_NO="P0",
                ),
            )
        )

        __REP_EPSI = POST_RELEVE_T(
            ACTION=(
                _F(
                    INTITULE="EPSILON",
                    RESULTAT=__EVOL,
                    NOM_CHAM=nomepsi,
                    TOUT_CMP="OUI",
                    OPERATION="EXTRACTION",
                    GROUP_NO="P0",
                ),
            )
        )

        __REP_SIGM = POST_RELEVE_T(
            ACTION=(
                _F(
                    INTITULE="SIGMA",
                    RESULTAT=__EVOL,
                    NOM_CHAM="SIGM_ELNO",
                    TOUT_CMP="OUI",
                    OPERATION="EXTRACTION",
                    GROUP_NO="P0",
                ),
            )
        )

        __REP_INV = POST_RELEVE_T(
            ACTION=(
                _F(
                    INTITULE="INV",
                    RESULTAT=__EVOL,
                    NOM_CHAM="SIGM_ELNO",
                    INVARIANT="OUI",
                    OPERATION="EXTRACTION",
                    GROUP_NO="P0",
                ),
            )
        )

        __REP_INV = CALC_TABLE(
            TABLE=__REP_INV,
            reuse=__REP_INV,
            ACTION=_F(OPERATION="EXTR", NOM_PARA=("INST", "TRACE", "VMIS")),
        )

        REPONSE = CALC_TABLE(
            TABLE=__REP_EPSI,
            TITRE="TABLE ",
            ACTION=(
                _F(OPERATION="COMB", TABLE=__REP_SIGM, NOM_PARA=("INST")),
                _F(OPERATION="COMB", TABLE=__REP_INV, NOM_PARA=("INST")),
                _F(OPERATION="COMB", TABLE=__REP_VARI, NOM_PARA=("INST")),
            ),
        )

    return REPONSE
