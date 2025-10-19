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

#              MACRO "TEST_COMPOR"
#           ----------------------------

import numpy as NP

from ..Messages import UTMESS

from ..Cata.Syntax import _F
from ..CodeCommands import (
    CALC_FONCTION,
    CALC_TABLE,
    DEBUG,
    DEFI_CONSTANTE,
    DEFI_FONCTION,
    DEFI_LIST_INST,
    DEFI_LIST_REEL,
    FORMULE,
    IMPR_FONCTION,
    IMPR_TABLE,
    SIMU_POINT_MAT,
    TEST_TABLE,
)
from .Utils.testcomp_utils import relative_error, vect_prod_rot
from .Utils.veri_matr_tang import VERI_MATR_TANG


def rename_components_tmp(i, N_pas, label_cal, ch_param, RESU, __RS_I):
    """On renomme les composantes en fonction de  l'ordre de discrétisation"""
    N = N_pas[i]
    chN = label_cal[i] + str(N)
    for ch in ch_param:
        j = ch_param.index(ch)
        chnew = ch + chN
        # Extraction par type de variable
        if __RS_I[j] is None:
            __RS_I[j] = CALC_TABLE(
                TABLE=RESU[i],
                TITRE=" ",
                ACTION=(
                    _F(OPERATION="EXTR", NOM_PARA=("INST", ch)),
                    _F(OPERATION="RENOMME", NOM_PARA=(ch, chnew)),
                ),
            )
        else:
            __TMP_S = CALC_TABLE(
                TABLE=RESU[i],
                TITRE=" ",
                ACTION=(
                    _F(OPERATION="EXTR", NOM_PARA=("INST", ch)),
                    _F(OPERATION="RENOMME", NOM_PARA=(ch, chnew)),
                ),
            )
            __RS_I[j] = CALC_TABLE(
                reuse=__RS_I[j],
                TABLE=__RS_I[j],
                TITRE=" ",
                ACTION=(_F(OPERATION="COMB", TABLE=__TMP_S, NOM_PARA="INST"),),
            )

    return __RS_I


def TEST_ECART(self, ch_param2, label_cal, N_pas, Ncal, ch_param, __RSI, prec_ecart, prec_zero):
    # Exploitations
    CH_V1 = ["INST"]
    C_Pa = 1.0e6
    ersi = []

    for ch in ch_param2:
        # CALCUL des ecarts relatifs
        i = ch_param2.index(ch)
        chref1 = ch + label_cal[4] + str(N_pas[4])
        chref2 = ch + label_cal[Ncal - 1] + str(N_pas[Ncal - 1])
        chref = [chref1, chref2]
        preczero = prec_zero[i]
        __ersi = None

        for j in range(Ncal):
            coef = 1.0
            ch_cal = ch + label_cal[j] + str(N_pas[j])
            ch_err = "ER_" + ch_cal
            if j < 4:
                if j == 0 and i > 0 and i < 9:
                    coef = 1 / C_Pa
                iref = 0
            else:
                iref = 1
                if i == 0:
                    CH_V1.append(ch_cal)
            #               calcul de l'erreur (ecart relatif)
            valfor = "relative_error(%s,%s,%f,%e)" % (ch_cal, chref[iref], coef, preczero)
            nompar1 = "%s" % (ch_cal)
            nompar2 = "%s" % (chref[iref])
            nompar = set([nompar1, nompar2])
            __errrel = FORMULE(NOM_PARA=list(nompar), VALE=valfor, relative_error=relative_error)
            if __ersi is None:
                __ersi = CALC_TABLE(
                    TABLE=__RSI[i],
                    TITRE="__RSI" + str(j),
                    ACTION=(_F(OPERATION="OPER", NOM_PARA=ch_err, FORMULE=__errrel),),
                )
            else:
                __ersi = CALC_TABLE(
                    TABLE=__ersi,
                    reuse=__ersi,
                    TITRE="__RSI" + str(j),
                    ACTION=(_F(OPERATION="OPER", NOM_PARA=ch_err, FORMULE=__errrel),),
                )
        ersi.append(__ersi)
    return ersi


#


def CHAR3D(self, POISSON, YOUNG, _tempsar, INFO):
    # definition du trajet de chargement 3D

    #
    # fonctions chargement
    calibrage = 3.5
    Ctau = (2 * calibrage) * ((1 + POISSON) / (2 * YOUNG))
    Ctrac = calibrage * (1 / YOUNG)
    __Ytrac = DEFI_LIST_REEL(VALE=(0.0, 150.0, 150.0, -50, 0.0, 50.0, -150.0, -150.0, 0.0))
    __Trace = DEFI_FONCTION(
        NOM_PARA="INST",
        VALE_PARA=_tempsar,
        VALE_FONC=__Ytrac,
        PROL_DROITE="EXCLU",
        PROL_GAUCHE="EXCLU",
    )

    __YDevi1 = DEFI_LIST_REEL(VALE=(0.0, 75.0, 150.0, 150, 0.0, -150.0, -150.0, -75.0, 0.0))
    __Devi1 = DEFI_FONCTION(
        NOM_PARA="INST",
        VALE_PARA=_tempsar,
        VALE_FONC=__YDevi1,
        PROL_DROITE="EXCLU",
        PROL_GAUCHE="EXCLU",
    )

    __YDevi2 = DEFI_LIST_REEL(VALE=(0.0, 75.0, -50.0, 100, 0.0, -100.0, 50.0, -75.0, 0.0))
    __Devi2 = DEFI_FONCTION(
        NOM_PARA="INST",
        VALE_PARA=_tempsar,
        VALE_FONC=__YDevi2,
        PROL_DROITE="EXCLU",
        PROL_GAUCHE="EXCLU",
    )

    __Ytauxy = DEFI_LIST_REEL(VALE=(0.0, 200.0, 100.0, 300, 0.0, -300.0, -100.0, -200.0, 0.0))
    __TAUxy = DEFI_FONCTION(
        NOM_PARA="INST",
        VALE_PARA=_tempsar,
        VALE_FONC=__Ytauxy,
        PROL_DROITE="EXCLU",
        PROL_GAUCHE="EXCLU",
    )

    __Ytauxz = DEFI_LIST_REEL(VALE=(0.0, -100.0, 100.0, 200, 0.0, -200.0, -100.0, 100.0, 0.0))
    __TAUxz = DEFI_FONCTION(
        NOM_PARA="INST",
        VALE_PARA=_tempsar,
        VALE_FONC=__Ytauxz,
        PROL_DROITE="EXCLU",
        PROL_GAUCHE="EXCLU",
    )

    __Ytauyz = DEFI_LIST_REEL(VALE=(0.0, 0.0, 200.0, -100, 0.0, 100.0, -200.0, 0.0, 0.0))
    __TAUyz = DEFI_FONCTION(
        NOM_PARA="INST",
        VALE_PARA=_tempsar,
        VALE_FONC=__Ytauyz,
        PROL_DROITE="EXCLU",
        PROL_GAUCHE="EXCLU",
    )

    __eps_xy = CALC_FONCTION(COMB=_F(FONCTION=__TAUxy, COEF=Ctau))

    __eps_xz = CALC_FONCTION(COMB=_F(FONCTION=__TAUxz, COEF=Ctau))

    __eps_yz = CALC_FONCTION(COMB=_F(FONCTION=__TAUyz, COEF=Ctau))

    __eps_xx = CALC_FONCTION(
        COMB=(_F(FONCTION=__Trace, COEF=Ctrac), _F(FONCTION=__Devi1, COEF=Ctrac))
    )

    __eps_yy = CALC_FONCTION(
        COMB=(
            _F(FONCTION=__Trace, COEF=Ctrac),
            _F(FONCTION=__Devi1, COEF=-(Ctrac)),
            _F(FONCTION=__Devi2, COEF=Ctrac),
        )
    )

    __eps_zz = CALC_FONCTION(
        COMB=(_F(FONCTION=__Trace, COEF=Ctrac), _F(FONCTION=__Devi2, COEF=-(Ctrac)))
    )
    eps_imp = [__eps_xx, __eps_yy, __eps_zz, __eps_xy, __eps_xz, __eps_yz]
    # rotation tenseur des def
    # angle de precession, nutation, rotation propre
    psi, teta, phi = 0.9, 0.7, 0.4
    cps, cte, cph = NP.cos(psi), NP.cos(teta), NP.cos(phi)
    sps, ste, sph = NP.sin(psi), NP.sin(teta), NP.sin(phi)
    # matrice de passage
    p11, p21, p31 = cph * cps - sph * cte * sps, cph * sps + sph * cte * cps, sph * ste
    p12, p22, p32 = -sph * cps - cph * cte * sps, -sph * sps + cph * cte * cps, cph * ste
    p13, p23, p33 = ste * sps, -ste * cps, cte
    V1 = [p11, p21, p31]
    V2 = [p12, p22, p32]
    V3 = [p13, p23, p33]
    # eps apres rotation
    VI = [[V1, V1], [V2, V2], [V3, V3], [V1, V2], [V1, V3], [V2, V3]]
    for vect_i in VI:
        V_COEF = vect_prod_rot(vect_i[0], vect_i[1])
        __epsrot = CALC_FONCTION(
            COMB=(
                _F(FONCTION=__eps_xx, COEF=V_COEF[0]),
                _F(FONCTION=__eps_yy, COEF=V_COEF[1]),
                _F(FONCTION=__eps_zz, COEF=V_COEF[2]),
                _F(FONCTION=__eps_xy, COEF=V_COEF[3]),
                _F(FONCTION=__eps_xz, COEF=V_COEF[4]),
                _F(FONCTION=__eps_yz, COEF=V_COEF[5]),
            )
        )

    # trace chargement
    if INFO == 2:
        IMPR_FONCTION(  # FORMAT='XMGRACE',PILOTE='INTERACTIF',
            COURBE=(
                _F(FONCTION=eps_imp[0]),
                _F(FONCTION=eps_imp[1]),
                _F(FONCTION=eps_imp[2]),
                _F(FONCTION=eps_imp[3]),
                _F(FONCTION=eps_imp[4]),
                _F(FONCTION=eps_imp[5]),
            ),
            UNITE=29,
        )
    return eps_imp


#


def CHAR2D(self, POISSON, YOUNG, _tempsar, INFO):
    # definition du trajet de chargement 2D

    # fonctions chargement
    calibrage = 4.5
    Ctau = (2 * calibrage) * ((1 + POISSON) / (2 * YOUNG))
    Ctrac = calibrage * (1 / YOUNG)
    __YTrace2d = DEFI_LIST_REEL(VALE=(0.0, 150.0, 200.0, 300.0, 0.0, -300.0, -200.0, -150.0, 0.0))
    __Trace2d = DEFI_FONCTION(
        NOM_PARA="INST",
        VALE_PARA=_tempsar,
        VALE_FONC=__YTrace2d,
        PROL_DROITE="EXCLU",
        PROL_GAUCHE="EXCLU",
    )

    __Y_Devi2 = DEFI_LIST_REEL(VALE=(0.0, 0.0, 100.0, 0.0, 0.0, 0.0, -100.0, 0.0, 0.0))
    __Devi1_2d = DEFI_FONCTION(
        NOM_PARA="INST",
        VALE_PARA=_tempsar,
        VALE_FONC=__Y_Devi2,
        PROL_DROITE="EXCLU",
        PROL_GAUCHE="EXCLU",
    )

    __Y_tauxy2 = DEFI_LIST_REEL(VALE=(0.0, 100.0, 40.0, 0.0, 0.0, 0.0, -40.0, -100.0, 0.0))
    __TAU_xy2 = DEFI_FONCTION(
        NOM_PARA="INST",
        VALE_PARA=_tempsar,
        VALE_FONC=__Y_tauxy2,
        PROL_DROITE="EXCLU",
        PROL_GAUCHE="EXCLU",
    )

    __eps_xy2 = CALC_FONCTION(COMB=_F(FONCTION=__TAU_xy2, COEF=Ctau))

    __eps_xx2 = CALC_FONCTION(
        COMB=(_F(FONCTION=__Trace2d, COEF=Ctrac), _F(FONCTION=__Devi1_2d, COEF=Ctrac))
    )

    __eps_yy2 = CALC_FONCTION(
        COMB=(_F(FONCTION=__Trace2d, COEF=Ctrac), _F(FONCTION=__Devi1_2d, COEF=-(Ctrac)))
    )
    # rotation tenseur des def
    c, s = NP.cos(0.9), NP.sin(0.9)
    c2, s2, cs = c * c, s * s, c * s
    __eps_rr2 = CALC_FONCTION(
        COMB=(
            _F(FONCTION=__eps_xx2, COEF=c2),
            _F(FONCTION=__eps_yy2, COEF=s2),
            _F(FONCTION=__eps_xy2, COEF=2 * cs),
        )
    )

    __eps_tt2 = CALC_FONCTION(
        COMB=(
            _F(FONCTION=__eps_xx2, COEF=s2),
            _F(FONCTION=__eps_yy2, COEF=c2),
            _F(FONCTION=__eps_xy2, COEF=-2 * cs),
        )
    )

    __eps_rt2 = CALC_FONCTION(
        COMB=(
            _F(FONCTION=__eps_xx2, COEF=-cs),
            _F(FONCTION=__eps_yy2, COEF=cs),
            _F(FONCTION=__eps_xy2, COEF=c2 - s2),
        )
    )
    eps_imp = [__eps_xx2, __eps_yy2, __eps_xy2, __eps_rr2, __eps_tt2, __eps_rt2]

    # Trace2d chargement
    if INFO == 2:
        IMPR_FONCTION(  # FORMAT='XMGRACE',PILOTE='INTERACTIF',
            COURBE=(
                _F(FONCTION=eps_imp[0]),
                _F(FONCTION=eps_imp[1]),
                _F(FONCTION=eps_imp[2]),
                _F(FONCTION=eps_imp[3]),
                _F(FONCTION=eps_imp[4]),
                _F(FONCTION=eps_imp[5]),
            ),
            UNITE=29,
        )
    return eps_imp


#


def test_compor_ops(self, **args):
    # seule l'option "THER", c'est à dire le test thermomecanique est programmé à ce jour
    # ajouter l'option MECA (tests comp001,002), l'option HYDR, etc..
    OPTION = args.get("OPTION")
    NEWTON = args.get("NEWTON")
    CONVERGENCE = args.get("CONVERGENCE")
    COMPORTEMENT = args.get("COMPORTEMENT")
    LIST_MATER = args.get("LIST_MATER")
    VARI_TEST = args.get("VARI_TEST")
    INFO = args.get("INFO")

    U = None

    motscles = {}
    if COMPORTEMENT:
        motscles["COMPORTEMENT"] = COMPORTEMENT

    motscles["CONVERGENCE"] = CONVERGENCE
    motscles["NEWTON"] = NEWTON

    if OPTION == "THER":
        epsi = 1.0e-10
        MATER = args.get("MATER")
        ALPHA = args.get("ALPHA")
        YOUNG = args.get("YOUNG")
        TEMP_INIT = args.get("TEMP_INIT")
        TEMP_FIN = args.get("TEMP_FIN")
        NB_VARI = args.get("NB_VARI")

        if args.get("INST_FIN") is not None:
            INST_FIN = args.get("INST_FIN")

        NCAL = len(LIST_MATER)

        __LINST0 = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=(_F(JUSQU_A=INST_FIN, NOMBRE=NCAL),))

        __LINST = DEFI_LIST_INST(
            METHODE="MANUEL",
            DEFI_LIST=_F(
                LIST_INST=__LINST0,
                # PAS_MINI=1.0E-12
            ),
            ECHEC=_F(EVENEMENT="ERREUR", ACTION="DECOUPE", SUBD_METHODE="MANUEL", SUBD_PAS=10),
            # ADAPTATION=_F(EVENEMENT='SEUIL'),
        )
        __TIMP = DEFI_FONCTION(
            NOM_PARA="INST", NOM_RESU="TEMP", VALE=(0.0, TEMP_INIT, INST_FIN, TEMP_FIN)
        )

        __zero = DEFI_CONSTANTE(VALE=0.0)

        if args.get("SUPPORT") is not None:
            motscles["SUPPORT"] = args.get("SUPPORT")

        __U = SIMU_POINT_MAT(
            MATER=MATER,
            INFO=INFO,
            AFFE_VARC=(_F(NOM_VARC="TEMP", VALE_FONC=__TIMP, VALE_REF=TEMP_INIT),),
            INCREMENT=_F(LIST_INST=__LINST),
            EPSI_IMPOSE=_F(EPXX=__zero),
            **motscles
        )
        # On ne garde pas les valeurs initiales (NUME_ORDRE = 0 exclu)

        __U = CALC_TABLE(
            reuse=__U,
            TABLE=__U,
            ACTION=(_F(OPERATION="FILTRE", NOM_PARA="INST", CRIT_COMP="GT", VALE=0.0),),
        )

        SXM = 0.0
        EXM = 0.0
        time = 0.0

        __RES = [None] * (NCAL)
        if NB_VARI > 0:
            Vim = NP.zeros(NB_VARI)

        for i in range(NCAL):
            timem = time

            time = timem + INST_FIN / NCAL

            Ti = TEMP_INIT + time / INST_FIN * (TEMP_FIN - TEMP_INIT)

            Tm = TEMP_INIT + timem / INST_FIN * (TEMP_FIN - TEMP_INIT)

            # deformation mecanique imposee correspondant a la deformation
            # thermique du premier calcul

            __epsimp = DEFI_CONSTANTE(VALE=-ALPHA(Ti) * (Ti - TEMP_INIT))

            # variation des coef du COMPORtement avec la temperature
            # correction eventuelle des valeurs initiales du temps ti

            if i > 0:
                SXM = SXM * (YOUNG(Ti) / YOUNG(Tm))
                # cas particuliers
                if COMPORTEMENT[0]["RELATION"] == "VMIS_CINE_LINE":
                    D_SIGM_EPSI = args["D_SIGM_EPSI"]
                    Vim[0:5] = Vim[0:5] * D_SIGM_EPSI(Ti) / D_SIGM_EPSI(Tm)

                if COMPORTEMENT[0]["RELATION"] == "VMIS_ECMI_LINE":
                    C_PRAG = args["C_PRAG"]
                    Vim[2:7] = Vim[2:7] * C_PRAG(Ti) / C_PRAG(Tm)

                if COMPORTEMENT[0]["RELATION"] == "VMIS_ECMI_TRAC":
                    C_PRAG = args["C_PRAG"]
                    Vim[2:7] = Vim[2:7] * C_PRAG(Ti) / C_PRAG(Tm)

            __list0 = DEFI_LIST_REEL(DEBUT=timem, INTERVALLE=(_F(JUSQU_A=time, NOMBRE=1),))

            __list = DEFI_LIST_INST(
                METHODE="MANUEL",
                DEFI_LIST=_F(
                    LIST_INST=__list0,
                    # PAS_MINI=1.0E-12
                ),
                ECHEC=_F(EVENEMENT="ERREUR", ACTION="DECOUPE", SUBD_METHODE="MANUEL"),
                # ADAPTATION=_F(EVENEMENT='SEUIL'),
            )
            if NB_VARI > 0:
                __RES[i] = SIMU_POINT_MAT(
                    INFO=INFO,
                    MATER=LIST_MATER[i],
                    INCREMENT=_F(LIST_INST=__list),
                    EPSI_IMPOSE=_F(EPXX=__epsimp),
                    SIGM_INIT=_F(SIXX=SXM),
                    VARI_INIT=_F(VALE=[Vim[j] for j in range(NB_VARI)]),
                    EPSI_INIT=_F(EPXX=EXM, EPYY=0.0, EPZZ=0.0, EPXY=0.0, EPXZ=0.0, EPYZ=0.0),
                    **motscles
                )

            else:
                __RES[i] = SIMU_POINT_MAT(
                    INFO=INFO,
                    MATER=LIST_MATER[i],
                    INCREMENT=_F(LIST_INST=__list),
                    EPSI_IMPOSE=_F(EPXX=__epsimp),
                    SIGM_INIT=_F(SIXX=SXM),
                    EPSI_INIT=_F(EPXX=EXM, EPYY=0.0, EPZZ=0.0, EPXY=0.0, EPXZ=0.0, EPYZ=0.0),
                    **motscles
                )

            # On ne teste pas les valeurs initiales (NUME_ORDRE = 0 exclu)
            __RES[i] = CALC_TABLE(
                reuse=__RES[i],
                TABLE=__RES[i],
                ACTION=(_F(OPERATION="FILTRE", NOM_PARA="INST", CRIT_COMP="EQ", VALE=time),),
            )

            # recuperation des valeurs initiales du futur pas de temps dans la
            # table resultat

            EXM = __RES[i]["EPXX", 1]

            SXM = __RES[i]["SIXX", 1]

            if NB_VARI > 0:
                for j in range(NB_VARI):
                    Vim[j] = __RES[i]["V" + str(j + 1), 1]

            # On ne peut pas faire de non régression puisqu'on ne connait pas ici
            # la valeur obtenue sur la machine de référence
            TEST_TABLE(
                TABLE=__RES[i],
                NOM_PARA="VMIS",
                VALE_CALC=0.0,
                VALE_REFE=__U["VMIS", i + 1],
                TOLE_MACHINE=1.0e-3,
                FILTRE=_F(NOM_PARA="INST", VALE=time),
                REFERENCE="AUTRE_ASTER",
            )

            TEST_TABLE(
                TABLE=__RES[i],
                NOM_PARA="TRACE",
                VALE_CALC=0.0,
                VALE_REFE=__U["TRACE", i + 1],
                TOLE_MACHINE=1.0e-3,
                FILTRE=_F(NOM_PARA="INST", VALE=time),
                REFERENCE="AUTRE_ASTER",
            )
            if NB_VARI > 0:
                if VARI_TEST is not None:
                    for j in range(len(VARI_TEST)):
                        nomvari = VARI_TEST[j]
                        if abs(__U[nomvari, i + 1]) > epsi:
                            TEST_TABLE(
                                TABLE=__RES[i],
                                NOM_PARA=nomvari,
                                VALE_CALC=0.0,
                                VALE_REFE=__U[nomvari, i + 1],
                                TOLE_MACHINE=1.0e-3,
                                FILTRE=_F(NOM_PARA="INST", VALE=time),
                                REFERENCE="AUTRE_ASTER",
                            )
                else:
                    for j in range(NB_VARI):
                        nomvari = "V" + str(j + 1)
                        if abs(__U[nomvari, i + 1]) > epsi:
                            TEST_TABLE(
                                TABLE=__RES[i],
                                NOM_PARA=nomvari,
                                VALE_CALC=0.0,
                                VALE_REFE=__U[nomvari, i + 1],
                                TOLE_MACHINE=1.0e-3,
                                FILTRE=_F(NOM_PARA="INST", VALE=time),
                                REFERENCE="AUTRE_ASTER",
                            )

    elif OPTION == "MECA":
        TEST_TANGENTE = args.get("TEST_TANGENTE")
        LIST_NPAS = args.get("LIST_NPAS")
        YOUNG = args.get("YOUNG")
        POISSON = args.get("POISSON")

        # Discretisation du calcul
        if args.get("LIST_NPAS") is not None:
            LIST_NPAS = args.get("LIST_NPAS")
        else:
            LIST_NPAS = 4 * [1] + [1, 5, 25]
        Ncal = len(LIST_NPAS)

        # les differents calculs effectues et les precisions sur chaque
        # TEST_RESU
        label_cal = ["_Pa_", "_Th_", "_sym_", "_rot_"] + (Ncal - 4) * ["_N"]
        if args.get("LIST_TOLE") is not None:
            LIST_TOLE = args.get("LIST_TOLE")
        else:
            LIST_TOLE = 4 * [1.0e-10] + [1.0e-1] + (Ncal - 5) * [1.0e-2] + [1.0e-8]

        if args.get("PREC_ZERO") is not None:
            PREC_ZERO = args["PREC_ZERO"]
            if len(PREC_ZERO) != len(VARI_TEST):
                UTMESS("F", "TESTCOMPOR_1")
        else:
            PREC_ZERO = len(VARI_TEST) * [1.0e-10]

        prec_tgt = LIST_TOLE[-1]

        # parametres vitesse de sollicitation
        t_0 = 1.0
        if COMPORTEMENT:
            if COMPORTEMENT[0]["RELATION"][0:4] == "VISC":
                vitesse = 1.0e-5
                t_0 = 5.0e-2 / (8.0 * vitesse)
        # liste d'archivage
        __tempsar = DEFI_LIST_REEL(VALE=[t_0 * i for i in range(9)])

        if args.get("MODELISATION") is not None:
            modelisation = args.get("MODELISATION")
        else:
            modelisation = "3D"
        if modelisation == "3D":
            #
            #  TEST 3D
            #

            eps_imp = CHAR3D(self, POISSON, YOUNG, __tempsar, INFO)

            ch_param2 = list(VARI_TEST)
            ch_param = ch_param2 + ["SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ"]

        elif modelisation == "C_PLAN":
            #
            #  TEST 2D C_PLAN
            #

            eps_imp = CHAR2D(self, POISSON, YOUNG, __tempsar, INFO)

            # les quantites (invariants...) sur lequels portent les calculs
            # d'erreur et les test_resu
            ch_param2 = list(VARI_TEST)
            ch_param = ch_param2 + ["SIXX", "SIYY", "SIZZ", "SIXY"]

        #
        # Discretisation et calcul
        #
        __RES = [None] * Ncal
        __RSI = [None] * len(ch_param)

        # pointeur materiau
        P_imat = [0] + [1] + (Ncal - 2) * [1]
        # pointeur deformation
        nomdef, symdef, rotdef = [0, 1, 2], [1, 0, 2], [3, 4, 5]
        P_idef = [nomdef, nomdef, symdef, rotdef]
        for i in range(Ncal - 4):
            P_idef.append(nomdef)

        ldicoeps = []
        dicoeps = {}
        dicoeps["EPXX"] = eps_imp[0]
        dicoeps["EPYY"] = eps_imp[1]
        dicoeps["EPZZ"] = eps_imp[2]
        dicoeps["EPXY"] = eps_imp[3]
        if modelisation == "3D":
            dicoeps["EPXZ"] = eps_imp[4]
            dicoeps["EPYZ"] = eps_imp[5]
        ldicoeps.append(dicoeps)
        motscles["EPSI_IMPOSE"] = ldicoeps

        if args.get("SUPPORT") is not None:
            if modelisation == "C_PLAN":
                motscles["MODELISATION"] = "C_PLAN"
                motscles["SUPPORT"] = "ELEMENT"
            else:
                motscles["SUPPORT"] = args.get("SUPPORT")

        # Boucle sur l'ensemble des calculs
        for i in range(Ncal):
            N = LIST_NPAS[i]
            imat = P_imat[i]
            __temps = DEFI_LIST_REEL(
                DEBUT=0.0,
                INTERVALLE=(
                    _F(JUSQU_A=t_0, NOMBRE=N),
                    _F(JUSQU_A=2.0 * t_0, NOMBRE=N),
                    _F(JUSQU_A=3.0 * t_0, NOMBRE=N),
                    _F(JUSQU_A=4.0 * t_0, NOMBRE=N),
                    _F(JUSQU_A=5.0 * t_0, NOMBRE=N),
                    _F(JUSQU_A=6.0 * t_0, NOMBRE=N),
                    _F(JUSQU_A=7.0 * t_0, NOMBRE=N),
                    _F(JUSQU_A=8.0 * t_0, NOMBRE=N),
                ),
            )

            __DEFLIST = DEFI_LIST_INST(
                DEFI_LIST=_F(LIST_INST=__temps),
                ECHEC=_F(
                    EVENEMENT="ERREUR",
                    ACTION="DECOUPE",
                    SUBD_METHODE="MANUEL",
                    SUBD_PAS=10,
                    SUBD_NIVEAU=10,
                ),
            )
            #       Resout le pb a deformation imposee
            __RES[i] = SIMU_POINT_MAT(
                INFO=INFO,
                MATER=LIST_MATER[imat],
                ARCHIVAGE=_F(LIST_INST=__tempsar),
                INCREMENT=_F(LIST_INST=__DEFLIST),
                **motscles
            )

            # On renomme les composantes en fonction de  l'ordre de discretisation
            rename_components_tmp(i, LIST_NPAS, label_cal, ch_param, __RES, __RSI)

        # TEST_RESU sur les erreurs relatives
        # les quantites (invariants...) sur lequels portent les calculs
        # d'erreur et les test_resu
        __RSI = TEST_ECART(
            self, ch_param2, label_cal, LIST_NPAS, Ncal, ch_param, __RSI, LIST_TOLE, PREC_ZERO
        )

        # On ne peut pas faire de non régression puisqu'on ne connait pas ici
        # la valeur obtenue sur la machine de référence
        for ch in ch_param2:
            i = ch_param2.index(ch)
            if INFO == 2:
                IMPR_TABLE(TABLE=__RSI[i])
            for j in range(Ncal):
                ch_cal = ch + label_cal[j] + str(LIST_NPAS[j])
                ch_err = "ER_" + ch_cal
                TEST_TABLE(
                    TABLE=__RSI[i],
                    NOM_PARA=ch_err,
                    TYPE_TEST="MAX",
                    VALE_CALC=0.0,
                    VALE_REFE=0.0,
                    CRITERE="ABSOLU",
                    TOLE_MACHINE=LIST_TOLE[j],
                    PRECISION=LIST_TOLE[j],
                    REFERENCE="ANALYTIQUE",
                )
        #
        # Test de la matrice tangente sur le calcul le plus fin
        #

        if TEST_TANGENTE == "OUI":
            N = LIST_NPAS[Ncal - 1]
            __Linst = DEFI_LIST_REEL(
                DEBUT=0.0,
                INTERVALLE=(
                    _F(JUSQU_A=t_0, NOMBRE=N),
                    _F(JUSQU_A=2.0 * t_0, NOMBRE=N),
                    _F(JUSQU_A=3.0 * t_0, NOMBRE=N),
                    _F(JUSQU_A=4.0 * t_0, NOMBRE=N),
                    _F(JUSQU_A=5.0 * t_0, NOMBRE=N),
                    _F(JUSQU_A=6.0 * t_0, NOMBRE=N),
                    _F(JUSQU_A=7.0 * t_0, NOMBRE=N),
                    _F(JUSQU_A=8.0 * t_0, NOMBRE=N),
                ),
            )

            if COMPORTEMENT:
                motscles["COMPORTEMENT"][0]["TYPE_MATR_TANG"] = "VERIFICATION"
                if args.get("VERI_MATR_OPTION") is not None:
                    motscles["COMPORTEMENT"][0]["VALE_PERT_RELA"] = args.get("VERI_MATR_OPTION")[0][
                        "VALE_PERT_RELA"
                    ]
                __DEFLIS2 = DEFI_LIST_INST(
                    DEFI_LIST=_F(LIST_INST=__Linst),
                    ECHEC=_F(
                        EVENEMENT="ERREUR",
                        ACTION="DECOUPE",
                        SUBD_METHODE="MANUEL",
                        SUBD_PAS=10,
                        SUBD_NIVEAU=10,
                    ),
                )
                #       Resout le pb a deformation imposee
                DEBUG(SDVERI="NON")
                U = SIMU_POINT_MAT(
                    INFO=INFO,
                    MATER=LIST_MATER[imat],
                    ARCHIVAGE=_F(LIST_INST=__tempsar),
                    INCREMENT=_F(LIST_INST=__DEFLIS2),
                    **motscles
                )
                DEBUG(SDVERI="OUI")
                motscles = {}
                if args.get("VERI_MATR_OPTION") is not None:
                    motscles["PRECISION"] = args.get("VERI_MATR_OPTION")[0]["PRECISION"]
                    motscles["PREC_ZERO"] = args.get("VERI_MATR_OPTION")[0]["PREC_ZERO"]

                __DIFFMAT = VERI_MATR_TANG(**motscles)

                # On ne peut pas faire de non régression puisqu'on ne connait pas ici
                # la valeur obtenue sur la machine de référence
                TEST_TABLE(
                    TABLE=__DIFFMAT,
                    NOM_PARA="MAT_DIFF",
                    TYPE_TEST="MAX",
                    VALE_CALC=0.0,
                    VALE_REFE=0.0,
                    CRITERE="ABSOLU",
                    TOLE_MACHINE=prec_tgt,
                    PRECISION=prec_tgt,
                    REFERENCE="ANALYTIQUE",
                )
                if INFO == 2:
                    IMPR_TABLE(TABLE=__DIFFMAT)
    if U is None:
        return
    else:
        return U
