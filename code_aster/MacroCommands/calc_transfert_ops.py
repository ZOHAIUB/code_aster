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


import numpy as np

from ..Cata.Syntax import _F
from ..CodeCommands import (
    CALC_FONC_INTERP,
    CALC_FONCTION,
    CREA_TABLE,
    DEFI_FONCTION,
    DEFI_LIST_REEL,
    RECU_FONCTION,
)
from ..Messages import UTMESS


def calc_transfert_ops(
    self,
    NOM_CHAM,
    ENTREE,
    SORTIE,
    RESULTAT_X,
    RESULTAT_Y,
    RESULTAT_Z=None,
    REPERE=None,
    SIGNAL=None,
    **args
):
    """
    Macro permettant le calcul de fonctions de transfert et de signaux deconvolues
    """

    # On importe les definitions des commandes a utiliser dans la macro
    # Le nom de la variable doit etre obligatoirement le nom de la commande

    # ......................................................
    # On cherche le type de resultat des calculs dynamiques
    # ......................................................
    l_resu = [RESULTAT_X, RESULTAT_Y]  # Liste des noms des resultats
    l_type = [
        "TRAN_GENE",
        "HARM_GENE",
        "DYNA_TRANS",
        "DYNA_HARMO",
    ]  # Liste des differents types de resultats possibles
    compo = ["X", "Y"]
    entrain = ["DX", "DY"]

    if RESULTAT_Z is not None:
        l_resu.append(RESULTAT_Z)
        compo.append("Z")
        entrain.append("DZ")

    for r_type in l_type:
        if RESULTAT_X.getType() == r_type:
            if RESULTAT_Y.getType() == r_type:
                if len(l_resu) == 3 and RESULTAT_Z.getType() == r_type:
                    typ_resu = r_type
                    break
                if len(l_resu) == 2:
                    typ_resu = r_type
                    break
                else:
                    UTMESS("F", "DYNAMIQUE_34")
            else:
                UTMESS("F", "DYNAMIQUE_33")

    # ....................................................................................
    # On extrait les resultats otenus pour ensuite traiter le cas de la liste de frequence
    # ....................................................................................

    # Recu fonction pour le noeud entree
    lst_entr = []
    if ENTREE is not None:
        motsentr = ENTREE.copy()
        for rr in l_resu:
            if typ_resu == "DYNA_TRANS" or typ_resu == "DYNA_HARMO":
                motsentr["RESULTAT"] = rr
            else:
                motsentr["RESU_GENE"] = rr
            for nomcmp in compo:
                _fonc = (RECU_FONCTION(NOM_CHAM=NOM_CHAM, NOM_CMP="D" + nomcmp, **motsentr),)
                lst_entr.append(_fonc)

    # Recu fonction pour le noeud sortie
    lst_sort = []
    for sort in SORTIE:
        motssort = SORTIE.copy()
        for rr in l_resu:
            if typ_resu == "DYNA_TRANS" or typ_resu == "DYNA_HARMO":
                motssort["RESULTAT"] = rr
            else:
                motssort["RESU_GENE"] = rr
            for nomcmp in compo:
                _fonc = (RECU_FONCTION(NOM_CHAM=NOM_CHAM, NOM_CMP="D" + nomcmp, **motssort),)
                lst_sort.append(_fonc)

    ##................................................................................................
    # Extraction des parties reelles et imaginaires des fonctions obtenues via les calculs dynamiques
    # ................................................................................................
    E_Lfreq = []
    E_L_Re = []
    E_L_Im = []
    for aa in lst_entr:
        if typ_resu == "TRAN_GENE" or typ_resu == "DYNA_TRANS":
            _aaa = CALC_FONCTION(COMB=_F(FONCTION=aa, COEF=1))
            (
                Tcal,
                Ordcal,
            ) = (
                _aaa.Valeurs()
            )  # on recupere la liste d'instant pour les signaux calcules par les operateurs de dynamique
            _Tcal = DEFI_LIST_REEL(
                VALE=Tcal
            )  # On cree la liste d'instant  associee au calcul dynamique
            _Ordcal = DEFI_LIST_REEL(
                VALE=Ordcal
            )  # On cree la liste des ordonnees associee au calcul dynamique
            _f0 = DEFI_FONCTION(
                NOM_PARA="INST", VALE_PARA=_Tcal, VALE_FONC=_Ordcal
            )  # On recree la fonction temporelle pour en faire sa fft
            _fonc = CALC_FONCTION(FFT=_F(FONCTION=_f0, METHODE="COMPLET"))
            __freq, __Re, __Im = _fonc.Valeurs()
            __freq = __freq[0 : int((len(__freq)) / 2)]
            __Re = __Re[0 : int((len(__Re)) / 2)]
            __Im = __Im[0 : int((len(__Im)) / 2)]
            E_Lfreq.append(__freq)
            E_L_Re.append(__Re)
            E_L_Im.append(__Im)
        else:
            _fonc = CALC_FONCTION(COMB_C=_F(FONCTION=aa, COEF_C=(1.0 + 0j)))
            __freq, __Re, __Im = _fonc.Valeurs()
            E_Lfreq.append(__freq)
            E_L_Re.append(__Re)
            E_L_Im.append(__Im)

    S_Lfreq = []
    S_L_Re = []
    S_L_Im = []
    for bb in lst_sort:
        if typ_resu == "TRAN_GENE" or typ_resu == "DYNA_TRANS":
            _bbb = CALC_FONCTION(COMB=_F(FONCTION=bb, COEF=1))
            (
                TcalS,
                OrdcalS,
            ) = (
                _bbb.Valeurs()
            )  # on recupere la liste d'instant pour les signaux calcules par les operateurs de dynamique
            _TcalS = DEFI_LIST_REEL(
                VALE=TcalS
            )  # On cree la liste d'instant  associee au calcul dynamique
            _OrdcalS = DEFI_LIST_REEL(
                VALE=OrdcalS
            )  # On cree la liste des ordonnees associee au calcul dynamique
            _f0 = DEFI_FONCTION(
                NOM_PARA="INST", VALE_PARA=_TcalS, VALE_FONC=_OrdcalS
            )  # On recree la fonction temporelle pour en faire sa fft
            _fonc = CALC_FONCTION(FFT=_F(FONCTION=_f0, METHODE="COMPLET"))
            __freq, __Re, __Im = _fonc.Valeurs()
            __freq = __freq[0 : int((len(__freq)) / 2)]
            __Re = __Re[0 : int((len(__Re)) / 2)]
            __Im = __Im[0 : int((len(__Im)) / 2)]
            S_Lfreq.append(__freq)
            S_L_Re.append(__Re)
            S_L_Im.append(__Im)
        else:
            _fonc = CALC_FONCTION(COMB_C=_F(FONCTION=bb, COEF_C=(1.0 + 0j)))
            _freq, _Re, _Im = _fonc.Valeurs()
            S_Lfreq.append(_freq)
            S_L_Re.append(_Re)
            S_L_Im.append(_Im)

    LISTFREQ = S_Lfreq[0]
    _LIST00 = DEFI_LIST_REEL(VALE=LISTFREQ)

    # On verifie que les calculs dynamiques ont ete faits sur les memes listes
    for ii in range(1, len(E_Lfreq)):
        if len(E_Lfreq[0]) != len(E_Lfreq[ii]):
            UTMESS("F", "DYNAMIQUE_35")
            break
        else:
            _p_0 = E_Lfreq[0][1] - E_Lfreq[0][0]
            _p_ii = E_Lfreq[ii][1] - E_Lfreq[ii][0]
            if (_p_0 - _p_ii) / _p_ii >= 1.0e-6:
                UTMESS("F", "DYNAMIQUE_35")
                break

    # ..............................................
    # On determine la matrice fonction de transfert
    # ..............................................

    Ab_Lfreq = []
    Ab_L_Re = []
    Ab_L_Im = []
    LTEST = []
    if REPERE == "RELATIF":
        entrainement = (args.get("ENTRAINEMENT"),)
        for mentr in entrainement:
            s_entr = mentr
            for mm in entrain:
                if s_entr[mm].getType() == "FONCTION_C":
                    _test = s_entr[mm]
                    Test_F, Test_Re, Test_Im = _test.Valeurs()
                    LTEST.append(Test_F)
                    if (Test_F[len(Test_F) - 1] - LISTFREQ[len(LISTFREQ) - 1]) / LISTFREQ[
                        len(LISTFREQ) - 1
                    ] > 1.0e-6:
                        # On n'a pas encore traite le cas ou la frequence finale des signaux d'entrainement est plus petite que celle des calculs dynamiques
                        UTMESS("F", "DYNAMIQUE_36")
                        break
                    else:
                        _interp = CALC_FONC_INTERP(LIST_PARA=_LIST00, FONCTION=s_entr[mm])
                        A, B, C = _interp.Valeurs()
                        Ab_Lfreq.append(A)
                        Ab_L_Re.append(B)
                        Ab_L_Im.append(C)
                else:
                    _test = s_entr[mm]
                    Test_T, Test_Ord = _test.Valeurs()
                    LTEST.append(Test_T)
                    if (Tcal[len(Tcal) - 1] - Test_T[len(Test_T) - 1]) / Test_T[
                        len(Test_T) - 1
                    ] > 1e-8:
                        UTMESS("F", "DYNAMIQUE_37")
                        # On n'a pas encore traite le cas ou l'instant final des signaux d'entrainement est plus petit que celui des calculs dynamiques
                        break
                    else:
                        _interp = CALC_FONC_INTERP(LIST_PARA=_Tcal, FONCTION=s_entr[mm])
                        _fonc = CALC_FONCTION(
                            FFT=_F(FONCTION=_interp, METHODE="COMPLET")
                        )  # on fait une FFT

                        A, B, C = _fonc.Valeurs()
                        Ab_Lfreq.append(A)
                        Ab_L_Re.append(B)
                        Ab_L_Im.append(
                            C
                        )  # On recupere la demi liste et on regarde si elle est compatible
    else:
        _A = []
        _A[0 : len(LISTFREQ)] = len(LISTFREQ) * [0]
        for i in range(3):
            Ab_L_Re.append(_A)
            Ab_L_Im.append(_A)

    # On verifie que les signaux d'entrainements sont discretises de la meme maniere
    for ii in range(1, len(LTEST)):
        if len(LTEST[0]) != len(LTEST[ii]):
            UTMESS("F", "DYNAMIQUE_40")
            break
        else:
            _p_0 = LTEST[0][1] - LTEST[0][0]
            _p_ii = LTEST[ii][1] - LTEST[ii][0]
            if (_p_0 - _p_ii) / _p_ii >= 1.0e-6:
                UTMESS("F", "DYNAMIQUE_40")
                break

    # les listes valent zero tout le temps
    dim_0 = len(l_resu)
    dim = len(l_resu) ** 2
    kk = nn = ll = 0
    motsfrf = {}
    mclist = []  # mot cle facteur FONCTION

    # On cree un dictionnaire pour ranger les fonctions de transferts. Pour chaque fonction de transfert sera associe les valeurs : (freq, Re,Im)
    d_frf = {}
    for ff in range(dim):
        d_frf["LRe_%d" % ff] = []
        d_frf["LIm_%d" % ff] = []

    for i in range(len(LISTFREQ)):
        A = np.zeros((dim, dim), complex)  # tableau
        B = np.zeros((dim), complex)  # vecteur ligne
        kk = nn = ll = 0
        for ii in range(dim):
            # On remplit le vecteur B
            if ll == kk:
                B[ii] = S_L_Re[ii][i] + 1j * S_L_Im[ii][i] + Ab_L_Re[nn][i] + 1j * Ab_L_Im[nn][i]
            else:
                B[ii] = S_L_Re[ii][i] + 1j * S_L_Im[ii][i]
            # On remplit la matrice A
            for jj in range(dim_0):
                if jj == nn:
                    A[ii, jj + kk] = (
                        E_L_Re[jj + ll][i]
                        + 1j * E_L_Im[jj + ll][i]
                        + Ab_L_Re[nn][i]
                        + 1j * Ab_L_Im[nn][i]
                    )
                else:
                    A[ii, jj + kk] = E_L_Re[jj + ll][i] + 1j * E_L_Im[jj + ll][i]
            kk = kk + dim_0
            if ii == dim_0 - 1 or ii == 2 * dim_0 - 1:
                kk = 0
                ll = ll + dim_0
                nn = nn + 1
        C = np.linalg.solve(A, B)
        for nb in range(len(C)):
            _frf = C[nb]
            d_frf["LRe_%d" % nb].append(_frf.real)
            d_frf["LIm_%d" % nb].append(_frf.imag)

    mclist.append(_F(LISTE_R=LISTFREQ, PARA="FREQ"))

    if dim_0 == 2:
        Lh = ["xx", "xy", "yx", "yy"]
    else:
        Lh = ["xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"]
    for rr in range(dim):
        _LRe0 = d_frf["LRe_%d" % rr]
        _LIm0 = d_frf["LIm_%d" % rr]
        mclist.append(_F(LISTE_R=_LRe0, PARA="Re_H%s" % (Lh[rr])))
        mclist.append(_F(LISTE_R=_LIm0, PARA="Im_H%s" % (Lh[rr])))

    motsfrf["LISTE"] = mclist
    tabfrf = CREA_TABLE(**motsfrf)

    # ...................................................
    # Cas ou l'utilisateur souhaite determiner un signal
    # ...................................................

    Signfreq = []
    Sign_Re = []
    Sign_Im = []
    motssign = {}
    mcsign = []  # mot cle facteur FONCTION
    d_signal = {}
    STEST = []
    # On recupere les signaux donnes par l'utilisateur
    if SIGNAL is not None:
        type_resu = SIGNAL["TYPE_RESU"]
        l_signal = ["MESURE_X", "MESURE_Y"]

        if RESULTAT_Z is not None:
            l_signal.append("MESURE_Z")

        # On cree un dictionnaire pour ranger les signaux calcules. Pour chaque signal sera associe les valeurs : (freq, Re,Im)
        for ff in range(dim_0):
            d_signal["LRe_%d" % ff] = []
            d_signal["LIm_%d" % ff] = []
            d_signal["Lff_%d" % ff] = []

        signal = (SIGNAL,)
        for sign in signal:
            s_sign = sign
            for ss in l_signal:
                if s_sign[ss].getType() == "FONCTION_C":
                    _test = s_sign[ss]
                    Test_F, Test_Re, Test_Im = _test.Valeurs()
                    STEST.append(Test_F)
                    if Test_F[len(Test_F) - 1] < LISTFREQ[len(LISTFREQ) - 1]:
                        UTMESS("F", "DYNAMIQUE_38")
                        # On n'a pas encore traite le cas ou la frequence finale des signaux mesures est plus petite que celle des calculs dynamiques
                        break
                    else:
                        _interp = CALC_FONC_INTERP(LIST_PARA=_LIST00, FONCTION=s_sign[ss])
                        A, B, C = _interp.Valeurs()
                        Signfreq.append(A)
                        Sign_Re.append(B)
                        Sign_Im.append(C)
                else:
                    _test = s_sign[ss]
                    Test_T, Test_Ord = _test.Valeurs()
                    STEST.append(Test_T)
                    if (Tcal[len(Tcal) - 1] - Test_T[len(Test_T) - 1]) / Test_T[
                        len(Test_T) - 1
                    ] > 1.0e-8:
                        UTMESS("F", "DYNAMIQUE_39")
                        # On n'a pas encore traite le cas ou l'instant final des signaux mesures est plus petit que celui des calculs dynamiques
                        break
                    else:
                        if typ_resu == "TRAN_GENE" or typ_resu == "DYNA_TRANS":
                            _interp = CALC_FONC_INTERP(LIST_PARA=_Tcal, FONCTION=s_sign[ss])
                            _fonc = CALC_FONCTION(FFT=_F(FONCTION=_interp, METHODE="COMPLET"))
                        else:
                            if np.remainder(len(Test_T), 2) == 0:  # Si pair
                                _fonc = CALC_FONCTION(
                                    FFT=_F(FONCTION=s_sign[ss], METHODE="COMPLET")
                                )
                            else:  # Si impair
                                _Test_T = DEFI_LIST_REEL(
                                    VALE=Test_T
                                )  # On cree la liste d'instant  associee au calcul dynamique
                                _Test_Ord = DEFI_LIST_REEL(
                                    VALE=Test_Ord
                                )  # On cree la liste des ordonnees associee au calcul dynamique
                                _f0 = DEFI_FONCTION(
                                    NOM_PARA="INST", VALE_PARA=_Test_T, VALE_FONC=_Test_Ord
                                )  # On recree la fonction temporelle pour en faire sa fft
                                _fonc = CALC_FONCTION(FFT=_F(FONCTION=_f0, METHODE="COMPLET"))

                        A, B, C = _fonc.Valeurs()
                        Signfreq.append(A)
                        Sign_Re.append(B)
                        Sign_Im.append(C)

        # On verifie que les signaux mesures sont discretises de la meme maniere
        for ii in range(1, len(STEST)):
            if len(STEST[0]) != len(STEST[ii]):
                UTMESS("F", "DYNAMIQUE_41")
                break
            else:
                _p_0 = STEST[0][1] - STEST[0][0]
                _p_ii = STEST[ii][1] - STEST[ii][0]
                if (_p_0 - _p_ii) / _p_ii >= 1.0e-6:
                    UTMESS("F", "DYNAMIQUE_41")
                    break

        # On determine le signal d'entree par inversion du systeme sur chaque frequence

        for i in range(len(LISTFREQ)):
            vv = 0
            A = np.zeros((dim_0, dim_0), complex)  # tableau
            B = np.zeros((dim_0), complex)  # vecteur ligne
            for ii in range(dim_0):
                B[ii] = Sign_Re[ii][i] + 1j * Sign_Im[ii][i]
                for jj in range(dim_0):
                    A[ii, jj] = d_frf["LRe_%d" % (jj + vv)][i] + 1j * d_frf["LIm_%d" % (jj + vv)][i]
                vv = vv + dim_0

            CC = np.linalg.solve(A, B)
            for nb in range(len(CC)):
                _signal = CC[nb]
                if type_resu == "HARMONIQUE":
                    d_signal["LRe_%d" % nb].append(_signal.real)
                    d_signal["LIm_%d" % nb].append(_signal.imag)
                else:
                    d_signal["Lff_%d" % nb].append(LISTFREQ[i])
                    d_signal["Lff_%d" % nb].append(_signal.real)
                    d_signal["Lff_%d" % nb].append(_signal.imag)

        # Choix de l'utilisateur sur la table de sortie : temporel ou harmonique
        if dim_0 == 2:
            Ls = ["X", "Y"]
        else:
            Ls = ["X", "Y", "Z"]

        for rr in range(dim_0):
            if type_resu == "HARMONIQUE":
                _LRe0 = d_signal["LRe_%d" % rr]
                _LIm0 = d_signal["LIm_%d" % rr]
                mcsign.append(_F(LISTE_R=_LRe0, PARA="Re_F%s" % (Ls[rr])))
                mcsign.append(_F(LISTE_R=_LIm0, PARA="Im_F%s" % (Ls[rr])))
            else:
                _L0 = d_signal["Lff_%d" % rr]
                _AXX = DEFI_FONCTION(NOM_PARA="FREQ", VALE_C=_L0)
                _ATT = CALC_FONCTION(FFT=_F(FONCTION=_AXX, METHODE="COMPLET", SYME="NON"))
                T_en, Ord_en = _ATT.Valeurs()
                mcsign.append(_F(LISTE_R=Ord_en, PARA="F%s" % (Ls[rr])))

        if type_resu == "HARMONIQUE":
            mcsign.insert(0, _F(LISTE_R=LISTFREQ, PARA="FREQ"))
        else:
            mcsign.insert(0, _F(LISTE_R=T_en, PARA="INST"))

        motssign["LISTE"] = mcsign
        table_s = CREA_TABLE(TYPE_TABLE="TABLE", **motssign)
        self.register_result(table_s, SIGNAL[0]["TABLE_RESU"])

    return tabfrf
