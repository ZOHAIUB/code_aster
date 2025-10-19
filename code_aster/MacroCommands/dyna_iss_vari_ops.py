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

# person_in_charge: irmela.zentner at edf.fr

"""Commande DYNA_ISS_VARI"""

import sys
import traceback
from math import ceil, floor, pi
import numpy as NP

import aster
from ..Cata.Syntax import _F
from ..CodeCommands import (
    CALC_FONCTION,
    COMB_MATR_ASSE,
    DEFI_FONCTION,
    DEFI_INTE_SPEC,
    DYNA_VIBRA,
    LIRE_FORC_MISS,
    LIRE_IMPE_MISS,
    REST_SPEC_TEMP,
)
from ..Messages import UTMESS
from .Miss.calc_miss_vari import compute_corr_vari, compute_mecmode, compute_pod
from .Utils.signal_correlation_utils import calc_coherency_matrix, calc_dist2, get_group_nom_coord


def dyna_iss_vari_ops(self, **kwargs):
    """Corps de la macro DYNA_ISS_VARI"""
    # conteneur des paramètres du calcul
    params = DynaISSParameters(**kwargs)
    # création de l'objet generator
    generator = Generator.factory(self, params)
    try:
        return generator.run()
    except Exception as err:
        trace = "".join(traceback.format_tb(sys.exc_info()[2]))
        UTMESS("F", "SUPERVIS2_3", valk=("DYNA_ISS_VARI", trace, str(err)))


class DynaISSParameters:
    def __init__(self, **kwargs):
        """Enregistrement des valeurs des mots-clés dans un dictionnaire.
        - On cree deux dictionnaires de parametres:
                      signal_process_keys et simulation_keys
        """
        others = list(kwargs.keys())
        cohekeys = kwargs.get("MATR_COHE")[0]
        genekeys = kwargs.get("MATR_GENE")[0]
        interfkeys = kwargs.get("INTERF")[0]
        self.cohe_keys = cohekeys.cree_dict_valeurs(cohekeys.mc_liste)
        self.mat_gene_keys = genekeys.cree_dict_valeurs(genekeys.mc_liste)
        self.interf_keys = interfkeys.cree_dict_valeurs(interfkeys.mc_liste)
        self.cas = "SPEC"
        if kwargs.get("EXCIT_SOL"):
            self.cas = "TRANS"
            excit_sol = kwargs["EXCIT_SOL"][0]
            self.excit_sol_keys = excit_sol.cree_dict_valeurs(excit_sol.mc_liste)
            others.remove("EXCIT_SOL")
        others.remove("MATR_GENE")
        others.remove("MATR_COHE")
        others.remove("INTERF")
        self.other_keys = {}
        for key in others:
            self.other_keys[key] = kwargs.get(key)
        self.other_keys["CASE"] = self.cas
        #  TEST: ERREUR SI ISSF ET INTERF AUTRE QUE QUELCONQUE!
        if self.other_keys["ISSF"] == "OUI":
            if self.interf_keys["MODE_INTERF"] != "QUELCONQUE":
                raise NotImplementedError("ISSF only available with INTERF='QUELCONQUE'")


class Generator:
    """Base class Generator"""

    @staticmethod
    def factory(macro, params):
        """create an instance of the appropriated type of Generator"""
        if params.cas == "TRANS":
            return GeneratorTRANS(macro, params)
        elif params.cas == "SPEC":
            return GeneratorSPEC(macro, params)
        else:
            raise ValueError("unknown configuration")

    def __init__(self, macro, params):
        """Constructor Base class"""
        self.macro = macro
        self.case = params.cas
        self.INFO = params.other_keys["INFO"]
        self.cohe_params = params.cohe_keys
        self.calc_params = params.other_keys
        self.mat_gene_params = params.mat_gene_keys
        self.interf_params = params.interf_keys
        self.NB_FREQ = None
        self.FREQ_PAS = None
        self.FREQ_INIT = None
        self.FMAX = None
        self.liste_freq = []
        self.liste_freq_sig = []
        self.cohe_params.update(params.interf_keys)
        self.NOM_CMP = None
        self.list_NOM_CMP = None

    def prepare_input(self):
        """run prepare data"""
        matr_rigi = self.mat_gene_params["MATR_RIGI"]
        bamo = matr_rigi.getModalBasis()
        ddlgene = matr_rigi.getGeneralizedDOFNumbering()
        self.mat_gene_params["NUME_DDL"] = ddlgene
        self.mat_gene_params["BASE"] = bamo
        mesh = bamo.getDOFNumbering().getMesh()
        nom_mail = mesh.getName()
        self.cohe_params["MAILLAGE"] = nom_mail
        nbmodt = bamo.getNumberOfIndexes()
        nbmodd = bamo.getNumberOfDynamicModes()
        nbmods = bamo.getNumberOfStaticModes()
        self.mat_gene_params["NBMODD"] = nbmodd
        self.mat_gene_params["NBMODS"] = nbmods
        self.mat_gene_params["NBMODT"] = nbmodt
        _, noe_interf = get_group_nom_coord(self.interf_params["GROUP_NO_INTERF"], mesh)
        self.cohe_params["DIST"] = calc_dist2(noe_interf)
        self.cohe_params["NOEUDS_INTERF"] = noe_interf
        self.interf_params["NBNO"] = len(noe_interf)

        if self.INFO == 2:
            texte1 = "NOMBRE DE MODES: " + str(nbmodt)
            texte2 = "MODES DYNAMIQUES: " + str(nbmodd)
            texte3 = "MODES STATIQUES: " + str(nbmods)
            aster.affiche("MESSAGE", texte1)
            aster.affiche("MESSAGE", texte2 + " - " + texte3)
            if self.case != "TRANS":
                aster.affiche("MESSAGE", "COMPOSANTE " + self.NOM_CMP)
            else:
                aster.affiche("MESSAGE", "COMPOSANTES " + ",".join(self.list_NOM_CMP))
            aster.affiche("MESSAGE", "NBNO INTERFACE: " + str(len(noe_interf)))

    def sampling(self):
        """run build_harm_gene"""
        raise NotImplementedError("must be implemented in a subclass")

    def build_result(self):
        """specific to each method: output table"""
        raise NotImplementedError("must be implemented in a subclass")

    def run(self):
        """compute SSI with spatial variability"""
        self.prepare_input()
        self.sampling()
        return self.build_result()

    def calc_miss_vari(self):
        """Compute SSI analysis with spatial variability"""

        NB_FREQ = 1 + int((self.FMAX - self.FREQ_INIT) / self.FREQ_PAS)
        RESU = [None] * len(self.list_NOM_CMP)
        list_NOM_CMP = self.list_NOM_CMP
        self.NOM_CMP = None

        frequencies = self.FREQ_INIT + self.FREQ_PAS * NP.arange(NB_FREQ)
        kwargs = self.cohe_params

        coherency_matrix = calc_coherency_matrix(
            frequencies=frequencies,
            model=kwargs["TYPE"],
            nom_mail=kwargs["MAILLAGE"],
            nom_group_inter=kwargs["GROUP_NO_INTERF"],
            **kwargs,
        )

        PVEC = compute_pod(
            coherency_matrix=coherency_matrix,
            precision=self.calc_params["PRECISION"],
            frequencies=frequencies,
            info=self.INFO,
        )

        # RECUPERATION DES MODES MECA (STATIQUES)
        nbmodd = self.mat_gene_params["NBMODD"]
        nbmods = self.mat_gene_params["NBMODS"]
        GROUP_NO_INTER = self.interf_params["GROUP_NO_INTERF"]

        # BOUCLE SUR LES DIRECTIONS
        for i_cmp, nom_cmp in enumerate(list_NOM_CMP):
            dict_modes = compute_mecmode(
                nom_cmp, GROUP_NO_INTER, self.mat_gene_params["BASE"], nbmods, nbmodd
            )
            dict_modes["nbmods"] = nbmods
            dict_modes["nbno"] = self.interf_params["NBNO"]
            dict_modes["NOM_CMP"] = nom_cmp
            # BOUCLE SUR LES FREQUENCES
            for k in range(NB_FREQ):
                dict_modes["nbpod"] = len(PVEC[k])
                # CALCUL ISS VARI
                if self.interf_params["MODE_INTERF"] != "QUELCONQUE":
                    RESU[i_cmp] = self.compute_freqk(k, RESU[i_cmp], PVEC[k], dict_modes)
                else:  # MODE_INTERF =='QUELCONQUE'
                    RESU[i_cmp] = self.compute_freqk_quelconque(k, RESU[i_cmp], PVEC[k], dict_modes)
        return RESU

    def compute_freqk(self):
        raise NotImplementedError("must be implemented in a subclass")

    def get_information_freqk(self, k, dict_modes, issf):
        nbmodt = self.mat_gene_params["NBMODT"]
        nbmodd = self.mat_gene_params["NBMODD"]
        freqk = self.FREQ_INIT + self.FREQ_PAS * k

        __impe = LIRE_IMPE_MISS(
            BASE=self.mat_gene_params["BASE"],
            TYPE=self.calc_params.get("TYPE"),
            NUME_DDL_GENE=self.mat_gene_params["NUME_DDL"],
            UNITE_RESU_IMPE=self.calc_params.get("UNITE_RESU_IMPE"),
            ISSF=issf,
            FREQ_EXTR=freqk,
        )

        __fosi = LIRE_FORC_MISS(
            BASE=self.mat_gene_params["BASE"],
            NUME_DDL_GENE=self.mat_gene_params["NUME_DDL"],
            NOM_CMP=dict_modes["NOM_CMP"],
            NOM_CHAM="DEPL",
            UNITE_RESU_FORC=self.calc_params["UNITE_RESU_FORC"],
            ISSF=issf,
            FREQ_EXTR=freqk,
        )

        __rito = COMB_MATR_ASSE(
            COMB_C=(
                _F(MATR_ASSE=__impe, COEF_C=1.0 + 0.0j),
                _F(MATR_ASSE=self.mat_gene_params["MATR_RIGI"], COEF_C=1.0 + 0.0j),
            ),
            SANS_CMP="LAGR",
        )

        return nbmodt, nbmodd, freqk, __impe, __fosi, __rito

    def compute_freqk_quelconque(self, k, RESU, VEC, dict_modes):
        """compute response for freqk - quelconque (trans and spec case)"""

        nbmodt, nbmodd, freqk, __impe, __fosi, __rito = self.get_information_freqk(
            k, dict_modes, self.calc_params["ISSF"]
        )
        MIMPE = __impe.toNumpy()
        #  extraction de la partie modes interface
        KRS = MIMPE[nbmodd:nbmodt, nbmodd:nbmodt]
        FSISM = __fosi.EXTR_VECT_GENE()
        FS0 = FSISM[nbmodd:nbmodt][:]
        FS = compute_corr_vari(dict_modes, VEC, KRS, FS0)
        FSISM[nbmodd:nbmodt][:] = FS
        __fosi.RECU_VECT_GENE(FSISM)
        if self.mat_gene_params["MATR_AMOR"] is not None:
            __dyge = DYNA_VIBRA(
                TYPE_CALCUL="HARM",
                BASE_CALCUL="GENE",
                MATR_MASS=self.mat_gene_params["MATR_MASS"],
                MATR_RIGI=__rito,
                FREQ=freqk,
                MATR_AMOR=self.mat_gene_params["MATR_AMOR"],
                EXCIT=_F(VECT_ASSE_GENE=__fosi, COEF_MULT=1.0),
            )
        else:
            __dyge = DYNA_VIBRA(
                TYPE_CALCUL="HARM",
                BASE_CALCUL="GENE",
                MATR_MASS=self.mat_gene_params["MATR_MASS"],
                MATR_RIGI=__rito,
                FREQ=freqk,
                EXCIT=_F(VECT_ASSE_GENE=__fosi, COEF_MULT=1.0),
            )
        #  recuperer le vecteur modal depl calcule par dyge
        RS = NP.array(__dyge.getDisplacement())
        VECRES = self.append_Vec(RS, k, RESU)
        return VECRES

    def compute_force_vari(self, dict_modes, VEC, *KRS):
        """compute seismic force with variability"""
        # CALCUL DE XO PAR PROJECTION-------------
        XO = NP.zeros(dict_modes["nbmods"])

        for mods in range(dict_modes["nbmods"]):
            MCMP = dict_modes["MCMP"][mods]
            som = dict_modes["som"][mods]
            maxm = dict_modes["maxm"][mods]
            is_som_bigger_than_zero = True
            if abs(som) < 10.0e-6:
                XO[mods] = 0.0
                is_som_bigger_than_zero = False

            #  CAS 1: MODES DE CORPS RIGIDE
            if self.interf_params["MODE_INTERF"] == "CORP_RIGI":
                # modes de translation
                if mods + 1 <= 3 and is_som_bigger_than_zero:
                    fact = 1.0 / som
                    XO[mods] = fact * abs(NP.inner(MCMP, VEC))
                # modes de rotation
                else:
                    if maxm < 10.0e-6:
                        if is_som_bigger_than_zero:
                            UTMESS("F", "ALGORITH6_86")
                    else:
                        fact = 1.0 / dict_modes["nbno"]
                        XO[mods] = 1.0 / (maxm**2.0) * fact * abs(NP.inner(MCMP, VEC))

            # CAS 2: MODES EF
            elif self.interf_params["MODE_INTERF"] == "TOUT":
                if is_som_bigger_than_zero:
                    fact = 1.0 / som
                    XO[mods] = fact * abs(NP.inner(MCMP, VEC))
                else:
                    if maxm < 10.0e-6:
                        UTMESS("F", "UTILITAI5_89")

        return NP.dot(KRS, NP.array(XO))
        # CAS 3: QUELCONQUE -> ROUTINE compute_corr_vari


#     -----------------------------------------------------------------
#         CLASSE TRANS
#     -----------------------------------------------------------------
class GeneratorTRANS(Generator):
    """TRANS class"""

    def __init__(self, macro, params):
        super().__init__(macro, params)
        self.excit_params = params.excit_sol_keys
        for dire in list(self.excit_params.keys()):
            if self.excit_params[dire] is None:
                del self.excit_params[dire]
        nom_cmp = [q.replace("ACCE_", "D") for q in list(self.excit_params.keys())]
        self.list_NOM_CMP = nom_cmp

    def sampling(self):
        """sampling for trans"""
        __foint = [None] * 3

        # verification que les abscisses des différents signaux sont les mêmes
        dire0 = list(self.excit_params.keys())[0]
        tt0, _ = self.excit_params[dire0].Valeurs()
        for idi, dire in enumerate(list(self.excit_params.keys())[1:]):
            tt, _ = self.excit_params[dire].Valeurs()
            if len(tt0) != len(tt):
                UTMESS("F", "SEISME_80", valk=[dire0, dire])
            max_abs = max(abs(NP.array(tt0)))
            prec = 1.0e-6
            for i in range(len(tt0)):
                if abs(tt0[i] - tt[i]) / max_abs > prec:
                    UTMESS("F", "SEISME_81", vali=[i + 1], valk=[dire0, dire])

            __foint[idi + 1] = CALC_FONCTION(
                FFT=_F(FONCTION=self.excit_params[dire], METHODE="PROL_ZERO")
            )
            self.excit_params[dire] = __foint[idi + 1]

        # utilisation du premier signal uniquement pour calcul des différents paramètres
        DT = tt0[1] - tt0[0]
        __foint[0] = CALC_FONCTION(FFT=_F(FONCTION=self.excit_params[dire0], METHODE="PROL_ZERO"))
        self.excit_params[dire0] = __foint[0]
        vale_fre, _, _ = __foint[0].Valeurs()
        FREQ_PAS = 1.0 / (len(vale_fre) * DT)
        NB_FREQ = int(floor(len(vale_fre) / 2))
        FREQ_INIT = 0.0
        FREQ_COUP = (NB_FREQ - 1) * FREQ_PAS
        FMAX = FREQ_COUP
        # signal nombre impair (floor(N/2)) ou signal nombre
        #       pair avec REST_SPEC_TEMP (prend N/2 pour N pair)
        # NB_FREQ= int(floor(len(vale_fre)/2)+1)
        # signal nombre pair: N/2+1
        # liste des frequences complete
        for k in range(NB_FREQ):
            self.liste_freq_sig.append(FREQ_INIT + FREQ_PAS * k)
        # si frequences de calcul donnees par l utilisateur
        FREQ_FIN = self.calc_params.get("FREQ_MAX")
        if FREQ_FIN is not None:
            if FREQ_FIN < FREQ_COUP:
                text = (
                    "FREQ_FIN ="
                    + str(FREQ_FIN)
                    + "Hz < "
                    + "FREQUENCE DE COUPURE ="
                    + str(FREQ_COUP)
                    + "Hz on complete par zero"
                )
                aster.affiche("MESSAGE", text)
            FREQ_PAS = self.calc_params.get("FREQ_PAS")
            NB_FREQ = int(ceil(FREQ_FIN / FREQ_PAS)) + 1
            FREQ_INIT = 0.0
            FMAX = (NB_FREQ - 1) * FREQ_PAS
        self.FREQ_PAS = FREQ_PAS
        self.NB_FREQ = NB_FREQ
        self.FREQ_INIT = FREQ_INIT
        self.FMAX = FMAX
        for k in range(self.NB_FREQ):
            self.liste_freq.append(self.FREQ_INIT + self.FREQ_PAS * k)
        if self.INFO == 2:
            text = (
                "DISCRETISATION UTILISATEUR : NB_FREQ, PAS, FREQ_COUP : "
                + str(self.NB_FREQ)
                + " ,"
                + str(self.FREQ_PAS)
                + " ,"
                + str(FREQ_COUP)
            )
            aster.affiche("MESSAGE", text)

    def build_result(self):
        return self.compute_result()

    def compute_result(self):
        L_VEC = self.compute_harm_gene()
        __dyge0 = self.create_host_sd()

        tup_re = []
        tup_im = []
        tup_re1 = []
        tup_im1 = []
        tup_re2 = []
        tup_im2 = []
        # boucle sur les 3 directions : sommation directionnelle
        n_dire = len(list(self.excit_params.keys()))
        for i, dire in enumerate(self.excit_params.keys()):
            VEC = L_VEC[i]

            #   ATTENTION:  on sort le champ en déplacement:
            #               c'est équivalent au champ en acceleration
            #               car on a applique un signal en ACCE pour fosi en acce
            # > cela evite de diviser par w2 pour integrer l'acceleration
            #               (erreurs numeriques au point 0)
            # > on remplace donc le champ en acceleration

            # si tous les point on été calculés: pas d'interpolation
            _, vale_re, vale_im = self.excit_params[dire].Valeurs()
            if self.calc_params.get("FREQ_MAX") is None:
                inul = 0
                for k, freqk in enumerate(self.liste_freq_sig):
                    omegk = 2.0 * pi * freqk
                    coef_a = vale_re[k] + vale_im[k] * 1.0j
                    VEC_comp = VEC[k] * coef_a
                    if i == 0:
                        tup_re.append(VEC_comp.real)
                        tup_im.append(VEC_comp.imag)
                    else:
                        tup_re[k] += VEC_comp.real
                        tup_im[k] += VEC_comp.imag
                    # traitement des vitesses et déplacements apres sommation des directions pour l'acceleration
                    if i == n_dire - 1:
                        if freqk > 1.0e-6:
                            tup_re1.append((-1.0) * tup_re[k] / (omegk * omegk))
                            tup_im1.append((-1.0) * tup_im[k] / (omegk * omegk))
                            tup_re2.append(tup_im[k] / omegk)
                            tup_im2.append((-1.0) * tup_re[k] / omegk)

                            if inul == 1:
                                inul = 2
                            elif inul == 2:
                                tup_re1[k - 2] = (2.0) * tup_re1[k - 1] - tup_re1[k]
                                tup_im1[k - 2] = (2.0) * tup_im1[k - 1] - tup_im1[k]
                                tup_re2[k - 2] = (2.0) * tup_re2[k - 1] - tup_re2[k]
                                tup_im2[k - 2] = (2.0) * tup_im2[k - 1] - tup_im2[k]
                                inul = 0
                        else:
                            inul = 1
                            tup_re1.append(VEC_comp.real * 0.0)
                            tup_im1.append(VEC_comp.imag * 0.0)
                            tup_re2.append(VEC_comp.real * 0.0)
                            tup_im2.append(VEC_comp.imag * 0.0)
            else:
                inul = 0
                for k, freqk in enumerate(self.liste_freq_sig):
                    coef_a = vale_re[k] + vale_im[k] * 1.0j
                    omegk = 2.0 * pi * freqk
                    if freqk >= self.calc_params.get("FREQ_MAX"):
                        # interpolation du vecteur POD VEC(NB_FREQ, nbmodt)
                        if i == 0:
                            tup_re.append(VEC[-1].real * 0.0)
                            tup_im.append(VEC[-1].imag * 0.0)
                    else:
                        vale_i = NP.searchsorted(self.liste_freq, freqk)
                        if vale_i == 0:
                            VEC_comp = VEC[0] * coef_a
                            if i == 0:
                                tup_re.append(VEC_comp.real)
                                tup_im.append(VEC_comp.imag)
                            else:
                                tup_re[k] += VEC_comp.real
                                tup_im[k] += VEC_comp.imag
                        else:
                            dfp = (freqk - self.liste_freq[vale_i - 1]) / (
                                self.liste_freq[vale_i] - self.liste_freq[vale_i - 1]
                            )
                            VEC_comp = coef_a * (
                                VEC[vale_i - 1] + dfp * (VEC[vale_i] - VEC[vale_i - 1])
                            )
                            if i == 0:
                                tup_re.append(VEC_comp.real)
                                tup_im.append(VEC_comp.imag)
                            else:
                                tup_re[k] += VEC_comp.real
                                tup_im[k] += VEC_comp.imag

                    # traitement des vitesses et déplacements apres sommation des directions pour l'acceleration
                    if i == n_dire - 1:
                        if freqk > 1.0e-6:
                            tup_re1.append((-1.0) * tup_re[k] / (omegk * omegk))
                            tup_im1.append((-1.0) * tup_im[k] / (omegk * omegk))
                            tup_re2.append(tup_im[k] / omegk)
                            tup_im2.append((-1.0) * tup_re[k] / omegk)

                            if inul == 1:
                                inul = 2
                            elif inul == 2:
                                tup_re1[k - 2] = (2.0) * tup_re1[k - 1] - tup_re1[k]
                                tup_im1[k - 2] = (2.0) * tup_im1[k - 1] - tup_im1[k]
                                tup_re2[k - 2] = (2.0) * tup_re2[k - 1] - tup_re2[k]
                                tup_im2[k - 2] = (2.0) * tup_im2[k - 1] - tup_im2[k]
                                inul = 0
                        else:
                            inul = 1
                            tup_re1.append(VEC_comp.real * 0.0)
                            tup_im1.append(VEC_comp.imag * 0.0)
                            tup_re2.append(VEC_comp.real * 0.0)
                            tup_im2.append(VEC_comp.imag * 0.0)

        # affectation des valeurs
        __dyge0.setDisplacement(NP.concatenate(tup_re1) + NP.concatenate(tup_im1) * 1j)
        __dyge0.setVelocity(NP.concatenate(tup_re2) + NP.concatenate(tup_im2) * 1j)
        __dyge0.setAcceleration(NP.concatenate(tup_re) + NP.concatenate(tup_im) * 1j)

        if self.INFO == 2:
            aster.affiche("MESSAGE", "START REST_SPEC_TEMP")

        dyha = REST_SPEC_TEMP(
            RESU_GENE=__dyge0,
            SYMETRIE="NON",
            METHODE="PROL_ZERO",
            TOUT_CHAM="OUI",
            N_PUIS=0,
            # NOM_CHAM='ACCE'
        )
        return dyha

    def create_host_sd(self):
        # on cree la SD resultat - factice
        # (le champ ACCE sera remplace dans la suite par celui calcule)
        __impe = LIRE_IMPE_MISS(
            BASE=self.mat_gene_params["BASE"],
            TYPE=self.calc_params.get("TYPE"),
            NUME_DDL_GENE=self.mat_gene_params["NUME_DDL"],
            UNITE_RESU_IMPE=self.calc_params.get("UNITE_RESU_IMPE"),
            ISSF=self.calc_params.get("ISSF"),
            FREQ_EXTR=self.FREQ_PAS,
        )
        __rito = COMB_MATR_ASSE(
            COMB_C=(
                _F(MATR_ASSE=__impe, COEF_C=1.0 + 0.0j),
                _F(MATR_ASSE=self.mat_gene_params["MATR_RIGI"], COEF_C=1.0 + 0.0j),
            ),
            SANS_CMP="LAGR",
        )
        # on fixe la composante car sa valeur n'a pas d'importance
        __fosi = LIRE_FORC_MISS(
            BASE=self.mat_gene_params["BASE"],
            NUME_DDL_GENE=self.mat_gene_params["NUME_DDL"],
            NOM_CMP="DX",
            NOM_CHAM="DEPL",
            ISSF=self.calc_params.get("ISSF"),
            UNITE_RESU_FORC=self.calc_params.get("UNITE_RESU_FORC"),
            FREQ_EXTR=self.FREQ_PAS,
        )

        __dyge0 = DYNA_VIBRA(
            TYPE_CALCUL="HARM",
            BASE_CALCUL="GENE",
            MATR_MASS=self.mat_gene_params["MATR_MASS"],
            MATR_RIGI=__rito,
            TOUT_CHAM="OUI",
            FREQ=self.liste_freq_sig,
            EXCIT=_F(VECT_ASSE_GENE=__fosi, COEF_MULT_C=1.0),
        )
        return __dyge0

    def append_Vec(self, RS, k, VEC=None):
        if VEC is None:
            nbmodt = self.mat_gene_params["NBMODT"]
            VEC = NP.zeros((self.NB_FREQ, nbmodt)) + 0j
        VEC[k] = RS
        return VEC

    #
    # Ecriture de  tran_gene
    # 1) on cree un concept harm_gene (factice) et le remplit
    # 2) On interpole les valeurs non calculés (liste_freq_sig)
    # 3) puis on fait la FFT pour obtenir le signal temporel

    def compute_harm_gene(self):
        """compute harm_gene for spec"""
        VEC = self.calc_miss_vari()
        return VEC

    def compute_freqk(self, k, RESU, VEC, dict_modes):
        """compute response for freqk - trans"""
        nbmodt, nbmodd, freqk, __impe, __fosi, __rito = self.get_information_freqk(
            k, dict_modes, "NON"
        )

        # IMPEDANCE
        MIMPE = __impe.toNumpy()
        #  extraction de la partie modes interface
        KRS = MIMPE[nbmodd:nbmodt, nbmodd:nbmodt]
        #  CALCUL FORCE SISMIQUE
        FSISM = __fosi.EXTR_VECT_GENE()
        #  extraction de la partie modes interface
        FS = 0.0
        for k1 in range(dict_modes["nbpod"]):
            FS = FS + self.compute_force_vari(dict_modes, VEC[k1], KRS)
        FSISM[nbmodd:nbmodt][:] = FS
        __fosi.RECU_VECT_GENE(FSISM)
        # CALCUL ISS
        if self.mat_gene_params["MATR_AMOR"] is not None:
            __dyge = DYNA_VIBRA(
                TYPE_CALCUL="HARM",
                BASE_CALCUL="GENE",
                MATR_MASS=self.mat_gene_params["MATR_MASS"],
                MATR_RIGI=__rito,
                FREQ=freqk,
                MATR_AMOR=self.mat_gene_params["MATR_AMOR"],
                EXCIT=_F(VECT_ASSE_GENE=__fosi, COEF_MULT=1.0),
            )
        else:
            __dyge = DYNA_VIBRA(
                TYPE_CALCUL="HARM",
                BASE_CALCUL="GENE",
                MATR_MASS=self.mat_gene_params["MATR_MASS"],
                MATR_RIGI=__rito,
                FREQ=freqk,
                EXCIT=_F(VECT_ASSE_GENE=__fosi, COEF_MULT=1.0),
            )
        #  recuperer le vecteur modal depl calcule par dyge
        RS = NP.array(__dyge.getDisplacement())
        VECRES = self.append_Vec(RS, k, RESU)
        return VECRES


#     -----------------------------------------------------------------
#         CLASSE SPEC
#     -----------------------------------------------------------------
class GeneratorSPEC(Generator):
    """SPEC class"""

    def __init__(self, macro, params):
        super().__init__(macro, params)
        self.NOM_CMP = params.other_keys["NOM_CMP"]
        self.list_NOM_CMP = [self.NOM_CMP]

    def sampling(self):
        """sampling for spec"""
        # discretisation temps et freq
        self.FREQ_INIT = self.calc_params.get("FREQ_INIT")
        self.FREQ_PAS = self.calc_params.get("FREQ_PAS")
        self.NB_FREQ = self.calc_params.get("NB_FREQ")
        self.FMAX = self.FREQ_INIT + self.FREQ_PAS * (self.NB_FREQ - 1)
        for k in range(self.NB_FREQ):
            self.liste_freq.append(self.FREQ_INIT + self.FREQ_PAS * k)

    def compute_harm_gene(self):
        """compute harm_gene for spec"""
        SPEC = self.calc_miss_vari()
        return SPEC[0]

    def compute_freqk(self, k, RESU, VEC, dict_modes):
        """compute response for freqk - spec"""
        nbmodt, nbmodd, freqk, __impe, __fosi, __rito = self.get_information_freqk(
            k, dict_modes, "NON"
        )
        # IMPEDANCE
        MIMPE = __impe.toNumpy()
        #  extraction de la partie modes interface
        KRS = MIMPE[nbmodd:nbmodt, nbmodd:nbmodt]
        #  CALCUL FORCE SISMIQUE AVEC VARIABILITE
        FSISM = __fosi.EXTR_VECT_GENE()
        SP = NP.zeros((nbmodt, nbmodt))
        for k1 in range(dict_modes["nbpod"]):
            #  calcul de la force sismique mode POD par mode POD
            FS = self.compute_force_vari(dict_modes, VEC[k1], KRS)
            FSISM[nbmodd:nbmodt][:] = FS
            #  Calcul harmonique
            __fosi.RECU_VECT_GENE(FSISM)
            if self.mat_gene_params["MATR_AMOR"] is not None:
                __dyge = DYNA_VIBRA(
                    TYPE_CALCUL="HARM",
                    BASE_CALCUL="GENE",
                    MATR_MASS=self.mat_gene_params["MATR_MASS"],
                    MATR_RIGI=__rito,
                    FREQ=freqk,
                    MATR_AMOR=self.mat_gene_params["MATR_AMOR"],
                    EXCIT=_F(VECT_ASSE_GENE=__fosi, COEF_MULT=1.0),
                )
            else:
                __dyge = DYNA_VIBRA(
                    TYPE_CALCUL="HARM",
                    BASE_CALCUL="GENE",
                    MATR_MASS=self.mat_gene_params["MATR_MASS"],
                    MATR_RIGI=__rito,
                    FREQ=freqk,
                    EXCIT=_F(VECT_ASSE_GENE=__fosi, COEF_MULT=1.0),
                )
            #  recuperer le vecteur modal depl calcule par dyge
            RS = NP.array(__dyge.getDisplacement())
            # stockage des matrices résultats: sum(s_q s_q* )
            SP = SP + RS * NP.conj(RS[:, NP.newaxis])
        SPEC = self.append_Vec(SP, k, RESU)
        return SPEC

    def compute_result(self):
        SPEC = self.compute_harm_gene()
        nbmodt = self.mat_gene_params["NBMODT"]
        mcfact = []
        for k2 in range(nbmodt):
            if self.calc_params.get("OPTION") == "DIAG":
                # on ecrit uniquement les termes diagonaux
                #        (autospectres) de la matrice
                foncc = []
                for k in range(self.NB_FREQ):
                    foncc.append(self.liste_freq[k])
                    foncc.append(SPEC[k][k2, k2].real)
                    foncc.append(SPEC[k][k2, k2].imag)
                _f = DEFI_FONCTION(NOM_PARA="FREQ", NOM_RESU="SPEC", VALE_C=foncc)
                mcfact.append(_F(NUME_ORDRE_I=k2 + 1, NUME_ORDRE_J=k2 + 1, FONCTION=_f))
            else:  # on ecrit tout
                for k1 in range(k2 + 1):
                    foncc = []
                    for k in range(self.NB_FREQ):
                        foncc.append(self.liste_freq[k])
                        foncc.append(SPEC[k][k1, k2].real)
                        foncc.append(SPEC[k][k1, k2].imag)
                    _f = DEFI_FONCTION(NOM_PARA="FREQ", NOM_RESU="SPEC", VALE_C=foncc)
                    mcfact.append(_F(NUME_ORDRE_I=k1 + 1, NUME_ORDRE_J=k2 + 1, FONCTION=_f))
        # Creation du concept en sortie
        dsp_out = DEFI_INTE_SPEC(PAR_FONCTION=mcfact, TITRE="DSP")
        return dsp_out

    def build_result(self):
        return self.compute_result()

    def append_Vec(self, RS, k, SPEC=None):
        if SPEC is None:
            nbmodt = self.mat_gene_params["NBMODT"]
            SPEC = NP.zeros((self.NB_FREQ, nbmodt, nbmodt)) + 0j
        if self.interf_params["MODE_INTERF"] == "QUELCONQUE":
            SPEC[k] = RS * NP.conj(RS[:, NP.newaxis])
        else:
            SPEC[k] = RS
        return SPEC
