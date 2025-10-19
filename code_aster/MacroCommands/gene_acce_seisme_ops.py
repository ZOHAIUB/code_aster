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


"""Commande GENE_ACCE_SEISME"""

import sys
import traceback
from datetime import datetime
from math import ceil, pi, sqrt

import numpy as NP

from ..Cata.Syntax import _F
from ..CodeCommands import CREA_TABLE, DEFI_FONCTION
from ..Messages import UTMESS
from ..Objects.function_py import t_fonction
from ..Objects.table_py import Table
from ..Utilities import disable_fpe
from .Utils.optimize import fmin
from .Utils.random_signal_utils import (
    DSP2ACCE1D,
    DSP2FR,
    RAND_DSP,
    RAND_VEC,
    SRO2DSP,
    Rice2,
    acce_filtre_CP,
    calc_dsp_KT,
    calc_phase_delay,
    corrcoefmodel,
    dsp_filtre_CP,
    f_ARIAS,
    f_ARIAS_TSM,
    f_opt1,
    f_opt2,
    f_opta,
    fonctm_gam,
    fonctm_JetH,
    gene_traj_gauss_evol1D,
    itersim_SRO,
    peak,
    zpa_match,
)
from .Utils.signal_correlation_utils import (
    CALC_CORRE,
    DSP2ACCE_ND,
    gene_traj_gauss_evol_ND,
    get_group_nom_coord,
    get_no_refe,
    itersimcor_SRO,
    itersimcortir_SRO,
)


def gene_acce_seisme_ops(self, **kwargs):
    """Corps de la macro GENE_ACCE_SEIMSE"""
    # conteneur des paramètres du calcul
    params = GeneAcceParameters(**kwargs)
    if not params.seed:
        now = datetime.now()
        params.seed = now.microsecond
        UTMESS("I", "SEISME_83", vali=params.seed)
    NP.random.seed(params.seed)
    # création de l'objet generator
    generator = Generator.factory(self, params)
    try:
        return generator.run()
    except Exception as err:
        trace = "".join(traceback.format_tb(sys.exc_info()[2]))
        UTMESS("F", "SUPERVIS2_3", valk=("GENE_ACCE_SEISME", trace, str(err)))


class GeneAcceParameters:
    def __init__(self, **kwargs):
        """Enregistrement des valeurs des mots-clés dans un dictionnaire.
        - On cree deux dictionnaires de parametres:
                      modulation_keys et simulation_keys
        """
        # GeneralKeys
        self.args = kwargs
        self.seed = kwargs.get("INIT_ALEA")
        self.norme = kwargs.get("PESANTEUR")
        #  # ModulationKeys
        modulation = kwargs.get("MODULATION")
        if modulation is not None and type(modulation) not in (list, tuple):
            modulation = (modulation,)
        modkeys = kwargs.get("MODULATION")[0]
        keys = self.modulation_keys = modkeys
        keys.update({"DUREE_PHASE_FORTE": kwargs.get("DUREE_PHASE_FORTE")})
        keys.update({"NORME": kwargs.get("PESANTEUR")})
        keys.update({"INFO": kwargs.get("INFO")})
        if "DSP" in kwargs:
            if keys["ECART_TYPE"]:
                keys["ECART_TYPE"] = keys["ECART_TYPE"] * self.norme
                if "ACCE_MAX" in keys:
                    del keys["ACCE_MAX"]
                if "INTE_ARIAS" in keys:
                    del keys["INTE_ARIAS"]
            elif keys["ACCE_MAX"]:
                keys["ACCE_MAX"] = keys["ACCE_MAX"] * self.norme
                if "ECART_TYPE" in keys:
                    del (keys["ECART_TYPE"],)
                if "INTE_ARIAS" in keys:
                    del keys["INTE_ARIAS"]
            elif keys["INTE_ARIAS"]:
                if "ECART_TYPE" in keys:
                    del (keys["ECART_TYPE"],)
                if "ACCE_MAX" in keys:
                    del keys["ACCE_MAX"]

        others = list(kwargs.keys())
        if "MODULATION" in others:
            others.remove("MODULATION")
        if "COEF_CORR" in others:
            others.remove("COEF_CORR")
        if "MATR_COHE" in others:
            others.remove("MATR_COHE")
        if "PHASE" in others:
            others.remove("PHASE")
        #  # SimulationKeys and MethodKeys
        if "COEF_CORR" in kwargs:
            corr_keys = {}
            corr_keys["TYPE"] = "COEF_CORR"
            corr_keys["COEF_CORR"] = kwargs.get("COEF_CORR")
            corr_keys["RATIO_HV"] = None
            if kwargs.get("SPEC_FRACTILE") is not None:
                corr_keys["RATIO_HV"] = kwargs.get("RATIO_HV")
        elif "MATR_COHE" in kwargs:
            ckeys = kwargs.get("MATR_COHE")[0]
            corr_keys = ckeys.cree_dict_valeurs(ckeys.mc_liste)
        elif "PHASE" in kwargs:
            ckeys = kwargs.get("PHASE")[0]
            corr_keys = ckeys.cree_dict_valeurs(ckeys.mc_liste)
            corr_keys["TYPE"] = "PHASE"
        else:
            corr_keys = {}
            corr_keys["TYPE"] = "SCALAR"
        self.simulation_keys = {"CORR_KEYS": corr_keys, "CORR_ZPA": False}

        if kwargs.get("DSP"):
            self.cas = "DSP"
            GeneratorKeys = kwargs.get("DSP")[0]
            method_keys = GeneratorKeys.cree_dict_valeurs(GeneratorKeys.mc_liste)
            if "DSP" in others:
                others.remove("DSP")
        else:
            self.cas = "SPECTRE"
            self.simulation_keys.update({"TYPE_ITER": "MOYENNE"})
            if kwargs.get("SPEC_FRACTILE"):
                GeneratorKeys = kwargs.get("SPEC_FRACTILE")[0]
                self.simulation_keys.update({"SPEC_METHODE": "SPEC_FRACTILE"})
                if "SPEC_FRACTILE" in others:
                    others.remove("SPEC_FRACTILE")
            elif kwargs.get("SPEC_MEDIANE"):
                GeneratorKeys = kwargs.get("SPEC_MEDIANE")[0]
                self.simulation_keys.update({"SPEC_METHODE": "SPEC_MEDIANE"})
                self.simulation_keys.update({"TYPE_ITER": "MEDIANE"})
                if "SPEC_MEDIANE" in others:
                    others.remove("SPEC_MEDIANE")
                if kwargs.get("NB_TIRAGE") == 1:
                    UTMESS("F", "SEISME_38")
            elif kwargs.get("SPEC_MOYENNE"):
                GeneratorKeys = kwargs.get("SPEC_MOYENNE")[0]
                self.simulation_keys.update({"SPEC_METHODE": "SPEC_MEDIANE"})
                if "SPEC_MOYENNE" in others:
                    others.remove("SPEC_MOYENNE")
                if kwargs.get("NB_TIRAGE") == 1:
                    UTMESS("F", "SEISME_38")
            elif kwargs.get("SPEC_UNIQUE"):
                GeneratorKeys = kwargs.get("SPEC_UNIQUE")[0]
                if "SPEC_UNIQUE" in others:
                    others.remove("SPEC_UNIQUE")
                self.simulation_keys.update({"SPEC_METHODE": "SPEC_UNIQUE"})

            method_keys = GeneratorKeys.cree_dict_valeurs(GeneratorKeys.mc_liste)

        self.method_keys = {}
        for key in method_keys:
            self.method_keys[key] = method_keys[key]
        if "NB_ITER" in self.method_keys:
            self.simulation_keys.update({"NB_ITER": self.method_keys["NB_ITER"]})
        if "CORR_ZPA" in self.method_keys:
            if self.method_keys["CORR_ZPA"] == "OUI":
                self.simulation_keys.update({"CORR_ZPA": True})
        # OtherKeys remplissage
        others_keys = {}
        for key in others:
            others_keys[key] = kwargs.get(key)
        self.simulation_keys.update(others_keys)
        self.simulation_keys.update({"CAS": self.cas})


class Generator:
    """Base class Generator"""

    @staticmethod
    def factory(macro, params):
        """create an instance of the appropriated type of Generator"""
        if params.cas == "DSP":
            return GeneratorDSP(macro, params)
        elif params.cas == "SPECTRE":
            return GeneratorSpectrum(macro, params)
        else:
            raise ValueError("unknown configuration")

    def __init__(self, macro, params):
        """Constructor Base class"""
        self.macro = macro
        self.norme = params.norme
        self.INFO = params.simulation_keys["INFO"]
        self.modul_params = params.modulation_keys
        self.method_params = params.method_keys
        self.simu_params = params.simulation_keys
        self.FREQ_FILTRE = params.simulation_keys.get("FREQ_FILTRE")
        self.FREQ_FILTRE_ZPA = params.simulation_keys.get("FREQ_FILTRE_ZPA")
        self.FREQ_CORNER = params.simulation_keys.get("FREQ_CORNER")
        self.FREQ_PENTE = params.simulation_keys.get("FREQ_PENTE")
        self.DSP_args = {}
        self.SRO_args = {"NORME": self.norme}
        self.tab = Table(titr="GENE_ACCE_SEISME")
        self.sampler = Sampler(params.modulation_keys, params.simulation_keys)
        # modulation indépendant de DSP/SPECTRE mais dépend de sampler:
        self.modulator = Modulator.factory(params.modulation_keys)
        # simulation depend de generator:
        self.simulator = Simulator.factory(params.simulation_keys)
        # parametres des t_fonctions  a creer
        # INTERPOL est modifié dans certains cas
        self.para_fonc_traj = {
            "NOM_PARA": "INST",
            "NOM_RESU": "ACCE",
            "PROL_DROITE": "EXCLU",
            "PROL_GAUCHE": "EXCLU",
            "TITRE": params.simulation_keys.get("TITRE"),
        }
        self.para_dsp = {
            "INTERPOL": ["LIN", "LIN"],
            "NOM_PARA": "FREQ",
            "PROL_DROITE": "CONSTANT",
            "PROL_GAUCHE": "EXCLU",
            "NOM_RESU": "ACCE",
        }
        self.para_sro = self.para_dsp

    def sampling(self):
        """run sampling"""
        self.sampler.run()

    def modulation(self):
        """run modulation"""
        self.modulator.run(self.sampler.liste_temps, self.sampler.DUREE_SIGNAL)

    def build_output(self):
        """run simulation"""
        self.simulator.run(self)

    def process_TimeHistory(self, Xt):
        """apply modulation and low pass filter if requested"""
        Xm = Xt * self.modulator.fonc_modul.vale_y
        if self.FREQ_FILTRE > 0.0:
            Xm = acce_filtre_CP(Xm, self.sampler.DT, self.FREQ_FILTRE)
        return Xm

    def prepare_data(self):
        """specific to each method"""
        raise NotImplementedError("must be implemented in a subclass")

    def build_DSP(self):
        """specific to each method"""
        raise NotImplementedError("must be implemented in a subclass")

    def build_result(self):
        """specific to each method: output table"""
        raise NotImplementedError("must be implemented in a subclass")

    def run(self):
        """Generate the signal"""
        self.sampling()
        self.modulation()
        self.prepare_data()
        self.build_DSP()
        return self.build_result()


class GeneratorDSP(Generator):
    """DSP class"""

    def prepare_data(self):
        """prepare data for DSP class"""
        self.DSP_args.update(
            {
                "FREQ_FOND": self.method_params.get("FREQ_FOND"),
                "AMORT": self.method_params.get("AMOR_REDUIT"),
            }
        )
        if self.FREQ_CORNER is None:
            self.FREQ_CORNER = 0.05 * self.DSP_args["FREQ_FOND"]
        # Il faut calculer le facteur de pic si la donnee = PGA
        # pour obtenir sigma et multiplier la modulation avec cette valeur
        if "ACCE_MAX" in self.modul_params:
            PeakFactor = self.calc_PeakFactor()
            sigma = self.modul_params["ACCE_MAX"] / PeakFactor
            self.modulator.sigma = sigma
            f_mod = t_fonction(
                self.sampler.liste_temps,
                self.modulator.fonc_modul.vale_y * sigma,
                para=self.modulator.para_fonc_modul,
            )
            self.modulator.fonc_modul = f_mod
            if self.INFO == 2:
                UTMESS("I", "PROBA0_12", valr=(PeakFactor))

    def calc_PeakFactor(self):
        """calculate peak factor"""
        spec = calc_dsp_KT(self, self.DSP_args["FREQ_FOND"], self.DSP_args["AMORT"])
        _, _, _, vop, delta = Rice2(self.sampler.liste_w2, spec)
        nup = peak(0.5, self.sampler.DUREE_PHASE_FORTE, vop, delta)
        return nup

    def build_DSP(self):
        """build DSP for DSP class"""
        if self.FREQ_PENTE is not None:
            self.DSP_args.update(
                {"FREQ_CORNER": self.FREQ_CORNER, "FREQ_PENTE": self.FREQ_PENTE, "TYPE_DSP": "KT"}
            )
        else:
            # calcul du facteur de normalisation
            dsp = calc_dsp_KT(self, self.DSP_args["FREQ_FOND"], self.DSP_args["AMORT"])
            # constante de normalisation pour que ecart_type=1:
            S_cst = 1.0 / (NP.trapz(dsp, self.sampler.liste_w2) * 2.0)
            # calcul de la DSP KT
            vale_dsp_KT = calc_dsp_KT(
                self, self.DSP_args["FREQ_FOND"], self.DSP_args["AMORT"], S_cst
            )
            fonc_dsp = t_fonction(self.sampler.liste_w2, vale_dsp_KT, para=self.para_dsp)
            self.DSP_args.update({"FONC_DSP": fonc_dsp, "TYPE_DSP": "KT"})

    def build_result(self):
        """Create the result function"""
        # --- construction des fonctions sortie
        self.build_output()
        # --- Creation du concept (table) en sortie
        dict_keywords = self.tab.dict_CREA_TABLE()
        tab_out = CREA_TABLE(TYPE_TABLE="TABLE_FONCTION", **dict_keywords)
        for func in self.tab.referenceToDataStructure:
            tab_out.addDependency(func)
        return tab_out


class GeneratorSpectrum(Generator):
    """Response Spectra class"""

    def prepare_data(self):
        """prepare data for Spectrum class"""
        if self.FREQ_CORNER is None:
            self.FREQ_CORNER = 0.0

        self.SRO_args.update({"FREQ_FILTRE_ZPA": self.FREQ_FILTRE_ZPA})

        if "NB_ITER" in self.method_params:
            dico_err = {
                "ERRE_ZPA": list(self.method_params.get("ERRE_ZPA")),
                "ERRE_MAX": list(self.method_params.get("ERRE_MAX")),
                "ERRE_RMS": list(self.method_params.get("ERRE_RMS")),
            }
            err_def = 0.2
            for keys in dico_err:
                if len(dico_err[keys]) < 2:
                    dico_err[keys].append(err_def)
            self.SRO_args.update({"DICO_ERR": dico_err, "NB_ITER": self.simu_params["NB_ITER"]})
            if self.simu_params["TYPE_ITER"] == "MEDIANE":
                self.SRO_args.update({"TYPE_ITER": "SPEC_MEDIANE"})
            elif self.simu_params["TYPE_ITER"] == "MOYENNE":
                self.SRO_args.update({"TYPE_ITER": "SPEC_MOYENNE"})

        spec_osci = self.method_params.get("SPEC_OSCI")
        l_freq_sro, sro_ref = spec_osci.Valeurs()
        ZPA = sro_ref[-1]
        self.SRO_args.update({"ZPA": ZPA})

        F_MIN = l_freq_sro[0]
        if self.sampler.FREQ_COUP > l_freq_sro[-1]:
            sro_ref = NP.append(sro_ref, ZPA)
            l_freq_sro = NP.append(l_freq_sro, self.sampler.FREQ_COUP)

        para_sro = dict(self.para_sro)
        para_sro["INTERPOL"] = spec_osci.Parametres()["INTERPOL"]
        f_spec = t_fonction(l_freq_sro, sro_ref, para=para_sro)

        self.SRO_args.update(
            {"FONC_SPEC": f_spec, "FMIN": F_MIN, "AMORT": self.method_params.get("AMOR_REDUIT")}
        )
        if "METHODE" in self.method_params:
            self.SRO_args.update({"METHODE_SRO": self.method_params.get("METHODE")})
        if "SPEC_1_SIGMA" in self.method_params:
            spec_sigma = self.method_params.get("SPEC_1_SIGMA")
            para_sro = dict(self.para_sro)
            para_sro["INTERPOL"] = spec_sigma.Parametres()["INTERPOL"]
            f_spec_sigma = t_fonction(spec_sigma.Absc(), spec_sigma.Ordo(), para=para_sro)
            f_spec_sigma = f_spec_sigma.evalfonc(l_freq_sro)
            sro_beta = NP.log(f_spec_sigma.vale_y / sro_ref)
            f_beta = t_fonction(l_freq_sro, sro_beta, para=para_sro)
            self.SRO_args.update({"FONC_BETA": f_beta})
        if "FREQ_PAS" in self.method_params:
            self.SRO_args.update({"PAS": self.method_params.get("FREQ_PAS")})
        elif "LIST_FREQ" in self.method_params:
            L_FREQ = self.method_params.get("LIST_FREQ").getValues()
            assert L_FREQ[0] > 0.0, "LIST_FREQ: il faut des valeurs >0.0"
            self.SRO_args.update({"LIST_FREQ": L_FREQ})
        else:
            self.SRO_args.update({"PAS": self.sampler.DW / 2.0 / pi})

    def build_DSP(self):
        """build DSP for Spectrum class"""
        #  CALCUL DE LA DSP SPECTRUM-COMPATIBLE

        f_dsp, f_spec_ref = SRO2DSP(
            self.sampler.FREQ_COUP, self.sampler.DUREE_PHASE_FORTE, **self.SRO_args
        )
        if self.FREQ_CORNER > 0.0:
            f_dsp = dsp_filtre_CP(f_dsp, self.FREQ_CORNER)
        #        if self.FREQ_FILTRE_ZPA > 0.0:
        #            f_dsp = butterfilter(self.FREQ_FILTRE_ZPA, f_dsp)

        fonc_dsp = f_dsp.evalfonc(self.sampler.liste_w2)
        self.DSP_args.update({"FONC_DSP": fonc_dsp, "TYPE_DSP": "SC", "FC": 0.05})
        self.SRO_args["FONC_SPEC"] = f_spec_ref

        if self.FREQ_PENTE is not None:
            self.DSP_args["TYPE_DSP"] = "FR"
            vop, amo, R0, R2, f_FIT = DSP2FR(self.DSP_args["FONC_DSP"], self.DSP_args["FC"])
            self.DSP_args.update(
                {
                    "FREQ_PENTE": self.FREQ_PENTE,
                    "FREQ_FOND": vop,
                    "AMORT": amo,
                    "para_R0": R0,
                    "para_R2": R2,
                    "fonc_FIT": f_FIT,
                    "TYPE_DSP": "FR",
                }
            )
        if self.simu_params["SPEC_METHODE"] == "SPEC_FRACTILE":
            Periods = 1.0 / (self.sampler.liste_w2 / (2.0 * pi))
            Periods, MAT_COVC = corrcoefmodel(Periods, self.SRO_args["FONC_BETA"])
            self.DSP_args.update({"PERIODS": Periods, "MAT_COVC": MAT_COVC})

    def build_result(self):
        """Create the result function"""
        # --- construction des fonctions sortie
        self.build_output()
        # --- Creation du concept (table) en sortie
        dict_keywords = self.tab.dict_CREA_TABLE()
        tab_out = CREA_TABLE(TYPE_TABLE="TABLE_FONCTION", **dict_keywords)
        for func in self.tab.referenceToDataStructure:
            tab_out.addDependency(func)
        return tab_out


#  ------------------------------------------------------------------
#  ECHANTILLONNAGE
#  ------------------------------------------------------------------


class Sampler:
    """class Sampling: common task for all cases"""

    def __init__(self, modul_params, method_params):
        self.FREQ_FILTRE = method_params.get("FREQ_FILTRE")
        self.INFO = modul_params["INFO"]
        self.DT = method_params.get("PAS_INST")
        self.modulation_type = modul_params["TYPE"]
        self.DUREE_PHASE_FORTE = modul_params["DUREE_PHASE_FORTE"]
        self.INST_INI = 0.0
        if "INST_INI" in modul_params:
            self.INST_INI = modul_params["INST_INI"]
        self.FREQ_COUP = 1.0 / (2.0 * self.DT)
        self.liste_temps = None
        self.DUREE_SIGNAL = None
        self.liste_w = None
        self.liste_w2 = None
        self.DW = None
        if method_params.get("NB_POIN") is None:
            self.NB_POIN = None
        else:
            self.NB_POIN = int(method_params.get("NB_POIN"))

    def run(self):
        """compute sampling"""
        # discretisation temps et freq
        OM = pi / self.DT

        if self.modulation_type == "CONSTANT":
            # on simule uniquement la phase forte si CONSTANT
            TTS = self.DUREE_PHASE_FORTE
            # on prend NB_POIN pair uniquement
            self.NB_POIN = int(ceil((TTS / self.DT + 1) / 2.0) * 2.0)
            DW = 2.0 * OM / self.NB_POIN
            TT = (self.NB_POIN - 1) * self.DT
        # on calcule la duree de simulation si NB_POIN donne
        elif self.NB_POIN is not None:
            if self.NB_POIN % 2 != 0:
                self.NB_POIN = self.NB_POIN + 1
            TT = (self.NB_POIN - 1) * self.DT
            DW = 2.0 * OM / self.NB_POIN
            if TT < self.DUREE_PHASE_FORTE * 1.5:
                UTMESS("A", "SEISME_39", valk=(str(TT)))

        else:  # on prend 3* phase forte comme duree de simulation
            TTS = self.INST_INI + 3.0 * self.DUREE_PHASE_FORTE
            # on prend NB_POIN pair uniquement
            self.NB_POIN = int(ceil((TTS / self.DT + 1) / 2.0) * 2.0)
            DW = 2.0 * OM / self.NB_POIN
            TT = (self.NB_POIN - 1) * self.DT

        liste_temps = NP.arange(0.0, self.NB_POIN * self.DT, self.DT)
        l_w = NP.arange(-OM + DW / 2.0, OM + DW / 2.0, DW)
        l_w2 = NP.arange(DW / 2.0, OM + DW / 2.0, DW)
        # parfois les listes ne sont pas bien construites
        # pour cause d'erreur num si valeurs reeles
        liste_temps = liste_temps[0 : self.NB_POIN]
        l_w = l_w[0 : self.NB_POIN]
        l_w2 = l_w2[0 : self.NB_POIN // 2]
        nbfreq = 2 * len(l_w2)
        assert self.NB_POIN == nbfreq
        assert len(liste_temps) == self.NB_POIN
        assert len(l_w) == self.NB_POIN
        self.liste_temps = liste_temps
        self.liste_w = l_w
        self.liste_w2 = l_w2
        self.DUREE_SIGNAL = TT
        self.DW = DW
        if self.INFO == 2:
            if self.FREQ_FILTRE > 0.0:
                vale_filtre = str(self.FREQ_FILTRE) + " Hz"
            else:
                vale_filtre = "None"
            valargs = _F(
                vali=self.NB_POIN,
                valr=(self.FREQ_COUP, self.DW / 2.0 / pi, self.DT, self.DUREE_SIGNAL),
                valk=vale_filtre,
            )
            UTMESS("I", "SEISME_9", **valargs)


#     -----------------------------------------------------------------
#          MODULATION   Gamma, JH, Constant
#     -----------------------------------------------------------------


class Modulator:
    """class Modulator: common task for all cases"""

    @staticmethod
    def factory(modul_params):
        """create an instance of the Modulator"""
        if modul_params["TYPE"] == "GAMMA":
            return ModulatorGamma(modul_params)
        elif modul_params["TYPE"] == "JENNINGS_HOUSNER":
            return ModulatorJH(modul_params)
        elif modul_params["TYPE"] == "CONSTANT":
            return ModulatorConstant(modul_params)
        else:
            raise ValueError("unknown configuration")

    def __init__(self, modul_params):
        self.para_fonc_modul = {
            "NOM_PARA": "INST",
            "NOM_RESU": "ACCE",
            "INTERPOL": ["LIN", "LIN"],
            "PROL_DROITE": "EXCLU",
            "PROL_GAUCHE": "EXCLU",
        }
        self.modul_params = modul_params
        self.DUREE_PHASE_FORTE = modul_params["DUREE_PHASE_FORTE"]
        self.norme = modul_params["NORME"]
        self.sigma = 1.0  # > il faut le calculer plus tard (GeneratorDSP) avec le facteur de pic.
        self.fonc_modul = None
        self.T1 = 0.0
        self.T2 = self.DUREE_PHASE_FORTE

    def run(self, sample_time, DUREE_SIGNAL):
        """Generate Modulating function: specific to each method"""
        raise NotImplementedError("must be implemented in a subclass")

    def calc_fonc_modul(self, sample_time, N1, N2, fqt):
        """determine amplitude of modulating function fqt"""
        if "INTE_ARIAS" in self.modul_params:
            vale_arias = f_ARIAS(sample_time, fqt, self.norme)
            fqt = fqt * sqrt(self.modul_params["INTE_ARIAS"] / vale_arias)
        elif "ECART_TYPE" in self.modul_params:
            int12 = NP.trapz((fqt[N1:N2]) ** 2, sample_time[N1:N2])
            fqt = fqt * self.modul_params["ECART_TYPE"] * sqrt(self.DUREE_PHASE_FORTE / int12)
        elif "ACCE_MAX" in self.modul_params:
            int12 = NP.trapz(fqt[N1:N2] ** 2, sample_time[N1:N2])
            fqt = fqt * self.sigma * sqrt(self.DUREE_PHASE_FORTE / int12)
        else:
            # equivalence energie totale avec signal module par CONSTANT sur
            # DUREE
            int12 = NP.trapz(fqt**2, sample_time)
            fqt = fqt * sqrt(self.DUREE_PHASE_FORTE / int12)
        f_mod = t_fonction(sample_time, fqt, para=self.para_fonc_modul)
        self.fonc_modul = f_mod


class ModulatorGamma(Modulator):
    """Modulator type Gamma"""

    def run(self, sample_time, DUREE_SIGNAL):
        T1 = self.modul_params["INST_INI"]
        T2 = T1 + self.DUREE_PHASE_FORTE
        x0 = [1.3, 0.25]
        liste_t = NP.arange(0.0, DUREE_SIGNAL, 0.01)
        N1 = NP.searchsorted(liste_t, T1)
        N2 = NP.searchsorted(liste_t, T2)
        fqt_ini = fonctm_gam(liste_t, 1.0, x0[0], x0[1])
        _, TSM, _, _ = f_ARIAS_TSM(liste_t, fqt_ini, self.norme)
        x_opt = fmin(f_opta, x0, args=(liste_t, N1, N2))
        a2 = x_opt[0]
        a3 = x_opt[1]
        fqt = fonctm_gam(sample_time, 1.0, a2, a3)
        _, TSM, self.T1, self.T2 = f_ARIAS_TSM(sample_time, fqt, self.norme)
        self.calc_fonc_modul(sample_time, N1, N2, fqt)
        if self.modul_params["INFO"] == 2:
            valargs = _F(valk=("GAMMA", str(a2) + " " + str(a3)), valr=(TSM, self.T1, self.T2))
            UTMESS("I", "SEISME_44", **valargs)


class ModulatorJH(Modulator):
    """Modulator type JH"""

    def run(self, sample_time, DUREE_SIGNAL):
        x0 = [0.5, 1.0]
        t1_ini = [2.0]
        #      identify T1 such that Inta1=0.05*Arias
        liste_t = NP.arange(0.0, DUREE_SIGNAL, 0.01)
        t_opt = fmin(f_opt1, t1_ini, args=(liste_t, self.DUREE_PHASE_FORTE, 0.5, 1.0))
        T1 = t_opt[0]
        x_opt = fmin(f_opt2, x0, args=(liste_t, T1, self.DUREE_PHASE_FORTE))
        alpha = x_opt[0]
        beta = x_opt[1]
        T2 = T1 + self.DUREE_PHASE_FORTE
        N1 = NP.searchsorted(liste_t, T1)
        N2 = NP.searchsorted(liste_t, T2)
        fqt = fonctm_JetH(sample_time, T1, T2, alpha, beta)
        _, TSM, self.T1, self.T2 = f_ARIAS_TSM(sample_time, fqt, self.norme)
        self.calc_fonc_modul(sample_time, N1, N2, fqt)
        if self.modul_params["INFO"] == 2:
            valargs = _F(
                valk=("JENNINGS & HOUSNER", str(alpha) + " " + str(beta)),
                valr=(TSM, self.T1, self.T2),
            )
            UTMESS("I", "SEISME_44", **valargs)


class ModulatorConstant(Modulator):
    """Modulator type Constant"""

    def run(self, sample_time, DUREE_SIGNAL):
        if "INTE_ARIAS" in self.modul_params:
            vale_arias = self.DUREE_PHASE_FORTE * pi / (2.0 * self.norme)
            fq = sqrt(self.modul_params["INTE_ARIAS"] / vale_arias)
        elif "ECART_TYPE" in self.modul_params:
            fq = self.modul_params["ECART_TYPE"]
        elif "ACCE_MAX" in self.modul_params:
            fq = self.sigma
        else:
            fq = 1.0
        fqt = NP.array([fq] * len(sample_time))
        f_mod = t_fonction(sample_time, fqt, para=self.para_fonc_modul)
        self.fonc_modul = f_mod
        if self.modul_params["INFO"] == 2:
            valargs = _F(
                valk=("CONSTANTE", "None"), valr=(self.DUREE_PHASE_FORTE, self.T1, self.T2)
            )
            UTMESS("I", "SEISME_44", **valargs)


#     -----------------------------------------------------------------
#          SIMULATION
#     -----------------------------------------------------------------


class Simulator:
    """class Simulation"""

    @staticmethod
    def factory(simu_params):
        """create an instance of the simulator"""
        if simu_params["CAS"] == "DSP":
            if simu_params["CORR_KEYS"]["TYPE"] == "PHASE":
                return SimulatorDSPPhase(simu_params)
            elif simu_params["CORR_KEYS"]["TYPE"] == "SCALAR":
                return SimulatorDSPScalar(simu_params)
            else:
                return SimulatorDSPVector(simu_params)
        #            else:
        #                raise ValueError('unknown configuration')
        elif simu_params["CAS"] == "SPECTRE":
            if simu_params["CORR_KEYS"]["TYPE"] == "PHASE":
                return SimulatorSPECPhase(simu_params)
            elif simu_params["CORR_KEYS"]["TYPE"] == "SCALAR":
                return SimulatorSPECScalar(simu_params)
            else:
                return SimulatorSPECVector(simu_params)
        #            else:
        #                raise ValueError('unknown configuration')
        else:
            raise ValueError("unknown configuration")

    def __init__(self, simu_params):
        self.simu_params = simu_params
        self.ntir = 0
        self.TYPE = simu_params["CORR_KEYS"]["TYPE"]
        self.DEFI_COHE = simu_params["CORR_KEYS"]
        self.INFO = simu_params["INFO"]
        self.nbtirage = simu_params["NB_TIRAGE"]
        self.FREQ_FILTRE = simu_params["FREQ_FILTRE"]
        self.para_fonc_traj = {
            "NOM_PARA": "INST",
            "NOM_RESU": "ACCE",
            "PROL_DROITE": "EXCLU",
            "PROL_GAUCHE": "EXCLU",
            "TITRE": simu_params.get("TITRE"),
        }

    def process_TimeHistory(self, generator, Xt):
        """apply modulation and low pass filter if requested"""
        Xm = Xt * generator.modulator.fonc_modul.vale_y

        if self.FREQ_FILTRE > 0.0:
            Xm = acce_filtre_CP(Xm, generator.sampler.DT, self.FREQ_FILTRE)

        if self.simu_params["CORR_ZPA"]:
            zpa = generator.SRO_args["ZPA"] * generator.SRO_args["NORME"]
            Xm = zpa_match((generator.sampler.liste_temps, Xm), zpa)

            if self.FREQ_FILTRE > 0.0:
                Xm = acce_filtre_CP(Xm, generator.sampler.DT, self.FREQ_FILTRE)

        return Xm

    def build_TimeHistory(self):
        """build_TimeHistory: specific to scalar or vector case"""
        raise NotImplementedError("must be implemented in a subclass")

    def run(self, generator):
        """run simulator: specific to scalar or vector case"""
        raise NotImplementedError("must be implemented in a subclass")


class SimulatorDSPScalar(Simulator):
    """Construct scalar signal for DSP class"""

    def build_TimeHistory(self, generator):
        """build scalar Time History for DSP class"""
        if self.INFO == 2:
            UTMESS("I", "PROBA0_13", vali=self.ntir + 1)
        if "FREQ_PENTE" in generator.DSP_args:
            Xt = gene_traj_gauss_evol1D(generator, **generator.DSP_args)
        else:
            Xt = DSP2ACCE1D(generator.DSP_args["FONC_DSP"])
        return Xt

    def run(self, generator):
        """Create the result table of functions"""
        for _ in range(self.nbtirage):
            Xt = self.build_TimeHistory(generator)
            Xt = self.process_TimeHistory(generator, NP.array(Xt))
            _f_out = DEFI_FONCTION(
                ABSCISSE=tuple(generator.sampler.liste_temps),
                ORDONNEE=tuple(Xt),
                **self.para_fonc_traj
            )
            generator.tab.append({"NUME_ORDRE": self.ntir + 1, "FONCTION": _f_out.getName()})
            generator.tab.referenceToDataStructure.append(_f_out)
            self.ntir = self.ntir + 1


class SimulatorDSPVector(Simulator):
    """Construct vector valued signal with correlation matrix for DSP class"""

    def build_TimeHistory(self, generator):
        """build vector valued Time History for DSP class"""
        if self.TYPE == "COEF_CORR":
            rho = self.DEFI_COHE["COEF_CORR"]
            dim = 2
            with disable_fpe():
                Mat_cohe = NP.linalg.cholesky(CALC_CORRE(rho, dim))
            Data_cohe = self.DEFI_COHE
            Data_cohe.update({"MAT_COHE": Mat_cohe})
        else:
            Data_cohe = self.DEFI_COHE

        if self.INFO == 2:
            UTMESS("I", "PROBA0_13", vali=self.ntir + 1)
        if "FREQ_PENTE" in generator.DSP_args:
            Xt = gene_traj_gauss_evol_ND(generator, Data_cohe, **generator.DSP_args)
        else:
            Xt = DSP2ACCE_ND(generator.DSP_args["FONC_DSP"], Data_cohe)
        return Xt

    def run(self, generator):
        """build result for vector DSP class"""
        if self.TYPE != "COEF_CORR":
            liste_nom, l2 = get_group_nom_coord(
                self.DEFI_COHE["GROUP_NO_INTERF"], self.DEFI_COHE["MAILLAGE"]
            )
            self.DEFI_COHE.update({"NOEUDS_INTERF": l2})
        for _ in range(self.nbtirage):
            Xt = self.build_TimeHistory(generator)
            nba = 1
            for acce in Xt:
                accef = self.process_TimeHistory(generator, NP.array(acce))
                _f_out = DEFI_FONCTION(
                    ABSCISSE=tuple(generator.sampler.liste_temps),
                    ORDONNEE=tuple(accef),
                    **self.para_fonc_traj
                )
                generator.tab.referenceToDataStructure.append(_f_out)
                if self.TYPE == "COEF_CORR":
                    nom_acce = "ACCE" + str(nba)
                    generator.tab.append(
                        {
                            "NUME_ORDRE": self.ntir + 1,
                            "FONCTION": _f_out.getName(),
                            "NOM_PARA": nom_acce,
                        }
                    )
                else:
                    nom_no = liste_nom[nba - 1]
                    generator.tab.append(
                        {"NUME_ORDRE": self.ntir + 1, "FONCTION": _f_out.getName(), "NOEUD": nom_no}
                    )
                nba = nba + 1
            self.ntir = self.ntir + 1


class SimulatorSPECVector(Simulator):
    """Construct vector valued signal with correlation matrix for SPEC class"""

    def run(self, generator):
        """build result for vector SPEC class"""
        if self.TYPE != "COEF_CORR":
            self.liste_nom, l2 = get_group_nom_coord(
                self.DEFI_COHE["GROUP_NO_INTERF"], self.DEFI_COHE["MAILLAGE"]
            )
            self.DEFI_COHE.update({"DIM": len(self.liste_nom)})
            self.DEFI_COHE.update({"NOEUDS_INTERF": l2})

        if self.simu_params["SPEC_METHODE"] == "SPEC_MEDIANE" and "NB_ITER" in self.simu_params:
            self.build_TimeHistories(generator)
        else:
            for _ in range(self.nbtirage):
                Xt = self.build_TimeHistory(generator)
                nba = 1
                for acce in Xt:
                    accef = self.process_TimeHistory(generator, NP.array(acce))
                    _f_out = DEFI_FONCTION(
                        ABSCISSE=tuple(generator.sampler.liste_temps),
                        ORDONNEE=tuple(accef),
                        **self.para_fonc_traj
                    )
                    generator.tab.referenceToDataStructure.append(_f_out)
                    if self.TYPE == "COEF_CORR":
                        nom_acce = "ACCE" + str(nba)
                        generator.tab.append(
                            {
                                "NUME_ORDRE": self.ntir + 1,
                                "FONCTION": _f_out.getName(),
                                "NOM_PARA": nom_acce,
                            }
                        )
                    else:
                        nom_no = self.liste_nom[nba - 1]
                        generator.tab.append(
                            {
                                "NUME_ORDRE": self.ntir + 1,
                                "FONCTION": _f_out.getName(),
                                "NOEUD": nom_no,
                            }
                        )
                    nba = nba + 1
                self.ntir = self.ntir + 1

    def build_TimeHistory(self, generator):
        """build vector valued Time History for Spectrum class"""
        specmethode = self.simu_params["SPEC_METHODE"]
        DSP_args = generator.DSP_args

        if self.TYPE == "COEF_CORR":
            rho = self.DEFI_COHE["COEF_CORR"]
            with disable_fpe():
                if self.DEFI_COHE["RATIO_HV"] is not None:
                    dim = 3
                    RS = (1.0 / self.DEFI_COHE["RATIO_HV"]) ** 2
                    Mat_cohe = NP.linalg.cholesky(CALC_CORRE(rho, dim, RS))
                else:
                    dim = 2
                    Mat_cohe = NP.linalg.cholesky(CALC_CORRE(rho, dim))
            Data_cohe = self.DEFI_COHE
            Data_cohe.update({"MAT_COHE": Mat_cohe, "DIM": dim})
        else:
            Data_cohe = self.DEFI_COHE

        if self.INFO == 2:
            UTMESS("I", "PROBA0_13", vali=self.ntir + 1)
        if specmethode == "SPEC_UNIQUE":
            if "NB_ITER" not in self.simu_params:
                if "FREQ_PENTE" in DSP_args:
                    Xt = gene_traj_gauss_evol_ND(generator, Data_cohe, **DSP_args)
                else:
                    Xt = DSP2ACCE_ND(DSP_args["FONC_DSP"], Data_cohe)
            else:  #'NB_ITER' in self.method_params
                Xt = []
                if "FREQ_PENTE" in DSP_args:
                    fonc_dsp_opt, liste_rv = itersimcor_SRO(
                        generator, DSP_args["FONC_DSP"], Data_cohe, **generator.SRO_args
                    )
                    vop, amo, R0, R2, f_FIT = DSP2FR(fonc_dsp_opt, DSP_args["FC"])
                    DSP_args.update(
                        {
                            "FREQ_FOND": vop,
                            "AMORT": amo,
                            "para_R0": R0,
                            "para_R2": R2,
                            "fonc_FIT": f_FIT,
                        }
                    )
                    Xt = gene_traj_gauss_evol_ND(generator, Data_cohe, rv=liste_rv, **DSP_args)
                else:
                    fonc_dsp_opt, liste_rv = itersimcor_SRO(
                        generator, DSP_args["FONC_DSP"], Data_cohe, **generator.SRO_args
                    )
                    Xt = DSP2ACCE_ND(fonc_dsp_opt, Data_cohe, rv=liste_rv)

        if specmethode == "SPEC_MEDIANE":
            if "NB_ITER" not in self.simu_params:
                if "FREQ_PENTE" in DSP_args:
                    Xt = gene_traj_gauss_evol_ND(generator, Data_cohe, **DSP_args)
                else:
                    Xt = DSP2ACCE_ND(DSP_args["FONC_DSP"], Data_cohe)

        if specmethode == "SPEC_FRACTILE":
            if "FREQ_PENTE" in DSP_args:
                alpha2 = RAND_VEC(DSP_args["MAT_COVC"], len(generator.sampler.liste_w2), para=2.0)
                DSP_args.update({"ALEA_DSP": alpha2})
                Xt = gene_traj_gauss_evol_ND(generator, Data_cohe, **DSP_args)
            else:
                fonc_dsp_rv = RAND_DSP(
                    DSP_args["MAT_COVC"], len(generator.sampler.liste_w2), DSP_args["FONC_DSP"]
                )
                Xt = DSP2ACCE_ND(fonc_dsp_rv, Data_cohe)
        return Xt

    def build_TimeHistories(self, generator):
        """build Time Histories for iterated median spec case"""
        DSP_args = generator.DSP_args
        if self.TYPE == "COEF_CORR":
            rho = self.DEFI_COHE["COEF_CORR"]
            dim = 2
            #        Mat_cor = NP.array([[1.0 , rho ],[rho ,1.0]])
            #             Mat_cor = CALC_CORRE(rho, dim)
            with disable_fpe():
                Mat_cohe = NP.linalg.cholesky(CALC_CORRE(rho, dim))
            Data_cohe = self.DEFI_COHE
            Data_cohe.update({"MAT_COHE": Mat_cohe, "DIM": dim})
        else:
            Data_cohe = self.DEFI_COHE

        if "FREQ_PENTE" in DSP_args:
            fonc_dsp_opt, liste_rv = itersimcortir_SRO(
                generator, DSP_args["FONC_DSP"], Data_cohe, self.nbtirage, **generator.SRO_args
            )
            vop, amo, R0, R2, f_FIT = DSP2FR(fonc_dsp_opt, DSP_args["FC"])
            DSP_args.update(
                {"FREQ_FOND": vop, "AMORT": amo, "para_R0": R0, "para_R2": R2, "fonc_FIT": f_FIT}
            )
            for rvtir in liste_rv:
                Xt = gene_traj_gauss_evol_ND(generator, Data_cohe, rv=rvtir, **DSP_args)
                nba = 1
                for acce in Xt:
                    accef = self.process_TimeHistory(generator, acce)
                    _f_out = DEFI_FONCTION(
                        ABSCISSE=tuple(generator.sampler.liste_temps),
                        ORDONNEE=tuple(accef),
                        **self.para_fonc_traj
                    )
                    if self.TYPE == "COEF_CORR":
                        nom_acce = "ACCE" + str(nba)
                        generator.tab.append(
                            {
                                "NUME_ORDRE": self.ntir + 1,
                                "FONCTION": _f_out.getName(),
                                "NOM_PARA": nom_acce,
                            }
                        )
                    else:
                        nom_no = self.liste_nom[nba - 1]
                        generator.tab.append(
                            {
                                "NUME_ORDRE": self.ntir + 1,
                                "FONCTION": _f_out.getName(),
                                "NOEUD": nom_no,
                            }
                        )
                    generator.tab.referenceToDataStructure.append(_f_out)
                    nba = nba + 1
                self.ntir = self.ntir + 1
        else:
            fonc_dsp_opt, liste_rv = itersimcortir_SRO(
                generator, DSP_args["FONC_DSP"], Data_cohe, self.nbtirage, **generator.SRO_args
            )
            for rvtir in liste_rv:
                Xt = DSP2ACCE_ND(fonc_dsp_opt, Data_cohe, rv=rvtir)
                nba = 1
                for acce in Xt:
                    accef = self.process_TimeHistory(generator, acce)
                    _f_out = DEFI_FONCTION(
                        ABSCISSE=tuple(generator.sampler.liste_temps),
                        ORDONNEE=tuple(accef),
                        **self.para_fonc_traj
                    )
                    generator.tab.referenceToDataStructure.append(_f_out)
                    if self.TYPE == "COEF_CORR":
                        nom_acce = "ACCE" + str(nba)
                        generator.tab.append(
                            {
                                "NUME_ORDRE": self.ntir + 1,
                                "FONCTION": _f_out.getName(),
                                "NOM_PARA": nom_acce,
                            }
                        )
                    else:
                        nom_no = self.liste_nom[nba - 1]
                        generator.tab.append(
                            {
                                "NUME_ORDRE": self.ntir + 1,
                                "FONCTION": _f_out.getName(),
                                "NOEUD": nom_no,
                            }
                        )
                    nba = nba + 1
                self.ntir = self.ntir + 1


class SimulatorSPECScalar(Simulator):
    """Construct scalar signal for SPEC class"""

    def run(self, generator):
        """Create the result table of functions"""
        if self.simu_params["SPEC_METHODE"] == "SPEC_MEDIANE" and "NB_ITER" in self.simu_params:
            self.build_TimeHistories(generator)
        else:
            for _ in range(self.nbtirage):
                Xt = self.build_TimeHistory(generator)
                Xt = self.process_TimeHistory(generator, NP.array(Xt))
                _f_out = DEFI_FONCTION(
                    ABSCISSE=tuple(generator.sampler.liste_temps),
                    ORDONNEE=tuple(Xt),
                    **self.para_fonc_traj
                )
                generator.tab.append({"NUME_ORDRE": self.ntir + 1, "FONCTION": _f_out.getName()})
                generator.tab.referenceToDataStructure.append(_f_out)
                self.ntir = self.ntir + 1

    def build_TimeHistory(self, generator):
        """build Time History for Spectrum class"""
        specmethode = self.simu_params["SPEC_METHODE"]
        DSP_args = generator.DSP_args
        if self.INFO == 2:
            UTMESS("I", "PROBA0_13", vali=self.ntir + 1)

        if specmethode == "SPEC_UNIQUE":
            if "NB_ITER" not in self.simu_params:
                if "FREQ_PENTE" in DSP_args:
                    Xt = gene_traj_gauss_evol1D(generator, **DSP_args)
                else:
                    Xt = DSP2ACCE1D(DSP_args["FONC_DSP"])
            else:  #'NB_ITER' in self.method_params
                if "FREQ_PENTE" in DSP_args:
                    fonc_dsp_opt, rv0 = itersim_SRO(
                        generator, DSP_args["FONC_DSP"], NB_TIRAGE=1, **generator.SRO_args
                    )
                    vop, amo, R0, R2, f_FIT = DSP2FR(fonc_dsp_opt, DSP_args["FC"])
                    DSP_args.update(
                        {
                            "FREQ_FOND": vop,
                            "AMORT": amo,
                            "para_R0": R0,
                            "para_R2": R2,
                            "fonc_FIT": f_FIT,
                        }
                    )
                    Xt = gene_traj_gauss_evol1D(generator, rv=rv0[0], **DSP_args)
                else:
                    fonc_dsp_opt, rv0 = itersim_SRO(
                        generator, DSP_args["FONC_DSP"], NB_TIRAGE=1, **generator.SRO_args
                    )
                    Xt = DSP2ACCE1D(fonc_dsp_opt, rv0[0])
        if specmethode == "SPEC_FRACTILE":
            if "FREQ_PENTE" in DSP_args:
                alpha2 = RAND_VEC(DSP_args["MAT_COVC"], len(generator.sampler.liste_w2), para=2.0)
                DSP_args.update({"ALEA_DSP": alpha2})
                Xt = gene_traj_gauss_evol1D(generator, **DSP_args)
            else:
                fonc_dsp_rv = RAND_DSP(
                    DSP_args["MAT_COVC"], len(generator.sampler.liste_w2), DSP_args["FONC_DSP"]
                )
                Xt = DSP2ACCE1D(fonc_dsp_rv)
        if specmethode == "SPEC_MEDIANE":
            if "NB_ITER" not in self.simu_params:
                if "FREQ_PENTE" in DSP_args:
                    Xt = gene_traj_gauss_evol1D(generator, **DSP_args)
                else:
                    Xt = DSP2ACCE1D(DSP_args["FONC_DSP"])
        return Xt

    def build_TimeHistories(self, generator):
        """build Time Histories for iterated median spec case"""
        DSP_args = generator.DSP_args
        if "FREQ_PENTE" in DSP_args:
            fonc_dsp_opt, liste_rv = itersim_SRO(
                generator, DSP_args["FONC_DSP"], NB_TIRAGE=self.nbtirage, **generator.SRO_args
            )
            vop, amo, R0, R2, f_FIT = DSP2FR(fonc_dsp_opt, DSP_args["FC"])
            DSP_args.update(
                {"FREQ_FOND": vop, "AMORT": amo, "para_R0": R0, "para_R2": R2, "fonc_FIT": f_FIT}
            )
            for rvtir in liste_rv:
                Xt = gene_traj_gauss_evol1D(generator, rv=rvtir, **DSP_args)
                Xt = self.process_TimeHistory(generator, Xt)
                _f_out = DEFI_FONCTION(
                    ABSCISSE=tuple(generator.sampler.liste_temps),
                    ORDONNEE=tuple(Xt),
                    **self.para_fonc_traj
                )
                generator.tab.referenceToDataStructure.append(_f_out)
                generator.tab.append({"NUME_ORDRE": self.ntir + 1, "FONCTION": _f_out.getName()})
                self.ntir = self.ntir + 1
        else:
            fonc_dsp_opt, liste_rv = itersim_SRO(
                generator, DSP_args["FONC_DSP"], NB_TIRAGE=self.nbtirage, **generator.SRO_args
            )
            for rvtir in liste_rv:
                Xt = DSP2ACCE1D(fonc_dsp_opt, rv=rvtir)
                Xt = self.process_TimeHistory(generator, Xt)
                _f_out = DEFI_FONCTION(
                    ABSCISSE=tuple(generator.sampler.liste_temps),
                    ORDONNEE=tuple(Xt),
                    **self.para_fonc_traj
                )
                generator.tab.referenceToDataStructure.append(_f_out)
                generator.tab.append({"NUME_ORDRE": self.ntir + 1, "FONCTION": _f_out.getName()})
                self.ntir = self.ntir + 1


class SimulatorSPECPhase(Simulator):
    """Construct vector valued signal with phase delay for SPEC class"""

    def run(self, generator):
        """build result for phase SPEC class"""
        self.liste_nom, l2 = get_group_nom_coord(
            self.DEFI_COHE["GROUP_NO_INTERF"], self.DEFI_COHE["MAILLAGE"]
        )
        self.DEFI_COHE.update({"NOEUDS_INTERF": l2})
        if self.DEFI_COHE["COOR_REFE"] is None:
            coord_ref = get_no_refe(self.DEFI_COHE)
            self.DEFI_COHE.update({"COOR_REFE": coord_ref})
            UTMESS("I", "SEISME_77", valr=(coord_ref[0], coord_ref[1], coord_ref[2]))

        if self.simu_params["SPEC_METHODE"] == "SPEC_MEDIANE" and "NB_ITER" in self.simu_params:
            self.build_TimeHistories(generator)
        #            raise NotImplementedError('must be implemented later on')
        else:
            for _ in range(self.nbtirage):
                Xt = self.build_TimeHistory(generator)
                nba = 1
                for accef in Xt:
                    _f_out = DEFI_FONCTION(
                        ABSCISSE=tuple(generator.sampler.liste_temps),
                        ORDONNEE=tuple(accef),
                        **self.para_fonc_traj
                    )
                    generator.tab.referenceToDataStructure.append(_f_out)
                    nom_no = self.liste_nom[nba - 1]
                    generator.tab.append(
                        {"NUME_ORDRE": self.ntir + 1, "FONCTION": _f_out.getName(), "NOEUD": nom_no}
                    )
                    nba = nba + 1
                self.ntir = self.ntir + 1

    def build_TimeHistory(self, generator):
        """build series of phase delayed Time History for Spectrum class"""
        specmethode = self.simu_params["SPEC_METHODE"]
        DSP_args = generator.DSP_args
        Data_phase = self.DEFI_COHE
        if self.INFO == 2:
            UTMESS("I", "PROBA0_13", vali=self.ntir + 1)

        if specmethode == "SPEC_UNIQUE":
            if "NB_ITER" not in self.simu_params:
                if "FREQ_PENTE" in DSP_args:
                    Xt = gene_traj_gauss_evol1D(generator, **DSP_args)
                else:
                    Xt = DSP2ACCE1D(DSP_args["FONC_DSP"])
            else:  #'NB_ITER' in self.method_params
                if "FREQ_PENTE" in DSP_args:
                    fonc_dsp_opt, rv0 = itersim_SRO(
                        generator, DSP_args["FONC_DSP"], NB_TIRAGE=1, **generator.SRO_args
                    )
                    vop, amo, R0, R2, f_FIT = DSP2FR(fonc_dsp_opt, DSP_args["FC"])
                    DSP_args.update(
                        {
                            "FREQ_FOND": vop,
                            "AMORT": amo,
                            "para_R0": R0,
                            "para_R2": R2,
                            "fonc_FIT": f_FIT,
                        }
                    )
                    Xt = gene_traj_gauss_evol1D(generator, rv=rv0[0], **DSP_args)
                else:
                    fonc_dsp_opt, rv0 = itersim_SRO(
                        generator, DSP_args["FONC_DSP"], NB_TIRAGE=1, **generator.SRO_args
                    )
                    Xt = DSP2ACCE1D(fonc_dsp_opt, rv0[0])
        if specmethode == "SPEC_FRACTILE":
            if "FREQ_PENTE" in DSP_args:
                alpha2 = RAND_VEC(DSP_args["MAT_COVC"], len(generator.sampler.liste_w2), para=2.0)
                DSP_args.update({"ALEA_DSP": alpha2})
                Xt = gene_traj_gauss_evol1D(generator, **DSP_args)
            else:
                fonc_dsp_rv = RAND_DSP(
                    DSP_args["MAT_COVC"], len(generator.sampler.liste_w2), DSP_args["FONC_DSP"]
                )
                Xt = DSP2ACCE1D(fonc_dsp_rv)
        if specmethode == "SPEC_MEDIANE":
            if "NB_ITER" not in self.simu_params:
                if "FREQ_PENTE" in DSP_args:
                    Xt = gene_traj_gauss_evol1D(generator, **DSP_args)
                else:
                    Xt = DSP2ACCE1D(DSP_args["FONC_DSP"])
        Xt = self.process_TimeHistory(generator, Xt)
        Xtl = calc_phase_delay(generator.sampler.liste_temps, Xt, Data_phase)
        return Xtl

    def build_TimeHistories(self, generator):
        """build Time Histories for iterated median spec case"""
        DSP_args = generator.DSP_args
        Data_phase = self.DEFI_COHE
        if "FREQ_PENTE" in DSP_args:
            fonc_dsp_opt, liste_rv = itersim_SRO(
                generator, DSP_args["FONC_DSP"], NB_TIRAGE=self.nbtirage, **generator.SRO_args
            )
            vop, amo, R0, R2, f_FIT = DSP2FR(fonc_dsp_opt, DSP_args["FC"])
            DSP_args.update(
                {"FREQ_FOND": vop, "AMORT": amo, "para_R0": R0, "para_R2": R2, "fonc_FIT": f_FIT}
            )
            for rvtir in liste_rv:
                Xt = gene_traj_gauss_evol1D(generator, rv=rvtir, **DSP_args)
                Xt = self.process_TimeHistory(generator, Xt)
                Xtl = calc_phase_delay(generator.sampler.liste_temps, Xt, Data_phase)
                nba = 1
                for accef in Xtl:
                    _f_out = DEFI_FONCTION(
                        ABSCISSE=tuple(generator.sampler.liste_temps),
                        ORDONNEE=tuple(accef),
                        **self.para_fonc_traj
                    )
                    generator.tab.referenceToDataStructure.append(_f_out)
                    nom_no = self.liste_nom[nba - 1]
                    generator.tab.append(
                        {"NUME_ORDRE": self.ntir + 1, "FONCTION": _f_out.getName(), "NOEUD": nom_no}
                    )
                    nba = nba + 1
                self.ntir = self.ntir + 1
        else:
            fonc_dsp_opt, liste_rv = itersim_SRO(
                generator, DSP_args["FONC_DSP"], NB_TIRAGE=self.nbtirage, **generator.SRO_args
            )
            for rvtir in liste_rv:
                Xt = DSP2ACCE1D(fonc_dsp_opt, rv=rvtir)
                Xt = self.process_TimeHistory(generator, Xt)
                Xtl = calc_phase_delay(generator.sampler.liste_temps, Xt, Data_phase)
                nba = 1
                for acce in Xtl:
                    _f_out = DEFI_FONCTION(
                        ABSCISSE=tuple(generator.sampler.liste_temps),
                        ORDONNEE=tuple(acce),
                        **self.para_fonc_traj
                    )
                    generator.tab.referenceToDataStructure.append(_f_out)
                    nom_no = self.liste_nom[nba - 1]
                    generator.tab.append(
                        {"NUME_ORDRE": self.ntir + 1, "FONCTION": _f_out.getName(), "NOEUD": nom_no}
                    )
                    nba = nba + 1
                self.ntir = self.ntir + 1


class SimulatorDSPPhase(Simulator):
    """Construct series of signals with phase delay for DSP class"""

    def run(self, generator):
        """build result for vector DSP class"""
        liste_nom, l2 = get_group_nom_coord(
            self.DEFI_COHE["GROUP_NO_INTERF"], self.DEFI_COHE["MAILLAGE"]
        )
        self.DEFI_COHE.update({"NOEUDS_INTERF": l2})
        if self.DEFI_COHE["COOR_REFE"] is None:
            coord_ref = get_no_refe(self.DEFI_COHE)
            UTMESS("I", "SEISME_77", valr=(coord_ref[0], coord_ref[1], coord_ref[2]))

        for _ in range(self.nbtirage):
            Xt = self.build_TimeHistory(generator)
            nba = 1
            for accef in Xt:
                _f_out = DEFI_FONCTION(
                    ABSCISSE=tuple(generator.sampler.liste_temps),
                    ORDONNEE=tuple(accef),
                    **self.para_fonc_traj
                )
                generator.tab.referenceToDataStructure.append(_f_out)
                nom_no = liste_nom[nba - 1]
                generator.tab.append(
                    {"NUME_ORDRE": self.ntir + 1, "FONCTION": _f_out.getName(), "NOEUD": nom_no}
                )
                nba = nba + 1
            self.ntir = self.ntir + 1

    def build_TimeHistory(self, generator):
        """build series of delayed Time History for DSP class"""
        Data_phase = self.DEFI_COHE
        if self.INFO == 2:
            UTMESS("I", "PROBA0_13", vali=self.ntir + 1)
        if "FREQ_PENTE" in generator.DSP_args:
            Xt = gene_traj_gauss_evol1D(generator, **generator.DSP_args)
        else:
            Xt = DSP2ACCE1D(generator.DSP_args["FONC_DSP"])
        Xt = self.process_TimeHistory(generator, Xt)
        Xtl = calc_phase_delay(generator.sampler.liste_temps, Xt, Data_phase)
        return Xtl
