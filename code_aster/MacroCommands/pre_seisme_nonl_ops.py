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

import pprint
import sys
import traceback

from libaster import AsterError

from ..Cata.Syntax import _F
from ..CodeCommands import (
    AFFE_CARA_ELEM,
    AFFE_CHAR_MECA,
    AFFE_MATERIAU,
    AFFE_MODELE,
    ASSE_MAILLAGE,
    ASSE_MATRICE,
    CALC_MATR_ELEM,
    CALC_MODES,
    COMB_MATR_ASSE,
    CREA_CHAMP,
    CREA_MAILLAGE,
    DEFI_BASE_MODALE,
    DEFI_GROUP,
    DEFI_INTERF_DYNA,
    DEFI_LIST_INST,
    DEFI_LIST_REEL,
    DEFI_MAILLAGE,
    DYNA_NON_LINE,
    LIRE_IMPE_MISS,
    MACR_ELEM_DYNA,
    MODE_STATIQUE,
    NUME_DDL,
    NUME_DDL_GENE,
    PROD_MATR_CHAM,
    STAT_NON_LINE,
)
from ..Helpers.LogicalUnit import LogicalUnitFile
from ..Messages import UTMESS
from ..Objects import Mesh as CA_Mesh
from ..Utilities import medcoupling as medc


def pre_seisme_nonl_ops(self, **args):
    """Corps de la macro PRE_SEISME_NONL"""
    # conteneur des paramètres du calcul
    param = Properties(**args)

    # création de l'objet PreSeismeNonL_xxx
    calcul = PreSeismeNonL.Factory(self, param)
    try:
        calcul.run()
    except AsterError:
        raise
    except Exception as err:
        trace = "".join(traceback.format_tb(sys.exc_info()[2]))
        UTMESS("F", "SUPERVIS2_3", valk=("PRE_SEISME_NONL", trace, str(err)))


class PreSeismeNonL:
    """Define a general methods for a PRE_SEISME_NONL calculation."""

    option_calcul = None

    @staticmethod
    def Factory(parent, param):
        """Factory that returns the calculation object"""
        if param["PRE_CALC_MISS"]:
            return PreCalcMiss(parent, param)
        elif param["POST_CALC_MISS"]:
            return PostCalcMiss(parent, param)
        elif param["STAT_DYNA"]:
            return StatDyna(parent, param)
        else:
            raise NotImplementedError("option calcul non défini")

    def __init__(self, parent, param):
        """initializations"""
        self.parent = parent
        self.param = param
        self.type_calcul = None
        self.typ_ISS = None
        self.type_IFS = None
        self.model = None
        self.mail = None
        self.mael = None
        self.bamo = None

    def run(self):
        """Execute calculation"""
        self.define_mesh()
        self.define_model()
        self.define_bamo()
        self.define_mael()

    def define_mesh(self):
        """Define the mesh"""
        self.mail = Mesh(self, self.param)

    def define_model(self):
        """Define the numerical model"""
        self.model = Model.factory(self, self.mail, self.param)
        self.model.DefineOut()

    def define_bamo(self):
        """Define the modal basis"""
        self.set_type()
        self.calc_base_modale()
        self.bamo.DefineOut()

    def define_mael(self):
        """Define the super-element"""
        self.mael = MacroElement(self, self.param, self.bamo)
        self.mael.DefineOut()

    def calc_base_modale(self):
        """Execute the eigenmodes calculation"""
        raise NotImplementedError("must be defined in a derivated class")

    def set_type(self):
        """Define different type of loads"""
        raise NotImplementedError("must be defined in a derivated class")

    def register_result(self, *args):
        """Define output depending on user choices"""
        self.parent.register_result(*args)


class PreCalcMiss(PreSeismeNonL):
    """Define the interface modal basis used for the soil impedance calculation (PRE_CALC_MISS)."""

    option_calcul = "PRE_CALC_MISS"

    def calc_base_modale(self):
        """Execute the eigenmodes calculation"""
        if self.type_calcul == "ISFS":
            nbmodes1 = self.param["PRE_CALC_MISS"]["NMAX_MODE_ISS"]
            nbmodes2 = self.param["PRE_CALC_MISS"]["NMAX_MODE_IFS"]
            bamoISS = BaseModale(self, self.param, self.model, nbmodes1, self.typ_ISS)
            bamoISFS = BaseModale(self, self.param, self.model, nbmodes2, self.typ_IFS)
            bamoISFS.defi_base_modale()
            bamoISFS.combine_base_modale(bamoISS)
            self.bamo = bamoISFS
        else:
            nbmodes1 = self.param["PRE_CALC_MISS"]["NMAX_MODE_ISS"]
            bamoISS = BaseModale(self, self.param, self.model, nbmodes1, self.typ_ISS)
            bamoISS.defi_base_modale()
            self.bamo = bamoISS

    def set_type(self):
        """Set the type of MISS calculation"""
        self.type_calcul = self.param["PRE_CALC_MISS"]["CALC_MISS_OPTION"]
        if self.param["PRE_CALC_MISS"]["REDUC_DYNA_ISS"] == "OUI":
            self.typ_ISS = "DYNA"
        else:
            self.typ_ISS = "STAT"
        if self.param["PRE_CALC_MISS"]["REDUC_DYNA_IFS"] == "OUI":
            self.typ_IFS = "DYNA"
        else:
            self.typ_IFS = "STAT"


class PostCalcMiss(PreSeismeNonL):
    """Define the numerical model necessary for a transient nonlinear calculation (POST_CALC_MISS)."""

    option_calcul = "POST_CALC_MISS"

    def define_mesh(self):
        """Factory that returns the Model object"""
        self.mael = MacroElement(self, self.param)
        self.mail = Mesh(self, self.param, self.mael)
        self.mail.build_mesh()
        self.mail.DefineOut()

    def define_mael(self):
        """Define the super-element"""

    def define_bamo(self):
        """Define the modal basis"""

    def calc_base_modale(self):
        """Execute the eigenmodes calculation"""

    def set_type(self):
        """Set the type of MISS calculation"""


class BaseModale:
    """Define a modal basis."""

    def __init__(self, parent, param, model=None, nb_modes=None, typ=None, bamo=None):
        """initializations"""
        self.parent = parent
        self.param = param
        self.model = model
        self.bamo = bamo
        self.cmd_bamo = None
        self.modes = None
        self.nb_modes = nb_modes
        self.interf_dyna = None
        self.nume_ddl = None
        if typ == "STAT":
            self.assemblage_stat()
            self.modes_statiques()
        else:
            self.assemblage_dyna()
            self.modes_dynamiques()
        self.grno_interf = self.from_grma_interf()

    def from_grma_interf(self):
        grma = self.param["PRE_CALC_MISS"]["GROUP_MA_INTERF"]
        self.model.mesh.add_group_no(grma)
        return grma

    def get_grno_interf(self):
        return self.grno_interf

    def get_num_modes(self):
        """Return the number of modes stored within self.modes"""
        return self.nb_modes

    def get_modes(self):
        """Return the calculated normal modes"""
        return self.modes

    def get_bamo(self):
        """Return the complete modal basis, which might contain other complementary modes"""
        return self.bamo

    def get_mass(self):
        """Return the FE mass matrix associated to the model self.model"""
        return self.matr_mass

    def get_rigi(self):
        """Return the FE stiffness matrix associated to the model self.model"""
        return self.matr_rigi

    def check_radier_rigide(self):
        """Check if solid rigid conditions are considered for the interface"""
        charge = self.param["AFFE_CHAR_MECA"]
        for key in charge:
            if key == "LIAISON_SOLIDE":
                msg_error = "\n\nLe mot-clé GROUP_NO_CENT est obligatoire lorsqu'une LIAISON_SOLIDE est définie"
                assert "GROUP_NO_CENT" in self.param["PRE_CALC_MISS"], msg_error
                return True
        return False

    def defi_interf_dyna(self):
        """Build the dynamic interface where the substructuring approach is applied"""
        if (
            self.param["PRE_CALC_MISS"]["REDUC_DYNA_ISS"] == "OUI"
            or self.param["PRE_CALC_MISS"]["REDUC_DYNA_IFS"] == "OUI"
        ):
            grno = self.get_grno_interf()
        else:
            if self.check_radier_rigide():
                grno = self.param["PRE_CALC_MISS"]["GROUP_NO_CENT"]
            else:
                grno = self.get_grno_interf()
        _C_LIM0 = AFFE_CHAR_MECA(
            MODELE=self.model.get_model(), DDL_IMPO=_F(GROUP_NO=grno, DX=0.0, DY=0.0, DZ=0.0)
        )
        _NUME0 = NUME_DDL(MODELE=self.model.get_model(), CHARGE=_C_LIM0)
        _INTERDY = DEFI_INTERF_DYNA(
            NUME_DDL=_NUME0, INTERFACE=_F(NOM="DROITE", TYPE="CRAIGB", GROUP_NO=grno)
        )
        self.interf_dyna = _INTERDY

    def assemblage_dyna(self):
        """Assemble FE matrices associated to self.model with non-static boundary conditions"""
        _MATMele = CALC_MATR_ELEM(
            MODELE=self.model.get_model(),
            CHAM_MATER=self.model.get_mate(),
            CARA_ELEM=self.model.get_cara_elem(),
            OPTION="MASS_MECA",
        )
        _MATRele = CALC_MATR_ELEM(
            MODELE=self.model.get_model(),
            CHAM_MATER=self.model.get_mate(),
            CARA_ELEM=self.model.get_cara_elem(),
            OPTION="RIGI_MECA",
        )
        _NUME = NUME_DDL(MATR_RIGI=_MATRele)
        _MATRIGI = ASSE_MATRICE(MATR_ELEM=_MATRele, NUME_DDL=_NUME)
        _MATMASS = ASSE_MATRICE(MATR_ELEM=_MATMele, NUME_DDL=_NUME)
        self.matr_rigi = _MATRIGI
        self.matr_mass = _MATMASS
        self.nume_ddl = _NUME

    def assemblage_stat(self):
        """Assemble FE matrices associated to self.model with static boundary conditions"""
        _MATMele = CALC_MATR_ELEM(
            MODELE=self.model.get_model(),
            CHAM_MATER=self.model.get_mate(),
            CARA_ELEM=self.model.get_cara_elem(),
            CHARGE=self.model.get_cond_lim(),
            OPTION="MASS_MECA",
        )
        _MATRele = CALC_MATR_ELEM(
            MODELE=self.model.get_model(),
            CHAM_MATER=self.model.get_mate(),
            CARA_ELEM=self.model.get_cara_elem(),
            CHARGE=self.model.get_cond_lim(),
            OPTION="RIGI_MECA",
        )
        _NUME = NUME_DDL(MATR_RIGI=_MATRele)
        _MATRIGI = ASSE_MATRICE(MATR_ELEM=_MATRele, NUME_DDL=_NUME)
        _MATMASS = ASSE_MATRICE(MATR_ELEM=_MATMele, NUME_DDL=_NUME)
        self.matr_rigi = _MATRIGI
        self.matr_mass = _MATMASS
        self.nume_ddl = _NUME

    def modes_statiques(self):
        """Calculate static/interface modes"""
        if self.check_radier_rigide():
            grno = self.param["PRE_CALC_MISS"]["GROUP_NO_CENT"]
        else:
            grno = self.get_grno_interf()
        _modsta = MODE_STATIQUE(
            MATR_RIGI=self.matr_rigi, MODE_STAT=_F(GROUP_NO=grno, TOUT_CMP="OUI")
        )
        self.modes = _modsta

    def modes_dynamiques(self):
        """Calculate dynamic modes or modal shapes"""
        _modyna = CALC_MODES(
            OPTION="PLUS_PETITE",
            MATR_RIGI=self.matr_rigi,
            MATR_MASS=self.matr_mass,
            CALC_FREQ=_F(NMAX_FREQ=self.nb_modes),
        )
        self.modes = _modyna

    def defi_base_modale(self):
        """Build the complete modal basis"""
        self.defi_interf_dyna()
        cmd_ritz_1 = _F(MODE_MECA=self.modes, NMAX_MODE=0)
        cmd_ritz_2 = _F(MODE_INTF=self.modes, NMAX_MODE=self.nb_modes)
        self.cmd_bamo = {
            "RITZ": (cmd_ritz_1, cmd_ritz_2),
            "INTERF_DYNA": self.interf_dyna,
            "NUME_REF": self.nume_ddl,
        }

    def combine_base_modale(self, bamo):
        """Combine two different modal basis"""
        _BAMO0 = DEFI_BASE_MODALE(**self.cmd_bamo)
        cmd_ritz_1 = _F(BASE_MODALE=_BAMO0)
        cmd_ritz_2 = _F(MODE_INTF=bamo.get_modes(), NMAX_MODE=bamo.get_num_modes())
        self.cmd_bamo = {
            "RITZ": (cmd_ritz_1, cmd_ritz_2),
            "INTERF_DYNA": self.interf_dyna,
            "NUME_REF": self.nume_ddl,
        }

    def DefineOut(self):
        """Define output depending on user choices"""
        if "BASE_MODALE" in self.param["RESULTAT"]:
            _BAMO = DEFI_BASE_MODALE(**self.cmd_bamo)
            self.parent.register_result(_BAMO, self.param["RESULTAT"]["BASE_MODALE"])
            self.bamo = _BAMO


class MacroElement:
    """Define a sub-structure, also known as super-element or macro-element."""

    def __init__(self, parent, param, BaMo=None):
        """initializations"""
        self.parent = parent
        self.param = param
        self.bamo = BaMo
        self.set_mael()

    def check_reduc_dyna(self):
        """Check if dynamic reduction is used within the super-element"""
        reduc_dyna = False

        if self.param["PRE_CALC_MISS"]:
            if (
                self.param["PRE_CALC_MISS"]["REDUC_DYNA_ISS"] == "OUI"
                or self.param["PRE_CALC_MISS"]["REDUC_DYNA_IFS"] == "OUI"
            ):
                reduc_dyna = True
        if self.param["POST_CALC_MISS"]:
            noeud_cmp = self.mael.getMechanicalMode().getAccessParameters()["NOEUD_CMP"]
            if noeud_cmp[0] is None:
                reduc_dyna = True
        return reduc_dyna

    def set_mael(self):
        """Build the super-element"""
        if self.param["POST_CALC_MISS"]:
            self.mael = self.param["POST_CALC_MISS"]["MACR_ELEM_DYNA"]
        else:
            self.mael = None

    def get_mael(self):
        """Return the super-element"""
        return self.mael

    def DefineOut(self):
        """Define output depending on user choices"""
        Bamo = self.bamo.get_bamo()
        _Mael = MACR_ELEM_DYNA(
            BASE_MODALE=Bamo,
            MATR_RIGI=self.bamo.get_rigi(),
            MATR_MASS=self.bamo.get_mass(),
            SANS_GROUP_NO=self.bamo.get_grno_interf(),
        )
        self.parent.register_result(_Mael, self.param["RESULTAT"]["MACR_ELEM_DYNA"])


class Properties:
    """Define a dictionary containing the keywords of the model properties."""

    def __init__(self, **kwargs):
        """initializations"""
        self._keywords = {}
        for key in list(kwargs.keys()):
            if hasattr(kwargs[key], "List_F"):
                self._keywords[key] = kwargs[key].List_F()[0]
            else:
                self._keywords[key] = kwargs[key]

    def __getitem__(self, key):
        return self._keywords.get(key)

    def get_nested_key(self, path):
        """Get a value from a nested dictionary"""
        diction = self._keywords
        for key in path:
            if key not in list(diction.keys()):
                diction[key] = {}
            diction = diction[key]
        return diction

    def set_key(self, path, value):
        """Set a value within nested dictionary"""
        self.get_nested_key(path[:-1])[path[-1]] = value

    def add_MCFACT(self, path, addedKey):
        """Add new Mot-Clé Facteur into a nested dictionary"""
        diction = self._keywords[path[0]]
        if path[1] not in list(diction.keys()):
            self.set_key(path, [])
        self._keywords[path[0]][path[1]] = list(self._keywords[path[0]][path[1]])
        self._keywords[path[0]][path[1]].insert(0, addedKey)

    def copy(self):
        """Return a new copy of the properties"""
        return self.__class__(**self._keywords.copy())

    def __repr__(self):
        """Print formatted dictionary of the attribut self._keywords"""
        return pprint.pformat(self._keywords)


class Model:
    """Define a numerical model."""

    option_calcul = None

    @staticmethod
    def factory(parent, mail, properties):
        """Factory that returns the Model object"""
        if properties["PRE_CALC_MISS"]:
            if mail.check_ficti_nodes():
                return ModelBaMoReduc(parent, mail, properties)
            else:
                return ModelBaseModale(parent, mail, properties)
        elif mail.check_ficti_nodes():
            return ModelDynaReduc(parent, mail, properties)
        else:
            return ModelMacrElem(parent, mail, properties)

    def __init__(self, parent, mail, properties):
        """initializations"""
        if not self.option_calcul:
            raise NotImplementedError("option_calcul non défini")
        self.parent = parent
        self.args = properties.copy()
        self.mesh = mail
        self.modele = None
        self.mate = None
        self.cara_elem = None
        self.cond_lim = None
        self.affe_modele()
        self.affe_materiau()
        self.affe_cara_elem()
        self.affe_char_meca()

    def get_properties(self):
        """Return the set of keywords related to the properties of the numerical model"""
        return self.args

    def modi_mesh(self, msh):
        """Modify the existing mesh of the current numerical model"""
        self.mesh = msh

    def get_model(self):
        """Return the numerical model"""
        return self.modele

    def get_mate(self):
        """Return the set of applied constitutive laws"""
        return self.mate

    def get_cara_elem(self):
        """Return the set of discrete properties"""
        return self.cara_elem

    def get_cond_lim(self):
        """Return the set of boundary conditions"""
        return self.cond_lim

    def affe_modele(self):
        """Define the physics of the numerical modelling"""
        raise NotImplementedError("must be defined in a derivated class")

    def affe_materiau(self):
        """Define constitutive discret properties"""
        raise NotImplementedError("must be defined in a derivated class")

    def affe_cara_elem(self):
        """Define discret properties"""
        raise NotImplementedError("must be defined in a derivated class")

    def affe_char_meca(self):
        """Define different type of loads"""
        raise NotImplementedError("must be defined in a derivated class")

    def DefineOut(self):
        """Define output depending on user choices"""

        _AMa = AFFE_MATERIAU(**self.args["AFFE_MATERIAU"])
        self.mate = _AMa
        if "CHAM_MATER" in self.args["RESULTAT"]:
            self.parent.register_result(_AMa, self.args["RESULTAT"]["CHAM_MATER"])

        _ACa = AFFE_CARA_ELEM(**self.args["AFFE_CARA_ELEM"])
        self.cara_elem = _ACa
        if "CARA_ELEM" in self.args["RESULTAT"]:
            self.parent.register_result(_ACa, self.args["RESULTAT"]["CARA_ELEM"])

        _ACh_CL = AFFE_CHAR_MECA(**self.args["AFFE_CHAR_MECA"])
        self.cond_lim = _ACh_CL
        if "CHARGE" in self.args["RESULTAT"]:
            for mcfact in self.args["RESULTAT"]["CHARGE"]:
                if mcfact["OPTION"] == "LAPL_TEMPS":
                    ACh_LT = AFFE_CHAR_MECA(**self.args["LAPL_TEMPS"])
                    self.parent.register_result(ACh_LT, mcfact["NOM"])
                elif mcfact["OPTION"] == "COND_LIM":
                    self.parent.register_result(_ACh_CL, mcfact["NOM"])


class StatDyna:
    def __init__(self, parent, properties):
        """initializations"""
        self.parent = parent
        self.args = properties
        self.resu_snl = properties["STAT_DYNA"]["RESULTAT"]

        self.modele = self.set_from_resu("model", self.resu_snl)
        self.maillage = self.set_from_resu("mesh", self.resu_snl)
        self.mater = self.set_from_resu("mater", self.resu_snl)
        self.cara_elem = self.set_from_resu("caraele", self.resu_snl)
        self.coef_amor = properties["STAT_DYNA"]["COEF_AMOR"]

        self.charges = list(properties["STAT_DYNA"]["EXCIT"])
        self.chsol = properties["STAT_DYNA"]["FORCE_SOL"]

        if properties["STAT_DYNA"]["COMPORTEMENT"]:
            self.comportement = properties["STAT_DYNA"]["COMPORTEMENT"]
        else:
            self.comportement = None

        if properties["STAT_DYNA"]["CONVERGENCE"]:
            self.converge = properties["STAT_DYNA"]["CONVERGENCE"]
        else:
            self.converge = None

        self.base_modale = properties["STAT_DYNA"]["BASE_MODALE"]
        self.UL_impe_freq = properties["STAT_DYNA"]["UNITE_IMPE_FREQ"]

        uniteImpeTemps = properties["STAT_DYNA"]["UNITE_IMPE_TEMPS"][0]
        if uniteImpeTemps["UNITE_RESU_RIGI"]:
            self.UL_impe_temps_K = uniteImpeTemps["UNITE_RESU_RIGI"]
            filename = LogicalUnitFile.filename_from_unit(self.UL_impe_temps_K)
            fid = open(filename, "r")
        if uniteImpeTemps["UNITE_RESU_AMOR"]:
            self.UL_impe_temps_C = uniteImpeTemps["UNITE_RESU_AMOR"]
            filename = LogicalUnitFile.filename_from_unit(self.UL_impe_temps_C)
            fid = open(filename, "r")
        if uniteImpeTemps["UNITE_RESU_MASS"]:
            self.UL_impe_temps_M = uniteImpeTemps["UNITE_RESU_MASS"]
            filename = LogicalUnitFile.filename_from_unit(self.UL_impe_temps_M)
            fid = open(filename, "r")

        data = fid.readline().split()
        fid.close()
        self.pas_inst_impe = float(data[1])
        self.nb_inst = properties["STAT_DYNA"]["NB_INST"]
        self.inst_init = self.resu_snl.LIST_PARA()["INST"][-1]

    def run(self):
        """Execute static-dynamic transition"""
        self.etapeStatique()
        self.etapeDynamique()

    def set_from_resu(self, what, resu):
        """Extract a parameter from a result"""
        assert what in ("mesh", "model", "caraele", "mater")
        if what == "mesh":
            return resu.getMesh()
        elif what == "model":
            return resu.getModel()
        elif what == "caraele":
            return resu.getElementaryCharacteristics()
        elif what == "mater":
            return resu.getMaterialField()

    def etapeStatique(self):
        """Execute static calculation"""
        _chsol = self.calc_chsol_equi()
        self.add_charge(_chsol)

        _lreel = DEFI_LIST_REEL(VALE=self.resu_snl.LIST_PARA()["INST"])
        _linst = DEFI_LIST_INST(METHODE="AUTO", DEFI_LIST=_F(LIST_INST=_lreel))

        _ResuSNL = STAT_NON_LINE(
            **self.non_line(EXCIT=self.charges, INCREMENT=_F(LIST_INST=_linst))
        )

        self.resu_snl = _ResuSNL

    def etapeDynamique(self):
        """Execute dynamic calculation"""
        _DEPF = CREA_CHAMP(
            TYPE_CHAM="NOEU_DEPL_R",
            OPERATION="EXTR",
            RESULTAT=self.resu_snl,
            NOM_CHAM="DEPL",
            INST=self.inst_init,
            INFO=1,
        )

        _SIGF = CREA_CHAMP(
            TYPE_CHAM="ELGA_SIEF_R",
            OPERATION="EXTR",
            RESULTAT=self.resu_snl,
            NOM_CHAM="SIEF_ELGA",
            INST=self.inst_init,
            INFO=1,
        )

        _VARF = CREA_CHAMP(
            TYPE_CHAM="ELGA_VARI_R",
            OPERATION="EXTR",
            RESULTAT=self.resu_snl,
            NOM_CHAM="VARI_ELGA",
            INST=self.inst_init,
            INFO=1,
        )

        _CNUL = CREA_CHAMP(
            TYPE_CHAM="NOEU_DEPL_R",
            OPERATION="ASSE",
            CHAM_NO=_DEPF,
            MODELE=self.modele,
            ASSE=(_F(TOUT="OUI", CHAM_GD=_DEPF, CUMUL="OUI", COEF_R=0.0),),
        )

        TFIN_TOTAL = self.pas_inst_impe * (self.nb_inst + 1)

        _larch = DEFI_LIST_REEL(
            DEBUT=self.inst_init, INTERVALLE=_F(JUSQU_A=TFIN_TOTAL, PAS=self.pas_inst_impe)
        )

        _linst2 = DEFI_LIST_INST(
            DEFI_LIST=_F(LIST_INST=_larch), ECHEC=_F(SUBD_PAS=2, SUBD_NIVEAU=10)
        )

        N_stab1 = int(0.75 * self.nb_inst)
        TFIN1 = N_stab1 * self.pas_inst_impe

        alpha_HHT = -7.0
        _ResuDNL = DYNA_NON_LINE(
            **self.non_line(
                SCHEMA_TEMPS=_F(
                    SCHEMA="HHT", ALPHA=alpha_HHT, MODI_EQUI="NON", FORMULATION="DEPLACEMENT"
                ),
                ETAT_INIT=_F(SIGM=_SIGF, VARI=_VARF, DEPL=_DEPF, VITE=_CNUL, ACCE=_CNUL),
                EXCIT=self.charges,
                INCREMENT=_F(LIST_INST=_linst2, INST_FIN=TFIN1),
                ARCHIVAGE=_F(LIST_INST=_larch),
            )
        )

        N_stab2 = self.nb_inst - N_stab1
        TFIN2 = TFIN1 + N_stab2 * self.pas_inst_impe
        TDEBUT2 = TFIN1

        self.modi_charge(self.chsol)
        _ResuDNL = DYNA_NON_LINE(
            **self.non_line(
                reuse=_ResuDNL,
                EXCIT=self.charges,
                SCHEMA_TEMPS=_F(
                    SCHEMA="HHT", ALPHA=alpha_HHT, MODI_EQUI="NON", FORMULATION="DEPLACEMENT"
                ),
                ETAT_INIT=_F(EVOL_NOLI=_ResuDNL, INST=TDEBUT2),
                INCREMENT=_F(LIST_INST=_linst2, INST_FIN=TFIN2),
                ARCHIVAGE=_F(LIST_INST=_larch),
            )
        )
        self.parent.register_result(_ResuDNL, self.args["RESULTAT"]["RESULTAT"])

    def non_line(self, **kwds):
        """Return the common keywords for STAT_NON_LINE and DYNA_NON_LINE
        All keywords can be overridden using `kwds`."""
        keywords = {
            "MODELE": self.modele,
            "CARA_ELEM": self.cara_elem,
            "CHAM_MATER": self.mater,
            "COMPORTEMENT": self.comportement,
            "CONVERGENCE": self.converge,
            "NEWTON": _F(
                MATRICE="TANGENTE",
                # PREDICTION='ELASTIQUE',
                REAC_ITER=1,
            ),
            "SOLVEUR": _F(METHODE="MUMPS", NPREC=20),
            "AFFICHAGE": _F(INFO_RESIDU="OUI"),
        }
        keywords.update(kwds)
        return keywords

    def calc_chsol_equi(self):
        """Compute the corrective force coming from the soil"""

        # compute impedence stiffness
        _NUMGEN = NUME_DDL_GENE(BASE=self.base_modale, STOCKAGE="PLEIN")

        _impeF = LIRE_IMPE_MISS(
            BASE=self.base_modale,
            UNITE_RESU_IMPE=self.UL_impe_freq,
            NUME_DDL_GENE=_NUMGEN,
            SYME="OUI",
            TYPE="ASCII",
            FREQ_EXTR=0.1,
        )

        _Ks = COMB_MATR_ASSE(COMB_C=(_F(MATR_ASSE=_impeF, COEF_C=1.0 + 0j),), SANS_CMP="LAGR")

        _Z0 = LIRE_IMPE_MISS(
            UNITE_RESU_IMPE=self.UL_impe_temps_K,
            SYME="OUI",
            INST_EXTR=0.0,
            BASE=self.base_modale,
            NUME_DDL_GENE=_NUMGEN,
        )

        _DIFFK = COMB_MATR_ASSE(
            COMB_C=(
                _F(MATR_ASSE=_Z0, COEF_C=1.0 + 0j),
                _F(MATR_ASSE=_Ks, COEF_C=-1.0 + 0j),  # Z0-KS
            ),
            SANS_CMP="LAGR",
        )

        # backup impedence stiffness
        mael = self.maillage.getDynamicMacroElements()[0]
        saved_matr_impe = mael.getImpedanceStiffnessMatrix()

        # set new impedence stiffness
        mael = MACR_ELEM_DYNA(
            reuse=mael, MACR_ELEM_DYNA=mael, BASE_MODALE=self.base_modale, MATR_IMPE_RIGI=_DIFFK
        )

        # compute stiffness with new impedence stiffness
        _rigiEle = CALC_MATR_ELEM(
            MODELE=self.modele,
            OPTION="RIGI_MECA",
            CALC_ELEM_MODELE="NON",
            CHAM_MATER=self.mater,
            CARA_ELEM=self.cara_elem,
            CHARGE=[elem["CHARGE"] for elem in self.charges],
        )

        _NUME = NUME_DDL(MATR_RIGI=_rigiEle)
        _MATKZ = ASSE_MATRICE(MATR_ELEM=_rigiEle, NUME_DDL=_NUME)

        # resore impedence stiffness
        mael = MACR_ELEM_DYNA(
            reuse=mael,
            MACR_ELEM_DYNA=mael,
            BASE_MODALE=self.base_modale,
            MATR_IMPE_RIGI=saved_matr_impe,
        )

        # compute load
        _DEPL0 = CREA_CHAMP(
            TYPE_CHAM="NOEU_DEPL_R",
            OPERATION="EXTR",
            RESULTAT=self.resu_snl,
            NOM_CHAM="DEPL",
            INST=self.inst_init,
            INFO=1,
        )
        _VFORC = PROD_MATR_CHAM(MATR_ASSE=_MATKZ, CHAM_NO=_DEPL0)
        _CHFORC = AFFE_CHAR_MECA(MODELE=self.modele, VECT_ASSE=_VFORC)

        return _CHFORC

    def add_charge(self, chsol):
        """Add the corrective force coming from the soil"""
        mcfact_chsol = _F(CHARGE=chsol)
        self.charges.append(mcfact_chsol)

    def modi_charge(self, chsol):
        """Replace the corrective force with the Laplace-Time force"""
        mcfact_chsol = _F(CHARGE=chsol)
        self.charges.pop()
        self.charges.append(mcfact_chsol)


class ModelMacrElem(Model):
    """Define a numerical model combined with superelements."""

    option_calcul = "Macro_Element"

    def affe_modele(self):
        """Define the physics of the numerical modelling"""
        self.link_macro_elem()
        self.fiction_model()
        self.parasol_model()
        self.args.set_key(("AFFE_MODELE", "MAILLAGE"), self.mesh.get_new_mesh())
        if "MODELE" in self.args["RESULTAT"]:
            _Modele = AFFE_MODELE(**self.args["AFFE_MODELE"])
            self.parent.register_result(_Modele, self.args["RESULTAT"]["MODELE"])
            self.modele = _Modele

    def fiction_model(self):
        """Define fictitious DoF's in the numerical modelling"""

    def fiction_cara_elem(self):
        """Define discret properties for the fictitious DoF's"""

    def link_macro_elem(self):
        """Link the supermesh to a physical phenomena or modelling"""
        super_ma = self.mesh.get_supermaille()
        mcfact_SupMa = _F(SUPER_MAILLE=super_ma, PHENOMENE="MECANIQUE")
        self.args.add_MCFACT(("AFFE_MODELE", "AFFE_SOUS_STRUC"), mcfact_SupMa)
        self.args.set_key(("AFFE_MODELE", "MAILLAGE"), self.mesh.get_new_mesh())

    def affe_materiau(self):
        """Define the constitutive law of the numerical modelling"""
        self.args.set_key(("AFFE_MATERIAU", "MAILLAGE"), self.mesh.get_new_mesh())

    def affe_cara_elem(self):
        """Define the constitutive law of discret elements"""
        self.fiction_cara_elem()
        self.parasol_cara_elem()
        self.args.set_key(("AFFE_CARA_ELEM", "MODELE"), self.modele)

    def affe_char_meca(self):
        """Define boundary conditions and external loads"""
        self.bound_conds()
        self.seismic_loads()
        self.other_loads()

    def bound_conds(self):
        """Define boundary conditions"""
        self.args.set_key(("AFFE_CHAR_MECA", "MODELE"), self.modele)

    def seismic_loads(self):
        """Define seismic loads"""
        mcfact = self.args.get_nested_key(("RESULTAT", "CHARGE"))
        for mm in mcfact:
            if mm["OPTION"] == "LAPL_TEMPS":
                cmd_charge = {"SUPER_MAILLE": "STAT1"}
                if "UNITE_RESU_RIGI" in self.args["POST_CALC_MISS"]:
                    UL_rigi = self.args["POST_CALC_MISS"]["UNITE_RESU_RIGI"]
                    cmd_charge["UNITE_RESU_RIGI"] = UL_rigi
                if "UNITE_RESU_MASS" in self.args["POST_CALC_MISS"]:
                    UL_mass = self.args["POST_CALC_MISS"]["UNITE_RESU_MASS"]
                    cmd_charge["UNITE_RESU_MASS"] = UL_mass
                if "UNITE_RESU_AMOR" in self.args["POST_CALC_MISS"]:
                    UL_amor = self.args["POST_CALC_MISS"]["UNITE_RESU_AMOR"]
                    cmd_charge["UNITE_RESU_AMOR"] = UL_amor
                charge_sol = _F(**cmd_charge)
                self.args.set_key(("LAPL_TEMPS", "MODELE"), self.modele)
                self.args.set_key(("LAPL_TEMPS", "FORCE_SOL"), charge_sol)

    def other_loads(self):
        """Define other loading entries"""

    def parasol_model(self):
        """Define the RIGI_PARASOL group within the model"""
        mcfact_DisTR = _F(GROUP_MA="PARA_SOL", PHENOMENE="MECANIQUE", MODELISATION="DIS_TR")
        self.args.add_MCFACT(("AFFE_MODELE", "AFFE"), mcfact_DisTR)

    def parasol_cara_elem(self):
        """Define the RIGI_PARASOL values of damping"""
        valC = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        mcfact_CEle = _F(
            GROUP_MA=self.args["POST_CALC_MISS"]["GROUP_MA_INTERF"],
            GROUP_MA_POI1="PARA_SOL",
            GROUP_NO_CENTRE=self.args["POST_CALC_MISS"]["GROUP_NO_CENT"],
            COEF_GROUP=1.0,
            CARA="A_TR_D_N",
            VALE=valC,
        )
        self.args.add_MCFACT(("AFFE_CARA_ELEM", "RIGI_PARASOL"), mcfact_CEle)
        valK = (1.0e-9, 1.0e-9, 1.0e-9, 1.0e-9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        valM = (1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3)
        mcfact_MEle = _F(GROUP_MA="PARA_SOL", CARA="M_TR_D_N", VALE=valK)
        mcfact_KEle = _F(GROUP_MA="PARA_SOL", CARA="K_TR_D_N", VALE=valM)
        self.args.add_MCFACT(("AFFE_CARA_ELEM", "DISCRET"), mcfact_MEle)
        self.args.add_MCFACT(("AFFE_CARA_ELEM", "DISCRET"), mcfact_KEle)


class ModelDynaReduc(ModelMacrElem):
    """Define a numerical model that uses dynamic reduction with superelements."""

    option_calcul = "Reduction_Dynamique"

    def fiction_model(self):
        """Define fictitious DoF's in the numerical modelling"""
        ma_fict_lst = []
        ma_fict_lst.append("MFICTIF")
        self.ma_fict = tuple(ma_fict_lst)
        mcfact_DisTR = _F(GROUP_MA=self.ma_fict, PHENOMENE="MECANIQUE", MODELISATION="DIS_TR")
        self.args.add_MCFACT(("AFFE_MODELE", "AFFE"), mcfact_DisTR)

    def fiction_cara_elem(self):
        """Define discret properties for the fictitious DoF's"""
        valK = (1.0e-9, 1.0e-9, 1.0e-9, 1.0e-9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        valM = (1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3)
        mcfact_MEle = _F(GROUP_MA=self.ma_fict, CARA="M_TR_D_N", VALE=valK)
        mcfact_KEle = _F(GROUP_MA=self.ma_fict, CARA="K_TR_D_N", VALE=valM)
        self.args.add_MCFACT(("AFFE_CARA_ELEM", "DISCRET"), mcfact_MEle)
        self.args.add_MCFACT(("AFFE_CARA_ELEM", "DISCRET"), mcfact_KEle)

    def other_loads(self):
        """Define the relation between the physical and generalized DoF's"""
        liaison_interf = _F(MACR_ELEM_DYNA=self.args["POST_CALC_MISS"]["MACR_ELEM_DYNA"])
        self.args.add_MCFACT(("AFFE_CHAR_MECA", "LIAISON_INTERF"), liaison_interf)


class ModelBaMoReduc(ModelDynaReduc):
    """Define a numerical model that uses dynamic reduction without superelements."""

    option_calcul = "Base_Modale"

    def link_macro_elem(self):
        """Assemble the mesh and the superelement"""

    def seismic_loads(self):
        """Define seismic loads"""

    def other_loads(self):
        """Define other loading entries"""

    def parasol_cara_elem(self):
        """D"""

    def parasol_model(self):
        """D"""


class ModelBaseModale(ModelMacrElem):
    """Define a numerical model that uses static reduction without superelements."""

    option_calcul = "Base_Modale"

    def link_macro_elem(self):
        """Assemble the mesh and the superelement"""

    def seismic_loads(self):
        """Define seismic loads"""

    def other_loads(self):
        """Define other loading entries"""

    def parasol_cara_elem(self):
        """D"""

    def parasol_model(self):
        """D"""


class Mesh:
    """Define the mesh of the numerical model."""

    def __init__(self, parent, param, mael=None):
        """initializations"""
        self.param = param
        self.parent = parent
        self.macro_elem = mael
        self.supermaille = []
        self.__builded = False
        self.new_mesh = param["AFFE_MODELE"]["MAILLAGE"]
        self.old_mesh = param["AFFE_MODELE"]["MAILLAGE"]
        self.add_fiction_mesh()

    def get_supermaille(self):
        """Return a set of supermeshes"""
        if self.__builded:
            return self.supermaille

    def build_mesh(self):
        """Assemble the mesh and the superelement"""
        self.old_mesh = self.new_mesh
        if self.macro_elem:
            list_SuperMa = self.__set_list_supermaille(self.macro_elem)
            _MAYADYN = DEFI_MAILLAGE(
                DEFI_SUPER_MAILLE=list_SuperMa,
                RECO_GLOBAL=_F(TOUT="OUI"),
                DEFI_NOEUD=_F(TOUT="OUI", INDEX=(1, 0, 1, 8)),
            )
            _MeshTmp = CREA_MAILLAGE(
                MAILLAGE=self.old_mesh,
                CREA_POI1=_F(
                    NOM_GROUP_MA="PARA_SOL",
                    GROUP_MA=self.param["POST_CALC_MISS"]["GROUP_MA_INTERF"],
                ),
            )
            self.new_mesh = ASSE_MAILLAGE(
                MAILLAGE_1=_MeshTmp, MAILLAGE_2=_MAYADYN, OPERATION="SOUS_STR"
            )
        self.__builded = True

    def DefineOut(self):
        """Define output depending on user choices"""
        if self.__builded:
            if "MAILLAGE" in self.param["RESULTAT"]:
                self.parent.register_result(self.new_mesh, self.param["RESULTAT"]["MAILLAGE"])

    def get_new_mesh(self):
        """Return the mesh concept which might contain a superelement"""
        return self.new_mesh

    def add_group_no(self, nom):
        lgrno = self.old_mesh.LIST_GROUP_NO()
        nomma = nom[0]
        check = 0
        for grp in lgrno:
            if grp[0] == nomma:
                check = 1
        if check == 0:
            DEFI_GROUP(
                reuse=self.old_mesh, MAILLAGE=self.old_mesh, CREA_GROUP_NO=(_F(GROUP_MA=nomma),)
            )

    def __set_list_supermaille(self, MacroElem):
        """Build a set of supermeshes"""
        list_sma = []
        key_index = 1
        list_me = [MacroElem.get_mael()]
        for mm in list_me:
            name = "STAT" + str(key_index)
            list_sma.append(_F(MACR_ELEM=mm, SUPER_MAILLE=name))
            self.supermaille.append(name)
            key_index += 1
        return list_sma

    def get_nb_ficti_no(self):
        """Count the number of fictitious cells and nodes to add to the mesh"""
        if self.macro_elem:
            mael = self.macro_elem.get_mael()
            Nb_no = mael.getNumberOfNodes()
        else:
            if self.param["PRE_CALC_MISS"]["NMAX_MODE_IFS"]:
                nb_modes_IFS = self.param["PRE_CALC_MISS"]["NMAX_MODE_IFS"]
            else:
                nb_modes_IFS = 0
            Nb_no = self.param["PRE_CALC_MISS"]["NMAX_MODE_ISS"] + nb_modes_IFS
        return Nb_no

    def check_ficti_nodes(self):
        """Check if fictitious cells and nodes should be added to the mesh"""
        if self.macro_elem:
            return self.macro_elem.check_reduc_dyna()
        elif (
            self.param["PRE_CALC_MISS"]["REDUC_DYNA_ISS"] == "OUI"
            or self.param["PRE_CALC_MISS"]["REDUC_DYNA_IFS"] == "OUI"
        ):
            return True
        else:
            return False

    def add_fiction_mesh(self):
        """Add fictitious cells and nodes to the mesh"""
        if not self.check_ficti_nodes():
            return None

        nb_new_nodes = self.get_nb_ficti_no()

        umesh = self.old_mesh.createMedCouplingMesh()

        coords = umesh.getCoords()
        nb_nodes = coords.getNumberOfTuples()
        coords.reAlloc(nb_nodes + nb_new_nodes)

        mesh0d = umesh.getMeshAtLevel(-3)
        nb_cells0d = mesh0d.getNumberOfCells()

        groups0d = [umesh.getGroupArr(-3, name) for name in umesh.getGroupsOnSpecifiedLev(-3)]

        cells = []
        for ii in range(nb_new_nodes):
            new_node = nb_nodes + ii
            new_cell = nb_cells0d + ii
            coords[new_node] = (0.0, 0.0, 1000000.0)
            mesh0d.insertNextCell(medc.NORM_POINT1, 1, [new_node])
            cells.append(new_cell)
        mesh0d.finishInsertingCells()
        mesh0d.checkConsistencyLight()

        group = medc.DataArrayInt(cells)
        group.setName("MFICTIF")
        groups0d.append(group)

        umesh.setMeshAtLevel(-3, mesh0d)
        umesh.setGroupsAtLevel(-3, groups0d)

        self.new_mesh = CA_Mesh()
        self.new_mesh.buildFromMedCouplingMesh(umesh)
