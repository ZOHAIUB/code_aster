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
    CALC_MATR_ELEM,
    CALC_MODES,
    DEFI_BASE_MODALE,
    DEFI_INTERF_DYNA,
    EXTR_MODE,
    MACR_ELEM_DYNA,
    MODE_STATIQUE,
    NUME_DDL,
)


def crea_elem_ssd_ops(self, **args):
    """
    Enchainement des commandes :
       CALC_MATR_ELEM + ASSE_MATRICE + CALC_MODES + MODE_STATIQUE
       DEFI_INTERF_DYNA + DEFI_BASE_MODALE + MACR_ELEM_DYNA
    """

    NUMEDDL = args.get("NUME_DDL")
    INTERFACE = args.get("INTERFACE")
    BASE_MODALE = args.get("BASE_MODALE")
    CALC_FREQ = args.get("CALC_FREQ")
    SOLVEUR = args.get("SOLVEUR")

    mSolveur = SOLVEUR[0].cree_dict_valeurs(SOLVEUR[0].mc_liste)

    _kelem = CALC_MATR_ELEM(
        CHARGE=args["CHARGE"],
        OPTION="RIGI_MECA",
        CARA_ELEM=args.get("CARA_ELEM"),
        MODELE=args["MODELE"],
        CHAM_MATER=args["CHAM_MATER"],
    )

    __melem = CALC_MATR_ELEM(
        CHARGE=args["CHARGE"],
        OPTION="MASS_MECA",
        CARA_ELEM=args.get("CARA_ELEM"),
        MODELE=args["MODELE"],
        CHAM_MATER=args["CHAM_MATER"],
    )

    _nume_ddl = NUME_DDL(MATR_RIGI=_kelem)

    if NUMEDDL:
        self.register_result(_nume_ddl, NUMEDDL)

    _matrigi = ASSE_MATRICE(NUME_DDL=_nume_ddl, MATR_ELEM=_kelem)

    _matmass = ASSE_MATRICE(NUME_DDL=_nume_ddl, MATR_ELEM=__melem)

    # recuperation des options de CALC_MODES
    motscit = {}
    motscfa = {}
    motscsm = {}

    NBANDE_OPTIONS = {"PLUS_PETITE": 1, "CENTRE": 1, "SANS": 0}

    nbande = NBANDE_OPTIONS.get(CALC_FREQ["OPTION"])

    if CALC_FREQ["OPTION"] == "BANDE":
        nbande = len(CALC_FREQ["FREQ"]) - 1

    if CALC_FREQ["DIM_SOUS_ESPACE"]:
        motscsm["DIM_SOUS_ESPACE"] = CALC_FREQ["DIM_SOUS_ESPACE"]

    motscit["VERI_MODE"] = _F(STOP_ERREUR=CALC_FREQ["STOP_ERREUR"])

    motfilt = {}
    motfilt["FILTRE_MODE"] = []
    for i in range(nbande):
        if CALC_FREQ["OPTION"] == "PLUS_PETITE":
            motscfa["NMAX_FREQ"] = CALC_FREQ["NMAX_FREQ"]

        if CALC_FREQ["OPTION"] == "CENTRE":
            motscfa["FREQ"] = CALC_FREQ["FREQ"]
            if CALC_FREQ["AMOR_REDUIT"]:
                motscfa["AMOR_REDUIT"] = CALC_FREQ["AMOR_REDUIT"]
            motscfa["NMAX_FREQ"] = CALC_FREQ["NMAX_FREQ"]

        if CALC_FREQ["OPTION"] == "BANDE":
            motscfa["FREQ"] = (CALC_FREQ["FREQ"][i], CALC_FREQ["FREQ"][i + 1])

        motscit["CALC_FREQ"] = _F(**motscfa)

        if CALC_FREQ["APPROCHE"] is not None:
            motscsm["APPROCHE"] = CALC_FREQ["APPROCHE"]

        if len(motscsm) > 0:
            motscit["SOLVEUR_MODAL"] = _F(**motscsm)

        __modes = CALC_MODES(
            MATR_RIGI=_matrigi,
            MATR_MASS=_matmass,
            OPTION=CALC_FREQ["OPTION"],
            INFO=args["INFO"],
            **motscit
        )

        motfilt["FILTRE_MODE"].append(_F(MODE=__modes, TOUT_ORDRE="OUI"))

    if nbande:
        _mode_meca = EXTR_MODE(**motfilt)

    if BASE_MODALE[0]["TYPE"] == "RITZ":
        mcfactc = []
        mcfactm = []
        mcfacti = []
        arg_no = []
        arg_grno = []
        for i in range(len(INTERFACE)):
            if BASE_MODALE[0]["TYPE_MODE"] == "INTERFACE":
                if INTERFACE[i]["TYPE"] == "CRAIGB":
                    if INTERFACE[i]["NOEUD"]:
                        if isinstance(INTERFACE[i]["NOEUD"], (list, tuple)):
                            for noeu in INTERFACE[i]["NOEUD"]:
                                arg_no.append(noeu)
                        else:
                            arg_no.append(INTERFACE[i]["NOEUD"])
                    if INTERFACE[i]["GROUP_NO"]:
                        if isinstance(INTERFACE[i]["GROUP_NO"], (list, tuple)):
                            for grno in INTERFACE[i]["GROUP_NO"]:
                                arg_grno.append(grno)
                        else:
                            arg_grno.append(INTERFACE[i]["GROUP_NO"])
            else:
                arg_int = {}
                if INTERFACE[i]["NOEUD"]:
                    arg_int["NOEUD"] = INTERFACE[i]["NOEUD"]
                if INTERFACE[i]["GROUP_NO"]:
                    arg_int["GROUP_NO"] = INTERFACE[i]["GROUP_NO"]
                arg_int["TOUT_CMP"] = "OUI"
                if INTERFACE[i]["TYPE"] == "CRAIGB":
                    mcfactc.append(_F(**arg_int))
                elif INTERFACE[i]["TYPE"] == "MNEAL":
                    mcfactm.append(_F(**arg_int))
        modstatc = {}
        modstatm = {}
        modstati = {}
        lmodint = []
        if mcfactc:
            modstatc["MODE_STAT"] = mcfactc
            _mode_intf = MODE_STATIQUE(MATR_RIGI=_matrigi, SOLVEUR=mSolveur, **modstatc)
            lmodint.append(_mode_intf)
        if mcfactm:
            modstatm["FORCE_NODALE"] = mcfactm
            _mode_intf = MODE_STATIQUE(MATR_RIGI=_matrigi, SOLVEUR=mSolveur, **modstatm)
            lmodint.append(_mode_intf)
        if BASE_MODALE[0]["TYPE_MODE"] == "INTERFACE":
            arg_int = {}
            if arg_no:
                arg_int["NOEUD"] = arg_no
            if arg_grno:
                arg_int["GROUP_NO"] = arg_grno
            arg_int["NB_MODE"] = BASE_MODALE[0]["NMAX_MODE_INTF"]
            arg_int["TOUT_CMP"] = "OUI"
            mcfacti.append(_F(**arg_int))
            modstati["MODE_INTERF"] = mcfacti
            _mode_intf = MODE_STATIQUE(
                MATR_RIGI=_matrigi, MATR_MASS=_matmass, SOLVEUR=mSolveur, **modstati
            )
            lmodint.append(_mode_intf)

    interface = {}
    mcfact = []
    freqnew = None
    ifreq = INTERFACE[0]["FREQ"]
    for i in range(len(INTERFACE)):
        arg_int = {}
        if INTERFACE[i]["NOEUD"]:
            arg_int["NOEUD"] = INTERFACE[i]["NOEUD"]
        if INTERFACE[i]["MASQUE"]:
            arg_int["MASQUE"] = INTERFACE[i]["MASQUE"]
        if INTERFACE[i]["GROUP_NO"]:
            arg_int["GROUP_NO"] = INTERFACE[i]["GROUP_NO"]
        mcfact.append(_F(NOM=INTERFACE[i]["NOM"], TYPE=INTERFACE[i]["TYPE"], **arg_int))
        ifreq_i = INTERFACE[i]["FREQ"]
        if ifreq != ifreq_i:
            freqnew = ifreq_i
    if freqnew:
        UTMESS("A", "SOUSTRUC2_12", valr=freqnew)
        ifreq = freqnew
    interface["INTERFACE"] = mcfact

    if args["INFO"]:
        interface["INFO"] = args["INFO"]
    if ifreq:
        interface["FREQ"] = ifreq

    _interf = DEFI_INTERF_DYNA(NUME_DDL=_nume_ddl, **interface)

    base = {}
    mcfact = []
    arg_base = {}
    if args["INFO"]:
        base["INFO"] = args["INFO"]

    if BASE_MODALE[0]["TYPE"] == "CLASSIQUE":
        type_base = "CLASSIQUE"
        arg_base["NMAX_MODE"] = CALC_FREQ[0]["NMAX_FREQ"]
        mcfact.append(_F(INTERF_DYNA=_interf, MODE_MECA=_mode_meca, **arg_base))

    if BASE_MODALE[0]["TYPE"] == "RITZ":

        type_base = "RITZ"
        if BASE_MODALE[0]["TYPE_MODE"] == "STATIQUE":
            mcfact.append(_F(MODE_MECA=_mode_intf))
        else:
            if nbande:
                if CALC_FREQ[0]["OPTION"] == "PLUS_PETITE" or CALC_FREQ[0]["OPTION"] == "CENTRE":
                    arg_base["NMAX_MODE"] = CALC_FREQ[0]["NMAX_FREQ"]
                mcfact.append(_F(MODE_MECA=_mode_meca, **arg_base))
            else:
                arg_base["NMAX_MODE"] = 0
                mcfact.append(_F(MODE_MECA=_mode_intf, **arg_base))

            arg_base = {}
            if BASE_MODALE[0]["NMAX_MODE_INTF"]:
                arg_base["NMAX_MODE"] = BASE_MODALE[0]["NMAX_MODE_INTF"]
            mcfact.append(_F(MODE_INTF=_mode_intf, **arg_base))

    if type_base == "CLASSIQUE":
        base["CLASSIQUE"] = mcfact
    elif type_base == "RITZ":
        base["RITZ"] = mcfact
        base["INTERF_DYNA"] = _interf
        base["NUME_REF"] = _nume_ddl

    _base_modale = DEFI_BASE_MODALE(**base)
    elem = {"MATR_RIGI": _matrigi, "MATR_MASS": _matmass}
    macr_elem = MACR_ELEM_DYNA(BASE_MODALE=_base_modale, **elem)

    return macr_elem
