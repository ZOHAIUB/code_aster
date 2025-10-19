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
# Contribution from Naval Group

from math import *

from ...Cata.DataStructure import *
from ...Cata.DataStructure import CO as typCO
from ...Cata.Syntax import *
from ...Cata.Syntax import _F, MACRO, SIMP
from ...CodeCommands import *
from ...CodeCommands import AFFE_MATERIAU, ASSE_MATRICE, CALC_MATR_ELEM, DEFI_MATERIAU, NUME_DDL
from ...Supervis import UserMacro


def calc_matr_ifs_ops(self, **args):
    MAILLAGE = args.get("MAILLAGE")
    MODELE = args.get("MODELE")
    CHAR_CINE = args.get("CHAR_CINE")
    NUME_DDL_OUT = args.get("NUME_DDL")
    GROUP_MA_ELAS = args.get("GROUP_MA_ELAS")
    GROUP_MA_FLUI = args.get("GROUP_MA_FLUI")
    GROUP_MA_VISC = args.get("GROUP_MA_VISC")
    RHO_ELAS = args.get("RHO_ELAS")
    NU_ELAS = args.get("NU_ELAS")
    RHO_FLUI = args.get("RHO_FLUI")
    R_FLUI = args.get("R_FLUI")
    C_FLUI = args.get("C_FLUI")
    RHO_VISC = args.get("RHO_VISC")
    NU_VISC = args.get("NU_VISC")
    MASS_E_OUT = args.get("MASS_E")
    MASS_F_OUT = args.get("MASS_F")
    RIGI_E_OUT = args.get("RIGI_E")
    RIGI_F_OUT = args.get("RIGI_F")
    RIGI_V_OUT = args.get("RIGI_V")
    IMPE_R_OUT = args.get("IMPE_R")

    # Commandes utilisées dans la macro

    # Définition des propriétés matériaux
    _FLUI = DEFI_MATERIAU(FLUIDE=_F(RHO=RHO_FLUI, CELE_R=C_FLUI, LONG_CARA=R_FLUI))
    _FLUI_0 = DEFI_MATERIAU(FLUIDE=_F(RHO=RHO_FLUI, CELE_R=0))
    _FLUI_00 = DEFI_MATERIAU(FLUIDE=_F(RHO=0.0, CELE_R=0))

    _ELAS_0 = DEFI_MATERIAU(ELAS=_F(E=0.0, NU=NU_ELAS, RHO=RHO_ELAS, AMOR_HYST=0.0))
    _ELAS_1 = DEFI_MATERIAU(ELAS=_F(E=1.0, NU=NU_ELAS, RHO=RHO_ELAS, AMOR_HYST=0.0))
    _ELAS_11 = DEFI_MATERIAU(ELAS=_F(E=1.0, NU=NU_ELAS, RHO=1.0, AMOR_HYST=0.0))
    _ELAS_00 = DEFI_MATERIAU(ELAS=_F(E=0.0, NU=NU_ELAS, RHO=0.0, AMOR_HYST=0.0))

    if GROUP_MA_VISC is not None:
        _VISC_0 = DEFI_MATERIAU(ELAS=_F(E=0.0, NU=NU_VISC, RHO=RHO_VISC, AMOR_HYST=0.0))
        _VISC_1 = DEFI_MATERIAU(ELAS=_F(E=1.0, NU=NU_VISC, RHO=RHO_VISC, AMOR_HYST=0.0))
        _VISC_00 = DEFI_MATERIAU(ELAS=_F(E=0.0, NU=NU_VISC, RHO=0.0, AMOR_HYST=0.0))

    # Materiaux pour matrices de masse
    if GROUP_MA_VISC is None:
        _MAT_MSE = AFFE_MATERIAU(
            MAILLAGE=MAILLAGE,
            AFFE=(
                _F(GROUP_MA=GROUP_MA_ELAS, MATER=_ELAS_11),
                _F(GROUP_MA=GROUP_MA_FLUI, MATER=_FLUI_00),
            ),
        )
        _MAT_MSF = AFFE_MATERIAU(
            MAILLAGE=MAILLAGE,
            AFFE=(
                _F(GROUP_MA=GROUP_MA_ELAS, MATER=_ELAS_00),
                _F(GROUP_MA=GROUP_MA_FLUI, MATER=_FLUI),
            ),
        )
    else:
        _MAT_MSE = AFFE_MATERIAU(
            MAILLAGE=MAILLAGE,
            AFFE=(
                _F(GROUP_MA=GROUP_MA_ELAS, MATER=_ELAS_11),
                _F(GROUP_MA=GROUP_MA_FLUI, MATER=_FLUI_00),
                _F(GROUP_MA=GROUP_MA_VISC, MATER=_VISC_00),
            ),
        )
        _MAT_MSF = AFFE_MATERIAU(
            MAILLAGE=MAILLAGE,
            AFFE=(
                _F(GROUP_MA=GROUP_MA_ELAS, MATER=_ELAS_00),
                _F(GROUP_MA=GROUP_MA_FLUI, MATER=_FLUI),
                _F(GROUP_MA=GROUP_MA_VISC, MATER=_VISC_0),
            ),
        )

    # Materiaux pour matrices de raideur
    if GROUP_MA_VISC is None:
        _MAT_MK1 = AFFE_MATERIAU(
            MAILLAGE=MAILLAGE,
            AFFE=(
                _F(GROUP_MA=GROUP_MA_ELAS, MATER=_ELAS_0),
                _F(GROUP_MA=GROUP_MA_FLUI, MATER=_FLUI),
            ),
        )

        _MAT_MK2 = AFFE_MATERIAU(
            MAILLAGE=MAILLAGE,
            AFFE=(
                _F(GROUP_MA=GROUP_MA_ELAS, MATER=_ELAS_1),
                _F(GROUP_MA=GROUP_MA_FLUI, MATER=_FLUI_0),
            ),
        )

        _MAT_MK3 = AFFE_MATERIAU(
            MAILLAGE=MAILLAGE,
            AFFE=(
                _F(GROUP_MA=GROUP_MA_ELAS, MATER=_ELAS_0),
                _F(GROUP_MA=GROUP_MA_FLUI, MATER=_FLUI_0),
            ),
        )
    else:
        _MAT_MK1 = AFFE_MATERIAU(
            MAILLAGE=MAILLAGE,
            AFFE=(
                _F(GROUP_MA=GROUP_MA_ELAS, MATER=_ELAS_0),
                _F(GROUP_MA=GROUP_MA_FLUI, MATER=_FLUI),
                _F(GROUP_MA=GROUP_MA_VISC, MATER=_VISC_0),
            ),
        )

        _MAT_MK2 = AFFE_MATERIAU(
            MAILLAGE=MAILLAGE,
            AFFE=(
                _F(GROUP_MA=GROUP_MA_ELAS, MATER=_ELAS_1),
                _F(GROUP_MA=GROUP_MA_FLUI, MATER=_FLUI_0),
                _F(GROUP_MA=GROUP_MA_VISC, MATER=_VISC_0),
            ),
        )

        _MAT_MK3 = AFFE_MATERIAU(
            MAILLAGE=MAILLAGE,
            AFFE=(
                _F(GROUP_MA=GROUP_MA_ELAS, MATER=_ELAS_0),
                _F(GROUP_MA=GROUP_MA_FLUI, MATER=_FLUI_0),
                _F(GROUP_MA=GROUP_MA_VISC, MATER=_VISC_1),
            ),
        )
    # Creation du NUME_DDL

    numeddl = NUME_DDL(MODELE=MODELE)
    self.register_result(numeddl, NUME_DDL_OUT)

    # Calcul des matrices de masse

    _ElMASSe = CALC_MATR_ELEM(OPTION="MASS_MECA", MODELE=MODELE, CHAM_MATER=_MAT_MSE)
    _ElMASSf = CALC_MATR_ELEM(OPTION="MASS_MECA", MODELE=MODELE, CHAM_MATER=_MAT_MSF)

    MASS_E = ASSE_MATRICE(MATR_ELEM=_ElMASSe, NUME_DDL=numeddl, CHAR_CINE=CHAR_CINE)
    self.register_result(MASS_E, MASS_E_OUT)

    MASS_F = ASSE_MATRICE(MATR_ELEM=_ElMASSf, NUME_DDL=numeddl, CHAR_CINE=CHAR_CINE)
    self.register_result(MASS_F, MASS_F_OUT)

    # Calcul des matrices de rigidité

    _ElRIGIf = CALC_MATR_ELEM(OPTION="RIGI_MECA", MODELE=MODELE, CHAM_MATER=_MAT_MK1)

    _ElRIGIe = CALC_MATR_ELEM(OPTION="RIGI_MECA", MODELE=MODELE, CHAM_MATER=_MAT_MK2)

    RIGI_F = ASSE_MATRICE(MATR_ELEM=_ElRIGIf, NUME_DDL=numeddl, CHAR_CINE=CHAR_CINE)
    self.register_result(RIGI_F, RIGI_F_OUT)
    RIGI_E = ASSE_MATRICE(MATR_ELEM=_ElRIGIe, NUME_DDL=numeddl, CHAR_CINE=CHAR_CINE)
    self.register_result(RIGI_E, RIGI_E_OUT)

    _ElIMPEr = CALC_MATR_ELEM(OPTION="IMPE_MECA", MODELE=MODELE, CHAM_MATER=_MAT_MK1)

    IMPE_R = ASSE_MATRICE(MATR_ELEM=_ElIMPEr, NUME_DDL=numeddl, CHAR_CINE=CHAR_CINE)
    self.register_result(IMPE_R, IMPE_R_OUT)

    if GROUP_MA_VISC is not None:
        _ElRIGf1 = CALC_MATR_ELEM(OPTION="RIGI_MECA", MODELE=MODELE, CHAM_MATER=_MAT_MK3)

        _ElRIGf2 = CALC_MATR_ELEM(
            OPTION="RIGI_MECA_HYST", MODELE=MODELE, CHAM_MATER=_MAT_MK3, RIGI_MECA=_ElRIGf1
        )
        RIGI_V = ASSE_MATRICE(MATR_ELEM=_ElRIGf2, NUME_DDL=numeddl, CHAR_CINE=CHAR_CINE)
        self.register_result(RIGI_V, RIGI_V_OUT)


def calc_matr_ifs_prod(self, **args):
    MASS_E = args.get("MASS_E")
    MASS_F = args.get("MASS_F")
    RIGI_E = args.get("RIGI_E")
    RIGI_F = args.get("RIGI_F")
    RIGI_V = args.get("RIGI_V")
    IMPE_R = args.get("IMPE_R")
    NUME_DDL = args.get("NUME_DDL")
    if args.get("__all__"):
        return (
            [None],
            [None, nume_ddl_sdaster],
            [None, matr_asse_depl_r],
            [None, matr_asse_depl_r],
            [None, matr_asse_depl_r],
            [None, matr_asse_depl_r],
            [None, matr_asse_depl_c],
            [None, matr_asse_depl_r],
        )

    if not NUME_DDL:
        raise CataError("Impossible de typer les concepts resultats")

    if NUME_DDL.is_typco():
        self.type_sdprod(NUME_DDL, nume_ddl_sdaster)

    if MASS_E:
        self.type_sdprod(MASS_E, matr_asse_depl_r)

    if MASS_F:
        self.type_sdprod(MASS_F, matr_asse_depl_r)

    if RIGI_E:
        self.type_sdprod(RIGI_E, matr_asse_depl_r)

    if RIGI_F:
        self.type_sdprod(RIGI_F, matr_asse_depl_r)

    if RIGI_V:
        self.type_sdprod(RIGI_V, matr_asse_depl_c)

    if IMPE_R:
        self.type_sdprod(IMPE_R, matr_asse_depl_r)

    return


CALC_MATR_IFS_CATA = MACRO(
    nom="CALC_MATR_IFS",
    op=calc_matr_ifs_ops,
    sd_prod=calc_matr_ifs_prod,
    reentrant="n",
    # INPUT GENERAL DATA
    MAILLAGE=SIMP(statut="o", typ=maillage_sdaster),
    MODELE=SIMP(statut="o", typ=modele_sdaster),
    CHAR_CINE=SIMP(statut="f", typ=char_cine_meca),
    INFO=SIMP(statut="f", typ="I", defaut=1, into=(1, 2)),
    NUME_DDL=SIMP(statut="o", typ=(nume_ddl_sdaster, typCO)),
    # Zones d'affectation
    GROUP_MA_ELAS=SIMP(statut="o", typ=grma, validators=NoRepeat(), max="**"),
    GROUP_MA_FLUI=SIMP(statut="o", typ=grma, validators=NoRepeat(), max="**"),
    GROUP_MA_VISC=SIMP(statut="f", typ=grma, validators=NoRepeat(), max="**"),
    # MATERIAL PROPERTIES
    RHO_ELAS=SIMP(statut="o", typ="R"),
    NU_ELAS=SIMP(statut="o", typ="R"),
    RHO_FLUI=SIMP(statut="o", typ="R"),
    R_FLUI=SIMP(statut="o", typ="R", default=0.0),
    C_FLUI=SIMP(statut="o", typ="R"),
    RHO_VISC=SIMP(statut="f", typ="R"),
    NU_VISC=SIMP(statut="f", typ="R"),
    # Sorties (matrices)
    MASS_E=SIMP(statut="f", typ=typCO),
    MASS_F=SIMP(statut="f", typ=typCO),
    RIGI_E=SIMP(statut="f", typ=typCO),
    RIGI_F=SIMP(statut="f", typ=typCO),
    RIGI_V=SIMP(statut="f", typ=typCO),
    IMPE_R=SIMP(statut="f", typ=typCO),
)

CALC_MATR_IFS = UserMacro("CALC_MATR_IFS", CALC_MATR_IFS_CATA, calc_matr_ifs_ops)
