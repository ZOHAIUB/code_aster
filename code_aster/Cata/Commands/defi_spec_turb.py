# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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


from ..Commons import *
from ..Language.DataStructure import *
from ..Language.Syntax import *

DEFI_SPEC_TURB = OPER(
    nom="DEFI_SPEC_TURB",
    op=145,
    sd_prod=spectre_sdaster,
    fr=tr("Definition d'un spectre d'excitation turbulente"),
    reentrant="n",
    regles=(
        UN_PARMI(
            "SPEC_LONG_COR_1",
            "SPEC_LONG_COR_2",
            "SPEC_LONG_COR_3",
            "SPEC_LONG_COR_4",
            "SPEC_CORR_CONV_1",
            "SPEC_CORR_CONV_2",
            "SPEC_CORR_CONV_3",
            "SPEC_FONC_FORME",
            "SPEC_EXCI_POINT",
        ),
    ),
    SPEC_LONG_COR_1=FACT(
        statut="f",
        LONG_COR=SIMP(statut="o", typ="R"),
        PROF_VITE_FLUI=SIMP(statut="o", typ=(fonction_sdaster, nappe_sdaster, formule)),
        VISC_CINE=SIMP(statut="o", typ="R"),
    ),
    SPEC_LONG_COR_2=FACT(
        statut="f",
        LONG_COR=SIMP(statut="o", typ="R"),
        PROF_VITE_FLUI=SIMP(statut="o", typ=(fonction_sdaster, nappe_sdaster, formule)),
        FREQ_COUP=SIMP(statut="f", typ="R", defaut=0.1),
        PHI0=SIMP(statut="f", typ="R", defaut=1.5e-3),
        BETA=SIMP(statut="f", typ="R", defaut=2.7),
    ),
    SPEC_LONG_COR_3=FACT(
        statut="f",
        LONG_COR=SIMP(statut="o", typ="R"),
        PROF_VITE_FLUI=SIMP(statut="o", typ=(fonction_sdaster, nappe_sdaster, formule)),
        FREQ_COUP=SIMP(statut="f", typ="R", defaut=0.2),
        PHI0_1=SIMP(statut="f", typ="R", defaut=5.0e-3),
        BETA_1=SIMP(statut="f", typ="R", defaut=0.5),
        PHI0_2=SIMP(statut="f", typ="R", defaut=4.0e-5),
        BETA_2=SIMP(statut="f", typ="R", defaut=3.5),
    ),
    SPEC_LONG_COR_4=FACT(
        statut="f",
        LONG_COR=SIMP(statut="o", typ="R"),
        PROF_VITE_FLUI=SIMP(statut="o", typ=(fonction_sdaster, nappe_sdaster, formule)),
        TAUX_VIDE=SIMP(statut="o", typ="R"),
        BETA=SIMP(statut="f", typ="R", defaut=2.0),
        GAMMA=SIMP(statut="f", typ="R", defaut=4.0),
    ),
    SPEC_CORR_CONV_1=FACT(
        statut="f",
        LONG_COR_1=SIMP(statut="f", typ="R", defaut=999),
        LONG_COR_2=SIMP(statut="f", typ="R"),
        VITE_FLUI=SIMP(statut="o", typ="R"),
        RHO_FLUI=SIMP(statut="o", typ="R"),
        FREQ_COUP=SIMP(statut="f", typ="R", defaut=1000),
        K=SIMP(statut="f", typ="R", defaut=5.8e-3),
        D_FLUI=SIMP(statut="o", typ="R"),
        COEF_VITE_FLUI_A=SIMP(statut="f", typ="R", defaut=0.65),
        COEF_VITE_FLUI_O=SIMP(statut="f", typ="R", defaut=999),
        METHODE=SIMP(statut="o", typ="TXM", into=("AU_YANG", "CORCOS")),
    ),
    SPEC_CORR_CONV_2=FACT(
        statut="f",
        FONCTION=SIMP(statut="o", typ=(fonction_sdaster, nappe_sdaster, formule)),
        LONG1_F=SIMP(statut="o", typ=(fonction_sdaster, nappe_sdaster, formule)),
        LONG2_F=SIMP(statut="o", typ=(fonction_sdaster, nappe_sdaster, formule)),
        COEF_VITE_FLUI_A=SIMP(statut="f", typ="R", defaut=0.65),
        VITE_FLUI=SIMP(statut="o", typ="R"),
        COEF_VITE_FLUI_O=SIMP(statut="f", typ="R", defaut=999.0),
        FREQ_COUP=SIMP(statut="f", typ="R", defaut=1000.0),
        METHODE=SIMP(statut="o", typ="TXM", into=("AU_YANG", "CORCOS")),
    ),
    SPEC_CORR_CONV_3=FACT(statut="f", TABLE_FONCTION=SIMP(statut="o", typ=(table_fonction))),
    SPEC_FONC_FORME=FACT(
        statut="f",
        regles=(
            UN_PARMI("INTE_SPEC", "GRAPPE_1"),
            ENSEMBLE("INTE_SPEC", "FONCTION"),
            UN_PARMI("NOEUD", "GROUP_NO"),
            EXCLUS("NOEUD", "GROUP_NO"),
        ),
        INTE_SPEC=SIMP(statut="f", typ=interspectre),
        FONCTION=SIMP(statut="f", typ=(table_fonction), max="**"),
        GRAPPE_1=SIMP(statut="f", typ="TXM", into=("DEBIT_180", "DEBIT_300")),
        NOEUD=SIMP(statut="c", typ=no),
        GROUP_NO=SIMP(statut="f", typ=grno),
        CARA_ELEM=SIMP(statut="o", typ=cara_elem),
        MODELE=SIMP(statut="o", typ=modele_sdaster),
    ),
    SPEC_EXCI_POINT=FACT(
        statut="f",
        regles=(UN_PARMI("INTE_SPEC", "GRAPPE_2"),),
        INTE_SPEC=SIMP(statut="f", typ=interspectre),
        GRAPPE_2=SIMP(statut="f", typ="TXM", into=("ASC_CEN", "ASC_EXC", "DES_CEN", "DES_EXC")),
        #  Quels sont les statuts des mots cles a l interieur des deux blocs qui suivent
        b_inte_spec=BLOC(
            condition="""exists("INTE_SPEC")""",
            regles=(UN_PARMI("NOEUD", "GROUP_NO"), EXCLUS("NOEUD", "GROUP_NO")),
            NATURE=SIMP(statut="o", typ="TXM", max="**", into=("FORCE", "MOMENT")),
            ANGLE=SIMP(statut="o", typ="R", max="**"),
            NOEUD=SIMP(statut="c", typ=no, max="**"),
            GROUP_NO=SIMP(statut="f", typ=grno, max="**"),
        ),
        b_grappe_2=BLOC(
            condition="""exists("GRAPPE_2")""",
            regles=(UN_PARMI("NOEUD", "GROUP_NO"), EXCLUS("NOEUD", "GROUP_NO")),
            RHO_FLUI=SIMP(statut="o", typ="R"),
            NOEUD=SIMP(statut="c", typ=no),
            GROUP_NO=SIMP(statut="f", typ=grno),
        ),
        CARA_ELEM=SIMP(statut="o", typ=cara_elem),
        MODELE=SIMP(statut="o", typ=modele_sdaster),
    ),
    TITRE=SIMP(statut="f", typ="TXM"),
)
