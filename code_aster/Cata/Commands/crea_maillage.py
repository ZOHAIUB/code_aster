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

# person_in_charge: jacques.pellet at edf.fr

from ..Commons import *
from ..Language.DataStructure import *
from ..Language.Syntax import *


def cream_prod(self, **args):
    if args.get("__all__"):
        return (maillage_sdaster, maillage_p)

    if "MAILLAGE" in args:
        if AsType(args["MAILLAGE"]) == "maillage_p":
            return maillage_p

    return maillage_sdaster


CREA_MAILLAGE = OPER(
    nom="CREA_MAILLAGE",
    op=167,
    sd_prod=cream_prod,
    reentrant="n",
    fr=tr("Crée un maillage à partir d'un maillage existant"),
    regles=(
        UN_PARMI(
            "COQU_VOLU",
            "CREA_FISS",
            "CREA_MAILLE",
            "CREA_POI1",
            "ECLA_PG",
            "HEXA20_27",
            "LINE_QUAD",
            "MODI_MAILLE",
            "QUAD_LINE",
            "REPERE",
            "RESTREINT",
            "PENTA15_18",
            "GEOM_FIBRE",
            "DECOUPE_LAC",
            "MODI_HHO",
            "COQUE_SOLIDE",
            "RAFFINEMENT",
        ),
    ),
    # le MAILLAGE est inutile si ECLA_PG et GEOM_FIBRE
    MAILLAGE=SIMP(statut="f", typ=(maillage_sdaster, maillage_p)),
    GEOM_FIBRE=SIMP(statut="f", max=1, typ=gfibre_sdaster),
    CREA_POI1=FACT(
        statut="f",
        max="**",
        fr=tr("Création de mailles de type POI1 à partir de noeuds"),
        regles=(AU_MOINS_UN("TOUT", "GROUP_MA", "GROUP_NO", "NOEUD"),),
        NOM_GROUP_MA=SIMP(statut="f", typ=grma, validators=NoRepeat(), max="**"),
        TOUT=SIMP(statut="f", typ="TXM", into=("OUI",)),
        GROUP_MA=SIMP(statut="f", typ=grma, validators=NoRepeat(), max="**"),
        GROUP_NO=SIMP(statut="f", typ=grno, validators=NoRepeat(), max="**"),
        NOEUD=SIMP(statut="c", typ=no, validators=NoRepeat(), max="**"),
    ),
    CREA_MAILLE=FACT(
        statut="f",
        max="**",
        fr=tr("Duplication de mailles"),
        regles=(UN_PARMI("TOUT", "GROUP_MA", TOUT="OUI"),),
        NOM=SIMP(statut="o", typ="TXM"),
        TOUT=SIMP(statut="f", typ="TXM", into=("OUI",)),
        GROUP_MA=SIMP(statut="f", typ=grma, validators=NoRepeat(), max="**"),
    ),
    RESTREINT=FACT(
        statut="f",
        fr=tr("Restreindre un maillage à des groupes de mailles"),
        max=1,
        GROUP_MA=SIMP(statut="o", typ=grma, validators=NoRepeat(), max="**"),
        TOUT_GROUP_MA=SIMP(statut="f", typ="TXM", defaut="NON", into=("OUI", "NON")),
        GROUP_NO=SIMP(statut="f", typ=grno, validators=NoRepeat(), max="**"),
        TOUT_GROUP_NO=SIMP(statut="f", typ="TXM", defaut="NON", into=("OUI", "NON")),
    ),
    COQU_VOLU=FACT(
        statut="f",
        fr=tr("Creation de mailles volumiques à partir de mailles surfaciques"),
        NOM=SIMP(statut="o", typ="TXM"),
        GROUP_MA=SIMP(statut="o", typ=grma, validators=NoRepeat(), max="**"),
        EPAIS=SIMP(statut="o", typ="R"),
        PLAN=SIMP(statut="o", typ="TXM", into=("SUP", "MOY", "INF")),
        b_MOY=BLOC(
            condition="""equal_to("PLAN", 'MOY')""",
            TRANSLATION=SIMP(statut="o", typ="TXM", into=("SUP", "INF")),
        ),
    ),
    MODI_MAILLE=FACT(
        statut="f",
        max="**",
        fr=tr("Modification du type de mailles"),
        regles=(AU_MOINS_UN("TOUT", "MAILLE", "GROUP_MA"),),
        TOUT=SIMP(statut="f", typ="TXM", into=("OUI",)),
        MAILLE=SIMP(statut="c", typ=ma, validators=NoRepeat(), max="**"),
        GROUP_MA=SIMP(statut="f", typ=grma, validators=NoRepeat(), max="**"),
        OPTION=SIMP(
            statut="o",
            typ="TXM",
            into=("TRIA6_7", "QUAD8_9", "SEG3_4", "QUAD_TRIA3"),
            validators=NoRepeat(),
            fr=tr("Choix de la transformation"),
        ),
    ),
    RAFFINEMENT=FACT(
        statut="f",
        fr=tr("Raffinement uniforme du maillage"),
        TOUT=SIMP(statut="f", typ="TXM", into=("OUI",), defaut="OUI"),
        NIVEAU=SIMP(statut="f", typ="I", defaut=1),
    ),
    CREA_FISS=FACT(
        statut="f",
        max="**",
        fr=tr("Creation d'une fissure potentielle avec elts de joint ou elts à disc"),
        NOM=SIMP(statut="o", typ="TXM"),
        GROUP_NO_1=SIMP(statut="o", typ=grno),
        GROUP_NO_2=SIMP(statut="o", typ=grno),
    ),
    LINE_QUAD=FACT(
        statut="f",
        fr=tr("Passage linéaire -> quadratique"),
        regles=(UN_PARMI("TOUT", "GROUP_MA", TOUT="OUI"),),
        TOUT=SIMP(statut="f", typ="TXM", into=("OUI",)),
        GROUP_MA=SIMP(statut="f", typ=grma, validators=NoRepeat(), max="**"),
    ),
    HEXA20_27=FACT(
        statut="f",
        fr=tr("Passage HEXA20 -> HEXA27"),
        regles=(UN_PARMI("TOUT", "GROUP_MA", TOUT="OUI"),),
        TOUT=SIMP(statut="f", typ="TXM", into=("OUI",)),
        GROUP_MA=SIMP(statut="f", typ=grma, validators=NoRepeat(), max="**"),
    ),
    PENTA15_18=FACT(
        statut="f",
        fr=tr("Passage PENTA15 -> PENTA18"),
        regles=(UN_PARMI("TOUT", "GROUP_MA", TOUT="OUI"),),
        TOUT=SIMP(statut="f", typ="TXM", into=("OUI",)),
        GROUP_MA=SIMP(statut="f", typ=grma, validators=NoRepeat(), max="**"),
    ),
    MODI_HHO=FACT(
        statut="f",
        fr=tr("Création maillage pour HHO"),
        TOUT=SIMP(statut="f", typ="TXM", into=("OUI",)),
    ),
    COQUE_SOLIDE=FACT(
        statut="f",
        max="**",
        fr=tr("Création maillage pour COQUE_SOLIDE"),
        GROUP_MA=SIMP(statut="o", typ=grma, validators=NoRepeat(), max="**"),
        GROUP_MA_SURF=SIMP(statut="f", typ=grma, validators=NoRepeat(), max="**"),
    ),
    QUAD_LINE=FACT(
        statut="f",
        fr=tr("Passage quadratique -> linéaire"),
        regles=(UN_PARMI("TOUT", "GROUP_MA", TOUT="OUI"),),
        TOUT=SIMP(statut="f", typ="TXM", into=("OUI",)),
        GROUP_MA=SIMP(statut="f", typ=grma, validators=NoRepeat(), max="**"),
    ),
    REPERE=FACT(
        statut="f",
        fr=tr(
            "changement de repère servant à déterminer les caractéristiques d'une section de poutre"
        ),
        TABLE=SIMP(
            statut="o",
            typ=table_sdaster,
            fr=tr("Nom de la table contenant les caractéristiques de la section de poutre"),
        ),
        NOM_ORIG=SIMP(
            statut="f", typ="TXM", into=("CDG", "TORSION"), fr=tr("Origine du nouveau repère")
        ),
        NOM_ROTA=SIMP(statut="f", typ="TXM", into=("INERTIE",), fr=tr("Direction du repére")),
        b_cdg=BLOC(
            condition="""equal_to("NOM_ORIG", 'CDG')""",
            GROUP_MA=SIMP(
                statut="f",
                typ=grma,
                fr=tr(
                    "Nom du groupe de mailles dont le centre de gravité sera l origine du nouveau repère"
                ),
            ),
        ),
    ),
    DECOUPE_LAC=FACT(
        statut="f",
        fr="creation des patchs LAC pour le  contact",
        regles=(AU_MOINS_UN("GROUP_MA_ESCL"),),
        GROUP_MA_ESCL=SIMP(statut="o", typ=grma, validators=NoRepeat(), max="**"),
        DECOUPE_HEXA=SIMP(
            statut="f", typ="TXM", into=("PYRA", "HEXA"), defaut="PYRA", validators=NoRepeat()
        ),
    ),
    ECLA_PG=FACT(
        statut="f",
        fr=tr("Eclatement des mailles en petites mailles contenant chacune un seul point de gauss"),
        MODELE=SIMP(statut="o", typ=modele_sdaster),
        TOUT=SIMP(statut="f", typ="TXM", into=("OUI",)),
        GROUP_MA=SIMP(statut="f", typ=grma, validators=NoRepeat(), max="**"),
        SHRINK=SIMP(statut="f", typ="R", defaut=0.9, fr=tr("Facteur de réduction")),
        TAILLE_MIN=SIMP(statut="f", typ="R", defaut=0.0, fr=tr("Taille minimale d'un coté")),
        NOM_CHAM=SIMP(
            statut="o",
            typ="TXM",
            validators=NoRepeat(),
            max="**",
            into=C_NOM_CHAM_INTO(type_cham="ELGA"),
        ),
    ),
    TITRE=SIMP(statut="f", typ="TXM"),
    #
    INFO=SIMP(statut="f", typ="I", defaut=1, into=(1, 2)),
)
