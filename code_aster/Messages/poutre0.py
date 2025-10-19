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

from ..Utilities import _

cata_msg = {
    1: _(
        """
Problème lors de l'utilisation de MACR_CARA_POUTRE

GROUP_MA et GROUP_MA_BORD doivent avoir le même nombre de groupes de mailles.
À un GROUP_MA doit correspondre un GROUP_MA_BORD.
"""
    ),
    3: _(
        """
Problème lors de l'utilisation de MACR_CARA_POUTRE

Vous avez renseigné NOEUD ou GROUP_NO avec GROUP_MA_BORD.
Il faut que NOEUD ou GROUP_NO contienne un noeud unique.

Si vous avez des groupes de mailles disjoints, il faut également renseigner GROUP_MA.
Et dans ce cas à un NOEUD doit correspondre un GROUP_MA et un GROUP_MA_BORD.
"""
    ),
    4: _(
        """
Poutre circulaire à variation de section homothétique.

Le rapport d'homothétie est assez différent entre les rayons et les épaisseurs :
    - rayon 1 = %(r1)16.9g
    - rayon 2 = %(r2)16.9g
        `- rapport d'homothétie = %(r5)16.9g

    - épaisseur 1 = %(r3)16.9g
    - épaisseur 2 = %(r4)16.9g
        `- rapport d'homothétie = %(r6)16.9g

La différence entre les rapports d'homothétie est supérieure à 1%%.
Les hypothèses du modèle de poutre à variation homothétique ne sont donc pas
respectées (consultez la documentation de référence des éléments poutre).


Risques et conseil:
    - Les calculs seront inexacts.
    - Raffiner le maillage permet de minimiser les écarts.
"""
    ),
    5: _(
        """
Problème lors de l'utilisation de MACR_CARA_POUTRE

Vous avez renseigné GROUP_NO avec GROUP_MA et GROUP_MA_BORD.
Le nombre de noeud dans GROUP_NO doit être identique au nombre de groupe de mailles
donné sous GROUP_MA et GROUP_MA_BORD.
"""
    ),
    6: _(
        """
Problème lors de l'utilisation de MACR_CARA_POUTRE

Vous avez renseigné GROUP_MA et GROUP_MA_BORD avec plusieurs groupes de mailles.
Vous devez renseigner 'LONGUEUR' 'MATERIAU' 'LIAISON', pour que l'on puisse calculer les
caractéristiques de torsion en prenant en compte les facteurs de participation de chaque
groupe de mailles.
"""
    ),
    7: _(
        """
Vous avez renseigné GROUP_MA et GROUP_MA_BORD avec un seul groupe de mailles par
mot clef.
Vous avez également renseigné 'LONGUEUR' 'MATERIAU' 'LIAISON', cela est inutile.
Les valeurs sont ignorées.
"""
    ),
    8: _(
        """
Problème lors de l'utilisation de MACR_CARA_POUTRE

Vous avez renseigné GROUP_NO. Le GROUP_NO %(k1)s n'existe pas dans le maillage.
"""
    ),
    10: _(
        """La caractéristique %(k1)8s est négative ou nulle %(r1)e
"""
    ),
    11: _(
        """
Problème lors de l'utilisation de MACR_CARA_POUTRE

La section présente des caractéristiques mécaniques négatives.

Vous avez renseigné l'option TABLE_CARA='OUI'. Si vous utilisez les résultats dans la commande
AFFE_CARA_ELEM / POUTRE avec TABLE_CARA, vous risquez d'avoir des résultats faux.

Conseil : La discrétisation de la section a un impact sur la qualité des résultats.
          Raffiner le maillage devrait permettre de lever cette erreur.
"""
    ),
    12: _(
        """
Les coordonnées du centre de gravité de la section sont G=(%(r1)e, %(r2)e)

Si vous utilisez des MULTIFIBRES et que le maillage de description des fibres est le même que
celui utilisé dans cette commande, il faut dans la commande DEFI_GEOM_FIBRE renseigner soit
le mot clef COOR_AXE_POUTRE soit le mot clef TABLE_CARA et NOM_SEC de façon à faire
correspondre le centre de gravité des fibres à l'axe neutre de la poutre :
    COOR_AXE_POUTRE = (%(r1)e, %(r2)e)
ou
    TABLE_CARA = %(k1)s, NOM_SEC = "%(k2)s",

Vous avez renseigné l'option TABLE_CARA='OUI'. Si vous utilisez les résultats dans la
commande : AFFE_CARA_ELEM / POUTRE avec TABLE_CARA et NOM_SEC,
                           / GEOM_FIBRE
                           / MULTIFIBRE

Vous risquez d'avoir des résultats inattendus, si vous ne renseignez pas :
    COOR_AXE_POUTRE ou TABLE_CARA et NOM_SEC
"""
    ),
    13: _(
        """
Le repère principal d'inertie est tourné d'un angle de %(r1)f° par rapport aux axes du maillage.

Si vous utilisez des MULTIFIBRES et que le maillage de description des fibres est le même que
celui utilisé dans cette commande, il faut dans la commande DEFI_GEOM_FIBRE renseigner soit
le mot clef ANGLE soit le mot clef TABLE_CARA et NOM_SEC de façon à faire correspondre le
repère principal d'inertie aux axes du maillage de la section :
    ANGLE = %(r2)e
ou
    TABLE_CARA = %(k1)s, NOM_SEC = "%(k2)s",

Vous avez renseigné l'option TABLE_CARA='OUI'. Si vous utilisez les résultats dans la commande
AFFE_CARA_ELEM / POUTRE avec TABLE_CARA, il est peut-être nécessaire dans la commande
AFFE_CARA_ELEM de renseigner le mot clef ORIENTATION de façon à définir la position du repère
principal d'inertie par rapport au repère global.

ORIENTATION=(
    _F(GROUP_MA='....', CARA='ANGL_VRIL', VALE= ??? )
)
Par défaut ANGL_VRIL= 0

Vous risquez d'avoir des résultats inattendus, si vous ne renseignez ni :
 - ANGLE ou  TABLE_CARA et NOM_SEC dans DEFI_GEOM_FIBRE
 - ORIENTATION dans AFFE_CARA_ELEM
"""
    ),
    14: _(
        """Vous utilisez un élément de type multifibre. Il faut que sous COMPORTEMENT le mot clef RELATION soit 'MULTIFIBRE'."""
    ),
    16: _(
        """Vous utilisez un élément de type multifibre avec DEFORMATION='%(k1)s'. L'option RIGI_GEOM='OUI' n'est pas autorisée dans ce cas."""
    ),
    17: _("""L'élément n'est utilisable qu'en élasticité."""),
    20: _(
        """
Problème lors de l'utilisation de MACR_CARA_POUTRE

Vous avez renseigné %(k2)s. Le GROUP_MA %(k1)s n'existe pas dans le maillage.
"""
    ),
    40: _("""Les éléments de poutre ne peuvent pas utiliser une déformation de type %(k1)s."""),
    41: _(
        """Ces éléments de poutre ne peuvent être utilisés qu'en grandes transformations (GROT_GDEP)."""
    ),
    42: _("""Il manque la masse volumique RHO qui est nécessaire en dynamique."""),
    43: _(
        """La matrice sécante n'est pas disponible avec le comportement %(k1)s sur les éléments barre."""
    ),
    49: _(
        """La méthode IMPLEX ne peut pas être utilisée avec la loi de comportement %(k1)s sur les éléments barre."""
    ),
    50: _(
        """On ne trouve pas de comportement élastique sur la poutre. C'est nécessaire pour récupérer la masse volumique."""
    ),
    59: _(
        """Le coefficient de poisson n'est pas constant. Les éléments de poutre n'en tiennent pas compte."""
    ),
    61: _(
        """La loi de comportement %(k1)s n'est pas disponible pour les poutres classiques. Utilisez les poutres multifibres."""
    ),
    62: _("""Les noeuds sont confondus pour un élément de barre."""),
    64: _(
        """Avec l'option GROT_GDEP, les coefficients de flexibilité ne sont pas pris en compte dans la matrice de raideur géométrique.
   Coefficient de flexibilité suivant y : %(r1)f
   Coefficient de flexibilité suivant z : %(r2)f"""
    ),
    81: _(
        """
Pour les éléments de poutre à section variable affine : seule une section rectangle plein est possible.
"""
    ),
    82: _(
        """
Pour les éléments de poutre à section variable homothétique : l'aire initiale est nulle.
"""
    ),
    90: _(
        """Le seul comportement élastique valide est ELAS pour les éléments de poutre squelette."""
    ),
}
