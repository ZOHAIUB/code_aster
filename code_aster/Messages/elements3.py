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
 On ne peut pas affecter des moments répartis sur des éléments de type %(k1)s.
"""
    ),
    2: _(
        """
 Le comportement ELAS_ISTR est interdit en 2D. Il faut utiliser ELAS_ORTH pour
 prendre en compte l'anisotropie.
"""
    ),
    3: _(
        """

 L'option MASS_FLUI_STRU n'est pas disponible pour les POU_D_TGM en multi-matériaux.
 Vous ne pouvez définir qu'un seul groupe de fibres.

"""
    ),
    13: _(
        """
  Les composantes SIXZ et SIYZ du champs de contraintes sont nulles pour les
  éléments de plaque, de coque et les tuyaux. Le calcul des composantes EPXZ et EPYZ du champs de déformation
  anélastiques donnerait des valeurs fausses. Ces valeurs sont donc mises à zéro et ne doivent pas être prises en compte.
"""
    ),
    16: _(
        """
 Comportement: %(k1)s non implanté
"""
    ),
    17: _(
        """
 Le matériau  %(k1)s  n'est pas connu
 Seuls sont admis les matériaux  'THER' et 'THER_COQUE' pour les coques thermiques
"""
    ),
    18: _(
        """
 Le matériau  %(k1)s  n'est pas connu
 Seuls sont admis les matériaux  'THER' et 'THER_COQUE' pour le calcul des flux pour les coques thermiques
"""
    ),
    26: _(
        """
 Mauvaise définition des caractéristiques de la section
"""
    ),
    30: _(
        """
 Section non tubulaire pour MASS_FLUI_STRU
"""
    ),
    31: _(
        """
 Pas de valeur utilisateur pour RHO
"""
    ),
    34: _(
        """
 Seules les forces suiveuses de type vent définies par un EVOL_CHAR sont autorisées
"""
    ),
    35: _(
        """
 Un champ de vitesse de vent est imposé sans donner un Cx dépendant de la vitesse sur une des barres
"""
    ),
    37: _(
        """
  Relation :  %(k1)s  non implantée sur les câbles
"""
    ),
    38: _(
        """
  Déformation :  %(k1)s  non implantée sur les câbles
"""
    ),
    39: _(
        """
 un champ de vitesse de vent est impose sans donner un Cx dépendant de la vitesse sur un des câbles.
"""
    ),
    49: _(
        """
 anisotropie non prévue pour coque1d
"""
    ),
    50: _(
        """
 nombre de couches limite a 30 pour les coques 1d
"""
    ),
    51: _(
        """
 Le nombre de couches défini dans DEFI_COMPOSITE et dans AFFE_CARA_ELEM dans n'est pas cohérent.
 Nombre de couches dans DEFI_COMPOSITE: %(i1)d
 Nombre de couches dans AFFE_CARA_ELEM: %(i2)d
"""
    ),
    52: _(
        """
 L'épaisseur totale des couches définie dans DEFI_COMPOSITE et celle définie dans AFFE_CARA_ELEM ne sont pas cohérentes.
 Épaisseur totale des couches dans DEFI_COMPOSITE: %(r1)f
 Épaisseur dans AFFE_CARA_ELEM: %(r2)f
"""
    ),
    56: _(
        """
 valeur utilisateur de RHO nulle
"""
    ),
    63: _(
        """
 ne pas utiliser THER_LINEAIRE avec des éléments de Fourier mais les commandes développées
"""
    ),
    67: _(
        """
 Élément dégénéré :
 revoir le maillage
"""
    ),
    74: _(
        """
 pour l'option "RICE_TRACEY", la relation " %(k1)s " n'est pas admise
"""
    ),
    91: _(
        """
  calcul de la masse non implanté pour les éléments COQUE_3D en grandes rotations, déformation : GROT_GDEP
"""
    ),
    93: _(
        """
  déformation :  %(k1)s  non implantée sur les éléments COQUE_3D en grandes rotations
  déformation : GROT_GDEP obligatoirement
"""
    ),
    98: _(
        """
 comportement coeur homogénéise inexistant
"""
    ),
    99: _(
        """
Seules les poutres à section constante sont admises !
"""
    ),
}
