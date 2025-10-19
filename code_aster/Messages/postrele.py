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

from ..Utilities import _

cata_msg = {
    1: _(
        """
 création/extension de la table %(k1)s
"""
    ),
    2: _(
        """
 post-traitement numéro :  %(i1)d
 l'instant demandé n'a pas été calculé
 pas de post-traitement
 """
    ),
    3: _(
        """
 post-traitement numéro :  %(i1)d
 aucune maille ne correspond aux critères demandés
 pas de post-traitement
"""
    ),
    4: _(
        """
Les données d'abscisses curvilignes ne sont pas disponibles pour la maille %(k1)s.
Si vous souhaitez les calculer il faut intégrer cette maille lors de l'appel à MODI_MAILLAGE/ABSC_CURV.
"""
    ),
    5: _(
        """
 il manque le vecteur des composantes
Ce message est un message d'erreur développeur.
Contactez le support technique.
"""
    ),
    6: _(
        """
 chemin nul ou défini en un noeud
"""
    ),
    7: _(
        """
 le nombre de composantes à traiter est limité à 6 pour opération "MOYENNE".
 utiliser "NOM_CMP" avec au plus 6 composantes.
"""
    ),
    8: _(
        """
 initialisation de la table %(k1)s
"""
    ),
    9: _(
        """
 pas de champ trouvé pour l'option %(k1)s
"""
    ),
    10: _(
        """
 paramètre %(k1)s de type %(k2)s
"""
    ),
    11: _(
        """
 on ne traite que les champs complexes
"""
    ),
    12: _(
        """
 tableau de travail limité, réduire le nombre de composantes à traiter
"""
    ),
    13: _(
        """
 plus de 3000 composantes.
 Contactez le support
"""
    ),
    14: _(
        """
 en repère local, on ne traite pas le champ %(k1)s
"""
    ),
    15: _(
        """
Ce message est un message d'erreur développeur.
Contactez le support technique.
"""
    ),
    16: _(
        """
Les données d'abscisses curvilignes ne sont pas présentes dans le maillage.
Pour les créer, il faut utiliser MODI_MAILLAGE/ABSC_CURV.
"""
    ),
    17: _(
        """
 on ne traite que des champs de type "DEPL_R" pour un changement de repère
"""
    ),
    18: _(
        """
Le champ d'intérêt n'a de valeurs pour aucun des noeuds renseignés.
Vérifiez qu'ils appartiennent au modèle.
"""
    ),
    19: _(
        """
Le champ d'intérêt n'a de valeurs pour aucune des mailles renseignées.
Vérifiez qu'elles appartiennent au modèle.
"""
    ),
    20: _(
        """
 on ne traite pas ce cas
Ce message est un message d'erreur développeur.
Contactez le support technique.
"""
    ),
    21: _(
        """
 avec VECT_Y, le groupe de noeuds doit contenir au moins 2 noeuds
"""
    ),
    22: _(
        """
 avec VECT_Y, il faut préciser
   - soit un seul groupe de noeuds
   - soit plusieurs noeuds
"""
    ),
    23: _(
        """
Le champ fourni doit être issu de l'option EFGE_ELNO et non %(k1)s.
"""
    ),
    24: _(
        """
Les opérations 'MOYENNE', 'EXTREMA' et 'MOYENNE_ARITH' ne sont pas autorisées
en parallélisme distribué.
"""
    ),
    26: _(
        """
 changement de repère:
 champ non traité %(k1)s
 option de calcul %(k2)s
"""
    ),
    27: _(
        """
 noeud sur l'AXE_Z
 maille      : %(k1)s
 noeud       : %(k2)s
 coordonnées : %(r1)f
"""
    ),
    28: _(
        """
 les noeuds du maillage ne sont pas tous dans un même plan Z = CONSTANTE
 changement de repère non traité
"""
    ),
    30: _(
        """
 le noeud %(k1)s est confondu avec l'origine
"""
    ),
    31: _(
        """
 le noeud %(k1)s est sur l'AXE_Z
"""
    ),
    32: _(
        """
 les noeuds du maillage ne sont pas tous dans un même plan Z = CSTE
 option TRAC_NOR non traitée
 utiliser l'option TRAC_DIR
"""
    ),
    33: _(
        """
 option non traitée: %(k1)s, post-traitement: %(i1)d
 les invariants tensoriels sont calculés
   pour les options :  %(k2)s
                       %(k3)s
                       %(k4)s
                       %(k5)s
                       %(k6)s
                       %(k7)s
                       %(k8)s
                       %(k9)s
                       %(k10)s
                       %(k11)s
                       %(k12)s
                       %(k13)s
                       %(k14)s
                       %(k15)s
                       %(k16)s
                       %(k17)s
                       %(k18)s
                       %(k19)s
                       %(k20)s
"""
    ),
    34: _(
        """
 option non traitée: %(k1)s, post-traitement: %(i1)d
 les traces normales sont calculées
   pour les options :  %(k2)s
                       %(k3)s
                       %(k4)s
                       %(k5)s
                       %(k6)s
                       %(k7)s
                       %(k8)s
                       %(k9)s
                       %(k10)s
                       %(k11)s
                       %(k12)s
                       %(k13)s
                       %(k14)s
                       %(k15)s
                       %(k16)s
                       %(k17)s
                       %(k18)s
   ou pour les grandeurs %(k19)s
                         %(k20)s
"""
    ),
    35: _(
        """
 option non traitée: %(k1)s, post-traitement: %(i1)d
 les traces directionnelles sont calculées
   pour les options :  %(k2)s
                       %(k3)s
                       %(k4)s
                       %(k5)s
                       %(k6)s
                       %(k7)s
                       %(k8)s
                       %(k9)s
                       %(k10)s
                       %(k11)s
                       %(k12)s
                       %(k13)s
                       %(k14)s
                       %(k15)s
                       %(k16)s
                       %(k17)s
                       %(k18)s
                       %(k19)s
                       %(k20)s
   ou pour les grandeurs %(k21)s
                         %(k22)s
                         %(k23)s
"""
    ),
    36: _(
        """
 trace directionnelle, post-traitement: %(i1)d
 direction nulle, pas de calcul
"""
    ),
    37: _(
        """
 attention post-traitement %(i1)d
 seules les composantes du tenseur des contraintes sont traitées
"""
    ),
    38: _(
        """
 Post-traitement %(i1)d
 Composante non traitée dans un changement de repère
 Contactez le support
"""
    ),
    39: _(
        """
 Post-traitement %(i1)d
 Grandeur %(k1)s non traitée dans un changement de repère
 Les changements de repère sont possibles
   pour la grandeur %(k2)s  option: %(k3)s
   pour la grandeur %(k4)s  options: %(k5)s %(k6)s
   pour les grandeurs %(k7)s  %(k8)s
"""
    ),
    40: _(
        """
 Le noeud numéro %(i1)d n'est pas connecté à la maille de nom %(k1)s
"""
    ),
    41: _(
        """
 Champ inexistant NOM_CHAM: %(k1)s  NUME_ORDRE: %(i1)d
"""
    ),
    42: _(
        """
 Occurrence %(i1)d du mot clé facteur ACTION :
 Les listes arguments des mots clés RESULTANTE et MOMENT doivent être de même longueur
 cette longueur doit être de 2 ou 3
"""
    ),
    43: _(
        """
 Occurrence %(i1)d du mot clé facteur ACTION :
 La liste arguments du mot clé POINT doit être de longueur 2 ou 3
"""
    ),
    44: _(
        """
 Occurrence %(i1)d du mot clé facteur ACTION :
 On ne peut accéder au RESULTAT de nom %(k1)s et de type %(k2)s par %(k3)s ou par %(k4)s
"""
    ),
    45: _(
        """
 Occurrence %(i1)d du mot clé facteur ACTION :
 On ne peut accéder au RESULTAT de nom %(k1)s et de type %(k2)s par %(k3)s
"""
    ),
    46: _(
        """
 Occurrence %(i1)d du mot clé facteur ACTION :
 Le NOM_CHAM %(k1)s n'est pas autorisé pour le RESULTAT %(k2)s de type %(k3)s
 ou le NOM_CHAM est autorisé mais aucun champ effectif n'existe.
"""
    ),
    47: _(
        """
 Occurrence %(i1)d du mot clé facteur ACTION :
 Le ou les champs élémentaires mis en jeu est ou sont donnés aux points de gauss
 ce ou ces champs ne sont pas traités.
"""
    ),
    48: _(
        """
 Occurrence %(i1)d du mot clé facteur ACTION :
 La composante %(k1)s n'est pas présente au catalogue des grandeurs.
"""
    ),
    49: _(
        """
 L'instant %(r1)f ne correspond à aucun numéro d'ordre du résultat fourni.
"""
    ),
    50: _(
        """
 Occurrence %(i1)d du mot clé facteur ACTION :
 Le groupe de noeuds %(k1)s ne fait pas partie du maillage sous-jacent au champ à traiter.
"""
    ),
    51: _(
        """
 Occurrence %(i1)d du mot clé facteur ACTION :
 Le noeud %(k1)s ne fait pas partie du maillage sous-jacent au champ à traiter.
"""
    ),
    52: _(
        """
 On ne traite pas le FORMAT_C %(k1)s
"""
    ),
    53: _(
        """
Ce message est un message d'erreur développeur.
Contactez le support technique.
"""
    ),
    54: _(
        """
 Occurrence %(i1)d du mot clé facteur ACTION :
 Impossible de récupérer les composantes du champ.
"""
    ),
    55: _(
        """
 la composante %(k1)s est en double.
"""
    ),
    56: _(
        """
 la composante %(k1)s n'est pas une composante de %(k2)s
"""
    ),
    57: _(
        """
 la grandeur %(k1)s est inconnue au catalogue.
"""
    ),
    58: _(
        """
 Occurrence %(i1)d du mot clé facteur ACTION :
 Le groupe de mailles %(k1)s ne fait pas partie du maillage sous-jacent au champ à traiter.
"""
    ),
    59: _(
        """
 Le contenu de la table n'est pas celui attendu !
 Contactez le support
"""
    ),
    60: _(
        """
 arrêt sur erreurs
 Contactez le support
"""
    ),
    61: _(
        """
 Nombre de cycles admissibles négatif,
 vérifiez la courbe de WOHLER
 contrainte calculée =  %(r1)f   Nombre admissible =  %(r2)f

"""
    ),
    62: _(
        """
 Attention lors de la définition de votre liste de noeuds,
 %(i1)d noeud est hors de la matière

"""
    ),
    63: _(
        """
 Attention lors de la définition de votre liste de noeuds,
 %(i1)d noeuds sont hors de la matière

"""
    ),
    64: _(
        """
Vous avez demandé à extraire les valeurs d'un champ sur des noeuds.
Or, ces noeuds ne sont pas présents dans le modèle qui a servi à calculer le champ.

Conseil :
    - Vérifiez qu'une modélisation a bien été affectée sur ces noeuds.
    - Si vous utilisez PROPA_FISS, peut-être s'agit-il de noeuds créés par cet opérateur.
"""
    ),
    65: _(
        """
La composante %(k1)s n'existe pas pour le champ de type %(k2)s
du résultat %(k3)s
"""
    ),
    66: _(
        """
Dans le cas d'un champ de type ELEM, l'utilisation des mots clés NOEUD, GROUP_NO
n'a pas de sens, seul le mot clé MAILLE ou GROUP_MA est autorisé.
"""
    ),
    67: _(
        """
Erreur utilisateur dans la commande POST_RELEVE_T :
   Le mot clé RESULTANTE n'est autorisé que pour OPERATION='EXTRACTION'.
"""
    ),
    68: _(
        """
Erreur utilisateur dans la commande POST_RELEVE_T :
  Pour les mots clés ACTION / RESULTANTE (et MOMENT), on ne peut utiliser que le repère
 'GLOBAL' ou le repère 'UTILISATEUR' (avec ANGL_NAUT).
"""
    ),
    70: _(
        """
Commande POST_RELEVE_T :
Votre modèle contient des éléments de structures et :
    - vous avez demandé le post-traitement du champ %(k1)s, qui n'a de signification que dans
      le repère LOCAL des éléments.
    - vous avez demandé le post-traitement dans un repère différent de LOCAL.

Le mot clef REPERE est ignoré.
Le post-traitement sera réalisé dans le repère LOCAL des éléments.

Conseil :
    Si vous avez précédemment utilisé la commande CALC_CHAMP et que les repères locaux des
    éléments connectés à un même noeud sont différents, une alarme a été émise.
    Dans ce cas les résultats obtenus aux noeuds ne sont pas pertinents.
"""
    ),
}
