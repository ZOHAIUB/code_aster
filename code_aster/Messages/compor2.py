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
        """Le type de DEFORMATION <%(k1)s>, est actuellement incompatible avec SUPPORT=POINT. On utilise donc SUPPORT=ELEMENT."""
    ),
    2: _("""On ne peut avoir à la fois SIGM et EPSI imposés sur la composante <%(k1)s>."""),
    3: _("""On doit avoir une seule composante donnée parmi <%(k1)s>."""),
    4: _(
        """Problème lors de l'inversion de la matrice jacobienne. On tente de subdiviser le pas de temps."""
    ),
    5: _(
        """Le nombre maximum d'itérations a été atteint. On tente de subdiviser le pas de temps."""
    ),
    6: _("""Le nombre de phases est trop grand."""),
    7: _("""On ne peut pas utiliser les KIT_THM dans SIMU_POINT_MAT."""),
    8: _(
        """La somme des fractions volumiques est très différente de 1.0, elle vaut <%(r1).15E>.
Vérifiez FRAC_VOL pour toutes les occurrences du mot clé POLYCRISTAL."""
    ),
    10: _(
        """
Le redécoupage local du pas de temps n'est pas compatible avec <%(k1)s>
"""
    ),
    11: _(
        """
La rotation de réseau n'est pas compatible avec RUNGE_KUTTA. Utiliser l'intégration IMPLICITE.
"""
    ),
    12: _(
        """
   On ne peut pas utiliser les variables de commandes intrinsèques dans SIMU_POINT_MAT en mode POINT.
   Votre loi de comportement en utilise. Par exemple, la taille de l'élément. Le calcul est donc impossible. Il faut demander une évolution.
"""
    ),
    13: _(
        """
  LA MODELISATION GRAD_SIGM N'EST COMPATIBLE QU'AVEC LA LOI ENDO_HETEROGENE.
"""
    ),
    14: _(
        """
 ENDO_HETEROGENE : Les critères entre KI et SY ne sont pas respectés ; baissez KI ou augmentez SY
"""
    ),
    15: _(
        """
 MONOCRISTAL : la matrice d'interaction fournie n'est pas carrée : nombre lignes = <%(r1)E>, nombre colonnes = <%(r2)E>.
"""
    ),
    16: _(
        """
 POLYCRISTAL : Il faut au maximum 5 mono cristaux différents sur l'ensemble des phases.  Ici,il y en a : <%(i1)i>.
"""
    ),
    17: _(
        """
 MONOCRISTAL : la matrice d'interaction fournie ne comporte pas le bon nombre de systèmes.
 il en faut : <%(i1)i>.
"""
    ),
    18: _(
        """
 MONOCRISTAL : la matrice d'interaction fournie n'est pas symétrique.
"""
    ),
    19: _(
        """
 MONOCRISTAL : le nombre de composantes de n et m n'est pas correct :  <%(r1)E> au lieu de 6.
"""
    ),
    20: _(
        """
 MONOCRISTAL : comme il y a  plusieurs familles de systèmes de glissement,
 il faut fournir une matrice d'interaction entre tous ces systèmes, de dimension  <%(i1)i>
"""
    ),
    21: _(
        """
 MONOCRISTAL : pas de matrice jacobienne programmée actuellement pour  MONO_DD_FAT.
 Utiliser ALGO_INTE='NEWTON_PERT'
"""
    ),
    22: _(
        """
   SIMU_POINT_MAT : Le type de DEFORMATION choisi,  <%(k1)s>, est incompatible avec GRAD_IMPOSE.
   GRAD_IMPOSE n'est utilisable qu'avec DEFORMATION='SIMO_MIEHE'.
"""
    ),
    23: _(
        """
   Il y a incohérence entre le champ des variables internes et le comportement affecté sur les mailles :
   si vous ne l'avez pas précisé (sans reuse), on suppose par défaut que le comportement est élastique.

   On ne peut pas modifier le champ des variables internes de manière automatique.
   Mais comme vous semblez être dans un cas autorisé (voir message précédent), la reprise de ce champ
   dans un calcul non-linéaire devrait bien se passer (correction automatique dans l'opérateur).
"""
    ),
    24: _(
        """
Erreur lors de la vérification de la cohérence entre les champs de variables internes.
Le maillage sur lequel s'appuie le modèle et le maillage du champ des variables internes ne sont pas les mêmes.
"""
    ),
    25: _(
        """On ne peut pas utiliser la variable de commande intrinsèque <%(k1)s> dans MFront car ce n'est pas un scalaire."""
    ),
    26: _(
        """La variable de commande correspondant au gradient de vélocité ne peut pas être utilisée pour ce comportement."""
    ),
    31: _(
        """
  CALC_ESSAI_GEOMECA : Erreur dans la saisie du mot clef facteur <%(k1)s> (occurrence %(i1)d).
  Vous n'avez pas renseigné le même nombre d'éléments pour les mots clefs simple %(k2)s.
"""
    ),
    32: _(
        """
  CALC_ESSAI_GEOMECA : Erreur dans la saisie du mot clef facteur <%(k1)s> (occurrence %(i1)d).
  La valeur <%(r1)E> a été renseignée pour le mot clef simple <%(k2)s>, au lieu d'une valeur strictement positive.
"""
    ),
    33: _(
        """
  CALC_ESSAI_GEOMECA : Erreur dans la saisie du mot clef facteur <%(k1)s> (occurrence %(i1)d).
  La valeur du mot clef simple <BIOT_COEF> doit être comprise dans l'intervalle ]0,1].
  Or vous avez renseigné la valeur <%(r1)E>
"""
    ),
    34: _(
        """
  CALC_ESSAI_GEOMECA : Erreur dans la saisie du mot clef facteur <%(k1)s> (occurrence %(i1)d).
  La valeur <%(r1)E> a été renseignée pour le mot clef simple <%(k2)s>, au lieu d'une valeur strictement positive.
"""
    ),
    35: _(
        """
  CALC_ESSAI_GEOMECA : Erreur dans la saisie du mot clef facteur <%(k1)s> (occurrence %(i1)d).
  Il y a n = %(i2)d occurrences du mot clef simple <PRES_CONF>, il faut donc renseigner n+1 = %(i3)d
  éléments pour le mot clef simple <TABLE_RESU>. Or vous en avez renseigné %(i4)d
"""
    ),
    36: _(
        """
  CALC_ESSAI_GEOMECA : Pour l'essai <%(k1)s>.
  Erreur lors du calcul du module de %(k2)s sécant maximal : pour les valeurs de paramètres matériau que vous
  avez choisies, la valeur par défaut du mot clef simple <%(k3)s> conduit à sortir du domaine d'élasticité du matériau.
  Il faut donc renseigner une valeur strictement inférieure à  <%(r1)E> pour <%(k3)s>
"""
    ),
    38: _(
        """
  CALC_ESSAI_GEOMECA : Erreur dans la saisie du mot clef facteur <%(k1)s> (occurrence %(i1)d).
  La liste de valeurs renseignées pour le mot clef simple <%(k2)s> doit être %(k4)s .
  Or vous avez renseigné la liste suivante :
  %(k3)s
"""
    ),
    39: _(
        """SIMU_POINT_MAT utilise une loi de comportement avec le gradient de vélocité comme paramètre. Ce paramètre a été initialisé à zéro car il n'est pas calculable en l'état dans la commande.
Assurez-vous que cette valeur par défaut correspond bien au cas que vous voulez traiter."""
    ),
    40: _(
        """
  CALC_ESSAI_GEOMECA : %(k1)s numéro %(k2)s.

  Pour le chargement :
    %(k3)s
    %(k4)s
  avec un nombre total de cycles :
    %(k5)s

  Le calcul a été mené jusqu'à %(i1)d cycles.

  Le critère de liquéfaction n'a pas été atteint.
"""
    ),
    41: _(
        """
  CALC_ESSAI_GEOMECA : %(k1)s numéro %(k2)s.

  Pour le chargement :
    %(k3)s
    %(k4)s
  avec un nombre total de cycles :
    %(k5)s

  Le calcul s'est arrêté au cours du cycle numéro %(i1)d pour cause d'absence de convergence,
  les pas archivés avant l'arrêt ont été sauvegardés.

  Le critère de liquéfaction a été atteint au cours du cycle numéro %(i2)d.
"""
    ),
    42: _(
        """
  CALC_ESSAI_GEOMECA : %(k1)s numéro %(k2)s.

  Pour le chargement :
    %(k3)s
    %(k4)s
  avec un nombre total de cycles :
    %(k5)s

  Le calcul s'est arrêté au cours du cycle numéro %(i1)d pour cause d'absence de convergence,
  les pas archivés avant l'arrêt ont été sauvegardés.

  Le critère de liquéfaction n'a pas été atteint.
"""
    ),
    43: _(
        """
  CALC_ESSAI_GEOMECA : Erreur dans la saisie du mot clef simple <TABLE_REF> pour la table <%(k1)s>.
  Vous avez indiqué dans cette table le TYPE suivant :
    <'%(k2)s'>
  Or on ne peut indiquer qu'un TYPE figurant parmi la liste de GRAPHIQUE demandée en sortie pour cet essai:
    <%(k3)s>
"""
    ),
    44: _(
        """
  CALC_ESSAI_GEOMECA : Erreur dans la saisie du mot clef simple <TABLE_REF> pour la table <%(k1)s>.
  La liste des paramètres d'une TABLE_REF doit nécessairement être <'TYPE','LEGENDE','ABSCISSE','ORDONNEE'>
  Or liste des paramètres de la table renseignée est :
  <%(k2)s>
"""
    ),
    45: _(
        """
  CALC_ESSAI_GEOMECA : Erreur dans la saisie du mot clef simple <TABLE_REF> pour la table <%(k1)s>.
  La colonne <%(k2)s> d'une TABLE_REF ne doit contenir qu'un élément de type chaîne de caractères.
"""
    ),
    46: _(
        """
  CALC_ESSAI_GEOMECA : Erreur dans la saisie du mot clef simple <TABLE_REF> pour la table <%(k1)s>.
  Les colonnes ABSCISSE et ORDONNEE d'une TABLE_REF doivent avoir même cardinal et contenir des réels.
"""
    ),
    47: _(
        """
  CALC_ESSAI_GEOMECA : Erreur dans la saisie du mot clef facteur <%(k1)s> (occurrence %(i1)d).
  Incohérence entre les valeurs saisies pour les mot clef simples <PRES_CONF>, <SIGM_IMPO> et <SIGM_décharge>.
  On doit toujours avoir PRES_CONF + SIGM_IMPO <= SIGM_décharge .
  Or vous avez renseigné < PRES_CONF = %(r1)E> et <SIGM_IMPO = %(r2)E>, soit PRES_CONF + SIGM_IMPO = %(r3)E supérieur à %(r4)E
"""
    ),
    52: _(
        """Maille: %(k1)-8s - Nombre de points précédent: %(i1)d - Nombre de points courant: %(i2)d"""
    ),
    54: _(
        """
Erreur lors de la vérification de la cohérence entre les champs de contrainte.

On a affiché ci-dessus la liste des mailles pour lesquelles le nombre de sous-points est différent.
"""
    ),
    91: _(
        """
   SIMU_POINT_MAT : META_LEMA_ANI est interdit à cause de la matrice de Hill définie en repère cylindrique.
   Si la matrice de Hill est isotrope, utiliser SIMU_POINT_MAT avec SUPPORT='ELEMENT'.
   Sinon, utiliser STAT_NON_LINE sur un maillage adéquat.
"""
    ),
    92: _(
        """
   SIMU_POINT_MAT : META_LEMA_ANI définit une matrice de Hill définie en repère cylindrique.
   Il faut donc vérifier que  la matrice de Hill est isotrope.
   Sinon, utiliser STAT_NON_LINE sur un maillage adéquat.
"""
    ),
    93: _(
        """
   Quand on utilise PETIT_REAC dans SIMU_POINT_MAT, il n'est pas possible d'imposer exactement une déformation donnée.
   Les déformations à la fin du calcul seront différentes de celles que l'on a donné, sauf en petites déformations.
   Si vous voulez simuler des grandes déformations, utilisez plutôt GDEF_LOG.
"""
    ),
    94: _(
        """
Il faut fournir le mot-clé BETON_DESORP à DEFI_MATERIAU pour calculer l'hygrométrie.
"""
    ),
    95: _(
        """
On ne trouve pas la variable de commande %(k1)s qui est nécessaire pour calculer
l'hygrométrie par isotherme de LEVERETT.
Il faut renseigner cette variable de commande dans le mot-clé AFFE_VARC d'AFFE_MATERIAU.
"""
    ),
    96: _(
        """
Pour le calcul de l'hygrométrie par isotherme de LEVERETT,
la valeur de la saturation en eau doit être comprise entre 0 et 1.

La saturation est obtenue en divisant le séchage (ou concentration en eau)
par la porosité et par le facteur 1000.

Teneur en eau : %(r1)E
"""
    ),
    97: _(
        """
Dans le cas du calcul de l'hygrométrie par isotherme de LEVERETT,
la température doit être inférieure à 300 degrés.
"""
    ),
    98: _(
        """
Hygrométrie par isotherme de LEVERETT :
La valeur de la saturation en eau est trop proche des bornes limites (0. et 1.).
On a translaté d'un epsilon cette valeur vers l'intérieur de l'intervalle.

La teneur ou saturation en eau est obtenue en divisant le séchage (concentration en eau)
par la porosité et par le facteur 1000.
"""
    ),
}
