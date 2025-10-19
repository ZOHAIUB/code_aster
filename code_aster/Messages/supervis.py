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
    1: _("""<INFO> Activation du mode parallélisme distribué."""),
    2: _(
        """
Les commandes DEBUT et POURSUITE doivent être appelées une fois et une seule.
"""
    ),
    3: _(
        """
  Erreur programmeur : %(k1)s non appariés.
"""
    ),
    4: _(
        """
Erreur de syntaxe dans %(k1)s

%(k2)s

Exception détaillée ci-dessous.
"""
    ),
    5: _("""Suppression de la référence : '%(k1)s'"""),
    6: _(
        """
Il n'y a pas d'objet référencé avec cette clé: '%(k1)s'.
    """
    ),
    7: _(
        """
La commande %(k1)s n'est pas disponible en mode parallélisme distribué.
    """
    ),
    8: _(
        """
  Un nom de concept intermédiaire doit commencer par '.' ou '_' et non :  %(k1)s
"""
    ),
    9: _(
        """
Fonctionnalité obsolète : %(k1)s

La fonctionnalité mentionnée ci-dessus est obsolète et son utilisation
est déconseillée dans la version du code que vous utilisez.
Elle sera supprimée dans la version %(i1)d.

Conseils :
- pour assurer la pérennité de votre étude, modifiez votre fichier de commandes ;
- si vous avez un besoin impérieux de cette fonctionnalité,
  rapprochez-vous de votre correspondant Club utilisateur ou de l'équipe de développement.

"""
    ),
    10: _(
        """
    Si le mot-clé TOUT est donné sous MASSIF, il ne peut y avoir qu'une seule occurrence du mot-clé facteur MASSIF.
"""
    ),
    11: _("""Les alarmes ne sont pas **encore** aggravées en Erreur sur cette plate-forme."""),
    12: _(
        """
  Exécution de JEVEUX en mode DEBUG
"""
    ),
    13: _(
        """
  %(k1)s  nom de base déjà définie
"""
    ),
    14: _(
        """
  %(k1)s  statut impossible pour la base globale
"""
    ),
    15: _(
        """
  Problème d'allocation des bases de données
"""
    ),
    16: _(
        """
  Écriture des catalogues des éléments faite.
"""
    ),
    17: _(
        """
  Relecture des catalogues des éléments faite.
"""
    ),
    18: _(
        """
  Trop de catalogues (maximum = 10)
"""
    ),
    20: _(
        """
  "%(k1)s" argument invalide du mot clé "FICHIER" du mot clé facteur "CATALOGUE"
"""
    ),
    21: _(
        """
  Erreur(s) fatale(s) lors de la lecture des catalogues
"""
    ),
    22: {
        "message": _(
            """
Les mots-clés CODE et DEBUG dans DEBUT/POURSUITE sont réservés aux cas-tests.
De même pour l'option, "%(k1)s" de %(k2)s.
Il ne faut pas les utiliser dans les études car ils modifient certaines valeurs par
défaut des commandes DEBUT/POURSUITE qui ont des conséquences sur le comportement
en cas d'erreur ou sur les performances.
"""
        ),
        "flags": "DECORATED",
    },
    23: _(
        """
  Débogage JXVERI demandé
"""
    ),
    24: _(
        """
  Débogage SDVERI demandé
"""
    ),
    25: _(
        """
L'usage de la commande INCLUDE n'est pas conseillé.

Une étude utilisant la commande INCLUDE ne peut pas être éditée en mode
graphique dans %(k1)s. Il est alors conseillé de remplacer INCLUDE par
un "Stage" supplémentaire.

De plus, le fonctionnement de INCLUDE est différent de ce qu'il était dans les versions
antérieures à la version 15. Il est à la fois plus simple et plus robuste.
S'il s'agit d'exécuter des instructions Python, il serait préférable d'utiliser "%(k2)s".
"""
    ),
    31: _(
        """
 Valeur invalide pour le mot clé RESERVE_CPU
"""
    ),
    38: _(
        """
 Il n'y a plus de temps pour continuer
"""
    ),
    39: _(
        """
Arrêt de l'exécution suite à la réception du signal utilisateur %(k1)s.
Fermeture des bases jeveux afin de permettre la POURSUITE ultérieure du calcul.
"""
    ),
    41: {
        "message": _(
            """La version %(k1)s a été modifiée par %(i1)d révisions.
"""
        ),
        "flags": "DECORATED",
    },
    42: _(
        """Les fichiers suivants ont été modifiés par rapport à la dernière révision %(k1)s :

%(k2)s
"""
    ),
    43: _(
        """
  Débogage %(k1)s suspendu
"""
    ),
    44: _(
        """
  Débogage %(k1)s demandé
"""
    ),
    50: _(
        """
 La commande a un numéro non appelable dans cette version.
 Le numéro erroné est  %(i1)d
"""
    ),
    52: _(
        """
  Fin de lecture (durée  %(r1)f  s.) %(k1)s
"""
    ),
    56: _(
        """
  Incohérence entre le catalogue et le corps de la macro-commande.
"""
    ),
    61: _(
        """
  La commande a un numéro non appelable dans cette version
  Le numéro erroné est : %(i1)d
"""
    ),
    63: _(
        """
     ARRET PAR MANQUE DE TEMPS CPU
     Les commandes suivantes sont ignorées, on passe directement dans FIN
     La base globale est sauvegardée
     Temps consommé de la réserve CPU        :  %(r1).2f s\n
"""
    ),
    64: _(
        """
  Valeur initiale du temps CPU maximum =   %(i1)d secondes
  Valeur du temps CPU maximum passé aux commandes =   %(i2)d secondes
  Réserve CPU prévue = %(i3)d secondes
"""
    ),
    81: _(
        """
 %(k1)s nom symbolique inconnu
  - nombre de valeurs attendues %(i1)d
  - valeurs attendues : %(k1)s, %(k2)s,...
"""
    ),
    82: _(
        """
 L'argument du mot clé "CAS" est erroné.
 Valeur lue %(k1)s
 nombre de valeurs attendues %(i1)d
 valeurs attendues : %(k1)s,%(k2)s, ...
"""
    ),
    83: _(
        """

 Le nombre d'enregistrements (NMAX_ENRE) et leurs longueurs (LONG_ENRE) conduisent à un
 fichier dont la taille maximale en Mo (%(i1)d) est supérieure à limite autorisée :  %(i2)d

 Vous pouvez augmenter cette limite en utilisant l'argument "--max_base" sur la ligne
 de commande suivi d'une valeur en Mo.

"""
    ),
    96: {
        "message": _(
            """

    Réception du signal USR1. Interruption du calcul demandée...

"""
        ),
        "flags": "DECORATED",
    },
    99: _("""Une erreur s'est produite."""),
}
