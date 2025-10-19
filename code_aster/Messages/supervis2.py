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
    2: _(
        """
    Vous utilisez une vieille version de Code_Aster.

    En mettant à jour votre version, vous bénéficierez des dernières améliorations
    apportées au code depuis 15 mois.
    Si vous avez des développements privés, vous risquez d'avoir un travail
    important de portage si vous ne suivez pas les mises à jour.
"""
    ),
    3: _(
        """
Erreur inattendue lors de l'exécution de la commande '%(k1)s'.
Merci de signaler cette anomalie.

Erreur :

%(k2)s
%(k3)s
"""
    ),
    5: _(
        """
Erreur lors du chargement du catalogue du matériau '%(k1)s'.
"""
    ),
    6: _(
        """
Erreur dans le catalogue du matériau '%(k1)s'.

Il n'est pas possible d'avoir plusieurs occurrences pour le
mot-clé facteur '%(k2)s'.
"""
    ),
    # Lignes d'entete
    7: {
        "message": _(
            """

-- CODE_ASTER -- VERSION : %(k1)s --
"""
        ),
        "flags": "CENTER",
    },
    8: {
        "message": _(
            """Version %(k1)s modifiée le %(k2)s
révision %(k3)s - branche '%(k4)s'"""
        ),
        "flags": "CENTER",
    },
    9: {
        "message": _(
            """Copyright EDF R&D %(k1)s - %(k2)s

Exécution du : %(k3)s
Nom de la machine : %(k4)s
Architecture : %(k5)s
Type de processeur : %(k6)s
Système d'exploitation : %(k7)s
Langue des messages : %(k8)s
"""
        ),
        "flags": "CENTER",
    },
    10: {
        "message": _(
            """Version de Python : %(k1)s
Version de NumPy : %(k2)s"""
        ),
        "flags": "CENTER",
    },
    # fin Lignes d'entete
    11: {
        "message": _(
            """Parallélisme MPI : actif
Rang du processeur courant : %(i1)d
Nombre de processeurs MPI utilisés : %(i2)d"""
        ),
        "flags": "CENTER",
    },
    12: {"message": _("""Parallélisme MPI : inactif"""), "flags": "CENTER"},
    13: {
        "message": _(
            """Parallélisme OpenMP : actif
Nombre de processus OpenMP utilisés : %(i1)d"""
        ),
        "flags": "CENTER",
    },
    14: {
        "message": _("""Version de la librairie HDF5 : %(i1)d.%(i2)d.%(i3)d"""),
        "flags": "CENTER",
    },
    15: {"message": _("""Librairie HDF5 : non disponible"""), "flags": "CENTER"},
    16: {"message": _("""Version de la librairie MED : %(i1)d.%(i2)d.%(i3)d"""), "flags": "CENTER"},
    17: {"message": _("""Librairie MED : non disponible"""), "flags": "CENTER"},
    18: {"message": _("""Version de la librairie MUMPS : %(k1)s"""), "flags": "CENTER"},
    19: {"message": _("""Librairie MUMPS : non disponible"""), "flags": "CENTER"},
    20: {
        "message": _("""Version de la librairie SCOTCH : %(i1)d.%(i2)d.%(i3)d"""),
        "flags": "CENTER",
    },
    21: {"message": _("""Librairie SCOTCH : non disponible"""), "flags": "CENTER"},
    22: {
        "message": _(
            """Mémoire limite pour l'allocation dynamique : %(r1).2f Mo
ajouté à l'initialisation du processus : %(r3).2f Mo
Limite cible du processus : %(r2).2f Mo"""
        ),
        "flags": "CENTER",
    },
    24: {
        "message": _(
            """Taille limite des fichiers d'échange : %(r1).2f Go
"""
        ),
        "flags": "CENTER",
    },
    25: {"message": _("""Version de la librairie PETSc : %(k1)s"""), "flags": "CENTER"},
    26: {"message": _("""Librairie PETSc : non disponible"""), "flags": "CENTER"},
    27: {"message": _("""Version de la librairie MFront : %(k1)s"""), "flags": "CENTER"},
    28: {"message": _("""Librairie MFront : non disponible"""), "flags": "CENTER"},
    29: {
        "message": _(
            """Mémoire limite pour l'allocation dynamique : %(r1).2f Mo
Limite cible du processus : %(r2).2f Mo"""
        ),
        "flags": "CENTER",
    },
    # Affichage des commandes
    69: _("""Le résultat intermédiaire '%(k1)s' sera disponible sous '%(k2)s'."""),
    70: """# ----------------------------------------------------------------------------------------------""",
    71: _("""# Commande #%(i1)04d de %(k1)s, ligne %(i2)d"""),
    72: _("""# Résultat commande #%(i1)04d (%(k1)s): %(k2)s %(k3)s"""),
    73: _(
        """# Mémoire (Mo) : %(r1)8.2f / %(r2)8.2f / %(r3)8.2f / %(r4)8.2f (VmPeak / VmSize / Optimum / Minimum)"""
    ),
    # attention au décalage
    74: _("""# Mémoire (Mo) : %(r3)8.2f / %(r4)8.2f (Optimum / Minimum)"""),
    75: _(
        """# Fin commande #%(i1)04d   user+syst:%(r1)12.2fs (syst:%(r2)12.2fs, elaps:%(r3)12.2fs)"""
    ),
    76: _(""" de type <%(k1)s>"""),
    # sans formatage pour STAT_NON_LINE (impmem)
    77: _(
        """
  Mémoire (Mo) : %(r1)8.2f / %(r2)8.2f / %(r3)8.2f / %(r4)8.2f (VmPeak / VmSize / Optimum / Minimum)
"""
    ),
    # attention au décalage
    78: _(
        """
  Mémoire (Mo) : %(r3)8.2f / %(r4)8.2f (Optimum / Minimum)
"""
    ),
    79: _(
        """
La commande est ré-entrante. Il faut renseigner le mot-clé '%(k1)s' pour
indiquer l'objet qui est modifié.
"""
    ),
    80: _("""# Dépend de :"""),
    98: _(
        """
 <INFO> Démarrage de l'exécution.
"""
    ),
    99: _("""<I> <FIN> ARRET NORMAL DANS "FIN" PAR APPEL A "JEFINI"."""),
}
