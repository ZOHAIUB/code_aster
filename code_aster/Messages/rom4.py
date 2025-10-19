# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

# person_in_charge: mickael.abbas at edf.fr

from ..Utilities import _

# Messages for DEFI_DOMAINE_REDUIT

cata_msg = {
    7: _("""Échec lors du calcul des modes pour l'estimation du domaine réduit."""),
    9: _("""Les bases ne sont pas définies sur le même maillage."""),
    10: _("""Les bases ne sont pas définies sur le même maillage que le résultat."""),
    11: _(
        """On ne peut pas utiliser une base avec des modes définis sur un autre support que les noeuds."""
    ),
    12: _("""Le GROUP_NO %(k1)s fait déjà partie du maillage."""),
    13: _("""Le GROUP_MA %(k1)s fait déjà partie du maillage."""),
    15: _(
        """Le nombre de couches du domaine est inférieur au nombre de couches de l'interface, on aura peut-être un problème de convergence"""
    ),
    16: _("""Les modes ne sont pas des champs du type attendu."""),
    17: _("""Les modes ne sont pas des champs gradients du type attendu (on attend %(k1)s)."""),
    20: _("""Calcul du domaine réduit. """),
    21: _("""Création des groupes dans le maillage pour l'estimation du domaine réduit."""),
    22: _("""Nombre d'éléments dans le domaine réduit: %(i1)d"""),
    23: _("""Nombre de noeuds sur l'interface du domaine réduit: %(i1)d"""),
    25: _("""Nombre total de points magiques pour l'estimation du domaine réduit: %(i1)d"""),
    26: _("""Nombre de noeuds dans la zone rigide pour le correcteur éléments finis : %(i1)d"""),
    27: _("""Il n'y a aucun noeud dans la zone d'interface."""),
    28: _("""Il n'y a aucun noeud dans la zone rigide pour le correcteur éléments finis."""),
    29: _("""Nombre de noeuds dans le domaine réduit: %(i1)d"""),
}
