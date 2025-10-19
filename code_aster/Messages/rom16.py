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

# Messages for REST_REDUIT_COMPLET

cata_msg = {
    1: _("""Lecture des paramètres de la commande."""),
    2: _("""Initialisations pour le post-traitement du calcul réduit."""),
    3: _("""Reconstruction des champs sur tout le domaine."""),
    4: _("""Initialisations pour la reconstruction des %(i1)d champs."""),
    5: _("""Copie des paramètres de la structure de données."""),
    20: _("""Les bases ne sont pas définies sur le même maillage."""),
    21: _("""Les bases ne sont pas définies sur le même modèle."""),
    22: _(
        """Le modèle d'une des bases n'est pas celui du modèle complet. Vérifiez que vous n'utilisez pas la base tronquée."""
    ),
    23: _("""Le modèle est le même pour la reconstruction que le modèle réduit d'origine."""),
    24: _(
        """Vous demandez à calculer un champ de type %(k1)s par REST_REDUIT_COMPLET alors que ce champ n'existe pas dans le résultat réduit./
Conseil: utilisez CALC_CHAMP pour calculer ce champ."""
    ),
    25: _("""Le résultat réduit n'est pas défini sur le même maillage que le résultat complet."""),
    26: _(
        """Vous demandez à calculer un champ avec un support de type %(k1)s en réduisant le domaine. Ce n'est pas autorisé."""
    ),
    30: _("""Liste des champs initialement présents dans les résultats réduits: """),
    31: _(""" Type du champ: %(k1)s. Ce champ sera reconstruit."""),
    32: _(""" Type du champ: %(k1)s. Ce champ ne sera pas reconstruit."""),
    50: _("""Le résultat sur le modèle complet sera de type %(k1)s."""),
    51: _("""Le résultat sur le modèle réduit contient %(i1)d numéros d'ordre."""),
}
