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

from ..Utilities import _

cata_msg = {
    1: _(
        """
 Impossible d'importer le catalogue de la structure de données '%(k1)s'
"""
    ),
    4: _(
        """
 Arguments incohérents :
      Nom des paramètres : %(k1)s
   Valeur des paramètres : %(k2)s
"""
    ),
    5: _(
        """
 Arguments invalide
Ce message est un message d'erreur développeur.
Contactez le support technique.
"""
    ),
}
