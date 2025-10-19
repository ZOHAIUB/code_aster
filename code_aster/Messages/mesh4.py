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
    1: _("""Les surfaces à apparier n'ont pas de cellules."""),
    # 2: _("""Les cellules de la surface ne sont pas toutes des mailles de peau."""),
    3: _(
        """
         Une erreur non fatale est détectée lors de l'appariement de deux cellules.
         Ce couple est ignoré, il y a potentiellement interpénétration entre ces deux cellules.
         Conseil :
             Vérifiez la qualité de l'appariement visuellement .
"""
    ),
    4: _(
        """On ne trouve pas les objets de l'appariement. Vérifiez que vous avez bien effectué le calcul."""
    ),
    5: _("""La paire demandée n'existe pas à cet index."""),
}
