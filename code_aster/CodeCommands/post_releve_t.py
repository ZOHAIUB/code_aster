# coding: utf-8

# Copyright (C) 1991 - 2023  EDF R&D                www.code-aster.org
#
# This file is part of Code_Aster.
#
# Code_Aster is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Code_Aster is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Code_Aster.  If not, see <http://www.gnu.org/licenses/>.

# person_in_charge: nicolas.sellenet@edf.fr

from ..Objects import Table
from ..Supervis import ExecuteCommand
from ..Messages import UTMESS


class PostReleveT(ExecuteCommand):
    """Command that defines :class:`~code_aster.Objects.Table`."""

    command_name = "POST_RELEVE_T"

    def create_result(self, keywords):
        """Initialize the result.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """
        self._result = Table()

    def exec_(self, keywords):
        """Execute the command.

        Arguments:
            keywords (dict): User's keywords.
        """
        for i, act in enumerate(keywords["ACTION"]):
            mesh = (act.get("RESULTAT") or act.get("CHAM_GD")).getMesh()
            if not mesh or not mesh.isParallel():
                break
            for grp in act.get("GROUP_MA", []):
                if not mesh.hasGroupOfCells(grp):
                    UTMESS("F", "POSTRELE_58", valk=grp, vali=i + 1)
            for grp in act.get("GROUP_NO", []):
                if not mesh.hasGroupOfNodes(grp):
                    UTMESS("F", "POSTRELE_50", valk=grp, vali=i + 1)

        super().exec_(keywords)

    def add_dependencies(self, keywords):
        """Do not keep any references to original objects.

        Arguments:
            keywords (dict): User's keywords.
        """

    def post_exec(self, keywords):
        """
        Arguments:
            keywords (dict): Keywords arguments of user's keywords, changed
                in place.
        """
        self._result.build()


POST_RELEVE_T = PostReleveT.run
