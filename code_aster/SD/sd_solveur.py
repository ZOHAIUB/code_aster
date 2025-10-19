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

from . import *


class sd_solveur(AsBase):
    nomj = SDNom(fin=19)
    SLVK = AsVK24(SDNom(debut=19), lonmax=14)
    SLVR = AsVR(SDNom(debut=19), lonmax=6)
    SLVI = AsVI(SDNom(debut=19), lonmax=9)
    SLVO = AsVK80(SDNom(debut=19), lonmax=31)

    def check_SLVK(self, checker):
        # ---------------------------------------------
        slvk = self.SLVK.get_stripped()
        method = slvk[0]
        assert slvk[12] in ("OUI", "NON", "")
        assert slvk[9] in ("XXXX",)
        if method == "MUMPS":
            assert slvk[1] in ("AUTO", "SANS")
            assert slvk[2] in ("NONSYM", "SYMGEN", "SYMDEF", "AUTO")
            assert slvk[3] in (
                "AMD",
                "AMF",
                "PORD",
                "METIS",
                "QAMD",
                "AUTO",
                "SCOTCH",
                "PARMETIS",
                "PTSCOTCH",
            )
            assert slvk[4] in ("AUTO", "FR", "FR+", "FR++", "LR", "LR+", "LR++")
            assert slvk[5] in ("LAGR2", "NON")
            assert slvk[7] in ("OUI", "NON", "XXXX")
            assert slvk[8] in ("IN_CORE", "OUT_OF_CORE", "AUTO", "EVAL", "XXXX")
            assert slvk[10] in ("SANS", "AUTO", "FORCE", "XXXX", "MINI")
            assert slvk[11] in ("XXXX",)  # unused
        elif method == "MULT_FRONT":
            assert slvk[1] in ("XXXX",)
            assert slvk[2] in ("XXXX",)
            assert slvk[3] in ("MD", "MDA")
            assert slvk[5] in ("XXXX",)
            assert slvk[6] in ("XXXX",)
            assert slvk[7] in ("XXXX",)
            assert slvk[8] in ("XXXX",)
            assert slvk[10] in ("XXXX",)
            assert slvk[11] in ("XXXX",)
        elif method == "LDLT":
            assert slvk[1] in ("XXXX",)
            assert slvk[2] in ("XXXX",)
            assert slvk[3] in ("RCMK",)
            assert slvk[5] in ("XXXX",)
            assert slvk[6] in ("XXXX",)
            assert slvk[7] in ("XXXX",)
            assert slvk[8] in ("XXXX",)
            assert slvk[10] in ("XXXX",)
            assert slvk[11] in ("XXXX",)
        elif method == "GCPC":
            assert slvk[1] in ("LDLT_INC", "LDLT_SP", "LDLT_DP", "SANS")
            assert slvk[3] in ("SANS", "RCMK")
            assert slvk[5] in ("XXXX",)
            assert slvk[6] in ("XXXX",)
            assert slvk[7] in ("XXXX",)
            assert slvk[8] in ("IN_CORE", "AUTO", "XXXX")
            assert slvk[10] in ("XXXX",)
            assert slvk[11] in ("XXXX",)
        elif method == "PETSC":
            assert slvk[1] in (
                "LDLT_INC",
                "LDLT_SP",
                "LDLT_DP",
                "JACOBI",
                "SOR",
                "ML",
                "BOOMER",
                "SANS",
            )
            assert slvk[3] in ("SANS", "RCMK")
            assert slvk[5] in ("CG", "CR", "GMRES", "GCR")
            assert slvk[6] in ("XXXX",)
            assert slvk[7] in ("XXXX",)
            assert slvk[8] in ("XXXX",)
            assert slvk[10] in ("XXXX",)
            assert slvk[11] in ("XXXX",)
        else:
            assert False, method
