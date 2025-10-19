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

from code_aster import CA

CA.init("--test", "--continue")

test = CA.TestCase()

test.assertIsNotNone(ch0.getDescription(), msg="chs has a description")

# same as done for Y to check the pickling of:
# - chs (SimpleFieldOnCellsReal),
# - vol (ndarray),
# - weight (ComponentOnCells)

# \int_{0}^{1} z.dz = 0.5
vz = chs.Z
vz.restrict(vol)
wvol = weight.onSupportOf(vz)
integr = (vz * wvol).sum()
test.assertAlmostEqual(integr, 0.5, msg="integr")

test.printSummary()

CA.close()
