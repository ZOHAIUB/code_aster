# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

"""
This module defines objects for the testing feature.
"""

import re
from glob import glob

import aster

from ..Messages import UTMESS

_trans = str.maketrans("e", "E")


def _fortran(srepr):
    """for fortran look"""
    return srepr.translate(_trans)


class TestResult:
    """This class provides the feature to print the testcase results.
    A singleton object is created to avoid to repeat some global tasks.
    """

    _isVerif = None

    @classmethod
    def isVerif(cls):
        """Tell if the testcase is a verification one"""
        if cls._isVerif is None:
            cls._isVerif = cls._checkVerif()
            if not cls._isVerif:
                UTMESS("I", "TEST0_19")
        return cls._isVerif

    @staticmethod
    def _printLine(text):
        aster.affiche("RESULTAT", text)

    @classmethod
    def write(cls, width, *args):
        """shortcut to print in the RESULTAT file"""
        fmtval = "%%-%ds" % width
        fmtcols = ["%-4s ", "%-16s", "%-16s", fmtval, fmtval, "%-16s", "%-16s"]
        assert len(args) <= 7, args
        fmt = " ".join(fmtcols[: len(args)])
        line = fmt % args
        cls._printLine(line)
        return line

    @classmethod
    def showResult(cls, type_ref, legend, label, skip, relative, tole, ref, val, compare=1.0):
        """Print a table for TEST_RESU

        type_ref : ANALYTIQUE, NON_REGRESSION, AUTRE_ASTER...
        legend : component name or XXXX
        label : boolean to print or not the labels
        skip : boolean to skip the test and print an empty line
        relative : boolean, True if for relative, False for absolute comparison
        tole : maximum error tolerated
        ref : reference value (integer, real or complex)
        val : computed value (same type as ref)
        compare : order of magnitude
        """
        # ignore NON_REGRESSION tests for validation testcases
        isNonRegr = type_ref.strip() == "NON_REGRESSION"
        isValidIgn = isNonRegr and not cls.isVerif()
        lines = ["pass in showResult"]
        # compute
        diag = "SKIP"
        error = "-"
        if not skip:
            error = abs(1.0 * ref - val)
            tole = 1.0 * tole
            if relative:
                ok = error <= abs((tole * ref))
                tole = tole * 100.0
                if ref != 0.0:
                    error = error / abs(ref) * 100.0
                elif ok:
                    error = 0.0
                else:
                    error = 999.999999
            else:
                tole = abs(tole * compare)
                ok = error <= tole
            diag = " OK " if ok else "NOOK"
        else:
            # do not warn if validation testcase
            if not isValidIgn:
                UTMESS("I", "TEST0_12")
        # formatting
        sref = "%s" % ref
        sval = "%s" % val
        width = max([16, len(sref), len(sval)]) + 2
        serr = "%s" % error
        if len(serr) > 15:
            serr = "%13.6e" % error
        stol = "%s" % tole
        if relative:
            serr += "%"
            stol += "%"
        sref, sval, serr, stol = [_fortran(i) for i in [sref, sval, serr, stol]]
        if diag == "SKIP":
            legend = sref = sval = serr = stol = "-"
        # printing
        if compare != 1.0:
            lines.append(cls.write(width, " ", "ORDRE DE GRANDEUR :", compare))
        if label:
            lines.append(
                cls.write(
                    width, " ", "REFERENCE", "LEGENDE", "VALE_REFE", "VALE_CALC", "ERREUR", "TOLE"
                )
            )
        if isValidIgn:
            lines.append(cls.write(width, "-", type_ref, legend, sref, sval, serr, "-"))
        else:
            lines.append(cls.write(width, diag, type_ref, legend, sref, sval, serr, stol))
        return lines

    @staticmethod
    def _checkVerif():
        """Check if the current execution is for a verification testcase
        (and not a validation one)."""
        exports = glob("*.export")
        if not exports:
            # export file not found, return "verification" that is more strict!
            return True
        with open(exports[0], "r") as f:
            text = f.read()
        expr = re.compile("^P +testlist.*validation", re.M)
        isVerif = expr.search(text) is None
        return isVerif


def testresu_print(type_ref, legend, label, skip, relative, tole, ref, val, compare=1.0):
    """Print a table for TEST_RESU

    type_ref : ANALYTIQUE, NON_REGRESSION, AUTRE_ASTER...
    legend : component name or XXXX
    label : boolean to print or not the labels
    skip : boolean to skip the test and print an empty line
    relative : boolean, True if for relative, False for absolute comparison
    tole : maximum error tolerated
    ref : reference value (integer, real or complex)
    val : computed value (same type as ref)
    compare : order of magnitude
    """
    lines = TestResult.showResult(type_ref, legend, label, skip, relative, tole, ref, val, compare)
    return lines


if __name__ == "__main__":
    testresu_print("NON_REGRESSION", "DX", True, False, False, 1.0e-6, 1.123e-6, 0.0, compare=275.0)
    testresu_print("AUTRE_ASTER", "DX", False, False, False, 1.0e-6, 1.123e-6, 0.0)
    print()

    testresu_print("NON_REGRESSION", "DX", True, True, False, 1.0e-6, 1.123e-6, 0.0)
    testresu_print(
        "NON_REGRESSION", "XXXXX", True, False, False, 1.0e-6, 1.123e-3, 0.0, compare=275.0
    )
    print()

    testresu_print("NON_REGRESSION", "XXXXX", True, False, True, 1.0e-6, 1.123e-2, 0.0)
    print()

    testresu_print("NON_REGRESSION", "XXXXX", True, False, True, 0.02, 456, 458)
    print()

    testresu_print("ANALYTIQUE", "DEPL_C", True, False, True, 1.0e-4, 1.0 + 1.0j, -0.5 + 0.99j)
    print()
