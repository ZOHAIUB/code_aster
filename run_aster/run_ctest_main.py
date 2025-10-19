#!/usr/bin/env python3
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

"""
``bin/run_ctest`` --- Script to execute code_aster testcases using ``ctest``
----------------------------------------------------------------------------

``bin/run_ctest`` executes code_aster testcases using ``ctest``.

Usage:

.. code-block:: sh

    bin/run_ctest [options] [ctest-options] [other arguments...]

`ctest-options` and `other arguments` are passed to ``ctest``.

The list of the testcases to be executed is built from ``--testlist`` argument
and a filter on the labels (taken from the ``.export`` files).

The *sequential* label is automatically added for a sequential version.

To show the list of labels, use:

.. code-block:: sh

    bin/run_ctest --print-labels

.. note::

  Difference from ``ctest``: all values passed to ``-L`` option are sorted and
  joined as a unique regular expression.

  Example:

    Using:

    .. code-block:: sh

        bin/run_ctest -L verification -L ci -N

    the ctest command will be:

    .. code-block:: sh

        ctest -N -j 6 -L 'ci.*verification'

See ``bin/run_ctest --help`` for the available options.

"""

import argparse
import os
import os.path as osp
import re
import sys
import tempfile
from glob import glob
from subprocess import run
from pathlib import Path

from .config import CFG
from .ctest2junit import XUnitReport
from .run import get_nbcores
from .utils import RUNASTER_ROOT, RUNASTER_PLATFORM

USAGE = """
    run_ctest [options] [ctest-options] [other arguments...]

Execute testcases using 'ctest'.
'ctest-options' and other arguments are passed to 'ctest'.
Use 'ctest --help' for details.

The list of the testcases to be executed is built from '--testlist' argument
and a filter on the labels (taken from the '.export' files).

The label 'sequential' is automatically added for a sequential version.

To show the list of labels, use:

    run_ctest --print-labels

To execute all sslp01* and zzzz100a/f testcases, use:

    run_ctest -R 'sslp01|zzzz100[af]'

The '--nlist=N' option creates N directories that can be run separately
(in different jobs under a batch scheduler for example). The lists may not be
well balanced if '-R' or '-L' filters are used.

Note:
  Difference from 'ctest': all values passed to '-L' option are sorted and joined
  as a unique regular expression.

  Example:

    Using:
        run_ctest -L verification -L ci -N

    the ctest command will be:

        ctest -N -j 8 -L 'ci.*verification'
"""


def parse_args(argv):
    """Parse command line arguments.

    Arguments:
        argv (list): List of command line arguments.
    """
    # command arguments parser
    parser = argparse.ArgumentParser(usage=USAGE)
    jobs = max(1, get_nbcores() - 2)
    parser.add_argument(
        "-j",
        "--jobs",
        action="store",
        type=int,
        default=jobs,
        help="run the tests in parallel using the given " f"number of jobs (default: {jobs})",
    )
    parser.add_argument(
        "--testlist", action="store", metavar="FILE", help="list of testcases to run"
    )
    parser.add_argument(
        "--exclude-testlist",
        action="store",
        metavar="FILE",
        help="list of testcases to be excluded",
    )
    parser.add_argument(
        "--resutest",
        action="store",
        metavar="DIR",
        help="directory to write the results of the testcases "
        "(relative to the current directory)",
    )
    parser.add_argument(
        "--no-resutest",
        action="store_const",
        dest="resutest",
        const="none",
        help="do not keep result files",
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        default="auto",
        help="remove the content of 'resutest' directory before starting",
    )
    parser.add_argument(
        "--no-clean",
        action="store_false",
        default="auto",
        dest="clean",
        help="do not remove the content of 'resutest' directory at startup",
    )
    parser.add_argument(
        "--timefactor",
        action="store",
        type=float,
        default=1.0,
        help="multiplicative factor applied to the time limit, "
        "passed through environment to run_aster "
        "(default: 1.0)",
    )
    parser.add_argument(
        "--only-failed-results",
        action="store_true",
        default=False,
        help="keep only the results of tests in failure",
    )
    run_sbatch = osp.join(RUNASTER_ROOT, "bin", "run_sbatch")
    parser.add_argument(
        "--sbatch",
        action="store_true",
        default=False,
        help=f"call run_sbatch instead of run_aster (see '{run_sbatch} --help' for specific options)",
    )
    parser.add_argument(
        "--nlist",
        metavar="N",
        action="store",
        type=int,
        default=None,
        help="create several CTestTestfile.cmake files, do not run testcases, show command lines",
    )
    parser.add_argument(
        "--report",
        action="store_true",
        default=None,
        help="create a 'run_testcases.xml' report file",
    )
    parser.add_argument(
        "--testdir",
        metavar="DIR",
        action="store",
        help="alternative directory containing the test cases (default: <installdir>/share/aster/tests)",
    )
    group = parser.add_argument_group("ctest options")
    group.add_argument(
        "--rerun-failed", action="store_true", help="run only the tests that failed previously"
    )
    group.add_argument(
        "-L",
        "--label-regex",
        action="append",
        metavar="regex",
        default=[],
        help="run tests with labels matching regular expression.",
    )
    group.add_argument(
        "-LE",
        "--label-exclude",
        action="append",
        metavar="regex",
        default=[],
        help="exclude tests with labels matching regular expression.",
    )
    group.add_argument(
        "--print-labels", action="store_true", help="print all available test labels"
    )

    args, others = parser.parse_known_args(argv)
    # resutest needed or not?
    if args.print_labels or "-N" in others or "--show-only" in others:
        args.resutest = "none"
    if not args.resutest:
        parser.error("'--resutest' argument is required")

    # args to be re-injected for ctest
    if args.print_labels:
        others.append("--print-labels")
    if args.rerun_failed:
        others.append("--rerun-failed")
    others.extend(["-j", str(args.jobs)])

    return args, others


def _run(cmd, shell=False):
    print("execute:", " ".join(cmd))
    return run(cmd, shell=shell)


def _rmtree(dirpath):
    if RUNASTER_PLATFORM == "linux":
        _run(["rm", "-rf", dirpath])
    else:
        _run(["rd", "/s", "/q", dirpath], shell=True)


def main(argv=None):
    """Entry point for testcases runner.

    Arguments:
        argv (list): List of command line arguments.
    """
    args, ctest_args = parse_args(argv or sys.argv[1:])
    # options passed through environment
    os.environ["FACMTPS"] = str(args.timefactor)
    if args.only_failed_results:
        os.environ["ASTER_ONLY_FAILED_RESULTS"] = "1"
        print("only the results files of testcases in failure will be kept!")

    use_tmp = args.resutest.lower() == "none"
    if use_tmp:
        resutest = tempfile.mkdtemp(prefix="resutest_")
    else:
        resutest = osp.join(os.getcwd(), args.resutest)

    if args.report and args.nlist:
        report = XUnitReport(resutest, "")
        for icount in range(args.nlist):
            resdir = osp.join(resutest, f"{icount + 1:03d}")
            report.read_ctest(resdir)
        report.write_xml("run_testcases.xml")
        return 0

    if not use_tmp and args.clean and not args.rerun_failed:  # clean = True or 'auto'
        if args.clean == "auto" and osp.exists(resutest):
            print(f"{resutest} will be removed.")
            answ = input("do you want to continue (y/n) ?")
            if answ.lower() not in ("y", "o"):
                print("interrupt by user")
                sys.exit(1)
        _rmtree(resutest)

    if not osp.exists(resutest):
        os.makedirs(resutest, exist_ok=True)

    testlist = osp.abspath(args.testlist) if args.testlist else ""
    excl = osp.abspath(args.exclude_testlist) if args.exclude_testlist else ""
    opts = "--sbatch" if args.sbatch else ""
    if not args.rerun_failed:
        # create CTestTestfile.cmake
        create_ctest_file(testlist, excl, resutest, opts, args.nlist, testdir=args.testdir)
    parallel = CFG.get("parallel", 0)
    labels = set()
    if not parallel:
        labels.add("sequential")
    if labels or args.label_regex:
        labels.update(args.label_regex)
        ctest_args.extend(["-L", ".*".join(sorted(labels))])
    if args.label_exclude:
        ctest_args.extend(["-LE", "|".join(sorted(args.label_exclude))])

    os.chdir(resutest)
    # show command lines to be run
    if args.nlist:
        print("Command lines for each batch job:")
        for icount in range(args.nlist):
            wrkdir = osp.join(resutest, f"{icount + 1:03d}")
            print(f"  cd {wrkdir} && ctest " + " ".join([f"'{i}'" for i in ctest_args]))
        print("Build consolidated report:")
        print(f"  cd {resutest}")
        bindir = osp.normpath(osp.join(RUNASTER_ROOT, "bin"))
        print("  " + osp.join(bindir, "run_ctest") + " " + " ".join(sys.argv[1:]) + " --report")
        return 0

    # execute ctest
    proc = _run(["ctest"] + ctest_args)
    if not use_tmp:
        legend = ""
        if testlist:
            legend += " from " + osp.join(*testlist.split(osp.sep)[-2:])
        report = XUnitReport(resutest, legend)
        report.read_ctest()
        report.write_xml("run_testcases.xml")
    else:
        _rmtree(resutest)
    return proc.returncode


def create_ctest_file(testlist, exclude, destdir, options, nlist=None, testdir=None):
    """Create the CTestTestfile.cmake file.

    Arguments:
        testlist (str): file containing a list of testcases.
        exclude (str): file containing a list of testcases to be excluded.
        destdir (str): Destination directory for the 'ctest' file(s).
        options (str): Additional command line options.
        nlist (int, optional): Number of files to be created.
        testdir (str, optional): directory containing the testcases.
    """
    datadir = Path(osp.normpath(osp.join(RUNASTER_ROOT, "share", "aster"))).as_posix()
    if testdir is None:
        testdir = osp.join(datadir, "tests")
    testdir = Path(testdir).absolute().as_posix()

    assert osp.isdir(testdir), f"no such directory {testdir}"
    re_comment = re.compile("^ *#.*$", re.M)
    if osp.isfile(testlist):
        with open(testlist, "r") as fobj:
            text = re_comment.sub("", fobj.read())
        ltests = set(text.split())
        lexport = [osp.join(testdir, tst + ".export") for tst in ltests]
    else:
        lexport = glob(osp.join(testdir, "*.export"))
    if osp.isfile(exclude):
        with open(exclude, "r") as fobj:
            text = re_comment.sub("", fobj.read())
        ltests = set(text.split())
        excl = [osp.join(testdir, tst + ".export") for tst in ltests]
        lexport = list(set(lexport).difference(excl))

    tag = CFG.get("version_tag", "")
    size = len(lexport) // (nlist or 1)
    if len(lexport) % (nlist or 1):
        size += 1
    icount = 0
    while lexport:
        icount += 1
        text = [
            f"set(COMPONENT_NAME ASTER_{tag})",
            _build_def(datadir, lexport[:size], options, testdir),
        ]
        lexport = lexport[size:]
        filename = osp.join(destdir, "CTestTestfile.cmake")
        if nlist:
            filename = osp.join(destdir, f"{icount:03d}", "CTestTestfile.cmake")
        os.makedirs(osp.dirname(filename), exist_ok=True)
        with open(filename, "w") as fobj:
            fobj.write("\n".join(text))


CTEST_DEF = """
set(TEST_NAME ${{COMPONENT_NAME}}_{testname})
add_test(${{TEST_NAME}} {ASTERDATADIR}/run_aster_for_ctest{ext} {options} {TESTDIR}/{testname}.export)
set_tests_properties(${{TEST_NAME}} PROPERTIES
                     LABELS "${{COMPONENT_NAME}} {labels}"
                     PROCESSORS {processors}
                     TIMEOUT {timeout}
                     COST {timeout})
"""

TEST_FILES_INTEGR = """
forma02a
forma01c
mumps01a
mfron01a
zzzz151a
zzzz200b
zzzz218a
zzzz401a
"""


def _build_def(datadir, lexport, options, testdir):
    re_list = re.compile("P +testlist +(.*)$", re.M)
    re_nod = re.compile("P +mpi_nbnoeud +([0-9]+)", re.M)
    re_mpi = re.compile("P +mpi_nbcpu +([0-9]+)", re.M)
    re_thr = re.compile("P +ncpus +([0-9]+)", re.M)
    re_time = re.compile("P +time_limit +([0-9]+)", re.M)
    text = []
    for exp in lexport:
        if not osp.isfile(exp):
            print(f"no such file: {exp}")
            continue
        testname = osp.splitext(osp.basename(exp))[0]
        lab = []
        nod = 1
        mpi = 1
        thr = 1
        tim = 86400
        with open(exp, "r") as fobj:
            export = fobj.read()
        mat = re_list.search(export)
        if mat:
            lab = mat.group(1).split()
        mat = re_nod.search(export)
        if mat:
            nod = int(mat.group(1))
        mat = re_mpi.search(export)
        if mat:
            mpi = int(mat.group(1))
        mat = re_thr.search(export)
        if mat:
            thr = int(mat.group(1))
        mat = re_time.search(export)
        if mat:
            tim = int(mat.group(1))
        lab.append(f"nodes={nod:02d}")
        if testname in TEST_FILES_INTEGR:
            lab.append("SMECA_INTEGR")
        procs = mpi * thr
        timeout = int(tim * 1.1 * float(os.environ["FACMTPS"]))
        if "sbatch" in options:
            timeout = 12 * 3600
            procs = 1
        text.append(
            CTEST_DEF.format(
                testname=testname,
                labels=" ".join(sorted(lab)),
                processors=procs,
                timeout=timeout,
                options=options,
                ASTERDATADIR=datadir,
                TESTDIR=testdir,
                ext=".bat" if RUNASTER_PLATFORM == "win" else "",
            )
        )
    return "\n".join(text)


if __name__ == "__main__":
    sys.exit(main())
