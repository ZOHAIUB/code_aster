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

# person_in_charge: mathieu.courtois at edf.fr

"""
Module fournissant quelques fonctions utilitaires.
"""

import contextlib
import os
import os.path as osp
import re
import socket
import tempfile
import time
import warnings
from subprocess import Popen

try:
    from asrun import create_run_instance

    HAS_ASRUN = True
except ImportError:
    HAS_ASRUN = False

import aster

from .logger import logger
from .mpi_utils import MPI
from .strfunc import convert, maximize_lines
from .version import get_version

DEBUG = False
TMPDIR = os.environ.get("ASTER_TMPDIR") or os.environ.get("TMPDIR") or "/tmp"


def set_debug(value):
    """Positionne la variable de déboggage"""
    global DEBUG
    DEBUG = value


def _print(*args):
    """Fonction 'print'."""
    l_str = []
    for arg in args:
        if type(arg) not in (str, str):
            arg = repr(arg)
        l_str.append(arg)
    text = convert(" ".join(l_str))
    aster.affiche("MESSAGE", text)


def _printDBG(*args):
    """Print debug informations."""
    if not DEBUG:
        return
    _print(*args)


# les commandes fortran pourraient appeler cette fonction
def get_titre_concept(co=None):
    """Retourne un titre automatique."""
    # ASTER 10.01.25 CONCEPT tab0 CALCULE LE 21/05/2010 A 17:58:50 DE TYPE
    # TABLE_SDASTER
    fmt = {
        "version": "ASTER %(version)s",
        "nomco": "CONCEPT %(nom_concept)s",
        "etatco": "CALCULE",
        "dateheure": "%(dateheure)s",
        "typeco": "DE TYPE %(type_concept)s",
    }
    format = [fmt["version"]]
    dfmt = {"version": get_version(), "dateheure": time.strftime("LE %m/%d/%Y A %H:%M:%S")}
    if co:
        dfmt["nom_concept"] = co.getName()
        format.append(fmt["nomco"])
        format.append(fmt["etatco"])
    format.append(fmt["dateheure"])
    if co:
        dfmt["type_concept"] = co.getType().upper()
        format.append(fmt["typeco"])
    globfmt = " ".join(format)
    titre = globfmt % dfmt
    ltit = titre.split()
    ltit = maximize_lines(ltit, 80, " ")
    return os.linesep.join(ltit)


def fmtF2PY(fformat):
    """Convertit un format Fortran en format Python (printf style).
    Gère uniquement les fortrans réels, par exemple : E12.5, 1PE13.6, D12.5...
    """
    fmt = ""
    matP = re.search("([0-9]+)P", fformat)
    if matP:
        fmt += " " * int(matP.group(1))
    matR = re.search(r"([eEdDfFgG]{1})([\.0-9]+)", fformat)
    if matR:
        fmt += "%" + matR.group(2) + re.sub("[dD]+", "E", matR.group(1))
    try:
        _ = fmt % -0.123
    except (ValueError, TypeError) as msg:
        fmt = "%12.5E"
        print("Error :", msg)
        print("Format par défaut utilisé :", fmt)
    return fmt


def encode_str(string):
    """Convert a string in an array of int"""
    return [ord(i) for i in string]


def decode_str(array):
    """Convert an array of int in a string"""
    return "".join([chr(i) for i in array])


def send_file(fname, dest):
    """Send a file into an existing remote destination directory using scp.

    Arguments:
        fname (str): local file name.
        dest (str): destination path '[host:]path'.
    """
    local = socket.gethostname()
    host, path = dest.split(":", maxsplit=1)
    if host in ("", local):
        dest = path
    dst = osp.join(dest, osp.basename(fname))
    proc = Popen(["scp", "-rBCq", "-o", "StrictHostKeyChecking=no", fname, dst])
    return proc.wait()


def get_time():
    """Return the current time with milliseconds"""
    ct = time.time()
    msec = (ct - int(ct)) * 1000
    return time.strftime("%H:%M:%S") + ".%03d" % msec


def get_shared_tmpdir(prefix, default_dir=TMPDIR):
    """Return a shared temporary directory.

    If asrun shared tmpdir is not known, use *default_dir*.
    """
    shared_tmp = default_dir
    if HAS_ASRUN:
        if getattr(get_shared_tmpdir, "cache_run", None) is None:
            get_shared_tmpdir.cache_run = create_run_instance(
                debug_stderr=False, log_progress="asrun.log"
            )
        run = get_shared_tmpdir.cache_run
        shared_tmp = run.get("shared_tmp")

    shared_tmp = shared_tmp or os.getcwd()
    tmpdir = tempfile.mkdtemp(dir=shared_tmp, prefix=prefix)
    return tmpdir


class SharedTmpdir(contextlib.AbstractContextManager):
    """Create a shared temporary directory with automatic cleanup.

    It can be used as a context manager::

        with SharedTmpdir("foo") as tmpdir:
            # Then use 'tmpdir.path'
            # The directory is removed when exiting the 'with' block.

    Or::

        tmpdir = SharedTmpdir("bar")
        # Then use 'tmpdir.path'
        # The directory will be removed when 'tmpdir' will be deleted,
        # may be forced with:
        del tmpdir

    Arguments:
        prefix (str): Prefix to be used for the temporary directory.
        default_dir (str): Local pathname that will be used if *asrun* is not
            available (default: /tmp).
    """

    def __init__(self, prefix, default_dir=TMPDIR):
        self._prefix = prefix
        self._dir = default_dir
        self._path = None

    @property
    def path(self):
        """str: Attribute containing the path to the temporary directory."""
        if not self._path:
            self._mkdtemp()
        return self._path

    def __enter__(self):
        self._mkdtemp()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self._rmdtemp()

    def __del__(self):
        self._rmdtemp()

    def _mkdtemp(self):
        rank = MPI.ASTER_COMM_WORLD.rank
        tmpdir = ""
        if rank == 0:
            # choose a non existing directory on proc #0
            # as the directory is shared, it should not exist on others
            tmpdir = get_shared_tmpdir(self._prefix, self._dir)
        # wait for #0 to create the temporary directory
        tmpdir = MPI.ASTER_COMM_WORLD.bcast(tmpdir)
        assert osp.isdir(tmpdir), f"Can not create directory on #{rank}: {tmpdir}"
        self._path = tmpdir
        logger.debug("created %s", self._path)
        MPI.ASTER_COMM_WORLD.Barrier()

    def _rmdtemp(self):
        MPI.ASTER_COMM_WORLD.Barrier()
        if MPI.ASTER_COMM_WORLD.rank == 0:
            logger.debug("removing %s...", self._path)
            os.system(f"rm -rf {self._path}")
        MPI.ASTER_COMM_WORLD.Barrier()


@contextlib.contextmanager
def shared_tmpdir(prefix, default_dir=TMPDIR):
    """Return a shared temporary directory with automatic cleanup, to be used
    as a context manager.

    Arguments:
        prefix (str): Prefix to be used for the temporary directory.
        default_dir (str): Local pathname that will be used if *asrun* is not
            available (default: /tmp).

    Returns:
        str: Path of the temporary directory.
    """
    warnings.warn("Prefer use SharedTmpdir object", DeprecationWarning)
    tmpdir = SharedTmpdir(prefix, default_dir)
    yield tmpdir.path
