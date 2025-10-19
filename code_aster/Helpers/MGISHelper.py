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

import os
import os.path as osp
import shutil
import tempfile
from subprocess import call

from ..Behaviours import catalc
from ..Messages import UTMESS
from ..Objects import MGISBehaviour
from ..Utilities import ExecutionParameter, config, logger


class MGISBuilder:
    """Object to easily create MGISBehaviour objects from '.mfront' user file,
    precompiled library or official behaviours.
    """

    _mgisnames = []

    @classmethod
    def from_embedded(cls, behaviour_name):
        """Create MGISBehaviour for an embedded behaviour.

        Arguments:
            behaviour_name (str): Name of the behaviour.
        """
        libpath = osp.join(
            os.environ["ASTER_LIBDIR"], "lib" + config["ASTER_BEHAVIOUR_LIB"] + ".so"
        )
        return MGISBuilder.from_library(behaviour_name, libpath)

    @classmethod
    def from_library(cls, behaviour_name, lib):
        """Create MGISBehaviour for an embedded behaviour.

        Arguments:
            behaviour_name (str): Name of the behaviour.
            lib (str): Path to a MGIS precompiled library.
        """
        builder = MGISBuilder(behaviour_name, lib=lib)
        return builder.build()

    @classmethod
    def from_source(cls, behaviour_name, src, lib=None, flags=None):
        """Create MGISBehaviour for an embedded behaviour.

        Arguments:
            behaviour_name (str): Name of the behaviour.
            src (str): Path to a MFront user source file.
            lib (str, optional): Path of the compiled MGIS library.
            flags (list[str]): Additional compilation flags.
        """
        builder = MGISBuilder(behaviour_name, src=src, lib=lib)
        builder.compile(flags)
        return builder.build()

    @classmethod
    def get_mgis_names(cls):
        """Return the list of MGIS behaviours names, embedded + prototype.

        Returns:
            list[str]: List of names accepted under RELATION.
        """
        if not cls._mgisnames:
            cls._mgisnames = [law.nom for law in catalc if law.ldctype == "mfront"]
            cls._mgisnames.append("MFRONT")
        return cls._mgisnames

    def __init__(self, behaviour_name, src=None, lib=None):
        self._name = behaviour_name
        self._src = src
        self._lib = lib

    def build(self):
        """Build and return the behaviour object.

        Returns:
            MGISBehaviour: Result object.
        """
        mgb = MGISBehaviour()
        mgb.setLibPath(self._lib)
        mgb.setBehaviourName(self._name)
        return mgb

    def compile(self, flags=None):
        """Compile a '.mfront' source file."""
        if not self._lib:
            self._lib = tempfile.NamedTemporaryFile(prefix="libMGIS", suffix=".so", dir=".").name
        cmd = [ExecutionParameter().get_option("prog:mfront"), "--build", "--interface=generic"]
        cmd.extend(flags or [])
        cmd.append(self._src)
        logger.debug("Execute command %r", cmd)
        try:
            call(cmd)
            libname = "libBehaviour.so"
            filename = osp.join("src", libname)
            if not osp.exists(filename):
                UTMESS("F", "MFRONT_4", valk=libname)
            shutil.copyfile(filename, self._lib)
        finally:
            for dname in ("src", "include"):
                if osp.isdir(dname):
                    shutil.rmtree(dname)


def adapt_for_mgis_behaviour(self, keywords):
    """Hook to adapt syntax *after* syntax checking.

    Arguments:
        keywords (dict): Keywords arguments of user's keywords, changed in place.
    """
    mgis_names = set(MGISBuilder.get_mgis_names())
    for occ in keywords.get("COMPORTEMENT", []):
        if occ.get("COMPOR_MFRONT"):
            continue
        rela = occ["RELATION"]
        assert rela != "MFRONT"
        kits = occ.get("RELATION_KIT", [])
        assert "MFRONT" not in kits
        inter = mgis_names.intersection(kits)
        if inter:
            assert len(inter) == 1, inter
            rela = inter.pop()
        if rela in mgis_names:
            occ["COMPOR_MFRONT"] = MGISBuilder.from_embedded(rela)
