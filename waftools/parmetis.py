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

import re
from functools import partial
from waflib import Configure, Utils, Errors


def options(self):
    group = self.add_option_group("Parmetis library options")
    group.add_option(
        "--disable-parmetis",
        action="store_false",
        default=None,
        dest="enable_parmetis",
        help="Disable PARMETIS support",
    )
    group.add_option(
        "--enable-parmetis",
        action="store_true",
        default=None,
        dest="enable_parmetis",
        help="Force PARMETIS support",
    )
    group.add_option(
        "--parmetis-libs",
        type=str,
        dest="parmetis_libs",
        default=None,
        help="parmetis librairies to use when linking",
    )
    group.add_option(
        "--embed-parmetis",
        dest="embed_parmetis",
        default=False,
        action="store_true",
        help="Embed PARMETIS libraries as static library",
    )


def configure(self):
    if not self.env.BUILD_MPI:
        self.undefine("ASTER_HAVE_PARMETIS")
        return
    try:
        self.env.stash()
        self.check_parmetis()
    except Errors.ConfigurationError:
        self.reset_msg()
        self.env.revert()
        self.undefine("ASTER_HAVE_PARMETIS")
        if self.options.enable_parmetis:
            raise
    else:
        self.define("ASTER_HAVE_PARMETIS", 1)


###############################################################################
@Configure.conf
def check_parmetis(self):
    opts = self.options
    if opts.enable_parmetis is False:
        raise Errors.ConfigurationError("PARMETIS disabled")
    if opts.parmetis_libs is None:
        opts.parmetis_libs = "parmetis"
    if opts.parmetis_libs:
        self.check_parmetis_libs()
    self.check_parmetis_idx()
    self.check_parmetis_real()
    self.check_parmetis_version()


@Configure.conf
def check_parmetis_libs(self):
    opts = self.options
    check_ = partial(self.check_cc, uselib_store="PARMETIS", use="PARMETIS METIS M", mandatory=True)
    if opts.embed_all or opts.embed_parmetis:
        check = lambda lib: check_(stlib=lib)
    else:
        check = lambda lib: check_(lib=lib)
    list(map(check, Utils.to_list(opts.parmetis_libs)))


@Configure.conf
def check_parmetis_idx(self):
    fragment = r"""
#include <stdio.h>
#include <parmetis.h>
int main(void){
    idx_t idx;
    printf("%d", (int)sizeof(idx));
    return 0;
}"""
    self.code_checker(
        "ASTER_PARMETIS_IDX_SIZE",
        self.check_cc,
        fragment,
        "Checking size of parmetis idx_t",
        "unexpected value for sizeof(idx_t): %(size)s",
        into=(4, 8),
        use="PARMETIS METIS M",
    )


@Configure.conf
def check_parmetis_real(self):
    fragment = r"""
#include <stdio.h>
#include <parmetis.h>
int main(void){
    real_t real;
    printf("%d", (int)(sizeof(real)));
    return 0;
}"""
    self.code_checker(
        "ASTER_PARMETIS_REAL_SIZE",
        self.check_cc,
        fragment,
        "Checking size of parmetis real_t",
        "unexpected value for sizeof(real_t): %(size)s",
        into=(4, 8),
        use="PARMETIS METIS M",
    )


@Configure.conf
def check_parmetis_version(self):
    fragment = r"""
#include <stdio.h>
#include <parmetis.h>
int main(void) {
    printf("PARMETISVER: %d.%d.%d", PARMETIS_MAJOR_VERSION, PARMETIS_MINOR_VERSION,
           PARMETIS_SUBMINOR_VERSION);
    return 0;
}"""
    self.start_msg("Checking parmetis version")
    try:
        ret = self.check_cc(
            fragment=fragment, use="PARMETIS METIS M", mandatory=True, execute=True, define_ret=True
        )
        mat5 = re.search(r"PARMETISVER: *(?P<vers>[0-9]+\.[0-9]+\.\w+)", ret)
        vers = mat5 and mat5.group("vers")
        major = int(vers.split(".")[0])
        if major != 4:
            self.end_msg("unsupported parmetis version: %s" % vers, "RED")
            raise Errors.ConfigurationError
    except:
        self.end_msg("can not get version", "RED")
        raise
    else:
        self.end_msg(vers)
