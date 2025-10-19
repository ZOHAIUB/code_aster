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

import os.path as osp
from functools import partial
from waflib import Configure, Utils, Errors


def options(self):
    group = self.add_option_group("HDF5/Med libraries options")
    group.add_option(
        "--med-libs",
        type=str,
        dest="med_libs",
        default=None,
        help="MED librairies to link against med",
    )
    group.add_option(
        "--embed-med",
        dest="embed_med",
        default=False,
        action="store_true",
        help="Embed MED libraries as static library",
    )
    group.add_option(
        "--disable-med",
        dest="enable_med",
        default=None,
        action="store_false",
        help="disable the MED support",
    )
    group.add_option(
        "--enable-med",
        dest="enable_med",
        default=None,
        action="store_true",
        help="force the MED support",
    )

    group.add_option(
        "--hdf5-libs", type=str, dest="hdf5_libs", default=None, help="HDF5 librairies to link with"
    )
    group.add_option(
        "--embed-hdf5",
        dest="embed_hdf5",
        default=None,
        action="store_true",
        help="Embed HDF5 libraries as static library",
    )
    group.add_option(
        "--disable-hdf5",
        dest="enable_hdf5",
        default=None,
        action="store_false",
        help="disable the HDF5 and MED support",
    )
    group.add_option(
        "--enable-hdf5",
        dest="enable_hdf5",
        default=None,
        action="store_true",
        help="force the HDF5 and MED support",
    )


def configure(self):
    opts = self.options
    try:
        self.env.stash()
        self.check_hdf5()
    except Errors.ConfigurationError:
        self.reset_msg()
        self.env.revert()
        if opts.enable_hdf5 is True:
            raise
        self.undefine("ASTER_HAVE_HDF5")
        self.undefine("ASTER_HAVE_MED")
    else:
        self.define("ASTER_HAVE_HDF5", 1)

    try:
        self.env.stash()
        self.check_med()
    except Errors.ConfigurationError:
        self.reset_msg()
        self.env.revert()
        if opts.enable_med is True:
            raise
        self.undefine("ASTER_HAVE_MED")
    else:
        self.define("ASTER_HAVE_MED", 1)
        self.env.BUILD_MED = True

    try:
        self.env.stash()
        self.check_medcoupling()
    except Errors.ConfigurationError:
        self.reset_msg()
        self.env.revert()


###############################################################################


@Configure.conf
def check_hdf5(self):
    opts = self.options
    if opts.enable_hdf5 is False:
        raise Errors.ConfigurationError("HDF5 disabled")

    if opts.hdf5_libs is None:
        opts.hdf5_libs = "hdf5"

    if opts.hdf5_libs:
        self.check_hdf5_libs()
    self.check_hdf5_headers()
    self.check_hdf5_version()
    self.check_hdf5_api()
    self.check_sizeof_hid()


@Configure.conf
def check_hdf5_libs(self):
    opts = self.options
    check_hdf5 = partial(self.check_cc, mandatory=True, uselib_store="HDF5", use="HDF5 Z")
    if opts.embed_all or opts.embed_hdf5:
        check_lib = lambda lib: check_hdf5(stlib=lib)
    else:
        check_lib = lambda lib: check_hdf5(lib=lib)
    list(map(check_lib, Utils.to_list(opts.hdf5_libs)))


@Configure.conf
def check_hdf5_headers(self):
    check = partial(self.check_cc, header_name="hdf5.h", uselib_store="HDF5", use="HDF5 Z")
    self.start_msg("Checking for header hdf5.h")
    try:
        if not check(mandatory=False):
            if not check(includes=[osp.join(self.env.INCLUDEDIR, "hdf5")], mandatory=False):
                check(includes=[osp.join(self.env.OLDINCLUDEDIR, "hdf5")], mandatory=True)
    except:
        self.end_msg("no", "YELLOW")
        raise
    else:
        self.end_msg("yes")


@Configure.conf
def check_hdf5_version(self):
    fragment = r"""
#include <stdio.h>
#include <hdf5.h>
int main(void){
    int ier;
    unsigned int n1=0, n2=0, n3=0;
    ier = (int)H5get_libversion( &n1, &n2, &n3 );
    printf("%d.%d.%d", n1, n2, n3);
    return 0;
}"""
    self.start_msg("Checking hdf5 version")
    try:
        ret = self.check_cc(
            fragment=fragment, use="HDF5 Z", mandatory=True, execute=True, define_ret=True
        )
    except:
        self.end_msg("no", "YELLOW")
        raise
    else:
        self.end_msg(ret)


@Configure.conf
def check_hdf5_api(self):
    fragv18 = "#include <hdf5.h>\nint main(){hid_t st=0;H5Eclear(st);return 0;}"

    self.start_msg("Checking for API hdf5 v18")
    try:
        self.to_log("check the v18 api and set H5_NO_DEPRECATED_SYMBOLS if it fails")
        check = partial(self.check_cc, execute=True, uselib_store="HDF5", use="HDF5 Z")
        v18 = check(fragment=fragv18, mandatory=False)
        if not v18:
            self.define("H5_NO_DEPRECATED_SYMBOLS", 1)
        self.to_log("try again by using H5_NO_DEPRECATED_SYMBOLS")
        check(fragment=fragv18, mandatory=True)
    except:
        self.end_msg("no", "RED")
        raise
    else:
        self.end_msg(v18 and "default v18" or "-DH5_NO_DEPRECATED_SYMBOLS")


@Configure.conf
def check_sizeof_hid(self):
    fragment = r"""
#include <stdio.h>
#include <hdf5.h>
int main(void){
    hid_t integer;
    printf("%d", (int)sizeof(integer));
    return 0;
}"""
    self.code_checker(
        "ASTER_HDF5_HID_SIZE",
        self.check_cc,
        fragment,
        "Checking size of hid_t integers",
        "unexpected value for sizeof(hid_t): %(size)s",
        into=(4, 8),
        use="HDF5 Z",
    )


@Configure.conf
def check_med(self):
    opts = self.options
    if opts.enable_med is False:
        raise Errors.ConfigurationError("MED disabled")

    self.check_med_libs()
    self.check_med_headers()
    self.check_sizeof_med_int()
    self.check_sizeof_med_idt()
    self.check_med_python()


@Configure.conf
def check_med_libs(self):
    opts = self.options
    candidates = ["med", "medfwrap medC"]
    if opts.med_libs is not None:
        candidates = [opts.med_libs]

    def do_check(libs):
        check_med = partial(self.check_cc, mandatory=True, uselib_store="MED", use="MED HDF5 Z")
        kwd = "stlib" if opts.embed_all or opts.embed_med else "lib"
        for lib in Utils.to_list(libs):
            check_med(**{kwd: lib})

    success = False
    while not success and candidates:
        libsset = candidates.pop(0)
        try:
            self.env.stash()
            do_check(libsset)
            # using MEDlibraryNumVersion symbol is sufficient
            self.check_med_version()
            success = True
        except Errors.ConfigurationError:
            self.env.revert()


@Configure.conf
def check_med_headers(self):
    check = partial(self.check_cc, header_name="med.h", uselib_store="MED", use="MED HDF5 Z")
    self.start_msg("Checking for header med.h")
    try:
        if not check(mandatory=False):
            if not check(includes=[osp.join(self.env.INCLUDEDIR, "med")], mandatory=False):
                check(includes=[osp.join(self.env.OLDINCLUDEDIR, "med")], mandatory=True)
    except:
        self.end_msg("no", "YELLOW")
        raise
    else:
        self.end_msg("yes")


@Configure.conf
def check_med_version(self):
    fragment = r"""
#include <stdio.h>
#include <hdf5.h>
#include <med.h>
int main(void){
    med_bool test; /* med >3.0 */
    med_int n1=0, n2=0, n3=0;
    MEDlibraryNumVersion( &n1, &n2, &n3 );
    printf("%d.%d.%d", n1, n2, n3);
    return 0;
}"""
    self.start_msg("Checking med version")
    try:
        ret = self.check_cc(
            fragment=fragment,
            use="MED HDF5 Z",
            uselib_store="MED",
            mandatory=True,
            execute=True,
            define_ret=True,
        )
    except:
        self.end_msg("no", "YELLOW")
        raise
    else:
        major, minor, release = ret.strip().split(".")
        self.define("ASTER_MED_VERSION_MAJOR", int(major))
        self.define("ASTER_MED_VERSION_MINOR", int(minor))
        self.define("ASTER_MED_VERSION_RELEASE", int(release))
        self.end_msg(ret)


@Configure.conf
def check_med_python(self):
    if not self.env["PYTHON"]:
        self.fatal("load python tool first")
    try:
        self.env.stash()
        self.check_python_module("med")
    except Errors.ConfigurationError:
        self.env.revert()
        if self.env.BUILD_MPI and self.options.with_py_med:
            raise


@Configure.conf
def check_medcoupling(self):
    if not self.env["PYTHON"]:
        self.fatal("load python tool first")
    self.check_python_module("medcoupling")


@Configure.conf
def check_sizeof_med_int(self):
    fragment = r"""
#include <stdio.h>
#include <hdf5.h>
#include <med.h>
int main(void){
    med_int integer;
    printf("%d", (int)sizeof(integer));
    return 0;
}"""
    self.code_checker(
        "ASTER_MED_INT_SIZE",
        self.check_cc,
        fragment,
        "Checking size of med_int integers",
        "unexpected value for sizeof(med_int): %(size)s",
        into=(4, 8),
        use="MED HDF5 Z",
    )


@Configure.conf
def check_sizeof_med_idt(self):
    fragment = r"""
#include <stdio.h>
#include <hdf5.h>
#include <med.h>
int main(void){
    med_idt integer;
    printf("%d", (int)sizeof(integer));
    return 0;
}"""
    self.code_checker(
        "ASTER_MED_IDT_SIZE",
        self.check_cc,
        fragment,
        "Checking size of med_idt integers",
        "unexpected value for sizeof(med_idt): %(size)s",
        into=(4, 8),
        use="MED HDF5 Z",
    )
