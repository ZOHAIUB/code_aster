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

import numpy as NP

from ..Messages import UTMESS

from ..Cata.Syntax import _F
from ..CodeCommands import (
    AFFE_MATERIAU,
    ASSE_MATRICE,
    CALC_MATR_ELEM,
    COMB_MATR_ASSE,
    DEFI_MATERIAU,
    NUME_DDL,
)
from .Utils.dyna_visco_harm import dyna_visco_harm
from .Utils.dyna_visco_modes import dyna_visco_modes


def dyna_visco_ops(
    self,
    MODELE,
    EXCIT,
    MATER_ELAS_FO,
    MATER_ELAS=None,
    CARA_ELEM=None,
    FREQ=None,
    LIST_FREQ=None,
    TYPE_MODE=None,
    RESI_RELA=None,
    INFO=None,
    TYPE_RESU=None,
    **args
):
    """
    Macro-command DYNA_VISCO, main file
    """

    coef_fmax = 1.0
    if "COEF_FREQ_MAX" in args:
        coef_fmax = args["COEF_FREQ_MAX"]

    list_FREQ = []

    if FREQ:
        list_FREQ = FREQ
    else:
        list_FREQ = LIST_FREQ.getValues()

    # check on the number of values of the frequencies list
    n_f = len(list_FREQ)
    if (TYPE_RESU == "MODE") & (n_f != 2):
        UTMESS("F", "DYNAVISCO_1")
    if (TYPE_RESU == "HARM") & (n_f < 2):
        UTMESS("F", "DYNAVISCO_2")

    fmax = coef_fmax * list_FREQ[-1]

    motscles = {}
    motscles["AFFE"] = []
    __listKv = {}
    ltrv = {}
    e0 = {}
    eta0 = {}
    ny = 0
    trKg = 0

    # ASSEMBLY OF THE MATRICES OF THE ELASTIC AND VISCOELASTIC PARTS

    if MATER_ELAS:
        for n in MATER_ELAS:
            if n["MATER"] is None:
                __new_matnv = DEFI_MATERIAU(
                    ELAS=_F(E=n["E"], NU=n["NU"], RHO=n["RHO"], AMOR_HYST=n["AMOR_HYST"])
                )

                motscles["AFFE"].append(_F(GROUP_MA=n["GROUP_MA"], MATER=__new_matnv))

            else:
                motscles["AFFE"].append(_F(GROUP_MA=n["GROUP_MA"], MATER=n["MATER"]))

    for y in MATER_ELAS_FO:
        e0[ny] = float(y["E"](list_FREQ[0]))
        eta0[ny] = float(y["AMOR_HYST"](list_FREQ[0]))

        __new_matv = DEFI_MATERIAU(ELAS=_F(E=e0[ny], NU=y["NU"], RHO=y["RHO"], AMOR_HYST=eta0[ny]))

        motscles["AFFE"].append(_F(GROUP_MA=y["GROUP_MA"], MATER=__new_matv))

        ny = ny + 1

    __affemat = AFFE_MATERIAU(MODELE=MODELE, **motscles)

    for y in MATER_ELAS_FO:
        del motscles["AFFE"][-1]

    charge = []
    for l_char in EXCIT:
        for l_char2 in l_char["CHARGE"]:
            charge.append(l_char2)

    __Mg = CALC_MATR_ELEM(
        OPTION="MASS_MECA", MODELE=MODELE, CHAM_MATER=__affemat, CARA_ELEM=CARA_ELEM, CHARGE=charge
    )

    __Kgr = CALC_MATR_ELEM(
        OPTION="RIGI_MECA", MODELE=MODELE, CHAM_MATER=__affemat, CARA_ELEM=CARA_ELEM, CHARGE=charge
    )

    __Kg = CALC_MATR_ELEM(
        OPTION="RIGI_MECA_HYST",
        MODELE=MODELE,
        CHAM_MATER=__affemat,
        CARA_ELEM=CARA_ELEM,
        CHARGE=charge,
        RIGI_MECA=__Kgr,
    )

    __num = NUME_DDL(MATR_RIGI=__Kg)

    __asseMg = ASSE_MATRICE(MATR_ELEM=__Mg, NUME_DDL=__num)

    __asseKgr = ASSE_MATRICE(MATR_ELEM=__Kgr, NUME_DDL=__num)

    __asseKg = ASSE_MATRICE(MATR_ELEM=__Kg, NUME_DDL=__num)

    if TYPE_MODE == "BETA_REEL":
        rig = extr_matr_elim_lagr(self, __asseKg)
        trKg = NP.trace(rig)

    # ONE SEPARATES THE CONTRIBUTION OF EACH VISCOELASTIC PART IN THE STIFFNESS MATRIX

    ny = 0

    for y in MATER_ELAS_FO:
        nv = 0
        for v in MATER_ELAS_FO:
            if nv == ny:
                e = 2 * e0[nv]
            else:
                e = e0[nv]

            __new_matv = DEFI_MATERIAU(ELAS=_F(E=e, NU=v["NU"], RHO=v["RHO"]))

            motscles["AFFE"].append(_F(GROUP_MA=v["GROUP_MA"], MATER=__new_matv))

            nv = nv + 1

        __affemat = AFFE_MATERIAU(MODELE=MODELE, **motscles)

        for v in MATER_ELAS_FO:
            del motscles["AFFE"][-1]

        __Kvr = CALC_MATR_ELEM(
            OPTION="RIGI_MECA",
            MODELE=MODELE,
            CHAM_MATER=__affemat,
            CARA_ELEM=CARA_ELEM,
            CHARGE=charge,
        )
        __asseKvr = ASSE_MATRICE(MATR_ELEM=__Kvr, NUME_DDL=__num)

        # copy of __asseKgr, but without the Lagrange DoF
        __asseKgr_sansLagr = COMB_MATR_ASSE(
            COMB_R=_F(MATR_ASSE=__asseKgr, COEF_R=1.0), SANS_CMP="LAGR"
        )

        __listKv[ny] = COMB_MATR_ASSE(
            COMB_R=(
                _F(MATR_ASSE=__asseKvr, COEF_R=1.0),
                _F(MATR_ASSE=__asseKgr_sansLagr, COEF_R=-1.0),
            )
        )

        ################################################################

        if TYPE_MODE == "BETA_REEL":
            rig = extr_matr_elim_lagr(self, __listKv[ny])
            ltrv[ny] = NP.trace(rig)

        ny = ny + 1

    # EIGENMODES COMPUTATION
    _modes = dyna_visco_modes(
        self,
        TYPE_RESU,
        TYPE_MODE,
        list_FREQ,
        fmax,
        RESI_RELA,
        MATER_ELAS_FO,
        __asseKg,
        __asseKgr,
        __asseMg,
        trKg,
        __listKv,
        e0,
        eta0,
        ltrv,
        **args
    )

    # FREQUENCY RESPONSE COMPUTATION
    if TYPE_RESU == "HARM":
        if args.get("MODE_MECA") is not None:
            self.register_result(_modes, args["MODE_MECA"])

        dyna_harm = dyna_visco_harm(
            self,
            EXCIT,
            list_FREQ,
            _modes,
            MATER_ELAS_FO,
            __asseKg,
            __asseKgr,
            __asseMg,
            __listKv,
            e0,
            eta0,
            __num,
            **args
        )

        return dyna_harm

    return _modes


###################################################################
def extr_matr_elim_lagr(self, matr_asse):
    matr_lagr = matr_asse.toNumpy()  # function EXTR_MATR available in the official source code
    # -----------------------------------------------------#
    # --                                                 --#
    # -- Elimination des degres de libertes de Lagranges --#
    # --                                                 --#
    # --        c.f. doc R4.06.02   - section 2.5        --#
    # --        + commentaires dans sdll123a.comm        --#
    # --                                                 --#
    # -----------------------------------------------------#

    dof_num = matr_asse.getDOFNumbering()
    ind_lag1 = dof_num.getLagrangeDOFs()
    ind_nolag = dof_num.getPhysicalDOFs()

    nlag1 = len(ind_lag1)
    nnolag = len(ind_nolag)

    Z = NP.zeros((nnolag - nlag1, nnolag))
    C = NP.vstack((matr_lagr[ind_lag1][:, ind_nolag], Z))
    Q, R = NP.linalg.qr(NP.transpose(C))

    dR = []
    for i1 in range(len(R)):
        dR.append(NP.abs(R[i1, i1]))

    mdR = NP.max(dR)
    indz = []
    for i1 in range(len(R)):
        if NP.abs(R[i1, i1]) <= mdR * 1.0e-16:
            indz.append(i1)

    matr_sans_lagr = NP.dot(
        NP.transpose(Q[:][:, indz]), NP.dot(matr_lagr[ind_nolag][:, ind_nolag], Q[:][:, indz])
    )

    # -- Fin elimination

    return matr_sans_lagr
