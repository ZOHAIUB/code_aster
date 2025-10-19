# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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

import aster
from ...Messages import UTMESS

from ...Cata.Syntax import _F
from ...CodeCommands import DEFI_BASE_MODALE, EXTR_MODE
from .dyna_visco_modes_calc import dyna_visco_modes_calc


def dyna_visco_modes(
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
):
    """
    Macro-command DYNA_VISCO,
    function to compute the eigenmodes of the structure
    """

    i = 0
    nmode = 5

    # EIGENMODES COMPUTATION

    freq1 = list_FREQ[0]

    # search for the 1st eigenfrequency
    j = 0
    [_modes, freq1, nmode] = dyna_visco_modes_calc(
        self,
        TYPE_MODE,
        freq1,
        nmode,
        RESI_RELA,
        i,
        j,
        MATER_ELAS_FO,
        e0,
        eta0,
        __asseMg,
        __asseKgr,
        __asseKg,
        __listKv,
        trKg,
        ltrv,
        TYPE_RESU,
        **args
    )

    # search for following eigenfrequencies
    while freq1 <= fmax:
        j = j + 1

        [_modes, freq1, nmode] = dyna_visco_modes_calc(
            self,
            TYPE_MODE,
            freq1,
            nmode,
            RESI_RELA,
            i,
            j,
            MATER_ELAS_FO,
            e0,
            eta0,
            __asseMg,
            __asseKgr,
            __asseKg,
            __listKv,
            trKg,
            ltrv,
            TYPE_RESU,
            reuse="oui",
            co_reuse=_modes,
            **args
        )

    if TYPE_MODE in ["REEL", "BETA_REEL"]:
        __modes_temp = EXTR_MODE(FILTRE_MODE=_F(MODE=_modes, TOUT_ORDRE="OUI"))
        _modes = DEFI_BASE_MODALE(ORTHO_BASE=_F(BASE=__modes_temp, MATRICE=__asseMg))

    #########################################################
    # PRINTING OF THE EIGENMODES
    UTMESS("I", "DYNAVISCO_3")
    eigenfreq = aster.GetResu(_modes.getName(), "VARI_ACCES")["FREQ"]

    if TYPE_MODE in ["REEL", "BETA_REEL"]:
        UTMESS("I", "DYNAVISCO_4")
    if TYPE_MODE == "COMPLEXE":
        eigendamping = aster.GetResu(_modes.getName(), "PARAMETRES")["AMOR_REDUIT"]
        UTMESS("I", "DYNAVISCO_6")
    for k in range(0, len(eigenfreq)):
        if TYPE_MODE in ["REEL", "BETA_REEL"]:
            UTMESS("I", "DYNAVISCO_5", vali=k + 1, valr=eigenfreq[k])
        if TYPE_MODE == "COMPLEXE":
            UTMESS("I", "DYNAVISCO_7", vali=k + 1, valr=(eigenfreq[k], eigendamping[k]))

    return _modes
