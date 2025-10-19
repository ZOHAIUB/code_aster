! --------------------------------------------------------------------
! Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! code_aster is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------
!
#include "asterf_types.h"
!
interface
    subroutine nxnewt(model, mateco, caraElem, listLoad, nume_dof, &
                      solver, tpsthe, timeMap, matass, cn2mbr, &
                      maprec, cnchci, varc_curr, temp_prev, temp_iter, &
                      vtempp, vec2nd, mediri, conver, hydr_prev, &
                      hydr_curr, comporTher, cnvabt, &
                      cnresi, ther_crit_i, ther_crit_r, reasma, ds_algorom, &
                      ds_print, sddisc, iter_newt, l_stat)
        use NonLin_Datastructure_type
        use Rom_Datastructure_type
        character(len=8), intent(in) :: model, caraElem
        character(len=24), intent(in) :: mateco
        character(len=24), intent(in) :: listLoad
        character(len=24), intent(in) :: nume_dof
        character(len=19), intent(in) :: solver
        real(kind=8) :: tpsthe(6)
        character(len=24), intent(in) :: timeMap
        character(len=19), intent(in) :: varc_curr, sddisc
        integer(kind=8), intent(in) :: iter_newt
        character(len=24) :: matass
        character(len=19) :: maprec
        character(len=24) :: cnchci
        character(len=24) :: cn2mbr
        character(len=24) :: temp_prev
        character(len=24) :: temp_iter
        character(len=24) :: vtempp
        character(len=24) :: vec2nd
        character(len=24) :: mediri
        aster_logical, intent(out) :: conver
        character(len=24) :: hydr_prev
        character(len=24) :: hydr_curr
        character(len=24) :: comporTher
        character(len=24) :: cnvabt
        character(len=24) :: cnresi
        integer(kind=8) :: ther_crit_i(*)
        real(kind=8) :: ther_crit_r(*)
        aster_logical :: reasma
        type(ROM_DS_AlgoPara), intent(in) :: ds_algorom
        type(NL_DS_Print), intent(inout) :: ds_print
        aster_logical, intent(in) :: l_stat
    end subroutine nxnewt
end interface
