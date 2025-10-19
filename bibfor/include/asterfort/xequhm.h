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
interface 
    subroutine xequhm(ds_thm,&
                      imate, option, ta, ta1, ndim,&
                      kpi, npg, dimenr,&
                      enrmec, dimdef, dimcon, nbvari, defgem,&
                      congem, vintm, defgep, congep, vintp,&
                      mecani, press1, press2, tempe,&
                      rinstp, dt, r, drds,&
                      dsde, retcom, angmas, enrhyd, nfh)
        use THM_type
        type(THM_DS), intent(inout) :: ds_thm
        integer(kind=8) :: nbvari
        integer(kind=8) :: dimcon
        integer(kind=8) :: dimdef
        integer(kind=8) :: dimenr
        integer(kind=8) :: imate
        character(len=16) :: option
        real(kind=8) :: ta
        real(kind=8) :: ta1
        integer(kind=8) :: ndim
        integer(kind=8) :: kpi
        integer(kind=8) :: npg
        integer(kind=8) :: enrmec(3)
        real(kind=8) :: defgem(dimdef)
        real(kind=8) :: congem(dimcon)
        real(kind=8) :: vintm(nbvari)
        real(kind=8) :: defgep(dimdef)
        real(kind=8) :: congep(dimcon)
        real(kind=8) :: vintp(nbvari)
        integer(kind=8) :: mecani(5)
        integer(kind=8) :: press1(7)
        integer(kind=8) :: press2(7)
        integer(kind=8) :: tempe(5)
        real(kind=8) :: rinstp
        real(kind=8) :: dt
        real(kind=8) :: r(dimenr)
        real(kind=8) :: drds(dimenr, dimcon)
        real(kind=8) :: dsde(dimcon, dimenr)
        integer(kind=8) :: retcom
        real(kind=8) :: angmas(3)
        integer(kind=8) :: enrhyd(3)
        integer(kind=8) :: nfh
    end subroutine xequhm
end interface 
