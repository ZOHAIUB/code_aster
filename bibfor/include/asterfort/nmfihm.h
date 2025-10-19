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
    subroutine nmfihm(ndim, nddl, nno1, nno2, npg,&
                      lgpg, ipg, wref, vff1, vff2,&
                      idf2, dffr2, mate, option, geom,&
                      ddlm, ddld, iu, ip, sigm,&
                      sigp, vect, matr, vim, vip,&
                      tm, tp, carcri, compor, typmod,&
                      lVect, lMatr, lSigm, codret)
        integer(kind=8) :: lgpg
        integer(kind=8) :: npg
        integer(kind=8) :: nno2
        integer(kind=8) :: nno1
        integer(kind=8) :: nddl
        integer(kind=8) :: ndim
        integer(kind=8) :: ipg
        real(kind=8) :: wref(npg)
        real(kind=8) :: vff1(nno1, npg)
        real(kind=8) :: vff2(nno2, npg)
        integer(kind=8) :: idf2
        real(kind=8) :: dffr2(ndim-1, nno2, npg)
        integer(kind=8) :: mate
        real(kind=8) :: geom(ndim, nno2)
        real(kind=8) :: ddlm(nddl)
        real(kind=8) :: ddld(nddl)
        integer(kind=8) :: iu(3, 16)
        integer(kind=8) :: ip(8)
        real(kind=8) :: sigm(2*ndim-1, npg)
        real(kind=8) :: sigp(2*ndim-1, npg)
        real(kind=8) :: vect(nddl)
        real(kind=8) :: matr(nddl*nddl)
        real(kind=8) :: vim(lgpg, npg)
        real(kind=8) :: vip(lgpg, npg)
        real(kind=8) :: tm
        real(kind=8) :: tp
        character(len=16), intent(in) :: option
        aster_logical, intent(in) :: lSigm, lVect, lMatr
        real(kind=8), intent(in) :: carcri(*)
        character(len=16), intent(in) :: compor(*)
        character(len=8), intent(in) :: typmod(*)
        integer(kind=8), intent(out) :: codret
    end subroutine nmfihm
end interface
