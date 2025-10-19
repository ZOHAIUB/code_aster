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
    subroutine nmplgs(ndim, nno1, vff1, idfde1, nno2,&
                      vff2, idfde2, npg, iw, geom,&
                      typmod, option, mate, compor, carcri,&
                      instam, instap, angmas, ddlm, ddld,&
                      sigm, lgpg, vim, sigp, vip,&
                      matr, vect, codret, livois,&
                      nbvois, numa, lisoco, nbsoco,&
                      lVari, lSigm, lMatr, lVect)
        integer(kind=8) :: lgpg
        integer(kind=8) :: npg
        integer(kind=8) :: nno2
        integer(kind=8) :: nno1
        integer(kind=8) :: ndim
        real(kind=8) :: vff1(nno1, npg)
        integer(kind=8) :: idfde1
        real(kind=8) :: vff2(nno2, npg)
        integer(kind=8) :: idfde2
        integer(kind=8) :: iw
        real(kind=8) :: geom(ndim, nno1)
        character(len=8) :: typmod(*)
        character(len=16) :: option
        integer(kind=8) :: mate
        character(len=16) :: compor(*)
        real(kind=8) :: carcri(*)
        real(kind=8) :: instam
        real(kind=8) :: instap
        real(kind=8) :: angmas(3)
        real(kind=8) :: ddlm(*)
        real(kind=8) :: ddld(*)
        real(kind=8) :: sigm(2*ndim, npg)
        real(kind=8) :: vim(lgpg, npg)
        real(kind=8) :: sigp(2*ndim, npg)
        real(kind=8) :: vip(lgpg, npg)
        real(kind=8) :: matr(*)
        real(kind=8) :: vect(*)
        integer(kind=8) :: codret
        integer(kind=8), parameter :: nvoima=12, nscoma=4
        integer(kind=8) :: livois(1:nvoima)
        integer(kind=8) :: nbvois
        integer(kind=8) :: numa
        integer(kind=8) :: lisoco(1:nvoima, 1:nscoma, 1:2)
        integer(kind=8) :: nbsoco(1:nvoima)
        aster_logical, intent(in) :: lVari, lSigm, lMatr, lVect
    end subroutine nmplgs
end interface
