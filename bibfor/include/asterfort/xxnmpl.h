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
#include "asterf_types.h"
!
interface
    subroutine xxnmpl(elrefp, elrese, ndim, coorse, igeom,&
                      he, nfh, ddlc, ddlm, nfe,&
                      instam, instap, ideplp, sigm, vip,&
                      basloc, nnop, npg, typmod, option,&
                      imate, compor, lgpg, carcri, idepl,&
                      lsn, lst, idecpg, sig, vi,&
                      matuu, ivectu, codret, nfiss, heavn, jstno,&
                      lMatr, lVect, lSigm)
        aster_logical, intent(in) :: lMatr, lVect, lSigm
        integer(kind=8) :: nfiss
        integer(kind=8) :: lgpg
        integer(kind=8) :: npg
        integer(kind=8) :: nnop
        integer(kind=8) :: nfe
        integer(kind=8) :: nfh
        integer(kind=8) :: ndim
        character(len=8) :: elrefp
        character(len=8) :: elrese
        real(kind=8) :: coorse(*)
        integer(kind=8) :: igeom
        real(kind=8) :: he(nfiss)
        integer(kind=8) :: ddlc
        integer(kind=8) :: ddlm
        real(kind=8) :: instam
        real(kind=8) :: instap
        integer(kind=8) :: ideplp
        real(kind=8) :: sigm(2*ndim, npg)
        real(kind=8) :: vip(lgpg, npg)
        real(kind=8) :: basloc(3*ndim*nnop)
        character(len=8) :: typmod(*)
        character(len=16) :: option
        integer(kind=8) :: imate
        character(len=16) :: compor(*)
        real(kind=8) :: carcri(*)
        integer(kind=8) :: idepl
        real(kind=8) :: lsn(nnop)
        real(kind=8) :: lst(nnop)
        integer(kind=8) :: idecpg
        real(kind=8) :: sig(2*ndim, npg)
        real(kind=8) :: vi(lgpg, npg)
        real(kind=8) :: matuu(*)
        integer(kind=8) :: ivectu
        integer(kind=8) :: codret
        integer(kind=8) :: heavn(nnop, 5)
        integer(kind=8) :: jstno
    end subroutine xxnmpl
end interface
