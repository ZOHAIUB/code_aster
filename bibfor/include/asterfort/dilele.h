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
#include "asterfort/Behaviour_type.h"
!
interface
    subroutine dilele(option, typmod, ds_dil, ndim, nnos, &
                      nnom, npg, nddl, dimdef, iw, vff, &
                      vffb, idff, idffb, geomi, compor, &
                      mate, lgpg, carcri, instam, instap, &
                      ddlm, ddld, siefm, vim, &
                      siefp, vip, fint, matr, &
                      lMatr, lVect, lSigm, codret)
        use dil_type
        aster_logical :: lSigm, lMatr, lVect
        character(len=8), intent(in) :: typmod(2)
        character(len=16), intent(in):: option, compor(COMPOR_SIZE)
        type(dil_modelisation) :: ds_dil
        integer(kind=8), intent(in)          :: ndim, nnos, nnom, npg, nddl, lgpg, dimdef
        integer(kind=8), intent(in)          :: mate, iw, idff, idffb
        real(kind=8), intent(in)     :: carcri(CARCRI_SIZE), instam, instap
        real(kind=8), intent(in)     :: geomi(ndim, nnos+nnom)
        real(kind=8), intent(in)     :: vff(nnos+nnom, npg), vffb(nnos, npg)
        real(kind=8), intent(in)     :: ddlm(nddl), ddld(nddl)
        real(kind=8), intent(in)     :: siefm(dimdef*npg), vim(lgpg*npg)
        real(kind=8), intent(inout)  :: siefp(dimdef*npg), vip(lgpg*npg)
        real(kind=8), intent(inout)  :: fint(nddl), matr(nddl, nddl)
        integer(kind=8), intent(inout)         :: codret
    end subroutine dilele
end interface
