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
! aslint: disable=C1505
#include "asterfort/Behaviour_type.h"

interface
    subroutine lc0058(BEHinteg, &
                      fami, kpg, ksp, ndim, typmod, &
                      imate, compor, carcri, instam, instap, &
                      neps, epsm, deps, nsig, sigm, &
                      nvi, vim, option, angmas, &
                      sigp, vip, ndsde, dsidep, codret)
        use Behaviour_type
        type(Behaviour_Integ), intent(in) :: BEHinteg
        character(len=*), intent(in) :: fami
        integer(kind=8), intent(in) :: kpg, ksp, ndim
        character(len=8), intent(in) :: typmod(*)
        integer(kind=8), intent(in) :: imate
        character(len=16), intent(in) :: compor(COMPOR_SIZE)
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        real(kind=8), intent(in) :: instam, instap
        integer(kind=8), intent(in) :: neps
        real(kind=8), intent(in) :: epsm(neps), deps(neps)
        integer(kind=8), intent(in) :: nsig
        real(kind=8), intent(in) :: sigm(nsig)
        integer(kind=8), intent(in) :: nvi
        real(kind=8), intent(in) :: vim(nvi)
        character(len=16), intent(in) :: option
        real(kind=8), intent(in) :: angmas(*)
        real(kind=8), intent(out) :: sigp(nsig)
        real(kind=8), intent(out) :: vip(nvi)
        integer(kind=8), intent(in) :: ndsde
        real(kind=8), intent(out) :: dsidep(merge(nsig, 6, nsig*neps .eq. ndsde), &
                                            merge(neps, 6, nsig*neps .eq. ndsde))
        integer(kind=8), intent(out) :: codret
    end subroutine lc0058
end interface
