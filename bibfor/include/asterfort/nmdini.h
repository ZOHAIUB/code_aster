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
    subroutine nmdini(factorKeyword, listInstJv, tole, &
                      nbInst, numeInstInit, instInit)
        character(len=16), intent(in) :: factorKeyword
        character(len=19), intent(in) :: listInstJv
        real(kind=8), intent(in) :: tole
        integer(kind=8), intent(in) :: nbInst
        integer(kind=8), intent(out) :: numeInstInit
        real(kind=8), intent(out) :: instInit
    end subroutine nmdini
end interface
