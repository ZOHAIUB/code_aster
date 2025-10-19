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
    subroutine compMecaChckStrain(iComp, &
                                  model, fullElemField, &
                                  lAllCellAffe, cellAffe, nbCellAffe, &
                                  lTotalStrain, lMfront, exteDefo, &
                                  defoComp, defoCompPY, &
                                  relaComp, relaCompPY)
        integer(kind=8), intent(in) :: iComp
        character(len=8), intent(in) :: model
        character(len=19), intent(in) :: fullElemField
        aster_logical, intent(in) :: lAllCellAffe
        character(len=24), intent(in) :: cellAffe
        integer(kind=8), intent(in) :: nbCellAffe
        aster_logical, intent(in) :: lTotalStrain, lMfront
        integer(kind=8), intent(in) :: exteDefo
        character(len=16), intent(in) :: defoComp, defoCompPY
        character(len=16), intent(in) :: relaComp, relaCompPY
    end subroutine compMecaChckStrain
end interface
