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
!
interface
    subroutine reliem(mo, ma, typem, motfaz, iocc,&
                      nbmocl, limocl, tymocl, litroz, nbtrou, l_keep_propz, l_allz)
        integer(kind=8) :: nbmocl
        character(len=*) :: mo
        character(len=8) :: ma
        character(len=*) :: typem
        character(len=*) :: motfaz
        integer(kind=8) :: iocc
        character(len=*) :: limocl(nbmocl)
        character(len=*) :: tymocl(nbmocl)
        character(len=*) :: litroz
        integer(kind=8) :: nbtrou
        aster_logical, optional, intent(in) :: l_keep_propz
    aster_logical, optional, intent(in) :: l_allz
    end subroutine reliem
end interface
