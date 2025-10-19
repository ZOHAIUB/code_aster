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

interface
    subroutine matass2petsc(matasz, local, petscMatz, iret)
        character(len=*), intent(in) :: matasz
        integer(kind=8), intent(in) :: local
#ifdef ASTER_HAVE_PETSC
        PetscErrorCode, intent(out) :: iret
        Mat, intent(out) :: petscMatz
#else
        integer(kind=8), intent(out) :: petscMatz, iret
#endif
    end subroutine matass2petsc
end interface
