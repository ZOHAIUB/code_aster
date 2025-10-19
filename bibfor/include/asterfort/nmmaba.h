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
    subroutine nmmaba(icodma, rela_comp, e, dsde, sigy,&
                      ncstpm, cstpm)
        integer(kind=8) :: ncstpm
        integer(kind=8) :: icodma
        character(len=16) :: rela_comp
        real(kind=8) :: e
        real(kind=8) :: dsde
        real(kind=8) :: sigy
        real(kind=8) :: cstpm(ncstpm)
    end subroutine nmmaba
end interface
