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
!
interface
    subroutine cvmcvg(dy, ddy, nr, itmax, toler,&
                      iter, intg, typess, essai, icomp,&
                      irteti)
        real(kind=8) :: dy(*)
        real(kind=8) :: ddy(*)
        integer(kind=8) :: nr
        integer(kind=8) :: itmax
        real(kind=8) :: toler
        integer(kind=8) :: iter
        integer(kind=8) :: intg
        integer(kind=8) :: typess
        real(kind=8) :: essai
        integer(kind=8) :: icomp
        integer(kind=8) :: irteti
    end subroutine cvmcvg
end interface
