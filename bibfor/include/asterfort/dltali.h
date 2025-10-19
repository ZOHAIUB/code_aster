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
    subroutine dltali(neq, result, imat, masse, rigid,&
                      liad, lifo, nchar, nveca, lcrea,&
                      lprem, lamort, t0, mate, mateco, carele,&
                      charge, infoch, fomult, modele, numedd,&
                      nume, solveu, criter, dep0, vit0,&
                      acc0, fexte0, famor0, fliai0, &
                      tabwk, force0, force1, ds_energy, kineLoad)
        use NonLin_Datastructure_type
        integer(kind=8) :: neq
        integer(kind=8) :: imat(3)
        character(len=8) :: masse
        character(len=8) :: rigid
        integer(kind=8) :: liad(*)
        character(len=24) :: lifo(*)
        integer(kind=8) :: nchar
        integer(kind=8) :: nveca
        aster_logical :: lcrea
        aster_logical :: lprem
        aster_logical :: lamort
        real(kind=8) :: t0
        character(len=24) :: mate, mateco, kineLoad
        character(len=24) :: carele
        character(len=24) :: charge
        character(len=24) :: infoch
        character(len=24) :: fomult
        character(len=24) :: modele
        character(len=24) :: numedd
        integer(kind=8) :: nume
        character(len=19) :: solveu
        character(len=24) :: criter
        real(kind=8) :: dep0(*)
        real(kind=8) :: vit0(*)
        real(kind=8) :: acc0(*)
        real(kind=8) :: fexte0(*)
        real(kind=8) :: famor0(*)
        real(kind=8) :: fliai0(*)
        real(kind=8) :: tabwk(*)
        character(len=19) :: force0
        character(len=19) :: force1
        character(len=8), intent(in) :: result
        type(NL_DS_Energy), intent(out) :: ds_energy
    end subroutine dltali
end interface
