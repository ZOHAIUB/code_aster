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
    subroutine ndxprm(modelz, ds_material, carele    , ds_constitutive, ds_algopara   ,&
                      lischa, numedd, solveu , ds_system     ,sddisc,&
                      sddyna, nlDynaDamping,&
                      ds_measure, nume_inst      , list_func_acti,&
                      valinc, solalg     , meelem    , measse     ,&
                      maprec, matass     , faccvg    , ldccvg)
        use NonLin_Datastructure_type
        use NonLinearDyna_type
        type(NL_DS_AlgoPara), intent(in) :: ds_algopara
        integer(kind=8), intent(in) :: list_func_acti(*), nume_inst
        character(len=*) :: modelz
        type(NL_DS_Material), intent(in) :: ds_material
        character(len=24) :: carele
        type(NL_DS_Measure), intent(inout) :: ds_measure
        character(len=24) :: numedd
        type(NL_DS_Constitutive), intent(in) :: ds_constitutive
        type(NL_DS_System), intent(in) :: ds_system
        character(len=19) :: sddisc, lischa, solveu
        character(len=19), intent(in) :: sddyna
        type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
        character(len=19) :: solalg(*), valinc(*)
        character(len=19) :: meelem(*), measse(*)
        character(len=19) :: maprec, matass
        integer(kind=8) :: faccvg, ldccvg
    end subroutine ndxprm
end interface
