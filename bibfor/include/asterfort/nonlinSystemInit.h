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
    subroutine nonlinSystemInit(list_func_acti, nume_dof, ds_algopara, ds_contact, ds_system)
        use NonLin_Datastructure_type
        integer(kind=8), intent(in) :: list_func_acti(*)
        character(len=24), intent(in) :: nume_dof
        type(NL_DS_AlgoPara), intent(in) :: ds_algopara
        type(NL_DS_Contact), intent(in) :: ds_contact
        type(NL_DS_System), intent(inout) :: ds_system
    end subroutine nonlinSystemInit
end interface
