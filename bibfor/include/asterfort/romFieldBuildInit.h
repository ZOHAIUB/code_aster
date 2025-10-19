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
    subroutine romFieldBuildInit(mesh         , nbNodeMesh , listNode       ,&
                                 nbFieldResult, resultField, resultFieldNume,&
                                 resultRom    , modelRom   , tablReduCoor   ,&
                                 fieldBuild)
        use Rom_Datastructure_type
        character(len=8), intent(in) :: mesh
        integer(kind=8), intent(in) :: nbNodeMesh
        integer(kind=8), pointer  :: listNode(:)
        integer(kind=8), intent(in)  :: nbFieldResult
        character(len=16), pointer :: resultField(:)
        integer(kind=8), pointer :: resultFieldNume(:)
        type(ROM_DS_Result), intent(in) :: resultRom
        character(len=8), intent(in) :: modelRom
        type(ROM_DS_TablReduCoor), intent(in) :: tablReduCoor
        type(ROM_DS_FieldBuild), intent(inout) :: fieldBuild
    end subroutine romFieldBuildInit
end interface
