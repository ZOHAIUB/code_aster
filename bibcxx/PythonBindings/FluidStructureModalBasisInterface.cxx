/**
 * @file FluidStructureModalBasisInterface.cxx
 * @brief Interface python de FluidStructureModalBasis
 * @author Natacha Béreux
 * @section LICENCE
 *   Copyright (C) 1991 - 2023  EDF R&D                www.code-aster.org
 *
 *   This file is part of Code_Aster.
 *
 *   Code_Aster is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Code_Aster is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Code_Aster.  If not, see <http://www.gnu.org/licenses/>.
 */

/* person_in_charge: natacha.bereux at edf.fr */

#include "PythonBindings/FluidStructureModalBasisInterface.h"

#include "aster_pybind.h"

void exportFluidStructureModalBasisToPython( py::module_ &mod ) {

    py::class_< FluidStructureModalBasis, FluidStructureModalBasis::FluidStructureModalBasisPtr,
                DataStructure >( mod, "FluidStructureModalBasis" )
        .def( py::init( &initFactoryPtr< FluidStructureModalBasis > ) )
        .def( py::init( &initFactoryPtr< FluidStructureModalBasis, std::string > ) )

        .def( "getTable", &ListOfTables::getTable, R"(
Extract a Table from the datastructure.

Arguments:
    identifier (str): Table identifier.

Returns:
    Table: Table stored with the given identifier.
        )",
              py::arg( "identifier" ) );
};
