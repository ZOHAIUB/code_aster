/**
 * @file GeneralizedAssemblyVectorInterface.cxx
 * @brief Interface python de GeneralizedAssemblyVector
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2024  EDF R&D                www.code-aster.org
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include "PythonBindings/GeneralizedAssemblyVectorInterface.h"

#include "aster_pybind.h"

void exportGeneralizedAssemblyVectorToPython( py::module_ &mod ) {

    py::class_< GenericGeneralizedAssemblyVector, GenericGeneralizedAssemblyVectorPtr,
                DataStructure >( mod, "GeneralizedAssemblyVector" );
    // fake initFactoryPtr: created by subclasses
    // fake initFactoryPtr: created by subclasses

    py::class_< GeneralizedAssemblyVectorReal, GeneralizedAssemblyVectorRealPtr,
                GenericGeneralizedAssemblyVector >( mod, "GeneralizedAssemblyVectorReal" )
        .def( py::init( &initFactoryPtr< GeneralizedAssemblyVectorReal > ) )
        .def( py::init( &initFactoryPtr< GeneralizedAssemblyVectorReal, std::string > ) )
        .def( "getValues",
              py::overload_cast<>( &GeneralizedAssemblyVectorReal::getValues, py::const_ ),
              R"(
            Return a list of values as (x1, y1, z1, x2, y2, z2...)

            Returns:
                list[float]: List of values.
            )" )
        .def( "setValues", &GeneralizedAssemblyVectorReal::setValues, R"(
            Set values of vector.
        
            Arguments:
                values (list[float]): set vector.
            )" );

    py::class_< GeneralizedAssemblyVectorComplex, GeneralizedAssemblyVectorComplexPtr,
                GenericGeneralizedAssemblyVector >( mod, "GeneralizedAssemblyVectorComplex" )
        .def( py::init( &initFactoryPtr< GeneralizedAssemblyVectorComplex > ) )
        .def( py::init( &initFactoryPtr< GeneralizedAssemblyVectorComplex, std::string > ) )
        .def( "getValues",
              py::overload_cast<>( &GeneralizedAssemblyVectorComplex::getValues, py::const_ ),
              R"(
            Return a list of values as (x1, y1, z1, x2, y2, z2...)

            Returns:
                list[complex]: List of values.
            )" )
        .def( "setValues", &GeneralizedAssemblyVectorComplex::setValues, R"(
            Set values of vector.
        
            Arguments:
                values (list[complex]): set vector.
            )" );
};
