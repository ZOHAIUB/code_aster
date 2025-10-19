/**
 * @file Result.cxx
 * @brief Implementation de Result
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2025  EDF R&D                www.code-aster.org
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

#include "DataFields/FieldBuilder.h"

std::set< std::string > FieldBuilder::_setGlobNume;
std::set< std::string > FieldBuilder::_setLigrel;

EquationNumberingPtr FieldBuilder::newEquationNumbering( const std::string &name,
                                                         const BaseMeshPtr mesh ) {
    if ( _setGlobNume.count( strip( name ) ) > 0 ) {
        raiseAsterError( "NUME_EQUA already exists: " + name );
    }

    EquationNumberingPtr curDesc;
    if ( mesh->isParallel() ) {
#ifdef ASTER_HAVE_MPI
        curDesc = std::make_shared< ParallelEquationNumbering >( name );
#endif
    } else {
        curDesc = std::make_shared< EquationNumbering >( name );
    }

    curDesc->setMesh( mesh );
    addEquationNumbering( curDesc );

    return curDesc;
};
