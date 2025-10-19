/**
 * @file DataStructureNaming.cxx
 * @brief Implementation de DataStructureNaming
 * @author Nicolas Sellenet
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
/* person_in_charge: nicolas.sellenet at edf.fr */

#include "DataStructures/DataStructureNaming.h"

#include "DataStructures/TemporaryDataStructureNaming.h"
#include "MemoryManager/JeveuxObject.h"
#include "Supervis/ResultNaming.h"

#include <sstream>
#include <string>

std::string DataStructureNaming::getNewName( JeveuxMemory memoryType, const int lengthName ) {
    if ( memoryType == Permanent ) {
        std::string tmpName = ResultNaming::getNewResultName();
        return std::string( tmpName + "                        ", 0, lengthName );
    } else if ( memoryType == Temporary )
        return TemporaryDataStructureNaming::getNewTemporaryName( lengthName );
    else
        throw std::runtime_error( "Programming error" );
    return "";
};
