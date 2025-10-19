/**
 * @file StructureInterface.cxx
 * @brief
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

#include "Interfaces/StructureInterface.h"

#include "aster_fort_superv.h"

#include "Supervis/CommandSyntax.h"

const std::vector< InterfaceTypeEnum > allInterfaceType = { MacNeal, CraigBampton,
                                                            HarmonicCraigBampton, NoInterfaceType };
const VectorString allInterfaceTypeNames = { "MNEAL", "CRAIGB", "CB_HARMO", "AUCUN" };

bool StructureInterface::build() {
    CommandSyntax cmdSt( "DEFI_INTERF_DYNA" );
    cmdSt.setResult( getName(), "INTERF_DYNA_CLAS" );

    CapyConvertibleSyntax syntax;
    syntax.setSimpleKeywordValues( _container );
    for ( const auto &iter : _interfDefs )
        syntax.addCapyConvertibleContainer( iter._container );

    cmdSt.define( syntax );

    ASTERINTEGER op = 98;
    CALL_EXECOP( &op );

    _isBuilt = true;
    return true;
};
