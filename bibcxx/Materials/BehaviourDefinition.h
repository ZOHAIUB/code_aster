#ifndef BEHAVIOURDEFINITION_H_
#define BEHAVIOURDEFINITION_H_

/**
 * @file BehaviourDefinition.h
 * @brief Fichier entete de la classe BehaviourDefinition
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

#include "astercxx.h"

#include "DataStructures/DataStructure.h"
#include "MemoryManager/JeveuxVector.h"
#include "Supervis/ResultNaming.h"

/**
 * @class BehaviourDefinition
 * @brief produit une sd identique a celle produite par DEFI_COMPOR
 * @author Nicolas Sellenet
 */
class BehaviourDefinition : public DataStructure {
  private:
    /** @brief Objet '.CPRK' */
    JeveuxVectorChar24 _cprk;
    /** @brief Objet '.CPRR' */
    JeveuxVectorReal _cprr;
    /** @brief Objet '.CPRI' */
    JeveuxVectorLong _cpri;

  public:
    /**
     * @typedef BehaviourDefinitionPtr
     * @brief Pointeur intelligent vers un BehaviourDefinition
     */
    typedef std::shared_ptr< BehaviourDefinition > BehaviourDefinitionPtr;

    /**
     * @brief Constructeur
     */
    BehaviourDefinition() : BehaviourDefinition( ResultNaming::getNewResultName() ) {};

    /**
     * @brief Constructeur
     */
    BehaviourDefinition( const std::string &name )
        : DataStructure( name, 8, "COMPOR" ),
          _cprk( JeveuxVectorChar24( getName() + ".CPRK" ) ),
          _cprr( JeveuxVectorReal( getName() + ".CPRR" ) ),
          _cpri( JeveuxVectorLong( getName() + ".CPRI" ) ) {};
};

/**
 * @typedef BehaviourDefinitionPtr
 * @brief Pointeur intelligent vers un BehaviourDefinition
 */
typedef std::shared_ptr< BehaviourDefinition > BehaviourDefinitionPtr;

#endif /* BEHAVIOURDEFINITION_H_ */
