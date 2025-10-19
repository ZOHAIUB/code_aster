#ifndef CPPTOFORTRANGLOSSARY_H_
#define CPPTOFORTRANGLOSSARY_H_

/**
 * @file CppToFortranGlossary.h
 * @brief Fichier entete de la struct CppToFortranGlossary
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

#include "astercxx.h"

#include "Loads/PhysicalQuantity.h"
#include "Modeling/PhysicsAndModelings.h"

#include <stdexcept>

/**
 * @class Glossary
 * @brief Classe de correspondance entre les noms Aster vers des types C++
 * @author Nicolas Sellenet
 */
class Glossary {
  private:
    /**
     * @typedef std::map d'une chaine et d'un entier
     * @todo Ce sera à changer quand on se passera des enums pour PhysicalQuantity, etc.
     */
    typedef std::map< std::string, int > MapStrInt;
    typedef MapStrInt::const_iterator MapStrIntIter;

    MapStrInt _strToInt;

  public:
    /**
     * @brief Constructeur vide
     */
    Glossary();

    /**
     * @brief Destructeur vide
     */
    ~Glossary() {};

    /**
     * @brief getPhysics
     * @param searchPhysics Nom d'un nom de composante dans le fichier de commande
     * @return une valeur dans l'enum PhysicalQuantityComponent
     */
    PhysicalQuantityComponent getComponent( std::string searchComp ) const {
        MapStrIntIter curIter = _strToInt.find( searchComp );
        if ( curIter == _strToInt.end() )
            throw std::runtime_error( "Unknown component" );
        return (PhysicalQuantityComponent)( curIter->second );
    };

    /**
     * @brief getModeling
     * @param searchMod Nom d'une physique dans le fichier de commande
     * @return une valeur dans l'enum Modelings
     */
    Modelings getModeling( std::string searchMod ) const {
        MapStrIntIter curIter = _strToInt.find( searchMod );
        if ( curIter == _strToInt.end() )
            throw std::runtime_error( "Unknown modeling" );
        return (Modelings)( curIter->second );
    };

    /**
     * @brief getPhysics
     * @param searchPhysics Nom d'une physique dans le fichier de commande
     * @return une valeur dans l'enum Physics
     */
    Physics getPhysics( std::string searchPhysics ) const {
        MapStrIntIter curIter = _strToInt.find( searchPhysics );
        if ( curIter == _strToInt.end() )
            throw std::runtime_error( "Unknown physics" );
        return (Physics)( curIter->second );
    };
};

extern Glossary *fortranGlossary2;

Glossary *getGlossary();

const Glossary &getReferenceToGlossary();

#endif /* CPPTOFORTRANGLOSSARY_H_ */
