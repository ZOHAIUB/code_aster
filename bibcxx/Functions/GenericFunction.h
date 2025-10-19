#ifndef GENERICFUNCTION_H_
#define GENERICFUNCTION_H_

/**
 * @file GenericFunction.h
 * @brief Fichier entete de la classe GenericFunction
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

/**
 * @class GenericFunction
 * @brief Base class of Function, Formula and Table
 * @author Nicolas Sellenet
 */
class GenericFunction : public DataStructure {
  protected:
    // Vecteur Jeveux '.PROL'
    JeveuxVectorChar24 _property;
    // Type of Function
    std::string _funct_type;

    void propertyAllocate() {
        // Create Jeveux vector ".PROL"
        _property->allocate( 6 );
        ( *_property )[0] = _funct_type;
        ( *_property )[1] = "LIN LIN";
        ( *_property )[2] = "";
        ( *_property )[3] = "TOUTRESU";
        ( *_property )[4] = "EE";
        ( *_property )[5] = getName();
    };

  public:
    /**
     * @typedef GenericFunctionPtr
     * @brief Pointeur intelligent vers un GenericFunction
     */
    typedef std::shared_ptr< GenericFunction > GenericFunctionPtr;

    /**
     * @brief Constructeur
     */
    GenericFunction( const std::string &name, const std::string &type,
                     const std::string &functType );

    /**
     * @brief Allocate function
     */
    virtual void allocate( ASTERINTEGER size ) {};

    /**
     * @brief Return the properties of the function
     * @return vector of strings
     */
    VectorString getProperties() const;

    /**
     * @brief Get the result name
     * @return  name of the result
     */
    virtual std::string getResultName() { return ""; };

    /**
     * @brief Return the number of points of the function
     */
    virtual ASTERINTEGER maximumSize() const { return 0; };

    /**
     * @brief Definition of the type of extrapolation
     * @param extrapolation type of extrapolation
     * @type  extrapolation string
     * @todo checking
     */
    void setExtrapolation( const std::string type );

    /**
     * @brief Return the number of points of the function
     */
    virtual ASTERINTEGER size() const { return 0; };
};

/**
 * @typedef GenericFunctionPtr
 * @brief Pointeur intelligent vers un GenericFunction
 */
typedef std::shared_ptr< GenericFunction > GenericFunctionPtr;

#endif /* GENERICFUNCTION_H_ */
