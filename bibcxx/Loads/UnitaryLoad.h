#ifndef UNITARYLOAD_H_
#define UNITARYLOAD_H_

/**
 * @file UnitaryLoad.h
 * @brief Fichier definissant la classe template UnitaryLoad
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

#include "astercxx.h"

#include "Loads/PhysicalQuantity.h"

#include <stdexcept>

/**
 * @class UnitaryLoad
 * @brief Classe template definissant une paire d'un MeshEntity et d'une valeur imposee sur une
 * composante
 * @author Nicolas Sellenet
 */
template < class PhysicalQuantityType >
class UnitaryLoad {
  private:
    /** @typedef Definition d'un pointeur intelligent sur un VirtualMeshEntity */
    typedef std::shared_ptr< VirtualMeshEntity > MeshEntityPtr;

    /** @typedef Definition du type informatique de la grandeur sur laquelle porte la charge */
    typedef typename PhysicalQuantityType::QuantityType ValueType;

    /** @brief MeshEntity sur laquelle repose le "blocage" */
    MeshEntityPtr _meshEntity;
    /** @brief "Numero" de la composante à imposer */
    PhysicalQuantityComponent _loadCoordinate;
    /** @brief Valeur a imposer */
    ValueType _value;

  public:
    /**
     * @brief Constructeur
     * @param meshEntity MeshEntity sur laquelle repose la chargement
     * @param curCoord Coordonnee de la grandeur sur laquelle on impose le chargement
     * @param value Valeur du chargement
     */
    UnitaryLoad( MeshEntityPtr meshEntity, PhysicalQuantityComponent curCoord, ValueType value )
        : _meshEntity( meshEntity ), _loadCoordinate( curCoord ), _value( value ) {
        if ( !PhysicalQuantityType::hasComponent( curCoord ) )
            throw std::runtime_error( ComponentNames.find( curCoord )->second + " not allowed" );
    };

    /**
     * @brief Obtenir le nom de la coordonnee
     * @return Chaine representant le nom Aster de la coordonnee
     */
    const std::string getAsterCoordinateName() const {
        return std::string( ComponentNames.find( _loadCoordinate )->second );
    };

    /**
     * @brief Obtenir l'entite du maillage sur laquelle repose le chargement
     * @return Un pointeur vers l'entite
     */
    const MeshEntityPtr &getMeshEntityPtr() const { return _meshEntity; };

    /**
     * @brief Obtenir la valeur du chargement
     * @return La valeur de type ValueType
     */
    const ValueType getValue() const { return _value; };
};

/** @typedef Definition d'un chargement sur DEPL_R */
typedef UnitaryLoad< DisplacementReal > DisplacementLoadReal;
/** @typedef Definition d'un chargement sur TEMP_R */
typedef UnitaryLoad< TemperatureReal > TemperatureLoadReal;
typedef UnitaryLoad< TemperatureFunction > TemperatureLoadFunction;
/** @typedef Definition d'un chargement pression */
typedef UnitaryLoad< PressureReal > PressureLoadReal;
#endif /* UNITARYLOAD_H_ */
