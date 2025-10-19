#ifndef CYCLICSYMMETRYMODE_H_
#define CYCLICSYMMETRYMODE_H_

/**
 * @file CyclicSymmetryMode.h
 * @brief Fichier entete de la classe CyclicSymmetryMode
 * @author Natacha Bereux
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

/* person_in_charge: natacha.bereux at edf.fr */

#include "astercxx.h"

#include "DataStructures/DataStructure.h"
#include "Interfaces/StructureInterface.h"
#include "MemoryManager/JeveuxVector.h"
#include "Meshes/Mesh.h"
#include "Modal/ModalBasis.h"
#include "Supervis/ResultNaming.h"

/**
 * @class CyclicSymmetryMode
 * @brief Cette classe correspond a un MOD_CYCL
 * @author Natacha Bereux
 */
class CyclicSymmetryMode : public DataStructure {
  private:
    /** @brief Objet Jeveux '.CYCL_TYPE */
    JeveuxVectorChar8 _type;
    /** @brief Objet Jeveux '.CYCL_DESC' */
    JeveuxVectorLong _desc;
    /** @brief Objet Jeveux '.CYCL_DIAM' */
    JeveuxVectorLong _diam;
    /** @brief Objet Jeveux '.CYCL_NBSC' */
    JeveuxVectorLong _numberOfSectors;
    /** @brief Objet Jeveux '.CYCL_REFE' */
    JeveuxVectorChar24 _refe;
    /** @brief Objet Jeveux '.CYCL_NUIN' */
    JeveuxVectorLong _interfaceIndices;
    /** @brief Objet Jeveux '.CYCL_CMODE' */
    JeveuxVectorComplex _cMode;
    /** @brief Objet Jeveux .CYCL_FREQ */
    JeveuxVectorComplex _cFreq;
    /** @brief Mesh  */
    MeshPtr _mesh;
    /** @brief Structure Interface */
    StructureInterfacePtr _structInterf;
    /** @brief Modal Basis */
    StandardModalBasisPtr _modalBasis;

  public:
    /**
     * @typedef CyclicSymmetryModePtr
     * @brief Pointeur intelligent vers un CyclicSymmetryMode
     */
    typedef std::shared_ptr< CyclicSymmetryMode > CyclicSymmetryModePtr;

    /**
     * @brief Constructeur
     */
    CyclicSymmetryMode( const std::string name = ResultNaming::getNewResultName() )
        : DataStructure( name, 8, "MODE_CYCL" ),
          _type( JeveuxVectorChar8( getName() + ".CYCL_TYPE" ) ),
          _desc( JeveuxVectorLong( getName() + ".CYCL_DESC" ) ),
          _diam( JeveuxVectorLong( getName() + "CYCL_DIAM" ) ),
          _numberOfSectors( JeveuxVectorLong( getName() + ".CYCL_NBSC" ) ),
          _refe( JeveuxVectorChar24( getName() + ".CYCL_REFE" ) ),
          _cFreq( JeveuxVectorComplex( getName() + ".CYCL_FREQ" ) ),
          _interfaceIndices( JeveuxVectorLong( getName() + ".CYCL_NUIN" ) ),
          _cMode( JeveuxVectorComplex( getName() + ".CYCL_MODE" ) ),
          /** Pointers to objects listed in _refe */
          _mesh( MeshPtr() ),
          _structInterf( StructureInterfacePtr() ),
          _modalBasis( StandardModalBasisPtr() )

              {};

    bool setMesh( MeshPtr &currentMesh ) {
        if ( currentMesh->isEmpty() ) {
            AS_ABORT( "Mesh is empty" );
        }
        _mesh = currentMesh;
        return true;
    };
    bool setStructureInterface( StructureInterfacePtr &currentStructInterf ) {
        _structInterf = currentStructInterf;
        return true;
    };
    bool setModalBasis( StandardModalBasisPtr &currentModalBasis ) {
        _modalBasis = currentModalBasis;
        return true;
    };
};

/**
 * @typedef CyclicSymmetryModePtr
 * @brief Pointeur intelligent vers un CyclicSymmetryMode
 */
typedef std::shared_ptr< CyclicSymmetryMode > CyclicSymmetryModePtr;

/**
 * @brief Definition du maillage
 * @param currentMesh objet MeshPtr sur lequel le Mode Cyclique reposera
 */

#endif /* CyclicSymmetryMode_H_ */
