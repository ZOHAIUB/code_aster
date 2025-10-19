#ifndef XFEMCRACK_H_
#define XFEMCRACK_H_

/**
 * @file XfemCrack.h
 * @brief Fichier entete de la classe XfemCrack
 * @author Nicolas Tardieu
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

/* person_in_charge: nicolas.tardieu at edf.fr */

#include "astercxx.h"

#include "Crack/CrackShape.h"
#include "DataFields/FieldOnNodes.h"
#include "DataFields/ListOfTables.h"
#include "DataStructures/DataStructure.h"
#include "Functions/Function.h"
#include "Meshes/Mesh.h"
#include "Modeling/Model.h"

/**
 * @class XfemCrack
 * @brief generates a data structure identical to DEFI_FISS_XFEM
 * @author Nicolas Tardieu
 */
class XfemCrack : public DataStructure, public ListOfTables {
  public:
    /**
     * @brief kind of forward declaration
     */
    typedef std::shared_ptr< XfemCrack > XfemCrackPtr;

  private:
    /** @typedef Definition of a smart pointer on VirtualMeshEntity */
    typedef std::shared_ptr< VirtualMeshEntity > MeshEntityPtr;

    /** @brief Mesh supporting the crack */
    MeshPtr _mesh;
    /** @brief Mesh of the auxiliary grid */
    MeshPtr _auxiliaryGrid;
    /** @brief Mesh of the previous auxiliary grid */
    XfemCrackPtr _existingCrackWithGrid;
    /** @brief Type de discontinuité */
    std::string _discontinuityType;
    /** @brief MeshEntity defining the lips of the crack */
    VectorString _crackLipsEntity;
    /** @brief MeshEntity defining the tip of the crack */
    VectorString _crackTipEntity;
    /** @brief Function defining the normal level set */
    FunctionPtr _normalLevelSetFunction;
    /** @brief Function defining the tangential level set */
    FunctionPtr _tangentialLevelSetFunction;
    /** @brief Crack Shape */
    CrackShapePtr _crackShape;
    /** @brief List of group of cells that define the crack tip in case of propagation in the
     * cohesive case */
    VectorString _cohesiveCrackTipForPropagation;
    /** @brief Field defining the normal level set */
    FieldOnNodesRealPtr _tangentialLevelSetField;
    /** @brief MeshEntity defining the enriched elements */
    FieldOnNodesRealPtr _normalLevelSetField;
    /** @brief Field defining the tangential level set */
    VectorString _enrichedCells;
    /** @brief Name of the discontinuous field */
    std::string _discontinuousField;
    /** @brief Type of the enrichment */
    std::string _enrichmentType;
    /** @brief Radius of th enriched zone */
    ASTERDOUBLE _enrichmentRadiusZone;
    /** @brief Number of enriched elements layers */
    ASTERINTEGER _enrichedLayersNumber;
    /** @brief List of juncting cracks */
    std::vector< XfemCrackPtr > _junctingCracks;
    /** @brief Point to define juncting cracks */
    VectorReal _pointForJunctingCracks;

    /** @brief Nodal Field named '.GRLTNO' */
    FieldOnNodesRealPtr _tangentialLevelSetGradient;
    /** @brief Nodal Field  named '.GRLNNO' */
    FieldOnNodesRealPtr _normalLevelSetGradient;
    /** @brief Nodal Field named '.BASLOC' */
    FieldOnNodesRealPtr _localBasis;
    /** @brief JeveuxVectorReal  named '.FONDFISS' */
    JeveuxVectorReal _crackTipCoords;
    /** @brief JeveuxVectorReal  named '.BASEFOND' */
    JeveuxVectorReal _crackTipBasis;
    /** @brief JeveuxVectorLong  named '.FONDMULT' */
    JeveuxVectorLong _crackTipMultiplicity;
    /** @brief JeveuxVectorReal  named '.CARAFOND' */
    JeveuxVectorReal _crackTipCharacteristics;
    /** @brief JeveuxVectorReal  named '.FOND.TAILLE_R' */
    JeveuxVectorReal _elementSize;
    /** @brief JeveuxVectorLong  named '.GROUP_NO_ENRI' */
    JeveuxVectorLong _enrichedNodes;
    /** @brief JeveuxVectorLong  named '.MAILFISS.CTIP' */
    JeveuxVectorLong _crackTipCells;
    /** @brief JeveuxVectorLong  named '.MAILFISS.HEAV' */
    JeveuxVectorLong _heavisideCells;
    /** @brief JeveuxVectorLong  named '.MAILFISS.HECT' */
    JeveuxVectorLong _crackTipAndHeavisideCells;
    /** @brief Field defining the status nodes */
    FieldOnNodesLongPtr _nodeStatusField;
    /** @brief JeveuxVectorLong  named '.MAILFISS.MAFOND' */
    JeveuxVectorLong _crackTipCellGroup;
    /** @brief JeveuxVectorLong  named '.GROUP_MA_ENRI' */
    JeveuxVectorLong _enrichedCellGroup;
    /** @brief JeveuxVectorChar8 named '.MAILLAGE' */
    JeveuxVectorChar8 _meshName;
    /** @brief JeveuxVectorChar8 named '.INFO' */
    JeveuxVectorChar16 _info;
    /** @brief JeveuxVectorLong  named '.NOFACPTFON' */
    JeveuxVectorLong _crackTipNodeFacesField;

  public:
    /**
     * @brief Constructeur
     */
    XfemCrack( MeshPtr mesh );

    /**
     * @brief Constructeur
     */
    XfemCrack( const std::string name, MeshPtr mesh );

    /**
     * @brief Construction du XfemCrack
     * @return Booleen indiquant que la construction s'est bien deroulee
     */
    bool build();

    /**
     * @brief Enrichissement d'un ModelPtr avec la fissure XFEM
     * @param baseModel modèle à enrichir
     * @return Modèle enrichi
     */
    ModelPtr enrichModelWithXfem( ModelPtr &baseModel );

    /**
     * @brief Series of getters and setters that check the syntax rules
     */

    const MeshPtr getMesh() const { return _mesh; };

    void setMesh( const MeshPtr &mesh ) { _mesh = mesh; }

    MeshPtr getAuxiliaryGrid() const { return _auxiliaryGrid; }

    void setAuxiliaryGrid( const MeshPtr &auxiliaryGrid ) {
        if ( _existingCrackWithGrid )
            throw std::runtime_error( "Cannot define auxiliary grid with already "
                                      "assigned Crack definition " );
        _auxiliaryGrid = auxiliaryGrid;
    }

    const XfemCrackPtr getExistingCrackWithGrid() const { return _existingCrackWithGrid; }

    void setExistingCrackWithGrid( const XfemCrackPtr &existingCrackWithGrid ) {
        if ( _auxiliaryGrid )
            throw std::runtime_error( "Cannot define existing Crack definition "
                                      "with already assigned auxiliary grid" );
        _existingCrackWithGrid = existingCrackWithGrid;
    }

    std::string getDiscontinuityType() const { return _discontinuityType; }

    void setDiscontinuityType( const std::string &discontinuityType ) {
        if ( discontinuityType != "Crack" && discontinuityType != "Interface" &&
             discontinuityType != "Cohesive" ) {
            throw std::runtime_error( "discontinuityType can only be Crack, Interface, Cohesive" );
        } else {
            _discontinuityType = discontinuityType;
        }
    }

    VectorString getCrackLipsEntity() const { return _crackLipsEntity; }

    void setCrackLipsEntity( const VectorString &crackLips ) {
        if ( !( _crackShape && _normalLevelSetField && _normalLevelSetFunction ) ) {
            _crackLipsEntity = crackLips;
        } else {
            throw std::runtime_error( "Only one to be assigned into crackShape, "
                                      "normalLevelSetFunction, crackLips, normalLevelSetField" );
        }
    }

    VectorString getCrackTipEntity() const { return _crackTipEntity; }

    void setCrackTipEntity( const VectorString &crackTip ) { _crackTipEntity = crackTip; }

    VectorString getCohesiveCrackTipForPropagation() const {
        return _cohesiveCrackTipForPropagation;
    }

    void setCohesiveCrackTipForPropagation( const VectorString &crackTipEntity ) {
        _cohesiveCrackTipForPropagation = crackTipEntity;
    }

    FunctionPtr getNormalLevelSetFunction() const { return _normalLevelSetFunction; }

    void setNormalLevelSetFunction( const FunctionPtr &normalLevelSetFunction ) {
        if ( !( _crackShape && _normalLevelSetField && _crackLipsEntity.size() != 0 ) ) {
            _normalLevelSetFunction = normalLevelSetFunction;
        } else {
            throw std::runtime_error( "Only one to be assigned into crackShape, "
                                      "normalLevelSetFunction, crackLips, normalLevelSetField" );
        }
    }

    FunctionPtr getTangentialLevelSetFunction() const { return _tangentialLevelSetFunction; }

    void setTangentialLevelSetFunction( const FunctionPtr &tangentialLevelSetFunction ) {
        _tangentialLevelSetFunction = tangentialLevelSetFunction;
    }

    CrackShapePtr getCrackShape() const { return _crackShape; }

    void setCrackShape( const CrackShapePtr &crackShape ) {
        if ( !( _normalLevelSetFunction && _normalLevelSetField &&
                _crackLipsEntity.size() != 0 ) ) {
            _crackShape = crackShape;
        } else {
            throw std::runtime_error( "Only one to be assigned into crackShape, "
                                      "normalLevelSetFunction, crackLips, normalLevelSetField" );
        }
    }

    FieldOnNodesRealPtr getNormalLevelSetField() const {
        _normalLevelSetField->updateValuePointers();
        return _normalLevelSetField;
    }

    void setNormalLevelSetField( const FieldOnNodesRealPtr &normalLevelSetField ) {
        if ( !( _crackShape && _crackLipsEntity.size() != 0 && _normalLevelSetFunction ) ) {
            _normalLevelSetField = normalLevelSetField;
        } else {
            throw std::runtime_error( "Only one to be assigned into crackShape, "
                                      "normalLevelSetFunction, crackLips, normalLevelSetField" );
        }
    }

    FieldOnNodesRealPtr getTangentialLevelSetField() const {
        _tangentialLevelSetField->updateValuePointers();
        return _tangentialLevelSetField;
    }

    void setTangentialLevelSetField( const FieldOnNodesRealPtr &tangentialLevelSet ) {
        _tangentialLevelSetField = tangentialLevelSet;
    }

    VectorString getEnrichedCells() const { return _enrichedCells; }

    void setEnrichedCells( const VectorString &enrichedCells ) { _enrichedCells = enrichedCells; }

    std::string getDiscontinuousField() const { return _discontinuousField; }

    void setDiscontinuousField( const std::string &discontinuousField ) {
        _discontinuousField = discontinuousField;
    }

    std::string getEnrichmentType() const { return _enrichmentType; }

    void setEnrichmentType( const std::string &enrichmentType ) {
        _enrichmentType = enrichmentType;
    }

    ASTERDOUBLE getEnrichmentRadiusZone() const { return _enrichmentRadiusZone; }

    void setEnrichmentRadiusZone( ASTERDOUBLE enrichmentRadiusZone ) {
        _enrichmentRadiusZone = enrichmentRadiusZone;
    }

    ASTERINTEGER getEnrichedLayersNumber() const { return _enrichedLayersNumber; }

    void setEnrichedLayersNumber( long enrichedLayersNumber ) {
        _enrichedLayersNumber = enrichedLayersNumber;
    }

    std::vector< XfemCrackPtr > getJunctingCracks() const { return _junctingCracks; }

    void insertJunctingCracks( const XfemCrackPtr &junctingCracks ) {
        _junctingCracks.push_back( junctingCracks );
    }

    void setPointForJunction( const VectorReal point ) { _pointForJunctingCracks = point; }

    const JeveuxVectorReal getCrackTipCoords();

    const JeveuxVectorReal getCrackTipBasis();

    const JeveuxVectorLong getCrackTipMultiplicity();

    const std::string getTipType();

    const JeveuxVectorLong getCrackTipNodeFacesField();

    const JeveuxVectorReal getCrackFrontRadius();
};

/**
 * @typedef MaterialPtr
 * @brief Pointeur intelligent vers un XfemCrack
 */
typedef std::shared_ptr< XfemCrack > XfemCrackPtr;

#endif /* XFEMCRACK_H_ */
