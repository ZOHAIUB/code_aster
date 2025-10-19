
#include "astercxx.h"

#ifdef ASTER_HAVE_MPI

#ifndef PARALLELMECHANICALLOAD_H_
#define PARALLELMECHANICALLOAD_H_

/**
 * @file ParallelMechanicalLoad.h
 * @brief Fichier entete de la classe ParallelMechanicalLoad
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

#include "DataFields/ConstantFieldOnCells.h"
#include "DataStructures/DataStructure.h"
#include "Loads/MechanicalLoad.h"
#include "Meshes/ConnectionMesh.h"
#include "Meshes/MeshExplorer.h"
#include "Modeling/FiniteElementDescriptor.h"
#include "Modeling/Model.h"
#include "Modeling/ParallelFiniteElementDescriptor.h"
#include "ParallelUtilities/AsterMPI.h"
#include "Supervis/ResultNaming.h"
#include "Utilities/SyntaxSaver.h"

#include <fstream>

/**
 * @class ParallelMechanicalLoad
 * @brief Classe definissant une charge dualisée parallèle
 */
template < typename ConstantFieldOnCellsType >
class ParallelMechanicalLoad : public DataStructure {
  public:
    using ConstantFieldOnCellsTypePtr = std::shared_ptr< ConstantFieldOnCellsType >;

    using ParallelMechanicalLoadPtr = std::shared_ptr< ParallelMechanicalLoad >;

  private:
    template < typename ConstantFieldOnCellsType2Ptr >
    void transferConstantFieldOnCells( const ConstantFieldOnCellsType2Ptr &fieldIn,
                                       ConstantFieldOnCellsType2Ptr &fieldOut ) {
        const auto &toKeep = _FEDesc->getVirtualCellsToKeep();

        std::string savedName( "" );
        fieldIn->build();
        fieldOut->allocate( fieldIn );
        const auto sizeFieldIn = ( *fieldIn ).size();
        const auto vect_resu = fieldIn->getValues();

        fieldIn->updateValuePointers();
        for ( int pos = 0; pos < sizeFieldIn; ++pos ) {
            const auto &zone = fieldIn->getZoneDescription( pos );
            const auto &curFEDesc = zone.getFiniteElementDescriptor();
            if ( curFEDesc->getName() != savedName && savedName != "" ) {
                std::string a( "Different FiniteElementDescriptor in one ConstantFieldOnCells is "
                               "not allowed" );
                throw std::runtime_error( a );
            }
            savedName = curFEDesc->getName();

            const auto &listCells = zone.getListOfCells();
            VectorLong toCopy;
            toCopy.reserve( listCells.size() );
            for ( const auto &num : listCells ) {
                if ( toKeep[-num - 1] != 1 )
                    toCopy.push_back( toKeep[-num - 1] );
            }

            if ( toCopy.size() != 0 ) {
                const auto newZone =
                    ConstantFieldOnZone( zone.getFiniteElementDescriptor(), toCopy );
                const auto resu = vect_resu[pos];
                fieldOut->setValueOnZone( newZone, resu );
            }
        }

        fieldOut->build();
    };

    void allocateFields( const MechanicalLoadPtr< ConstantFieldOnCellsType > &load ) {
        auto typeLoadStd = load->getType();
        typeLoadStd->updateValuePointer();

        _type->allocate( 1 );
        ( *_type )[0] = ( *typeLoadStd )[0];
        _modelName->allocate( 1 );
        ( *_modelName )[0] = _model->getName();

        transferConstantFieldOnCells( load->getMechanicalLoadDescription()->getImposedField(),
                                      _cimpo );
        transferConstantFieldOnCells(
            load->getMechanicalLoadDescription()->getMultiplicativeField(), _cmult );
    };

  protected:
    /** @brief Modele */
    ModelPtr _model;
    /** @brief Vecteur Jeveux '.LIGRE' */
    ParallelFiniteElementDescriptorPtr _FEDesc;
    /** @brief Carte '.CIMPO' */
    ConstantFieldOnCellsTypePtr _cimpo;
    /** @brief Carte '.CMULT' */
    ConstantFieldOnCellsRealPtr _cmult;
    /** @brief Vecteur Jeveux '.TYPE' */
    JeveuxVectorChar8 _type;
    /** @brief Vecteur Jeveux '.MODEL.NOMO' */
    JeveuxVectorChar8 _modelName;
    /** @brief Syntax used to build object */
    SyntaxSaverPtr _syntax;
    /** @brief node and cell support groups of ParallelMechnicalLoad */
    VectorString _grpNo, _grpMa;

  public:
    /**
     * @brief Constructeur
     */
    ParallelMechanicalLoad( void ) = delete;

    /**
     *
     * @brief Constructeur
     */
    ParallelMechanicalLoad( const MechanicalLoadPtr< ConstantFieldOnCellsType > &load,
                            const ModelPtr &model )
        : ParallelMechanicalLoad( ResultNaming::getNewResultName(), load, model ) {};

    /** @brief Constructor */
    ParallelMechanicalLoad( const ParallelMechanicalLoadPtr load, const ModelPtr &model )
        : DataStructure( ResultNaming::getNewResultName(), 8, "CHAR_MECA" ),
          _type( getName() + ".TYPE" ),
          _modelName( getName() + ".CHME.MODEL.NOMO" ),
          _model( model ) {
        _syntax = load->_syntax;
        _grpNo = load->_grpNo;
        _grpMa = load->_grpMa;
        auto pMesh = std::static_pointer_cast< ParallelMesh >( model->getMesh() );
        ConnectionMeshPtr connectionMesh( new ConnectionMesh( pMesh, _grpNo, _grpMa ) );
        ModelPtr connectionModel( new Model( connectionMesh ) );
        connectionModel->setFrom( model );
        MechanicalLoadPtr< ConstantFieldOnCellsType > partialMechanicalLoad(
            new MechanicalLoad< ConstantFieldOnCellsType >( connectionModel ) );
        partialMechanicalLoad->buildFromSyntax( _syntax );

        _FEDesc = std::make_shared< ParallelFiniteElementDescriptor >(
            getName() + ".CHME.LIGRE", partialMechanicalLoad->getFiniteElementDescriptor(),
            partialMechanicalLoad->getModel()->getConnectionMesh(), model );
        _cimpo = std::make_shared< ConstantFieldOnCellsType >( getName() + ".CHME.CIMPO", _FEDesc );
        _cmult = std::make_shared< ConstantFieldOnCellsReal >( getName() + ".CHME.CMULT", _FEDesc );

        allocateFields( partialMechanicalLoad );
    };

    /**
     * @brief Constructeur
     */
    ParallelMechanicalLoad( const std::string &name,
                            const MechanicalLoadPtr< ConstantFieldOnCellsType > &load,
                            const ModelPtr &model )
        : DataStructure( name, 8, "CHAR_MECA" ),
          _FEDesc( std::make_shared< ParallelFiniteElementDescriptor >(
              getName() + ".CHME.LIGRE", load->getFiniteElementDescriptor(),
              load->getModel()->getConnectionMesh(), model ) ),
          _cimpo(
              std::make_shared< ConstantFieldOnCellsType >( getName() + ".CHME.CIMPO", _FEDesc ) ),
          _cmult(
              std::make_shared< ConstantFieldOnCellsReal >( getName() + ".CHME.CMULT", _FEDesc ) ),
          _type( getName() + ".TYPE" ),
          _modelName( getName() + ".CHME.MODEL.NOMO" ),
          _model( model ) {
        allocateFields( load );
    };

    /**
     * @brief Constructeur
     */
    ParallelMechanicalLoad( const std::string &name,
                            const ParallelFiniteElementDescriptorPtr &fEDesc,
                            const ModelPtr &model )
        : DataStructure( name, 8, "CHAR_MECA" ),
          _FEDesc( fEDesc ),
          _cimpo(
              std::make_shared< ConstantFieldOnCellsType >( getName() + ".CHME.CIMPO", _FEDesc ) ),
          _cmult(
              std::make_shared< ConstantFieldOnCellsReal >( getName() + ".CHME.CMULT", _FEDesc ) ),
          _type( getName() + ".TYPE" ),
          _modelName( getName() + ".CHME.MODEL.NOMO" ),
          _model( model ) {};

    /**
     * @brief Function membre debugPrint
     * @param logicalUnit Unite logique d'impression
     */
    void debugPrint( const int logicalUnit = 6 ) const {
        this->DataStructure::debugPrint( logicalUnit );
        const auto &explorer = _FEDesc->getVirtualCellsExplorer();
        const auto &mesh = _FEDesc->getMesh();
        auto &LToGmapMesh = mesh->getLocalToGlobalNodeIds();
        auto &LToGmapFE = _FEDesc->getLocalToGlobalMapping();
        auto e1 = LToGmapMesh.exists();
        auto e2 = LToGmapFE.exists();
        if ( !e1 || !e2 )
            return;
        LToGmapMesh->updateValuePointer();
        LToGmapFE->updateValuePointer();

        std::ofstream outFile;
        outFile.open( "fort." + std::to_string( logicalUnit ), std::ios::app );
        outFile << "\nEcriture de la connectivité des mailles fantômes en numérotation globale\n";
        for ( const auto meshElem : explorer ) {
            const auto &numElem = meshElem.getId() + 1;
            outFile << numElem << " : ";
            for ( auto numNode : meshElem ) {
                if ( numNode > 0 ) {
                    outFile << ( *LToGmapMesh )[numNode - 1] << " ";
                } else {
                    outFile << ( *LToGmapFE )[-numNode - 1] << " ";
                }
            }
            outFile << std::endl;
        }
        outFile.close();
    };

    /**
     * @brief Get the finite element descriptor
     */
    ParallelFiniteElementDescriptorPtr getFiniteElementDescriptor() const { return _FEDesc; };

    /**
     * @brief Get the model
     */
    ModelPtr getModel() const {
        if ( ( !_model ) || _model->isEmpty() )
            throw std::runtime_error( "Model of current load is empty" );
        return _model;
    };

    ConstantFieldOnCellsRealPtr getMultiplicativeField() const { return _cmult; };

    ConstantFieldOnCellsTypePtr getImposedField() const { return _cimpo; };

    /**
     * @brief Function to set rebuild parameters (used for balancing)
     * @param syntax aster syntax
     * @param grpNo cell groups
     * @param syntax node groups
     */
    void setRebuildParameters( SyntaxSaverPtr syntax, const VectorString &grpNo,
                               const VectorString &grpMa ) {
        _syntax = syntax;
        _grpNo = grpNo;
        _grpMa = grpMa;
    };
};

/**
 * @typedef ParallelMechanicalLoadPtr
 * @brief Pointeur intelligent vers un ParallelMechanicalLoad
 */
using ParallelMechanicalLoadReal = ParallelMechanicalLoad< ConstantFieldOnCellsReal >;

using ParallelMechanicalLoadFunction = ParallelMechanicalLoad< ConstantFieldOnCellsChar24 >;

using ParallelMechanicalLoadRealPtr = std::shared_ptr< ParallelMechanicalLoadReal >;

using ParallelMechanicalLoadFunctionPtr = std::shared_ptr< ParallelMechanicalLoadFunction >;

/** @typedef std::list de ParallelMechanicalLoad */
using ListParaMecaLoadReal = std::list< ParallelMechanicalLoadRealPtr >;
/** @typedef Iterateur sur une std::list de ParallelMechanicalLoad */
using ListParaMecaLoadRealIter = ListParaMecaLoadReal::iterator;
/** @typedef Iterateur constant sur une std::list de ParallelMechanicalLoad */
using ListParaMecaLoadRealCIter = ListParaMecaLoadReal::const_iterator;

/** @typedef std::list de ParallelMechanicalLoad */
using ListParaMecaLoadFunction = std::list< ParallelMechanicalLoadFunctionPtr >;
/** @typedef Iterateur sur une std::list de ParallelMechanicalLoad */
using ListParaMecaLoadFunctionIter = ListParaMecaLoadFunction::iterator;
/** @typedef Iterateur constant sur une std::list de ParallelMechanicalLoad */
using ListParaMecaLoadFunctionCIter = ListParaMecaLoadFunction::const_iterator;

#endif /* PARALLELMECHANICALLOAD_H_ */

#endif /* ASTER_HAVE_MPI */
