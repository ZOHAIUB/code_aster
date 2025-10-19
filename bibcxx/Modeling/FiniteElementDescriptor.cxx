/**
 * @file FiniteElementDescriptor.cxx
 * @brief Implementation de FiniteElementDescriptor
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

#include "Modeling/FiniteElementDescriptor.h"

#include "aster_fort_ds.h"

#include "Meshes/BaseMesh.h"
#include "Meshes/ConnectionMesh.h"
#include "Modeling/PhysicalQuantityManager.h"
#include "Modeling/PhysicsAndModelings.h"
#include "ParallelUtilities/AsterMPI.h"
#include "Utilities/Tools.h"

FiniteElementDescriptor::FiniteElementDescriptor( const std::string &name, const std::string &type,
                                                  const BaseMeshPtr mesh )
    : DataStructure( name, 19, type ),
      _numberOfDelayedNumberedConstraintNodes( getName() + ".NBNO" ),
      _parameters( getName() + ".LGRF" ),
      _dofDescriptor( getName() + ".PRNM" ),
      _listOfGroupsOfElements( getName() + ".LIEL" ),
      _groupsOfCellsNumberByElement( getName() + ".REPE" ),
      _virtualCellsDescriptor( getName() + ".NEMA" ),
      _dofOfDelayedNumberedConstraintNodes( getName() + ".PRNS" ),
      _virtualNodesNumbering( getName() + ".LGNS" ),
      _superElementsDescriptor( getName() + ".SSSA" ),
      _nameOfNeighborhoodStructure( getName() + ".NVGE" ),
      _typeOfCells( getName() + ".TYFE" ),
      _mesh( mesh ),
      _partition( nullptr ),
      _explorer(
          FiniteElementDescriptor::ConnectivityVirtualCellsExplorer( _virtualCellsDescriptor ) ),
      _explorer2(
          FiniteElementDescriptor::ConnectivityVirtualCellsExplorer( _listOfGroupsOfElements ) ) {
    if ( _parameters->exists() ) {
        _parameters->updateValuePointer();
        if ( strip( ( *_parameters )[1] ) != "" ) {
            _partition = std::make_shared< Partition >( ( *_parameters )[1] );
        }
    };
};

FiniteElementDescriptor::FiniteElementDescriptor( const std::string &name, const BaseMeshPtr mesh )
    : FiniteElementDescriptor( name, "LIGREL", mesh ) {};

FiniteElementDescriptor::FiniteElementDescriptor( const BaseMeshPtr mesh )
    : FiniteElementDescriptor( DataStructureNaming::getNewName(), mesh ) {};

FiniteElementDescriptor::FiniteElementDescriptor( const FiniteElementDescriptorPtr FEDesc,
                                                  const VectorString &groupOfCells )
    : FiniteElementDescriptor( FEDesc->getMesh() ) {

    VectorLong listOfCells = _mesh->getCells( groupOfCells );

    std::string base( "G" );
    ASTERINTEGER nbCells = listOfCells.size();
    for ( auto &cell : listOfCells )
        cell += 1;
    CALL_EXLIM2( listOfCells.data(), &nbCells, FEDesc->getName(), base, getName() );
    this->build();
};

FiniteElementDescriptor::FiniteElementDescriptorPtr
    FiniteElementDescriptor::restrict( const VectorString &groupsOfCells ) const {

    VectorLong listOfCells = _mesh->getCells( groupsOfCells );

    return this->restrict( listOfCells );
};

FiniteElementDescriptor::FiniteElementDescriptorPtr
    FiniteElementDescriptor::restrict( const VectorLong &cells ) const {

    auto fed = std::make_shared< FiniteElementDescriptor >( getMesh() );

    VectorLong listOfCells = cells;

    std::string base( "G" );
    ASTERINTEGER nbCells = listOfCells.size();
    for ( auto &cell : listOfCells )
        cell += 1;
    CALL_EXLIM2( listOfCells.data(), &nbCells, getName(), base, fed->getName() );
    fed->build();
    return fed;
};

const FiniteElementDescriptor::ConnectivityVirtualCellsExplorer &
FiniteElementDescriptor::getVirtualCellsExplorer() const {
    _virtualCellsDescriptor->build();
    return _explorer;
};

const JeveuxVectorLong &FiniteElementDescriptor::getVirtualNodesComponentDescriptor() const {
    return _dofOfDelayedNumberedConstraintNodes;
};

const JeveuxVectorLong &FiniteElementDescriptor::getVirtualNodesNumbering() const {
    _virtualNodesNumbering->updateValuePointer();
    return _virtualNodesNumbering;
};

const FiniteElementDescriptor::ConnectivityVirtualCellsExplorer &
FiniteElementDescriptor::getListOfGroupsOfElementsExplorer() const {
    _listOfGroupsOfElements->build();
    return _explorer2;
};

const JeveuxContiguousCollectionLong &FiniteElementDescriptor::getListOfGroupsOfElements() const {
    return _listOfGroupsOfElements;
};

const JeveuxContiguousCollectionLong &FiniteElementDescriptor::getVirtualCellsDescriptor() const {
    return _virtualCellsDescriptor;
};

ASTERINTEGER FiniteElementDescriptor::getNumberOfVirtualNodes() const {
    _numberOfDelayedNumberedConstraintNodes->updateValuePointer();
    return ( *_numberOfDelayedNumberedConstraintNodes )[0];
};

void FiniteElementDescriptor::setNumberOfVirtualNodes( const ASTERINTEGER nbNodes ) {

    if ( !_numberOfDelayedNumberedConstraintNodes.exists() ) {
        _numberOfDelayedNumberedConstraintNodes->allocate( 1 );
    } else {
        _numberOfDelayedNumberedConstraintNodes->updateValuePointer();
    }
    ( *_numberOfDelayedNumberedConstraintNodes )[0] = nbNodes;
};

JeveuxVectorLong FiniteElementDescriptor::getNumberOfVirtualNodesDescriptor() {
    return _numberOfDelayedNumberedConstraintNodes;
};

JeveuxVectorChar8 FiniteElementDescriptor::getParameters() const { return _parameters; };

const JeveuxVectorLong &FiniteElementDescriptor::getPhysicalNodesComponentDescriptor() const {
    _dofDescriptor->updateValuePointer();
    return _dofDescriptor;
};

const JeveuxVectorLong &FiniteElementDescriptor::getListOfGroupsOfElementsbyElement() const {
    _groupsOfCellsNumberByElement->updateValuePointer();
    return _groupsOfCellsNumberByElement;
};

JeveuxVectorLong FiniteElementDescriptor::getFiniteElementType() const { return _typeOfCells; }

ASTERINTEGER FiniteElementDescriptor::getNumberOfCells() const {
    _groupsOfCellsNumberByElement->updateValuePointer();

    auto size = _groupsOfCellsNumberByElement->size() / 2;
    ASTERINTEGER nbCells = 0;

    for ( int i = 0; i < size; i++ ) {
        if ( ( *_groupsOfCellsNumberByElement )[2 * i] != 0 &&
             ( *_groupsOfCellsNumberByElement )[2 * i + 1] != 0 ) {
            nbCells += 1;
        }
    }

    return nbCells;
};

const BaseMeshPtr FiniteElementDescriptor::getMesh() const { return _mesh; };

void FiniteElementDescriptor::setMesh( const BaseMeshPtr &currentMesh ) { _mesh = currentMesh; };

int FiniteElementDescriptor::getPhysics( void ) const {
    const std::string docu = strip( _parameters->getInformationParameter() );

    if ( docu == "MECA" )
        return Physics::Mechanics;
    else if ( docu == "THER" )
        return Physics::Thermal;
    else if ( docu == "ACOU" )
        return Physics::Acoustic;
    else
        throw std::runtime_error( "Unknown physics" );

    return -1;
};

ASTERINTEGER FiniteElementDescriptor::numberOfSuperElement() {
    const std::string typeco( "LIGREL" );
    ASTERINTEGER repi = 0, ier = 0;
    JeveuxChar32 repk( " " );
    const std::string arret( "C" );
    const std::string questi( "NB_SS_ACTI" );

    CALLO_DISMOI( questi, getName(), typeco, &repi, repk, arret, &ier );

    return repi;
};

bool FiniteElementDescriptor::existsFiniteElement() {
    const std::string typeco( "LIGREL" );
    ASTERINTEGER repi = 0, ier = 0;
    JeveuxChar32 repk( " " );
    const std::string arret( "C" );
    const std::string questi( "EXI_ELEM" );

    CALLO_DISMOI( questi, getName(), typeco, &repi, repk, arret, &ier );
    auto retour = strip( repk.toString() );
    if ( retour == "OUI" )
        return true;
    return false;
};

bool FiniteElementDescriptor::existsSuperElement() { return ( numberOfSuperElement() > 0 ); }

bool FiniteElementDescriptor::exists() const {
    if ( _parameters.exists() && _numberOfDelayedNumberedConstraintNodes.exists() )
        return true;

    return false;
};

bool FiniteElementDescriptor::build() {
    // too costly in affe_char_meca
    // _virtualCellsDescriptor->build();
    _listOfGroupsOfElements->build();

    return true;
};

ASTERINTEGER FiniteElementDescriptor::getElemTypeNume( const std::string elemTypeName ) const {

    ASTERINTEGER elemTypeNume;
    JeveuxChar32 objName( " " );
    std::string name = "&CATA.TE.NOMTE";
    CALLO_JEXNOM( objName, name, elemTypeName );
    CALLO_JENONU( objName, &elemTypeNume );
    return elemTypeNume;
};

const std::string FiniteElementDescriptor::getPartitionMethod() const {
    if ( _partition ) {
        return _partition->getMethod();
    }

    return ModelSplitingMethodNames[(int)Centralized];
};

#ifdef ASTER_HAVE_MPI
void FiniteElementDescriptor::transferDofDescriptorFrom( FiniteElementDescriptorPtr &other ) {
    // "the mesh associated to finiteElementDescriptor is not a partial mesh"
    AS_ASSERT( getMesh()->isConnection() );
    const ConnectionMeshPtr connectionMesh =
        std::static_pointer_cast< ConnectionMesh >( getMesh() );

    // "parallel mesh associated to partial mesh of FiniteElementDescriptor \n"
    //        "does not correspond to other FiniteElementDescriptor mesh"
    AS_ASSERT( connectionMesh->getParallelMesh() == other->getMesh() );

    const JeveuxVectorLong &otherDofDescriptor = other->getPhysicalNodesComponentDescriptor();

    const int rank = getMPIRank();
    const int size = getMPISize();
    int nbNodes = connectionMesh->getNumberOfNodes();
    int nec = otherDofDescriptor->size() / other->getMesh()->getNumberOfNodes();

    const JeveuxVectorLong &localNumbering = connectionMesh->getNodesLocalNumbering();
    const JeveuxVectorLong &owner = connectionMesh->getNodesOwner();

    int nbNodesLoc = 0;
    for ( int i = 0; i < nbNodes; ++i ) {
        if ( ( *owner )[i] == rank )
            nbNodesLoc++;
    }

    VectorLong buffer( nec * nbNodesLoc, 0 );
    nbNodesLoc = 0;
    for ( int i = 0; i < nbNodes; ++i ) {
        int nodeNum = ( *localNumbering )[i] - 1;
        if ( ( *owner )[i] == rank ) {
            for ( int j = 0; j < nec; ++j )
                buffer[nbNodesLoc * nec + j] = ( *otherDofDescriptor )[nodeNum * nec + j];

            nbNodesLoc++;
        }
    }

    std::vector< VectorLong > gathered;
    AsterMPI::all_gather( buffer, gathered );
    buffer.clear();

    _dofDescriptor->allocate( nbNodes * nec );
    VectorLong nbNodesProc( size, 0 );
    for ( int i = 0; i < nbNodes; ++i ) {
        auto rowner = ( *owner )[i];
        for ( int j = 0; j < nec; ++j )
            ( *_dofDescriptor )[i * nec + j] = gathered[rowner][nbNodesProc[rowner] * nec + j];
        nbNodesProc[rowner]++;
    }
};

void FiniteElementDescriptor::transferListOfGroupOfCellFrom( FiniteElementDescriptorPtr &other ) {
    // "the mesh associated to finiteElementDescriptor is not a partial mesh"
    AS_ASSERT( getMesh()->isConnection() );
    const ConnectionMeshPtr connectionMesh =
        std::static_pointer_cast< ConnectionMesh >( getMesh() );

    // "parallel mesh associated to partial mesh of FiniteElementDescriptor \n"
    //        "does not correspond to other FiniteElementDescriptor mesh"
    AS_ASSERT( connectionMesh->getParallelMesh() == other->getMesh() );

    const int rank = getMPIRank();
    const int size = getMPISize();

    auto &otherRepe = other->getListOfGroupsOfElementsbyElement();
    auto &otherLiel = other->getListOfGroupsOfElements();

    const auto nbCells = connectionMesh->getNumberOfCells();
    auto &cellsLocNum = connectionMesh->getCellsLocalNumbering();
    auto &cellsOwner = connectionMesh->getCellsOwner();

    int nbCellsLoc = 0;
    for ( int i = 0; i < nbCells; ++i ) {
        if ( ( *cellsOwner )[i] == rank )
            nbCellsLoc++;
    }

    VectorLong typeCellFE;
    typeCellFE.reserve( nbCellsLoc );
    for ( int i = 0; i < nbCells; ++i ) {
        if ( ( *cellsOwner )[i] == rank ) {
            auto cellId = ( *cellsLocNum )[i] - 1;
            auto numGrel = ( *otherRepe )[2 * cellId];
            if ( numGrel > 0 ) {
                auto &grel = ( *otherLiel )[numGrel];
                grel->updateValuePointer();
                auto typeFE = ( *grel )[grel->size() - 1];
                typeCellFE.push_back( typeFE );
            } else {
                typeCellFE.push_back( 0 );
            }
        }
    }

    std::vector< VectorLong > typeFEGathered;
    AsterMPI::all_gather( typeCellFE, typeFEGathered );
    typeCellFE.clear();

    std::vector< VectorLong > listOfGrel;
    MapLong listOfGrelNume, listOfGrelNumeInv;

    _groupsOfCellsNumberByElement->allocate( 2 * nbCells );

    VectorLong nbCellsProc( size, 0 );
    for ( int i = 0; i < nbCells; ++i ) {
        auto rowner = ( *cellsOwner )[i];
        auto typeEF = typeFEGathered[rowner][nbCellsProc[rowner]];
        nbCellsProc[rowner]++;

        if ( typeEF > 0 ) {
            if ( listOfGrelNume.count( typeEF ) == 0 ) {
                listOfGrelNume[typeEF] = listOfGrelNume.size();
                listOfGrelNumeInv[listOfGrelNumeInv.size()] = typeEF;
                listOfGrel.push_back( VectorLong() );
                listOfGrel[listOfGrelNume[typeEF]].reserve( nbCells );
            }

            listOfGrel[listOfGrelNume[typeEF]].push_back( i + 1 );

            ( *_groupsOfCellsNumberByElement )[2 * i] = listOfGrelNume[typeEF] + 1;
            ( *_groupsOfCellsNumberByElement )[2 * i + 1] =
                listOfGrel[listOfGrelNume[typeEF]].size();
        } else {
            ( *_groupsOfCellsNumberByElement )[2 * i] = 0;
            ( *_groupsOfCellsNumberByElement )[2 * i + 1] = 0;
        }
    }

    int totalSize = listOfGrelNume.size();
    int pos = 0;
    for ( auto &grel : listOfGrel ) {
        // il faut rajouter le nom de l'EF à la fin
        grel.push_back( listOfGrelNumeInv[pos++] );
        totalSize += grel.size();
    }

    _listOfGroupsOfElements->allocate( listOfGrel.size(), totalSize );
    for ( auto &grel : listOfGrel ) {
        _listOfGroupsOfElements->push_back( grel );
    }
};

void FiniteElementDescriptor::setFrom( FiniteElementDescriptorPtr &other ) {
    // "the mesh associated to finiteElementDescriptor is not a partial mesh"
    AS_ASSERT( getMesh()->isConnection() );
    const ConnectionMeshPtr connectionMesh =
        std::static_pointer_cast< ConnectionMesh >( getMesh() );

    // "parallel mesh associated to partial mesh of FiniteElementDescriptor \n"
    //        "does not correspond to other FiniteElementDescriptor mesh"
    AS_ASSERT( connectionMesh->getParallelMesh() == other->getMesh() );

    // Fill '.LGRF'
    _parameters->allocate( 3 );
    ( *_parameters )[0] = getMesh()->getName();
    _parameters->setInformationParameter( other->getParameters()->getInformationParameter() );

    // Fill '.NBNO'
    _numberOfDelayedNumberedConstraintNodes->allocate( 1 );
    ( *_numberOfDelayedNumberedConstraintNodes )[0] = 0;

    // Fill '.PRNM'
    transferDofDescriptorFrom( other );

    // Fill 'LIEL'
    transferListOfGroupOfCellFrom( other );

    //   tranfer .TYFE
    const auto nbCells = connectionMesh->getNumberOfCells();
    _typeOfCells->allocate( nbCells );
    auto &cellsLocNum = connectionMesh->getCellsLocalNumbering();
    auto &cellsOwner = connectionMesh->getCellsOwner();

    const auto rank = getMPIRank();
    const auto size = getMPISize();

    int nbCellsLoc = 0;
    for ( int i = 0; i < nbCells; ++i ) {
        if ( ( *cellsOwner )[i] == rank )
            nbCellsLoc++;
    }

    auto typeCellsOther = other->getFiniteElementType();
    typeCellsOther->updateValuePointer();

    VectorLong buffer;
    buffer.reserve( nbCellsLoc );
    for ( int i = 0; i < nbCells; ++i ) {
        if ( ( *cellsOwner )[i] == rank )
            buffer.push_back( ( *typeCellsOther )[( *cellsLocNum )[i] - 1] );
    }

    std::vector< VectorLong > gathered;
    AsterMPI::all_gather( buffer, gathered );

    VectorLong nbCellsProc( size, 0 );
    for ( int i = 0; i < nbCells; ++i ) {
        auto rowner = ( *cellsOwner )[i];
        ( *_typeOfCells )[i] = gathered[rowner][nbCellsProc[rowner]];
        nbCellsProc[rowner]++;
    }
}

#endif /* ASTER_HAVE_MPI */
