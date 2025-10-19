/**
 * @file MPIGroup.cxx
 * @brief Implementation de MPIGroup
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

#include "ParallelUtilities/MPIGroup.h"

#ifdef ASTER_HAVE_MPI
#include "ParallelUtilities/AsterMPI.h"

void MPIGroup::buildFromProcsVector( const MPI_Comm &parentComm, const VectorInt &procIdVector ) {
    MPI_Comm_group( parentComm, &_parentGroup );
    int nbProc = procIdVector.size();
    MPI_Comm_size( parentComm, &nbProc );

    const int size = procIdVector.size();
    if ( size < nbProc ) {
        MPI_Group_incl( _parentGroup, size, procIdVector.data(), &_currentGroup );
        MPI_Comm_create_group( parentComm, _currentGroup, 0, &_groupComm );
        groupCreated = true;
        newGroupCreated = true;
        auto commWorld = aster_get_comm_world();
        asterComm.id = _groupComm;
        asterComm.parent = commWorld;
        asterComm.level = 1;
        asterComm.nbchild = 0;
        AS_ASSERT( commWorld->nbchild < MAX_CHILDS );
        commWorld->childs[commWorld->nbchild] = &asterComm;
        commWorld->nbchild++;
    } else if ( size == nbProc ) {
        _currentGroup = _parentGroup;
        _groupComm = parentComm;
        groupCreated = true;
        newGroupCreated = false;
    } else {
        throw std::runtime_error( "MPI parent group must be bigger than child one" );
    }
};

MPI_Comm MPIGroup::getCommunicator() const {
    if ( groupCreated )
        return _groupComm;
    else
        throw std::runtime_error( "MPI group not created" );
};

MPIGroup::~MPIGroup() {
    if ( newGroupCreated ) {
        int i = 0, j;
        aster_comm_t *node = &asterComm;
        auto parent = node->parent;
        auto nb = parent->nbchild;
        while ( i < nb && parent->childs[i] != node ) {
            i++;
        }
        AS_ASSERT( i < nb );
        for ( j = i + 1; j < nb; j++ ) {
            parent->childs[j - 1] = parent->childs[j];
        }
        parent->childs[nb - 1] = NULL;
        parent->nbchild = nb - 1;

        MPI_Group_free( &_currentGroup );
    }
    groupCreated = false;
    newGroupCreated = false;
};

#endif /* ASTER_HAVE_MPI */
