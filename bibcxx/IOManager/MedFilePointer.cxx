/**
 * @file MedProfle.cxx
 * @brief Implementation de MedFilePointer
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include "IOManager/MedFilePointer.h"

#include "ParallelUtilities/AsterMPI.h"

#ifdef ASTER_HAVE_MED

int MedFilePointer::close() {
    MEDfileClose( _fileId );
    _fileId = -1;
    _isOpen = false;
    return 0;
};

med_idt MedFilePointer::getFileId() const {
    if ( !_isOpen )
        throw std::runtime_error( "Med file not open" );
    return _fileId;
};

int MedFilePointer::open( const std::filesystem::path &filename,
                          const MedFileAccessType &openType ) {
    med_access_mode medAccessMode = MED_ACC_UNDEF;
    if ( openType == MedReadOnly ) {
        medAccessMode = MED_ACC_RDONLY;
    } else if ( openType == MedReadWrite ) {
        medAccessMode = MED_ACC_RDEXT;
    } else if ( openType == MedCreate ) {
        medAccessMode = MED_ACC_CREAT;
    } else {
        throw std::runtime_error( "Med file access type not allowed" );
    }
    _fileId = MEDfileOpen( filename.string().c_str(), medAccessMode );
    _isOpen = true;
    _parallelOpen = false;
    return 0;
};

int MedFilePointer::openParallel( const std::filesystem::path &filename,
                                  const MedFileAccessType &openType ) {
#ifdef ASTER_HAVE_MPI
    MPI_Info info = MPI_INFO_NULL;
    MPI_Comm comm = aster_get_comm_world()->id;
    med_access_mode medAccessMode = MED_ACC_UNDEF;
    if ( openType == MedReadOnly ) {
        medAccessMode = MED_ACC_RDONLY;
    } else if ( openType == MedReadWrite ) {
        medAccessMode = MED_ACC_RDEXT;
    } else if ( openType == MedCreate ) {
        medAccessMode = MED_ACC_CREAT;
    } else {
        throw std::runtime_error( "Med file access type not allowed" );
    }
    _fileId = MEDparFileOpen( filename.string().c_str(), medAccessMode, comm, MPI_INFO_NULL );
    _isOpen = true;
    _parallelOpen = true;
    return 0;
#else
    throw std::runtime_error( "Parallel opening not available in sequential" );
    return 0;
#endif /* ASTER_HAVE_MPI */
};
#endif
