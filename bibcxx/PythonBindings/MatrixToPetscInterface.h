#ifndef MATRIXTOPETSCINTERFACE_H_
#define MATRIXTOPETSCINTERFACE_H_

/**
 * @file FortranInterface.cxx
 * @brief Python bindings for Fortran interface.
 * @author Mathieu Courtois
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

/* person_in_charge: mathieu.courtois@edf.fr */
#include "aster_pybind.h"

void exportMatrixToPetscToPython( py::module_ &mod );

#endif /* MATRIXTOPETSCINTERFACE_H_ */
