/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org             */
/* This file is part of code_aster.                                     */
/*                                                                      */
/* code_aster is free software: you can redistribute it and/or modify   */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */
/*                                                                      */
/* code_aster is distributed in the hope that it will be useful,        */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */
/*                                                                      */
/* You should have received a copy of the GNU General Public License    */
/* along with code_aster.  If not, see <http://www.gnu.org/licenses/>.  */
/* -------------------------------------------------------------------- */

/* person_in_charge: mathieu.courtois at edf.fr */

#ifndef LOGICAL_UNIT_H_
#define LOGICAL_UNIT_H_

#include "aster.h"
/*
 *   PUBLIC FUNCTIONS
 *
 */

#ifdef __cplusplus
extern "C" {
#endif

extern int openLogicalUnitFile( const char *name, const int type, const int access );
extern void releaseLogicalUnitFile( const int logicalUnit );

/* FIN LOGICAL_UNIT_H */
#ifdef __cplusplus
}
#endif

#endif
