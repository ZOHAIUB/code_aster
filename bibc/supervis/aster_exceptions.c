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

#include "aster.h"

#include "aster_exceptions.h"

#include <setjmp.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * Emulate the behavior of exceptions using the system functions 'setjmp/longjmp'.
 *
 * The fortran subroutines can raise an exception by calling interruptTry via XFINI
 * or UEXCEP (usually UEXCEP if called through UTEXCP or U2MESS subroutines).
 * XFINI is called at the end of a normal execution and raise the EOFError exception.
 *
 * try {                                            if ((gExcNumb = setjmp(env)) == 0) { <--
 *      ...                                                 ...                     |
 *      interruptTry( code );                               longjmp(env, code);   ---
 *      ...                                                 ...
 * }                                --->            }
 * except( code ) {                                 else if (gExcEnv == code ) {
 *      ...                                                 ...
 * }                                                }
 * exceptAll {                                      else {
 *      ...                                                 ...
 * }                                                }
 * endTry();
 *
 * NB: there are two differences with the Python behavior/syntax.
 *     1. There is no finally clause.
 *     2. An additional statement endTry() to decrement the counter level.
 *        Do not forget endTry() if there is a return statement in a block.
 *
 * except( code ) : will probably be not very usefull in C.
 *
 * Global variables:
 *  gExcNumb: code of the exception to raise
 *  gExcEnv : array to store the stack environment
 *  gExcArgs: arguments passed to the exception raised
 *
 */
int gExcLvl = 0;
int gExcNumb = -1;
jmp_buf gExcEnv[NIVMAX + 1];

static PyObject *gExcArgs = NULL;
static PyObject *exc_module = NULL;

/*
 *   PRIVATE FUNCTIONS
 *
 */
void _new_try() {
    /* Begin of try : `gExcLvl` incremented
     */
    gExcLvl += 1;
    if ( NIVMAX < gExcLvl ) {
        printf( "AssertionError: too many nested try/except statements: %d\n", gExcLvl );
        abort();
    }
}

void _end_try() {
    /* End if try : `gExcLvl` is decremented
     */
    gExcLvl -= 1;
}
