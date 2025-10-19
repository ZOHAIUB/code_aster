/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org             */
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
/* ------------------------------------------------------------------ */
#include "aster.h"

#include "aster_module.h"

#include "aster_core_module.h"
#include "aster_exceptions.h"
#include "aster_fort_ds.h"
#include "aster_fort_jeveux.h"
#include "aster_fort_material.h"
#include "aster_fort_mesh.h"
#include "aster_fort_python.h"
#include "aster_fort_superv.h"
#include "aster_fort_utils.h"
#include "aster_mpi.h"
#include "aster_utils.h"
#include "shared_vars.h"

#include <ctype.h>
#include <signal.h>
#include <stdlib.h>

#ifdef ASTER_HAVE_PETSC
#include "petsc.h"

void charToPetscArgcArgv( char *buffer, char *delim, char **Output, int *index ) {

    int partcount = 0;
    *index = 0;
    Output[partcount++] = "petsc_aster";
    *index += 1;
    Output[partcount++] = buffer;

    char *ptr = buffer;

    while ( ptr != 0 ) { // check if the string is over
        ptr = strstr( ptr, delim );
        if ( ptr != NULL ) {
            *ptr = 0;
            Output[partcount++] = ptr + strlen( delim );
            ptr = ptr + strlen( delim );
            *index += 1;
            if ( *index > 2499 )
                MYABORT( "Erreur dans charToPetscArgcArgv" );
        }
    }
    Output[partcount++] = NULL;
    *index += 1;
}

void DEFSP( ASTER_PETSC_INITIALIZE, aster_petsc_initialize, _IN char *options,
            _IN STRING_SIZE loptions, ASTERINTEGER *ierr ) {
    int ierr2 = 0;
    char **myargs;
    char *options2;
    int myargc;

    options2 = MakeCStrFromFStr( options, loptions );
    myargs = (char **)malloc( 2500 * sizeof( *myargs ) );

    charToPetscArgcArgv( options2, " ", myargs, &myargc );

    *ierr = (ASTERINTEGER)PetscInitialize( &myargc, &myargs, NULL, NULL );
    ierr2 = (ASTERINTEGER)PetscInitializeFortran();
    *ierr += ierr2;
    ierr2 = (ASTERINTEGER)PetscLogDefaultBegin();
    *ierr += ierr2;

    FreeStr( options2 );
    free( myargs );
};

void DEF0( ASTER_PETSC_FINALIZE, aster_petsc_finalize ) { PetscFinalize(); };

#endif

/*
 *   PRIVATE FUNCTIONS
 *
 */

void TraiteMessageErreur( _IN char * );

ASTERINTEGER DEF0( ISJVUP, isjvup ) {
    /* "is jeveux up ?" : retourne 1 si jeveux est démarré/initialisé, sinon 0. */
    return (ASTERINTEGER)get_sh_jeveux_status();
}

void DEFP( XFINI, xfini, _IN ASTERINTEGER *code ) {
    /* XFINI est n'appelé que par JEFINI avec code=19 (=EOFError) */
    /* jeveux est fermé */
    register_sh_jeveux_status( 0 );

    // Do not raise EOFError when using new language description.
    // interruptTry(*code);
}

/*
 *   Ce module crée de nombreux objets Python. Il doit respecter les règles
 *   générales de création des objets et en particulier les règles sur le
 *   compteur de références associé à chaque objet.
 *   Tous les objets sont partagés. Seules des références à des objets peuvent
 *   etre acquises.
 *   Si une fonction a acquis une référence sur un objet elle doit la traiter
 *   proprement, soit en la transférant (habituellement à l'appelant), soit en
 *   la relachant (par appel à Py_DECREF ou Py_XDECREF).
 *   Quand une fonction transfere la propriété d'une référence, l'appelant recoit
 *   une nouvelle référence. Quand la propriété n'est pas transférée, l'appelant
 *   emprunte la référence.
 *   Dans l'autre sens, quand un appelant passe une référence à une fonction, il y a
 *   deux possibilités : la fonction vole une référence à l'objet ou elle ne le fait
 *   pas. Peu de fonctions (de l'API Python) volent des références : les deux exceptions
 *   les plus notables sont PyList_SetItem() et PyTuple_SetItem() qui volent une
 *   référence à l'item qui est inséré dans la liste ou dans le tuple.
 *   Ces fonctions qui volent des références existent, en général, pour alléger
 *   la programmation.
 */
/* ------------------------------------------------------------------ */

void TraiteMessageErreur( _IN char *message ) {
    ASTERINTEGER ier = SIGABRT;
    printf( "%s\n", message );
    if ( PyErr_Occurred() )
        PyErr_Print();
    CALL_ASABRT( &ier );
}

/* ------------------------------------------------------------------ */
void PRE_myabort( _IN const char *nomFichier, _IN const int numeroLigne, _IN const char *message ) {
    /*
    Procedure : PRE_myabort
    Intention
            Cette procedure prepare la chaine de caracteres affichee par TraiteMessageErreur()
            en ajoutant devant cette chaine, le nom du fichier source et le numero
            de la ligne a partir desquels PRE_myabort a ete appelee.
            Puis elle appelle elle-meme TraiteMessageErreur().
            Voir aussi la macro MYABORT qui permet de generer automatiquement le nom
            du fichier et le numero de la ligne.
    */
    char *chaine = (char *)0;
    int longueur = 0;
    longueur += strlen( nomFichier );
    longueur += 1; /* pour le blanc de separation */
    longueur += 5; /* pour le numero de la ligne */
    longueur += 3; /* pour les deux points entre deux blancs */
    longueur += ( message != (const char *)0 ) ? strlen( message ) : 0;
    longueur += 1; /* pour le caractere de fin de chaine */

    chaine = (char *)( malloc( longueur * sizeof( char ) ) );
    sprintf( chaine, "%s %u : %s", nomFichier, numeroLigne, message );
    TraiteMessageErreur( chaine );

    free( chaine );
    chaine = (char *)0;
    longueur = 0;
}

/* ------------------------------------------------------------------ */
#define CALL_GETLTX( a, b, c, d, e, f, g ) CALLSSPPPPP( GETLTX, getltx, a, b, c, d, e, f, g )

void DEFSSPPPPP( GETLTX, getltx, _IN char *motfac, _IN STRING_SIZE lfac, _IN char *motcle,
                 _IN STRING_SIZE lcle, _IN ASTERINTEGER *iocc, _IN ASTERINTEGER *taille,
                 _IN ASTERINTEGER *mxval, _OUT ASTERINTEGER *isval, _OUT ASTERINTEGER *nbval ) {
    /*
    Procedure : getltx_ (appelee par le fortran sous le nom GETLTX)
    */
    PyObject *res = (PyObject *)0;
    PyObject *tup = (PyObject *)0;
    char *mfc = (char *)0;
    char *mcs = (char *)0;
    int ok = 0;
    int nval = 0;
    int ioc = 0;

    mfc = MakeCStrFromFStr( motfac, lfac );
    DEBUG_ASSERT( mfc != (char *)0 );
    mcs = MakeCStrFromFStr( motcle, lcle );
    DEBUG_ASSERT( mcs != (char *)0 );
    ioc = (int)*iocc;
    ioc = ioc - 1;
    res = PyObject_CallMethod( get_sh_etape(), "getltx", "ssiii", mfc, mcs, ioc, (int)*mxval,
                               (int)*taille );

    /*  si le retour est NULL : exception Python a transferer
        normalement a l appelant mais FORTRAN ??? */
    if ( res == NULL )
        MYABORT( "erreur dans la partie Python" );

    ok = PyArg_ParseTuple( res, "iO", &nval, &tup );
    if ( !ok )
        MYABORT( "erreur dans la partie Python" );

    *nbval = (ASTERINTEGER)nval;
    if ( nval < 0 )
        nval = (int)*mxval;
    if ( nval > 0 ) {
        convert( nval, tup, isval );
    }
    Py_DECREF( res ); /*  decrement sur le refcount du retour */
    FreeStr( mfc );
    FreeStr( mcs );
    return;
}

/* ------------------------------------------------------------------ */
void DEFSP( GETFAC, getfac, _IN char *nomfac, _IN STRING_SIZE lfac, _OUT ASTERINTEGER *occu ) {
    /*
      Procedure GETFAC pour le FORTRAN : emule le fonctionnement
      de la procedure equivalente ASTER
      Entrees :
        le nom d un mot cle facteur : nomfac (string)
      Retourne :
        le nombre d occurence de ce mot cle dans les args : occu (entier)
        dans l'etape (ou la commande) courante
    */
    PyObject *res = (PyObject *)0;
    char *mfc;
    mfc = MakeCStrFromFStr( nomfac, lfac );
    res = PyObject_CallMethod( get_sh_etape(), "getfac", "s", mfc );

    /*  si le retour est NULL : exception Python a transferer
        normalement a l appelant mais FORTRAN ??? */
    if ( res == NULL )
        MYABORT( "erreur dans la partie Python" );

    *occu = (ASTERINTEGER)PyLong_AsLong( res );

    Py_DECREF( res ); /*  decrement sur le refcount du retour */
    FreeStr( mfc );
    return;
}

/* ------------------------------------------------------------------ */
void DEFP( GETRAN, getran, _OUT ASTERDOUBLE *rval ) {
    /*
      Procedure GETRAN pour le FORTRAN : recupere un réel aleatoire (loi uniforme 0-1)
      du module python Random
      Entrees :
        neant
      Retourne :
        un reel tiré au hasard
    */
    PyObject *res = (PyObject *)0;
    PyObject *val = (PyObject *)0;
    PyObject *func = GetJdcAttr( "getran" );
    int ok = 0;

    res = PyObject_CallFunction( func, "" );
    /*  si le retour est NULL : exception Python a transferer
        normalement a l appelant mais FORTRAN ??? */
    if ( res == NULL )
        MYABORT( "erreur dans la partie Python" );

    ok = PyArg_ParseTuple( res, "O", &val );
    if ( !ok )
        MYABORT( "erreur dans la partie Python" );

    *rval = (ASTERDOUBLE)PyFloat_AsDouble( val );

    Py_DECREF( func );
    Py_DECREF( res ); /*  decrement sur le refcount du retour */
    return;
}

/* ------------------------------------------------------------------ */
void DEFP( INIRAN, iniran, _IN ASTERINTEGER *jump ) {
    /*
      Procedure INIRAN pour le FORTRAN : recupere un réel aleatoire (loi uniforme 0-1)
      du module python Random
      avec un shift eventuel de jump termes
    */
    PyObject *res = (PyObject *)0;
    PyObject *func = GetJdcAttr( "iniran" );

    res = PyObject_CallFunction( func, "i", (int)*jump );
    Py_DECREF( func );
    Py_DECREF( res ); /*  decrement sur le refcount du retour */
    return;
}

/* ------------------------------------------------------------------ */
void DEFPS( GETMAT, getmat, _INOUT ASTERINTEGER *nbarg, _OUT char *motcle, _IN STRING_SIZE lcle ) {
    /*
      Procedure GETMAT pour le FORTRAN
      Routine a l usage de DEFI_MATERIAU : consultation du catalogue (et non de l etape)
      Retourne :
        le nombre de mots cles facteur sous la commande, y compris en eliminant les blocs
        la liste de leur noms
    */
    PyObject *res = (PyObject *)0;
    PyObject *lnom = (PyObject *)0; /* liste python des noms */
    int nval = 0;
    int k = 0;
    DEBUG_ASSERT( lcle > 0 );
    for ( k = 0; k < lcle; k++ )
        motcle[k] = ' ';
    res = PyObject_CallMethod( get_sh_etape(), "getmat", "" );
    /*  si le retour est NULL : exception Python a transferer
        normalement a l appelant mais FORTRAN ??? */
    if ( res == NULL )
        MYABORT( "erreur dans la partie Python" );
    /*  si non impression du retour */

    if ( !PyArg_ParseTuple( res, "O", &lnom ) )
        MYABORT( "erreur dans la partie Python" );
    nval = PyList_Size( lnom );

    if ( nval > 0 && *nbarg > 0 ) {
        converltx( nval, lnom, motcle, lcle ); /* conversion  */
    }
    *nbarg = (ASTERINTEGER)nval;

    Py_DECREF( res ); /*  decrement sur le refcount du retour */
    return;
}

/* ------------------------------------------------------------------ */
void DEFSPPSSP( GETMJM, getmjm, _IN char *nomfac, _IN STRING_SIZE lfac, _IN ASTERINTEGER *iocc,
                _IN ASTERINTEGER *nbval, _OUT char *motcle, _IN STRING_SIZE lcle, _OUT char *type,
                _IN STRING_SIZE ltyp, _OUT ASTERINTEGER *nbarg ) {
    /*
      Procedure GETMJM : emule la procedure equivalente ASTER
       Retourne les nbval premiers mots cles du mot cle facteur nomfac du catalogue
       de la commande en cours
      Entrees :
       nomfac : nom du mot cle facteur
       iocc   : numero d occurence du mot cle facteur
       nbval  : nombre de mots cles facteurs demandes
      Retourne :
       motcle : liste des mots cles du mot cle facteur demande
       type   : liste des types des mots cles du mot cle facteur demande
                R8 , R8L : un reel ou une liste de reels ;
                C8 , C8L : un complexe ou une liste de complexes ;
                 ...
                CO , COL : un concept ou une liste de concepts.
       nbarg  : nombre d arguments des mots cles du mot cle facteur
    */

    PyObject *res = (PyObject *)0;
    PyObject *lnom = (PyObject *)0;
    PyObject *lty = (PyObject *)0; /* liste python des noms */
    int nval = 0;
    int k = 0;
    int ioc = 0;
    char *mfc;
    DEBUG_ASSERT( ltyp > 0 );
    for ( k = 0; k < ltyp; k++ )
        type[k] = ' ';
    ioc = (int)*iocc;
    ioc = ioc - 1;
    mfc = MakeCStrFromFStr( nomfac, lfac );
    res = PyObject_CallMethod( get_sh_etape(), "getmjm", "sii", mfc, ioc, (int)*nbval );
    /*  si le retour est NULL : exception Python a transferer
        normalement a l appelant mais FORTRAN ??? */
    if ( res == NULL )
        MYABORT( "erreur dans la partie Python" );
    /*  si non impression du retour */

    if ( !PyArg_ParseTuple( res, "OO", &lnom, &lty ) )
        MYABORT( "erreur dans la partie Python" );
    nval = (int)PyList_Size( lnom );
    *nbarg = nval > *nbval ? (ASTERINTEGER)-nval : (ASTERINTEGER)nval;
    DEBUG_ASSERT( ( ( nval <= *nbval ) && ( *nbarg == nval ) ) || ( *nbarg == -nval ) );
    if ( *nbarg < 0 )
        nval = (int)*nbval;

    if ( nval > 0 ) {
        converltx( nval, lnom, motcle, lcle ); /* conversion  */
        converltx( nval, lty, type, ltyp );
    }

    /*
    A la demande des developpeurs (J. Pellet), le nom des concepts retourne par
    la methode EXECUTION.getmjm (par exemple grma) est ici remplace par
    la chaine CO (pour COncept).
    les types retournes sont donc parmi les valeurs : R8 , C8 , IS , TX et CO.
    */
    for ( k = 0; k < nval * ltyp; k += ltyp ) {
        char *mot = (char *)0;
        mot = type + k;
        if ( strncmp( mot, "R8", 2 ) != 0 && strncmp( mot, "IS", 2 ) != 0 &&
             strncmp( mot, "TX", 2 ) != 0 && strncmp( mot, "C8", 2 ) != 0 ) {
            int j = 0;

            DEBUG_ASSERT( ltyp > 2 );
            mot[0] = 'C';
            mot[1] = 'O';
            for ( j = 2; j < ltyp; j++ )
                mot[j] = ' ';
        }
    }
    Py_DECREF( res ); /*  decrement sur le refcount du retour */
    FreeStr( mfc );
    return;
}

/* ------------------------------------------------------------------ */
ASTERINTEGER DEFSS( GETEXM, getexm, _IN char *motfac, _IN STRING_SIZE lfac, _IN char *motcle,
                    _IN STRING_SIZE lcle ) {
    /*
      Procedure GETEXM pour le FORTRAN : emule le fonctionnement
      de la procedure equivalente ASTER
      Entrees :
        le nom d un mot cle facteur : motfac (string)
        le nom d un mot cle simple ou sous mot cle : motcle (string)
      Retourne :
        0 si n existe pas 1 si existe
    */
    PyObject *res = (PyObject *)0;
    char *mfc, *mcs;
    ASTERINTEGER presence;
    if ( get_sh_etape() == Py_None ) {
        return (ASTERINTEGER)0;
    }
    DEBUG_ASSERT( motcle != (char *)0 );
    mfc = MakeCStrFromFStr( motfac, lfac );
    mcs = MakeCStrFromFStr( motcle, lcle );
    res = PyObject_CallMethod( get_sh_etape(), "getexm", "ss", mfc, mcs );
    /*  si le retour est NULL : exception Python a transferer
        normalement a l appelant mais FORTRAN ??? */
    if ( res == NULL )
        MYABORT( "erreur dans la partie Python" );
    presence = (ASTERINTEGER)PyLong_AsLong( res );
    /*  decrement sur le refcount du retour */
    Py_DECREF( res ); /*  decrement sur le refcount du retour */
    FreeStr( mfc );
    FreeStr( mcs );
    return presence;
}

/* ------------------------------------------------------------------ */
void DEFSSS( GETRES, getres, _OUT char *nomres, _IN STRING_SIZE lres, _OUT char *concep,
             _IN STRING_SIZE lconc, _OUT char *nomcmd, _IN STRING_SIZE lcmd ) {
    /*
      Procedure GETRES pour le FORTRAN : emule le fonctionnement
      de la procedure equivalente ASTER
      Retourne
        le nom utilisateur du resultat : nomres (string)
        le nom du concept resultat     : concep (string)
        le nom de la commande          : nomcmd (string)
    */
    PyObject *res = (PyObject *)0;
    PyObject *etape;
    int ok;
    char *ss1, *ss2, *ss3;

    /* (MC) le 1er test ne me semble pas suffisant car entre deux commandes,
       commande n'est pas remis à (PyObject*)0... */
    etape = get_sh_etape();
    if ( etape == (PyObject *)0 || PyObject_HasAttrString( etape, "getres" ) == 0 ) {
        /* Aucune commande n'est active on retourne des chaines blanches */
        BlankStr( nomres, lres );
        BlankStr( concep, lconc );
        BlankStr( nomcmd, lcmd );
        return;
    }
    res = PyObject_CallMethod( etape, "getres", "" );

    ok = PyArg_ParseTuple( res, "sss", &ss1, &ss2, &ss3 );
    if ( !ok )
        MYABORT( "erreur dans la partie Python" );

    /* le fortran attend des chaines de caracteres completees par des blancs */
    CopyCStrToFStr( nomres, ss1, lres );
    CopyCStrToFStr( concep, ss2, lconc );
    CopyCStrToFStr( nomcmd, ss3, lcmd );

    Py_DECREF( res ); /*  decrement sur le refcount du retour */
    return;
}

/* ------------------------------------------------------------------ */
void DEFSSPPPP( GETVC8_WRAP, getvc8_wrap, _IN char *motfac, _IN STRING_SIZE lfac, _IN char *motcle,
                _IN STRING_SIZE lcle, _IN ASTERINTEGER *iocc, _IN ASTERINTEGER *mxval,
                _INOUT ASTERDOUBLE *val, _OUT ASTERINTEGER *nbval ) {
    /*
      Procedure GETVC8 pour le FORTRAN : emule le fonctionnement
      de la procedure equivalente ASTER
      Entrees :
        le nom d un mot cle facteur : motfac (string)
        le nom d un mot cle simple ou sous mot cle : motcle (string)
        le numero de l occurence du mot cle facteur : iocc (entier)
        le nombre max de valeur attendues dans val : mxval (entier)
      Retourne :
        le tableau des valeurs attendues : val (2 reels (double) par complexe)
        le nombre de valeurs effectivement retournees : nbval (entier)
           si pas de valeur nbval =0
           si plus de valeur que mxval nbval <0 et valeur abs = nbre valeurs
           si moins de valeurs que mxval nbval>0 et egal au nombre retourne
    */
    PyObject *res = (PyObject *)0;
    PyObject *tup = (PyObject *)0;
    int ok = 0;
    int nval = 0;
    int ioc = 0;
    char *mfc = (char *)0;
    char *mcs = (char *)0;
    mfc = MakeCStrFromFStr( motfac, lfac );
    DEBUG_ASSERT( mfc != (char *)0 );
    mcs = MakeCStrFromFStr( motcle, lcle );
    /*
            VERIFICATION
            Si le mot-cle simple est recherche sous un mot-cle facteur et uniquement dans ce
            cas, le numero d'occurrence (*iocc) doit etre strictement positif.
            Si le mot-cle simple est recherche sans un mot-cle facteur iocc n'est pas utilise

    */
    if ( isalpha( mfc[0] ) && ( *iocc <= 0 ) ) {
        printf( "<F> GETVC8 : le numero d'occurence (IOCC=%ld) est invalide\n", *iocc );
        printf( "             commande : %s\n",
                PyUnicode_AsUTF8( PyObject_CallMethod( get_sh_etape(), "getName", "" ) ) );
        printf( "             mot-cle facteur : %s\n", mfc );
        printf( "             mot-cle simple  : %s\n", mcs );
        MYABORT( "erreur d'utilisation detectee" );
    }

    ioc = (int)*iocc;
    ioc = ioc - 1;
    res = PyObject_CallMethod( get_sh_etape(), "getvc8", "ssii", mfc, mcs, ioc, (int)*mxval );

    /*  si le retour est NULL : exception Python a transferer
        normalement a l appelant mais FORTRAN ??? */
    if ( res == NULL )
        MYABORT( "erreur dans la partie Python" );
    DEBUG_ASSERT( PyTuple_Check( res ) );
    ok = PyArg_ParseTuple( res, "iO", &nval, &tup );
    if ( !ok )
        MYABORT( "erreur dans la partie Python" );

    *nbval = (ASTERINTEGER)nval;
    if ( nval < 0 )
        nval = (int)*mxval;
    if ( nval > 0 ) {
        convc8( nval, tup, val );
    }
    Py_DECREF( res );
    FreeStr( mfc );
    FreeStr( mcs );
    return;
}

/* ------------------------------------------------------------------ */
void DEFSSPPPP( GETVR8_WRAP, getvr8_wrap, _IN char *motfac, _IN STRING_SIZE lfac, _IN char *motcle,
                _IN STRING_SIZE lcle, _IN ASTERINTEGER *iocc, _IN ASTERINTEGER *mxval,
                _INOUT ASTERDOUBLE *val, _OUT ASTERINTEGER *nbval ) {
    /*
      Procedure GETVR8 pour le FORTRAN : emule le fonctionnement
      de la procedure equivalente ASTER
      Entrees :
        le nom d un mot cle facteur : motfac (string)
        le nom d un mot cle simple ou sous mot cle : motcle (string)
        le numero de l occurence du mot cle facteur : iocc (entier)
        le nombre max de valeur attendues dans val : mxval (entier)
      Retourne :
        le tableau des valeurs attendues : val (tableau de R8    )
        le nombre de valeurs effectivement retournees : nbval (entier)
           si pas de valeur nbval =0
           si plus de valeur que mxval nbval <0 et valeur abs = nbre valeurs
           si moins de valeurs que mxval nbval>0 et egal au nombre retourne
    */
    PyObject *res = (PyObject *)0;
    PyObject *tup = (PyObject *)0;
    int ok = 0;
    int nval = 0;
    int ioc = 0;
    char *mfc = (char *)0;
    char *mcs = (char *)0;
    mfc = MakeCStrFromFStr( motfac, lfac ); /* conversion chaine fortran en chaine C */
    DEBUG_ASSERT( mfc != (char *)0 );
    mcs = MakeCStrFromFStr( motcle, lcle );
    /*
     * VERIFICATION
     * Si le mot-cle simple est recherche sous un mot-cle facteur et uniquement dans ce cas,
     * le numero d'occurrence (*iocc) doit etre strictement positif.
     * Si le mot-cle simple est recherche sans un mot-cle facteur iocc n'est pas utilise
     */
    if ( isalpha( mfc[0] ) && ( *iocc <= 0 ) ) {
        printf( "<F> GETVR8 : le numero d'occurence (IOCC=%ld) est invalide\n", *iocc );
        printf( "             commande : %s\n",
                PyUnicode_AsUTF8( PyObject_CallMethod( get_sh_etape(), "getName", "" ) ) );
        printf( "             mot-cle facteur : %s\n", mfc );
        printf( "             mot-cle simple  : %s\n", mcs );
        MYABORT( "erreur d'utilisation detectee" );
    }
    ioc = (int)*iocc;
    ioc = ioc - 1;
    res = PyObject_CallMethod( get_sh_etape(), "getvr8", "ssii", mfc, mcs, ioc, (int)*mxval );
    /*  si le retour est NULL : exception Python a transferer
        normalement a l appelant mais FORTRAN ??? */
    if ( res == NULL )
        MYABORT( "erreur dans la partie Python" );
    DEBUG_ASSERT( PyTuple_Check( res ) );
    ok = PyArg_ParseTuple( res, "iO", &nval, &tup );
    if ( !ok )
        MYABORT( "erreur dans la partie Python" );

    *nbval = (ASTERINTEGER)nval;
    if ( nval < 0 )
        nval = (int)*mxval;
    if ( nval > 0 ) {
        convr8( nval, tup, val );
    }
    Py_DECREF( res ); /*  decrement sur le refcount du retour */
    FreeStr( mfc );
    FreeStr( mcs );
    return;
}

/* ------------------------------------------------------------------ */
void DEFSSPPPP( GETVIS_WRAP, getvis_wrap, _IN char *motfac, _IN STRING_SIZE lfac, _IN char *motcle,
                _IN STRING_SIZE lcle, _IN ASTERINTEGER *iocc, _IN ASTERINTEGER *mxval,
                _INOUT ASTERINTEGER *val, _OUT ASTERINTEGER *nbval ) {
    /*
      Procedure GETVIS pour le FORTRAN : emule le fonctionnement
      de la procedure equivalente ASTER
      Entrees :
        le nom d un mot cle facteur : motfac (string)
        le nom d un mot cle simple ou sous mot cle : motcle (string)
        le numero de l occurence du mot cle facteur : iocc (entier)
        le nombre max de valeur attendues dans val : mxval (entier)
      Retourne :
        le tableau des valeurs attendues : val (tableau d entier )
        le nombre de valeurs effectivement retournees : nbval (entier)
           si pas de valeur nbval =0
           si plus de valeur que mxval nbval <0 et valeur abs = nbre valeurs
           si moins de valeurs que mxval nbval>0 et egal au nombre retourne
    */
    PyObject *res = (PyObject *)0;
    PyObject *tup = (PyObject *)0;
    int ok = 0;
    int nval = 0;
    int ioc = 0;
    char *mfc = (char *)0;
    char *mcs = (char *)0;
    DEBUG_ASSERT( ( *iocc > 0 ) || ( FStrlen( motfac, lfac ) == 0 ) );
    mfc = MakeCStrFromFStr( motfac, lfac ); /* conversion chaine fortran en chaine C */
    DEBUG_ASSERT( mfc != (char *)0 );
    mcs = MakeCStrFromFStr( motcle, lcle );
    /*
            VERIFICATION
            Si le mot-cle simple est recherche sous un mot-cle facteur et uniquement dans ce
            cas, le numero d'occurrence (*iocc) doit etre strictement positif.
            Si le mot-cle simple est recherche sans un mot-cle facteur iocc n'est pas utilise

    */
    if ( isalpha( mfc[0] ) && ( *iocc <= 0 ) ) {
        printf( "<F> GETVIS : le numero d'occurence (IOCC=%ld) est invalide\n", *iocc );
        printf( "             commande : %s\n",
                PyUnicode_AsUTF8( PyObject_CallMethod( get_sh_etape(), "getName", "" ) ) );
        printf( "             mot-cle facteur : %s\n", mfc );
        printf( "             mot-cle simple  : %s\n", mcs );
        MYABORT( "erreur d'utilisation detectee" );
    }
    ioc = (int)*iocc;
    ioc = ioc - 1;
    res = PyObject_CallMethod( get_sh_etape(), "getvis", "ssii", mfc, mcs, ioc, (int)*mxval );

    /*  si le retour est NULL : exception Python a transferer
        normalement a l appelant mais FORTRAN ??? */
    if ( res == NULL )
        MYABORT( "erreur dans la partie Python" );

    ok = PyArg_ParseTuple( res, "iO", &nval, &tup );
    if ( !ok )
        MYABORT( "erreur dans la partie Python" );

    *nbval = (ASTERINTEGER)nval;
    if ( nval < 0 )
        nval = (int)*mxval;
    if ( nval > 0 ) {
        convert( nval, tup, val );
    }
    Py_DECREF( res ); /*  decrement sur le refcount du retour */
    FreeStr( mfc );
    FreeStr( mcs );
    return;
}

/* ------------------------------------------------------------------ */
void DEFSSPPSP( GETVTX_WRAP, getvtx_wrap, _IN char *motfac, _IN STRING_SIZE lfac, _IN char *motcle,
                _IN STRING_SIZE lcle, _IN ASTERINTEGER *iocc, _IN ASTERINTEGER *mxval,
                _INOUT char *txval, _IN STRING_SIZE ltx, _OUT ASTERINTEGER *nbval ) {
    /*
      Procedure GETVTX pour le FORTRAN : emule le fonctionnement
      de la procedure equivalente ASTER
      Entrees :
        le nom d un mot cle facteur : motfac (string)
        le nom d un mot cle simple ou sous mot cle : motcle (string)
        le numero de l occurence du mot cle facteur : iocc (entier)
        le nombre max de valeur attendues dans val : mxval (entier)
      Retourne :
        le tableau des valeurs attendues : txval (tableau de string)
        ATTENTION : txval arrive avec une valeur par defaut
        le nombre de valeurs effectivement retournees : nbval (entier)
           si pas de valeur nbval =0
           si plus de valeur que mxval nbval <0 et valeur abs = nbre valeurs
           si moins de valeurs que mxval nbval>0 et egal au nombre retourne

    */
    PyObject *res = (PyObject *)0;
    PyObject *tup = (PyObject *)0;
    int ok = 0;
    int nval = 0;
    int ioc = 0;
    char *mfc = (char *)0;
    char *mcs = (char *)0;

    mfc = MakeCStrFromFStr( motfac, lfac );
    mcs = MakeCStrFromFStr( motcle, lcle );
    /*
            VERIFICATION
            Si le mot-cle simple est recherche sous un mot-cle facteur et uniquement dans ce
            cas, le numero d'occurrence (*iocc) doit etre strictement positif.
            Si le mot-cle simple est recherche sans un mot-cle facteur iocc n'est pas utilise

    */
    if ( isalpha( mfc[0] ) && ( *iocc <= 0 ) ) {
        printf( "<F> GETVTX : le numero d'occurence (IOCC=%ld) est invalide\n", *iocc );
        printf( "             commande : %s\n",
                PyUnicode_AsUTF8( PyObject_CallMethod( get_sh_etape(), "getName", "" ) ) );
        printf( "             mot-cle facteur : %s\n", mfc );
        printf( "             mot-cle simple  : %s\n", mcs );
        MYABORT( "erreur d'utilisation detectee" );
    }
    ioc = (int)*iocc;
    ioc = ioc - 1;
    res = PyObject_CallMethod( get_sh_etape(), "getvtx", "ssii", mfc, mcs, ioc, (int)*mxval );

    /*  si le retour est NULL : exception Python a transferer
        normalement a l appelant mais FORTRAN ??? */
    if ( res == NULL ) {
        printf( "<F> GETVTX : numero d'occurence (IOCC=%ld) \n", *iocc );
        printf( "             commande : %s\n",
                PyUnicode_AsUTF8( PyObject_CallMethod( get_sh_etape(), "getName", "" ) ) );
        printf( "             mot-cle facteur : %s\n", mfc );
        printf( "             mot-cle simple  : %s\n", mcs );
        MYABORT( "erreur dans la partie Python" );
    }

    ok = PyArg_ParseTuple( res, "iO", &nval, &tup );
    if ( !ok )
        MYABORT( "erreur au decodage d'une chaine dans le module C aster.getvtx" );

    *nbval = (ASTERINTEGER)nval;
    if ( nval < 0 )
        nval = (int)*mxval;
    if ( nval > 0 ) {
        convertxt( nval, tup, txval, ltx );
    }
    /* ATTENTION : il ne faut decrementer le compteur de references de res
     *             qu'apres en avoir fini avec l'utilisation de tup.
     *             NE PAS decrementer le compteur de references de tup car
     *             la disparition de res entrainera un decrement automatique
     *             du compteur de tup (res=(nbval,tup))
     */
    Py_DECREF( res ); /*  decrement sur le refcount du retour */
    FreeStr( mfc );
    FreeStr( mcs );
    return;
}

/* ------------------------------------------------------------------ */
void DEFSSPPSP( GETVID_WRAP, getvid_wrap, _IN char *motfac, _IN STRING_SIZE lfac, _IN char *motcle,
                _IN STRING_SIZE lcle, _IN ASTERINTEGER *iocc, _IN ASTERINTEGER *mxval,
                _INOUT char *txval, _IN STRING_SIZE ltx, _OUT ASTERINTEGER *nbval ) {
    /*
      Procedure GETVID pour le FORTRAN : emule le fonctionnement
      de la procedure equivalente ASTER
      Entrees :
        le nom d un mot cle facteur : motfac (string)
        le nom d un mot cle simple ou sous mot cle : motcle (string)
        le numero de l occurence du mot cle facteur : iocc (entier)
        le nombre max de valeur attendues dans val : mxval (entier)
      Retourne :
        le tableau des valeurs attendues : val (tableau de string)
        le nombre de valeurs effectivement retournees : nbval (entier)
           si pas de valeur nbval =0
           si plus de valeur que mxval nbval <0 et valeur abs = nbre valeurs
           si moins de valeurs que mxval nbval>0 et egal au nombre retourne
    */
    PyObject *res = (PyObject *)0;
    PyObject *tup = (PyObject *)0;
    int ok, nval, ioc;
    char *mfc;
    char *mcs;
    DEBUG_ASSERT( ( *iocc > 0 ) || ( FStrlen( motfac, lfac ) == 0 ) );
    mfc = MakeCStrFromFStr( motfac, lfac ); /* conversion chaine fortran en chaine C */
    DEBUG_ASSERT( mfc != (char *)0 );
    mcs = MakeCStrFromFStr( motcle, lcle );
    /*
            VERIFICATION
            Si le mot-cle simple est recherche sous un mot-cle facteur et uniquement dans ce
            cas, le numero d'occurrence (*iocc) doit etre strictement positif.
            Si le mot-cle simple est recherche sans un mot-cle facteur iocc n'est pas utilise

    */
    if ( isalpha( mfc[0] ) && ( *iocc <= 0 ) ) {
        printf( "<F> GETVID : le numero d'occurence (IOCC=%ld) est invalide\n", *iocc );
        printf( "             commande : %s\n",
                PyUnicode_AsUTF8( PyObject_CallMethod( get_sh_etape(), "getName", "" ) ) );
        printf( "             mot-cle facteur : %s\n", mfc );
        printf( "             mot-cle simple  : %s\n", mcs );
        MYABORT( "erreur d'utilisation detectee" );
    }
    ioc = (int)*iocc;
    ioc = ioc - 1;
    res = PyObject_CallMethod( get_sh_etape(), "getvid", "ssii", mfc, mcs, ioc, (int)*mxval );

    /*  si le retour est NULL : exception Python a transferer
        normalement a l appelant mais FORTRAN ??? */
    if ( res == NULL )
        MYABORT( "erreur dans la partie Python" );

    ok = PyArg_ParseTuple( res, "iO", &nval, &tup );
    if ( !ok )
        MYABORT( "erreur dans la partie Python" );

    *nbval = (ASTERINTEGER)nval;
    if ( nval < 0 )
        nval = (int)*mxval;
    if ( nval > 0 ) {
        convertxt( nval, tup, txval, ltx );
    }
    Py_DECREF( res ); /*  decrement sur le refcount du retour */
    FreeStr( mfc );
    FreeStr( mcs );
    return;
}

/* ------------------------------------------------------------------ */
void getvpy( _IN const char *mfc, _IN const char *mcs, _IN const int ioc, _OUT PyObject **rhs ) {
    /*
      Return an array of Python objects for a user keyword.

      Inputs:
          motfac (string): factor keyword or " "
          motcle (string): simple keyword
          iocc (int): number of occurrence (1-based)

      Outputs:
          rhs (PyObject**): pointer on the tuple

      The caller must free the array of 'PyObject*'!
    */
    PyObject *res = (PyObject *)0;
    PyObject *tup = (PyObject *)0;
    int nbval = 0;
    int ok = 0;
    /*
            VERIFICATION
            Si le mot-cle simple est recherche sous un mot-cle facteur et uniquement dans ce
            cas, le numero d'occurrence (ioc) doit etre strictement positif.
            Si le mot-cle simple est recherche sans un mot-cle facteur ioc n'est pas utilise
    */
    if ( isalpha( mfc[0] ) && ( ioc <= 0 ) ) {
        printf( "<F> GETVPY : le numero d'occurence (IOCC=%d) est invalide\n", ioc );
        printf( "             commande : %s\n",
                PyUnicode_AsUTF8( PyObject_CallMethod( get_sh_etape(), "getName", "" ) ) );
        printf( "             mot-cle facteur : %s\n", mfc );
        printf( "             mot-cle simple  : %s\n", mcs );
        MYABORT( "erreur d'utilisation detectee" );
    }
    res = PyObject_CallMethod( get_sh_etape(), "getvpy", "ssi", mfc, mcs, ioc - 1 );

    if ( res == NULL )
        MYABORT( "erreur dans la partie Python" );

    ok = PyArg_ParseTuple( res, "iO", &nbval, &tup );
    if ( !ok )
        MYABORT( "erreur dans la partie Python" );

    *rhs = tup;
    Py_INCREF( tup );
    Py_DECREF( res ); /*  decrement sur le refcount du retour */
    return;
}

/* ------------------------------------------------------------------ */
void DEFP( PUTVIR, putvir, _IN ASTERINTEGER *ival ) {
    /*
       Entrees:
          ival entier à affecter
       Fonction:
          renseigner l'attribut valeur associé à la sd
          n'est utile que pour DEFI_FICHIER, EXTR_TABLE
          cet attribut est ensuite évalué par la méthode traite_value
          de B_ETAPE.py
    */
    PyObject *res = (PyObject *)0;

    res = PyObject_CallMethod( get_sh_etape(), "setres", "i", (int)*ival );
    if ( res == NULL )
        MYABORT( "erreur a l appel de setres dans la partie Python" );

    Py_DECREF( res );
}

void DEFP( PUTVRR, putvrr, _IN ASTERDOUBLE *rval ) {
    /*
       Entrees:
          rval réel à affecter
       Fonction:
          renseigner l'attribut valeur associé à la sd
          n'est utile que pour EXTR_TABLE
          cet attribut est ensuite évalué par la méthode traite_value
          de B_ETAPE.py
    */
    PyObject *res = (PyObject *)0;

    res = PyObject_CallMethod( get_sh_etape(), "setres", "d", (double)*rval );
    if ( res == NULL )
        MYABORT( "erreur a l appel de setres dans la partie Python" );

    Py_DECREF( res );
}

/* ------------------------------------------------------------------ */
void DEFP( GCECDU, gcecdu, ASTERINTEGER *numint ) {
    /*
      Sortie :
        numint  numero de l operateur de la commande
      Fonction:
         Recuperation du numero de l operateur
    */
    PyObject *res = (PyObject *)0;
    res = PyObject_CallMethod( get_sh_etape(), "getoper", "" );
    /*
                Si le retour est NULL : une exception a ete levee dans le code Python appele
                Cette exception est a transferer normalement a l appelant mais FORTRAN ???
                On produit donc un abort en ecrivant des messages sur la stdout
    */
    if ( res == NULL )
        MYABORT( "erreur a l appel de gcecdu dans la partie Python" );

    *numint = (ASTERINTEGER)PyLong_AsLong( res );
    Py_DECREF( res );
}

/* ------------------------------------------------------------------ */
void gcncon2_( char *type, char *resul, STRING_SIZE ltype, int lresul ) {
    /* CCAR : cette fonction devrait s appeler gcncon mais elle est utilisee par
              tous les operateurs (???) et pas seulement dans les macros
              Pour le moment il a ete decide de ne pas l'emuler dans le superviseur
              Python mais d'utiliser les fonctions FORTRAN existantes
              Ceci a l avantage d'assurer la coherence entre tous les operateurs
              et de conserver les fonctionnalites de poursuite pour les macros
    */
    /*
      Entrees:
        type vaut soit
                '.' : le concept sera detruit en fin de job
                '_' : le concept ne sera pas detruit

      Sorties:
        resul  nom d'un concept delivre par le superviseur
               Ce nom est de la forme type // '000ijkl' ou ijkl est un nombre
               incremente a chaque appel pour garantir l unicite des noms

      Fonction:
        Delivrer un nom de concept non encore utilise et unique
    */
    MYABORT( "Cette procedure n est pas implementee" );
}

/* ------------------------------------------------------------------ */
static char getvectjev_doc[] = "\
Retourne la valeur du concept nomsd \n\
dans un tuple.";

static PyObject *aster_getvectjev( PyObject *self, PyObject *args ) {
    char *nomsd, *nomsd32;
    char *nomob;
    ASTERDOUBLE *f;
    ASTERINTEGER *l;
    ASTERINTEGER4 *i4;
    char *kvar;
    PyObject *tup = NULL;
    ASTERINTEGER lcon, iob;
    int ishf = 0, ilng = 0;
    ASTERINTEGER shf;
    ASTERINTEGER lng;
    ASTERINTEGER ctype = 0;
    int i, ksize = 0;
    char *iaddr;

    if ( !PyArg_ParseTuple( args, "s|ii:getvectjev", &nomsd, &ishf, &ilng ) )
        return NULL;
    shf = (ASTERINTEGER)ishf;
    lng = (ASTERINTEGER)ilng;
    iob = 0;
    nomsd32 = MakeFStrFromCStr( nomsd, 32 );
    nomob = MakeBlankFStr( 24 );

    try {
        CALL_JEMARQ();
        CALL_GETCON( nomsd32, &iob, &shf, &lng, &ctype, &lcon, &iaddr, nomob );
        FreeStr( nomsd32 );
        FreeStr( nomob );
        if ( ctype < 0 ) {
            /* Erreur : vecteur jeveux inexistant, on retourne None */
            CALL_JEDEMA();
            endTry();
            Py_INCREF( Py_None );
            return Py_None;
        } else if ( ctype == 0 ) {
            /* Liste vide */
            tup = PyTuple_New( 0 );
        } else if ( ctype == 1 ) {
            /* REEL */
            f = (ASTERDOUBLE *)iaddr;
            tup = PyTuple_New( (Py_ssize_t)lcon );
            for ( i = 0; i < lcon; i++ ) {
                PyTuple_SetItem( tup, i, PyFloat_FromDouble( (double)f[i] ) );
            }
        } else if ( ctype == 2 ) {
            /* ENTIER */
            l = (ASTERINTEGER *)iaddr;
            tup = PyTuple_New( (Py_ssize_t)lcon );
            for ( i = 0; i < lcon; i++ ) {
                PyTuple_SetItem( tup, i, PyLong_FromLong( (long)l[i] ) );
            }
        } else if ( ctype == 9 ) {
            /* ENTIER COURT */
            i4 = (ASTERINTEGER4 *)iaddr;
            tup = PyTuple_New( (Py_ssize_t)lcon );
            for ( i = 0; i < lcon; i++ ) {
                PyTuple_SetItem( tup, i, PyLong_FromLong( (long)i4[i] ) );
            }
        } else if ( ctype == 3 ) {
            /* COMPLEXE */
            f = (ASTERDOUBLE *)iaddr;
            tup = PyTuple_New( (Py_ssize_t)lcon );
            for ( i = 0; i < lcon; i++ ) {
                PyTuple_SetItem( tup, i,
                                 PyComplex_FromDoubles( (double)f[2 * i], (double)f[2 * i + 1] ) );
            }
        } else if ( ctype == 4 || ctype == 5 || ctype == 6 || ctype == 7 || ctype == 8 ) {
            switch ( ctype ) {
            case 4:
                ksize = 8;
                break;
            case 5:
                ksize = 16;
                break;
            case 6:
                ksize = 24;
                break;
            case 7:
                ksize = 32;
                break;
            case 8:
                ksize = 80;
                break;
            }
            /* CHAINE DE CARACTERES */
            tup = PyTuple_New( (Py_ssize_t)lcon );
            for ( i = 0; i < lcon; i++ ) {
                kvar = iaddr + i * ksize;
                PyTuple_SetItem( tup, i, PyUnicode_FromStringAndSize( kvar, ksize ) );
            }
        }
        CALL_JEDETR( "&&GETCON.PTEUR_NOM" );
        CALL_JEDEMA();
    }
    exceptAll {
        CALL_JEDEMA();
        FreeStr( nomsd32 );
        FreeStr( nomob );
        raiseException();
    }
    endTry();
    return tup;
}

static char getcolljev_doc[] = "\
\n\
Retourne la valeur du concept nomsd \n\
dans un tuple.";

/* ------------------------------------------------------------------ */
static PyObject *aster_getcolljev( PyObject *self, PyObject *args ) {
    char *nomsd, *nom, *nomsd32;
    char *nomob;
    ASTERDOUBLE *f;
    ASTERINTEGER *l;
    ASTERINTEGER4 *i4;
    char *kvar;
    PyObject *tup = NULL, *dico, *key;
    ASTERINTEGER iob, j, ishf, ilng;
    ASTERINTEGER lcon;
    ASTERINTEGER ctype = 0;
    ASTERINTEGER *val, nbval;
    int i, ksize = 0;
    char *iaddr;

    if ( !PyArg_ParseTuple( args, "s:getcolljev", &nomsd ) )
        return NULL;

    /* Taille de la collection */
    nbval = 1;
    nomsd32 = MakeFStrFromCStr( nomsd, 32 );
    nomob = MakeBlankFStr( 24 );
    val = (ASTERINTEGER *)malloc( ( nbval ) * sizeof( ASTERINTEGER ) );
    nom = MakeFStrFromCStr( "LIST_COLLECTION", 24 );
    CALL_JEMARQ();
    CALL_TAILSD( nom, nomsd32, val, &nbval );
    iob = val[0];
#define DictSetAndDecRef( dico, key, item )                                                        \
    PyDict_SetItem( dico, key, item );                                                             \
    Py_DECREF( key );                                                                              \
    Py_DECREF( item );
    dico = PyDict_New();
    try {
        for ( j = 1; j < iob + 1; j++ ) {
            ishf = 0;
            ilng = 0;
            CALL_GETCON( nomsd32, &j, &ishf, &ilng, &ctype, &lcon, &iaddr, nomob );
            if ( nomob[0] == ' ' ) {
                key = PyLong_FromLong( (long)j );
            } else {
                key = PyUnicode_FromStringAndSize( nomob, 24 );
            }
            switch ( ctype ) {
            case 0:
                Py_INCREF( Py_None );
                PyDict_SetItem( dico, key, Py_None );
                Py_DECREF( key );
                break;
            case 1:
                /* REEL */
                f = (ASTERDOUBLE *)iaddr;
                tup = PyTuple_New( (Py_ssize_t)lcon );
                for ( i = 0; i < lcon; i++ ) {
                    PyTuple_SetItem( tup, i, PyFloat_FromDouble( (double)f[i] ) );
                }
                DictSetAndDecRef( dico, key, tup );
                break;
            case 2:
                /* ENTIER */
                l = (ASTERINTEGER *)iaddr;
                tup = PyTuple_New( (Py_ssize_t)lcon );
                for ( i = 0; i < lcon; i++ ) {
                    PyTuple_SetItem( tup, i, PyLong_FromLong( (long)l[i] ) );
                }
                DictSetAndDecRef( dico, key, tup );
                break;
            case 9:
                /* ENTIER COURT */
                i4 = (ASTERINTEGER4 *)iaddr;
                tup = PyTuple_New( (Py_ssize_t)lcon );
                for ( i = 0; i < lcon; i++ ) {
                    PyTuple_SetItem( tup, i, PyLong_FromLong( (long)i4[i] ) );
                }
                DictSetAndDecRef( dico, key, tup );
                break;
            case 3:
                /* COMPLEXE */
                f = (ASTERDOUBLE *)iaddr;
                tup = PyTuple_New( (Py_ssize_t)lcon );
                for ( i = 0; i < lcon; i++ ) {
                    PyTuple_SetItem(
                        tup, i, PyComplex_FromDoubles( (double)f[2 * i], (double)f[2 * i + 1] ) );
                }
                DictSetAndDecRef( dico, key, tup );
                break;
            case 4:
            case 5:
            case 6:
            case 7:
            case 8:
                switch ( ctype ) {
                case 4:
                    ksize = 8;
                    break;
                case 5:
                    ksize = 16;
                    break;
                case 6:
                    ksize = 24;
                    break;
                case 7:
                    ksize = 32;
                    break;
                case 8:
                    ksize = 80;
                    break;
                }
                /* CHAINE DE CARACTERES */
                tup = PyTuple_New( (Py_ssize_t)lcon );
                for ( i = 0; i < lcon; i++ ) {
                    kvar = iaddr + i * ksize;
                    PyTuple_SetItem( tup, i, PyUnicode_FromStringAndSize( kvar, ksize ) );
                }
                DictSetAndDecRef( dico, key, tup );
                break;
            default:
                /* Erreur */
                FreeStr( nom );
                FreeStr( nomob );
                FreeStr( nomsd32 );
                free( val );
                raiseExceptionString( PyExc_KeyError, "Concept inexistant, type inconnu" );
            }
        }
        CALL_JEDETR( "&&GETCON.PTEUR_NOM" );
        CALL_JEDEMA();
        FreeStr( nom );
        FreeStr( nomob );
        FreeStr( nomsd32 );
        free( val );
    }
    exceptAll {
        CALL_JEDEMA();
        FreeStr( nom );
        FreeStr( nomob );
        FreeStr( nomsd32 );
        free( val );
        raiseException();
    }
    endTry();
    return dico;
}

/* ------------------------------------------------------------------ */
static PyObject *aster_GetResu( PyObject *self, PyObject *args )

/* Construit sous forme d'un dictionnaire Python l'architecture d'une SD resultat

   Arguments :
   IN Nom de la SD resultat
   IN Nature des informations recherchees
   CHAMPS      -> Champs de resultats
   COMPOSANTES -> Liste des composantes des champs
   VARI_ACCES  -> Variables d'acces
   PARAMETRES  -> Parametres


   OUT dico
   Si 'CHAMPS'
   dico['NOM_CHAM'] -> [] si le champ n'est pas calcule
   -> Liste des numeros d'ordre ou le champ est calcule

   Si 'COMPOSANTES'
   dico['NOM_CHAM'] -> [] si le champ n'est pas calcule
   -> Liste des composantes du champ (enveloppe sur tous les instants)

   Si 'VARI_ACCES'
   dico['NOM_VA']   -> Liste des valeurs de la variable d'acces

   Si 'PARAMETRES'
   dico['NOM_VA']   -> Liste des valeurs du parametre

*/
{
    ASTERINTEGER nbchmx, nbpamx, nbord, numch, numva, ier, nbcmp;
    ASTERINTEGER *liord, *ival;
    ASTERINTEGER *val, nbval;
    ASTERDOUBLE *rval;
    char *nomsd, *mode, *liscmp, *nom, *nomsd32, *cmp;
    char *kval, *kvar;
    char *nomch, *nomva;
    int i, lo, ksize = 0, ksizemax = 80, inbord;
    ASTERINTEGER icode, ctype;
    PyObject *dico = NULL, *liste, *key, *value;
    char blanc[80];

    BlankStr( blanc, 80 );

    if ( !PyArg_ParseTuple( args, "ss", &nomsd, &mode ) )
        return NULL;
    nomsd32 = MakeFStrFromCStr( nomsd, 32 );

    /* Identifiant de la SD resultat */
    nbval = 3;
    val = (ASTERINTEGER *)malloc( ( nbval ) * sizeof( ASTERINTEGER ) );
    nom = MakeFStrFromCStr( "LIST_RESULTAT", 24 );

    /* Taille de la SD resultat : nbr champs, nbr paras, nbr numeros d'ordre */
    CALL_JEMARQ();
    try {
        CALL_TAILSD( nom, nomsd32, val, &nbval );
    }
    exceptAll {
        FreeStr( nomsd32 );
        FreeStr( nom );
        free( val );
        raiseException();
    }
    endTry();
    nbchmx = val[0];
    nbpamx = val[1];
    nbord = val[2];
    inbord = (int)nbord;

    if ( strcmp( mode, "CHAMPS" ) == 0 || strcmp( mode, "COMPOSANTES" ) == 0 ) {
        /* Construction du dictionnaire : cle d'acces = nom du champ */
        liord = (ASTERINTEGER *)malloc( inbord * sizeof( ASTERINTEGER ) );
        liscmp = MakeTabFStr( 500, 8 );
        dico = PyDict_New();
        for ( numch = 1; numch <= nbchmx; numch++ ) {
            nomch = MakeBlankFStr( 16 );
            try {
                CALL_RSACCH( nomsd32, &numch, nomch, &nbord, liord, &nbcmp, liscmp );
                inbord = (int)nbord;
                lo = FStrlen( nomch, 16 );
                key = PyUnicode_FromStringAndSize( nomch, lo );
                liste = PyList_New( 0 );
                if ( strcmp( mode, "CHAMPS" ) == 0 ) {
                    for ( i = 0; i < inbord; i++ ) {
                        value = PyLong_FromLong( (long)liord[i] );
                        PyList_Append( liste, value );
                        Py_DECREF( value );
                    }
                }
                if ( strcmp( mode, "COMPOSANTES" ) == 0 ) {
                    for ( i = 0; i < nbcmp; i++ ) {
                        cmp = &( liscmp[i * 8] );
                        lo = FStrlen( cmp, 8 );
                        value = PyUnicode_FromStringAndSize( cmp, lo );
                        PyList_Append( liste, value );
                        Py_DECREF( value );
                    }
                }
                PyDict_SetItem( dico, key, liste );
                Py_DECREF( key );
                Py_DECREF( liste );
                FreeStr( nomch );
            }
            exceptAll {
                FreeStr( nomch );
                raiseException();
            }
            endTry();
        }
        free( liord );
        FreeStr( liscmp );
    } else if ( strcmp( mode, "VARI_ACCES" ) == 0 || strcmp( mode, "PARAMETRES" ) == 0 ) {
        icode = 2;
        if ( strcmp( mode, "VARI_ACCES" ) == 0 ) {
            icode = 0;
        }
        /* Extraction des paramètres ou variables d'accès */
        ival = (ASTERINTEGER *)malloc( inbord * sizeof( ASTERINTEGER ) );
        rval = (ASTERDOUBLE *)malloc( inbord * sizeof( ASTERDOUBLE ) );
        kval = MakeTabFStr( inbord, ksizemax );

        dico = PyDict_New();
        for ( numva = 0; numva <= nbpamx; numva++ ) {
            nomva = MakeBlankFStr( 16 );
            CALL_RSACPA( nomsd32, &numva, &icode, nomva, &ctype, ival, rval, kval, &ier );
            if ( ier != 0 ) {
                FreeStr( nomva );
                continue;
            }

            lo = FStrlen( nomva, 16 );
            key = PyUnicode_FromStringAndSize( nomva, lo );

            liste = PyList_New( 0 );
            if ( ctype < 0 ) {
                /* Erreur */
                PyErr_SetString( PyExc_KeyError, "Type incorrect" );
                return NULL;
            } else if ( ctype == 1 ) {
                for ( i = 0; i < inbord; i++ ) {
                    if ( rval[i] != CALL_R8VIDE() ) {
                        value = PyFloat_FromDouble( (double)rval[i] );
                        PyList_Append( liste, value );
                        Py_DECREF( value );
                    } else {
                        PyList_Append( liste, Py_None );
                    }
                }
            } else if ( ctype == 2 ) {
                for ( i = 0; i < inbord; i++ ) {
                    if ( ival[i] != CALL_ISNNEM() ) {
                        value = PyLong_FromLong( (long)ival[i] );
                        PyList_Append( liste, value );
                        Py_DECREF( value );
                    } else {
                        PyList_Append( liste, Py_None );
                    }
                }
            } else if ( ctype == 4 || ctype == 5 || ctype == 6 || ctype == 7 || ctype == 8 ) {
                switch ( ctype ) {
                case 4:
                    ksize = 8;
                    break;
                case 5:
                    ksize = 16;
                    break;
                case 6:
                    ksize = 24;
                    break;
                case 7:
                    ksize = 32;
                    break;
                case 8:
                    ksize = 80;
                    break;
                }
                for ( i = 0; i < inbord; i++ ) {
                    kvar = kval + i * ksizemax;
                    if ( strncmp( kvar, blanc, ksize ) != 0 ) {
                        value = PyUnicode_FromStringAndSize( kvar, ksize );
                        PyList_Append( liste, value );
                        Py_DECREF( value );
                    } else {
                        PyList_Append( liste, Py_None );
                    }
                }
            }
            PyDict_SetItem( dico, key, liste );
            Py_DECREF( key );
            Py_DECREF( liste );
            FreeStr( nomva );
        }

        free( ival );
        free( rval );
        FreeStr( kval );
    } else {
        PyErr_SetString( PyExc_KeyError, "Mode incorrect" );
        return NULL;
    }
    CALL_JEDEMA();
    FreeStr( nom );
    FreeStr( nomsd32 );
    free( val );
    return dico;
}

/* ------------------------------------------------------------------ */
static PyObject *aster_impers( PyObject *self, PyObject *args ) {
    CALL_IMPERS();
    Py_INCREF( Py_None );
    return Py_None;
}

/* ------------------------------------------------------------------ */
static PyObject *aster_affich( PyObject *self, PyObject *args ) {
    char *texte;
    char *nomfic;

    if ( !PyArg_ParseTuple( args, "ss:affiche", &nomfic, &texte ) )
        return NULL;
    CALLSS( AFFICH, affich, nomfic, texte );

    Py_INCREF( Py_None );
    return Py_None;
}

/* ------------------------------------------------------------------ */
static PyObject *aster_ulopen( PyObject *self, PyObject *args ) {
    char *fichie;
    char *name;
    char *acces;
    char *autor;
    int iunit = 0;
    ASTERINTEGER unit;

    if ( !PyArg_ParseTuple( args, "ssssi:ulopen", &fichie, &name, &acces, &autor, &iunit ) )
        return NULL;
    unit = (ASTERINTEGER)iunit;
    CALL_ULOPEN( &unit, fichie, name, acces, autor );

    Py_INCREF( Py_None );
    return Py_None;
}

/* ------------------------------------------------------------------ */
static PyObject *aster_fclose( PyObject *self, PyObject *args ) {
    int iunit = 0;
    ASTERINTEGER unit;

    if ( !PyArg_ParseTuple( args, "i:fclose", &iunit ) )
        return NULL;
    unit = (ASTERINTEGER)iunit;
    CALL_FCLOSE( &unit );

    Py_INCREF( Py_None );
    return Py_None;
}

/* ------------------------------------------------------------------ */
static PyObject *aster_gcncon( PyObject *self, PyObject *args ) {
    PyObject *res;
    char *type, *Fty, *result;

    if ( !PyArg_ParseTuple( args, "s", &type ) )
        return NULL;
    result = MakeBlankFStr( 8 );
    Fty = MakeFStrFromCStr( type, 1 );
    if ( CALL_ISJVUP() == 1 ) {
        try {
            CALL_GCNCON( Fty, result );
        }
        exceptAll {
            FreeStr( result );
            FreeStr( Fty );
            raiseException();
        }
        endTry();
    }
    res = PyUnicode_FromStringAndSize( result, FStrlen( result, 8 ) );
    FreeStr( result );
    FreeStr( Fty );
    return res;
}

/* ---------------------------------------------------------------------- */
static char rcvale_doc[] =
    "Interface d'appel a la routine fortran RCVALE.\n"
    "   Arguments : nommat, phenomene, nompar, valpar, nomres, stop\n"
    "   Retourne  : valres, codret (tuples)\n"
    " Aucune verification n'est faite sur les arguments d'entree (c'est l'appelant,\n"
    " a priori mater_sdaster.rcvale, qui le fait)";

static PyObject *aster_rcvale( PyObject *self, PyObject *args ) {
    char *nommat, *phenom;
    int istop;
    PyObject *t_nompar, *t_valpar, *t_nomres;
    PyObject *t_valres, *t_codret;
    PyObject *t_res;
    int inbres, inbpar;
    ASTERINTEGER nbpar, nbres, stop;
    char *nompar, *nomres;
    ASTERINTEGER *codret;
    ASTERDOUBLE *valpar, *valres;
    int long_nompar = 8;  /* doivent impérativement correspondre aux  */
    int long_nomres = 16; /* longueurs des chaines de caractères      */
    void *malloc( size_t size );

    if ( !PyArg_ParseTuple( args, "ssOOOi", &nommat, &phenom, &t_nompar, &t_valpar, &t_nomres,
                            &istop ) )
        return NULL;

    /* Conversion en tableaux de chaines et réels */
    inbpar = PyTuple_Size( t_nompar );
    nbpar = (ASTERINTEGER)inbpar;
    nompar = MakeTabFStr( inbpar, long_nompar );
    convertxt( inbpar, t_nompar, nompar, long_nompar );

    valpar = (ASTERDOUBLE *)malloc( inbpar * sizeof( ASTERDOUBLE ) );
    convr8( inbpar, t_valpar, valpar );

    inbres = PyTuple_Size( t_nomres );
    nbres = (ASTERINTEGER)inbres;
    stop = (ASTERINTEGER)istop;
    nomres = MakeTabFStr( inbres, long_nomres );
    convertxt( inbres, t_nomres, nomres, long_nomres );

    /* allocation des variables de sortie */
    valres = (ASTERDOUBLE *)malloc( inbres * sizeof( ASTERDOUBLE ) );
    codret = (ASTERINTEGER *)malloc( inbres * sizeof( ASTERINTEGER ) );

    CALL_RCVALE_WRAP( nommat, phenom, &nbpar, nompar, valpar, &nbres, nomres, valres, codret,
                      &stop );

    /* création des tuples de sortie */
    t_valres = MakeTupleFloat( (long)inbres, valres );
    t_codret = MakeTupleInt( (long)inbres, codret );

    /* retour de la fonction */
    t_res = PyTuple_New( 2 );
    PyTuple_SetItem( t_res, 0, t_valres );
    PyTuple_SetItem( t_res, 1, t_codret );

    FreeStr( nompar );
    free( valpar );
    FreeStr( nomres );
    free( valres );
    free( codret );

    return t_res;
}

/* ---------------------------------------------------------------------- */
static PyObject *aster_gmardm( PyObject *self, PyObject *args ) {
    char *nomgrm, *modele;
    char *Fnom, *Fmod;
    ASTERINTEGER iret;
    PyObject *res;

    if ( !PyArg_ParseTuple( args, "ss", &nomgrm, &modele ) )
        return NULL;

    Fnom = MakeFStrFromCStr( nomgrm, 32 );
    Fmod = MakeFStrFromCStr( modele, 32 );
    CALL_GMARDM( Fnom, Fmod, &iret );

    res = Py_BuildValue( "i", (int)iret );

    FreeStr( Fnom );
    FreeStr( Fmod );

    return res;
}

/* ---------------------------------------------------------------------- */
static char postkutil_doc[] =
    "Interface d'appel a la routine fortran postkutil.\n"
    "   usage: materiau, modelisa = aster.postkutil(imater,resu, fiss) \n\n"
    "     imater : 1 si on veut chercher un materiau, sinon 0"
    "     resu : nom d'une sd_resultat\n"
    "     fiss : nom d'une sd_fond_fiss ou d'une sd_fiss_xfem\n"
    "   Retourne :\n"
    "     materiau : nom d'une sd_mater\n"
    "     modelisa : une chaine parmi 3D, AXIS, D_PLAN_C_PLAN \n"
    "Voir la description de la routine fortran postkutil pour plus de detail\n";

static PyObject *aster_postkutil( PyObject *self, PyObject *args ) {
    int ilmater;
    char *nomres, *nomfis, *repmod, *nommat;
    ASTERINTEGER lmater;
    char *Fres, *Ffis, *Fmod, *Fmat;
    PyObject *res;

    repmod = MakeBlankFStr( 8 );
    nommat = MakeBlankFStr( 8 );
    if ( !PyArg_ParseTuple( args, "iss", &ilmater, &nomres, &nomfis ) )
        return NULL;

    Fres = MakeFStrFromCStr( nomres, 8 );
    Ffis = MakeFStrFromCStr( nomfis, 8 );
    lmater = (ASTERINTEGER)ilmater;
    CALL_POSTKUTIL( &lmater, Fres, Ffis, nommat, repmod );
    Fmat = MakeCStrFromFStr( nommat, 8 );
    Fmod = MakeCStrFromFStr( repmod, 8 );

    res = Py_BuildValue( "ss", Fmat, Fmod );

    FreeStr( Fres );
    FreeStr( Ffis );
    FreeStr( Fmod );
    FreeStr( Fmat );
    return res;
}

/* ---------------------------------------------------------------------- */
static PyObject *aster_mdnoma( PyObject *self, PyObject *args ) {
    PyObject *temp = (PyObject *)0;
    ASTERINTEGER lnomam = 0;
    ASTERINTEGER codret = 0;
    char *nomast, *Fnom;
    char *nomamd;

    if ( !PyArg_ParseTuple( args, "s", &nomast ) )
        return NULL;
    nomamd = MakeBlankFStr( 64 );
    Fnom = MakeFStrFromCStr( nomast, 8 );
    CALL_MDNOMA( nomamd, &lnomam, Fnom, &codret );

    temp = PyUnicode_FromStringAndSize( nomamd, FStrlen( nomamd, (Py_ssize_t)lnomam ) );
    FreeStr( nomamd );
    FreeStr( Fnom );
    return temp;
}

/* ------------------------------------------------------------------ */
static PyObject *aster_mdnoch( PyObject *self, PyObject *args ) {
    PyObject *temp = (PyObject *)0;
    ASTERINTEGER lnochm = 0;
    ASTERINTEGER lresu;
    int ilresu;
    char *noresu;
    char *nomsym;
    char *nopase;
    ASTERINTEGER codret = 0;
    char *nochmd, *n1, *n2, *n3;

    if ( !PyArg_ParseTuple( args, "isss", &ilresu, &noresu, &nomsym, &nopase ) )
        return NULL;
    nochmd = MakeBlankFStr( 64 );
    n1 = MakeFStrFromCStr( noresu, 32 );
    n2 = MakeFStrFromCStr( nomsym, 16 );
    n3 = MakeFStrFromCStr( nopase, 8 );
    lresu = (ASTERINTEGER)ilresu;
    CALL_MDNOCH( nochmd, &lnochm, &lresu, n1, n2, n3, &codret );
    temp = PyUnicode_FromStringAndSize( nochmd, FStrlen( nochmd, (Py_ssize_t)lnochm ) );
    FreeStr( nochmd );
    FreeStr( n1 );
    FreeStr( n2 );
    FreeStr( n3 );
    return temp;
}

/* ------------------------------------------------------------------ */
static PyObject *jeveux_getobjects( PyObject *self, PyObject *args ) {
    ASTERINTEGER nmax, total;
    char *base;
    PyObject *the_list, *pystr;
    char *dummy;
    char *tmp, *ptr;
    int i;

    if ( !PyArg_ParseTuple( args, "s", &base ) )
        return NULL;

    if ( strlen( base ) != 1 ) {
        MYABORT( "le type de base doit etre 1 caractere" );
    }

    dummy = MakeBlankFStr( 24 );
    nmax = 0;
    /* premier appel avec nmax==0 pour connaitre le total */
    CALL_JELST3( base, dummy, &nmax, &total );
    FreeStr( dummy );
    tmp = MakeTabFStr( total, 24 );
    nmax = total;
    /* second appel après allocation mémoire */
    CALL_JELST3( base, tmp, &nmax, &total );

    the_list = PyList_New( (Py_ssize_t)total );
    for ( i = 0, ptr = tmp; i < total; ++i, ptr += 24 ) {
        pystr = PyUnicode_FromStringAndSize( ptr, 24 );
        PyList_SetItem( the_list, i, pystr );
    }
    FreeStr( tmp );
    return the_list;
}

/* ------------------------------------------------------------------ */
static PyObject *jeveux_getattr( PyObject *self, PyObject *args ) {
    PyObject *res;
    char *nomobj, *attr;
    char *charval;
    ASTERINTEGER intval = 0;

    charval = MakeBlankFStr( 32 );
    if ( !PyArg_ParseTuple( args, "ss", &nomobj, &attr ) )
        return NULL;
    CALL_JELIRA( nomobj, attr, &intval, charval );

    res = Py_BuildValue( "is", (int)intval, charval );
    FreeStr( charval );
    return res;
}

static PyObject *jeveux_exists( PyObject *self, PyObject *args ) {
    char *nomobj;
    char *tmpbuf;
    ASTERINTEGER intval = 0;

    if ( !PyArg_ParseTuple( args, "s", &nomobj ) )
        return NULL;
    tmpbuf = MakeFStrFromCStr( nomobj, 32 );
    try {
        CALL_JEEXIN( tmpbuf, &intval );
        FreeStr( tmpbuf );
    }
    exceptAll {
        FreeStr( tmpbuf );
        raiseException();
    }
    endTry();

    if ( intval == 0 ) {
        Py_INCREF( Py_False );
        return Py_False;
    } else {
        Py_INCREF( Py_True );
        return Py_True;
    }
}

/* ------------------------------------------------------------------ */
/*   Routines d'interface pour le catalogue de loi de comportement    */
/* ------------------------------------------------------------------ */
void DEFPSS( LCCREE, lccree, _IN ASTERINTEGER *nbkit, _IN char *lkit, STRING_SIZE llkit,
             _OUT char *compor, STRING_SIZE lcompor ) {
    /*
       Créer un assemblage de LC composé des comportements listés dans 'list_kit'
       et retourne le nom attribué automatiquement à ce comportement.

          CALL LCCREE(NBKIT, LKIT, COMPOR)
          ==> comport = catalc.create(*list_kit)
    */
    PyObject *catalc, *res, *tup_kit;
    const char *scomp;

    catalc = GetJdcAttr( "catalc" );
    /* transforme le tableau de chaines fortran en tuple */
    tup_kit = MakeTupleString( (long)*nbkit, lkit, llkit, NULL );

    res = PyObject_CallMethod( catalc, "create", "O", tup_kit );
    if ( res == NULL ) {
        MYABORT( "Echec lors de la creation du comportement (lccree/create) !" );
    }

    scomp = PyUnicode_AsUTF8( res );
    CopyCStrToFStr( compor, scomp, lcompor );

    Py_XDECREF( res );
    Py_XDECREF( tup_kit );
    Py_XDECREF( catalc );
}

void DEFSS( LCALGO, lcalgo, _IN char *compor, STRING_SIZE lcompor, _OUT char *algo,
            STRING_SIZE lalgo ) {
    /*
       Retourne le premier algorithme d'intégration

          CALL LCALGO(COMPOR, ALGO)
          ==> algo_inte = catalc.get_algo(COMPOR)
    */
    PyObject *catalc, *res;

    catalc = GetJdcAttr( "catalc" );
    res = PyObject_CallMethod( catalc, "get_algo", "s#", compor, lcompor );
    if ( res == NULL ) {
        MYABORT( "Echec lors de la recuperation du premier algorithme "
                 "d'integration (lcalgo/get_algo) !" );
    }

    convertxt( 1, res, algo, lalgo );

    Py_XDECREF( res );
    Py_XDECREF( catalc );
}

void DEFSPPP( LCINFO, lcinfo, _IN char *compor, STRING_SIZE lcompor, _OUT ASTERINTEGER *numlc,
              _OUT ASTERINTEGER *nbvari, _OUT ASTERINTEGER *nbvari_exte ) {
    /*
       Retourne le numéro de routine et le nbre de variables internes

          CALL LCINFO(COMPOR, NUMLC, NBVARI, NBVARI_EXTE)
          ==> num_lc, nb_vari, nb_vari_exte = catalc.get_info(COMPOR)
    */
    PyObject *catalc, *res;

    catalc = GetJdcAttr( "catalc" );
    res = PyObject_CallMethod( catalc, "get_info", "s#", compor, lcompor );
    if ( res == NULL ) {
        MYABORT( "Echec lors de la recuperation des informations sur le "
                 "comportement (lcinfo/get_info) !" );
    }

    *numlc = (ASTERINTEGER)PyLong_AsLong( PyTuple_GetItem( res, 0 ) );
    *nbvari = (ASTERINTEGER)PyLong_AsLong( PyTuple_GetItem( res, 1 ) );
    *nbvari_exte = (ASTERINTEGER)PyLong_AsLong( PyTuple_GetItem( res, 2 ) );

    Py_XDECREF( res );
    Py_XDECREF( catalc );
}

void DEFSPS( LCVARI, lcvari, _IN char *compor, STRING_SIZE lcompor, _IN ASTERINTEGER *nbvari,
             _OUT char *nomvar, STRING_SIZE lnomvar ) {
    /*
       Retourne la liste des variables internes

          CALL LCVARI(COMPOR, NBVARI, LVARI)
          ==> nom_vari = catalc.get_vari(COMPOR)
    */
    PyObject *catalc, *res;

    catalc = GetJdcAttr( "catalc" );
    res = PyObject_CallMethod( catalc, "get_vari", "s#", compor, lcompor );
    if ( res == NULL ) {
        MYABORT( "Echec lors de la recuperation des noms des variables internes du "
                 "comportement (lcvari/get_vari) !" );
    }

    convertxt( (int)*nbvari, res, nomvar, lnomvar );

    Py_XDECREF( res );
    Py_XDECREF( catalc );
}

void DEFSPS( LCEXTEVARI, lcextevari, _IN char *compor, STRING_SIZE lcompor,
             _IN ASTERINTEGER *nbvari, _OUT char *nomvar, STRING_SIZE lnomvar ) {
    /*
       Retourne la liste des variables externes

          CALL LCEXTEVARI(COMPOR, NBVARI, LVARI)
          ==> exte_vari = catalc.get_variexte(COMPOR)
    */
    PyObject *catalc, *res;

    catalc = GetJdcAttr( "catalc" );
    res = PyObject_CallMethod( catalc, "get_variexte", "s#", compor, lcompor );
    if ( res == NULL ) {
        MYABORT( "Echec lors de la recuperation des noms des variables externes du "
                 "comportement (lcvariexte/get_variexte) !" );
    }

    convertxt( (int)*nbvari, res, nomvar, lnomvar );

    Py_XDECREF( res );
    Py_XDECREF( catalc );
}

void DEFSSSP( LCTEST, lctest, _IN char *compor, STRING_SIZE lcompor, _IN char *prop,
              STRING_SIZE lprop, _IN char *valeur, STRING_SIZE lvaleur, _OUT ASTERINTEGER *iret ) {
    /*
       Est-ce que VALEUR est un valeur autorisée de PROPRIETE ?
             CALL LCTEST(COMPOR, PROPRIETE, VALEUR, IRET)
             ==> iret = catalc.query(COMPOR, PROPRIETE, VALEUR)
    */
    PyObject *catalc, *res;

    catalc = GetJdcAttr( "catalc" );
    res = PyObject_CallMethod( catalc, "query", "s#s#s#", compor, lcompor, prop, lprop, valeur,
                               lvaleur );
    if ( res == NULL ) {
        MYABORT( "Echec lors du test d'une propriete du comportement (lctest/query) !" );
    }

    *iret = (ASTERINTEGER)PyLong_AsLong( res );

    Py_XDECREF( res );
    Py_XDECREF( catalc );
}

void DEFSS( LCTYPE, lctype, _IN char *compor, STRING_SIZE lcompor, _OUT char *typ,
            STRING_SIZE ltyp ) {
    /*
       Retourne le type de comportement
             CALL LCTYPE(COMPOR, TYPE)
             ==> ldctype = catalc.get_type(COMPOR)
    */
    PyObject *catalc, *res;
    const char *styp;

    catalc = GetJdcAttr( "catalc" );
    res = PyObject_CallMethod( catalc, "get_type", "s#", compor, lcompor );
    if ( res == NULL ) {
        MYABORT( "Echec lors du test d'une propriete du comportement (lctype/get_type) !" );
    }

    styp = PyUnicode_AsUTF8( res );
    CopyCStrToFStr( typ, styp, ltyp );

    Py_XDECREF( res );
    Py_XDECREF( catalc );
}

void DEFS( LCDISCARD, lcdiscard, _IN char *compor, STRING_SIZE lcompor ) {
    /*
       Supprime la loi "de travail"
       Si compor=' ', on supprime toutes les lois de travail

           CALL LCDISCARD(COMPOR)
           ==> catalc.discard(COMPOR)
    */
    PyObject *catalc, *res;

    catalc = GetJdcAttr( "catalc" );
    if ( compor[0] == ' ' ) {
        res = PyObject_CallMethod( catalc, "discard", NULL );
    } else {
        res = PyObject_CallMethod( catalc, "discard", "s#", compor, lcompor );
    }
    if ( res == NULL ) {
        MYABORT( "Echec lors de la suppression des comportements (lcdiscard) !" );
    }

    Py_XDECREF( res );
    Py_XDECREF( catalc );
}

void DEFSS( LCSYMB, lcsymb, _IN char *compor, STRING_SIZE lcompor, _OUT char *name,
            STRING_SIZE lname ) {
    /*
       Retourne le nom de la fonction dans la bibliothèque MFront

          CALL LCFUNC(COMPOR, NAME)
          ==> name = catalc.get_symbol(COMPOR)
    */
    PyObject *catalc, *res;

    catalc = GetJdcAttr( "catalc" );
    res = PyObject_CallMethod( catalc, "get_symbol", "s#", compor, lcompor );
    if ( res == NULL ) {
        MYABORT( "Echec lors de la recuperation du nom de la fonction "
                 "d'integration dans MFront (lcsymb/get_symbol) !" );
    }
    convertxt( 1, res, name, lname );

    Py_XDECREF( res );
    Py_XDECREF( catalc );
}

void DEFSS( LCSYMM, lcsymm, _IN char *compor, STRING_SIZE lcompor, _OUT char *symm,
            STRING_SIZE lsymm ) {
    /*
       Retourne le type de symétrie de la matrice

          CALL LCSYMM(COMPOR, SYMMETRY)
          ==> name = catalc.get_symmetry(COMPOR)
    */
    PyObject *catalc, *res;

    catalc = GetJdcAttr( "catalc" );
    res = PyObject_CallMethod( catalc, "get_symmetry", "s#", compor, lcompor );
    if ( res == NULL ) {
        MYABORT( "Echec lors de la récuperation de la symétrie de la matrice !" );
    }
    convertxt( 1, res, symm, lsymm );

    Py_XDECREF( res );
    Py_XDECREF( catalc );
}

void DEFSS( LCDEFORMLDC, lcdeformldc, _IN char *compor, STRING_SIZE lcompor, _OUT char *deformldc,
            STRING_SIZE ldeformldc ) {
    /*
       Retourne la nature de la déformation en entrée de la ldc

          CALL LCDEFORMLDC(COMPOR, DEFORM_LDC)
          ==> deform_ldc = catalc.get_deformldc(COMPOR)
    */
    PyObject *catalc, *res;

    catalc = GetJdcAttr( "catalc" );
    res = PyObject_CallMethod( catalc, "get_deformldc", "s#", compor, lcompor );
    if ( res == NULL ) {
        MYABORT( "Echec lors de la récuperation de la nature de la  "
                 "déformation en entrée de la ldc (lcdeformldc/get_deformldc) !" );
    }
    convertxt( 1, res, deformldc, ldeformldc );

    Py_XDECREF( res );
    Py_XDECREF( catalc );
}

void DEFSS( LCREGUVISC, lcreguvisc, _IN char *compor, STRING_SIZE lcompor, _OUT char *reguvisc,
            STRING_SIZE lreguvisc ) {
    /*
       Retourne la nature de la régularisation visqueuse

          CALL LCREGUVISC(COMPOR, REGU_VISC)
          ==> regu_visc = catalc.get_reguvisc(COMPOR)
    */
    PyObject *catalc, *res;

    catalc = GetJdcAttr( "catalc" );
    res = PyObject_CallMethod( catalc, "get_reguvisc", "s#", compor, lcompor );
    if ( res == NULL ) {
        MYABORT( "Echec lors de la récuperation de la nature de la  "
                 "régularisation visqueuse !" );
    }
    convertxt( 1, res, reguvisc, lreguvisc );

    Py_XDECREF( res );
    Py_XDECREF( catalc );
}

void DEFPSS( LCKITREAD, lckitread, _IN ASTERINTEGER *nbkit, _IN char *lkit, STRING_SIZE llkit,
             _OUT char *lrela, STRING_SIZE llrela ) {
    PyObject *catalc, *tup_lkit, *res;
    int nb_tuple;
    catalc = GetJdcAttr( "catalc" );

    // Input in tuple: first => name of kit, next => components of kit
    tup_lkit = MakeTupleString( (long)*nbkit, lkit, llkit, NULL );

    // Call C function
    res = PyObject_CallMethod( catalc, "get_kit", "O", tup_lkit );

    if ( res == NULL ) {
        MYABORT( "Echec lors de la lecture du kit (lckitread/get_kit) !" );
        DEBUG_ASSERT( PyTuple_Check( res ) );
    }
    nb_tuple = PyTuple_Size( tup_lkit );
    if ( nb_tuple != 5 ) {
        MYABORT( "Echec lors de la lecture du kit (lckitread/get_kit) !" );
        DEBUG_ASSERT( PyTuple_Check( res ) );
    }

    // Convert output
    convertxt( 4, res, lrela, llrela );

    Py_XDECREF( res );
    Py_XDECREF( tup_lkit );
    Py_XDECREF( catalc );
}

void DEFPS( INT_TO_STRING_CONVERSION, int_to_string_conversion, ASTERINTEGER *value, char *out,
            STRING_SIZE sout ) {
    char tmp[9];
    sprintf( tmp, "%-8ld", *value );
    CopyCStrToFStr( out, tmp, 8 );
    return;
}

/* ----------   FIN catalogue de loi de comportement   -------------- */
/* ------------------------------------------------------------------ */

/* ------------------------------------------------------------------ */
static PyObject *aster_argv( _UNUSED PyObject *self, _IN PyObject *args ) {
    Py_INCREF( Py_None );
    return Py_None;
}

/* List of functions defined in the module */
static PyMethodDef aster_methods[] = {
    // {"onFatalError", aster_onFatalError, METH_VARARGS},
    { "fclose", aster_fclose, METH_VARARGS },
    { "ulopen", aster_ulopen, METH_VARARGS },
    { "affiche", aster_affich, METH_VARARGS },
    { "impers", aster_impers, METH_VARARGS },
    { "mdnoma", aster_mdnoma, METH_VARARGS },
    { "mdnoch", aster_mdnoch, METH_VARARGS },
    { "rcvale", aster_rcvale, METH_VARARGS, rcvale_doc },
    { "gmardm", aster_gmardm, METH_VARARGS },
    { "postkutil", aster_postkutil, METH_VARARGS, postkutil_doc },
    { "argv", aster_argv, METH_VARARGS },
    { "getvectjev", aster_getvectjev, METH_VARARGS, getvectjev_doc },
    { "getcolljev", aster_getcolljev, METH_VARARGS, getcolljev_doc },
    { "GetResu", aster_GetResu, METH_VARARGS },
    { "jeveux_getobjects", jeveux_getobjects, METH_VARARGS },
    { "jeveux_getattr", jeveux_getattr, METH_VARARGS },
    { "jeveux_exists", jeveux_exists, METH_VARARGS },
    { "get_nom_concept_unique", aster_gcncon, METH_VARARGS },
    { NULL, NULL } /* sentinel */
};

#ifndef ASTER_WITHOUT_PYMOD
/* Initialization function for the module (*must* be called initaster) */
static char aster_module_documentation[] = "C implementation of the Python aster module\n"
                                           "\n";

static struct PyModuleDef aster_def = { PyModuleDef_HEAD_INIT,
                                        "aster",
                                        aster_module_documentation,
                                        -1,
                                        aster_methods,
                                        NULL,
                                        NULL,
                                        NULL,
                                        NULL };

PyObject *PyInit_aster( void ) {
    PyObject *aster = (PyObject *)0;

    /* Create the module and add the functions */
    aster = PyModule_Create( &aster_def );

    init_etape_stack();
    return aster;
}
#endif
