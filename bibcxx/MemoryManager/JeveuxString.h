#ifndef JEVEUXSTRING_H_
#define JEVEUXSTRING_H_

/**
 * @file JeveuxString.h
 * @brief Definition d'une chaine a la maniere Fortran (sans \0 a la fin)
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

#include "astercxx.h"

#include <stdexcept>
#include <string>

#include <string.h>

/**
 * @class JeveuxString
 * @brief Cette classe template permet de definir une chaine fortran rapidement manipulable
 * @author Nicolas Sellenet
 */
template < int lengthT >
class JeveuxString {
  private:
    /** @brief Pointeur vers la chaine de caractere */
    char currentValue[lengthT];

    /**
     * @brief Recopie securisee d'un char*
     * @param chaine char*
     * @param chaine Taille de la chaine a recopier
     */
    inline void safeCopyFromChar( const char *chaine, const int size ) {
        if ( size < 1 )
            return;
        if ( size < lengthT ) {
            memset( &currentValue, ' ', sizeof( char ) * lengthT );
            memcpy( &currentValue, chaine, sizeof( char ) * size );
        } else {
            memcpy( &currentValue, chaine, sizeof( char ) * lengthT );
        }
    };

  public:
    /**
     * @brief Constructeur par defaut
     */
    inline JeveuxString() {};

    /**
     * @brief Constructeur a partir d'un char*
     * @param chaine Chaine a recopier en style fortran
     */
    inline JeveuxString( const JeveuxString< lengthT > &chaine ) {
        memcpy( &currentValue, chaine.c_str(), sizeof( char ) * lengthT );
    };

    /**
     * @brief Constructeur rapide a partir d'un char*
     * @param chaine Chaine a recopier en style fortran
     * @param size Taille de la chaine a recopier
     */
    inline JeveuxString( const char *chaine, const int size ) { safeCopyFromChar( chaine, size ); };

    /**
     * @brief Constructeur a partir d'un char*
     * @param chaine Chaine a recopier en style fortran
     */
    inline JeveuxString( const char *chaine ) { safeCopyFromChar( chaine, strlen( chaine ) ); };

    /**
     * @brief Constructeur a partir d'un char*
     * @param chaine Chaine a recopier en style fortran
     */
    inline JeveuxString( const std::string str ) { safeCopyFromChar( str.c_str(), str.length() ); };

    /**
     * @brief Surcharge de l'operateur = pour une affectation rapide
     * @param chaine Recopie a partir d'un JeveuxString
     * @return Reference vers la chaine recopiee
     */
    inline JeveuxString &operator=( const JeveuxString< lengthT > &chaine ) {
        memcpy( &currentValue, &( chaine.currentValue ), sizeof( char ) * lengthT );
        return *this;
    };

    /**
     * @brief Surcharge de l'operateur = pour une affectation a partir d'un char*
     * @param chaine Recopie a partir d'un char*
     * @return reference vers la chaine recopiee
     */
    inline JeveuxString &operator=( const char *chaine ) {
        safeCopyFromChar( chaine, strlen( chaine ) );
        return *this;
    };

    /**
     * @brief Surcharge de l'operateur = pour une affectation a partir d'une string
     * @param chaine Recopie a partir d'une string
     * @return reference vers la chaine recopiee
     */
    inline JeveuxString &operator=( const std::string &chaine ) {
        safeCopyFromChar( chaine.c_str(), chaine.size() );
        return *this;
    };

    /**
     * @brief Operator ==
     */
    inline bool operator==( const JeveuxString< lengthT > &chaine ) {
        int ret = strncmp( this->c_str(), chaine.c_str(), lengthT );
        if ( ret == 0 )
            return true;
        return false;
    };

    /**
     * @brief Fonction permettant d'obtenir le pointeur vers le debut de la chaine
     * @return Pointeur vers le debut de la chaine. Attention pas de \0 a la fin !!
     */
    inline const char *c_str() { return currentValue; };

    /**
     * @brief Fonction permettant d'obtenir le pointeur vers le debut de la chaine
     * @return Pointeur vers le debut de la chaine. Attention pas de \0 a la fin !!
     */
    inline const char *c_str() const { return currentValue; };

    std::string rstrip() const {
        std::string buff( currentValue, lengthT );
        std::string whitespaces( " \t\f\v\n\r" );
        std::size_t found = buff.find_last_not_of( whitespaces );
        if ( found != std::string::npos )
            buff.erase( found + 1 );
        else
            buff.clear();
        return buff;
    };

    /**
     * @brief Fonction renvoyant la taille de la chaine
     * @return un entier
     */
    inline int size() const { return lengthT; };

    /**
     * @brief Fonction renvoyant la taille de la chaine (équivalent de 'size').
     * @return un entier
     */
    inline int length() const { return lengthT; };

    /**
     * @brief Fonction renvoyant une string correspondant a une recopie de l'objet
     * @return string contenant la chaine contenue dans l'objet
     */
    inline std::string toString() const { return std::string( currentValue, lengthT ); };

    /**
     * @brief Unsafe fast copy from a char*
     * @param chaine String to copy
     */
    inline void unsafeFastCopy( const char *chaine ) {
#ifndef NDEBUG
        if ( strlen( chaine ) < lengthT )
            throw std::runtime_error( "String size error" );
#endif
        memcpy( &currentValue, chaine, sizeof( char ) * lengthT );
    };

    inline operator std::string() const { return toString(); };

    /** @brief overload << operator */
    friend std::ostream &operator<<( std::ostream &os, const JeveuxString< lengthT > &toPrint ) {
        os << toPrint.toString();

        return os;
    };

    inline const char &operator[]( const ASTERINTEGER &i ) const {
#ifdef ASTER_DEBUG_CXX_OBJECTS
        if ( i < 0 && i >= this->size() ) {
            std::string error = "Out of range of JeveuxString, index = " + std::to_string( i ) +
                                " ( size = " + std::to_string( this->size() ) + " )";
            AS_ABORT( error );
        }
#endif
        return currentValue[i];
    };

    /**
     * @brief Surcharge de l'operateur [] sans const (pour les lvalue)
     * @param i Indice dans le tableau Jeveux
     * @return la valeur du tableau Jeveux a la position i
     */
    inline char &operator[]( const ASTERINTEGER &i ) {
#ifdef ASTER_DEBUG_CXX_OBJECTS
        if ( i < 0 && i >= this->size() ) {
            std::string error = "Out of range of JeveuxString, index = " + std::to_string( i ) +
                                " ( size = " + std::to_string( this->size() ) + " )";
            AS_ABORT( error );
        }
#endif
        return currentValue[i];
    };
};

/** @typedef Definition d'une chaine Jeveux de longueur 8 */
typedef JeveuxString< 8 > JeveuxChar8;
/** @typedef Definition d'une chaine Jeveux de longueur 16 */
typedef JeveuxString< 16 > JeveuxChar16;
/** @typedef Definition d'une chaine Jeveux de longueur 24 */
typedef JeveuxString< 24 > JeveuxChar24;
/** @typedef Definition d'une chaine Jeveux de longueur 32 */
typedef JeveuxString< 32 > JeveuxChar32;
/** @typedef Definition d'une chaine Jeveux de longueur 80 */
typedef JeveuxString< 80 > JeveuxChar80;

#endif
