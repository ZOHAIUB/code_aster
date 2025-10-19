#ifndef JEVEUXALLOWEDTYPES_H_
#define JEVEUXALLOWEDTYPES_H_

/**
 * @file JeveuxAllowedTypes.h
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

#include "astercxx.h"

#include "asterc_config.h"

#include "MemoryManager/JeveuxString.h"

#include <complex>

/**
 * @enum JeveuxTypes
 * @brief Fournit tous les types autorises dans le gestionnaire memoire Jeveux
 */
enum JeveuxTypes {
    Integer,
    Integer4,
    Real,
    Complex,
    Char8,
    Char16,
    Char24,
    Char32,
    Char80,
    Logical
};
/**
 * @def JeveuxTypesNames
 * @brief Fournit sous forme de chaine les types Jeveux existant
 */
static const std::string JeveuxTypesNames[10] = { "I",   "S",   "R",   "C",   "K8",
                                                  "K16", "K24", "K32", "K80", "L" };

/**
 * @struct AllowedJeveuxType
 * @brief Structure template permettant de limiter le type instanciable de JeveuxVectorClass
 * @tparam T Type autorise
 */
template < typename T >
struct AllowedJeveuxType; // undefined for bad types!

template <>
struct AllowedJeveuxType< ASTERINTEGER > {
    static const unsigned short numTypeJeveux = Integer;
    typedef ASTERINTEGER type;
};

template <>
struct AllowedJeveuxType< ASTERINTEGER4 > {
    static const unsigned short numTypeJeveux = Integer4;
    typedef ASTERINTEGER4 type;
};

template <>
struct AllowedJeveuxType< ASTERDOUBLE > {
    static const unsigned short numTypeJeveux = Real;
    typedef ASTERDOUBLE type;
};

template <>
struct AllowedJeveuxType< ASTERCOMPLEX > {
    static const unsigned short numTypeJeveux = Complex;
    typedef ASTERCOMPLEX type;
};

template <>
struct AllowedJeveuxType< JeveuxChar8 > {
    static const unsigned short numTypeJeveux = Char8;
    typedef JeveuxChar8 type;
};

template <>
struct AllowedJeveuxType< JeveuxChar16 > {
    static const unsigned short numTypeJeveux = Char16;
    typedef JeveuxChar16 type;
};

template <>
struct AllowedJeveuxType< JeveuxChar24 > {
    static const unsigned short numTypeJeveux = Char24;
    typedef JeveuxChar24 type;
};

template <>
struct AllowedJeveuxType< JeveuxChar32 > {
    static const unsigned short numTypeJeveux = Char32;
    typedef JeveuxChar32 type;
};

template <>
struct AllowedJeveuxType< JeveuxChar80 > {
    static const unsigned short numTypeJeveux = Char80;
    typedef JeveuxChar80 type;
};

#if ASTER_LOGICAL_SIZE != 1
#error Size for ASTER_LOGICAL_SIZE (!= 1) not allowed
#endif

template <>
struct AllowedJeveuxType< bool > {
    static const unsigned short numTypeJeveux = Logical;
    typedef bool type;
};

#endif /* JEVEUXALLOWEDTYPES_H_ */
