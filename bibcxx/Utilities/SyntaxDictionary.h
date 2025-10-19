#ifndef SYNTAXDICTIONARY_H_
#define SYNTAXDICTIONARY_H_

/**
 * @file SyntaxDictionary.h
 * @brief Fichier entete de la struct SyntaxMapContainer
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

#include "astercxx.h"

#include "aster_pybind.h"

#include <typeinfo>

struct SyntaxMapContainer;

/**
 * @typedef ListSyntaxMapContainer
 * @brief Type permettant de definir les occurence des mots-cles facteurs
 */
typedef std::list< SyntaxMapContainer > ListSyntaxMapContainer;
typedef ListSyntaxMapContainer::iterator ListSyntaxMapContainerIter;
typedef ListSyntaxMapContainer::const_iterator ListSyntaxMapContainerCIter;

/**
 * @class SyntaxMapContainer
 * @brief Cette struct decrit un dictionnaire permettant de contenir la syntaxe des commandes Aster
 * @author Nicolas Sellenet
 */
class SyntaxMapContainer {
  public:
    /** @brief Typedef definissant divers types */

    using VariantSyntaxType =
        std::variant< ASTERINTEGER, ASTERDOUBLE, ASTERCOMPLEX, std::string, VectorLong,
                      VectorString, VectorReal, VectorComplex, ListSyntaxMapContainer >;

    /** @brief Typedef definissant un map associant une chaine a divers types */
    typedef std::map< std::string, VariantSyntaxType > SyntaxMap;
    typedef SyntaxMap::iterator SyntaxMapIter;
    typedef SyntaxMap::const_iterator SyntaxMapCIter;

    /** @brief Conteneur a proprement parler */
    SyntaxMap container;

    /**
     * @brief Opérateur +=
     * @param toAdd SyntaxMapContainer à ajouter
     * @return reference to the current object
     */
    SyntaxMapContainer &operator+=( const SyntaxMapContainer &toAdd ) {
        container.insert( toAdd.container.begin(), toAdd.container.end() );
        return *this;
    };

    /**
     * @brief Convertisseur du conteneur en dictionnaire python
     * @return un dict python contenant la syntaxe valorisable par l'objet CommandSyntax
     * @todo Probleme de refcounting : ajouter un objet wrapper qui se chargera de la destruction
     * @todo seul CommandSyntax devrait pouvoir appeler cette fonction ?
     * @todo ajouter un const pour this
     */
    PyObject *convertToPythonDictionnary( PyObject *returnDict = NULL ) const;
    py::dict convertToPyDict() const;

  private:
    friend class CommandSyntax;
    friend class MaterialField;
    friend SyntaxMapContainer operator+( const SyntaxMapContainer &, const SyntaxMapContainer & );
};

/**
 * @brief Opérateur +
 * @param toAdd1 SyntaxMapContainer à ajouter
 * @param toAdd2 SyntaxMapContainer à ajouter
 * @return SyntaxMapContainer résultat
 */
SyntaxMapContainer operator+( const SyntaxMapContainer &toAdd1, const SyntaxMapContainer &toAdd2 );

#endif /* SYNTAXDICTIONARY_H_ */
