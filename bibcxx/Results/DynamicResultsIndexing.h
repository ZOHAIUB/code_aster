#ifndef DYNAMICRESULTSINDEXING_H_
#define DYNAMICRESULTSINDEXING_H_

/**
 * @file DynamicResultsIndexing.h
 * @brief Fichier entete de la classe DynamicResultsIndexing
 * @author Natacha Béreux
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

/* person_in_charge: natacha.bereux at edf.fr */

#include "astercxx.h"

#include "MemoryManager/JeveuxCollection.h"
#include "MemoryManager/JeveuxVector.h"
/**
 * @class DynamicResultsIndexing
 * @brief Cette classe correspond aux objets Aster permettant d'indexer les résultats "généralisés"
 * i.e. les résultats sur base modale.
 * @author Natacha Béreux
 */
class DynamicResultsIndexing {
  private:
    /** @brief Collection '.REFD' */
    /** Collection de vecteurs contenant les informations sur les matrices de masse, amortissement
        et rigidité mais également d’autres informations telles que la numérotation (NUME_DDL)
        ou encore la liste d’interfaces dynamiques et statiques ayant servi au calcul produisant le
       concept*/
    JeveuxCollectionChar24 _refd;

    /** @brief Vecteur Jeveux '.INDI' */
    /** Cet objet définit l’indirection entre le numéro d’ordre des champs enregistrés et les objets
        enregistrés dans  la collection .REFD */
    JeveuxVectorLong _indi;

  public:
    /**
     * @brief Constructeur
     */
    DynamicResultsIndexing( const std::string &name )
        : _refd( JeveuxCollectionChar24( name + ".REFD" ) ),
          _indi( JeveuxVectorLong( name + ".INDI" ) ) {};
};

#endif /* DYNAMICRESULTSINDEXING_H_ */
