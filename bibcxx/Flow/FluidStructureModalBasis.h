#ifndef FLUIDSTRUCTMODALBASIS_H_
#define FLUIDSTRUCTMODALBASIS_H_

/**
 * @file FluidStructureModalBasis.h
 * @brief Fichier entete de la classe FluidStructureModalBasis
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include "astercxx.h"

#include "DataFields/ListOfTables.h"
#include "DataStructures/DataStructure.h"
#include "MemoryManager/JeveuxVector.h"
#include "Supervis/ResultNaming.h"

/**
 * @class FluidStructureModalBasis
 * @brief Cette classe correspond a une sd_melasflu
 * @author Natacha Béreux
 */
class FluidStructureModalBasis : public DataStructure, public ListOfTables {
  private:
    /** @brief Objet Jeveux '.REMF' */
    JeveuxVectorChar8 _remf;
    /** @brief Objet Jeveux '.DESC' */
    JeveuxVectorChar16 _desc;
    /** @brief Objet Jeveux '.FACT' */
    JeveuxVectorReal _fact;
    /** @brief Objet Jeveux '.FREQ' */
    JeveuxVectorReal _freq;
    /** @brief Objet Jeveux '.MASG' */
    JeveuxVectorReal _masg;
    /** @brief Objet Jeveux '.NUMO' */
    JeveuxVectorLong _numo;
    /** @brief Objet Jeveux '.VITE' */
    JeveuxVectorReal _vite;

  public:
    /**
     * @typedef FluidStructureModalBasisPtr
     * @brief Pointeur intelligent vers un FluidStructureModalBasis
     */
    typedef std::shared_ptr< FluidStructureModalBasis > FluidStructureModalBasisPtr;

    /**
     * @brief Constructeur
     */
    FluidStructureModalBasis() : FluidStructureModalBasis( ResultNaming::getNewResultName() ) {};
    /**
     * @brief Constructeur
     */
    FluidStructureModalBasis( const std::string name )
        : DataStructure( name, 8, "MELASFLU" ),
          ListOfTables( name ),
          _remf( JeveuxVectorChar8( getName() + ".REMF" ) ),
          _desc( JeveuxVectorChar16( getName() + ".DESC" ) ),
          _fact( JeveuxVectorReal( getName() + ".FACT" ) ),
          _freq( JeveuxVectorReal( getName() + ".FREQ" ) ),
          _masg( JeveuxVectorReal( getName() + ".MASG" ) ),
          _numo( JeveuxVectorLong( getName() + ".NUMO" ) ),
          _vite( JeveuxVectorReal( getName() + ".VITE" ) ) {};
};

/**
 * @typedef FluidStructureModalBasisPtr
 * @brief Pointeur intelligent vers un FluidStructureModalBasis
 */
typedef std::shared_ptr< FluidStructureModalBasis > FluidStructureModalBasisPtr;

#endif /* FLUIDSTRUCTMODALBASIS_H_ */
