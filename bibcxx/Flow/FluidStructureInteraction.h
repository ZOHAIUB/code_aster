#ifndef FLUIDSTRUCTUREINTERACTION_H_
#define FLUIDSTRUCTUREINTERACTION_H_

/**
 * @file FluidStructureInteraction.h
 * @brief Fichier entete de la classe FluidStructureInteraction
 * @author Nicolas Sellenet
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

#include "DataStructures/DataStructure.h"
#include "MemoryManager/JeveuxVector.h"
#include "Supervis/ResultNaming.h"

/**
 * @class FluidStructureInteraction
 * @brief Cette classe correspond a une sd spectre
 * @author Nicolas Sellenet
 */
class FluidStructureInteraction : public DataStructure {
  private:
    /** @brief Objet Jeveux '.FSIC' */
    JeveuxVectorLong _fsic;
    /** @brief Objet Jeveux '.FSVI' */
    JeveuxVectorLong _fsvi;
    /** @brief Objet Jeveux '.FSVK' */
    JeveuxVectorChar8 _fsvk;
    /** @brief Objet Jeveux '.FSVR' */
    JeveuxVectorReal _fsvr;
    /** @brief Objet Jeveux '.FSGM' */
    JeveuxVectorChar24 _fsgm;
    /** @brief Objet Jeveux '.FSGR' */
    JeveuxVectorReal _fsgr;
    /** @brief Objet Jeveux '.FSCR' */
    JeveuxVectorReal _fscr;
    /** @brief Objet Jeveux '.UNIT_FAISCEAU' */
    JeveuxVectorLong _unitFaisceau;
    /** @brief Objet Jeveux '.UNIT_GRAPPES' */
    JeveuxVectorLong _unitGrappes;

  public:
    /**
     * @typedef FluidStructureInteractionPtr
     * @brief Pointeur intelligent vers un FluidStructureInteraction
     */
    typedef std::shared_ptr< FluidStructureInteraction > FluidStructureInteractionPtr;

    /**
     * @brief Constructeur
     */
    FluidStructureInteraction() : FluidStructureInteraction( ResultNaming::getNewResultName() ) {};
    /**
     * @brief Constructeur
     */
    FluidStructureInteraction( const std::string name )
        : DataStructure( name, 8, "TYPE_FLUI_STRU" ),
          _fsic( JeveuxVectorLong( getName() + "           .FSIC" ) ),
          _fsvi( JeveuxVectorLong( getName() + "           .FSVI" ) ),
          _fsvk( JeveuxVectorChar8( getName() + "           .FSVK" ) ),
          _fsvr( JeveuxVectorReal( getName() + "           .FSVR" ) ),
          _fsgm( JeveuxVectorChar24( getName() + "           .FSGM" ) ),
          _fsgr( JeveuxVectorReal( getName() + "           .FSGR" ) ),
          _fscr( JeveuxVectorReal( getName() + "           .FSCR" ) ),
          _unitFaisceau( JeveuxVectorLong( getName() + ".UNIT_FAISCEAU" ) ),
          _unitGrappes( JeveuxVectorLong( getName() + ".UNIT_GRAPPES" ) ) {};
};

/**
 * @typedef FluidStructureInteractionPtr
 * @brief Pointeur intelligent vers un FluidStructureInteraction
 */
typedef std::shared_ptr< FluidStructureInteraction > FluidStructureInteractionPtr;

#endif /* FLUIDSTRUCTUREINTERACTION_H_ */
