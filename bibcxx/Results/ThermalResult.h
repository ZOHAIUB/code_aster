#ifndef EVOLUTIVETHERMALLOAD_H_
#define EVOLUTIVETHERMALLOAD_H_

/**
 * @file ThermalResult.h
 * @brief Fichier entete de la classe ThermalResult
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

#include "astercxx.h"

#include "Results/TransientResult.h"
#include "Supervis/ResultNaming.h"

/**
 * @class ThermalResult
 * @brief Cette classe correspond a un evol_ther
 * @author Nicolas Sellenet
 */
class ThermalResult : public TransientResult {
  public:
    /**
     * @brief Constructeur
     */
    ThermalResult( const std::string name = ResultNaming::getNewResultName() )
        : TransientResult( name, "EVOL_THER" ) {};
};

/**
 * @typedef ThermalResultPtr
 * @brief Pointeur intelligent vers un ThermalResult
 */
typedef std::shared_ptr< ThermalResult > ThermalResultPtr;

#endif /* EVOLUTIVETHERMALLOAD_H_ */
