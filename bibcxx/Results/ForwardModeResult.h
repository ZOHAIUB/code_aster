#ifndef FORWARDMECHANICALMODERESULT_H_
#define FORWARDMECHANICALMODERESULT_H_

/**
 * @file ForwardModeResult.h
 * @brief Fichier entete de la classe ForwardModeResult
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

class ModeResult;
typedef std::shared_ptr< ModeResult > ModeResultPtr;

/**
 * @class ForwardModeResultPtr
 * @brief Forward definition of ModeResultPtr
 * @author Nicolas Sellenet
 */
class ForwardModeResultPtr {
  private:
    /** @brief Pointer to ModeResult */
    ModeResultPtr _ptr;
    bool _isSet;

  public:
    /**
     * @brief Constructor
     */
    ForwardModeResultPtr();

    /**
     * @brief Constructor
     */
    ForwardModeResultPtr( const ModeResultPtr &ptr );

    void operator=( const ModeResultPtr &ptr );

    ModeResultPtr getPointer();

    bool isSet() const;

    void setPointer( const ModeResultPtr &ptr );
};

#endif /* FORWARDMECHANICALMODERESULT_H_ */
