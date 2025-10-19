# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
# This file is part of code_aster.
#
# code_aster is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# code_aster is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
# --------------------------------------------------------------------

from cataelem.Tools.base_objects import LocatedComponents, ArrayOfComponents, SetOfNodes, ElrefeLoc
from cataelem.Tools.base_objects import Calcul, Element
import cataelem.Commons.physical_quantities as PHY
import cataelem.Commons.located_components as LC
import cataelem.Commons.parameters as SP
import cataelem.Commons.mesh_types as MT
from cataelem.Options.options import OP
import cataelem.Commons.attributes as AT

# --------------------------------------------------------------------------------------------------
# Located components
# --------------------------------------------------------------------------------------------------
DDL_THER = LocatedComponents(
    phys=PHY.TEMP_R, type="ELNO", diff=True, components=(("EN1", ("HHO_FT[4]",)), ("EN2", ()))
)

CHHOBS = LocatedComponents(
    phys=PHY.N3600R, type="ELNO", diff=True, components=(("EN1", ("X[10]",)), ("EN2", ()))
)

MVECTTR = ArrayOfComponents(phys=PHY.VTEM_R, locatedComponents=DDL_THER)

MMATTTR = ArrayOfComponents(phys=PHY.MTEM_R, locatedComponents=DDL_THER)


# --------------------------------------------------------------------------------------------------
class THER_2D_HHO3_F(Element):
    """Thermics - Skin element HHO_QUAD - PLAN - SEG"""

    meshType = MT.SEG3
    nodes = (SetOfNodes("EN1", (3,)), SetOfNodes("EN2", (1, 2)))
    attrs = ((AT.BORD_ISO, "OUI"),)
    elrefe = (ElrefeLoc(MT.SE3, gauss=("RIGI=FPG4",), mater=("RIGI",)),)
    calculs = (
        OP.CHAR_THER_FLUN_F(
            te=461,
            para_in=(
                (SP.PFLUXNF, LC.CFLUXNF),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PINSTR, LC.CTIMETR),
                (OP.CHAR_THER_FLUN_F.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_FLUN_R(
            te=461,
            para_in=(
                (SP.PFLUXNR, LC.CFLUXNR),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PINSTR, LC.CTIMETR),
                (OP.CHAR_THER_FLUN_R.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_FLUNL(
            te=461,
            para_in=(
                (SP.PFLUXNL, LC.CFLUXNF),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
                (OP.CHAR_THER_FLUNL.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_ECHA_F(
            te=461,
            para_in=(
                (SP.PCOEFHF, LC.CHECHPF),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
                (SP.PT_EXTF, LC.CTEMPEF),
                (OP.CHAR_THER_ECHA_F.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_ECHA_R(
            te=461,
            para_in=(
                (SP.PCOEFHR, LC.CHECHPR),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
                (SP.PT_EXTR, LC.ET_EXTR),
                (OP.CHAR_THER_ECHA_R.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_RAYO_F(
            te=461,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PRAYONF, LC.CRAYONF),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
                (OP.CHAR_THER_RAYO_F.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_RAYO_R(
            te=461,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PRAYONR, LC.CRAYONR),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
                (OP.CHAR_THER_RAYO_R.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.RIGI_THER_ECHA_F(
            te=457,
            para_in=(
                (SP.PCOEFHF, LC.CHECHPF),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PINSTR, LC.CTIMETR),
                (OP.RIGI_THER_ECHA_F.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PMATTTR, MMATTTR),),
        ),
        OP.RIGI_THER_ECHA_R(
            te=457,
            para_in=(
                (SP.PCOEFHR, LC.CHECHPR),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PINSTR, LC.CTIMETR),
                (OP.RIGI_THER_ECHA_R.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PMATTTR, MMATTTR),),
        ),
        OP.MTAN_THER_FLUXNL(
            te=457,
            para_in=(
                (SP.PFLUXNL, LC.CFLUXNF),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PTEMPEI, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
                (OP.MTAN_THER_FLUXNL.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PMATTTR, MMATTTR),),
        ),
        OP.MTAN_THER_RAYO_F(
            te=457,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PRAYONF, LC.CRAYONF),
                (SP.PTEMPEI, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
                (OP.MTAN_THER_RAYO_F.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PMATTTR, MMATTTR),),
        ),
        OP.MTAN_THER_RAYO_R(
            te=457,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PRAYONR, LC.CRAYONR),
                (SP.PTEMPEI, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
                (OP.MTAN_THER_RAYO_R.PCHHOBS, CHHOBS),
            ),
            para_out=((OP.MTAN_THER_RAYO_R.PMATTTR, MMATTTR),),
        ),
        OP.TOU_INI_ELGA(
            te=99,
            para_out=(
                (OP.TOU_INI_ELGA.PFLUX_R, LC.EFLUX2R),
                (OP.TOU_INI_ELGA.PGEOM_R, LC.EGGEO2D),
                (OP.TOU_INI_ELGA.PNEUT_F, LC.EGTINIF),
                (OP.TOU_INI_ELGA.PNEUT_R, LC.EGTINIR),
                (SP.PTEMP_R, LC.ETEMPPG),
            ),
        ),
        OP.TOU_INI_ELEM(
            te=99,
            para_out=((OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM2D), (OP.TOU_INI_ELEM.PFLUN_R, LC.CFLUXNR)),
        ),
        OP.TOU_INI_ELNO(
            te=99,
            para_out=(
                (OP.TOU_INI_ELNO.PFLUX_R, LC.NFLUX2R),
                (OP.TOU_INI_ELNO.PGEOM_R, LC.EGEOM2D),
                (OP.TOU_INI_ELNO.PHYDR_R, LC.EHYDRNO),
                (OP.TOU_INI_ELNO.PINST_R, LC.ENINST_R),
                (OP.TOU_INI_ELNO.PNEUT_F, LC.ENNEUT_F),
                (OP.TOU_INI_ELNO.PNEUT_R, LC.ENNEUT_R),
                (OP.TOU_INI_ELNO.PVARI_R, LC.EPHASES),
            ),
        ),
    )


# --------------------------------------------------------------------------------------------------
class THER_AX_HHO3_F(THER_2D_HHO3_F):
    """Thermics - Skin element HHO_QUAD - AXIS - SEG"""
