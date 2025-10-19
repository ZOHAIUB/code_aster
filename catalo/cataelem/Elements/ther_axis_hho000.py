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

# --------------------------------------------------------------------------------------------------
# Located components
# --------------------------------------------------------------------------------------------------
DDL_THER = LocatedComponents(
    phys=PHY.TEMP_R,
    type="ELNO",
    diff=True,
    components=(("EN1", ("HHO_FT[1]",)), ("EN2", ()), ("EN3", ("HHO_CT[1]"))),
)

TEMPHHO = LocatedComponents(phys=PHY.TEMP_R, type="ELNO", components=("TEMP",))

CHHOGT = LocatedComponents(phys=PHY.N1920R, type="ELEM", components=("X[10]",))

CHHOST = LocatedComponents(phys=PHY.N1360R, type="ELEM", components=("X[25]",))

CHHOBS = LocatedComponents(
    phys=PHY.N3600R,
    type="ELNO",
    diff=True,
    components=(("EN1", ("X[1]",)), ("EN2", ()), ("EN3", ("X[6]"))),
)

PFONC = LocatedComponents(phys=PHY.NEUT_K8, type="ELEM", components=("Z[2]",))

MVECTTR = ArrayOfComponents(phys=PHY.VTEM_R, locatedComponents=DDL_THER)

MMATTTR = ArrayOfComponents(phys=PHY.MTEM_R, locatedComponents=DDL_THER)

MMATTSR = ArrayOfComponents(phys=PHY.MTNS_R, locatedComponents=DDL_THER)


# --------------------------------------------------------------------------------------------------
class THERAXQ9_HHO000(Element):
    """Thermics - HHO_CSTE - AXIS - QUAD"""

    meshType = MT.QUAD9
    nodes = (
        SetOfNodes("EN1", (5, 6, 7, 8)),
        SetOfNodes("EN2", (1, 2, 3, 4)),
        SetOfNodes("EN3", (9,)),
    )
    elrefe = (
        ElrefeLoc(
            MT.QU9, gauss=("RIGI=FPG1", "FPG1=FPG1", "MASS=FPG1"), mater=("RIGI", "FPG1", "MASS")
        ),
    )
    calculs = (
        OP.CARA_CISA(te=-1),
        OP.CARA_GAUCHI(te=-1),
        OP.CARA_TORSION(te=-1),
        OP.CHAR_THER_EVOL(
            te=445,
            para_in=(
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
                (OP.CHAR_THER_EVOL.PVARCPR, LC.ZVARCPG),
                (OP.CHAR_THER_EVOL.PCHHOGT, CHHOGT),
                (OP.CHAR_THER_EVOL.PCHHOST, CHHOST),
                (OP.CHAR_THER_EVOL.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_SOURNL(
            te=465,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PSOURNL, LC.CSOURCF),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
                (OP.CHAR_THER_SOURNL.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_SOUR_F(
            te=465,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PSOURCF, LC.CSOURCF),
                (SP.PINSTR, LC.CTIMETR),
                (OP.CHAR_THER_SOUR_F.PVARCPR, LC.ZVARCPG),
                (OP.CHAR_THER_SOUR_F.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_SOUR_R(
            te=465,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PSOURCR, LC.ESOURCR),
                (SP.PINSTR, LC.CTIMETR),
                (OP.CHAR_THER_SOUR_R.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.COOR_ELGA(
            te=488,
            para_in=((SP.PGEOMER, LC.EGEOM2D),),
            para_out=((OP.COOR_ELGA.PCOORPG, LC.EGGAU2D),),
        ),
        OP.FLUX_ELGA(
            te=487,
            para_in=(
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
                (OP.FLUX_ELGA.PVARCPR, LC.ZVARCPG),
                (OP.FLUX_ELGA.PCHHOGT, CHHOGT),
                (OP.FLUX_ELGA.PCHHOST, CHHOST),
            ),
            para_out=((OP.FLUX_ELGA.PFLUXPG, LC.EFLUX2R),),
        ),
        OP.FLUX_ELNO(
            te=4,
            para_in=((OP.FLUX_ELNO.PFLUXPG, LC.EFLUX2R),),
            para_out=((SP.PFLUXNO, LC.NFLUX2R),),
        ),
        OP.HHO_PRECALC_BS(
            te=494,
            para_in=((SP.PGEOMER, LC.EGEOM2D),),
            para_out=((OP.HHO_PRECALC_BS.PCHHOBO, CHHOBS),),
        ),
        OP.HHO_PRECALC_OP(
            te=460,
            para_in=((SP.PGEOMER, LC.EGEOM2D), (OP.HHO_PRECALC_OP.PCHHOBS, CHHOBS)),
            para_out=((OP.HHO_PRECALC_OP.PCHHOGT, CHHOGT), (OP.HHO_PRECALC_OP.PCHHOST, CHHOST)),
        ),
        OP.HHO_PROJ_THER(
            te=473,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (OP.HHO_PROJ_THER.PFUNC_R, PFONC),
                (SP.PINSTPR, LC.MTEMPSR),
                (OP.HHO_PROJ_THER.PCHHOBS, CHHOBS),
            ),
            para_out=((OP.HHO_PROJ_THER.PTEMP_R, DDL_THER),),
        ),
        OP.HHO_PROJ2_THER(
            te=484,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (OP.HHO_PROJ2_THER.PH1TP_R, TEMPHHO),
                (OP.HHO_PROJ2_THER.PCHHOBS, CHHOBS),
            ),
            para_out=((OP.HHO_PROJ_THER.PTEMP_R, DDL_THER),),
        ),
        OP.HHO_PROJ3_THER(
            te=484,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (OP.HHO_PROJ3_THER.PQPTP_R, LC.ETEMPPG),
                (OP.HHO_PROJ3_THER.PCHHOBS, CHHOBS),
            ),
            para_out=((OP.HHO_PROJ3_THER.PTEMP_R, DDL_THER),),
        ),
        OP.HHO_TEMP_THER(
            te=456,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PTMPCHF, DDL_THER),
                (OP.HHO_TEMP_THER.PCHHOBS, CHHOBS),
            ),
            para_out=((OP.HHO_TEMP_THER.PTEMP_R, TEMPHHO),),
        ),
        OP.HHO_CINE_R_THER(
            te=492,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (OP.HHO_CINE_R_THER.PCMPVALE, TEMPHHO),
                (OP.HHO_CINE_R_THER.PCHHOBS, CHHOBS),
            ),
            para_out=((OP.HHO_CINE_R_THER.PCINE, DDL_THER),),
        ),
        OP.INIT_VARC(te=99, para_out=((OP.INIT_VARC.PVARCPR, LC.ZVARCPG),)),
        OP.MASS_THER(
            te=449,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.CTIMETR),
                (OP.MASS_THER.PVARCPR, LC.ZVARCPG),
                (OP.MASS_THER.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PMATTTR, MMATTTR),),
        ),
        OP.MASS_THER_TANG(
            te=429,
            para_in=(
                (OP.MASS_THER_TANG.PCOMPOR, LC.CCOMPOT),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTEMPEI, DDL_THER),
                (OP.MASS_THER_TANG.PVARCPR, LC.ZVARCPG),
                (OP.MASS_THER_TANG.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PMATTTR, MMATTTR),),
        ),
        OP.MASS_THER_RESI(
            te=429,
            para_in=(
                (OP.MASS_THER_RESI.PCOMPOR, LC.CCOMPOT),
                (SP.PGEOMER, LC.EGEOM2D),
                (OP.MASS_THER_RESI.PHYDRPR, LC.EHYDRR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTEMPEI, DDL_THER),
                (OP.MASS_THER_RESI.PVARCPR, LC.ZVARCPG),
                (OP.MASS_THER_RESI.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PRESIDU, MVECTTR),),
        ),
        OP.MTAN_THER_SOURNL(
            te=437,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PSOURNL, LC.CSOURCF),
                (SP.PTEMPEI, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
                (OP.MTAN_THER_SOURNL.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PMATTTR, MMATTTR),),
        ),
        OP.NSPG_NBVA(
            te=496,
            para_in=((OP.NSPG_NBVA.PCOMPOR, LC.CCOMPO2),),
            para_out=((SP.PDCEL_I, LC.EDCEL_I),),
        ),
        OP.RAPH_THER(
            te=477,
            para_in=(
                (OP.RAPH_THER.PCOMPOR, LC.CCOMPOT),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTEMPEI, DDL_THER),
                (OP.RAPH_THER.PVARCPR, LC.ZVARCPG),
                (OP.RAPH_THER.PCHHOGT, CHHOGT),
                (OP.RAPH_THER.PCHHOST, CHHOST),
                (OP.RAPH_THER.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PRESIDU, MVECTTR), (OP.RAPH_THER.PFLUXPR, LC.EFLUX2R)),
        ),
        OP.RIGI_THER(
            te=454,
            para_in=(
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.CTIMETR),
                (OP.RIGI_THER.PVARCPR, LC.ZVARCPG),
                (OP.RIGI_THER.PCHHOGT, CHHOGT),
                (OP.RIGI_THER.PCHHOST, CHHOST),
                (OP.RIGI_THER.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PMATTTR, MMATTTR),),
        ),
        OP.RIGI_THER_TANG(
            te=454,
            para_in=(
                (OP.RIGI_THER_TANG.PCOMPOR, LC.CCOMPOT),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTEMPEI, DDL_THER),
                (OP.RIGI_THER_TANG.PVARCPR, LC.ZVARCPG),
                (OP.RIGI_THER_TANG.PCHHOGT, CHHOGT),
                (OP.RIGI_THER_TANG.PCHHOST, CHHOST),
                (OP.RIGI_THER_TANG.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PMATTTR, MMATTTR), (OP.RIGI_THER_TANG.PMATTSR, MMATTSR)),
        ),
        OP.TEMP_ELGA(
            te=456,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PTEMPER, DDL_THER),
                (OP.TEMP_ELGA.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PTEMP_R, LC.ETEMPPG),),
        ),
        OP.TOU_INI_ELEM(te=99, para_out=((OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM3D),)),
        OP.TOU_INI_ELEM(
            te=99,
            para_out=((OP.TOU_INI_ELEM.PCOEH_R, LC.EHECHPR), (OP.TOU_INI_ELEM.PSOUR_R, LC.CSOURCR)),
        ),
        OP.TOU_INI_ELGA(
            te=99,
            para_out=(
                (OP.TOU_INI_ELGA.PFLUX_R, LC.EFLUX2R),
                (OP.TOU_INI_ELGA.PGEOM_R, LC.EGGEO2D),
                (OP.TOU_INI_ELGA.PNEUT_F, LC.EGTINIF),
                (OP.TOU_INI_ELGA.PNEUT_R, LC.EGTINIR),
                (OP.TOU_INI_ELGA.PSOUR_R, LC.ESOURCR),
                (OP.TOU_INI_ELGA.PVARI_R, LC.ZVARIPG),
                (SP.PTEMP_R, LC.ETEMPPG),
            ),
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
                (OP.TOU_INI_ELNO.PSOUR_R, LC.NSOURCR),
            ),
        ),
        OP.VERI_JACOBIEN(
            te=328, para_in=((SP.PGEOMER, LC.EGEOM2D),), para_out=((SP.PCODRET, LC.ECODRET),)
        ),
    )


# --------------------------------------------------------------------------------------------------
class THERAXT7_HHO000(THERAXQ9_HHO000):
    """Thermics - HHO_CSTE - AXIS - TRIA"""

    meshType = MT.TRIA7
    nodes = (SetOfNodes("EN1", (4, 5, 6)), SetOfNodes("EN2", (1, 2, 3)), SetOfNodes("EN3", (7,)))
    elrefe = (
        ElrefeLoc(
            MT.TR7, gauss=("RIGI=FPG1", "FPG1=FPG1", "MASS=FPG1"), mater=("RIGI", "FPG1", "MASS")
        ),
    )
