# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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

from code_aster.Commands import *
from code_aster import CA

test = CA.TestCase()

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))


pMesh2 = CA.Mesh()
pMesh2.readMedFile("zzzz504i.med")
DEFI_GROUP(reuse=pMesh2, MAILLAGE=pMesh2, CREA_GROUP_NO=_F(TOUT_GROUP_MA="OUI"))

model = AFFE_MODELE(
    MAILLAGE=pMesh2,
    AFFE=_F(MODELISATION="3D", PHENOMENE="MECANIQUE", TOUT="OUI"),
    DISTRIBUTION=_F(METHODE="CENTRALISE"),
)

char_cin = AFFE_CHAR_CINE(MODELE=model, MECA_IMPO=_F(GROUP_NO="COTE_H", DX=0.0, DY=0.0, DZ=0.0))

char_meca = AFFE_CHAR_MECA(
    MODELE=model,
    LIAISON_DDL=_F(GROUP_NO=("A", "B"), DDL=("DX", "DX"), COEF_MULT=(1.0, -1.0), COEF_IMPO=0),
    DDL_IMPO=_F(GROUP_NO="A", DX=1.0),
)

MATER1 = DEFI_MATERIAU(ELAS=_F(E=200000.0, NU=0.3, RHO=1.0e-3))

char_vol = AFFE_CHAR_MECA(MODELE=model, PESANTEUR=_F(GRAVITE=9.81, DIRECTION=(0.0, 0.0, -1.0)))


AFFMAT = AFFE_MATERIAU(MAILLAGE=pMesh2, AFFE=_F(TOUT="OUI", MATER=MATER1))

LI = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=_F(JUSQU_A=1.0, NOMBRE=1))

resu = STAT_NON_LINE(
    CHAM_MATER=AFFMAT,
    METHODE="NEWTON",
    COMPORTEMENT=_F(RELATION="ELAS"),
    CONVERGENCE=_F(RESI_GLOB_RELA=1.0e-8),
    EXCIT=(_F(CHARGE=char_cin), _F(CHARGE=char_vol), _F(CHARGE=char_meca)),
    INCREMENT=_F(LIST_INST=LI),
    MODELE=model,
    NEWTON=_F(MATRICE="TANGENTE", REAC_ITER=1),
    SOLVEUR=_F(METHODE="PETSC", RESI_RELA=1.0e-5, PRE_COND="JACOBI"),
)


# if (parallel):
# rank = MPI.ASTER_COMM_WORLD.Get_rank()
# myFile='par.txt'
# os.system("sed 's/Mat_.*\=/par\ \=/g' /tmp/par.txt > /tmp/par_clean.txt && mv /tmp/par_clean.txt /tmp/par.txt")
# if (rank==0): os.system( """grep -v %% /tmp/%s | grep -v zzz | grep -v \] | grep -v Mat | awk '{print $3}' | LANG=en_US.UTF-8  sort -g > /tmp/%s_sorted"""%(myFile,myFile) )
# if (rank==0): os.system( """grep -v Object /tmp/rhs_par.txt | grep -v type | grep -v Process | LANG=en_US.UTF-8  sort -g > /tmp/rhs_par_sorted.txt  """)
# if (rank==0): os.system( """grep -v Object /tmp/sol_par.txt | grep -v type | grep -v Process | LANG=en_US.UTF-8  sort -g > /tmp/sol_par_sorted.txt  """)
# else:
# myFile='seq.txt'
# os.system("sed 's/Mat_.*\=/seq\ \=/g' /tmp/seq.txt > /tmp/seq_clean.txt && mv /tmp/seq_clean.txt /tmp/seq.txt")
# os.system("grep -v %% /tmp/%s | grep -v zzz | grep -v \] | grep -v Mat | awk '{print $3}' | LANG=en_US.UTF-8  sort -g > /tmp/%s_sorted"%(myFile,myFile))
# os.system( """grep -v Object /tmp/rhs_seq.txt | grep -v type | grep -v Process | LANG=en_US.UTF-8  sort -g > /tmp/rhs_seq_sorted.txt  """)
# os.system( """grep -v Object /tmp/sol_seq.txt | grep -v type | grep -v Process | LANG=en_US.UTF-8  sort -g > /tmp/sol_seq_sorted.txt  """)

# if (parallel):
# rank = MPI.ASTER_COMM_WORLD.Get_rank()
# if (rank==0): os.system( """cp fort.11 /tmp/ddl0.txt """ )
# if (rank==1): os.system( """cp fort.12 /tmp/ddl1.txt """ )
# if (rank==2): os.system( """cp fort.13 /tmp/ddl2.txt """ )
# if (rank==3): os.system( """cp fort.14 /tmp/ddl3.txt """ )
# else:
# os.system( """cp fort.11 /tmp/ddl_seq.txt """ )

# if parallel:
#     rank = MPI.ASTER_COMM_WORLD.Get_rank()
#     resu.printMedFile('/tmp/par_%d.resu.med'%rank)
# else:
#     resu.printMedFile('/tmp/seq.resu.med')

MyFieldOnNodes = resu.getField("DEPL", 1)
sfon = MyFieldOnNodes.toSimpleFieldOnNodes()

test.assertAlmostEqual(sfon[0, 0], 0.0)
test.assertAlmostEqual(sfon[240, 0], 0.02832458511648515)

test.printSummary()

FIN()
