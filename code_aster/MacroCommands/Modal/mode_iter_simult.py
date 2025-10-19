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

#

from ...Cata.Commons import *
from ...Cata.DataStructure import *
from ...Cata.Syntax import *
from ...Supervis import ExecuteCommand
from ...Objects import (
    AssemblyMatrixDisplacementComplex,
    AssemblyMatrixDisplacementReal,
    AssemblyMatrixEliminatedReal,
    AssemblyMatrixPressureReal,
    AssemblyMatrixTemperatureReal,
    BucklingModeResult,
    GeneralizedAssemblyMatrixComplex,
    GeneralizedAssemblyMatrixReal,
    GeneralizedModeResult,
    ModeResultComplex,
    ModeResult,
    AcousticModeResult,
)


def mode_iter_simult_prod(TYPE_RESU, **args):
    if args.get("__all__"):
        return (mode_flamb, mode_meca_c, mode_meca, mode_acou, mode_gene, ASSD)

    if TYPE_RESU not in ["DYNAMIQUE", "MODE_FLAMB", "GENERAL"]:
        # on retourne un type fictif pour que le plantage aie lieu dans la lecture du catalogue
        return ASSD
    if TYPE_RESU == "MODE_FLAMB":
        return mode_flamb
    if TYPE_RESU == "GENERAL":
        return mode_flamb
    # sinon on est dans le cas 'DYNAMIQUE' donc **args doit contenir les mots-clés
    # MATR_RIGI et (faculativement) MATR_AMOR, et on peut y accéder
    vale_rigi = args["MATR_RIGI"]
    if vale_rigi is None:  # si MATR_RIGI non renseigné
        # on retourne un type fictif pour que le plantage aie lieu dans la lecture du catalogue
        return ASSD
    vale_amor = args.get("MATR_AMOR")
    if AsType(vale_amor) == matr_asse_depl_r:
        return mode_meca_c
    if AsType(vale_rigi) == matr_asse_depl_r:
        return mode_meca
    if AsType(vale_rigi) == matr_asse_elim_r:
        return mode_meca
    if AsType(vale_rigi) == matr_asse_temp_r:
        return mode_meca
    if AsType(vale_rigi) == matr_asse_depl_c:
        return mode_meca_c
    if AsType(vale_rigi) == matr_asse_pres_r:
        return mode_acou
    if AsType(vale_rigi) == matr_asse_gene_r:
        return mode_gene
    if AsType(vale_rigi) == matr_asse_gene_c:
        return mode_gene

    raise CataError("type de concept resultat non prevu")


MODE_ITER_SIMULT_CATA = OPER(
    nom="MODE_ITER_SIMULT",
    op=45,
    sd_prod=mode_iter_simult_prod,
    fr=tr(
        "Calcul des modes propres par itérations simultanées : valeurs propres et modes propres réels ou complexes"
    ),
    reentrant="n",
    METHODE=SIMP(
        statut="f", typ="TXM", defaut="SORENSEN", into=("TRI_DIAG", "JACOBI", "SORENSEN", "QZ")
    ),
    b_tri_diag=BLOC(
        condition="METHODE == 'TRI_DIAG'",
        PREC_ORTHO=SIMP(statut="f", typ="R", defaut=1.0e-12, val_min=0.0e0),
        NMAX_ITER_ORTHO=SIMP(statut="f", typ="I", defaut=5, val_min=0),
        PREC_LANCZOS=SIMP(statut="f", typ="R", defaut=1.0e-8, val_min=0.0e0),
        NMAX_ITER_QR=SIMP(statut="f", typ="I", defaut=30, val_min=0),
    ),
    b_jacobi=BLOC(
        condition="METHODE == 'JACOBI'",
        PREC_BATHE=SIMP(statut="f", typ="R", defaut=1.0e-10, val_min=0.0e0),
        NMAX_ITER_BATHE=SIMP(statut="f", typ="I", defaut=40, val_min=0),
        PREC_JACOBI=SIMP(statut="f", typ="R", defaut=1.0e-2, val_min=0.0e0),
        NMAX_ITER_JACOBI=SIMP(statut="f", typ="I", defaut=12, val_min=0),
    ),
    b_sorensen=BLOC(
        condition="METHODE == 'SORENSEN'",
        PREC_SOREN=SIMP(statut="f", typ="R", defaut=0.0e0, val_min=0.0e0),
        NMAX_ITER_SOREN=SIMP(statut="f", typ="I", defaut=20, val_min=1),
        PARA_ORTHO_SOREN=SIMP(statut="f", typ="R", defaut=0.717),
    ),
    b_qz=BLOC(
        condition="METHODE == 'QZ'",
        TYPE_QZ=SIMP(
            statut="f", typ="TXM", defaut="QZ_SIMPLE", into=("QZ_QR", "QZ_SIMPLE", "QZ_EQUI")
        ),
    ),
    CHAM_MATER=SIMP(statut="f", typ=cham_mater),
    CARA_ELEM=SIMP(statut="f", typ=cara_elem),
    TYPE_RESU=SIMP(
        statut="f",
        typ="TXM",
        defaut="DYNAMIQUE",
        into=("DYNAMIQUE", "MODE_FLAMB", "GENERAL"),
        fr=tr("Type d analyse"),
    ),
    OPTION=SIMP(
        statut="f",
        typ="TXM",
        defaut="SANS",
        into=("MODE_RIGIDE", "SANS"),
        fr=tr("Calcul des modes de corps rigide, uniquement pour la méthode TRI_DIAG"),
    ),
    b_dynam=BLOC(
        condition="TYPE_RESU == 'DYNAMIQUE'",
        MATR_RIGI=SIMP(
            statut="o",
            typ=(
                matr_asse_depl_r,
                matr_asse_depl_c,
                matr_asse_temp_r,
                matr_asse_gene_r,
                matr_asse_gene_c,
                matr_asse_pres_r,
                matr_asse_elim_r,
            ),
        ),
        MATR_MASS=SIMP(
            statut="o",
            typ=(
                matr_asse_depl_r,
                matr_asse_gene_r,
                matr_asse_pres_r,
                matr_asse_temp_r,
                matr_asse_elim_r,
            ),
        ),
        MATR_AMOR=SIMP(statut="f", typ=(matr_asse_depl_r, matr_asse_gene_r)),
        CALC_FREQ=FACT(
            statut="d",
            min=0,
            OPTION=SIMP(
                statut="f",
                typ="TXM",
                defaut="PLUS_PETITE",
                into=("PLUS_PETITE", "PLUS_GRANDE", "BANDE", "CENTRE", "TOUT"),
                fr=tr("Choix de l option et par conséquent du shift du problème modal"),
            ),
            b_plus_petite=BLOC(
                condition="OPTION == 'PLUS_PETITE'",
                fr=tr("Recherche des plus petites fréquences propres"),
                NMAX_FREQ=SIMP(statut="f", typ="I", defaut=10, val_min=0),
            ),
            b_plus_grande=BLOC(
                condition="OPTION == 'PLUS_GRANDE'",
                fr=tr("Recherche des plus grandes fréquences propres"),
                NMAX_FREQ=SIMP(statut="f", typ="I", defaut=1, val_min=0),
            ),
            b_centre=BLOC(
                condition="OPTION == 'CENTRE'",
                fr=tr("Recherche des fréquences propres les plus proches d'une valeur donnée"),
                FREQ=SIMP(
                    statut="o",
                    typ="R",
                    fr=tr("Fréquence autour de laquelle on cherche les fréquences propres"),
                ),
                AMOR_REDUIT=SIMP(statut="f", typ="R"),
                NMAX_FREQ=SIMP(statut="f", typ="I", defaut=10, val_min=0),
            ),
            b_bande=BLOC(
                condition="(OPTION == 'BANDE')",
                fr=tr("Recherche des fréquences propres dans une bande donnée"),
                FREQ=SIMP(
                    statut="o",
                    typ="R",
                    min=2,
                    max=2,
                    validators=AndVal((OrdList("croissant"), NoRepeat())),
                    fr=tr("Valeur des deux fréquences délimitant la bande de recherche"),
                ),
                TABLE_FREQ=SIMP(statut="f", typ=table_sdaster),
            ),
            APPROCHE=SIMP(
                statut="f",
                typ="TXM",
                defaut="REEL",
                into=("REEL", "IMAG", "COMPLEXE"),
                fr=tr(
                    "Choix du pseudo-produit scalaire pour la résolution du problème quadratique"
                ),
            ),
            regles=(EXCLUS("DIM_SOUS_ESPACE", "COEF_DIM_ESPACE"),),
            DIM_SOUS_ESPACE=SIMP(statut="f", typ="I"),
            COEF_DIM_ESPACE=SIMP(statut="f", typ="I"),
            NMAX_ITER_SHIFT=SIMP(statut="f", typ="I", defaut=3, val_min=0),
            PREC_SHIFT=SIMP(statut="f", typ="R", defaut=5.0e-2, val_min=0.0e0),
            SEUIL_FREQ=SIMP(statut="f", typ="R", defaut=1.0e-2, val_min=0.0e0),
        ),
    ),
    b_general=BLOC(
        condition="TYPE_RESU == 'GENERAL'",
        MATR_A=SIMP(statut="o", typ=(matr_asse_depl_r, matr_asse_gene_r, matr_asse_pres_r)),
        MATR_B=SIMP(statut="o", typ=(matr_asse_depl_r, matr_asse_gene_r, matr_asse_pres_r)),
    ),
    b_flamb=BLOC(
        condition="TYPE_RESU == 'MODE_FLAMB'",
        MATR_RIGI=SIMP(statut="o", typ=(matr_asse_depl_r, matr_asse_gene_r, matr_asse_pres_r)),
        MATR_RIGI_GEOM=SIMP(statut="o", typ=(matr_asse_depl_r, matr_asse_gene_r, matr_asse_pres_r)),
    ),
    b_flamb_general=BLOC(
        condition="(TYPE_RESU == 'MODE_FLAMB') or (TYPE_RESU == 'GENERAL')",
        CALC_CHAR_CRIT=FACT(
            statut="d",
            min=0,
            OPTION=SIMP(
                statut="f",
                typ="TXM",
                defaut="PLUS_PETITE",
                into=("PLUS_PETITE", "BANDE", "CENTRE", "TOUT"),
                fr=tr("Choix de l option et par conséquent du shift du problème modal"),
            ),
            b_plus_petite=BLOC(
                condition="OPTION == 'PLUS_PETITE'",
                fr=tr("Recherche des plus petites valeurs propres"),
                NMAX_CHAR_CRIT=SIMP(statut="f", typ="I", defaut=10, val_min=0),
            ),
            b_centre=BLOC(
                condition="OPTION == 'CENTRE'",
                fr=tr("Recherche des valeurs propres les plus proches d une valeur donnée"),
                CHAR_CRIT=SIMP(
                    statut="o",
                    typ="R",
                    fr=tr(
                        "Charge critique autour de laquelle on cherche les charges critiques propres"
                    ),
                ),
                NMAX_CHAR_CRIT=SIMP(statut="f", typ="I", defaut=10, val_min=0),
            ),
            b_bande=BLOC(
                condition="(OPTION == 'BANDE')",
                fr=tr("Recherche des valeurs propres dans une bande donnée"),
                CHAR_CRIT=SIMP(
                    statut="o",
                    typ="R",
                    min=2,
                    max=2,
                    validators=AndVal((OrdList("croissant"), NoRepeat())),
                    fr=tr("Valeur des deux charges critiques délimitant la bande de recherche"),
                ),
                TABLE_CHAR_CRIT=SIMP(statut="f", typ=table_sdaster),
            ),
            APPROCHE=SIMP(
                statut="f",
                typ="TXM",
                defaut="REEL",
                into=("REEL", "IMAG"),
                fr=tr(
                    "Choix du pseudo-produit scalaire pour la résolution du problème quadratique"
                ),
            ),
            regles=(EXCLUS("DIM_SOUS_ESPACE", "COEF_DIM_ESPACE"),),
            DIM_SOUS_ESPACE=SIMP(statut="f", typ="I"),
            COEF_DIM_ESPACE=SIMP(statut="f", typ="I"),
            NMAX_ITER_SHIFT=SIMP(statut="f", typ="I", defaut=3, val_min=0),
            PREC_SHIFT=SIMP(statut="f", typ="R", defaut=5.0e-2, val_min=0.0e0),
            SEUIL_CHAR_CRIT=SIMP(statut="f", typ="R", defaut=1.0e-2, val_min=0.0e0),
        ),
    ),
    # -------------------------------------------------------------------
    #        Catalogue commun SOLVEUR
    SOLVEUR=C_SOLVEUR("MODE_ITER_SIMULT"),
    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    #  Mot-cles caches pour activer le parallelisme au sein d'une macro-commande
    PARALLELISME_MACRO=FACT(
        statut="d",
        min=0,
        TYPE_COM=SIMP(
            statut="c", typ="I", defaut=-999, into=(-999, 1), fr=tr("Type de communication")
        ),
        IPARA1_COM=SIMP(
            statut="c", typ="I", defaut=-999, fr=tr("Parametre entier n 1 de la communication")
        ),
        IPARA2_COM=SIMP(
            statut="c", typ="I", defaut=-999, fr=tr("Parametre entier n 2 de la communication")
        ),
    ),
    # -------------------------------------------------------------------
    VERI_MODE=FACT(
        statut="d",
        min=0,
        STOP_ERREUR=SIMP(statut="f", typ="TXM", defaut="OUI", into=("OUI", "NON")),
        PREC_SHIFT=SIMP(statut="f", typ="R", defaut=5.0e-3, val_min=0.0e0),
        SEUIL=SIMP(
            statut="f",
            typ="R",
            defaut=1.0e-6,
            val_min=0.0e0,
            fr=tr("Valeur limite admise pour l'erreur a posteriori des modes"),
        ),
        STURM=SIMP(statut="f", typ="TXM", defaut="OUI", into=("OUI", "NON")),
    ),
    STOP_BANDE_VIDE=SIMP(statut="f", typ="TXM", defaut="OUI", into=("OUI", "NON")),
    INFO=SIMP(statut="f", typ="I", defaut=1, into=(1, 2)),
    TITRE=SIMP(statut="f", typ="TXM", max="**"),
)


class ModalCalculationSimult(ExecuteCommand):
    """Internal (non public) command to call the underlying operator."""

    command_name = "MODE_ITER_SIMULT"
    command_cata = MODE_ITER_SIMULT_CATA

    def create_result(self, keywords):
        """Initialize the result.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """
        TYPE_RESU = keywords.get("TYPE_RESU")
        if TYPE_RESU in ("MODE_FLAMB", "GENERAL"):
            self._result = BucklingModeResult()
            return

        vale_rigi = keywords.get("MATR_RIGI")
        vale_amor = keywords.get("MATR_AMOR")
        if vale_amor is not None and isinstance(vale_amor, AssemblyMatrixDisplacementReal):
            self._result = ModeResultComplex()
        elif isinstance(vale_rigi, AssemblyMatrixEliminatedReal):
            self._result = ModeResult()
        elif isinstance(vale_rigi, AssemblyMatrixDisplacementReal):
            self._result = ModeResult()
        elif isinstance(vale_rigi, AssemblyMatrixTemperatureReal):
            self._result = ModeResult()
        elif isinstance(vale_rigi, AssemblyMatrixDisplacementComplex):
            self._result = ModeResultComplex()
        elif isinstance(vale_rigi, AssemblyMatrixPressureReal):
            self._result = AcousticModeResult()
        elif isinstance(vale_rigi, GeneralizedAssemblyMatrixReal):
            self._result = GeneralizedModeResult()
        elif isinstance(vale_rigi, GeneralizedAssemblyMatrixComplex):
            self._result = GeneralizedModeResult()

    def post_exec(self, keywords):
        """Execute the command.

        Arguments:
            keywords (dict): User's keywords.
        """
        if self.exception:
            return
        matrRigi = keywords.get("MATR_RIGI")
        if matrRigi is not None:
            if isinstance(self._result, GeneralizedModeResult):
                nume_ddl_gene = matrRigi.getGeneralizedDOFNumbering()
                self._result.setGeneralizedDOFNumbering(nume_ddl_gene)
                basis = nume_ddl_gene.getModalBasis()
                if basis is not None:
                    mesh = basis.getMesh()
                    if mesh is not None:
                        self._result.setMesh(mesh)
            else:
                self._result.setDOFNumbering(matrRigi.getDOFNumbering())
            self._result.setStiffnessMatrix(matrRigi)
        matrAmor = keywords.get("MATR_AMOR")
        if matrAmor is not None:
            self._result.setDampingMatrix(matrAmor)

        matrA = keywords.get("MATR_A")
        if matrA is not None:
            if self._result.getMesh() is None:
                self._result.setMesh(matrA.getMesh())
            self._result.setDOFNumbering(matrA.getDOFNumbering())

        self._result.build()


MODE_ITER_SIMULT = ModalCalculationSimult.run
