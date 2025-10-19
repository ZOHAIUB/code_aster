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

# person_in_charge: albert.alarcon at edf.fr

# Fichier comprenant une procédure de test de CALC_ESSAI avec les différentes
# options de calcul.

from .ce_calcul_modifstruct import CalcEssaiModifStruct
from ...Cata.Syntax import _F

from .ce_calcul_expansion import CalcEssaiExpansion
from .ce_calcul_identification import CalcEssaiIdentification


def TestCalcEssai(
    macro,
    mess,
    out_identification,
    out_modifstru,
    objects,
    EXPANSION,
    IDENTIFICATION,
    MODIFSTRUCT,
    GROUP_NO_CAPTEURS,
    GROUP_NO_EXTERIEUR,
):

    #
    # Transformation des objets Aster en classes python #
    # Test de ces classes                  #
    #

    # classe CalcEssaiObject
    # La classe MO recupere tous les concepts aster disponibles et
    # les range dans des classes python telles que Resultat, InterSpectre
    # Modele, CaraElem...

    #
    # 1) Test des fonctions de correlation #
    #
    if EXPANSION:
        assert len(EXPANSION) == 1, "'**' but not supported here!"
        EXPANSION = EXPANSION[0]
        calc_essai = CalcEssaiExpansion(macro, mess, objects)

        # Transformation des objets Aster en classes Python
        resu_num = objects.resultats[EXPANSION["CALCUL"].getName()]
        resu_exp = objects.resultats[EXPANSION["MESURE"].getName()]

        # Set default indices for the modes
        if EXPANSION["NUME_MODE_CALCUL"] != 0:
            nume_mode_calcul = list(EXPANSION["NUME_MODE_CALCUL"])
        else:
            nume_mode_calcul = None
        if EXPANSION["NUME_MODE_MESURE"] != 0:
            nume_mode_mesure = list(EXPANSION["NUME_MODE_MESURE"])
        else:
            nume_mode_mesure = None

        parametres = {}
        if EXPANSION["RESOLUTION"] == "LU":
            parametres["METHODE"] = "LU"
        elif EXPANSION["RESOLUTION"] == "SVD":
            parametres["METHODE"] = "SVD"
            parametres["EPS"] = EXPANSION["EPS"]

        calc_essai.setup(resu_num, nume_mode_calcul, resu_exp, nume_mode_mesure, parametres)
        # Les concepts resultats sont srockes sous le nom RESU_+type de resu
        result = calc_essai.calc_proj_resu()
        calc_essai.calc_mac_mode(result.ET, result.NX, resu_num.mass)

    #
    # 4) Test des fonctions de l'onglet identification #
    #
    if IDENTIFICATION:
        assert len(IDENTIFICATION) == 1, "'**' but supported here!"
        IDENTIFICATION = IDENTIFICATION[0]
        calcturb = CalcEssaiIdentification(objects, mess)

        inter_spec_name = IDENTIFICATION["INTE_SPEC"].getName()
        obs_name = IDENTIFICATION["OBSERVABILITE"].getName()
        com_name = IDENTIFICATION["COMMANDABILITE"].getName()
        res_base = IDENTIFICATION["BASE"].getName()

        calcturb.set_interspectre(objects.get_inter_spec(inter_spec_name))
        calcturb.set_observabilite(objects.get_mode_meca(obs_name))
        calcturb.set_commandabilite(objects.get_mode_meca(com_name))
        calcturb.set_base(objects.get_mode_meca(res_base))

        calcturb.set_alpha(IDENTIFICATION["ALPHA"])
        calcturb.set_epsilon(IDENTIFICATION["EPS"])
        calcturb.set_mcoeff(0.0)

        calcturb.calculate_force()

        calcturb.Sff.make_inte_spec(
            titre="Resultat identification : efforts", paras_out=out_identification
        )
        calcturb.Syy.make_inte_spec(
            titre="Resultat identification : mesure", paras_out=out_identification
        )
        calcturb.Syy_S.make_inte_spec(
            titre="Resultat identification : synthese", paras_out=out_identification
        )

    #
    # 5) Test des fonctions de l'onglet modif struct #
    #
    if MODIFSTRUCT:
        lance_modif_struct_calcul(
            macro, objects, MODIFSTRUCT, GROUP_NO_CAPTEURS, GROUP_NO_EXTERIEUR, out_modifstru
        )


class MessageBox:
    """!Classe qui permet d'ecrire dans un .mess separe"""

    def __init__(self, unite):
        self.unite = unite  # unite d'ecriture
        self.mess_file = open("fort." + str(unite), "w")

    def disp_mess(self, new_mess):
        """!Ecriture des messages dans le fichier sortie
        s'il existe et dans la fenetre de message"""
        self.mess_file.writelines(new_mess + "\n")

    def close_file(self):
        """Ferme le fichier de message a la fin de l'utilisation"""
        self.mess_file.close()


def to_dict_lst(groups):
    """Transforme la liste renvoyée par le superviseur
    en une liste de dictionaires. Il est important que
    les degrés de liberté soient dans une liste."""
    res_lst = []
    for grp in groups:
        rdict = {"NOM": grp["GROUP_NO"], "NOM_CMP": list(grp["NOM_CMP"])}
        res_lst.append(rdict)
    return res_lst


def lance_modif_struct_calcul(
    macro, ce_objects, MODIFSTRUCT, GROUP_NO_CAPTEURS, GROUP_NO_EXTERIEUR, out_modifstru
):
    """Démarre le calcul CALC_ESSAI sur la structure modifiée.

    :param macro: la macro étape CALC_ESSAI.

    :param ce_objects: les objects, utilisés par le module CALC_ESSAI,
                           présents dans la mémoire JEVEUX au moment
                           de la macro étape.

    :param MODIFSTRUCT: la macro étape permettant le calcul de la structure
                        modifiée depuis le fichier de commande Aster.

    :param CAPTEURS: dictionaire (ou FACT) contenant les choix
                     de groupes de noeuds et leurs degrés de liberté
                     pour les capteurs.

    :param EXTERIEUR: dictionaire (ou FACT) contenant les choix
                       de groupes de noeuds et leurs degrés de liberté
                       pour les ddl exterieurs.

    :param out_modifstru: dictionaire (ou FACT) utilisé pour les résultats."""

    modif_struct = CalcEssaiModifStruct(macro, ce_objects, ce_objects.mess, out_modifstru)

    modif_struct.find_experimental_result_from(MODIFSTRUCT["MESURE"].getName())
    modif_struct.find_support_modele_from(MODIFSTRUCT["MODELE_SUP"].getName())
    modif_struct.set_stiffness_matrix(MODIFSTRUCT["MATR_RIGI"])
    modif_struct.set_method_name(MODIFSTRUCT["RESOLUTION"])

    modif_struct.set_sensor_groups(to_dict_lst(GROUP_NO_CAPTEURS))
    modif_struct.set_interface_groups(to_dict_lst(GROUP_NO_EXTERIEUR))

    modif_struct.set_modes_ide(MODIFSTRUCT["NUME_MODE_MESU"])
    modif_struct.set_modes_expansion(MODIFSTRUCT["NUME_MODE_CALCUL"])

    modif_struct.set_sumail_name("SUMAIL")  # nom super-maille par defaut
    modif_struct.find_modele_modif_from(MODIFSTRUCT["MODELE_MODIF"].getName())
    modif_struct.find_maillage_modif_from(MODIFSTRUCT["MODELE_MODIF"].getName())

    modif_struct.calc_base_proj()

    modif_struct.set_param_condens({"METHODE": "SVD", "EPS": 1.0e-5, "REGUL": "NON"})
    modif_struct.creation_modele_couple()

    mode_simult = 1  # methode utilisee : MODE_ITER_SIMULT
    calc_freq = _F(OPTION="PLUS_PETITE", NMAX_FREQ=20, SEUIL_FREQ=1.0e-4)
    calc_freq["SEUIL_FREQ"] = 1e-4

    modif_struct.calc_modes_modele_couple(mode_simult, calc_freq)
