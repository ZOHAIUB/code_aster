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

# person_in_charge: harinaivo.andriambololona at edf.fr

import copy

import numpy

from ..Cata.Syntax import _F
from ..CodeCommands import CREA_CHAMP, CREA_RESU, DEFI_GROUP
from ..CodeCommands import MODI_REPERE as MODI_REPERE_CMD
from ..CodeCommands import POST_RELEVE_T, PROJ_CHAMP, RECU_TABLE
from ..Messages import UTMESS


def observation_ops(
    self,
    PROJECTION=None,
    MODELE_1=None,
    MODELE_2=None,
    RESULTAT=None,
    MATR_RIGI=None,
    MATR_MASS=None,
    MODI_REPERE=None,
    NOM_CHAM=None,
    FILTRE=None,
    EPSI_MOYENNE=None,
    **args
):
    """
    Ecriture de la macro MACRO_OBSERVATION
    """
    # dans **args, on range les options de PROJ_CHAMP facultatives, et dont on
    # ne sert pas par la suite
    mcfact = args

    if not isinstance(NOM_CHAM, tuple) and not isinstance(NOM_CHAM, list):
        NOM_CHAM = [NOM_CHAM]

    TYPE_CHAM = None
    RESU = None

    TYPE_RESU = RESULTAT.getType()

    # recuperation du maillage associe au modele numerique
    mayanum = MODELE_1.getMesh()

    # vis_a_vis
    if "VIS_A_VIS" in mcfact.keys():
        for arg in mcfact["VIS_A_VIS"]:
            if arg["GROUP_MA_1"] != arg["GROUP_MA_2"]:
                UTMESS("F", "CALCESSAI0_11")

    # modele numerique 2D ou 3D
    typmod = mayanum.getDimension()

    # recuperation du maillage associe au modele experimental
    mayaexp = MODELE_2.getMesh()

    # cham_mater et cara_elem pour le resultat a projeter
    if RESULTAT.getNumberOfIndexes() > 0:
        cara_elems = []
        for j in RESULTAT.getIndexes():
            if RESULTAT.hasElementaryCharacteristics(j):
                cara_elems += [RESULTAT.getElementaryCharacteristics(j)]
        assert len(cara_elems) <= 1
        if len(cara_elems):
            cara_elem = cara_elems[0]
        else:
            cara_elem = None
    else:
        cara_elem = None

    if RESULTAT.getNumberOfIndexes() > 0:
        cham_maters = []
        for j in RESULTAT.getIndexes():
            if RESULTAT.hasMaterialField(j):
                cham_maters += [RESULTAT.getMaterialField(j)]

        cham_maters = list(set(cham_maters))
        assert len(cham_maters) <= 1
        if len(cham_maters):
            cham_mater = cham_maters[0]
        else:
            cham_mater = None
    else:
        cham_mater = None

    # afreq pour les frequences propres
    # if isinstance(RESULTAT, mode_meca):
    if TYPE_RESU == "MODE_MECA":
        # frequences propres
        __freq = RECU_TABLE(CO=RESULTAT, NOM_PARA="FREQ")
        afreq = __freq.EXTR_TABLE().Array("NUME_ORDRE", "FREQ")
        # noms des matrices
        if MATR_RIGI is not None or MATR_MASS is not None:
            # recherche du nume_ddl associe
            NUME_DDL = MATR_RIGI.getDOFNumbering()
        else:
            UTMESS("A", "CALCESSAI0_9")
            NUME_DDL = None

    else:
        afreq = None
        NUME_DDL = None

    indice = list(range(RESULTAT.getNumberOfIndexes()))

    # ***********************************************
    #  PHASE DE CALCUL DE LA DEFORMATION MOYENNE AUX NOEUDS
    #  CHAMP CALCULE SUR LE MODELE NUMERIQUE
    # ***********************************************

    resu_epsi = None
    if EPSI_MOYENNE is not None:
        for nomcham in NOM_CHAM:
            if nomcham == "EPSI_NOEU":
                if RESULTAT.getType() == "DYNA_HARMO":
                    TYPE_CHAM = "NOEU_EPSI_C"
                else:
                    TYPE_CHAM = "NOEU_EPSI_R"

        if TYPE_CHAM is None:
            UTMESS("F", "UTILITAI8_24", valk=["NOEU_EPSI", nomcham])
        else:
            num_ordr = RESULTAT.getIndexes()

            if RESULTAT.getType() == "EVOL_ELAS":
                list_inst = RESULTAT.LIST_VARI_ACCES()["INST"]
            if RESULTAT.getType() == "DYNA_TRANS":
                list_inst = RESULTAT.LIST_VARI_ACCES()["INST"]
            if RESULTAT.getType() == "DYNA_HARMO":
                list_freq = RESULTAT.LIST_VARI_ACCES()["FREQ"]

            liste = []

            # il faut calculer le champ complet
            if typmod == 2:
                nom_cmps = ["EPXX", "EPYY", "EPZZ", "EPXY"]
            else:
                nom_cmps = ["EPXX", "EPYY", "EPZZ", "EPXY", "EPXZ", "EPYZ"]

            argsi = {"ACTION": []}
            lnoeuds = {}
            nb_mcfact = 0
            seuil = []
            masque = []
            for epsi_moye in EPSI_MOYENNE:
                mcfactr = {}
                mcfacti = {}
                l_noeud = None
                val_masque = []
                seuil_lu = epsi_moye["SEUIL_VARI"]
                if type(seuil_lu) == tuple:
                    val_seuil = seuil_lu[0]
                else:
                    val_seuil = seuil_lu
                seuil.append(val_seuil)
                masque_lu = epsi_moye["MASQUE"]
                if type(masque_lu) != tuple:
                    val_masque.append(masque_lu)
                else:
                    val_masque = masque_lu
                masque.append(val_masque)
                for typ in ["NOEUD", "GROUP_NO", "MAILLE", "GROUP_MA"]:
                    if epsi_moye[typ] is not None:
                        l_noeud = find_no(mayanum, {typ: epsi_moye[typ]})
                        nb_mcfact = nb_mcfact + 1
                        for i in range(len(l_noeud)):
                            l_noeud[i] = l_noeud[i].strip()
                        lnoeuds[str(nb_mcfact)] = l_noeud

                if l_noeud is None:
                    UTMESS("F", "MODELISA3_13", valk=["EPSI_MOYENNE"])

                if TYPE_CHAM[-1:] == "C":
                    mcfactr = {
                        "NOM_CMP": nom_cmps,
                        "OPERATION": "EXTRACTION",
                        "INTITULE": str("R" + str(nb_mcfact)),
                        "FORMAT_C": "REEL",
                        "NOEUD": l_noeud,
                        "NOM_CHAM": "EPSI_NOEU",
                        "RESULTAT": RESULTAT,
                        "NUME_ORDRE": num_ordr,
                    }
                    argsi["ACTION"].append(mcfactr)
                    mcfacti = {
                        "NOM_CMP": nom_cmps,
                        "OPERATION": "EXTRACTION",
                        "INTITULE": str("I" + str(nb_mcfact)),
                        "FORMAT_C": "IMAG",
                        "NOEUD": l_noeud,
                        "NOM_CHAM": "EPSI_NOEU",
                        "RESULTAT": RESULTAT,
                        "NUME_ORDRE": num_ordr,
                    }
                    argsi["ACTION"].append(mcfacti)
                else:
                    mcfactr = {
                        "NOM_CMP": nom_cmps,
                        "OPERATION": "EXTRACTION",
                        "INTITULE": str(nb_mcfact),
                        "NOEUD": l_noeud,
                        "NOM_CHAM": "EPSI_NOEU",
                        "RESULTAT": RESULTAT,
                        "NUME_ORDRE": num_ordr,
                    }
                    argsi["ACTION"].append(mcfactr)

            _tepsi = POST_RELEVE_T(**argsi)

            table = _tepsi.EXTR_TABLE()

            mcfact2 = {}
            __chame = [None] * len(indice)
            for ind in indice:
                argsa = {"AFFE": []}
                for mcfacta in range(nb_mcfact):
                    l_noeud_mcfact = lnoeuds[str(mcfacta + 1)]
                    l_vmoye = []
                    l_cmp_vari = []
                    seuil_mc = seuil[mcfacta]
                    masque_mc = masque[mcfacta]
                    for cmp in nom_cmps:
                        lur = 0
                        lui = 0
                        l_valr = []
                        l_vali = []
                        l_valc = []
                        l_val = []
                        for row in table.rows:
                            if TYPE_CHAM[-1:] == "C":
                                if (
                                    row["INTITULE"].strip() == str("R" + str(mcfacta + 1))
                                    and row["NUME_ORDRE"] == num_ordr[ind]
                                ):
                                    l_valr.append(row[cmp])
                                    lur = 1
                                elif (
                                    row["INTITULE"].strip() == str("I" + str(mcfacta + 1))
                                    and row["NUME_ORDRE"] == num_ordr[ind]
                                ):
                                    l_vali.append(row[cmp])
                                    lui = 1

                            else:
                                if (
                                    row["INTITULE"].strip() == str(mcfacta + 1)
                                    and row["NUME_ORDRE"] == num_ordr[ind]
                                ):
                                    l_val.append(row[cmp])

                        if TYPE_CHAM[-1:] == "C":
                            if lur and lui:
                                if len(l_valr) != len(l_vali):
                                    UTMESS("F", "POSTRELE_59")
                                for i in range(len(l_valr)):
                                    l_valc.append(complex(l_valr[i], l_vali[i]))

                                lur = 0
                                lui = 0
                            else:
                                UTMESS("F", "POSTRELE_59")

                            # on regarde a la fois la partie reelle et la partie imag
                            # pour les complexes
                            vmoyer = sum(l_valr) / len(l_valr)
                            vmoyei = sum(l_vali) / len(l_vali)
                            vmoye = sum(l_valc) / len(l_valc)
                            vminr = min(l_valr)
                            vmini = min(l_vali)
                            vmaxr = max(l_valr)
                            vmaxi = max(l_vali)
                            if vmoyer > 0:
                                if (vmaxr > vmoyer * (1.0 + seuil_mc)) or (
                                    vminr < vmoyer * (1 - seuil_mc)
                                ):
                                    l_cmp_vari.append(cmp)
                            if vmoyei > 0:
                                if (vmaxi > vmoyei * (1.0 + seuil_mc)) or (
                                    vmini < vmoyei * (1 - seuil_mc)
                                ):
                                    l_cmp_vari.append(cmp)
                            if vmoyer < 0:
                                if (vminr > vmoyer * (1.0 - seuil_mc)) or (
                                    vmaxr < vmoyer * (1 + seuil_mc)
                                ):
                                    l_cmp_vari.append(cmp)
                            if vmoyei < 0:
                                if (vmini > vmoyei * (1.0 - seuil_mc)) or (
                                    vmaxi < vmoyei * (1 + seuil_mc)
                                ):
                                    l_cmp_vari.append(cmp)
                        else:
                            vmoye = sum(l_val) / len(l_val)
                            vmin = min(l_val)
                            vmax = max(l_val)
                            if vmoye > 0:
                                if (vmax > vmoye * (1.0 + seuil_mc)) or (
                                    vmin < vmoye * (1 - seuil_mc)
                                ):
                                    l_cmp_vari.append(cmp)
                            if vmoye < 0:
                                if (vmin > vmoye * (1.0 - seuil_mc)) or (
                                    vmax < vmoye * (1 + seuil_mc)
                                ):
                                    l_cmp_vari.append(cmp)

                        l_vmoye.append(vmoye)

                    if len(l_cmp_vari) > 0:
                        for cmp in nom_cmps:
                            vu = 0
                            for cmp_vari in l_cmp_vari:
                                if cmp_vari not in masque_mc:
                                    if cmp == cmp_vari and not vu:
                                        if EPSI_MOYENNE[mcfacta]["MAILLE"] is not None:
                                            entite = str(
                                                "MAILLE : " + str(EPSI_MOYENNE[mcfacta]["MAILLE"])
                                            )
                                        if EPSI_MOYENNE[mcfacta]["GROUP_MA"] is not None:
                                            entite = str(
                                                "GROUP_MA : "
                                                + str(EPSI_MOYENNE[mcfacta]["GROUP_MA"])
                                            )
                                        valargs = _F(
                                            vali=[num_ordr[ind]],
                                            valr=[seuil_mc],
                                            valk=[entite, cmp],
                                        )
                                        UTMESS("A", "OBSERVATION_8", **valargs)
                                        vu = 1

                    if TYPE_CHAM[-1:] == "C":
                        mcfactc = {"NOM_CMP": nom_cmps, "NOEUD": l_noeud_mcfact, "VALE_C": l_vmoye}
                    else:
                        mcfactc = {"NOM_CMP": nom_cmps, "NOEUD": l_noeud_mcfact, "VALE": l_vmoye}

                    argsa["AFFE"].append(mcfactc)

                __chame[ind] = CREA_CHAMP(
                    OPERATION="AFFE",
                    MODELE=MODELE_1,
                    TYPE_CHAM=TYPE_CHAM,
                    OPTION="EPSI_NOEU",
                    **argsa
                )

                if RESULTAT.getType() == "MODE_MECA":
                    mcfact2 = {
                        "CHAM_GD": __chame[ind],
                        "MODELE": MODELE_1,
                        "NUME_MODE": int(afreq[ind, 0]),
                        "FREQ": afreq[ind, 1],
                    }

                if RESULTAT.getType() == "EVOL_ELAS":
                    mcfact2 = {"CHAM_GD": __chame[ind], "MODELE": MODELE_1, "INST": list_inst[ind]}

                if RESULTAT.getType() == "DYNA_TRANS":
                    mcfact2 = {"CHAM_GD": __chame[ind], "MODELE": MODELE_1, "INST": list_inst[ind]}

                if RESULTAT.getType() == "DYNA_HARMO":
                    mcfact2 = {"CHAM_GD": __chame[ind], "MODELE": MODELE_1, "FREQ": list_freq[ind]}

                if cham_mater is not None:
                    mcfact2["CHAM_MATER"] = cham_mater
                if cara_elem is not None:
                    mcfact2["CARA_ELEM"] = cara_elem
                mcfact2["NOM_CHAM"] = "EPSI_NOEU"

                liste.append(mcfact2)

            resu_epsi = "EPSI_NOEU"
            RESU = CREA_RESU(OPERATION="AFFE", TYPE_RESU=TYPE_RESU, AFFE=liste)

    # ***********************************************
    #  BOUCLE SUR LES NOM_CHAM
    # ***********************************************
    for nomcham in NOM_CHAM:
        if nomcham == "DEPL" or nomcham == "VITE" or nomcham == "ACCE":
            if RESULTAT.getType() == "DYNA_HARMO":
                TYPE_CHAM = "NOEU_DEPL_C"
            else:
                TYPE_CHAM = "NOEU_DEPL_R"
        elif nomcham == "EPSI_NOEU":
            if RESULTAT.getType() == "DYNA_HARMO":
                TYPE_CHAM = "NOEU_EPSI_C"
            else:
                TYPE_CHAM = "NOEU_EPSI_R"
        else:
            UTMESS("F", "ELEMENTS4_48", valk=[nomcham])

        # ***********************************************
        #  PHASE DE PROJECTION
        # ***********************************************

        if PROJECTION == "OUI":
            if resu_epsi and nomcham == "EPSI_NOEU":
                __proj = PROJ_CHAMP(
                    METHODE="COLLOCATION",
                    RESULTAT=RESU,
                    MODELE_1=MODELE_1,
                    MODELE_2=MODELE_2,
                    NOM_CHAM=nomcham,
                    **mcfact
                )
            else:
                __proj = PROJ_CHAMP(
                    METHODE="COLLOCATION",
                    RESULTAT=RESULTAT,
                    MODELE_1=MODELE_1,
                    MODELE_2=MODELE_2,
                    NOM_CHAM=nomcham,
                    **mcfact
                )
            modele = MODELE_2
        else:
            if resu_epsi and nomcham == "EPSI_NOEU":
                __proj = RESU
            else:
                __proj = RESULTAT
            modele = MODELE_1

        # ***********************************************
        #  PHASE DE CHANGEMENT DE REPERE
        # ***********************************************
        # Le changement de repere se fait dans les routines exterieures
        # crea_normale et crea_repere
        # On range dans le mcfact MODI_REPERE une liste de modifications. ex :
        #     MODI_REPERE = ( _F( GROUP_NO  = toto,
        #                         REPERE    = 'NORMALE' ),
        #                         CONDITION = (1.0,0.0,0.0),
        #                         NOM_PARA  = 'X')
        #                     _F( NOEUD   = ('a','z','e'...),
        #                         REPERE  = 'CYLINDRIQUE',
        #                         ORIGINE = (0.,0.,0.),
        #                         AXE_Z   = (0.,1.,0.), ),
        #                     _F( GROUP_NO = titi,
        #                         REPERE   = 'UTILISATEUR',
        #                         ANGL_NAUT = (alpha, beta, gamma), )
        #                    )
        if MODI_REPERE is not None:
            num_ordr = __proj.getIndexes()

            for modif_rep in MODI_REPERE:
                modi_rep = modif_rep
                type_cham = modif_rep["TYPE_CHAM"]
                nom_cmp = modif_rep["NOM_CMP"]

                if type_cham == "TENS_2D":
                    nomchamx = "EPSI_NOEU"
                elif type_cham == "TENS_3D":
                    nomchamx = "EPSI_NOEU"
                else:
                    nomchamx = "DEPL"

                if nomcham == nomchamx:
                    mcfact1 = {"TYPE_CHAM": type_cham, "NOM_CHAM": nomcham}

                    mcfact2 = {}

                    if modi_rep["REPERE"] == "DIR_JAUGE":
                        vect_x = None
                        vect_y = None
                        if "VECT_X" in modi_rep:
                            vect_x = modi_rep["VECT_X"]
                        if "VECT_Y" in modi_rep:
                            vect_y = modi_rep["VECT_Y"]

                        # il faut des mailles pour les tenseurs d'ordre 2
                        for typ in ["MAILLE", "GROUP_MA"]:
                            if typ in modi_rep:
                                if PROJECTION == "OUI":
                                    maya = mayaexp
                                else:
                                    maya = mayanum
                                list_ma = find_ma(maya, {typ: modi_rep[typ]})

                        angl_naut = crea_repere_xy(vect_x, vect_y)

                        mcfact2.update({"MAILLE": list_ma, "ANGL_NAUT": angl_naut})

                        argsm = {"MODI_CHAM": mcfact1, "AFFE": mcfact2}

                        __bidon = MODI_REPERE_CMD(
                            RESULTAT=__proj, REPERE="UTILISATEUR", NUME_ORDRE=num_ordr, **argsm
                        )
                        __proj = __bidon

                    if modi_rep["REPERE"] == "NORMALE":
                        # Cas ou l'utilisateur choisit de creer les reperes locaux
                        # selon la normale. On fait un changement de repere local
                        # par noeud
                        for option in ["VECT_X", "VECT_Y", "CONDITION_X", "CONDITION_Y"]:
                            if option in modi_rep:
                                vect = {option: modi_rep[option]}
                        if len(vect) != 1:
                            UTMESS("E", "UTILITAI7_9")

                        chnorm = crea_normale(
                            self, MODELE_1, MODELE_2, NUME_DDL, cham_mater, cara_elem
                        )
                        modele = MODELE_2
                        chnormx, description = chnorm.getValuesWithDescription("DX")
                        ind_noeuds = description[0]
                        nom_allno = [mayaexp.getNodeName(k) for k in ind_noeuds]

                        # on met les noeuds conernes sous forme de liste et on va
                        # chercher les noeuds des mailles pour 'MAILLE' et
                        # 'GROUP_MA'
                        for typ in ["NOEUD", "GROUP_NO", "MAILLE", "GROUP_MA"]:
                            if typ in modi_rep:
                                list_no_exp = find_no(mayaexp, {typ: modi_rep[typ]})

                        # boucle sur les noeuds pour modifier les reperes.
                        __bid = [None] * (len(list_no_exp) + 1)
                        __bid[0] = __proj
                        k = 0
                        for nomnoe in list_no_exp:
                            ind_no = nom_allno.index(nomnoe)
                            angl_naut = crea_repere(chnorm, ind_no, vect)

                            mcfact2.update({"NOEUD": nomnoe, "ANGL_NAUT": angl_naut})
                            argsm = {"MODI_CHAM": mcfact1, "AFFE": mcfact2}
                            __bid[k + 1] = MODI_REPERE_CMD(
                                RESULTAT=__bid[k],
                                REPERE="UTILISATEUR",
                                TOUT_ORDRE="OUI",
                                CRITERE="RELATIF",
                                **argsm
                            )
                            k = k + 1

                        __proj = __bid[-1:][0]

                    if modi_rep["REPERE"] == "UTILISATEUR" or modi_rep["REPERE"] == "CYLINDRIQUE":
                        if type_cham == "TENS_2D" or type_cham == "TENS_3D":
                            for typ in ["MAILLE", "GROUP_MA"]:
                                if typ in modi_rep:
                                    mcfact2.update({typ: modi_rep[typ]})
                        else:
                            for typ in ["NOEUD", "GROUP_NO", "MAILLE", "GROUP_MA"]:
                                if typ in modi_rep:
                                    mcfact2.update({typ: modi_rep[typ]})

                        if modi_rep["REPERE"] == "CYLINDRIQUE":
                            origine = modi_rep["ORIGINE"]
                            axe_z = modi_rep["AXE_Z"]
                            mcfact2.update({"ORIGINE": origine, "AXE_Z": axe_z})

                        elif modi_rep["REPERE"] == "UTILISATEUR":
                            angl_naut = modi_rep["ANGL_NAUT"]
                            mcfact2.update({"ANGL_NAUT": angl_naut})

                        argsm = {
                            "MODI_CHAM": mcfact1,
                            "REPERE": modi_rep["REPERE"],
                            "AFFE": mcfact2,
                        }

                        __bidon = MODI_REPERE_CMD(RESULTAT=__proj, CRITERE="RELATIF", **argsm)
                        __proj = __bidon

        # *************************************************
        # Phase de selection des DDL de mesure
        # *************************************************
        resu_filtre = None
        if FILTRE is not None:
            num_ordr = __proj.getIndexes()

            __chamf = [None] * len(indice)
            if RESULTAT.getType() == "EVOL_ELAS":
                list_inst = __proj.LIST_VARI_ACCES()["INST"]
            if RESULTAT.getType() == "DYNA_TRANS":
                list_inst = __proj.LIST_VARI_ACCES()["INST"]
            if RESULTAT.getType() == "DYNA_HARMO":
                list_freq = __proj.LIST_VARI_ACCES()["FREQ"]

            liste = []

            for ind in indice:
                mcfact2 = {}
                filtres = []
                __chamex = CREA_CHAMP(
                    TYPE_CHAM=TYPE_CHAM,
                    OPERATION="EXTR",
                    RESULTAT=__proj,
                    NOM_CHAM=nomcham,
                    NUME_ORDRE=num_ordr[ind],
                )

                for filtre in FILTRE:
                    try:
                        nomchamx = filtre["NOM_CHAM"]
                    except KeyError:
                        nomchamx = None

                    if nomchamx is None or nomchamx == nomcham:
                        mcfact1 = {}

                        atraiter = None
                        if filtre["DDL_ACTIF"][0][0] == "E" and nomcham == "EPSI_NOEU":
                            atraiter = nomcham
                        elif filtre["DDL_ACTIF"][0][0] == "D" and nomcham == "DEPL":
                            atraiter = nomcham
                        elif filtre["DDL_ACTIF"][0][0] == "D" and nomcham == "VITE":
                            atraiter = nomcham
                        elif filtre["DDL_ACTIF"][0][0] == "D" and nomcham == "ACCE":
                            atraiter = nomcham

                        if atraiter:
                            for typ in ["NOEUD", "GROUP_NO", "MAILLE", "GROUP_MA"]:
                                if typ in filtre:
                                    mcfact1.update({typ: filtre[typ]})
                            mcfact1.update({"NOM_CMP": filtre["DDL_ACTIF"], "CHAM_GD": __chamex})
                            filtres.append(mcfact1)

                if len(filtres) > 0:
                    if nomcham == "DEPL" or nomcham == "VITE" or nomcham == "ACCE":
                        __chamf[ind] = CREA_CHAMP(
                            TYPE_CHAM=TYPE_CHAM, OPERATION="ASSE", MODELE=modele, ASSE=filtres
                        )

                    elif nomcham == "EPSI_NOEU":
                        __chamf[ind] = CREA_CHAMP(
                            TYPE_CHAM=TYPE_CHAM, OPERATION="ASSE", MODELE=modele, ASSE=filtres
                        )
                    else:
                        valk = []
                        valk.append(nomcham)
                        valk.append("DEPL VITE ACCE")
                        valk.append("EPSI_NOEU")
                        UTMESS("F", "OBSERVATION_6", valk)

                    argsr = {}
                    if RESULTAT.getType() == "MODE_MECA":
                        mcfact2 = {
                            "CHAM_GD": __chamf[ind],
                            "MODELE": modele,
                            "NUME_MODE": int(afreq[ind, 0]),
                            "FREQ": afreq[ind, 1],
                        }
                        # on recopie les matrices associees au MODELE_2 dans le resultat final
                        # NB : ce n'est peut-etre pas propre, car ces matrices ne
                        # veulent plus rien dire si on a filtre des DDL !!!!!
                        argsr = {"MATR_RIGI": MATR_RIGI, "MATR_MASS": MATR_MASS}

                    if RESULTAT.getType() == "EVOL_ELAS":
                        mcfact2 = {
                            "CHAM_GD": __chamf[ind],
                            "MODELE": modele,
                            "INST": list_inst[ind],
                        }

                    if RESULTAT.getType() == "DYNA_TRANS":
                        mcfact2 = {
                            "CHAM_GD": __chamf[ind],
                            "MODELE": modele,
                            "INST": list_inst[ind],
                        }

                    if RESULTAT.getType() == "DYNA_HARMO":
                        mcfact2 = {
                            "CHAM_GD": __chamf[ind],
                            "MODELE": MODELE_2,
                            "FREQ": list_freq[ind],
                        }

                    if cham_mater is not None:
                        mcfact2["CHAM_MATER"] = cham_mater
                    if cara_elem is not None:
                        mcfact2["CARA_ELEM"] = cara_elem
                    mcfact2["NOM_CHAM"] = nomcham

                    liste.append(mcfact2)

            if len(filtres) > 0 and len(liste) > 0:
                resu_filtre = nomcham
                if RESU:
                    RESU = CREA_RESU(
                        reuse=RESU,
                        RESULTAT=RESU,
                        OPERATION="AFFE",
                        TYPE_RESU=TYPE_RESU,
                        AFFE=liste,
                        **argsr
                    )
                else:
                    RESU = CREA_RESU(OPERATION="AFFE", TYPE_RESU=TYPE_RESU, AFFE=liste, **argsr)

        # *************************************************
        # Recopie de __proj dans RESU si celle-ci
        # n'a pas encore ete faite via FILTRE
        # *************************************************
        if resu_filtre is None:
            RESU = PROJ_CHAMP(
                METHODE="COLLOCATION",
                RESULTAT=__proj,
                MODELE_1=modele,
                MODELE_2=modele,
                NUME_DDL=NUME_DDL,
                NOM_CHAM=nomcham,
                **mcfact
            )

    return RESU


# **********************************************
# RECUPERATION DES NORMALES
# **********************************************
def crea_normale(self, modele_1, modele_2, nume_ddl, cham_mater=None, cara_elem=None):
    """Cree un champ de vecteurs normaux sur le maillage experimental, par
    projection du champ de normales cree sur le maillage numerique
    les mailles doivent etre des elements de <peau> (facettes)
    """
    # recherche du maillage associe au modele numerique
    mayanum = modele_1.getMesh()

    DEFI_GROUP(reuse=mayanum, MAILLAGE=mayanum, CREA_GROUP_MA=_F(NOM="&&TOUMAI", TOUT="OUI"))

    __norm1 = CREA_CHAMP(
        MODELE=modele_1, OPERATION="NORMALE", TYPE_CHAM="NOEU_GEOM_R", GROUP_MA="&&TOUMAI"
    )

    DEFI_GROUP(reuse=mayanum, MAILLAGE=mayanum, DETR_GROUP_MA=_F(NOM="&&TOUMAI"))

    __norm2 = CREA_CHAMP(
        OPERATION="ASSE",
        TYPE_CHAM="NOEU_DEPL_R",
        MODELE=modele_1,
        ASSE=_F(
            TOUT="OUI", CHAM_GD=__norm1, NOM_CMP=("X", "Y", "Z"), NOM_CMP_RESU=("DX", "DY", "DZ")
        ),
    )

    affe_dct = {"NOM_CHAM": "DEPL", "CHAM_GD": __norm2, "INST": 1, "MODELE": modele_1}
    if cham_mater is not None:
        affe_dct["CHAM_MATER"] = cham_mater
    if cara_elem is not None:
        affe_dct["CARA_ELEM"] = cara_elem

    __norm3 = CREA_RESU(OPERATION="AFFE", TYPE_RESU="EVOL_ELAS", AFFE=_F(**affe_dct))

    __norm4 = PROJ_CHAMP(
        METHODE="COLLOCATION",
        RESULTAT=__norm3,
        MODELE_1=modele_1,
        MODELE_2=modele_2,
        NOM_CHAM="DEPL",
        TOUT_ORDRE="OUI",
        NUME_DDL=nume_ddl,
    )

    # __norm5 : toutes les normales au maillage au niveau des capteurs
    __norm5 = CREA_CHAMP(
        RESULTAT=__norm4, OPERATION="EXTR", NUME_ORDRE=1, NOM_CHAM="DEPL", TYPE_CHAM="NOEU_DEPL_R"
    )

    return __norm5


# **********************************************************************
# Calcul des angles nautiques pour le repere local associe a la normale
# **********************************************************************


def crea_repere(chnorm, ind_no, vect):
    """Creation d'un repere orthonormal a partir du vecteur normale et
    d'une equation supplementaire donnee par l'utilisateur sous forme
    de trois parametres et du vecteur de base concerne.
    """
    nom_para = list(vect.keys())[0]  # nom_para = 'VECT_X' ou 'VECT_Y'

    # 1) pour tous les noeuds du maillage experimental, recuperer la normale
    #    calculee a partir du maillage numerique
    chnormx, _ = chnorm.getValuesWithDescription("DX")
    chnormy, _ = chnorm.getValuesWithDescription("DY")
    chnormz, _ = chnorm.getValuesWithDescription("DZ")

    normale = [chnormx[ind_no], chnormy[ind_no], chnormz[ind_no]]

    # 2.1) soit l'utilisateur a donne un deuxieme vecteur explicitement
    # (option VECT_X Ou VECT_Y). Dans ce cas la, le 3e est le produit
    # vectoriel des deux premiers.
    if nom_para == "VECT_X":
        vect1 = numpy.array(list(vect[nom_para]))  # vect x du reploc
        vect2 = numpy.cross(normale, vect1)
        reploc = numpy.array([vect1.tolist(), vect2.tolist(), normale])
        reploc = numpy.transpose(reploc)

    elif nom_para == "VECT_Y":
        vect2 = numpy.array(list(vect[nom_para]))  # vect y du reploc
        vect1 = numpy.cross(vect2, normale)
        reploc = numpy.array([vect1.tolist(), vect2.tolist(), normale])
        reploc = numpy.transpose(reploc)

    # 2.2) TODO : plutot que de donner explicitement un vecteur du repere
    # local avec VECT_X/Y, on devrait aussi pouvoir donner une condition
    # sous forme d'une equation sur un de ces vecteurs. Par exemple,
    # CONDITION_X = (0.,1.,0.) signifierait que le vecteur X1 verifie
    # x(X1) + y(X1) + z(X1) = 0
    elif nom_para == "CONDITION_X":
        pass
    elif nom_para == "CONDITION_Y":
        pass

    # 3) Calcul de l'angle nautique associe au repere local
    angl_naut = anglnaut(reploc)

    return angl_naut


# **********************************************************************
# Calcul des angles nautiques pour le repere associe a VECT_X et VECT_Y
# **********************************************************************


def crea_repere_xy(vect_x, vect_y):
    """Calcul des angles nautiques a partir des directions vect_x et vect_y.
    Si vect_x is not None et vect_y is not None alors on impose le premier vecteur de base
        colineaire a vect_x et le deuxieme vecteur dans le plan (vect_x,vect_y)
    Si vect_x is not None et vect_y is None alors on impose le premier vecteur de base
        colineaire a vect_x
    Si vect_x is None et vect_y is not None alors on impose le deuxieme vecteur de base
        colineaire a vect_y
    Si vect_x is None et vect_y is None alors on ne fait rien
    """
    if vect_x is None and vect_y is None:
        angl_naut = (0.0, 0.0, 0.0)
    else:
        if vect_x and vect_y:
            vx = numpy.array(list(vect_x))
            vy = numpy.array(list(vect_y))
            vect1 = vx
            vect3 = numpy.cross(vx, vy)
            vect2 = numpy.cross(vect3, vx)

        elif vect_x:
            vx = numpy.array(list(vect_x))
            vy1 = numpy.cross((1.0, 0.0, 0.0), vx)
            vy2 = numpy.cross((0.0, 1.0, 0.0), vx)
            vy3 = numpy.cross((0.0, 0.0, 1.0), vx)
            n1 = norm(vy1)
            n2 = norm(vy2)
            n3 = norm(vy3)
            nmax = max(n1, n2, n3)
            if nmax == n1:
                vy = vy1
            elif nmax == n2:
                vy = vy2
            elif nmax == n3:
                vy = vy3
            else:
                UTMESS("F", "UTILITAI_7")
            vect3 = numpy.cross(vx, vy)
            vect1 = vx
            vect2 = numpy.cross(vect3, vect1)

        elif vect_y:
            vy = numpy.array(list(vect_y))
            vx1 = numpy.cross((1.0, 0.0, 0.0), vy)
            vx2 = numpy.cross((0.0, 1.0, 0.0), vy)
            vx3 = numpy.cross((0.0, 0.0, 1.0), vy)
            n1 = norm(vx1)
            n2 = norm(vx2)
            n3 = norm(vx3)
            nmax = max(n1, n2, n3)
            if nmax == n1:
                vx = vx1
            elif nmax == n2:
                vx = vx2
            elif nmax == n3:
                vx = vx3
            else:
                UTMESS("F", "UTILITAI_7")
            vect3 = numpy.cross(vx, vy)
            vect2 = vy
            vect1 = numpy.cross(vect2, vect3)

        norm12 = numpy.dot(vect1, vect1)
        norm22 = numpy.dot(vect2, vect2)
        norm32 = numpy.dot(vect3, vect3)
        if norm12 == 0 or norm22 == 0 or norm32 == 0:
            UTMESS("F", "UTILITAI_7")
        else:
            reploc = numpy.array([vect1.tolist(), vect2.tolist(), vect3.tolist()])
            reploc = numpy.transpose(reploc)
            angl_naut = anglnaut(reploc)

    return angl_naut


# *****************************************************************************
# Aller chercher une liste de noeuds pour un mot cle 'NOEUD', 'GROUP_NO'
# 'MAILLE' ou 'GROUP_MA'
# *****************************************************************************


def find_no(maya, mcsimp):
    """Si on demande une liste de noeuds, c'est simple, on retourne les noeuds
    Si on demande une liste de groupes de noeuds, on va chercher les noeuds
      dans ces groupes, en faisant attention a ne pas etre redondant
    Si on demande un liste de mailles, on va chercher dans la connectivité
      du maillage les indices, puis les noms des noeuds concernes
    etc...
    """
    list_no = []
    connect = maya.getConnectivity()
    if "GROUP_NO" in mcsimp and type(mcsimp["GROUP_NO"]) not in (tuple, list):
        mcsimp["GROUP_NO"] = [mcsimp["GROUP_NO"]]
    if "MAILLE" in mcsimp and type(mcsimp["MAILLE"]) not in (tuple, list):
        mcsimp["MAILLE"] = [mcsimp["MAILLE"]]
    if "GROUP_MA" in mcsimp and type(mcsimp["GROUP_MA"]) not in (tuple, list):
        mcsimp["GROUP_MA"] = [mcsimp["GROUP_MA"]]

    if "NOEUD" in mcsimp:
        list_no = list(mcsimp["NOEUD"])
    elif "GROUP_NO" in mcsimp:
        for group in mcsimp["GROUP_NO"]:
            list_ind_no = maya.getNodes(group)
            for ind_no in list_ind_no:
                nomnoe = maya.getNodeName(ind_no)
                if nomnoe not in list_no:
                    list_no.append(nomnoe)
    elif "MAILLE" in mcsimp:
        for nu_ma in maya.getCells():
            if maya.getCellName(nu_ma) in mcsimp["MAILLE"]:
                for ind_no in connect[nu_ma]:
                    nomnoe = maya.getNodeName(ind_no)
                    if nomnoe not in list_no:
                        list_no.append(nomnoe)
    elif "GROUP_MA" in mcsimp:
        for group in mcsimp["GROUP_MA"]:
            list_nu_ma = maya.getCells(group)
            for nu_ma in list_nu_ma:
                for ind_no in connect[nu_ma]:
                    nomnoe = maya.getNodeName(ind_no)
                    if nomnoe not in list_no:
                        list_no.append(nomnoe)

    return list_no


# *****************************************************************************
# Aller chercher une liste de mailles pour un mot cle 'MAILLE' ou 'GROUP_MA'
# *****************************************************************************


def find_ma(maya, mcsimp):
    """Si mot cle MAILLE, on retourne la liste des mailles
    Si mot cle GROUP_MA, on va chercher les mailles dans ces groupes
    """
    list_ma = []
    if "GROUP_MA" in mcsimp and type(mcsimp["GROUP_MA"]) != tuple:
        mcsimp["GROUP_MA"] = [mcsimp["GROUP_MA"]]
    if "MAILLE" in mcsimp and type(mcsimp["MAILLE"]) != tuple:
        mcsimp["MAILLE"] = [mcsimp["MAILLE"]]

    if "MAILLE" in mcsimp:
        for mail in mcsimp["MAILLE"]:
            list_ma.append(mail)
    elif "GROUP_MA" in mcsimp:
        for group in mcsimp["GROUP_MA"]:
            list_ind_ma = maya.getCells(group)
            for ind_ma in list_ind_ma:
                nommail = maya.getCellName(ind_ma)
                if nommail not in list_ma:
                    list_ma.append(nommail)

    return list_ma


# ****************************************************
# Quelques utilitaires de calculs d'angles nautiques
# ****************************************************


def norm(x):
    """Calcul de la norme euclidienne d'un vecteur"""
    tmp = numpy.sqrt(numpy.dot(x, x))
    return tmp


def anglnaut(P):
    """Calcule les angles nautiques correspondant a un repere local
    NB : seuls les deux premiers vecteurs de P (les images respectives
    de X et Y) sont utiles pour le calcul des angles
    """
    # expression des coordonnees globales des 3 vecteurs de base locale
    y = numpy.array([0.0, 1.0, 0.0])

    xg = P[:, 0]
    yg = P[:, 1]

    # calcul des angles nautiques
    x1 = copy.copy(xg)
    # x1 projection de xg sur la plan xy, non norme
    x1[2] = 0.0
    # produit scalaire X xg
    normx = norm(x1)
    if normx == 0.0:  # on impose alpha = 0 pour lever l'indetermination
        COSA = 1.0
        SINA = 0.0
    else:
        COSA = x1[0] / normx
        # produit vectoriel X xg
        SINA = x1[1] / normx
    ar = numpy.arctan2(SINA, COSA)
    alpha = ar * 180 / numpy.pi

    COSB = norm(x1)
    SINB = -xg[2]
    beta = numpy.arctan2(SINB, COSB) * 180 / numpy.pi

    P2 = numpy.zeros((3, 3))
    P2[0, 0] = numpy.cos(ar)
    P2[1, 0] = numpy.sin(ar)
    P2[1, 1] = numpy.cos(ar)
    P2[0, 1] = -numpy.sin(ar)
    y1 = numpy.dot(P2, y)
    y1n = y1 / norm(y1)

    # calcul de gamma
    COSG = numpy.dot(y1n, yg)
    SING = numpy.dot(xg, numpy.cross(y1n, yg))
    gamma = numpy.arctan2(SING, COSG) * 180 / numpy.pi

    return alpha, beta, gamma


# NB : Equations de passage : un vecteur de coordonnees globales (X,Y,Z) a pour
# coordonnees locales (X1,Y1,Z1) avec
# _                  _  _                   _  _                   _  _ _     _  _
# | 1     0      0     || cos(B) 0    -sin(B) ||  cos(A)  sin(A)   0 || X |   | X1 |
# | 0   cos(G)  sin(G) ||   0    1      0     || -sin(A)  cos(A)   0 || Y | = | Y1 |
# |_0  -sin(G)  cos(G)_||_sin(B) 0     cos(B)_||_   0       0      1_||_Z_|   |_Z1_|
#
# A (alpha), B(beta), gamma (G) sont les angle nautiques que l'on donne habituellemet
# dans les MODI_REPERE. Les equations a resoudre sont les suivantes :
# cos(A)cos(B)                      = reploc[0][0]
# -cos(G)sin(A) + sin(G)cos(A)sin(B) = reploc[0][1]
# sin(A)sin(G) + cos(A)sin(B)cos(G) = reploc[0][2]
#
# sin(A)cos(B)                      = reploc[1][0]
# cos(A)cos(G) + sin(A)sin(B)sin(G) = reploc[1][1]
# -cos(A)sin(G) + sin(A)sin(B)cos(G) = reploc[1][2]
#
# -sin(B) = reploc[2][0]
# cos(B)sin(G) = reploc[2][1]
# cos(B)cos(G) = reploc[2][2]
