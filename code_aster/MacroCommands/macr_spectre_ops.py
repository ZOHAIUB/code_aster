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

from ..Cata.Syntax import _F
from ..CodeCommands import CALC_FONCTION, CREA_TABLE, IMPR_FONCTION, RECU_FONCTION
from ..Messages import UTMESS


def macr_spectre_ops(self, **args):
    """
    Ecriture de la macro MACR_SPECTRE
    """
    MAILLAGE = args.get("MAILLAGE")
    PLANCHER = args.get("PLANCHER")
    NOM_CHAM = args.get("NOM_CHAM")
    ENVELOPPE = args.get("ENVELOPPE")
    CALCUL = args.get("CALCUL")
    RESU = args.get("RESU")
    IMPRESSION = args.get("IMPRESSION")
    FREQ = args.get("FREQ")
    LIST_FREQ = args.get("LIST_FREQ")
    LIST_INST = args.get("LIST_INST")
    AMOR_SPEC = args.get("AMOR_SPEC")

    EnumType = (list, tuple)

    # On importe les definitions des commandes a utiliser dans la macro

    # construction de la liste des noeuds à traiter
    planch_nodes = {}
    l_plancher = []
    l_batiment = []
    l_commentaire = []
    #
    dplancher = []
    for j in PLANCHER:
        dplancher.append(j.cree_dict_valeurs(j.mc_liste))
    #
    for plancher in dplancher:
        liste_no = []
        clefs = list(plancher.keys())
        if "NOEUD" in clefs:
            if plancher["NOEUD"] is not None:
                if type(plancher["NOEUD"]) == str:
                    liste_no.append(plancher["NOEUD"])
                else:
                    for noeud in plancher["NOEUD"]:
                        liste_no.append(noeud)
        if "GROUP_NO" in clefs:
            if plancher["GROUP_NO"] is not None:
                assert MAILLAGE is not None
                if type(plancher["GROUP_NO"]) == str:
                    noms_no = [
                        MAILLAGE.getNodeName(n) for n in MAILLAGE.getNodes(plancher["GROUP_NO"])
                    ]
                    liste_no = liste_no + noms_no
                else:
                    for group_no in plancher["GROUP_NO"]:
                        noms_no = [MAILLAGE.getNodeName(n) for n in MAILLAGE.getNodes(group_no)]
                        liste_no = liste_no + noms_no
        planch_nodes[plancher["NOM"]] = liste_no
        l_plancher.append(plancher["NOM"])
        l_batiment.append(plancher["BATIMENT"])
        l_commentaire.append(plancher["COMMENTAIRE"])

    if AMOR_SPEC is not None and type(AMOR_SPEC) not in EnumType:
        AMOR_SPEC = (AMOR_SPEC,)
    #
    if NOM_CHAM == "ACCE":
        dico_glob = {}
    if NOM_CHAM == "DEPL":
        if ENVELOPPE == "OUI":
            dico_glob = {"DX_max": [], "DY_max": [], "DZ_max": [], "DH_max": []}
        else:
            dico_glob = {"DX_max": [], "DY_max": [], "DZ_max": []}
    #
    # ---------------------------------------------------------------------------------------------
    # boucle 1 sur les planchers
    for plancher in l_plancher:
        #
        if NOM_CHAM == "ACCE":
            __moy_x = [None] * len(planch_nodes[plancher])
            __moy_y = [None] * len(planch_nodes[plancher])
            __moy_z = [None] * len(planch_nodes[plancher])
        if NOM_CHAM == "DEPL":
            dicDmax = {}
        # -----------------------------------------------------------------------------------------
        # boucle 2 sur les noeuds du plancher
        indexn = 0
        for node in planch_nodes[plancher]:
            # -------------------------------------------------------------------------------------
            # boucle 3 sur les directions (X,Y,Z)
            for dd in ("X", "Y", "Z"):
                # ---------------------------------------------------------------------------------
                # boucle 4 sur les résultats
                l_fonc = []
                for resu in RESU:
                    # Etape 1: Récupération des fonctions
                    motscles = {}
                    if resu["RESU_GENE"] is not None:
                        if CALCUL == "ABSOLU" and args.get("MULT_APPUI") is not None:
                            motscles["MULT_APPUI"] = args["MULT_APPUI"]
                        #
                        motscles["RESU_GENE"] = resu["RESU_GENE"]
                        __spo = RECU_FONCTION(
                            NOM_CHAM=NOM_CHAM,
                            TOUT_ORDRE="OUI",
                            NOM_CMP="D" + dd,
                            INTERPOL="LIN",
                            PROL_GAUCHE="CONSTANT",
                            PROL_DROITE="CONSTANT",
                            NOEUD=node,
                            **motscles
                        )
                    elif args.get("MULT_APPUI") is not None:
                        UTMESS("F", "SPECTRAL0_13")
                    #
                    if resu["RESULTAT"] is not None:
                        motscles["RESULTAT"] = resu["RESULTAT"]
                        __spo = RECU_FONCTION(
                            NOM_CHAM=NOM_CHAM,
                            TOUT_ORDRE="OUI",
                            NOM_CMP="D" + dd,
                            INTERPOL="LIN",
                            PROL_GAUCHE="CONSTANT",
                            PROL_DROITE="CONSTANT",
                            NOEUD=node,
                            **motscles
                        )
                    #
                    if resu["TABLE"] is not None:
                        # 2 formats de table possible. Avec les colonnes :
                        #   INST NOEUD NOM_CHAM NOM_CMP VALE
                        #   INST NOEUD NOM_CHAM DX DY DZ
                        # récupération du nom des colonnes de la table
                        nomcol = resu["TABLE"].get_nom_para()
                        #
                        lst1 = ("INST", "NOEUD", "NOM_CHAM", "NOM_CMP", "VALE")
                        ok1 = True
                        for para in lst1:
                            ok1 = ok1 and (para in nomcol)
                        #
                        lst2 = ("INST", "NOEUD", "NOM_CHAM", "D" + dd)
                        ok2 = True
                        for para in lst2:
                            ok2 = ok2 and (para in nomcol)
                        #
                        if not ok1 ^ ok2:
                            print(nomcol)
                            assert ok1 ^ ok2
                        #
                        if ok1:
                            __spo = RECU_FONCTION(
                                TABLE=resu["TABLE"],
                                PARA_X="INST",
                                PARA_Y="VALE",
                                INTERPOL="LIN",
                                FILTRE=(
                                    _F(NOM_PARA="NOEUD", VALE_K=node),
                                    _F(NOM_PARA="NOM_CHAM", VALE_K=NOM_CHAM),
                                    _F(NOM_PARA="NOM_CMP", VALE_K="D" + dd),
                                ),
                            )
                        #
                        if ok2:
                            __spo = RECU_FONCTION(
                                TABLE=resu["TABLE"],
                                PARA_X="INST",
                                PARA_Y="D" + dd,
                                INTERPOL="LIN",
                                FILTRE=(
                                    _F(NOM_PARA="NOEUD", VALE_K=node),
                                    _F(NOM_PARA="NOM_CHAM", VALE_K=NOM_CHAM),
                                ),
                            )
                    #
                    # Etape 2: Combinaisons
                    if NOM_CHAM == "ACCE":
                        # Accelerations relatives
                        if CALCUL == "RELATIF":
                            # Combinaison avec fonction d'accélération
                            motscles = {}
                            if LIST_INST is not None:
                                motscles["LIST_PARA"] = LIST_INST
                            __spo = CALC_FONCTION(
                                COMB=(
                                    _F(FONCTION=__spo, COEF=1.0),
                                    _F(FONCTION=resu["ACCE_" + dd], COEF=1.0),
                                ),
                                **motscles
                            )

                        # Etape 3: Calcul des spectres d'oscillateur
                        motscles = {}
                        if FREQ is not None:
                            motscles["FREQ"] = FREQ
                        if LIST_FREQ is not None:
                            motscles["LIST_FREQ"] = LIST_FREQ
                        __spo = CALC_FONCTION(
                            SPEC_OSCI=_F(
                                FONCTION=__spo,
                                AMOR_REDUIT=AMOR_SPEC,
                                NORME=args["NORME"],
                                **motscles
                            )
                        )
                        l_fonc.append(__spo)

                    elif NOM_CHAM == "DEPL":
                        if CALCUL == "ABSOLU":
                            # On retranche les deplacements d entrainement
                            motscles = {}
                            if LIST_INST is not None:
                                motscles["LIST_PARA"] = LIST_INST
                            __spo = CALC_FONCTION(
                                COMB=(
                                    _F(FONCTION=__spo, COEF=1.0),
                                    _F(FONCTION=resu["DEPL_" + dd], COEF=-1.0),
                                ),
                                **motscles
                            )
                        l_fonc.append(__spo)

                # fin boucle 4 sur les résultats
                # ---------------------------------------------------------------------------------
                #
                # calcul de la moyenne sur les resultats à noeud et direction
                # fixes
                nbresu = len(RESU)
                if NOM_CHAM == "ACCE":
                    mcfCMBx = []
                    mcfCMBy = []
                    mcfCMBz = []
                    for spo in l_fonc:
                        mcfCMBx.append(_F(FONCTION=spo, COEF=1.0 / float(nbresu)))
                        mcfCMBy.append(_F(FONCTION=spo, COEF=1.0 / float(nbresu)))
                        mcfCMBz.append(_F(FONCTION=spo, COEF=1.0 / float(nbresu)))
                    motscles = {}
                    if LIST_FREQ is not None:
                        motscles["LIST_PARA"] = LIST_FREQ
                    if dd == "X":
                        __moy_x[indexn] = CALC_FONCTION(COMB=mcfCMBx, **motscles)
                    if dd == "Y":
                        __moy_y[indexn] = CALC_FONCTION(COMB=mcfCMBy, **motscles)
                    if dd == "Z":
                        __moy_z[indexn] = CALC_FONCTION(COMB=mcfCMBz, **motscles)
                #
                elif NOM_CHAM == "DEPL":
                    moy = 0.0
                    for spo in l_fonc:
                        fspo = spo.convert()
                        aspo = fspo.abs()
                        vmax = aspo.extreme()["max"]
                        moy = moy + vmax[-1][-1]
                    dicDmax[(node, dd)] = moy / nbresu
            # fin boucle 3 sur les directions
            # -------------------------------------------------------------------------------------
            #
            # impressions en chaque noeud
            if NOM_CHAM == "ACCE" and IMPRESSION is not None:
                if IMPRESSION["TOUT"] == "OUI":
                    __moyxa = [None] * len(AMOR_SPEC)
                    __moyya = [None] * len(AMOR_SPEC)
                    __moyza = [None] * len(AMOR_SPEC)
                    for i in range(len(AMOR_SPEC)):
                        __moyxa[i] = RECU_FONCTION(
                            NAPPE=__moy_x[indexn], VALE_PARA_FONC=AMOR_SPEC[i]
                        )
                        __moyya[i] = RECU_FONCTION(
                            NAPPE=__moy_y[indexn], VALE_PARA_FONC=AMOR_SPEC[i]
                        )
                        __moyza[i] = RECU_FONCTION(
                            NAPPE=__moy_z[indexn], VALE_PARA_FONC=AMOR_SPEC[i]
                        )
                    motscles = {}
                    dI = IMPRESSION[0].cree_dict_valeurs(IMPRESSION[0].mc_liste)
                    if "PILOTE" in dI:
                        motscles["PILOTE"] = IMPRESSION["PILOTE"]
                    if IMPRESSION["FORMAT"] != "TABLEAU":
                        motscles["ECHELLE_X"] = "LOG"
                    if IMPRESSION["TRI"] == "AMOR_SPEC":
                        for i in range(len(AMOR_SPEC)):
                            TITRE = (
                                "Spectres / Plancher = "
                                + plancher
                                + " / amor="
                                + str(AMOR_SPEC[i])
                                + " / noeud="
                                + node
                            )
                            IMPR_FONCTION(
                                FORMAT=IMPRESSION["FORMAT"],
                                UNITE=IMPRESSION["UNITE"],
                                COURBE=(
                                    _F(FONCTION=__moyxa[i], LEGENDE="X"),
                                    _F(FONCTION=__moyya[i], LEGENDE="Y"),
                                    _F(FONCTION=__moyza[i], LEGENDE="Z"),
                                ),
                                TITRE=TITRE,
                                **motscles
                            )
                    elif IMPRESSION["TRI"] == "DIRECTION":
                        for dd in ("X", "Y", "Z"):
                            TITRE = (
                                "Spectres / Plancher = "
                                + plancher
                                + " / direction = "
                                + dd
                                + " / noeud = "
                                + node
                            )
                            legende = "amor=" + str(AMOR_SPEC[i])
                            if dd == "X":
                                l_fonc = [
                                    _F(FONCTION=__moyxa[i], LEGENDE=legende)
                                    for i in range(len(AMOR_SPEC))
                                ]
                            if dd == "Y":
                                l_fonc = [
                                    _F(FONCTION=__moyya[i], LEGENDE=legende)
                                    for i in range(len(AMOR_SPEC))
                                ]
                            if dd == "Z":
                                l_fonc = [
                                    _F(FONCTION=__moyza[i], LEGENDE=legende)
                                    for i in range(len(AMOR_SPEC))
                                ]
                            IMPR_FONCTION(
                                FORMAT=IMPRESSION["FORMAT"],
                                UNITE=IMPRESSION["UNITE"],
                                COURBE=l_fonc,
                                TITRE=TITRE,
                                **motscles
                            )

            # increment de l'indice de noeud
            indexn = indexn + 1
        # fin boucle 2 sur les noeuds du plancher
        # -----------------------------------------------------------------------------------------
        #
        # Etape 4: Calcul des enveloppes des spectres ou des deplacements max
        if NOM_CHAM == "ACCE":
            mcslx = []
            mcsly = []
            mcslz = []
            indexn = 0
            for node in planch_nodes[plancher]:
                mcslx.append(__moy_x[indexn])
                mcsly.append(__moy_y[indexn])
                mcslz.append(__moy_z[indexn])
                indexn = indexn + 1
            __snx = CALC_FONCTION(ENVELOPPE=_F(FONCTION=mcslx))
            __sny = CALC_FONCTION(ENVELOPPE=_F(FONCTION=mcsly))
            __snz = CALC_FONCTION(ENVELOPPE=_F(FONCTION=mcslz))
            if ENVELOPPE == "OUI":
                __snh = CALC_FONCTION(ENVELOPPE=_F(FONCTION=(__snx, __sny)))
        elif NOM_CHAM == "DEPL":
            DRmX = max([dicDmax[(node, "X")] for node in planch_nodes[plancher]])
            DRmY = max([dicDmax[(node, "Y")] for node in planch_nodes[plancher]])
            DRmZ = max([dicDmax[(node, "Z")] for node in planch_nodes[plancher]])
            if ENVELOPPE == "OUI":
                DRmH = max([DRmX, DRmY])
        #
        # Renseignement de la table finale des résultats
        if NOM_CHAM == "ACCE":
            nbind = len(AMOR_SPEC)
            for i in range(nbind):
                dico_glob["FREQ"] = __snx.Valeurs()[1][i][0]
                dico_glob["eX_%d_%s" % (i, plancher)] = __snx.Valeurs()[1][i][1]
                dico_glob["eY_%d_%s" % (i, plancher)] = __sny.Valeurs()[1][i][1]
                dico_glob["eZ_%d_%s" % (i, plancher)] = __snz.Valeurs()[1][i][1]
                if ENVELOPPE == "OUI":
                    dico_glob["eH_%d_%s" % (i, plancher)] = __snh.Valeurs()[1][i][1]
        elif NOM_CHAM == "DEPL":
            dico_glob["DX_max"].append(DRmX)
            dico_glob["DY_max"].append(DRmY)
            dico_glob["DZ_max"].append(DRmZ)
            if ENVELOPPE == "OUI":
                dico_glob["DH_max"].append(DRmH)
        #
        # Etape 5: Impression des courbes
        if NOM_CHAM == "ACCE" and IMPRESSION is not None:
            motscles = {}
            dI = IMPRESSION[0].cree_dict_valeurs(IMPRESSION[0].mc_liste)
            if "PILOTE" in dI:
                motscles["PILOTE"] = IMPRESSION["PILOTE"]
            if IMPRESSION["FORMAT"] != "TABLEAU":
                motscles["ECHELLE_X"] = "LOG"
            __snxa = [None] * len(AMOR_SPEC)
            __snya = [None] * len(AMOR_SPEC)
            __snza = [None] * len(AMOR_SPEC)
            if ENVELOPPE == "OUI":
                __snha = [None] * len(AMOR_SPEC)
            for i in range(nbind):
                __snxa[i] = RECU_FONCTION(NAPPE=__snx, VALE_PARA_FONC=AMOR_SPEC[i])
                __snya[i] = RECU_FONCTION(NAPPE=__sny, VALE_PARA_FONC=AMOR_SPEC[i])
                __snza[i] = RECU_FONCTION(NAPPE=__snz, VALE_PARA_FONC=AMOR_SPEC[i])
                if ENVELOPPE == "OUI":
                    __snha[i] = RECU_FONCTION(NAPPE=__snh, VALE_PARA_FONC=AMOR_SPEC[i])
            if IMPRESSION["TRI"] == "AMOR_SPEC":
                for i in range(nbind):
                    TITRE = (
                        "Spectres moyens / Plancher = " + plancher + " / amor=" + str(AMOR_SPEC[i])
                    )
                    if ENVELOPPE == "OUI":
                        IMPR_FONCTION(
                            FORMAT=IMPRESSION["FORMAT"],
                            UNITE=IMPRESSION["UNITE"],
                            COURBE=(
                                _F(FONCTION=__snxa[i], LEGENDE="X"),
                                _F(FONCTION=__snya[i], LEGENDE="Y"),
                                _F(FONCTION=__snza[i], LEGENDE="Z"),
                                _F(FONCTION=__snha[i], LEGENDE="H"),
                            ),
                            TITRE=TITRE,
                            **motscles
                        )
                    else:
                        IMPR_FONCTION(
                            FORMAT=IMPRESSION["FORMAT"],
                            UNITE=IMPRESSION["UNITE"],
                            COURBE=(
                                _F(FONCTION=__snxa[i], LEGENDE="X"),
                                _F(FONCTION=__snya[i], LEGENDE="Y"),
                                _F(FONCTION=__snza[i], LEGENDE="Z"),
                            ),
                            TITRE=TITRE,
                            **motscles
                        )
            elif IMPRESSION["TRI"] == "DIRECTION":
                if ENVELOPPE == "OUI":
                    liste_dir = ("X", "Y", "Z", "H")
                else:
                    liste_dir = ("X", "Y", "Z")

                for dd in liste_dir:
                    TITRE = "Spectres moyens / Plancher = " + plancher + " / direction = " + dd
                    legende = "amor=" + str(AMOR_SPEC[i])
                    l_fonc = []
                    if dd == "X":
                        l_fonc = [
                            _F(FONCTION=__snxa[i], LEGENDE=legende) for i in range(len(AMOR_SPEC))
                        ]
                    if dd == "Y":
                        l_fonc = [
                            _F(FONCTION=__snya[i], LEGENDE=legende) for i in range(len(AMOR_SPEC))
                        ]
                    if dd == "Z":
                        l_fonc = [
                            _F(FONCTION=__snza[i], LEGENDE=legende) for i in range(len(AMOR_SPEC))
                        ]
                    if ENVELOPPE == "OUI":
                        if dd == "H":
                            l_fonc = [
                                _F(FONCTION=__snha[i], LEGENDE=legende)
                                for i in range(len(AMOR_SPEC))
                            ]
                    IMPR_FONCTION(
                        FORMAT=IMPRESSION["FORMAT"],
                        UNITE=IMPRESSION["UNITE"],
                        COURBE=l_fonc,
                        TITRE=TITRE,
                        **motscles
                    )
    # fin boucle 1 sur les planchers
    # ---------------------------------------------------------------------------------------------
    #
    # Etape6 : Renseignement de la table finale des résultats
    lListe = []
    nb_amor = 0
    if NOM_CHAM == "DEPL":
        titre = "Calcul des spectres enveloppes"
    elif NOM_CHAM == "ACCE":
        titre = "Calcul des spectres enveloppes par planchers"
        infos_amor = {}
        infos_amor["NUME_AMOR"] = []
        infos_amor["AMOR"] = []
        nb_amor = len(AMOR_SPEC)
        for i in range(len(AMOR_SPEC)):
            infos_amor["NUME_AMOR"].append(i)
            infos_amor["AMOR"].append(AMOR_SPEC[i])

    nb_plancher = len(l_plancher)
    lkeys = list(dico_glob.keys())
    lkeys.sort()
    for key in lkeys:
        nb_lignes = len(dico_glob[key])
        lListe.append(
            _F(
                LISTE_R=dico_glob[key],
                PARA=key,
                NUME_LIGN=list(
                    range(nb_amor + nb_plancher + 1, nb_amor + nb_plancher + nb_lignes + 1)
                ),
            )
        )
    if NOM_CHAM == "ACCE":
        lListe.append(
            _F(
                LISTE_I=infos_amor["NUME_AMOR"],
                PARA="NUME_AMOR",
                NUME_LIGN=list(range(nb_plancher + 1, nb_plancher + nb_amor + 1)),
            )
        )
        lListe.append(
            _F(
                LISTE_R=infos_amor["AMOR"],
                PARA="AMOR",
                NUME_LIGN=list(range(nb_plancher + 1, nb_plancher + nb_amor + 1)),
            )
        )

    lListe.append(_F(LISTE_K=l_plancher, TYPE_K="K24", PARA="NOM"))
    l_bat = [i for i in l_batiment if i is not None]
    l_com = [i for i in l_commentaire if i is not None]
    if l_bat != []:
        l_bat2 = ["-" if i is None else i for i in l_batiment]
        lListe.append(_F(LISTE_K=l_bat2, TYPE_K="K24", PARA="BATIMENT"))
    if l_com != []:
        l_com2 = ["-" if i is None else i for i in l_commentaire]
        lListe.append(_F(LISTE_K=l_com2, TYPE_K="K24", PARA="COMMENTAIRE"))

    tab = CREA_TABLE(LISTE=lListe, TITRE=titre)
    return tab
