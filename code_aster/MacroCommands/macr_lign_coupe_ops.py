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

import os
from math import asin, atan2, cos, pi, sin, sqrt

from ..Cata.Syntax import _F
from ..CodeCommands import (
    COPIER,
    CREA_RESU,
    CREA_TABLE,
    DEFI_GROUP,
    LIRE_MAILLAGE,
    MODI_REPERE,
    POST_RELEVE_T,
    PROJ_CHAMP,
)
from ..Helpers import FileAccess, LogicalUnitFile
from ..Messages import UTMESS, MasquerAlarme, RetablirAlarme

#
# script PYTHON de creation du résultat local
#

#
# verification que les points de la ligne de coupe sont dans la matiere


def crea_grp_matiere(groupe, newgrp, iocc, m, __remodr, NOM_CHAM, __macou):
    motscles = {}
    if m["NOM_CMP"]:
        motscles["NOM_CMP"] = m["NOM_CMP"]
    else:
        motscles["TOUT_CMP"] = "OUI"
    motscles["OPERATION"] = "EXTRACTION"

    __tab = POST_RELEVE_T(
        ACTION=_F(
            INTITULE=newgrp, RESULTAT=__remodr, NOM_CHAM=NOM_CHAM, GROUP_NO=groupe, **motscles
        )
    )

    node_by_name = {__macou.getNodeName(node): node for node in __macou.getNodes(groupe)}

    # dictb=table initiale (contenant éventuellement des noeuds hors matière)
    dictb = __tab.EXTR_TABLE()
    # listenoe_b=liste ordonnee des noeuds de la ligne de coupe (avec doublons)
    listenoe_b = dictb.NOEUD.values()

    # dictc=table (extraite de dictb) contenant uniquement des noeuds dans la
    # matière
    if m["NOM_CMP"]:
        dictc = getattr(dictb, m["NOM_CMP"][0]).NON_VIDE()
        lno_c2 = set(dictc.NOEUD.values())
    else:  # TOUT_CMP='OUI'
        # on garde uniquement les composantes pour conserver les noeuds où il y
        # a des valeurs
        a_suppr = set(
            [
                "INTITULE",
                "RESU",
                "NOM_CHAM",
                "NUME_ORDRE",
                "INST",
                "ABSC_CURV",
                "COOR_X",
                "COOR_Y",
                "COOR_Z",
            ]
        )
        new_para = set(dictb.para)
        new_para.difference_update(a_suppr)

        lno_c2 = set()
        for comp in new_para.difference(["NOEUD"]):
            dictc = getattr(dictb, comp).NON_VIDE()
            lno_c2.update(dictc.NOEUD.values())

    # on réordonne la liste des noeuds de lno_c2 (selon leur position dans listenoe_b) => l_matiere
    # l_horsmat=liste des noeuds hors matière
    l_matiere = [j for j in listenoe_b if j in lno_c2]
    nderm = l_matiere.index(l_matiere[len(l_matiere) - 1])
    l_horsmat = [j for j in listenoe_b if j not in lno_c2]

    # si on est en présence de noeuds hors matière,
    # on emet une alarme pour informer l'utilisateur
    nbpoin = m["NB_POINTS"]
    reste = nbpoin - len(l_matiere)
    if len(l_horsmat) > 0:
        nderh = l_horsmat.index(l_horsmat[len(l_horsmat) - 1])
        coord = __macou.getCoordinates()
        indent = os.linesep + " " * 12
        l_surlig = []
        l_horslig = []
        for j in l_matiere[: nderm + 1]:
            node = node_by_name[j]
            text_coordo = "(%f, %f, %f)" % (
                coord.getNode(node).x(),
                coord.getNode(node).y(),
                coord.getNode(node).z(),
            )
            l_surlig.append(text_coordo)
        for j in l_horsmat[: nderh + 1]:
            node = node_by_name[j]
            text_coordo = "(%f, %f, %f)" % (
                coord.getNode(node).x(),
                coord.getNode(node).y(),
                coord.getNode(node).z(),
            )
            l_horslig.append(text_coordo)
        UTMESS("A", "POST0_8", valk=[indent.join(l_surlig), indent.join(l_horslig)])

    elif reste > 0:
        coord = __macou.getCoordinates()
        indent = os.linesep + " " * 12
        l_surlig = []
        for j in l_matiere[: nderm + 1]:
            node = node_by_name[j]
            text_coordo = "(%f, %f, %f)" % (
                coord.getNode(node).x(),
                coord.getNode(node).y(),
                coord.getNode(node).z(),
            )
            l_surlig.append(text_coordo)
        UTMESS("A", "POST0_24", vali=[iocc, reste], valk=[indent.join(l_surlig)])

    __macou = DEFI_GROUP(
        reuse=__macou, MAILLAGE=__macou, CREA_GROUP_NO=_F(NOM=newgrp, NOEUD=l_matiere[: nderm + 1])
    )

    return


def crea_resu_local(dime, NOM_CHAM, m, resin, mail, nomgrma):
    epsi = 0.00000001

    if NOM_CHAM == "DEPL":
        if dime == 2:
            TYPE_CHAM = "VECT_2D"
        elif dime == 3:
            TYPE_CHAM = "VECT_3D"
    elif NOM_CHAM in ("SIGM_NOEU", "SIEF_ELNO", "SIGM_ELNO"):
        if dime == 2:
            TYPE_CHAM = "TENS_2D"
        elif dime == 3:
            TYPE_CHAM = "TENS_3D"
    elif NOM_CHAM in ("FLUX_ELNO", "FLUX_NOEU"):
        if dime == 2:
            TYPE_CHAM = "VECT_2D"
        elif dime == 3:
            TYPE_CHAM = "VECT_3D"
        else:
            assert 0
    else:
        assert 0

    if m["TYPE"] == "SEGMENT" and m["REPERE"] != "CYLINDRIQUE":
        if m["REPERE"] == "LOCAL":
            # --- determination des angles nautiques
            cx1 = m["COOR_EXTR"][0] - m["COOR_ORIG"][0]
            cx2 = m["COOR_EXTR"][1] - m["COOR_ORIG"][1]
            cx3 = 0.0
            if dime == 3:
                cx3 = m["COOR_EXTR"][2] - m["COOR_ORIG"][2]
            nvx = sqrt(cx1**2 + cx2**2 + cx3**2)
            if abs(nvx) < epsi:
                UTMESS("F", "POST0_1")
            cx1 = cx1 / nvx
            cx2 = cx2 / nvx
            cx3 = cx3 / nvx
            if m["VECT_Y"]:
                cy1 = m["VECT_Y"][0]
                cy2 = m["VECT_Y"][1]
                cy3 = 0.0
                if dime == 3:
                    cy3 = m["VECT_Y"][2]
            else:
                UTMESS("F", "POST0_50")
            nvy = sqrt(cy1**2 + cy2**2 + cy3**2)
            if abs(nvy) < epsi:
                UTMESS("F", "POST0_2")
            cy1 = cy1 / nvy
            cy2 = cy2 / nvy
            cy3 = cy3 / nvy
            if (abs(cx1 - cy1) < epsi and abs(cx2 - cy2) < epsi and abs(cx3 - cy3) < epsi) or (
                abs(cx1 + cy1) < epsi and abs(cx2 + cy2) < epsi and abs(cx3 + cy3) < epsi
            ):
                UTMESS("F", "POST0_3")
            if abs(cx1 * cy1 + cx2 * cy2 + cx3 * cy3) > epsi:
                cz1 = cx2 * cy3 - cx3 * cy2
                cz2 = cx3 * cy1 - cx1 * cy3
                cz3 = cx1 * cy2 - cx2 * cy1
                nvz = sqrt(cz1**2 + cz2**2 + cz3**2)
                cz1 = cz1 / nvz
                cz2 = cz2 / nvz
                cz3 = cz3 / nvz
                cy1 = cz2 * cx3 - cz3 * cx2
                cy2 = cz3 * cx1 - cz1 * cx3
                cy3 = cz1 * cx2 - cz2 * cx1
                nvy = sqrt(cy1**2 + cy2**2 + cy3**2)
                cy1 = cy1 / nvy
                cy2 = cy2 / nvy
                cy3 = cy3 / nvy
                UTMESS("A", "POST0_4", valr=[cy1, cy2, cy3])
            else:
                cz1 = cx2 * cy3 - cx3 * cy2
                cz2 = cx3 * cy1 - cx1 * cy3
                cz3 = cx1 * cy2 - cx2 * cy1
            beta = 0.0
            gamma = 0.0
            if dime == 2:
                alpha = atan2(cx2, cx1)
            else:
                if cx1**2 + cx2**2 > epsi:
                    alpha = atan2(cx2, cx1)
                    beta = -asin(cx3)
                    gamma = atan2(cy3, cz3)
                else:
                    alpha = atan2(-cy1, cy2)
                    beta = -asin(cx3)
                    gamma = 0.0
            alpha = alpha * 180 / pi
            beta = beta * 180 / pi
            gamma = gamma * 180 / pi

        elif m["REPERE"] == "UTILISATEUR":
            alpha = m["ANGL_NAUT"][0]
            beta = m["ANGL_NAUT"][1]
            gamma = m["ANGL_NAUT"][2]

        motscles = {}
        motscles["MODI_CHAM"] = []
        motscles["AFFE"] = []
        motscles["MODI_CHAM"].append(_F(NOM_CHAM=NOM_CHAM, TYPE_CHAM=TYPE_CHAM))
        ANGL_NAUT = []
        ANGL_NAUT.append(alpha)
        if dime == 3:
            ANGL_NAUT.append(beta)
            ANGL_NAUT.append(gamma)
        motscles["AFFE"].append(_F(ANGL_NAUT=ANGL_NAUT, TOUT="OUI"))
        __remodr = MODI_REPERE(RESULTAT=resin, REPERE="UTILISATEUR", **motscles)

    if m["TYPE"] == "ARC":
        if m["REPERE"] == "CYLINDRIQUE":
            motscles = {}
            motscles["MODI_CHAM"] = []
            motscles["AFFE"] = []
            motscles["MODI_CHAM"].append(_F(NOM_CHAM=NOM_CHAM, TYPE_CHAM=TYPE_CHAM))
            ORIGINE = []
            ORIGINE.append(m["CENTRE"][0])
            ORIGINE.append(m["CENTRE"][1])
            if dime == 3:
                ORIGINE.append(m["CENTRE"][2])
                AXE_Z = []
                AXE_Z.append(m["DNOR"][0])
                AXE_Z.append(m["DNOR"][1])
                AXE_Z.append(m["DNOR"][2])
                motscles["AFFE"].append(_F(ORIGINE=ORIGINE, AXE_Z=AXE_Z, TOUT="OUI"))
            elif dime == 2:
                motscles["AFFE"].append(_F(ORIGINE=ORIGINE, TOUT="OUI"))
            __remodr = MODI_REPERE(RESULTAT=resin, REPERE="CYLINDRIQUE", **motscles)
        else:
            UTMESS("F", "POST0_5", valk=[m["TYPE"], m["REPERE"]])

    if m["TYPE"][:5] == "GROUP" or m["TYPE"] == "SEGMENT":
        if m["TYPE"][:5] == "GROUP" and m["REPERE"] == "LOCAL":
            # determination du repère local (v1,v2,v3)
            # ---------------------------------------
            connex = mail.getConnectivity()
            coord = mail.getCoordinates()

            cells = mail.getCells(nomgrma)
            dictu = {}
            #     initialisations
            for cell in cells:
                n1 = connex[cell][0]
                n2 = connex[cell][1]
                dictu[n1] = []
                dictu[n2] = []
            #     determination du vecteur tangent (v1) + normalisation
            for cell in cells:
                vectu1 = []
                vectu2 = []
                n1 = connex[cell][0]
                n2 = connex[cell][1]
                ux = coord.getNode(n2).x() - coord.getNode(n1).x()
                uy = coord.getNode(n2).y() - coord.getNode(n1).y()
                vectu1.append(ux)
                vectu1.append(uy)
                vectu2.append(ux)
                vectu2.append(uy)
                if dime == 3:
                    uz = coord.getNode(n2).z() - coord.getNode(n2).z()
                    vectu1.append(uz)
                    vectu2.append(uz)
                dictu[n1].append(vectu1)
                dictu[n2].append(vectu2)
            for i in dictu:
                if len(dictu[i]) == 2:
                    dictu[i][0][0] = dictu[i][0][0] + dictu[i][1][0]
                    dictu[i][0][1] = dictu[i][0][1] + dictu[i][1][1]
                    if dime == 3:
                        dictu[i][0][2] = dictu[i][0][2] + dictu[i][1][2]
                    del dictu[i][1]
            for i in dictu:
                if dime == 2:
                    norm = sqrt(dictu[i][0][0] ** 2 + dictu[i][0][1] ** 2)
                    dictu[i][0][0] = dictu[i][0][0] / norm
                    dictu[i][0][1] = dictu[i][0][1] / norm
                elif dime == 3:
                    norm = sqrt(dictu[i][0][0] ** 2 + dictu[i][0][1] ** 2 + dictu[i][0][2] ** 2)
                    dictu[i][0][0] = dictu[i][0][0] / norm
                    dictu[i][0][1] = dictu[i][0][1] / norm
                    dictu[i][0][2] = dictu[i][0][2] / norm
            #     determination du vecteur normal (v2):
            #     on projete VECT_Y sur le plan orthogonal au vecteur v1.
            #     (ce vecteur normal est obtenu par 2 produits vectoriels successifs en 3D)
            if dime == 3:
                norm = sqrt(m["VECT_Y"][0] ** 2 + m["VECT_Y"][1] ** 2 + m["VECT_Y"][2] ** 2)
                tmpy = [m["VECT_Y"][0] / norm, m["VECT_Y"][1] / norm, m["VECT_Y"][2] / norm]
            j = 0
            __resu = [None] * (len(dictu) + 1)
            __resu[0] = resin
            for i in dictu:
                j = j + 1
                vecty = []
                if dime == 2:
                    vecty.append(-dictu[i][0][1])
                    vecty.append(dictu[i][0][0])
                    dictu[i].append(vecty)
                elif dime == 3:
                    # v3= v1 vectoriel vect_y
                    vectz = []
                    vectz.append(dictu[i][0][1] * tmpy[2] - dictu[i][0][2] * tmpy[1])
                    vectz.append(dictu[i][0][2] * tmpy[0] - dictu[i][0][0] * tmpy[2])
                    vectz.append(dictu[i][0][0] * tmpy[1] - dictu[i][0][1] * tmpy[0])
                    normz = sqrt(vectz[0] ** 2 + vectz[1] ** 2 + vectz[2] ** 2)
                    vectz[0] = vectz[0] / normz
                    vectz[1] = vectz[1] / normz
                    vectz[2] = vectz[2] / normz
                    vecty.append(vectz[1] * dictu[i][0][2] - vectz[2] * dictu[i][0][1])
                    vecty.append(vectz[2] * dictu[i][0][0] - vectz[0] * dictu[i][0][2])
                    vecty.append(vectz[0] * dictu[i][0][1] - vectz[1] * dictu[i][0][0])
                    normy = sqrt(vecty[0] ** 2 + vecty[1] ** 2 + vecty[2] ** 2)
                    vecty[0] = vecty[0] / normy
                    vecty[1] = vecty[1] / normy
                    vecty[2] = vecty[2] / normy
                    dictu[i].append(vecty)
                    dictu[i].append(vectz)
                cx1 = dictu[i][0][0]
                cx2 = dictu[i][0][1]
                cy1 = dictu[i][1][0]
                cy2 = dictu[i][1][1]
                if dime == 3:
                    cx3 = dictu[i][0][2]
                    cy3 = dictu[i][1][2]
                    cz1 = dictu[i][2][0]
                    cz2 = dictu[i][2][1]
                    cz3 = dictu[i][2][2]

                # determination des angles nautiques (alpha,beta,gamma)
                # ----------------------------------------------------
                beta = 0.0
                gamma = 0.0
                if dime == 2:
                    alpha = atan2(cx2, cx1)
                else:
                    if cx1**2 + cx2**2 > epsi:
                        alpha = atan2(cx2, cx1)
                        beta = -asin(cx3)
                        gamma = atan2(cy3, cz3)
                    else:
                        alpha = atan2(-cy1, cy2)
                        beta = -asin(cx3)
                        gamma = 0.0
                alpha = alpha * 180 / pi
                beta = beta * 180 / pi
                gamma = gamma * 180 / pi
                motscles = {}
                motscles["MODI_CHAM"] = []
                motscles["AFFE"] = []
                motscles["MODI_CHAM"].append(_F(NOM_CHAM=NOM_CHAM, TYPE_CHAM=TYPE_CHAM))
                ANGL_NAUT = []
                ANGL_NAUT.append(alpha)
                if dime == 3:
                    ANGL_NAUT.append(beta)
                    ANGL_NAUT.append(gamma)
                motscles["AFFE"].append(_F(ANGL_NAUT=ANGL_NAUT, NOEUD=mail.getNodeName(i)))
                __resu[j] = MODI_REPERE(RESULTAT=__resu[j - 1], REPERE="UTILISATEUR", **motscles)
            __remodr = __resu[j]

        motscles = {}
        motscles["MODI_CHAM"] = []
        motscles["AFFE"] = []
        motscles["MODI_CHAM"].append(_F(NOM_CHAM=NOM_CHAM, TYPE_CHAM=TYPE_CHAM))
        if m["REPERE"] == "CYLINDRIQUE":
            if dime == 3:
                motscles["AFFE"].append(_F(ORIGINE=m["ORIGINE"], AXE_Z=m["AXE_Z"], TOUT="OUI"))
            elif dime == 2:
                motscles["AFFE"].append(_F(ORIGINE=m["ORIGINE"], TOUT="OUI"))
            __remodr = MODI_REPERE(RESULTAT=resin, REPERE="CYLINDRIQUE", **motscles)
        elif m["REPERE"] == "UTILISATEUR":
            alpha = m["ANGL_NAUT"][0]
            beta = m["ANGL_NAUT"][1]
            gamma = m["ANGL_NAUT"][2]
            ANGL_NAUT = []
            ANGL_NAUT.append(alpha)
            if dime == 3:
                ANGL_NAUT.append(beta)
                ANGL_NAUT.append(gamma)
            motscles["AFFE"].append(_F(ANGL_NAUT=ANGL_NAUT, TOUT="OUI"))
            __remodr = MODI_REPERE(RESULTAT=resin, REPERE="UTILISATEUR", **motscles)

    return __remodr


#
# script PYTHON de creation des noeuds d'une ligne de coupe 'arc'


def crea_noeu_lig_coup(dimension, pt1, pt2, anglj, dnor):
    a = pt1[0] - pt2[0]
    b = pt1[1] - pt2[1]
    eps = 0.00000001
    anglr = anglj * pi / 180.0
    if dimension == 2:
        r = sqrt(a**2 + b**2)
        if abs(r) < eps:
            UTMESS("F", "POST0_6")
        x = pt2[0] + a * cos(anglr) - b * sin(anglr)
        y = pt2[1] + b * cos(anglr) + a * sin(anglr)
        return x, y
    elif dimension == 3:
        c = pt1[2] - pt2[2]
        r = sqrt(a**2 + b**2 + c**2)
        if abs(r) < eps:
            UTMESS("F", "POST0_6")
        d1 = dnor[0]
        d2 = dnor[1]
        d3 = dnor[2]
        d = sqrt(d1**2 + d2**2 + d3**2)
        if abs(r) < eps:
            UTMESS("F", "POST0_7")
        x = pt2[0] + a * cos(anglr) + sin(anglr) * (c * d2 - b * d3) / d
        y = pt2[1] + b * cos(anglr) + sin(anglr) * (a * d3 - c * d1) / d
        z = pt2[2] + c * cos(anglr) + sin(anglr) * (b * d1 - a * d2) / d
        return x, y, z


#
# determination de la distance min entre 2 points consécutifs de la ligne
# de coupe


def dist_min_deux_points(mail):
    nno = mail.getNumberOfNodes()
    coordinates = mail.getCoordinates()
    l_coor1 = []
    l_coor2 = []
    for i in range(nno - 1):
        l_coor1 = coordinates[i]
        l_coor2 = coordinates[i + 1]
        d = sqrt(
            (l_coor1[0] - l_coor2[0]) ** 2
            + (l_coor1[1] - l_coor2[1]) ** 2
            + (l_coor1[2] - l_coor2[2]) ** 2
        )
        if i == 0:
            dist = d
        else:
            dist = min(d, dist)
    return dist


#
# script PYTHON de creation d un maillage de ligne de coupe


def crea_mail_lig_coup(dimension, lignes, groups, arcs):
    # construction du maillage au format Aster des segments de lignes de coupe

    resu = []
    nblig = len(lignes)
    nbngr = len(groups)
    nbarc = len(arcs)

    resu.append("TITRE")
    resu.append("FINSF")
    resu.append("COOR_%dD" % dimension)
    epsi = 0.00000001

    # creation des noeuds
    nbno = 0
    for i in range(nblig):
        pt1 = lignes[i][0]
        pt2 = lignes[i][1]
        nbp_lig_coupe = lignes[i][2]
        for j in range(nbp_lig_coupe):
            if dimension == 2:
                x = pt1[0] + j * (pt2[0] - pt1[0]) / (nbp_lig_coupe - 1)
                y = pt1[1] + j * (pt2[1] - pt1[1]) / (nbp_lig_coupe - 1)
                nbno = nbno + 1
                noeud = "N%d %21.14E %21.14E" % (nbno, x, y)
                resu.append(noeud)
            elif dimension == 3:
                x = pt1[0] + j * (pt2[0] - pt1[0]) / (nbp_lig_coupe - 1)
                y = pt1[1] + j * (pt2[1] - pt1[1]) / (nbp_lig_coupe - 1)
                z = pt1[2] + j * (pt2[2] - pt1[2]) / (nbp_lig_coupe - 1)
                nbno = nbno + 1
                noeud = "N%d %21.14E %21.14E %21.14E" % (nbno, x, y, z)
                resu.append(noeud)
    for i in range(nbngr):
        for pt in groups[i][1:]:
            if dimension == 2:
                nbno = nbno + 1
                noeud = "N%d %21.14E %21.14E" % (nbno, pt[0], pt[1])
                resu.append(noeud)
            elif dimension == 3:
                nbno = nbno + 1
                noeud = "N%d %21.14E %21.14E %21.14E" % (nbno, pt[0], pt[1], pt[2])
                resu.append(noeud)
    angles = [None] * nbarc
    for i in range(nbarc):
        pt1 = arcs[i][0]
        pt2 = arcs[i][1]
        nbp_lig_coupe = arcs[i][2]
        angle = arcs[i][3]
        if abs(angle - 360.0) < epsi:
            nbpt = nbp_lig_coupe + 1
        else:
            nbpt = nbp_lig_coupe
        if dimension == 3:
            dnor = arcs[i][4]
        angles[i] = []
        for j in range(nbp_lig_coupe):
            anglj = j * angle / (nbpt - 1)
            angles[i].append(anglj)
            if dimension == 2:
                nbno = nbno + 1
                x, y = crea_noeu_lig_coup(dimension, pt1, pt2, anglj, dnor=[])
                noeud = "N%d %21.14E %21.14E" % (nbno, x, y)
                resu.append(noeud)
            elif dimension == 3:
                nbno = nbno + 1
                x, y, z = crea_noeu_lig_coup(dimension, pt1, pt2, anglj, dnor)
                noeud = "N%d %21.14E %21.14E %21.14E" % (nbno, x, y, z)
                resu.append(noeud)
    resu.append("FINSF")

    # creation des mailles
    nbma = 0
    for i in range(nblig):
        nbp_lig_coupe = lignes[i][2]
        resu.append("SEG2")
        for j in range(nbp_lig_coupe - 1):
            nbma = nbma + 1
            maille = "M%d N%d N%d" % (nbma, nbma + i, nbma + 1 + i)
            resu.append(maille)
        resu.append("FINSF")
    for i in range(nbngr):
        resu.append("SEG2")
        for pt in groups[i][1:-1]:
            nbma = nbma + 1
            maille = "M%d N%d N%d" % (nbma, nbma + nblig + i, nbma + nblig + 1 + i)
            resu.append(maille)
        resu.append("FINSF")
    nprec = 0

    for i in range(nbarc):
        nbp_lig_coupe = arcs[i][2]
        angle = arcs[i][3]
        resu.append("SEG2")
        nbmai = nbma + nblig + nbngr + nprec + i + 1
        for j in range(nbp_lig_coupe - 1):
            nbma = nbma + 1
            maille = "M%d N%d N%d" % (
                nbma,
                nbma + nblig + nbngr + nprec + i,
                nbma + nblig + nbngr + nprec + 1 + i,
            )
            resu.append(maille)
        if abs(angle - 360.0) < epsi:
            nbma = nbma + 1
            maille = "M%d N%d N%d" % (nbma, nbma + nblig + nbngr + nprec + i, nbmai)
            nprec = nprec - 1
            resu.append(maille)
        resu.append("FINSF")

    # creation des groupes de mailles (1 par ligne de coupe)
    nbma = 0
    for i in range(nblig):
        resu.append("GROUP_MA")
        resu.append("  LICOU%d" % (i + 1))
        nbp_lig_coupe = lignes[i][2]
        for j in range(nbp_lig_coupe - 1):
            nbma = nbma + 1
            resu.append("  M%d" % nbma)
        resu.append("")
        resu.append("FINSF")
    for i in range(nbngr):
        resu.append("GROUP_MA")
        resu.append(groups[i][0])
        nbp_lig_coupe = len(groups[i]) - 1
        for j in range(nbp_lig_coupe - 1):
            nbma = nbma + 1
            resu.append("  M%d" % nbma)
        resu.append("")
        resu.append("FINSF")
    arcgma = []
    for i in range(nbarc):
        resu.append("GROUP_MA")
        k = nblig + i
        resu.append("  LICOU%d" % (k + 1))
        arcgma.append("LICOU%d" % (k + 1))
        nbp_lig_coupe = arcs[i][2]
        angle = arcs[i][3]
        if abs(angle - 360.0) < epsi:
            nbpt = nbp_lig_coupe + 1
        else:
            nbpt = nbp_lig_coupe
        for j in range(nbpt - 1):
            nbma = nbma + 1
            resu.append("  M%d" % nbma)
        resu.append("")
        resu.append("FINSF")
    resu.append("FIN")
    resu.append("")

    return resu, arcgma, angles, nbno


def get_coor(LIGN_COUPE, position, coord, mesh):
    """
    Extrait les coordonnées du noeud ORIG ou EXTR à partir des coordonnées
    ou bien d'un groupe de noeuds ne contenant qu'un seul noeud.
    """
    assert position in ("ORIG", "EXTR")
    if "GROUP_NO_" + position in LIGN_COUPE:
        group = LIGN_COUPE["GROUP_NO_" + position]
        if not mesh.hasGroupOfNodes(group):
            UTMESS("F", "POST0_13", valk=[group, mesh.getName()])
        nodes = mesh.getNodes(group)
        if len(nodes) != 1:
            UTMESS("F", "POST0_27", valk=group, vali=len(nodes))
        node = nodes[0]
        coor = [coord.getNode(node).x(), coord.getNode(node).y(), coord.getNode(node).z()]
        # only COOR_xxxx can be used further
        LIGN_COUPE["COOR_" + position] = coor
    elif "COOR_" + position in LIGN_COUPE:
        coor = LIGN_COUPE["COOR_" + position]
    else:
        # Utilisation impossible d'après le catalogue
        assert False
    return coor


#
def macr_lign_coupe_ops(
    self, LIGN_COUPE, RESULTAT=None, CHAM_GD=None, NOM_CHAM=None, MODELE=None, **args
):
    """
    Ecriture de la macro MACR_LIGN_COUPE
    """

    # La valeur par défaut n'est pas dans le catalogue, sinon le mot-clé devient
    # obligatoire dans AsterStudy
    UNITE_MAILLAGE = args.get("UNITE_MAILLAGE")
    if not UNITE_MAILLAGE:
        logical_unit = LogicalUnitFile.new_free(access=FileAccess.New)
        UNITE_MAILLAGE = logical_unit.unit

    # On importe les definitions des commandes a utiliser dans la macro

    #
    MasquerAlarme("ALGORITH12_43")
    MasquerAlarme("MODELISA5_53")
    MasquerAlarme("MODELE1_58")
    MasquerAlarme("MODELE1_63")
    MasquerAlarme("MODELE1_64")

    mcORDR = {}

    model = MODELE
    mesh = None

    l_mode_meca_sans_modele = False

    if RESULTAT:
        if "NUME_ORDRE" in args:
            mcORDR["NUME_ORDRE"] = args["NUME_ORDRE"]
        elif "NUME_MODE" in args:
            mcORDR["NUME_MODE"] = args["NUME_MODE"]
        elif "LIST_ORDRE" in args:
            mcORDR["LIST_ORDRE"] = args["LIST_ORDRE"]
        elif "INST" in args:
            mcORDR["INST"] = args["INST"]
            mcORDR["CRITERE"] = args["CRITERE"]
            mcORDR["PRECISION"] = args["PRECISION"]
        elif "LIST_INST" in args:
            mcORDR["LIST_INST"] = args["LIST_INST"]
            mcORDR["CRITERE"] = args["CRITERE"]
            mcORDR["PRECISION"] = args["PRECISION"]
        else:
            mcORDR["TOUT_ORDRE"] = "OUI"

        type_resu = RESULTAT.getType().lower()
        model_resu = RESULTAT.getModel()

        if (not model_resu) or (model_resu.getName() in ("", "#AUCUN")):
            if not MODELE:
                if type_resu != "mode_meca":
                    UTMESS("F", "POST0_9", valk=RESULTAT.getName())
                # si le résultat en entrée est un mode_meca et qu'il ne contient pas de modèle (il est obtenu par sous-structuration, par exemple)
                # on passe le message fatal et on récupérera directement le
                # maillage (ou squelette)
                else:
                    l_mode_meca_sans_modele = True
                    mesh = RESULTAT.getMesh()
                    UTMESS("I", "POST0_23", valk=RESULTAT.getName())
            else:
                mesh = MODELE.getMesh()
        else:
            mesh = model_resu.getMesh()
            if not MODELE:
                model = model_resu

        if not mesh:
            raise Exception("Empty mesh")

    elif CHAM_GD:
        mcORDR["TOUT_ORDRE"] = "OUI"
        if not MODELE:
            UTMESS("F", "POST0_10")

        # récupération de la grandeur du champ
        nomgd = CHAM_GD.getPhysicalQuantity()

        # détermination du type de résultat à créer
        if nomgd[:6] == "TEMP_R":
            TYPE_RESU = "EVOL_THER"
            if not NOM_CHAM:
                NOM_CHAM = "TEMP"
        elif nomgd[:6] == "DEPL_R":
            TYPE_RESU = "EVOL_ELAS"
            if not NOM_CHAM:
                NOM_CHAM = "DEPL"
        elif nomgd[:6] == "NEUT_R":
            TYPE_RESU = "EVOL_VARC"
            if not NOM_CHAM:
                NOM_CHAM = "NEUT"
        elif nomgd[:6] == "EPSI_R":
            TYPE_RESU = "EVOL_ELAS"
        elif nomgd[:6] == "VAR2_R":
            TYPE_RESU = "EVOL_NOLI"
        elif nomgd[:6] == "VARI_R":
            TYPE_RESU = "EVOL_NOLI"
        elif nomgd[:6] == "SIEF_R":
            if not NOM_CHAM:
                TYPE_RESU = "EVOL_ELAS"
                NOM_CHAM = "DEPL"
            elif NOM_CHAM[:4] == "SIGM":
                TYPE_RESU = "EVOL_ELAS"
            elif NOM_CHAM[:4] == "SIEF":
                TYPE_RESU = "EVOL_NOLI"
        else:
            assert 0, "grandeur imprevue : " + nomgd

        # création d'un concept résultat à partir du champ CHAM_GD
        __resuch = CREA_RESU(
            OPERATION="AFFE",
            TYPE_RESU=TYPE_RESU,
            AFFE=_F(NOM_CHAM=NOM_CHAM, CHAM_GD=CHAM_GD, INST=0.0, MODELE=MODELE),
        )
        RESULTAT = __resuch

    # Maillage sur lequel s'appuie le résultat à projeter
    if not mesh:
        mesh = model.getMesh()
    # le maillage est-il 2D ou 3D ?
    dime = mesh.getDimension()
    coord = mesh.getCoordinates()
    lignes = []
    groups = []
    arcs = []
    minidim = dime

    for m in LIGN_COUPE:
        if m["TYPE"] == "SEGMENT":
            coor_orig = get_coor(m, "ORIG", coord, mesh)
            coor_extr = get_coor(m, "EXTR", coord, mesh)
            lignes.append((coor_orig, coor_extr, m["NB_POINTS"]))
            minidim = min(minidim, len(coor_orig), len(coor_extr))
            if minidim != dime:
                UTMESS("F", "POST0_11")

        elif m["TYPE"] == "ARC":
            minidim = min(minidim, len(m["COOR_ORIG"]), len(m["CENTRE"]))
            if minidim != dime:
                UTMESS("F", "POST0_11")
            if dime == 2:
                arcs.append((m["COOR_ORIG"], m["CENTRE"], m["NB_POINTS"], m["ANGLE"]))
            elif dime == 3:
                if str(m["DNOR"]) == "None":
                    UTMESS("F", "POST0_12")
                arcs.append((m["COOR_ORIG"], m["CENTRE"], m["NB_POINTS"], m["ANGLE"], m["DNOR"]))

        elif m["TYPE"] == "GROUP_NO":
            group = m["GROUP_NO"]
            if not mesh.hasGroupOfNodes(group):
                UTMESS("F", "POST0_13", valk=[group, mesh.getName()])
            l_coor_group = [group]
            for node in mesh.getNodes(group):
                l_coor_group.append(
                    [coord.getNode(node).x(), coord.getNode(node).y(), coord.getNode(node).z()]
                )
            groups.append(l_coor_group)

        elif m["TYPE"] == "GROUP_MA":
            group = m["GROUP_MA"]
            if not mesh.hasGroupOfCells(group):
                UTMESS("F", "POST0_14", valk=[group, mesh.getName()])
            for cell in mesh.getCells(group):
                if mesh.getCellTypeName(cell)[:3] != "SEG":
                    UTMESS("F", "POST0_15", valk=[group, mesh.getCellName(cell)])
            __mailla = COPIER(CONCEPT=m["MAILLAGE"])

            m2 = m.cree_dict_valeurs(m.mc_liste)
            argsup = {}
            if m2.get("GROUP_NO_ORIG"):
                argsup["GROUP_NO_ORIG"] = m2.get("GROUP_NO_ORIG")
            if m2.get("GROUP_NO_EXTR"):
                argsup["GROUP_NO_EXTR"] = m2.get("GROUP_NO_EXTR")
            if m2.get("VECT_ORIE"):
                argsup["VECT_ORIE"] = m2.get("VECT_ORIE")

            __mailla = DEFI_GROUP(
                reuse=__mailla,
                MAILLAGE=__mailla,
                DETR_GROUP_NO=_F(NOM=str(m["GROUP_MA"])),
                CREA_GROUP_NO=_F(
                    OPTION="NOEUD_ORDO",
                    NOM=str(m["GROUP_MA"]),
                    GROUP_MA=m["GROUP_MA"],
                    ORIGINE="SANS",
                    **argsup,
                ),
            )

            l_coor_group = [group]
            for node in __mailla.getNodes(group):
                l_coor_group.append(
                    [coord.getNode(node).x(), coord.getNode(node).y(), coord.getNode(node).z()]
                )
            groups.append(l_coor_group)

    if arcs and (lignes or groups):
        UTMESS("F", "POST0_16")

    # Création du maillage des NB_POINTS segments entre COOR_ORIG et COOR_EXTR
    # ainsi que des segments reliant les noeuds issus des group_no demandés
    # par appel au script python crea_mail_lig_coup
    # le maillage est ensuite recopié dans l unité logique UNITE_MAILLAGE

    resu_mail, arcgma, angles, nbno = crea_mail_lig_coup(dime, lignes, groups, arcs)

    nomFichierSortie = LogicalUnitFile.filename_from_unit(UNITE_MAILLAGE)
    with open(nomFichierSortie, "w") as fproc:
        fproc.write(os.linesep.join(resu_mail))

    # Lecture du maillage de seg2 contenant toutes les lignes de coupe
    __macou = LIRE_MAILLAGE(FORMAT="ASTER", UNITE=UNITE_MAILLAGE)
    LogicalUnitFile.release_from_number(UNITE_MAILLAGE)

    motscles = {}
    iocc = 1
    motscles["CREA_GROUP_NO"] = []
    for m in LIGN_COUPE:
        if m["TYPE"] in ("GROUP_NO", "GROUP_MA"):
            motscles["CREA_GROUP_NO"].append(_F(GROUP_MA=m[m["TYPE"]].ljust(24)))
        else:
            motscles["CREA_GROUP_NO"].append(_F(GROUP_MA="LICOU" + str(iocc)))
            iocc = iocc + 1

    __macou = DEFI_GROUP(reuse=__macou, MAILLAGE=__macou, **motscles)

    # issue27543 : l'utilisation d'un modèle n'est pas utile et elle
    # engendrait dans le cas thermique l'émission des alarmes MODELE1_3 et MODELE1_53
    # Dans le cas où il y aurait besoin de réintroduire les modèles, j'ai remplacé
    # la modélisation PLAN par COQUE (laissée en commentaire) ce qui permet également
    # de supprimer les alarmes.

    motscles = {}
    motscles.update(mcORDR)
    motscles["VIS_A_VIS"] = []
    if "VIS_A_VIS" in args:
        for v in args["VIS_A_VIS"]:
            if v["GROUP_MA_1"]:
                motscles["VIS_A_VIS"].append(_F(GROUP_MA_1=v["GROUP_MA_1"], TOUT_2="OUI"))
            elif v["MAILLE_1"]:
                motscles["VIS_A_VIS"].append(_F(MAILLE_1=v["MAILLE_1"], TOUT_2="OUI"))

    if NOM_CHAM[5:9] == "ELGA":
        UTMESS("A", "POST0_18", valk=[NOM_CHAM])

    if l_mode_meca_sans_modele:
        motscles["MAILLAGE_1"] = mesh
    else:
        motscles["MODELE_1"] = model

    # On récupère les minimum des DISTANCE_MAX et DISTANCE_ALARME, ces valeurs
    # serviront pour le PROJ_CHAMP pour toutes les coupes

    all_distance_max = [m.get("DISTANCE_MAX") for m in LIGN_COUPE if m.get("DISTANCE_MAX")]
    all_distance_alarme = [m.get("DISTANCE_ALARME") for m in LIGN_COUPE if m.get("DISTANCE_ALARME")]

    min_distance_max = min(all_distance_max, default=0.0)
    min_distance_alarme = min(all_distance_alarme, default=0.0)

    # On émet une alarme dès qu'au moins une des valeurs des DISTANCE_MAX
    # ou DISTANCE_ALARME est différente de la valeur retenue

    if not all(dist == min_distance_max for dist in all_distance_max):
        UTMESS("A", "POST0_53", valr=[min_distance_max])

    if not all(dist == min_distance_alarme for dist in all_distance_alarme):
        UTMESS("A", "POST0_54", valr=[min_distance_alarme])

    __recou = PROJ_CHAMP(
        METHODE="COLLOCATION",
        RESULTAT=RESULTAT,
        DISTANCE_MAX=min_distance_max,
        DISTANCE_ALARME=min_distance_alarme,
        # issue27543 #MODELE_2=__mocou,
        MAILLAGE_2=__macou,
        TYPE_CHAM="NOEU",
        NOM_CHAM=NOM_CHAM,
        **motscles,
    )

    __remodr = __recou
    ioc2 = 0
    mcACTION = []

    if RESULTAT.getType().lower() in (
        "evol_ther",
        "evol_sech",
        "evol_elas",
        "evol_noli",
        "mode_meca",
        "evol_varc",
        "comb_fourier",
        "mult_elas",
        "fourier_elas",
        "dyna_trans",
    ):

        for iocc, m in enumerate(LIGN_COUPE):
            motscles = {}
            motscles["OPERATION"] = m["OPERATION"]
            if m["NOM_CMP"]:
                motscles["NOM_CMP"] = m["NOM_CMP"]
                if m["TRAC_NOR"]:
                    motscles["TRAC_NOR"] = m["TRAC_NOR"]
                elif m["TRAC_DIR"]:
                    motscles["TRAC_DIR"] = m["TRAC_DIR"]
                    motscles["DIRECTION"] = m["DIRECTION"]
            elif m["INVARIANT"]:
                motscles["INVARIANT"] = m["INVARIANT"]
            elif m["RESULTANTE"]:
                motscles["RESULTANTE"] = m["RESULTANTE"]
                if m["MOMENT"]:
                    motscles["MOMENT"] = m["MOMENT"]
                    motscles["POINT"] = m["POINT"]
            elif m["ELEM_PRINCIPAUX"]:
                motscles["ELEM_PRINCIPAUX"] = m["ELEM_PRINCIPAUX"]
            else:
                motscles["TOUT_CMP"] = "OUI"

            if m["TYPE"] in ("GROUP_NO", "GROUP_MA"):
                groupe = m[m["TYPE"]].ljust(8)
                nomgrma = groupe
            else:
                ioc2 += 1
                nomgrma = " "
                groupe = f"LICOF{ioc2}"
                crea_grp_matiere(f"LICOU{ioc2}", groupe, iocc, m, __remodr, NOM_CHAM, __macou)

            # on definit l'intitulé
            if m["INTITULE"]:
                intitl = m["INTITULE"]
            elif m["TYPE"] in ("GROUP_NO", "GROUP_MA"):
                intitl = groupe
            else:
                intitl = f"l.coupe{ioc2}"

            # Expression des contraintes aux noeuds ou des déplacements dans le
            # repere local
            if m["REPERE"] != "GLOBAL":
                if NOM_CHAM in (
                    "DEPL",
                    "SIEF_ELNO",
                    "SIGM_NOEU",
                    "SIGM_ELNO",
                    "FLUX_ELNO",
                    "FLUX_NOEU",
                ):
                    if m["REPERE"] == "POLAIRE":
                        mcACTION.append(
                            _F(
                                INTITULE=intitl,
                                RESULTAT=__remodr,
                                REPERE=m["REPERE"],
                                GROUP_NO=groupe,
                                NOM_CHAM=NOM_CHAM,
                                **motscles,
                            )
                        )
                    else:
                        __remodr = crea_resu_local(dime, NOM_CHAM, m, __recou, __macou, nomgrma)
                        mcACTION.append(
                            _F(
                                INTITULE=intitl,
                                RESULTAT=__remodr,
                                GROUP_NO=groupe,
                                NOM_CHAM=NOM_CHAM,
                                **motscles,
                            )
                        )

                else:
                    UTMESS("A", "POST0_17", valk=[NOM_CHAM, m["REPERE"]])
                    mcACTION.append(
                        _F(
                            INTITULE=intitl,
                            RESULTAT=__recou,
                            GROUP_NO=groupe,
                            NOM_CHAM=NOM_CHAM,
                            **motscles,
                        )
                    )

            # Expression des contraintes aux noeuds ou des déplacements dans le
            # repere global
            else:
                mcACTION.append(
                    _F(
                        INTITULE=intitl,
                        RESULTAT=__recou,
                        GROUP_NO=groupe,
                        NOM_CHAM=NOM_CHAM,
                        **motscles,
                    )
                )

    else:
        assert 0

    __tabitm = POST_RELEVE_T(ACTION=mcACTION)

    # on repasse par les tables python pour supprimer les paramètres inutiles
    # NOEUD (car il est propre au maillage de la ligne) et RESU

    dictab = __tabitm.EXTR_TABLE()
    # Ajout de la colonne theta
    if len(arcgma) > 0 and "ABSC_CURV" in dictab.para:
        val = dictab["ABSC_CURV"].values()["ABSC_CURV"]
        nbi = len(val) // nbno
        nba = len(angles)
        tmp = []
        for k in range(nba):
            for j in range(nbi):
                for i in range(len(angles[k])):
                    tmp.append(angles[k][i])
        dictab["ANGLE"] = tmp

    if "RESU" in dictab.para:
        del dictab["RESU"]
    if "NOEUD" in dictab.para:
        del dictab["NOEUD"]
    dprod = dictab.dict_CREA_TABLE()
    # Remove TITRE key from the dictionary so that the newly created table uses
    # the user-given concept name in its title
    dprod.pop("TITRE")

    nomres = CREA_TABLE(**dprod)

    RetablirAlarme("ALGORITH12_43")
    RetablirAlarme("MODELISA5_53")
    RetablirAlarme("MODELE1_58")
    RetablirAlarme("MODELE1_63")
    RetablirAlarme("MODELE1_64")
    return nomres
