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

from ..Messages import UTMESS

from ..Cata.Syntax import _F
from ..CodeCommands import DEFI_GROUP
from .Utils.defi_cable_op import DEFI_CABLE_OP


# ===========================================================================
#           CORPS DE LA MACRO "DEFI_CABLE_BP"
#           -------------------------------------
# USAGE :
# Entrée :
#  - MODELE
#  - CABLE
#  - CHAM_MATER
#  - CARA_ELEM
#  - GROUP_MA_BETON
#  - DEFI_CABLE
#  - MODI_CABLE_ETCC
#  - MODI_CABLE_RUPT
#  - ADHERENT
#  - TYPE_ANCRAGE
#  - TENSION_INIT
#  - RECUL_ANCRAGE
#  - RELAXATION
#  - CONE
#      RAYON
#      LONGUEUR
#      PRESENT          OUI ou NON deux fois
#  - TITRE
#  - INFO               1 / 2
#
# ===========================================================================
def verif_table(tabin):
    """verification que l'on a le bon nombre de colonnes et les bons noms"""

    nbcab = len(tabin)
    if nbcab < 1:
        UTMESS("F", "CABLE0_22")
    for nom_para in ("GROUP_MA", "GROUP_NO1", "GROUP_NO2"):
        if nom_para not in tabin.para:
            UTMESS("F", "CABLE0_23")
    UTMESS("I", "CABLE0_24", vali=nbcab)
    return nbcab


def defi_cable_bp_ops(
    self,
    MODELE,
    CHAM_MATER,
    CARA_ELEM,
    GROUP_MA_BETON,
    ADHERENT,
    TYPE_ANCRAGE,
    TYPE_RELAX,
    INFO,
    CONE=None,
    **args
):
    """
    Ecriture de la macro DEFI_CABLE_BP
    """
    motscles = {}
    args = _F(args)

    # RECUPERATION DES INFOS DONNEES PAR LE MOT-CLE "CONE"

    if CONE:
        dCONE = CONE[0].cree_dict_valeurs(CONE[0].mc_liste)
        for i in list(dCONE.keys()):
            if dCONE[i] is None:
                del dCONE[i]

        RAYON = dCONE["RAYON"]
        LONGUEUR = dCONE["LONGUEUR"]

        motscles["CONE"] = []
        motscles["CONE"].append(dCONE)

        # RECUPERATION DU MAILLAGE A PARTIR DU MODELE
        MAILLAGE = MODELE.getMesh()

        # DEFINITION DU NOM DES GROUP_NO
        __NOM = "AN__"
        __LGNO = MAILLAGE.LIST_GROUP_NO()
        __LGN1 = []
        for i in __LGNO:
            __LGN1.append(i[0][: len(__NOM)])

        __NB = __LGN1.count(__NOM)

    # FIN RECUPERATION DES INFOS DONNEES PAR LE MOT-CLE "CONE"

    # RECUPERATION DES INFOS DONNEES PAR LE MOT-CLE "DEFI_CABLE"
    dDEFI_CABLE = []
    if args["DEFI_CABLE"] is not None:
        ANALY = "DEFI"
        for j in args["DEFI_CABLE"]:
            dDEFI_CABLE.append(j)
            for i in list(dDEFI_CABLE[-1].keys()):
                if dDEFI_CABLE[-1][i] is None:
                    del dDEFI_CABLE[-1][i]
            if "TABL_CABLE" in j:
                __TAB = j["TABL_CABLE"].EXTR_TABLE()
                nbcab = verif_table(__TAB)
                for ic in range(nbcab):
                    __gma = __TAB.GROUP_MA.values()[ic]
                    __gno1 = __TAB.GROUP_NO1.values()[ic]
                    __gno2 = __TAB.GROUP_NO2.values()[ic]
                    #              print(__gma,__gno1,__gno2)
                    dDEFI_CABLE.append(_F(GROUP_MA=__gma, GROUP_NO_ANCRAGE=(__gno1, __gno2)))

    elif args["MODI_CABLE_ETCC"] is not None:
        ANALY = "ETCC"
        # RECUPERATION DES INFOS DONNEES PAR LE MOT-CLE "MODI_CABLE_ETCC"
        for j in args["MODI_CABLE_ETCC"]:
            dDEFI_CABLE.append(j)
            for i in list(dDEFI_CABLE[-1].keys()):
                if dDEFI_CABLE[-1][i] is None:
                    del dDEFI_CABLE[-1][i]

    # RECUPERATION DES INFOS DONNEES PAR LE MOT-CLE "MODI_CABLE_RUPT"
    elif args["MODI_CABLE_RUPT"] is not None:
        ANALY = "RUPT"
        for j in args["MODI_CABLE_RUPT"]:
            dDEFI_CABLE.append(j)
            for i in list(dDEFI_CABLE[-1].keys()):
                if dDEFI_CABLE[-1][i] is None:
                    del dDEFI_CABLE[-1][i]

    else:
        UTMESS("F", "CABLE0_25")

    # BOUCLE SUR LES FACTEURS DU MOT-CLE "DEFI_CABLE" et "MODI_CABLE"
    motscles["DEFI_CABLE"] = []

    for i in dDEFI_CABLE:
        #   CAS OU ON RENTRE UNE TENSION INITIALE DU CABLE (TYPE_RELAX='ETCC_REPRISE')
        motscle3 = {}
        if ("TENSION" in i) == 1:
            motscle3 = {"TENSION": i["TENSION"]}

        # CAS OU L'ON A DEFINI LE MOT-CLE "CONE"
        if CONE:
            # CREATION DU PREMIER TUNNEL

            if dCONE["PRESENT"][0] == "OUI":
                __NB = __NB + 1
                __NOM1 = __NOM + str(int(__NB))

                motscle2 = {}
                motscle2["CREA_GROUP_NO"] = []

                if ("GROUP_MA" in i) == 1:
                    __CAB = i["GROUP_MA"]

                    if type(GROUP_MA_BETON) in [tuple, list]:
                        gma = list(GROUP_MA_BETON)
                    else:
                        gma = [GROUP_MA_BETON]
                    gma.insert(0, __CAB)

                    motscle2 = {
                        "CREA_GROUP_NO": [
                            {
                                "LONGUEUR": LONGUEUR,
                                "RAYON": RAYON,
                                "OPTION": "TUNNEL",
                                "GROUP_MA": gma,
                                "GROUP_MA_AXE": __CAB,
                                "NOM": __NOM1,
                            }
                        ]
                    }
                if ("MAILLE" in i) == 1:
                    UTMESS("F", "CABLE0_2")
                if ("GROUP_NO_ANCRAGE" in i) == 1:
                    __PC1 = i["GROUP_NO_ANCRAGE"][0]
                    motscle2["CREA_GROUP_NO"][0]["GROUP_NO_ORIG"] = __PC1

                DEFI_GROUP(reuse=MAILLAGE, MAILLAGE=MAILLAGE, INFO=INFO, ALARME="NON", **motscle2)

            # CREATION DU DEUXIEME TUNNEL

            if dCONE["PRESENT"][1] == "OUI":
                __NB = __NB + 1
                __NOM2 = __NOM + str(int(__NB))

                motscle2 = {}
                motscle2["CREA_GROUP_NO"] = []

                if ("GROUP_MA" in i) == 1:
                    __CAB = i["GROUP_MA"]

                    if type(GROUP_MA_BETON) in [tuple, list]:
                        gma = list(GROUP_MA_BETON)
                    else:
                        gma = [GROUP_MA_BETON]
                    gma.insert(0, __CAB)

                    motscle2 = {
                        "CREA_GROUP_NO": [
                            {
                                "LONGUEUR": LONGUEUR,
                                "RAYON": RAYON,
                                "OPTION": "TUNNEL",
                                "GROUP_MA": gma,
                                "GROUP_MA_AXE": __CAB,
                                "NOM": __NOM2,
                            }
                        ]
                    }
                if ("MAILLE" in i) == 1:
                    UTMESS("F", "CABLE0_2")
                if ("GROUP_NO_ANCRAGE" in i) == 1:
                    __PC1 = i["GROUP_NO_ANCRAGE"][1]
                    motscle2["CREA_GROUP_NO"][0]["GROUP_NO_ORIG"] = __PC1

                DEFI_GROUP(reuse=MAILLAGE, MAILLAGE=MAILLAGE, INFO=INFO, ALARME="NON", **motscle2)

            # CREATION DES NOUVEAUX FACTEURS DU MOT-CLE "DEFI_CABLE" POUR
            # DEFI_CABLE_BP
            if dCONE["PRESENT"][0] == "OUI" and dCONE["PRESENT"][1] == "OUI":
                if ("GROUP_MA" in i) == 1 and ("GROUP_NO_ANCRAGE" in i) == 1:
                    motscles["DEFI_CABLE"].append(
                        _F(
                            GROUP_MA=i["GROUP_MA"],
                            GROUP_NO_ANCRAGE=i["GROUP_NO_ANCRAGE"],
                            GROUP_NO_FUT=(__NOM1, __NOM2),
                            **motscle3
                        )
                    )

            if dCONE["PRESENT"][0] == "OUI" and dCONE["PRESENT"][1] == "NON":
                if ("GROUP_MA" in i) == 1 and ("GROUP_NO_ANCRAGE" in i) == 1:
                    motscles["DEFI_CABLE"].append(
                        _F(
                            GROUP_MA=i["GROUP_MA"],
                            GROUP_NO_ANCRAGE=i["GROUP_NO_ANCRAGE"],
                            GROUP_NO_FUT=(__NOM1,),
                            **motscle3
                        )
                    )

            if dCONE["PRESENT"][0] == "NON" and dCONE["PRESENT"][1] == "OUI":
                if ("GROUP_MA" in i) == 1 and ("GROUP_NO_ANCRAGE" in i) == 1:
                    motscles["DEFI_CABLE"].append(
                        _F(
                            GROUP_MA=i["GROUP_MA"],
                            GROUP_NO_ANCRAGE=i["GROUP_NO_ANCRAGE"],
                            GROUP_NO_FUT=(__NOM2,),
                            **motscle3
                        )
                    )

            if dCONE["PRESENT"][0] == "NON" and dCONE["PRESENT"][1] == "NON":
                if ("GROUP_MA" in i) == 1 and ("GROUP_NO_ANCRAGE" in i) == 1:
                    motscles["DEFI_CABLE"].append(
                        _F(
                            GROUP_MA=i["GROUP_MA"],
                            GROUP_NO_ANCRAGE=i["GROUP_NO_ANCRAGE"],
                            **motscle3
                        )
                    )

        # CAS OU L'ON A PAS DEFINI LE MOT-CLE "CONE"
        else:
            if ("GROUP_MA" in i) == 1 and ("GROUP_NO_ANCRAGE" in i) == 1:
                motscles["DEFI_CABLE"].append(
                    _F(GROUP_MA=i["GROUP_MA"], GROUP_NO_ANCRAGE=i["GROUP_NO_ANCRAGE"], **motscle3)
                )

            if ("MAILLE" in i) == 1 and ("GROUP_NO_ANCRAGE" in i) == 1:
                motscles["DEFI_CABLE"].append(
                    _F(MAILLE=i["MAILLE"], GROUP_NO_ANCRAGE=i["GROUP_NO_ANCRAGE"], **motscle3)
                )

    # FIN BOUCLE sur i in DEFI_CABLE
    # LANCEMENT DE DEFI_CABLE_BP
    #    TRAITEMENT DE LA RELAXATION
    if TYPE_RELAX == "ETCC_DIRECT" or ANALY == "ETCC":
        motscles["NBH_RELAX"] = args["NBH_RELAX"]

    if TYPE_RELAX == "BPEL":
        motscles["R_J"] = args["R_J"]

    #  Traitement mot-cle facultatif
    tit = args.get("TITRE")
    if tit is not None:
        motscles["TITRE"] = tit
    #  if PERT_ELAS=='OUI':
    #    motscles['ESP_CABLE']=args['ESP_CABLE'] ;
    #    motscles['EP_BETON']=args['EP_BETON'] ;

    __DC = DEFI_CABLE_OP(
        MODELE=MODELE,
        CHAM_MATER=CHAM_MATER,
        CARA_ELEM=CARA_ELEM,
        GROUP_MA_BETON=GROUP_MA_BETON,
        ADHERENT=ADHERENT,
        TYPE_ANCRAGE=TYPE_ANCRAGE,
        TENSION_INIT=args.get("TENSION_INIT", None),
        RECUL_ANCRAGE=args.get("RECUL_ANCRAGE", None),
        TYPE_RELAX=TYPE_RELAX,
        ANALYSE=ANALY,
        INFO=INFO,
        **motscles
    )

    return __DC
