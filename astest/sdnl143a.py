# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

import numpy

##################################################################################################################################################################################
# CONDITIONS DE MER
H = 4.0
T = 10.0
proinf = 0
depth = 20.0
thetai = 0.0
lstrech = 0


def hunt(T, DEPTH, PROINF):
    PI = numpy.pi
    GRAVIT = 9.81
    OMEGA = 2.0 * PI / T
    if PROINF == 1:
        L = (GRAVIT / 2.0 * PI) * T**2
    else:
        #        HUNT (1979) C=C=C=C=C=C=C=C=C
        T2 = OMEGA**2 * DEPTH / GRAVIT
        T1 = 1.0 + T2 * (
            0.66667
            + T2
            * (
                0.3555
                + T2
                * (
                    0.16084
                    + T2
                    * (
                        0.0632
                        + T2
                        * (
                            0.02174
                            + T2 * (0.00654 + T2 * (0.00171 + T2 * (0.00039 + T2 * 0.00011)))
                        )
                    )
                )
            )
        )
        T1 = numpy.sqrt(T2**2 + (T2 / T1))
        KK = T1 / DEPTH
        L = 2.0 * PI / KK
    return L


L = hunt(T, depth, proinf)

# CINEMATIQUES ET FORCES
def cinelin(XX, YY, ZZ, TEMPS, H, T, L, THETA, UC, THETAC, PROINF, DEPTH, LSTRECH):
    PI = numpy.pi
    GRAVIT = 9.81
    OMEGA = 2.0 * PI / T
    AA = H / 2.0
    THETAR = THETA * PI / 180.0
    THETACR = THETAC * PI / 180.0

    KK = 2.0 * PI / L
    KX = KK * numpy.cos(THETAR)
    KY = KK * numpy.sin(THETAR)

    OMEGAR = OMEGA - KK * UC * numpy.cos(THETAR - THETACR)

    if PROINF == 1:
        ETA = AA * numpy.cos(KX * XX + KY * YY - OMEGA * TEMPS)
        #               -- MODIFICATION DU ZERO DE REFERENCE POUR LE STRECHING DE WHEELER
        if LSTRECH == 1:
            ETAJ = ETA
        else:
            ETAJ = 0.0
            #
            EKZ = numpy.exp(KK * (ZZ - ETAJ))
            UX = (AA * GRAVIT * KX / OMEGAR) * numpy.cos(KX * XX + KY * YY - OMEGA * TEMPS) * EKZ
            UY = (AA * GRAVIT * KY / OMEGAR) * numpy.cos(KX * XX + KY * YY - OMEGA * TEMPS) * EKZ
            UZ = (AA * GRAVIT * KK / OMEGAR) * numpy.sin(KX * XX + KY * YY - OMEGA * TEMPS) * EKZ
            #
            AX = (
                (AA * GRAVIT * KX)
                * (OMEGA / OMEGAR)
                * numpy.sin(KX * XX + KY * YY - OMEGA * TEMPS)
                * EKZ
            )
            AY = (
                (AA * GRAVIT * KY)
                * (OMEGA / OMEGAR)
                * numpy.sin(KX * XX + KY * YY - OMEGA * TEMPS)
                * EKZ
            )
            AZ = (
                -(AA * GRAVIT * KK)
                * (OMEGA / OMEGAR)
                * numpy.cos(KX * XX + KY * YY - OMEGA * TEMPS)
                * EKZ
            )
    else:
        ETA = AA * numpy.cos(KX * XX + KY * YY - OMEGA * TEMPS)
        #               -- MODIFICATION DU ZERO DE REFERENCE POUR LE STRECHING DE WHEELER
        if LSTRECH == 1:
            ETAJ = ETA
        else:
            ETAJ = 0.0
        CHKZ = numpy.cosh(KK * (ZZ + DEPTH) * DEPTH / (ETAJ + DEPTH))
        SHKZ = numpy.sinh(KK * (ZZ + DEPTH) * DEPTH / (ETAJ + DEPTH))
        CHKD = numpy.cosh(KK * DEPTH)

        UX = (
            (AA * GRAVIT * KX / OMEGAR) * numpy.cos(KX * XX + KY * YY - OMEGA * TEMPS) * CHKZ / CHKD
        )
        UY = (
            (AA * GRAVIT * KY / OMEGAR) * numpy.cos(KX * XX + KY * YY - OMEGA * TEMPS) * CHKZ / CHKD
        )
        UZ = (
            (AA * GRAVIT * KK / OMEGAR) * numpy.sin(KX * XX + KY * YY - OMEGA * TEMPS) * SHKZ / CHKD
        )

        AX = (
            (AA * GRAVIT * KX)
            * (OMEGA / OMEGAR)
            * numpy.sin(KX * XX + KY * YY - OMEGA * TEMPS)
            * CHKZ
            / CHKD
        )
        AY = (
            (AA * GRAVIT * KY)
            * (OMEGA / OMEGAR)
            * numpy.sin(KX * XX + KY * YY - OMEGA * TEMPS)
            * CHKZ
            / CHKD
        )
        AZ = (
            -(AA * GRAVIT * KK)
            * (OMEGA / OMEGAR)
            * numpy.cos(KX * XX + KY * YY - OMEGA * TEMPS)
            * SHKZ
            / CHKD
        )

    #   Prise en compte du courant
    UX = UX + UC * numpy.cos(THETACR)
    UY = UY + UC * numpy.sin(THETACR)
    return ETA, UX, UY, UZ, AX, AY, AZ


def projvect(VX, VY, VZ, EX, EY, EZ):
    VTX = VX - (VX * EX + VY * EY + VZ * EZ) * EX
    VTY = VY - (VX * EX + VY * EY + VZ * EZ) * EY
    VTZ = VZ - (VX * EX + VY * EY + VZ * EZ) * EZ
    return VTX, VTY, VTZ


def dragquad(VTX, VTY, VTZ, VT, D, RHO, CD):
    FXR = 0.5 * RHO * CD * D * VT * VTX
    FYR = 0.5 * RHO * CD * D * VT * VTY
    FZR = 0.5 * RHO * CD * D * VT * VTZ
    return FXR, FYR, FZR


def FP(
    rho,
    cd,
    ca,
    D,
    xx,
    yy,
    zz,
    tps,
    H,
    T,
    L,
    THETA,
    UC,
    THETAC,
    PROINF,
    dd,
    LSTRECH,
    vite_x,
    vite_y,
    vite_z,
    dx,
    dy,
    dz,
    dir,
    cinelin,
    projvect,
    dragquad,
):
    #   Cinematique de vague
    eta, ux, uy, uz, ax, ay, az = cinelin(
        xx, yy, zz, tps, H, T, L, THETA, UC, THETAC, PROINF, dd, LSTRECH
    )
    #   verif si dans l'eau et pas dans le sol !!
    if (eta > zz) and (zz > -dd):
        norme = numpy.sqrt(dx**2 + dy**2 + dz**2)
        ex = dx / norme
        ey = dy / norme
        ez = dz / norme
        # projection vitesse relative ux-vite_x,uy-vite_y,uz-vite_z |_ dx,dy,dz
        vtx, vty, vtz = projvect(ux - vite_x, uy - vite_y, uz - vite_z, ex, ey, ez)
        # projection ax,ay,az |_ dx,dy,dz
        atx, aty, atz = projvect(ax, ay, az, ex, ey, ez)
        # calcul 1/2 rho cd D V abs(V) sur 3 direction
        vt = numpy.sqrt(vtx**2 + vty**2 + vtz**2)
        FX, FY, FZ = dragquad(vtx, vty, vtz, vt, D, rho, cd)
        #       Direction
        if dir == 1:
            F = FX + numpy.pi * D**2 / 4.0 * (1.0 + ca) * rho * atx
        elif dir == 2:
            F = FY + numpy.pi * D**2 / 4.0 * (1.0 + ca) * rho * aty
        elif dir == 3:
            F = FZ + numpy.pi * D**2 / 4.0 * (1.0 + ca) * rho * atz
        else:
            F = 0
    else:
        F = 0.0
    return F


##################################################################################################################################################################################
