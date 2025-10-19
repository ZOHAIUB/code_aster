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

# from numpy import array, add, sqrt, dot, ones, linspace, outer, eye, reshape, sin, cos, log, zeros
from numpy import *
from numpy.linalg import norm
from scipy.integrate import odeint

MU = 80000.0  # MPa
NU = 0.3
E = MU * 2.0 * (1.0 + NU)
k = E / (3.0 * (1.0 - 2.0 * NU))
TAU_F = 105.0
RHO_0 = 1.0e6  # en mm**-2
RHOREF = 1.0e6  # en mm**-2
N = 5.0
GAMMA0 = 1.0e-3
ALPHA = 0.35
BETA = 2.54e-7  # 2.54 Angstrom
A = 0.13
B = 0.005
Y = 2.5e-7  # 2.5 Angstrom

rho = 12 * [RHO_0]
rho = array(rho)

list_cfc = [
    ["B4", ((1, 1, 1), (-1, 0, 1))],
    ["B2", ((1, 1, 1), (0, -1, 1))],
    ["B5", ((1, 1, 1), (-1, 1, 0))],
    ["D4", ((1, -1, 1), (-1, 0, 1))],
    ["D1", ((1, -1, 1), (0, 1, 1))],
    ["D6", ((1, -1, 1), (1, 1, 0))],
    ["A2", ((-1, 1, 1), (0, -1, 1))],
    ["A6", ((-1, 1, 1), (1, 1, 0))],
    ["A3", ((-1, 1, 1), (1, 0, 1))],
    ["C5", ((-1, -1, 1), (-1, 1, 0))],
    ["C3", ((-1, -1, 1), (1, 0, 1))],
    ["C1", ((-1, -1, 1), (0, 1, 1))],
]

HSB = array(
    [
        [0.124, 0.124, 0.124, 0.625, 0.137, 0.137, 0.070, 0.137, 0.122, 0.070, 0.122, 0.137],
        [0.124, 0.124, 0.124, 0.137, 0.070, 0.122, 0.137, 0.625, 0.137, 0.122, 0.070, 0.137],
        [0.124, 0.124, 0.124, 0.137, 0.122, 0.070, 0.122, 0.137, 0.070, 0.137, 0.137, 0.625],
        [0.625, 0.137, 0.137, 0.124, 0.124, 0.124, 0.070, 0.122, 0.137, 0.070, 0.137, 0.122],
        [0.137, 0.070, 0.122, 0.124, 0.124, 0.124, 0.122, 0.070, 0.137, 0.137, 0.625, 0.137],
        [0.137, 0.122, 0.070, 0.124, 0.124, 0.124, 0.137, 0.137, 0.625, 0.122, 0.137, 0.070],
        [0.070, 0.137, 0.122, 0.070, 0.122, 0.137, 0.124, 0.124, 0.124, 0.625, 0.137, 0.137],
        [0.137, 0.625, 0.137, 0.122, 0.070, 0.137, 0.124, 0.124, 0.124, 0.137, 0.070, 0.122],
        [0.122, 0.137, 0.070, 0.137, 0.137, 0.625, 0.124, 0.124, 0.124, 0.137, 0.122, 0.070],
        [0.070, 0.122, 0.137, 0.070, 0.137, 0.122, 0.625, 0.137, 0.137, 0.124, 0.124, 0.124],
        [0.122, 0.070, 0.137, 0.137, 0.625, 0.137, 0.137, 0.070, 0.122, 0.124, 0.124, 0.124],
        [0.137, 0.137, 0.625, 0.122, 0.137, 0.070, 0.137, 0.122, 0.070, 0.124, 0.124, 0.124],
    ]
)
NN = [6, 8, 7, 1, 0, 2, 11, 10, 9, 4, 3, 5]

H = zeros(144)
H.resize(12, 12)

for i in range(12):
    for j in range(12):
        H[NN[i], NN[j]] = HSB[i, j]

print("H=", H)

l0 = array([0, 0, 1], dtype=float)

# tenseurs d'orientation pour les systemes de glissements
tens = []
tens_mu = []

for i in range(0, len(list_cfc)):
    tens.append(
        outer(
            list_cfc[i][1][1] / norm(list_cfc[i][1][1], 2),
            list_cfc[i][1][0] / norm(list_cfc[i][1][0], 2),
        )
    )
    tens_mu.append(1.0 / 2.0 * (tens[i] + tens[i].transpose()))

tens_mu = array(tens_mu)

vect_mu = array(
    [
        tens_mu[:, 0, 0],
        tens_mu[:, 1, 1],
        tens_mu[:, 2, 2],
        sqrt(2) * tens_mu[:, 0, 1],
        sqrt(2) * tens_mu[:, 0, 2],
        sqrt(2) * tens_mu[:, 1, 2],
    ]
).transpose()


def facteur_schmid(l=l0, list_orient=list_cfc):
    """
    Cette fonction calcul les facteurs de schmid a partir d'une direction de sollicitation l et de la liste des systemes de glissement nommes suivant
    la convention de Schmid et Boas (normal puis direction de glissement)
    Elle retourne une liste contenant les facteurs de Schmid pour un systeme donne (ainsi qu'une valeur qui multipliee par sqrt(6) donne le denominateur et le numerateur)
    """
    schmid = []
    for i in range(0, len(list_orient)):
        schmid.append([])
        schmid[i].append(list_orient[i][0])
        schmid[i].append(
            [
                dot(list_orient[i][1][0] / norm(list_orient[i][1][0], 2), l / norm(l, 2))
                * dot(list_orient[i][1][1] / norm(list_orient[i][1][1], 2), l / norm(l, 2)),
                norm(l, 2) ** 2,
                dot(list_orient[i][1][0], l) * dot(list_orient[i][1][1], l),
            ]
        )
        # print schmid[i]
    return schmid


schmid = facteur_schmid()
schmid = array([schmid[i][1][0] for i in range(0, len(schmid))])

print("schmid", schmid)

# systemes coplanaires et forets pour l'evolution de rho
copla = array(
    [
        [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
    ]
)

forest = ones(12) - copla

sigma = 550.0
contr = sigma * outer(l0 / norm(l0, 2), l0 / norm(l0, 2))
vect_contr = array(
    [
        contr[0, 0],
        contr[1, 1],
        contr[2, 2],
        sqrt(2) * contr[0, 1],
        sqrt(2) * contr[0, 2],
        sqrt(2) * contr[1, 2],
    ]
)


def Ceff(Srho):
    return 0.2 + 0.8 * log(ALPHA * BETA * sqrt(Srho)) / log(ALPHA * BETA * sqrt(RHOREF))


# array
def taumu(rho):
    Srho = add.reduce(rho)
    Hrho = dot(H, rho)
    #     SHrho=add.reduce(Hrho)
    #     return MU*BETA*sqrt(SHrho)*Ceff(Srho)
    return MU * BETA * sqrt(Hrho) * Ceff(Srho)


# array
def gammap(rho, t):
    test = (schmid * sigma * t) // (TAU_F + taumu(rho))
    # print t, test
    return GAMMA0 * ((schmid * sigma * t / (TAU_F + taumu(rho))) ** N - 1.0) * test


# array
def rhopoint(rho, t):
    Srho = add.reduce(rho)
    # Hrho=dot(H,rho)
    Hrho = zeros(12)
    Hfrho = zeros(12)
    Hcrho = zeros(12)
    vfor = zeros(12)
    for i in range(12):
        for j in range(12):
            Hrho[i] += sqrt(H[i, j] * rho[j])
            Hfrho[i] += sqrt(H[i, j]) * rho[j] * forest[i, j]
            Hcrho[i] += sqrt(H[i, j] * rho[j]) * copla[i, j]
    tmp = A * Hfrho / Hrho + B * Ceff(Srho) * Hcrho - Y * rho
    return tmp * gammap(rho, t) / BETA


discretis = 100000
t = linspace(0.0, 1.0, discretis)

# pas robuste si N > 5 !
print("attention solution par odeint valable pour n<5")

rho_sol = odeint(rhopoint, rho, t)

rho_solb2 = rho_sol * BETA**2
print("rho final", rho_solb2[-1])
