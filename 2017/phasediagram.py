#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# 2022-03-24
# 一般化相図上に、メタンを位置付ける。
# 温度や圧力を変えた時に、点がどう動くかを見せる。

import vdwp.vdWP as vdWP
import vdwp.crystals as crystals
import CageIntegral.physconst as pc
import CageIntegral.molecule as molecule
import CageIntegral.chempot as chempot

import numpy as np
# import vdwp.interpolate as ip
import matplotlib.pyplot as plt
from itertools import combinations
from logging import getLogger, DEBUG, INFO, basicConfig

# basicConfig(level=DEBUG, format="%(levelname)s %(message)s")
basicConfig(level=INFO, format="%(levelname)s %(message)s")
logger = getLogger()
logger.debug("Debug mode.")




def drawLine(A, B, C):
    """
    Ax + By + C = 0
    """
    if A == 0:
        X = np.linspace(-0.4, 0.6, 100)
        Y = np.zeros_like(X) - C / B
        plt.plot(X, Y)
    elif B == 0:
        Y = np.linspace(-0.4, 0.6, 100)
        X = np.zeros_like(Y) - C / A
        plt.plot(X, Y)
    else:
        X = np.linspace(-0.4, 0.6, 100)
        Y = (-C - A * X) / B
        plt.plot(X, Y)

# def elimempty(values):
#   newvalues = []
#   for value in values:
#     if value != "":
#       newvalues.append(value)
#   return newvalues


#! Test case for methane hydrate


def LoadMoleculeDict(filename):
    file = open(filename, encoding="utf-8")
    moldict = molecule.loadInfo(file)
    # Supplementary info
    for key, mol in moldict.items():
        if sum(mol.moi) == 0:
            mol.symm = 1
            mol.dimen = 0
            mol.dof = 3
        elif np.product(mol.moi) == 0:
            mol.symm = 2  # not always correct
            mol.dimen = 1
            mol.dof = 5
        else:
            mol.symm = 2  # not always correct
            mol.dimen = 3
            mol.dof = 6
        mol.name = key
    # Aliases
    moldict["NEGTHF__"].name = "THF(invert)"
    moldict["CJTHF___"].name = "THF"
    moldict["CPENTANE"].name = "cPentane"
    moldict["CPENTAN+"].name = "cPentane (105%)"
    moldict["CPENTAN-"].name = "cPentane (95%)"
    moldict["CPEN+5L_"].name = "cPentane (105% L)"
    moldict["CPEN-5L_"].name = "cPentane (95% L)"
    moldict["CPENTA++"].name = "cPentane (110%)"
    moldict["CPENTA--"].name = "cPentane (90%)"
    moldict["LJME____"].name = "Methane"
    moldict["TYTHF___"].name = "THF(Yaga)"
    moldict["LJAR____"].name = "Ar"
    moldict["LJET____"].name = "Et"
    moldict["LJXE____"].name = "Xe"
    moldict["LJBR2___"].name = "Br2"
    return moldict


moldict = LoadMoleculeDict("data/DEFR")
#
# 変数の命名則
# f, mu: free energy/chemical potential
#h    : enthalpy
#
# _w   : of water
# _i   : of ice
# _e   : of the empty lattice
# _f   : of the filled lattice
# _g   : of the guest
# _c   : in the cage
#

# User variables
pressure = 101325.00 * 50  # Pa
temperatures = np.array(np.arange(270,301,10))
guest = "LJME____"
host = "TIP4PICE"
structures = ["CS2", "CS1", "HS1", "TS1"]  # , "TS1", "HS1"]

# #Winget's formula
# def VaporPressure(temp, FE, dens, mass):
#     #std pressure
#     p0 = temp * 22.4 / 273.15
#     pvap = p0*dens/mass*np.exp(FE/(0.008314*temp))
#     return pvap


####### chemical potential of gas ######################################
mol = moldict[guest]
logger.info(f"Calculating chemical potential of guest {mol.name}...")
logger.info(f"Regard the guest in gas state at {pressure} Pa.")

stericterm = chempot.StericFix(temperatures, mol.mass, mol.symm, mol.moi)


mu_g = (
    chempot.chempot(temperatures, pressure) +
    chempot.IntegrationFixMinus(temperatures, mol.dimen) + stericterm)

####### structure-dependent terms ######################################

# for hydrate structure types
mu_e = dict()
for structure in structures:
    logger.info(
        f"Calculating chemical potential of empty clathrate {structure}...")
    mu_e[structure] = crystals.U_e[structure] + \
        vdWP.FreeEnergyOfVibration(crystals.nma_file[structure], temperatures)


plt.xlim(-0.4, 0.6)
plt.ylim(-0.4, 0.6)
plt.xlabel(r"\Delta\mu_c^{CS1} - \Delta\mu_c^{HS1}")
plt.ylabel(r"\Delta\mu_c^{CS2} - \Delta\mu_c^{HS1}")

for t, temperature in enumerate(temperatures):
    for s1, s2 in combinations(structures, 2):
        # 共存線の方程式
        # Ax + By + C = 0
        # A = xA - xB
        # B = yA - yB
        # C = mu_eA- mu_eB
        A = crystals.ratios[s1][0] - crystals.ratios[s2][0]
        B = crystals.ratios[s1][1] - crystals.ratios[s2][1]
        C = mu_e[s1][t] - mu_e[s2][t]
        drawLine(A, B, C)

####### cage-dependent terms ###########################################

f_c = vdWP.EncagingFE(temperatures, guest, stericterm)
Deltamu = vdWP.ChemPotByOccupation(temperatures, f_c, mu_g, structures)

# position of methane


# 1. When the temperature changes.
X = Deltamu["CS1"] - Deltamu["HS1"]
Y = Deltamu["CS2"] - Deltamu["HS1"]
plt.plot(X, Y, "-")
# ほとんど動かない。

# 2. When the pressure changes.
temperatures = 273.15
stericterm = chempot.StericFix(temperatures, mol.mass, mol.symm, mol.moi)
f_c = vdWP.EncagingFE(temperatures, guest, stericterm)

for p0 in (50,500,5000,50000,500000):
    pressure = p0*101326
    mu_g = (
        chempot.chempot(temperatures, pressure) +
        chempot.IntegrationFixMinus(temperatures, mol.dimen) + stericterm)
    Deltamu = vdWP.ChemPotByOccupation(temperatures, f_c, mu_g, structures)

    X = Deltamu["CS1"] - Deltamu["HS1"]
    Y = Deltamu["CS2"] - Deltamu["HS1"]
    plt.plot(X, Y, ".")
    # こちらもほとんど動かない。

# 3. Meのsigmaを調節する。
temperatures = 273.15
pressure = 50*101326
stericterm = chempot.StericFix(temperatures, mol.mass, mol.symm, mol.moi)
mu_g = (
    chempot.chempot(temperatures, pressure) +
    chempot.IntegrationFixMinus(temperatures, mol.dimen) + stericterm)

X = []
Y = []
for guest in ("LJMEs2__", "LJMEs3__", "LJMEs4__", "LJMEs5__"):
    f_c = vdWP.EncagingFE(temperatures, guest, stericterm)
    Deltamu = vdWP.ChemPotByOccupation(temperatures, f_c, mu_g, structures)

    x = Deltamu["CS1"] - Deltamu["HS1"]
    y = Deltamu["CS2"] - Deltamu["HS1"]
    X.append(x)
    Y.append(y)

plt.plot(X, Y, "s-")


for guest in ("LJAR____", "LJME____", "LJXE____", "LJET____", "LJBR2___"):
    mol = moldict[guest]
    # 分子が変わると、質量が変わる。ただし、この項はΔμの差には効かないので、計算しなくてもいい。
    stericterm = chempot.StericFix(temperatures, mol.mass, mol.symm, mol.moi)
    # ゲスト分子の理想気体状態での化学ポテンシャル
    mu_g = (
        chempot.chempot(temperatures, pressure) +
        chempot.IntegrationFixMinus(temperatures, mol.dimen) + stericterm)
    # ゲスト分子がケージに閉じこめられた場合の自由エネルギー(化学ポテンシャル)
    f_c = vdWP.EncagingFE(temperatures, guest, stericterm)
    # 混合物の場合は?
    Deltamu = vdWP.ChemPotByOccupation(temperatures, f_c, mu_g, structures)

    x = Deltamu["CS1"] - Deltamu["HS1"]
    y = Deltamu["CS2"] - Deltamu["HS1"]
    plt.plot(x, y, "o")

    plt.annotate(mol.name, # this is the text
                 (x,y), # these are the coordinates to position the label
                 textcoords="offset points", # how to position the text
                 xytext=(0,10), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center



plt.show()
