#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# 2022-03-24
# 一般化相図上に、メタンを位置付ける。
# 温度や圧力を変えた時に、点がどう動くかを見せる。

from vdwp import vdWP, crystals, normalmode, general, cage
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


figure = plt.figure()


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
    moldict["LJKR____"].name = "Kr"
    moldict["LJCO2___"].name = "CO2"
    moldict["LJBR2___"].name = "Br2"
    moldict["EKABR2__"].name = "EkaBr2"
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
temperatures = np.array(np.arange(270, 301, 10))
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
        normalmode.FreeEnergyOfVibration(crystals.nma_file[structure], temperatures)


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
        general.drawLine(A, B, C)

####### cage-dependent terms ###########################################

f_c = cage.EncagingFE(temperatures, guest, stericterm)
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
f_c = cage.EncagingFE(temperatures, guest, stericterm)

for p0 in (50, 500, 5000, 50000, 500000):
    pressure = p0 * 101326
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
pressure = 50 * 101326
stericterm = chempot.StericFix(temperatures, mol.mass, mol.symm, mol.moi)
mu_g = (
    chempot.chempot(temperatures, pressure) +
    chempot.IntegrationFixMinus(temperatures, mol.dimen) + stericterm)

X = []
Y = []
for guest in ("LJMEs2__", "LJMEs3__", "LJMEs4__", "LJMEs5__"):
    f_c = cage.EncagingFE(temperatures, guest, stericterm)
    Deltamu = vdWP.ChemPotByOccupation(temperatures, f_c, mu_g, structures)

    x = Deltamu["CS1"] - Deltamu["HS1"]
    y = Deltamu["CS2"] - Deltamu["HS1"]
    X.append(x)
    Y.append(y)

plt.plot(X, Y, "s-")

# 4. いろんな分子。BR2のパラメータは手で調節した。

guests = """
LJME____ LJMEs3__ LJMEs4__ LJMEs5__ LJMEe0__ LJMEe2__ \
	LJAR____ LJXE____ LJET____ LJBR2___ LJCS2___ LJKR____ LJCO2___ EKABR2__ \
	LJCF4___ HTCO2___ \
	CH2CL2__ CHCL3___ CH3CL___ CH3BR___ C2H4F2__ 	C2H5F___ 	C2HF5___ \
		CCL2F2__	CCL3F___	CH2CLF__	CH2FCF3_	CH3CL___ \
		CHCL2F__	CHCLF2__	CL2FCCH3	CLF2CCH3	PROPANE_ \
		ISOBUTAN
        """.split()

# ("LJAR____", "LJME____", "LJXE____", "LJET____", "LJBR2___", "LJPR____", "LJKR____", "LJCO2___", "EKABR2__"):
for guest in guests:
    mol = moldict[guest]
    # 分子が変わると、質量が変わる。ただし、この項はΔμの差には効かないので、計算しなくてもいい。
    stericterm = chempot.StericFix(temperatures, mol.mass, mol.symm, mol.moi)
    # ゲスト分子の理想気体状態での化学ポテンシャル
    mu_g = (
        chempot.chempot(temperatures, pressure) +
        chempot.IntegrationFixMinus(temperatures, mol.dimen) + stericterm)
    # ゲスト分子がケージに閉じこめられた場合の自由エネルギー(化学ポテンシャル)
    f_c = cage.EncagingFE(temperatures, guest, stericterm)
    # 混合物の場合は?
    Deltamu = vdWP.ChemPotByOccupation(temperatures, f_c, mu_g, structures)

    x = Deltamu["CS1"] - Deltamu["HS1"]
    y = Deltamu["CS2"] - Deltamu["HS1"]
    plt.plot(x, y, "o")

    plt.annotate(mol.name,  # this is the text
                 (x, y),  # these are the coordinates to position the label
                 textcoords="offset points",  # how to position the text
                 xytext=(0, 10),  # distance from text to points (x,y)
                 ha='center')  # horizontal alignment can be left, right or center


# 5. Mixture of methane and ethane.
temperatures = 273.15
p0 = 101326 * 50  # 50 bar


pairs = (("LJME____", "LJET____"),
         ("LJXE____", "LJET____"),
         ("LJME____", "C2H4F2__"),
         ("LJXE____", "C2H4F2__"),
         ("LJCO2___", "C2H4F2__"),
         ("LJME____", "LJBR2___"),
         ("LJET____", "LJBR2___"),
         ("LJME____", "LJXE____"),
         ("LJME____", "LJCO2___"),
         ("LJCO2___", "LJBR2___"),
         ("LJCO2___", "C2H5F___"),
         ("LJME____", "C2H5F___"),
         ("LJCO2___", "EKABR2__"),
         ("LJET____", "LJCO2___"),
         ("LJXE____", "LJBR2___"))


def DoubleClathrate(me, et, ticks=np.linspace(0.0, 1.0, 100)):
    mol_me = moldict[me]
    stericterm_me = chempot.StericFix(
        temperatures, mol_me.mass, mol_me.symm, mol_me.moi)
    f_me = cage.EncagingFE(temperatures, me, stericterm_me)

    mol_et = moldict[et]
    stericterm_et = chempot.StericFix(
        temperatures, mol_et.mass, mol_et.symm, mol_et.moi)
    f_et = cage.EncagingFE(temperatures, et, stericterm_et)

    X = []
    Y = []
    for r in ticks:
        p_me = (1 - r) * p0
        p_et = r * p0

        if p_me == 0:
            mu_et = (
                chempot.chempot(
                    temperatures,
                    p_et) +
                chempot.IntegrationFixMinus(
                    temperatures,
                    mol_et.dimen) +
                stericterm_et)
            Deltamu = vdWP.ChemPotByOccupation(
                temperatures, f_et, mu_et, structures)
        elif p_et == 0:
            mu_me = (
                chempot.chempot(
                    temperatures,
                    p_me) +
                chempot.IntegrationFixMinus(
                    temperatures,
                    mol_me.dimen) +
                stericterm_me)
            Deltamu = vdWP.ChemPotByOccupation(
                temperatures, f_me, mu_me, structures)
        else:
            mu_me = (
                chempot.chempot(
                    temperatures,
                    p_me) +
                chempot.IntegrationFixMinus(
                    temperatures,
                    mol_me.dimen) +
                stericterm_me)
            mu_et = (
                chempot.chempot(
                    temperatures,
                    p_et) +
                chempot.IntegrationFixMinus(
                    temperatures,
                    mol_et.dimen) +
                stericterm_et)
            Deltamu = vdWP.ChemPotByOccupation(
                temperatures, (f_me, f_et), (mu_me, mu_et), structures)

        x = Deltamu["CS1"] - Deltamu["HS1"]
        y = Deltamu["CS2"] - Deltamu["HS1"]
        X.append(x)
        Y.append(y)
    return X, Y


for me, et in pairs:
    ticks = np.concatenate(
        [np.arange(0, 0.01, 0.0001), np.arange(0.01, 1, 0.01)])
    X, Y = DoubleClathrate(me, et, ticks)
    plt.plot(X, Y, "-")
    X, Y = DoubleClathrate(me, et, np.linspace(0.0, 1.0, 11))
    plt.plot(X, Y, ".")


plt.show()
figure.savefig("phasediagram.pdf")
