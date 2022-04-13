#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# phasediagram.pyからの変更点:
# ケージ積分の代わりにljdを用いる。
# 必然的に、単原子モデル以外は扱えない。
# 相境界の条件は50気圧、273 Kのみにする。
# moldictの利用をやめる。

from vdwp import vdWP, crystals, normalmode, general, chempot
from vdwp.physconst import NkB, NA
# import CageIntegral.molecule as molecule
from ljd.ljd import fvalue
from LJparam import gases, inter

import numpy as np
# import vdwp.interpolate as ip
import matplotlib.pyplot as plt
from itertools import combinations
from logging import getLogger, DEBUG, INFO, basicConfig

# basicConfig(level=DEBUG, format="%(levelname)s %(message)s")
basicConfig(level=INFO, format="%(levelname)s %(message)s")
logger = getLogger()
logger.debug("Debug mode.")


def determine_phase(mu_e, Deltamu, structures):
    mumin = 0
    stmin = None
    for s in structures:
        if mu_e[s] + Deltamu[s] < mumin:
            mumin = mu_e[s] + Deltamu[s]
            stmin = s
    return stmin


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


def MultipleClathrate(gases, pressures, temperatures, structures):
    beta = 1 / (NkB * temperatures)
    f = []
    mu = []
    for gas, pressure in zip(gases, pressures):
        if pressure > 0.0:
            mu.append(
                chempot.chempot(temperatures, pressure) +
                chempot.IntegrationFixMinus(temperatures, dimen=0))
            ff = dict()
            for cage, R in radii.items():
                # sigma and epsilon must be the intermolecular ones.
                ff[cage] = fvalue({R: nmemb[cage]}, gas.sig,
                                  gas.epsK * 8.314 / 1000, beta)
            f.append(ff)
    Deltamu = vdWP.ChemPotByOccupation(temperatures, f, mu, structures)

    X = Deltamu["CS1"] - Deltamu["HS1"]
    Y = Deltamu["CS2"] - Deltamu["HS1"]
    return X, Y


radii = {
    12: 3.894043995956962,
    14: 4.3269124645840025,
    15: 4.479370276434863,
    16: 4.694779956202102}
nmemb = {12: 20, 14: 24, 15: 26, 16: 28}

# User variables
pressure = 101325.00 * 50  # Pa
temperatures = 273.15
beta = 1.0 / (NkB * temperatures)


####### chemical potential of gas ######################################
# Here we do not care it because this term in mu_g and f_g cancels.
# chempot.StericFix(temperatures, mass=1.0, symm=1, moi=(0,0,0))
stericterm = 0.0

mu_g = (
    chempot.chempot(temperatures, pressure) +
    chempot.IntegrationFixMinus(temperatures, dimen=0) + stericterm)

# for hydrate structure types
mu_e = dict()
for structure in crystals.names:
    logger.info(
        f"Calculating chemical potential of empty clathrate {structure}...")
    mu_e[structure] = crystals.U_e[structure] + \
        normalmode.FreeEnergyOfVibration(crystals.nma_file[structure], temperatures)

#plt.rcParams['text.usetex'] = True
plt.rcParams["font.size"] = 14
plt.rcParams["font.family"] = "sans-serif"

figure = plt.figure(figsize=(5, 5))


plt.xlim(-0.4, 0.3)
plt.ylim(-0.1, 0.6)
plt.xlabel(
    r"$(\Delta\mu_c^\mathrm{CS1} - \Delta\mu_c^\mathrm{HS1}) / \mathrm{kJ~mol}^{-1}$")
plt.ylabel(
    r"$(\Delta\mu_c^\mathrm{CS2} - \Delta\mu_c^\mathrm{HS1}) / \mathrm{kJ~mol}^{-1}$")
plt.axis("square")

for s1, s2 in combinations(crystals.names, 2):
    # 共存線の方程式
    # Ax + By + C = 0
    # A = xA - xB
    # B = yA - yB
    # C = mu_eA- mu_eB
    if s1 == "HS1" or s2 == "HS1":
        # skip metastable coexistences
        continue
    A = crystals.ratios[s1][0] - crystals.ratios[s2][0]
    B = crystals.ratios[s1][1] - crystals.ratios[s2][1]
    C = mu_e[s1] - mu_e[s2]
    general.drawLine(A, B, C, style="-k")

# manually labelled
plt.annotate("I",  # this is the text
             xy=(0.4, 0.9),  # these are the coordinates to position the label
             xycoords="axes fraction",
             fontsize=24, )
plt.annotate("II",  # this is the text
             xy=(0.8, 0.4),  # these are the coordinates to position the label
             xycoords="axes fraction",
             fontsize=24, )
plt.annotate("III",  # this is the text
             xy=(0.8, 0.9),  # these are the coordinates to position the label
             xycoords="axes fraction",
             fontsize=24, )


####### cage-dependent terms ###########################################
for name, gas in inter.items():
    sigma = gas.sig
    epsilon = gas.epsK * 8.314 / 1000  # in kJ/mol
    f_c = dict()
    for cage, R in radii.items():
        f_c[cage] = fvalue({R: nmemb[cage]}, sigma, epsilon, beta) + stericterm
    Deltamu = vdWP.ChemPotByOccupation(temperatures, f_c, mu_g, crystals.names)

    X = Deltamu["CS1"] - Deltamu["HS1"]
    Y = Deltamu["CS2"] - Deltamu["HS1"]
    plt.plot(X, Y, "ok")
    if name in ("Methane", "Kr", "n-Butane"):
        ha = "right"
        xytext = (-2, -16)
    else:
        ha = "left"
        xytext = (2, 2)
    plt.annotate(gas.TeX,  # this is the text
                 (X, Y),  # these are the coordinates to position the label
                 textcoords="offset points",  # how to position the text
                 xytext=xytext,  # distance from text to points (x,y)
                 ha=ha)  # horizontal alignment can be left, right or center

# for guest, label in [("CPENTANE", "cPen")]:
#     f_c = vdWP.EncagingFE(temperatures, guest, stericterm)
#     Deltamu = vdWP.ChemPotByOccupation(temperatures, f_c, mu_g, crystals.names)

#     X = Deltamu["CS1"] - Deltamu["HS1"]
#     Y = Deltamu["CS2"] - Deltamu["HS1"]
#     plt.plot(X, Y, "ok")
#     ha = "left"
#     xytext=(2,2)
#     plt.annotate(label, # this is the text
#                 (X, Y), # these are the coordinates to position the label
#                 textcoords="offset points", # how to position the text
#                 xytext=xytext, # distance from text to points (x,y)
#                 ha=ha) # horizontal alignment can be left, right or center


for a, b, color in [("Methane", "Ethane", "red"), ("Methane", "C2H4",
                                                   "blue"), ("Xe", "Br2", "green"), ("Methane", "cC3H6", "#cc0")]:
    pressures = np.zeros(2)
    X = []
    Y = []
    for pressures[1] in np.concatenate(
            [np.linspace(0.000, 0.99, 1000), np.linspace(1.0, 50 - pressures[0], 100)]):
        pressures[0] = 50.0 - pressures[1]
        x, y = MultipleClathrate(
            (inter[a], inter[b]), pressures * 101326, temperatures, crystals.names)
        X.append(x)
        Y.append(y)
    plt.plot(X, Y, "-", color=color)
    X = []
    Y = []
    for pressures[0] in np.linspace(0.0, 50, 11):
        pressures[1] = 50.0 - pressures[0]
        x, y = MultipleClathrate(
            (inter[a], inter[b]), pressures * 101326, temperatures, crystals.names)
        X.append(x)
        Y.append(y)
    plt.plot(X, Y, ".", color=color)


for a, b in [("Ethane", "Br2"), ("CO2", "Br2"),
             ("N2O", "Br2"), ("C2H4", "Br2"), ("cC3H6", "Br2")]:
    pressures = np.zeros(2)
    X = []
    Y = []
    for pressures[1] in np.concatenate(
            [np.linspace(0.000, 0.99, 1000), np.linspace(1.0, 50 - pressures[0], 100)]):
        pressures[0] = 50.0 - pressures[1]
        x, y = MultipleClathrate(
            (inter[a], inter[b]), pressures * 101326, temperatures, crystals.names)
        X.append(x)
        Y.append(y)
    if a == "cC3H6":
        plt.plot(X, Y, "--k", linewidth=0.5, dashes=(10,5))
    else:
        plt.plot(X, Y, "-k", linewidth=0.5)


plt.tight_layout()
plt.show()
figure.savefig("Figure3.pdf")
figure.savefig("Figure3.png")
