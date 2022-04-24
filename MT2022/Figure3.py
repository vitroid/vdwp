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

# figure = plt.figure(figsize=(5, 5))
fig, axes = plt.subplots(
    nrows=1, ncols=2, figsize=(
        9, 5), sharey=True, gridspec_kw={
            'wspace': 0})

ax = axes[0]
ax.set_ylabel(
    r"$(\Delta\mu_c^\mathrm{CS2} - \Delta\mu_c^\mathrm{HS1}) / \mathrm{kJ~mol}^{-1}$")

for ax in axes:
    ax.set_xlim(-0.4, 0.2)
    ax.set_ylim(0.05, 0.65)
    ax.set_xlabel(
        r"$(\Delta\mu_c^\mathrm{CS1} - \Delta\mu_c^\mathrm{HS1}) / \mathrm{kJ~mol}^{-1}$")
    # ax.axis("square")

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
        general.drawLine(A, B, C, style="-k", ax=ax)

    # manually labelled
    ax.annotate("I",  # this is the text
                xy=(0.4, 0.9),  # these are the coordinates to position the label
                xycoords="axes fraction",
                fontsize=24, )
    ax.annotate("II",  # this is the text
                xy=(0.8, 0.05),  # these are the coordinates to position the label
                xycoords="axes fraction",
                fontsize=24, )
    ax.annotate("III",  # this is the text
                xy=(0.8, 0.9),  # these are the coordinates to position the label
                xycoords="axes fraction",
                fontsize=24, )

# plt.show()

for a, ax in enumerate(axes):
    ####### cage-dependent terms ###########################################
    for name, gas in inter.items():
        sigma = gas.sig
        epsilon = gas.epsK * 8.314 / 1000  # in kJ/mol
        f_c = dict()
        for cage, R in radii.items():
            f_c[cage] = fvalue({R: nmemb[cage]}, sigma, epsilon, beta) + stericterm
        pressures = (50,30,10)
        if name == "cC3H6":
            pressures = (50, 30, 10, 0.75)
        # elif name == "Br2":
        #     pressures = (50, 30, 10, 1.0)
        x = []
        y = []
        for pressure in pressures:
            if a == 1 and pressure < 50:
                continue
            if a == 1 and name not in ("Br2", "Methane", "Ethane", "C2H4", "cC3H6", "Xe"):
                continue

            mu_g0 = (
                chempot.chempot(temperatures, pressure*101326) +
                chempot.IntegrationFixMinus(temperatures, dimen=0) + stericterm)
            Deltamu = vdWP.ChemPotByOccupation(temperatures, f_c, mu_g0, crystals.names)

            X = Deltamu["CS1"] - Deltamu["HS1"]
            Y = Deltamu["CS2"] - Deltamu["HS1"]
            label = gas.TeX
            color = "#000"
            symbol = "o"
            if pressure <= 1.0:
                label += f" ({pressure} bar)"
                color = "#888"
            elif pressure < 50.0:
                label = ""
                if pressure == 30.0:
                    color = "green"
                    symbol = "."
                elif pressure == 10.0:
                    color = "blue"
                    symbol = "."

            # ax.plot(X, Y, symbol, color=color)
            x.append(X)
            y.append(Y)

            if name in ("Methane", "Kr", "n-Butane", "C2H4", "Xe") or (name == "cC3H6" and pressure == 50):
                ha = "right"
                xytext = (-2, -16)
            else:
                ha = "left"
                xytext = (5, 5)


            ax.annotate(label,  # this is the text
                        (X, Y),  # these are the coordinates to position the label
                        textcoords="offset points",  # how to position the text
                        xytext=xytext,  # distance from text to points (x,y)
                        color=color,
                        ha=ha)  # horizontal alignment can be left, right or center
        # if a == 0:
        if len(x) > 3:
            ax.plot(x[2:4], y[2:4], ":k", linewidth=0.5)
        ax.plot(x[:3], y[:3], "-k", linewidth=0.5)
        s = (40,30,30,20)
        markers = ("o", "^", "+", "o")
        facecolors=("white", "blue", "green", "white")
        edgecolors=("black", "blue", "green", "gray")
        for i, (X, Y) in enumerate(zip(x,y)):
            ax.scatter(X, Y, s=s[i], marker= markers[i], facecolors=facecolors[i], edgecolors=edgecolors[i])
        # ax.plot(x[0], y[0], "o", color="black", fillcolor="white")

ticks = np.concatenate(
    [np.array([0.0]), np.logspace(-5, 0.0, 300)])  # 1e-4 .. 1e0
ax = axes[1]

for a, b, color in [("Methane", "Ethane", "red"), ("Methane", "C2H4",
                                                   "blue"), ("Xe", "Br2", "green"), ("Methane", "cC3H6", "#cc0")]:
    pressures = np.zeros(2)
    X = []
    Y = []
    for frac in ticks:
        pressures[1] = 50.0*frac
        pressures[0] = 50.0*(1.0 - frac)
        x, y = MultipleClathrate(
            (inter[a], inter[b]), pressures * 101326, temperatures, crystals.names)
        X.append(x)
        Y.append(y)
    ax.plot(X, Y, "-", color=color)
    X = []
    Y = []
    for frac in np.linspace(0, 1.0, 11):
        pressures[1] = 50.0*frac
        pressures[0] = 50.0*(1.0 - frac)
        x, y = MultipleClathrate(
            (inter[a], inter[b]), pressures * 101326, temperatures, crystals.names)
        X.append(x)
        Y.append(y)
    ax.plot(X, Y, ".", color=color)


for a, b in [("Ethane", "Br2"), ("Methane", "Br2"),
             ("Xe", "Br2"), ("C2H4", "Br2"), ("cC3H6", "Br2")]:
    pressures = np.zeros(2)
    X = []
    Y = []
    for frac in ticks:
        pressures[1] = 50.0*frac
        pressures[0] = 50.0*(1.0 - frac)
        x, y = MultipleClathrate(
            (inter[a], inter[b]), pressures * 101326, temperatures, crystals.names)
        X.append(x)
        Y.append(y)
    ax.plot(X, Y, "-k", linewidth=0.5)


plt.tight_layout()
plt.show()
fig.savefig("Figure3.pdf")
fig.savefig("Figure3.png")
