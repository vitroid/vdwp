#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# phasediagram.pyからの変更点:
# ケージ積分の代わりにljdを用いる。
# 必然的に、単原子モデル以外は扱えない。
# 相境界の条件は50気圧、273 Kのみにする。
# moldictの利用をやめる。

import vdwp.vdWP as vdWP
import vdwp.crystals as crystals
from CageIntegral.physconst import NkB, NA
import CageIntegral.molecule as molecule
import CageIntegral.chempot as chempot
from ljd.ljd import fvalue

import numpy as np
# import vdwp.interpolate as ip
import matplotlib.pyplot as plt
from itertools import combinations
from logging import getLogger, DEBUG, INFO, basicConfig
from attrdict import AttrDict

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
                ff[cage] = fvalue({R: nmemb[cage]}, gas.sig, gas.epsK * 8.314/1000, beta)
            f.append(ff)
    Deltamu = vdWP.ChemPotByOccupation(temperatures, f, mu, structures)

    X = Deltamu["CS1"] - Deltamu["HS1"]
    Y = Deltamu["CS2"] - Deltamu["HS1"]
    return X, Y

gastable = """
Methane 3.758 148.6 $\mathrm{CH}_4$
Ethane 4.520 208.8
Ne 2.749 35.6
Ar 3.405 119.8
Kr 3.60 171.0
Xe 4.047 231.0
Br2 4.933 488 $\mathrm{Br}_2$
CO2 4.486 189.0 $\mathrm{CO}_2$
CS2 4.438 488 $\mathrm{CS_2}$
N2O 4.59 189 $\mathrm{N_2O}$
CF4 4.70 152.5 $\mathrm{CF_4}$
n-Butane 4.997 410 {\it n}-Butane
""".splitlines()

tip4pice = AttrDict({"sig": 3.1668, "epsK": 106.1})

gases = dict()
inter = dict()
for gas in gastable:
    cols = gas.rstrip().split(maxsplit=3)
    if len(cols) == 0:
        continue
    if len(cols) == 3:
        tex = r"$\mathrm{" + cols[0] + r"}$"
        cols.append(tex)
    # print(len(cols))
    name, sig, epsK, tex = cols
    gases[name] = AttrDict({"sig": float(sig), "epsK": float(epsK), "TeX": tex})
    inter[name] = AttrDict({"sig": (float(sig)+tip4pice.sig)/2,
                            "epsK": (float(epsK)*tip4pice.epsK)**0.5,
                            "TeX": tex})

# User variables
pressure = 101325.00 * 50  # Pa
temperatures = 273.15
# guest = "Methane"
structures = ["CS2", "CS1", "HS1", "TS1"]  # , "TS1", "HS1"]

####### chemical potential of gas ######################################
stericterm = chempot.StericFix(temperatures, mass=1.0, symm=1, moi=(0,0,0))

mu_g = (
    chempot.chempot(temperatures, pressure) +
    chempot.IntegrationFixMinus(temperatures, dimen=0) + stericterm)

####### structure-dependent terms ######################################

# for hydrate structure types
mu_e = dict()
for structure in structures:
    logger.info(
        f"Calculating chemical potential of empty clathrate {structure}...")
    mu_e[structure] = crystals.U_e[structure] + \
        vdWP.FreeEnergyOfVibration(crystals.nma_file[structure], temperatures)

figure = plt.figure(figsize=(5,5))

plt.rcParams['text.usetex'] = True

plt.xlim(-0.3, 0.3)
plt.ylim(-0.1, 0.5)
plt.xlabel(r"$(\Delta\mu_c^\mathrm{CS1} - \Delta\mu_c^\mathrm{HS1}) / \mathrm{kJ~mol}^{-1}$")
plt.ylabel(r"$(\Delta\mu_c^\mathrm{CS2} - \Delta\mu_c^\mathrm{HS1}) / \mathrm{kJ~mol}^{-1}$")
plt.axis("square")

for s1, s2 in combinations(structures, 2):
    # 共存線の方程式
    # Ax + By + C = 0
    # A = xA - xB
    # B = yA - yB
    # C = mu_eA- mu_eB
    A = crystals.ratios[s1][0] - crystals.ratios[s2][0]
    B = crystals.ratios[s1][1] - crystals.ratios[s2][1]
    C = mu_e[s1] - mu_e[s2]
    drawLine(A, B, C)

def determine_phase(mu_e, Deltamu, structures):
    mumin = 0
    stmin = None
    for s in structures:
        if mu_e[s] + Deltamu[s] < mumin:
            mumin = mu_e[s] + Deltamu[s]
            stmin = s
    return stmin

####### cage-dependent terms ###########################################
radii = {12: 3.894043995956962, 14: 4.3269124645840025, 15: 4.479370276434863, 16: 4.694779956202102}
nmemb = {12:20, 14:24, 15:26, 16:28}
beta = 1.0 / (NkB * temperatures)
for name, gas in inter.items():
    sigma = gas.sig
    epsilon = gas.epsK * 8.314 / 1000# in kJ/mol
    f_c = dict()
    for cage, R in radii.items():
        f_c[cage] = fvalue({R: nmemb[cage]}, sigma, epsilon, beta) + stericterm
    Deltamu = vdWP.ChemPotByOccupation(temperatures, f_c, mu_g, structures)

    X = Deltamu["CS1"] - Deltamu["HS1"]
    Y = Deltamu["CS2"] - Deltamu["HS1"]
    plt.plot(X, Y, "o")
    plt.annotate(gas.TeX, # this is the text
                (X, Y), # these are the coordinates to position the label
                textcoords="offset points", # how to position the text
                xytext=(0,10), # distance from text to points (x,y)
                ha='center') # horizontal alignment can be left, right or center

    phase = determine_phase(mu_e, Deltamu, structures)

for a, b in ("Methane", "Ethane"), ("Methane", "Xe"), ("Xe", "Br2"):
    pressures = np.zeros(2)
    X = []
    Y = []
    for pressures[1] in np.concatenate([np.linspace(0.000, 0.99, 1000), np.linspace(1.0, 50 - pressures[0],100)]):
        pressures[0] = 50.0 - pressures[1]
        x,y = MultipleClathrate((inter[a], inter[b]), pressures*101326, temperatures, structures)
        X.append(x)
        Y.append(y)
    plt.plot(X,Y, "-")
    X = []
    Y = []
    for pressures[0] in np.linspace(0.0, 50, 11):
        pressures[1] = 50.0 - pressures[0]
        x,y = MultipleClathrate((inter[a], inter[b]), pressures*101326, temperatures, structures)
        X.append(x)
        Y.append(y)
    plt.plot(X,Y, ".")


# pressures = np.zeros(3)
# for pressures[0] in np.linspace(0.0, 45.0, 10):
#     X = []
#     Y = []
#     for pressures[2] in np.concatenate([np.linspace(0.000, 0.99, 1000), np.linspace(1.0, 50 - pressures[0],100)]):
#         pressures[1] = 50.0 - pressures[0] - pressures[2]
#         x,y = MultipleClathrate(gases, pressures*101326, temperatures, structures)
#         X.append(x)
#         Y.append(y)
#     plt.plot(X,Y, "-")

plt.show()

figure.savefig("simple.pdf")










def DoubleClathrate(g1, g2, beta, pressure, ticks=np.linspace(0.0, 1.0, 100)):
    f1 = dict()
    f2 = dict()
    for cage, R in radii.items():
        # sigma and epsilon must be the intermolecular ones.
        f1[cage] = fvalue({R: nmemb[cage]}, g1.sig, g1.epsK*8.314/1000, beta)
        f2[cage] = fvalue({R: nmemb[cage]}, g2.sig, g2.epsK*8.314/1000, beta)

    phases = set()
    for r in ticks:
        p1 = (1-r)*pressure
        p2 = r*pressure

        if p1 == 0:
            mu2 = (
                chempot.chempot(temperatures, p2) +
                chempot.IntegrationFixMinus(temperatures, dimen=0))
            Deltamu = vdWP.ChemPotByOccupation(temperatures, f2, mu2, structures)
        elif p2 == 0:
            mu1 = (
                chempot.chempot(temperatures, p1) +
                chempot.IntegrationFixMinus(temperatures, dimen=0))
            Deltamu = vdWP.ChemPotByOccupation(temperatures, f1, mu1, structures)
        else:
            mu1 = (
                chempot.chempot(temperatures, p1) +
                chempot.IntegrationFixMinus(temperatures, dimen=0))
            mu2 = (
                chempot.chempot(temperatures, p2) +
                chempot.IntegrationFixMinus(temperatures, dimen=0))
            Deltamu = vdWP.ChemPotByOccupation(temperatures, (f1, f2), (mu1, mu2), structures)

        phase = determine_phase(mu_e, Deltamu, structures)
        phases.add(phase)
    return phases, phase


def mark(sig, epsK, phases, lastphase, ax):
    markers = {1:"o", 2:"+", 3:"^"}
    if lastphase == "CS1":
        color = "lightgreen"
    elif lastphase == "CS2":
        color = "#4666ff"
    else:
        color = "brown"
    # I-III-Iの場合を区別して描く。
    if lastphase == "CS1" and "TS1" in phases and len(phases) == 2:
        ax.plot(sig, epsK, color="orange", marker=".")
    else:
        ax.plot(sig, epsK, color=color, marker=markers[len(phases)])
        if len(phases) == 3:
            print(sig, epsK, phases, lastphase)

def cities(ax, gases):
    # 地図の上に、都市(分子)を描く。
    for name, gas in gases.items():
        sig = gas.sig
        eps = gas.epsK
        ax.plot(sig, eps, "ok")
        ax.annotate(gas.TeX, # this is the text
                    (sig, eps), # these are the coordinates to position the label
                    textcoords="offset points", # how to position the text
                    xytext=(0,10), # distance from text to points (x,y)
                    ha='center') # horizontal alignment can be left, right or center


fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8, 4), sharey=True, gridspec_kw={'wspace': 0})
# Meと混合すると、I-II-III-I転移する条件をさがす。(てあたりしだい)
guest = "Methane"
axes[0].set_xlabel(r"$\sigma_g / \mathrm\AA$")
axes[0].set_ylabel(r"$\epsilon_g / \mathrm K$")
axes[0].set_xlim(3.5, 5.0)
axes[0].set_ylim(120.0, 600)
xyticks = 61
# ticks = np.concatenate([np.arange(0, 0.01, 0.0001), np.arange(0.01, 1, 0.01)])
ticks = np.arange(0., 1, 0.001)
for sig in np.linspace(3.5, 5.0, xyticks):
    for epsK in np.linspace(120.0, 600, xyticks):
        inter2 = AttrDict({"sig": (tip4pice.sig+sig)/2,
                           "epsK": (tip4pice.epsK*epsK)**0.5})
        phases, lastphase = DoubleClathrate(inter[guest], inter2, beta, pressure, ticks = ticks)
        mark(sig, epsK, phases, lastphase, ax=axes[0])

cities(ax=axes[0], gases=gases)

# Xeと混合すると、I-II-III-I転移する条件をさがす。(てあたりしだい)
markers = {1:"o", 2:"+", 3:"^"}
guest = "Xe"
axes[1].set_xlabel(r"$\sigma_g / \mathrm\AA$")
# axes[1].set_ylabel(r"\epsilon_g / K")
axes[1].set_xlim(4.6, 5.0)
axes[1].set_ylim(120.0, 600)
xyticks = 61 # 60
# ticks = np.concatenate([np.arange(0, 0.01, 0.0001), np.arange(0.01, 1, 0.01)])
ticks = np.arange(0., 1, 0.001)
for sig in np.linspace(4.6, 5.0, xyticks):
    for epsK in np.linspace(120.0, 600., xyticks):
        inter2 = AttrDict({"sig": (tip4pice.sig+sig)/2,
                           "epsK": (tip4pice.epsK*epsK)**0.5})
        phases, lastphase = DoubleClathrate(inter[guest], inter2, beta, pressure, ticks = ticks)
        mark(sig, epsK, phases, lastphase, ax=axes[1])

cities(ax=axes[1], gases=gases)
plt.tight_layout()

plt.show()
fig.savefig("s-e.pdf")
