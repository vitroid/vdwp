#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# phasediagram.pyからの変更点:
# ケージ積分の代わりにljdを用いる。
# 必然的に、単原子モデル以外は扱えない。
# 相境界の条件は50気圧、273 Kのみにする。

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
    moldict["LJKR____"].name = "Kr"
    moldict["LJCO2___"].name = "CO2"
    moldict["LJBR2___"].name = "Br2"
    moldict["LJBR2LJD"].name = "Br2"
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


def MultipleClathrate(gases, pressures, temperatures, structures):
    beta = 1 / (NkB * temperatures)
    f = []
    mu = []
    for gas, pressure in zip(gases, pressures):
        if pressure > 0.0:
            mu.append(
                chempot.chempot(temperatures, pressure) +
                chempot.IntegrationFixMinus(temperatures, gas.dimen))
            ff = dict()
            for cage, R in radii.items():
                # sigma and epsilon must be the intermolecular ones.
                ff[cage] = fvalue({R: nmemb[cage]}, gas.sigma, gas.epsilon, beta)
            f.append(ff)
    Deltamu = vdWP.ChemPotByOccupation(temperatures, f, mu, structures)

    X = Deltamu["CS1"] - Deltamu["HS1"]
    Y = Deltamu["CS2"] - Deltamu["HS1"]
    return X, Y



# User variables
pressure = 101325.00 * 50  # Pa
temperatures = 273.15
guest = "LJME____"
host = "TIP4PICE"
structures = ["CS2", "CS1", "HS1", "TS1"]  # , "TS1", "HS1"]

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
for guest in ("LJME____", "LJET____", "LJXE____", "LJBR2LJD", "LJAR____", "LJKR____", "LJNE____", "LJCO2HT_"):
    sigma = (moldict[host].intr[0][1] + moldict[guest].intr[0][1])/2
    epsilon = (moldict[host].intr[0][0] * moldict[guest].intr[0][0])**0.5 # in kJ/mol
    f_c = dict()
    for cage, R in radii.items():
        f_c[cage] = fvalue({R: nmemb[cage]}, sigma, epsilon, beta) + stericterm
    Deltamu = vdWP.ChemPotByOccupation(temperatures, f_c, mu_g, structures)

    X = Deltamu["CS1"] - Deltamu["HS1"]
    Y = Deltamu["CS2"] - Deltamu["HS1"]
    plt.plot(X, Y, ".")
    plt.annotate(moldict[guest].name, # this is the text
                (X, Y), # these are the coordinates to position the label
                textcoords="offset points", # how to position the text
                xytext=(0,10), # distance from text to points (x,y)
                ha='center') # horizontal alignment can be left, right or center

    phase = determine_phase(mu_e, Deltamu, structures)
    print(f"{guest}: {phase}")


from attrdict import AttrDict
gases = []
for guest in ["LJME____", "LJET____", "LJBR2LJD", "LJXE____", ]:
    gases.append(AttrDict({"sigma": (moldict[host].intr[0][1] + moldict[guest].intr[0][1])/2,
                           "epsilon": (moldict[host].intr[0][0] * moldict[guest].intr[0][0])**0.5,
                           "dimen": 0,
                           "name": moldict[guest].name}))

for a, b in ((gases[0],gases[1]), (gases[0],gases[3]), (gases[2],gases[3]),):
    pressures = np.zeros(2)
    X = []
    Y = []
    for pressures[0] in np.concatenate([np.linspace(0.000, 0.99, 1000), np.linspace(1.0, 50 - pressures[0],100)]):
        pressures[1] = 50.0 - pressures[0]
        x,y = MultipleClathrate((a, b), pressures*101326, temperatures, structures)
        X.append(x)
        Y.append(y)
    plt.plot(X,Y, "-")
    X = []
    Y = []
    for pressures[0] in np.linspace(0.0, 50, 11):
        pressures[1] = 50.0 - pressures[0]
        x,y = MultipleClathrate((a, b), pressures*101326, temperatures, structures)
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
        f1[cage] = fvalue({R: nmemb[cage]}, g1.sigma, g1.epsilon, beta)
        f2[cage] = fvalue({R: nmemb[cage]}, g2.sigma, g2.epsilon, beta)

    phases = set()
    for r in ticks:
        p1 = (1-r)*pressure
        p2 = r*pressure

        if p1 == 0:
            mu2 = (
                chempot.chempot(temperatures, p2) +
                chempot.IntegrationFixMinus(temperatures, g2.dimen))
            Deltamu = vdWP.ChemPotByOccupation(temperatures, f2, mu2, structures)
        elif p2 == 0:
            mu1 = (
                chempot.chempot(temperatures, p1) +
                chempot.IntegrationFixMinus(temperatures, g1.dimen))
            Deltamu = vdWP.ChemPotByOccupation(temperatures, f1, mu1, structures)
        else:
            mu1 = (
                chempot.chempot(temperatures, p1) +
                chempot.IntegrationFixMinus(temperatures, g1.dimen))
            mu2 = (
                chempot.chempot(temperatures, p2) +
                chempot.IntegrationFixMinus(temperatures, g2.dimen))
            Deltamu = vdWP.ChemPotByOccupation(temperatures, (f1, f2), (mu1, mu2), structures)

        phase = determine_phase(mu_e, Deltamu, structures)
        phases.add(phase)
    return phases, phase


from attrdict import AttrDict


def cities(ax):
    # 地図の上に、都市(分子)を描く。
    # import fridge as fr
    # for row in fr.refrigerants().itertuples():
    #     index, name, sig, eps = row
    #     print(name, sig, eps)
    #     plt.plot(sig, eps, "ok")
    #     plt.annotate(name, # this is the text
    #                 (sig, eps), # these are the coordinates to position the label
    #                 textcoords="offset points", # how to position the text
    #                 xytext=(0,10), # distance from text to points (x,y)
    #                 ha='center') # horizontal alignment can be left, right or center

    for guest in ("LJME____", "LJET____", "LJXE____", "LJBR2LJD", "LJAR____", "LJNE____", "LJKR____", "LJCO2HT_"):
        sig = moldict[guest].intr[0][1]
        eps = moldict[guest].intr[0][0] / 8.314 * 1000 # K
        ax.plot(sig, eps, "ok")
        ax.annotate(moldict[guest].name, # this is the text
                    (sig, eps), # these are the coordinates to position the label
                    textcoords="offset points", # how to position the text
                    xytext=(0,10), # distance from text to points (x,y)
                    ha='center') # horizontal alignment can be left, right or center

    # # Critical point of Lennard-Jones fluid
    # Tc = 1.321 # T*, T*=kT/eps これで温度からepsを算出できる。
    # Pc = 0.129 # P*, P* = P sig**3 / eps

    # cp = {"Acetone": {"Tc": 508, # K
    #                   "Pc": 48*101326 #Pa
    #                   },
    #       "DME": {"Tc": 402, # K
    #               "Pc": 52.6*101326 #Pa
    #               },
    #       "EO": {"Tc": 469, # K
    #              "Pc": 72.3*101326 #Pa
    #              },}

    # for guest, spec in cp.items():
    #     tc, pc = spec["Tc"], spec["Pc"]
    #     eps = NkB*tc / Tc # in kJ/mol
    #     sig = (Pc * eps*1e3/NA / pc)**(1/3)*1e10  # epsをkJ/molからJ / moleculeに換算。これはそれらしい数値。
    #     print(guest, sig, eps)
    #     plt.plot(sig, eps, "ok")
    #     plt.annotate(guest, # this is the text
    #                 (sig, eps), # these are the coordinates to position the label
    #                 textcoords="offset points", # how to position the text
    #                 xytext=(0,10), # distance from text to points (x,y)
    #                 ha='center') # horizontal alignment can be left, right or center

    gases = [AttrDict({"sig":4.438, # AA
                       "eps":488, # K
                       "name":"CS2",}),
             AttrDict({"sig":4.997, # AA
                       "eps":410, # K
                       "name":"n-C4H10",}),
             AttrDict({"sig":4.70, # AA
                       "eps":152.5, # K
                       "name":"CF4",}),
             AttrDict({"sig":4.59, # AA
                       "eps":189, # K
                       "name":"N2O",}),
             AttrDict({"sig":5.061, # AA
                       "eps":254, # K
                       "name":"Propane",}),
    ]
# CF4 4.70 152.5
# N2O 4.59 189
    for gas in gases:
        sig = gas.sig
        eps = gas.eps # K
        name = gas.name
        ax.plot(sig, eps, "ok")
        ax.annotate(name, # this is the text
                    (sig, eps), # these are the coordinates to position the label
                    textcoords="offset points", # how to position the text
                    xytext=(0,10), # distance from text to points (x,y)
                    ha='center') # horizontal alignment can be left, right or center


fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8, 4), sharey=True, gridspec_kw={'wspace': 0})
# Meと混合すると、I-II-III-I転移する条件をさがす。(てあたりしだい)
markers = {1:"o", 2:"+", 3:"^"}
guest = "LJME____"
axes[0].set_xlabel(r"$\sigma_g / \mathrm\AA$")
axes[0].set_ylabel(r"$\epsilon_g / \mathrm K$")
axes[0].set_xlim(3.5, 5.0)
axes[0].set_ylim(120.0, 600)
xyticks = 40
# ticks = np.concatenate([np.arange(0, 0.01, 0.0001), np.arange(0.01, 1, 0.01)])
ticks = np.arange(0., 1, 0.001)
for sig in np.linspace(3.5, 5.0, xyticks):
    for eps in np.linspace(120.0, 600, xyticks):
        g1 = AttrDict({"sigma":(moldict[host].intr[0][1]+moldict[guest].intr[0][1])/2,
                       "epsilon":(moldict[host].intr[0][0]*moldict[guest].intr[0][0])**0.5,
                       "dimen":0})
        g2 = AttrDict({"sigma":(moldict[host].intr[0][1]+sig)/2,
                       "epsilon":(moldict[host].intr[0][0]*eps*8.314/1000)**0.5,
                       "dimen":0})
        phases, lastphase = DoubleClathrate(g1, g2, beta, pressure, ticks = ticks)
        if lastphase == "CS1":
            color = "lightgreen"
        elif lastphase == "CS2":
            color = "#4666ff"
        else:
            color = "brown"
        # I-III-Iの場合を区別して描く。
        if lastphase == "CS1" and "TS1" in phases and len(phases) == 2:
            axes[0].plot(sig, eps, color="orange", marker=".")
        else:
            axes[0].plot(sig, eps, color=color, marker=markers[len(phases)])
        print(sig, eps, phases, lastphase)

cities(ax=axes[0])

# Xeと混合すると、I-II-III-I転移する条件をさがす。(てあたりしだい)
markers = {1:"o", 2:"+", 3:"^"}
guest = "LJXE____"
axes[1].set_xlabel(r"$\sigma_g / \mathrm\AA$")
# axes[1].set_ylabel(r"\epsilon_g / K")
axes[1].set_xlim(4.6, 5.0)
axes[1].set_ylim(120.0, 600)
xyticks = 40 # 60
# ticks = np.concatenate([np.arange(0, 0.01, 0.0001), np.arange(0.01, 1, 0.01)])
ticks = np.arange(0., 1, 0.001)
for sig in np.linspace(4.6, 5.0, xyticks):
    for eps in np.linspace(120.0, 600., xyticks):
        g1 = AttrDict({"sigma":(moldict[host].intr[0][1]+moldict[guest].intr[0][1])/2,
                       "epsilon":(moldict[host].intr[0][0]*moldict[guest].intr[0][0])**0.5,
                       "dimen":0})
        g2 = AttrDict({"sigma":(moldict[host].intr[0][1]+sig)/2,
                       "epsilon":(moldict[host].intr[0][0]*eps*8.314/1000)**0.5,
                       "dimen":0})
        phases, lastphase = DoubleClathrate(g1, g2, beta, pressure, ticks = ticks)
        if lastphase == "CS1":
            color = "lightgreen"
        elif lastphase == "CS2":
            color = "#4666ff"
        else:
            color = "brown"
        # I-III-Iの場合を区別して描く。
        if lastphase == "CS1" and "TS1" in phases and len(phases) == 2:
            axes[1].plot(sig, eps, color="orange", marker=".")
        else:
            axes[1].plot(sig, eps, color=color, marker=markers[len(phases)])
        print(sig, eps, phases, lastphase)

cities(ax=axes[1])
plt.tight_layout()

plt.show()
fig.savefig("s-e.pdf")
