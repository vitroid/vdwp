#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# phasediagram.pyからの変更点:
# ケージ積分の代わりにljdを用いる。
# 必然的に、単原子モデル以外は扱えない。
# 相境界の条件は50気圧、273 Kのみにする。
# moldictの利用をやめる。

import vdwp.vdWP as vdWP
import vdwp.crystals as crystals
import vdwp.normalmode as normalmode
from vdwp.physconst import NkB, NA
# import CageIntegral.molecule as molecule
import vdwp.chempot as chempot
# from ljd.ljd import fvalue
from vdwp.general import drawLine

import numpy as np
# import vdwp.interpolate as ip
import matplotlib.pyplot as plt
from itertools import combinations
from logging import getLogger, DEBUG, INFO, basicConfig
# from attrdict import AttrDict

# basicConfig(level=DEBUG, format="%(levelname)s %(message)s")
basicConfig(level=INFO, format="%(levelname)s %(message)s")
logger = getLogger()
logger.debug("Debug mode.")


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


# plt.rcParams['text.usetex'] = True
plt.rcParams["font.size"] = 14
plt.rcParams["font.family"] = "sans-serif"


figure = plt.figure(figsize=(5, 5))


plt.xlim(-0.4, 0.6)
plt.ylim(-0.4, 0.6)
plt.xticks([-0.4, -0.2, 0, 0.2, 0.4, 0.6])
plt.xlabel(
    r"$(\Delta\mu_c^\mathrm{CS1} - \Delta\mu_c^\mathrm{HS1}) / \mathrm{kJ~mol}^{-1}$")
plt.ylabel(
    r"$(\Delta\mu_c^\mathrm{CS2} - \Delta\mu_c^\mathrm{HS1}) / \mathrm{kJ~mol}^{-1}$")
plt.axis("square")

temperatures = 273.15

####### chemical potential of gas ######################################
stericterm = chempot.StericFix(temperatures, mass=1.0, symm=1, moi=(0, 0, 0))

####### structure-dependent terms ######################################

# for hydrate structure types
mu_e = dict()
for structure in crystals.names:
    logger.info(
        f"Calculating chemical potential of empty clathrate {structure}...")
    mu_e[structure] = crystals.U_e[structure] + \
        normalmode.FreeEnergyOfVibration(crystals.nma_file[structure], temperatures)
# mu_e = crystals.mu_e

for s1, s2 in combinations(crystals.names, 2):
    # 共存線の方程式
    # Ax + By + C = 0
    # A = xA - xB
    # B = yA - yB
    # C = mu_eA- mu_eB
    A = crystals.ratios[s1][0] - crystals.ratios[s2][0]
    B = crystals.ratios[s1][1] - crystals.ratios[s2][1]
    C = mu_e[s1] - mu_e[s2]
    drawLine(A, B, C, style="-k")

# manually labelled
plt.annotate("I (CS-I)",  # this is the text
             xy=(0.2, 0.9),  # these are the coordinates to position the label
             xycoords="axes fraction",
             fontsize=24,
             ha="center")
plt.annotate("II (CS-II)",  # this is the text
             xy=(0.8, 0.4),  # these are the coordinates to position the label
             xycoords="axes fraction",
             fontsize=24,
             ha="center")
plt.annotate("III\n(TS-I)",  # this is the text
             xy=(0.57, 0.8),  # these are the coordinates to position the label
             xycoords="axes fraction",
             fontsize=24,
             ha="center")
plt.annotate("IV\n(HS-I)",  # this is the text
             xy=(0.85, 0.8),  # these are the coordinates to position the label
             xycoords="axes fraction",
             fontsize=24,
             ha="center")

plt.tight_layout()
plt.show()
figure.savefig("Figure1.pdf")
