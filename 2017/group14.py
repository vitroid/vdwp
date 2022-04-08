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




def drawLine(A, B, C, ax):
    """
    Ax + By + C = 0
    """
    if A == 0 and B == 0:
        return
    if A == 0:
        X = np.linspace(-0.4, 0.6, 100)
        Y = np.zeros_like(X) - C / B
        ax.plot(X, Y)
    elif B == 0:
        Y = np.linspace(-0.4, 0.6, 100)
        X = np.zeros_like(Y) - C / A
        ax.plot(X, Y)
    else:
        X = np.linspace(-0.4, 0.6, 100)
        Y = (-C - A * X) / B
        ax.plot(X, Y)


# Karttunen, A. J., Fässler, T. F., Linnolahti, M. & Pakkanen, T. A. Structural principles of semiconducting Group 14 clathrate frameworks. Inorg. Chem. 50, 1733–1742 (2011)
# Table 2
elements = ["C", "Si", "Ge", "Sn"]
karttunen = {element:dict() for element in elements}
density   = {element:dict() for element in elements}
for line in """
I       0.15 0.10 0.05 0.04   6.2 2.8 2.4 1.8 3.08 2.00 4.66 5.02
II      0.11 0.08 0.04 0.03   6.0 2.8 1.9 1.4 3.06 1.98 4.61 4.96
III     0.15 0.10 0.06 0.04   6.1 2.7 2.2 1.6 3.07 2.00 4.64 4.99
IV      0.17 0.12 0.07 0.05   6.1 2.7 2.0 1.4 3.05 1.99 4.61 4.97
V       0.12 0.09 0.04 0.03   6.0 2.8 1.9 1.4 3.05 1.98 4.61 4.96
II-4H   0.12 0.09 0.04 0.03   6.0 2.8 1.9 1.4 3.05 1.98 4.61 4.96
II+IV-a 0.20 0.13 0.08 0.06   5.9 2.6 1.5 1.0 3.04 1.98 4.60 4.96
II+IV-b 0.15 0.10 0.06 0.04   6.0 2.7 2.0 1.5 3.05 1.98 4.61 4.96
""".splitlines():
    cols = line.rstrip().split()
    if len(cols) > 0:
        phase = cols.pop(0)
        for i, element in enumerate(elements):
            karttunen[element][phase] = float(cols[i])
            density[element][phase] = float(cols[i+8])


print(karttunen)

ratios = {"I": (1.0, 0.0, 0.0),
          "II": (0.0, 1.0, 0.0),
          "IV": (0.0, 0.0, 1.0),
          "III": (23/43, 0.0, 20/43),
          "V": (0.0, 1.0, 0.0),
          "II-4H": (0.0, 1.0, 0.0),
          "II+IV-a": (0, 0.739130435, 0.260869565),
          "II+IV-b": (0, 0.459459459, 0.540540541),
          }

fig = plt.figure(figsize=(7,7))
gs = fig.add_gridspec(2, 2, hspace=0, wspace=0)
(ax1, ax2), (ax3, ax4) = gs.subplots(sharex='col', sharey='row')
# fig.suptitle('Sharing x per column, y per row')
# ax1.plot(x, y)
# ax2.plot(x, y**2, 'tab:orange')
# ax3.plot(x + 1, -y, 'tab:green')
# ax4.plot(x + 2, -y**2, 'tab:red')

# for ax in axs.flat:
#     ax.label_outer()


for element, panel in zip(elements, (ax1, ax2, ax3, ax4)):
    structures = karttunen[element].keys()
    colors = ("black", "red", "blue", "green", "pink", "orange", "purple", "brown")
    for x in np.linspace(-0.02, 0.06, 40):
        for y in np.linspace(0.0, 0.08, 40):
            smin = -1
            emin = 1e10
            for i, s in enumerate(structures):
                e = karttunen[element][s] + ratios[s][0]*x + ratios[s][1]*y
                if e < emin:
                    emin = e
                    smin = i
            panel.plot(x,y, "o", color=colors[smin])

for element, panel in zip(elements, (ax1, ax2, ax3, ax4)):
    panel.annotate(element, # this is the text
                   (0.8, 0.2), # these are the coordinates to position the label
                   xycoords="subfigure fraction", # how to position the text
                   ha='center') # horizontal alignment can be left, right or center
    structures = karttunen[element].keys()
    for s1, s2 in combinations(structures, 2):
        panel.set_xlim(-0.02, 0.06)
        panel.set_ylim(-0.0, 0.08)
        panel.set_xlabel(r"\Delta\mu_c^{I} - \Delta\mu_c^{IV}")
        panel.set_ylabel(r"\Delta\mu_c^{II} - \Delta\mu_c^{IV}")
        # 共存線の方程式
        # Ax + By + C = 0
        # A = xA - xB
        # B = yA - yB
        # C = mu_eA- mu_eB
        A = ratios[s1][0] - ratios[s2][0]
        B = ratios[s1][1] - ratios[s2][1]
        C = karttunen[element][s1] - karttunen[element][s2]
        drawLine(A, B, C, ax=panel)



plt.show()
# fig.savefig("group14.pdf")

fig = plt.figure(figsize=(7,7))
gs = fig.add_gridspec(2, 2, hspace=0, wspace=0)
(ax1, ax2), (ax3, ax4) = gs.subplots() #sharex='col', sharey='row')
for element, panel in zip(elements, (ax1, ax2, ax3, ax4)):
    for s in ratios:
        v0 = 1.0 / density[element][s]
        v1 = ratios[s][0] / density[element]["I"] + ratios[s][1] / density[element]["II"] + ratios[s][2] / density[element]["IV"]
        panel.plot(v0, v1, "o")

plt.show()
