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



figure = plt.figure()


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
    moldict["EKABR2__"].name = "EkaBr2"
    return moldict

def LoadAR3A(file):
    tag = ""
    for line in iter(file.readline, ""):
        if line[0] == "@":
            tag = line.rstrip()
        elif tag == "@BOX3":
            box = [float(x) for x in line.split()]
            tag = ""
        elif tag == "@NX4A":
            nmol = int(line)
            tag = "+"
            atoms = []
        elif tag == "+":
            atoms.append([float(x) for x in line.split()[:3]])
            nmol -= 1
            if nmol == 0:
                return box, np.array(atoms)

with open("data/12hedra.nx4a") as f:
    box12, cage12 = LoadAR3A(f)

cage12 /= box12
cage12 -= np.floor(cage12+0.5)
cage12 *= box12
cage12 = cage12[np.linalg.norm(cage12, axis=1)<6.0]

with open("data/16hedra.nx4a") as f:
    box16, cage16 = LoadAR3A(f)

cage16 /= box16
cage16 -= np.floor(cage16+0.5)
cage16 *= box16
cage16 = cage16[np.linalg.norm(cage16, axis=1)<6.0]

with open("data/15hedra.nx4a") as f:
    box15, cage15 = LoadAR3A(f)

cage15 /= box15
cage15 -= np.floor(cage15+0.5)
cage15 *= box15
cage15 = cage15[np.linalg.norm(cage15, axis=1)<6.0]


moldict = LoadMoleculeDict("data/DEFR")

# User variables
guest = "LJME____"
host = "TIP4PICE"
structures = ["CS2", "CS1", "HS1", "TS1"]  # , "TS1", "HS1"]

temperatures = 273.15

mol = moldict[guest]
stericterm = chempot.StericFix(temperatures, mol.mass, mol.symm, mol.moi)
f_c = vdWP.EncagingFE(temperatures, guest, stericterm)
print(f_c)

sigma = (moldict[host].intr[0][1] + moldict[guest].intr[0][1])/2
epsilon = (moldict[host].intr[0][0] * moldict[guest].intr[0][0])**0.5 # in kJ/mol
beta = 1.0 / (pc.NkB * temperatures)

from ljd.ljd import fvalue

plt.xlabel("r / AA")
# plt.ylim(0, 10)
plt.ylim(-10, 0)
plt.ylabel("w(r)/kT")
print(beta)

R12 = 3.902
z12 = 20
f12 = fvalue(R12, sigma, epsilon, beta, z12)
C12 = f12 + stericterm
print(12, f12, C12)

R14 = 4.326
z14 = 24
f14 = fvalue(R14, sigma, epsilon, beta, z14)
C14 = f14 + stericterm
print(14, f14, C14)

R16 = 4.683
z16 = 28
f16 = fvalue(R16, sigma, epsilon, beta, z16)
C16 = f16 + stericterm
print(16, f16, C16)

# direct calculation of the potential
r = np.linspace(0, 2.0, 100)
x = np.zeros([100,3])
x[:,0] = r
ep = 0
for w in cage12:
    sr = sigma / np.linalg.norm(x - w, axis=1)
    ep += 4*epsilon*(sr**12 - sr**6)
plt.plot(r, ep*beta)

ep = 0
for w in cage16:
    sr = sigma / np.linalg.norm(x - w, axis=1)
    ep += 4*epsilon*(sr**12 - sr**6)
plt.plot(r, ep*beta)
# きれいに一致するので、cage potentialは正しい。なぜ積分が全然違う値になるのか。
# もしかして、EncagingFEのfvalueにはde Broglie wave lengthが含まれている?

plt.show()

