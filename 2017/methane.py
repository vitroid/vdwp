#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# 2022-03-24
# 一般化相図上に、メタンを位置付ける。
# 温度や圧力を変えた時に、点がどう動くかを見せる。

import histo2f
import histo
import derivative
import physconst as pc
import math
import numpy as np
import crystals
# import sys
import interpolate as ip
import matplotlib.pyplot as plt
# import os
import chempot
import vdWP
from logging import getLogger
logger = getLogger()


def elimempty(values):
  newvalues = []
  for value in values:
    if value != "":
      newvalues.append(value)
  return newvalues


    
    

def prod(x):
    p = 1
    for v in x:
        p *= v
    return p

#! Test case for methane hydrate
import molecule
def LoadMoleculeDict(filename):
    file = open(filename,encoding="utf-8")
    moldict = molecule.loadInfo( file )
    #Supplementary info
    for key,mol in moldict.items():
        if sum(mol.moi) == 0:
            mol.symm = 1
            mol.dimen = 0
            mol.dof = 3
        elif prod(mol.moi) == 0:
            mol.symm = 2      #not always correct
            mol.dimen = 1
            mol.dof = 5
        else:
            mol.symm = 2      #not always correct
            mol.dimen = 3
            mol.dof = 6
        mol.name = key
    #Aliases
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
    return moldict


moldict = LoadMoleculeDict("data/DEFR")
#
#変数の命名則
#f, mu: free energy/chemical potential
#h    : enthalpy
#
#_w   : of water
#_i   : of ice
#_e   : of the empty lattice
#_f   : of the filled lattice
#_g   : of the guest
#_c   : in the cage
#

#User variables
pressure     = 101325.00 * 1# Pa
temperatures = np.array(np.arange(262,295,2))
#temperatures = np.array(np.arange(273,275,1))
mode = "Graphs"
# mode = "None"
guest = "LJME____"
host  = "TIP4P   "
structure = "sII"
debug = True

#分解状態での化学ポテンシャル(あるいは蒸気圧)
status = {"CJTHF___" : ("LIQUID-ANTOINE", 4.12118, 1202.942, -46.818),
          "NEGTHF__" : ("LIQUID-ANTOINE", 4.12118, 1202.942, -46.818),
          #"CPENTANE" : ("LIQUID-ANTOINE", 4.24714, 1235.305, -30.666),
          "CPENTANE" : ("LIQUID-LOGP"   ,   #進捗6より。
                        [ 252.1, 260., 270., 280., 290., 300. ],
                        [-2.49764352, -2.18736956, -1.72377413, -1.29552641, -0.80673168, -0.46915021]),
          "TIP4P   " : ("LIQUID-LOGP"   ,
                        [ 252.1, 260., 270., 280., 290., 300. ],
                        [-6.41141701, -5.95776395, -4.80618801, -4.00800805, -2.95803771, -2.85252098]),
          "CPENTAN+" : ("LIQUID-ANTOINE", 4.24714, 1235.305, -30.666),
          "CPENTAN-" : ("LIQUID-ANTOINE", 4.24714, 1235.305, -30.666),
          "CPENTA++" : ("LIQUID-ANTOINE", 4.24714, 1235.305, -30.666),
          "CPENTA--" : ("LIQUID-ANTOINE", 4.24714, 1235.305, -30.666),
          "LJME____" : ("IDEALGAS", pressure)}


    
#Winget's formula
def VaporPressure(temp, FE, dens, mass):
    #std pressure
    p0 = temp * 22.4 / 273.15
    pvap = p0*dens/mass*np.exp(FE/(0.008314*temp))
    return pvap

    

#generate the list of guests, structures, and cages
cages = crystals.cage_components[structure][1]


#gas terms (methane)
mol = moldict[guest]
logger.info(f"Calculating chemical potential of guest {mol.name}...")
logger.info(f"Regard the guest in gas state at {pressure} Pa.")
pressures = pressure

# mu_g ######################(15)
mu_g = chempot.chempot(temperatures, pressures) + chempot.IntegrationFixMinus(temperatures, mol.dimen) + chempot.StericFix(temperatures, mol.mass, mol.symm, mol.moi[0], mol.moi[1], mol.moi[2])
# h_g ###################### (28)
# h_g = dict()
mol = moldict[guest]
# logger.info("Calculating enthalpy of guest "+mol.name+"...")
# h_g = mu_g - temperatures* np.array(derivative.derivatives(temperatures,mu_g))


####### Clathrate terms ################################################

#for hydrate structure types
logger.info(f"Calculating chemical potential of empty clathrate {structure}...")
mu_e = crystals.U_e[structure] + vdWP.FreeEnergyOfVibration(crystals.nma_file[structure], temperatures)

#for gases in the cage
f_c = dict()
# h_k = dict()
for cage in cages:
    mol = moldict[guest]

    histofile   = "data/" + guest + "." + ("%d" % cage) + "hedra.histo"
    histogram   = histo.loadAHisto( open(histofile) )

    # f_c ######################
    if histogram != None:
        logger.info(f"Calculating free energy of guest {mol.name} in a {cage} hedral cage...")

        f_c[cage] = histo2f.fvalue(histogram, temperatures) + chempot.StericFix(temperatures, mol.mass, mol.symm, mol.moi[0], mol.moi[1], mol.moi[2])
        print(cage,":",f_c[cage])

#for gases in the lattice
mu_f = dict()
Deltamu = dict()
mol = moldict[guest]
logger.info(f"Calculating chemical potential of guest {mol.name} in {structure}...")

Nw        = crystals.cage_components[structure][0]

# mu_f ######################
sum = np.zeros_like(temperatures)
beta = 1.0/(pc.NkB*temperatures)
for cage in f_c:
    coeffs = -pc.NkB*temperatures*cages[cage]/Nw
    body  = np.log(1+np.exp(beta*(mu_g-f_c[cage])))
    add = coeffs*body
    sum = sum + add
Deltamu = sum
mu_f = mu_e + sum
# h_f ###################### (49)
# logger.info("Calculating enthalpy of guest "+mol.name+" in "+structure+"...")
# h_f = mu_f - temperatures* derivative.derivatives(temperatures,mu_f)

beta = 1.0/(pc.NkB*temperatures)
num_guest = np.zeros_like(temperatures)
for cage in f_c:
    e = np.exp(beta*(mu_g - f_c[cage]))
    y = e/(1.0+e)
#    print(y)
#    print(mu_g)
#    print(f_c[cage])
    #logger.info("occupancy in {0}:{1}".format(cage, y))
    num_guest = num_guest + y * cages[cage]/Nw
if num_guest[0] > 1e-4:
    if debug:
        print("Num of guest per a water molecule:", num_guest)






