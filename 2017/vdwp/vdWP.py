# coding: utf-8
import physconst as pc
import numpy as np
from logging import getLogger
from functools import lru_cache
import CageIntegral.histo as histo
import CageIntegral.histo2f as histo2f
import vdwp.crystals as crystals

#占有率を求める。
# mu: chemical potential of the gas
# f:  list of the free energy of cage occupation
def Occupancy( mu, f, Temp ):
    Beta = 1.0 / ( pc.NkB * Temp )
    occ = dict()
    for i in f.keys():
        occ[i] = np.exp(Beta*(mu - f[i])) / (1.0 + np.exp(Beta*(mu - f[i])))
    return occ


#Free Energy of Vibration / kJ mol^-1
def FreeEnergyOfVibration(nmafile, Temp):
    import os
    import os.path
    sum=0.0
    wsum = 0.0
    dof = 6.0
    file = open(nmafile,"r")
    Beta = 1.0 / ( pc.kB * Temp )
    while True:
        line = file.readline()
        if line == "":
            break
        (bin,weight) = line.split()
        omega = float(bin) #+0.5
        weight = float(weight)
        #if omega > 1.0:
        sum  += weight * pc.NkB * Temp * np.log(Beta * pc.h * omega * pc.cc)
        wsum += weight
    f = sum * dof / wsum
    return f

# @lru_cache
def EncagingFE(temperatures, guest, stericterm):
    logger = getLogger()
    f_c = dict()
    for cage in (12, 14, 15, 16):
        histofile = f"data/{guest}.ice{cage}.histo"
        histogram = histo.loadAHisto(open(histofile))

        # f_c ######################
        if histogram is not None:

            f_c[cage] = (histo2f.fvalue(histogram, temperatures) + stericterm)
            logger.debug(f"{cage}: {f_c[cage]}")
    return f_c


def ChemPotByOccupation(temperatures, f_c, mu_guest, structures):
    logger = getLogger()
    Deltamu = dict()

    for structure in structures:
        # number of water molecules and cages in a unit cell
        Nw, cages = crystals.cage_components[structure]
        # mu_f ######################
        sum = np.zeros_like(temperatures)
        beta = 1.0 / (pc.NkB * temperatures)
        for cage in cages:
            if cage in f_c:
                coeffs = -pc.NkB * temperatures * cages[cage] / Nw
                body = np.log(1 + np.exp(beta * (mu_guest - f_c[cage])))
                add = coeffs * body
                sum = sum + add
        Deltamu[structure] = sum

    return Deltamu
