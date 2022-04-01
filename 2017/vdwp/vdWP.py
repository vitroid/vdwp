# coding: utf-8
import CageIntegral.physconst as pc
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

    if type(f_c) == dict:
        f_c = [f_c]
        mu_guest = [mu_guest]

    for structure in structures:
        # number of water molecules and cages in a unit cell
        Nw, cages = crystals.cage_components[structure]
        # mu_f ######################
        sum = np.zeros_like(temperatures)
        beta = 1.0 / (pc.NkB * temperatures)
        for cage in cages:
            alpha = cages[cage] / Nw
            coeffs = -pc.NkB * temperatures * alpha
            sum0 = 1.0
            for f, mu in zip(f_c, mu_guest):
                # assume single occupation
                if cage in f:
                    sum0 += np.exp(beta * (mu - f[cage]))
            sum = sum + coeffs * np.log(sum0)
        Deltamu[structure] = sum

    return Deltamu
