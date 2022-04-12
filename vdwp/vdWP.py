# coding: utf-8
import vdwp.physconst as pc
import numpy as np
from logging import getLogger
import vdwp.crystals as crystals

# 占有率を求める。
# mu: chemical potential of the gas
# f:  list of the free energy of cage occupation


def Occupancy(mu, f, Temp):
    Beta = 1.0 / (pc.NkB * Temp)
    occ = dict()
    for i in f.keys():
        occ[i] = np.exp(Beta * (mu - f[i])) / \
            (1.0 + np.exp(Beta * (mu - f[i])))
    return occ


def ChemPotByOccupation(temperatures, f_c, mu_guest, structures):
    logger = getLogger()
    Deltamu = dict()

    if isinstance(f_c, dict):
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
