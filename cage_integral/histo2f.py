#!/usr/bin/env python

import sys
import numpy as np
import vdwp.physconst as pc


def fvalue(histo: dict, temperature: np.ndarray):
    Rgas = 0.0083144621  # kJ/mol
    beta = 1.0 / (Rgas * temperature)
    pf = partition_function(histo, temperature)
    return (-np.log(pf)) / beta


def partition_function(histo: dict, temperature: np.ndarray):
    Rgas = 0.0083144621  # kJ/mol
    beta = 1.0 / (Rgas * temperature)
    pf = 0.0
    for energy in histo:
        weight = histo[energy]
        pf += weight * np.exp(-beta * energy)
    return pf


def energy(histo: dict, temperature: np.ndarray):
    Rgas = 0.0083144621  # kJ/mol
    beta = 1.0 / (Rgas * temperature)
    num = 0.0
    for energy in histo:
        weight = histo[energy]
        num += weight * np.exp(-beta * energy) * energy
    return num / partition_function(histo, temperature)


def test():
    from loader import loadAHisto

    temp = 273.15
    if len(sys.argv) > 1:
        if sys.argv[1] == "-t":
            temp = float(sys.argv[2])
    histo = loadAHisto(sys.stdin)
    if histo != None:
        print(partition_function(histo, temp), fvalue(histo, temp), energy(histo, temp))


if __name__ == "__main__":
    test()
