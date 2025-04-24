#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import histo2f
import histo

# import math
import sys
import vdwp.chempot as chempot


def prod(x):
    p = 1
    for v in x:
        p *= v
    return p


import molecule


def LoadMoleculeDict(filename):
    file = open(filename, encoding="utf-8")
    moldict = molecule.loadInfo(file)
    # Supplementary info
    for key, mol in moldict.items():
        if sum(mol.moi) == 0:
            mol.symm = 1
            mol.dimen = 0
            mol.dof = 3
        elif prod(mol.moi) == 0:
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
    return moldict


moldict = LoadMoleculeDict("DEFR")
guest = sys.argv[1]
T = 273.15

for cage in 12, 14, 16:
    mol = moldict[guest]
    histofile = guest + "." + ("%d" % cage) + "hedra.histo"
    histogram = histo.loadAHisto(open(histofile))
    # f_c ######################
    if histogram != None:
        print(
            cage,
            histo2f.fvalue(histogram, T),
            histo2f.fvalue(histogram, T)
            + chempot.molecular_chemical_potential_corrections(
                T, mol.mass, mol.symm, mol.moi
            ),
        )


sys.exit(0)
