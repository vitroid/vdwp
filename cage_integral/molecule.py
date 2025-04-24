#!/usr/bin/env python

# handle @DEFR and @DEFP

import numpy as np


def loadDEFR(file):
    line = file.readline()
    id = line[0:8]
    line = file.readline()
    nsite = int(line)
    sites = []
    intr = []
    for site in range(nsite):
        line = file.readline()
        columns = line.split()  # x,y,z,mass,label
        columns[0:4] = map(float, columns[0:4])
        sites.append(columns)
    for site in range(nsite):
        line = file.readline()
        columns = line.split()  # eps, sig, charge
        columns = map(float, columns[0:3])
        intr.append(columns)
    return id, sites, intr


def loadDEFP(file):
    line = file.readline()
    id = line[0:8]
    line = file.readline()
    nsite = int(line)
    sites = []
    intr = []
    for site in range(nsite):
        line = file.readline()
        columns = line.split()  # mass,label
        columns = [0.0, 0.0, 0.0, float(columns[0]), columns[1]]
        sites.append(columns)
    for site in range(nsite):
        line = file.readline()
        columns = line.split()  # eps, sig, charge
        columns = map(float, columns[0:3])
        intr.append(columns)
    return id, sites, intr


class Molecule:
    def __init__(self):
        pass

    def load(self, file, tag):
        if tag == "@DEFR":
            self.id, self.sites, self.intr = loadDEFR(file)
        elif tag == "@DEFP":
            self.id, self.sites, self.intr = loadDEFP(file)
        mass = 0.0
        moi = [0.0, 0.0, 0.0]
        for site in self.sites:
            mass += site[3]
            moi[0] += site[3] * (site[1] ** 2 + site[2] ** 2)
            moi[1] += site[3] * (site[2] ** 2 + site[0] ** 2)
            moi[2] += site[3] * (site[0] ** 2 + site[1] ** 2)
        self.mass = mass
        self.moi = moi
        return self.id


def loadInfo(file):
    moldict = dict()
    while True:
        line = file.readline()
        if len(line) == 0:
            break
        columns = line.split()
        if len(columns) > 0 and columns[0] in ("@DEFR", "@DEFP"):
            tag = columns[0]
            mol = Molecule()
            id = mol.load(file, tag)
            moldict[id] = mol
    # symmetry
    for key, mol in moldict.items():
        if sum(mol.moi) == 0:
            mol.symm = 1
            mol.dimen = 0
            mol.dof = 3
        elif np.prod(mol.moi) == 0:
            mol.symm = 2  # not always correct
            mol.dimen = 1
            mol.dof = 5
        else:
            mol.symm = 2  # not always correct
            mol.dimen = 3
            mol.dof = 6
        mol.name = key
    # Aliases
    aliases = {
        "NEGTHF__": "THF(invert)",
        "CJTHF___": "THF",
        "CPENTANE": "cPentane",
        "CPENTAN+": "cPentane (105%)",
        "CPENTAN-": "cPentane (95%)",
        "CPEN+5L_": "cPentane (105% L)",
        "CPEN-5L_": "cPentane (95% L)",
        "CPENTA++": "cPentane (110%)",
        "CPENTA--": "cPentane (90%)",
        "LJME____": "Methane",
        "TYTHF___": "THF(Yaga)",
    }
    for key, mol in moldict.items():
        if key in aliases:
            mol.name = aliases[key]
    return moldict
