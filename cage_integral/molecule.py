#!/usr/bin/env python

# handle @DEFR and @DEFP
from logging import getLogger
import numpy as np
from pymatgen.symmetry.analyzer import PointGroupAnalyzer
from pymatgen.core import Molecule


def load_DEFR_and_DEFP(file, drop_virtual_sites=False):
    logger = getLogger()

    line = file.readline()
    id = line[0:8]
    line = file.readline()
    nsite = int(line)
    sites = []
    intr = []
    mass = []
    labels = []
    for site in range(nsite):
        line = file.readline()
        columns = line.split()  # x,y,z,mass,label
        if len(columns) == 5:
            if drop_virtual_sites and float(columns[3]) == 0.0:
                continue
            sites.append([float(x) for x in columns[:3]])
            mass.append(float(columns[3]))
            labels.append(columns[4])
        else:
            sites.append((0.0, 0.0, 0.0))
            mass.append(float(columns[0]))
            labels.append(columns[1])
    for site in range(nsite):
        line = file.readline()
        columns = line.split()  # eps, sig, charge
        columns = map(float, columns[0:3])
        intr.append(columns)
    molecule = dict()
    molecule["id"] = id
    molecule["sites"] = sites
    molecule["mass"] = mass
    molecule["labels"] = labels
    molecule["intr"] = intr
    return molecule


class RigidMolecule:
    def __init__(self, file, drop_virtual_sites=False):
        logger = getLogger()

        mol = load_DEFR_and_DEFP(file, drop_virtual_sites)
        self.id = mol["id"]
        self.sites = mol["sites"]
        self.site_mass = mol["mass"]
        self.labels = mol["labels"]
        self.intr = mol["intr"]
        inertia_moment = [0.0, 0.0, 0.0]
        # 慣性モーメントは対角化され、重心が原点にあることを仮定する。
        for site, site_mass in zip(self.sites, self.site_mass):
            inertia_moment[0] += site_mass * (site[1] ** 2 + site[2] ** 2)
            inertia_moment[1] += site_mass * (site[2] ** 2 + site[0] ** 2)
            inertia_moment[2] += site_mass * (site[0] ** 2 + site[1] ** 2)
        self.mass = sum(self.site_mass)
        self.moi = inertia_moment

        # monatomic
        if sum(self.moi) == 0:
            self.dimen = 0
            self.dof = 3
            self.symm = 1
            return

        # linear
        if np.prod(self.moi) == 0:
            self.dimen = 1
            self.dof = 5
        # nonlinear
        else:
            self.dimen = 3
            self.dof = 6

        molecule = Molecule(self.labels, self.sites)
        pga = PointGroupAnalyzer(molecule)
        point_group = pga.get_pointgroup()
        if point_group.sch_symbol in ("C2v", "D*h", "C*v", "Cs"):
            self.symm = 2
        elif point_group.sch_symbol in ("D5h",):
            self.symm = 10
        elif point_group.sch_symbol in ("D3h",):
            self.symm = 6
        elif point_group.sch_symbol in ("C3v",):
            self.symm = 3
        else:
            logger.error(f"Invalid point group: {point_group.sch_symbol}")
        logger.debug([self.id, point_group.sch_symbol, self.symm])


def loadInfo(file):
    moldict = dict()
    while True:
        line = file.readline()
        if len(line) == 0:
            break
        columns = line.split()
        if len(columns) > 0 and columns[0] in ("@DEFR", "@DEFP"):
            tag = columns[0]
            mol = RigidMolecule(file, drop_virtual_sites=True)
            moldict[mol.id] = mol
    # symmetry
    for key, mol in moldict.items():
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


def test():
    with open("DEFR", "r") as f:
        info = loadInfo(f)


if __name__ == "__main__":
    test()
