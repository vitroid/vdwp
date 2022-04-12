#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# 2022-04-01
# ケージの球殻近似

import numpy as np


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


def cageradii():
    radii = dict()
    for size in (12, 14, 15, 16):
        with open(f"data/{size}hedra.nx4a") as f:
            box, cage = LoadAR3A(f)
        cage /= box
        cage -= np.floor(cage + 0.5)
        cage *= box
        cage = cage[np.linalg.norm(cage, axis=1) < 6.0]
        cage -= np.average(cage, axis=0)
        radii[size] = np.average(np.linalg.norm(cage, axis=1))
    return radii


if __name__ == "__main__":
    print(cageradii())
