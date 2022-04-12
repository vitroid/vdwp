# coding: utf-8
import vdwp.physconst as pc
import numpy as np


# Free Energy of Vibration / kJ mol^-1
def FreeEnergyOfVibration(nmafile, Temp):
    import os
    import os.path
    sum = 0.0
    wsum = 0.0
    dof = 6.0
    file = open(nmafile, "r")
    Beta = 1.0 / (pc.kB * Temp)
    while True:
        line = file.readline()
        if line == "":
            break
        (bin, weight) = line.split()
        omega = float(bin)  # +0.5
        weight = float(weight)
        # if omega > 1.0:
        sum += weight * pc.NkB * Temp * np.log(Beta * pc.h * omega * pc.cc)
        wsum += weight
    f = sum * dof / wsum
    return f
