# coding: utf-8
import physconst as pc
import numpy as np

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
