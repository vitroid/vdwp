#!/usr/bin/env python

import physconst as pc
from math import *


#
# atomic mass in g/mol
# OBSOLETE
def IdealGas(Temp, pressure, atomicmass):
    A = 2.0 * pi * (atomicmass * 0.001 / pc.NA) * pc.kB * Temp / (pc.h**2)
    return -pc.NkB * Temp * (1.5 * log(A) + log(pc.kB * Temp / pressure))


# a_coeff == coeff comes from mass(g/mol)
def Avalue(T, mass):
    return (
        -pc.NkB
        * T
        * (3.0 / 2.0)
        * log(mass * 0.001 / pc.NA * 2.0 * pi * pc.kB * T / pc.h**2)
    )


# ixx,iyy,izz: moment of inertia in (atomic mass g/mol) x A**2 (Use output of momentofinertia.py)
# c_coeff == coeff comes from moment of inertia
def Cvalue(T, ixx, iyy, izz):
    ixx *= 0.001 / pc.NA * 1e-10**2
    iyy *= 0.001 / pc.NA * 1e-10**2
    izz *= 0.001 / pc.NA * 1e-10**2
    return (
        -pc.NkB
        * T
        * (
            (1.0 / 2.0) * log(ixx * iyy * izz)
            + (3.0 / 2.0) * log(2.0 * pi * pc.kB * T / pc.h**2)
        )
    )


def Cvalue2(T, ixx):
    ixx *= 0.001 / pc.NA * 1e-10**2
    return -pc.NkB * T * (log(ixx) + log(2.0 * pi * pc.kB * T / pc.h**2))


def SymmetryFix(T, symm):
    return pc.NkB * T * log(symm)


# In case dimen == 3, even if you change it,
# both water and guest share the same value and therefore they cancel.
# Correction for empty integration.
def IntegrationFixMinus(T, dimen):
    # 2014-12-24 Removed
    # return 0.0
    # 2015-2-5 revival.
    if dimen == 3:
        return -pc.NkB * T * log(8 * pi**2)
    elif dimen == 1:
        return -pc.NkB * T * log(4 * pi)
    else:
        return 0.0


# ixx,iyy,izz: moment of inertia in (atomic mass) x A**2 (Use output of momentofinertia.py)
# For use with ideal gas EOS.
def StericFix(T, mass, symm, ixx=0, iyy=0, izz=0, debug=False):
    value = Avalue(T, mass)
    if ixx == 0:
        # monatomic
        dimen = 0
        if debug:
            print(T, "\tTemperature")
            print(value, "\tAvalue")
        return value
    elif izz == 0:
        # rodlike
        dimen = 1
        if debug:
            print(T, "\tTemperature")
            # print IntegrationFix(T,dimen), "\tIntegration Fix"
            # print SymmetryFix(T, symm), "\tSymmetryFix"
            print(value, "\tAvalue")
            print(SymmetryFix(T, symm) + Cvalue2(T, ixx), "\tCvalue2")
        # return value + IntegrationFix(T,dimen) + SymmetryFix(T, symm) + Cvalue2(T,ixx)
        return value + SymmetryFix(T, symm) + Cvalue2(T, ixx)
    else:
        # polyatomic
        dimen = 3
        if debug:
            print(T, "\tTemperature")
            # print IntegrationFix(T,dimen), "\tIntegration Fix"
            # print SymmetryFix(T, symm), "\tSymmetryFix"
            print(value, "\tAvalue")
            print(SymmetryFix(T, symm) + Cvalue(T, ixx, iyy, izz), "\tCvalue")
        # return value + IntegrationFix(T,dimen) + SymmetryFix(T, symm) + Cvalue(T,ixx, iyy, izz)
        return value + SymmetryFix(T, symm) + Cvalue(T, ixx, iyy, izz)


# #ixx,iyy,izz: moment of inertia in (atomic mass) x A**2 (Use output of momentofinertia.py)
# def StericFix2(T,mass,symm,ixx=0,iyy=0,izz=0,debug=False):
#     value = Avalue(T,mass)
#     if ixx == 0:
#         #monatomic
#         dimen = 0
#         if debug:
#             print value, "\tAvalue"
#         return value
#     elif izz == 0:
#         #rodlike
#         dimen = 1
#         if debug:
#             print value, "\tAvalue"
#             print Cvalue2(T,ixx), "\tCvalue2"
#         return value + Cvalue2(T,ixx)
#     else:
#         #polyatomic
#         dimen = 3
#         if debug:
#             print value, "\tAvalue"
#             print Cvalue(T,ixx, iyy, izz), "\tCvalue"
#         return value + Cvalue(T,ixx, iyy, izz)


# #for use with histogram (== result of integration)
# def StericFix3(T,mass,symm,ixx=0,iyy=0,izz=0,debug=False):
#     value = Avalue(T,mass)
#     if ixx == 0:
#         #monatomic
#         dimen = 0
#         if debug:
#             print value, "\tAvalue"
#         return value
#     elif izz == 0:
#         #rodlike
#         dimen = 1
#         if debug:
#             print SymmetryFix(T, symm), "\tSymmetryFix"
#             print value, "\tAvalue"
#             print Cvalue2(T,ixx), "\tCvalue2"
#         return value + SymmetryFix(T, symm) + Cvalue2(T,ixx)
#     else:
#         #polyatomic
#         dimen = 3
#         if debug:
#             print SymmetryFix(T, symm), "\tSymmetryFix"
#             print value, "\tAvalue"
#             print Cvalue(T,ixx, iyy, izz), "\tCvalue"
#         return value + SymmetryFix(T, symm) + Cvalue(T,ixx, iyy, izz)


def chempot(T, p):
    return -pc.NkB * T * (log(pc.kB * T / p))


def test():
    import antoine

    A = 5.40221
    B = 1838.675
    C = -31.737
    for T in range(273, 293):
        print(
            T,
            chempot(T, 101326) + Avalue(T, 16),
            chempot(T, antoine.VaporPressure(T, A, B, C)),
        )


if __name__ == "__main__":
    test()
