#!/usr/bin/env python

import vdwp.physconst as pc
import numpy as np
from logging import getLogger


# atomic mass in g/mol
# OBSOLETE
def _IdealGas(Temp, pressure, atomicmass):
    A = 2.0 * np.pi * (atomicmass * 0.001 / pc.NA) * pc.kB * Temp / (pc.h**2)
    return -pc.NkB * Temp * (1.5 * np.log(A) + np.log(pc.kB * Temp / pressure))


# a_coeff == coeff comes from mass(g/mol)
def Avalue(Temp, mass):
    return (
        -pc.NkB
        * Temp
        * (3.0 / 2.0)
        * np.log(mass * 0.001 / pc.NA * 2.0 * np.pi * pc.kB * Temp / pc.h**2)
    )


# ixx,iyy,izz: moment of inertia in (atomic mass g/mol) x A**2 (Use output of momentofinertia.py)
# c_coeff == coeff comes from moment of inertia
def Cvalue(T, ixx, iyy, izz):
    ixx *= 0.001 / pc.NA * 1e-10**2
    iyy *= 0.001 / pc.NA * 1e-10**2
    izz *= 0.001 / pc.NA * 1e-10**2
    return (
        -pc.NkB
        * Temp
        * (
            (1.0 / 2.0) * np.log(ixx * iyy * izz)
            + (3.0 / 2.0) * np.log(2.0 * np.pi * pc.kB * Temp / pc.h**2)
        )
    )


def Cvalue2(Temp, ixx):
    ixx *= 0.001 / pc.NA * 1e-10**2
    return -pc.NkB * Temp * (np.log(ixx) + np.log(2.0 * np.pi * pc.kB * Temp / pc.h**2))


def SymmetryFix(Temp, symm):
    return pc.NkB * Temp * np.log(symm)


# In case dimen == 3, even if you change it,
# both water and guest share the same value and therefore they cancel.
# Correction for empty integration.
def IntegrationFixMinus(Temp, dimen):
    # 2014-12-24 Removed
    # return 0.0
    # 2015-2-5 revival.
    if dimen == 3:
        return -pc.NkB * Temp * np.log(8 * np.pi**2)
    elif dimen == 1:
        return -pc.NkB * Temp * np.log(4 * np.pi)
    else:
        return 0.0


# ixx,iyy,izz: moment of inertia in (atomic mass) x A**2 (Use output of momentofinertia.py)
# For use with ideal gas EOS.
def StericFix(Temp, mass, symm, moi):
    logger = getLogger()
    ixx, iyy, izz = moi
    value = Avalue(Temp, mass)
    if ixx == 0:
        # monatomic
        dimen = 0
        logger.debug(f"{Temp}\tTemperature")
        logger.debug(f"{value}\tAvalue")
        return value
    elif izz == 0:
        # rodlike
        dimen = 1
        logger.debug(f"{Temp}\tTemperature")
        # logger.debug IntegrationFix(Temp,dimen), "\tIntegration Fix"
        # logger.debug SymmetryFix(Temp, symm), "\tSymmetryFix"
        logger.debug(f"{value}\tAvalue")
        logger.debug(SymmetryFix(Temp, symm) + Cvalue2(Temp, ixx), "\tCvalue2")
        # return value + IntegrationFix(Temp,dimen) + SymmetryFix(Temp, symm) + Cvalue2(Temp,ixx)
        return value + SymmetryFix(Temp, symm) + Cvalue2(Temp, ixx)
    else:
        # polyatomic
        dimen = 3
        logger.debug(f"{Temp}\tTemperature")
        # logger.debug IntegrationFix(Temp,dimen), "\tIntegration Fix"
        # logger.debug SymmetryFix(Temp, symm), "\tSymmetryFix"
        logger.debug(f"{value}\tAvalue")
        logger.debug(SymmetryFix(Temp, symm) + Cvalue(Temp, ixx, iyy, izz), "\tCvalue")
        # return value + IntegrationFix(Temp,dimen) + SymmetryFix(Temp, symm) + Cvalue(Temp,ixx, iyy, izz)
        return value + SymmetryFix(Temp, symm) + Cvalue(Temp, ixx, iyy, izz)


# #ixx,iyy,izz: moment of inertia in (atomic mass) x A**2 (Use output of momentofinertia.py)
# def StericFix2(T,mass,symm,ixx=0,iyy=0,izz=0,debug=False):
#     value = Avalue(T,mass)
#     if ixx == 0:
#         #monatomic
#         dimen = 0
#         if debug:
#             logger.debug value, "\tAvalue"
#         return value
#     elif izz == 0:
#         #rodlike
#         dimen = 1
#         if debug:
#             logger.debug value, "\tAvalue"
#             logger.debug Cvalue2(T,ixx), "\tCvalue2"
#         return value + Cvalue2(T,ixx)
#     else:
#         #polyatomic
#         dimen = 3
#         if debug:
#             logger.debug value, "\tAvalue"
#             logger.debug Cvalue(T,ixx, iyy, izz), "\tCvalue"
#         return value + Cvalue(T,ixx, iyy, izz)


# #for use with histogram (== result of integration)
# def StericFix3(T,mass,symm,ixx=0,iyy=0,izz=0,debug=False):
#     value = Avalue(T,mass)
#     if ixx == 0:
#         #monatomic
#         dimen = 0
#         if debug:
#             logger.debug value, "\tAvalue"
#         return value
#     elif izz == 0:
#         #rodlike
#         dimen = 1
#         if debug:
#             logger.debug SymmetryFix(T, symm), "\tSymmetryFix"
#             logger.debug value, "\tAvalue"
#             logger.debug Cvalue2(T,ixx), "\tCvalue2"
#         return value + SymmetryFix(T, symm) + Cvalue2(T,ixx)
#     else:
#         #polyatomic
#         dimen = 3
#         if debug:
#             logger.debug SymmetryFix(T, symm), "\tSymmetryFix"
#             logger.debug value, "\tAvalue"
#             logger.debug Cvalue(T,ixx, iyy, izz), "\tCvalue"
#         return value + SymmetryFix(T, symm) + Cvalue(T,ixx, iyy, izz)


def chempot(Temp, p):
    return -pc.NkB * Temp * (np.log(pc.kB * Temp / p))


def test():
    import antoine

    A = 5.40221
    B = 1838.675
    C = -31.737
    for T in range(273, 293):
        logger.debug(
            T,
            chempot(T, 101326) + Avalue(T, 16),
            chempot(T, antoine.VaporPressure(T, A, B, C)),
        )


if __name__ == "__main__":
    test()
