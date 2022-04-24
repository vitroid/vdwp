#!/usr/bin/env python
# coding: utf-8

from math import *
import physconst as pc

# use ABC values in NIST database
# http://webbook.nist.gov/cgi/cbook.cgi?ID=C287923&Mask=4&Type=ANTOINE&Plot=on


def VaporPressure(T, A, B, C):  # in Pa
    return 10.0**(A - B / (C + T)) * 101325.0

# for Relative value. See Antoine.pages


def ChemicalPotentialPolyatomic(T, A, B, C):
    return -pc.NkB * T * (4 * log(T) + (-A + B / (C + T)) * log(10))


def ChemicalPotentialOfWater(T):
    # from NIST webbook, correct between 273..303K
    A = 5.40221
    B = 1838.675
    C = -31.737
    # from NIST webbook, correct between 255.9..373K
    #A = 4.6543
    #B = 1435.264
    #C = -64.848
    return ChemicalPotentialPolyatomic(T, A, B, C)

