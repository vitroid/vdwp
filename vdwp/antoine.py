#!/usr/bin/env python
# coding: utf-8

from math import *
import physconst as pc
from vdwp.decorators import deprecated_alias

# use ABC values in NIST database
# http://webbook.nist.gov/cgi/cbook.cgi?ID=C287923&Mask=4&Type=ANTOINE&Plot=on


def calculate_vapor_pressure(T, A, B, C):  # in Pa
    """Calculate vapor pressure using Antoine equation.

    Args:
        T: Temperature in Kelvin
        A, B, C: Antoine equation coefficients
    Returns:
        Vapor pressure in Pascal
    """
    return 10.0 ** (A - B / (C + T)) * 101325.0


@deprecated_alias("calculate_vapor_pressure")
def VaporPressure(T, A, B, C):
    pass


# for Relative value. See Antoine.pages
def calculate_chemical_potential_for_polyatomic(T, A, B, C):
    """Calculate chemical potential for polyatomic molecules using Antoine equation.

    Args:
        T: Temperature in Kelvin
        A, B, C: Antoine equation coefficients
    Returns:
        Chemical potential in J/mol
    """
    return -pc.NkB * T * (4 * log(T) + (-A + B / (C + T)) * log(10))


@deprecated_alias("calculate_chemical_potential_for_polyatomic")
def ChemicalPotentialPolyatomic(T, A, B, C):
    pass


def calculate_chemical_potential_of_water(T):
    """Calculate chemical potential of water using Antoine equation.

    Args:
        T: Temperature in Kelvin
    Returns:
        Chemical potential of water in J/mol
    """
    # from NIST webbook, correct between 273..303K
    A = 5.40221
    B = 1838.675
    C = -31.737
    # from NIST webbook, correct between 255.9..373K
    # A = 4.6543
    # B = 1435.264
    # C = -64.848
    return calculate_chemical_potential_for_polyatomic(T, A, B, C)


@deprecated_alias("calculate_chemical_potential_of_water")
def ChemicalPotentialOfWater(T):
    pass
