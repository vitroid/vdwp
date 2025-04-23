#!/usr/bin/env python

import vdwp.physconst as pc
import numpy as np
from logging import getLogger, DEBUG, basicConfig
from vdwp.decorators import deprecated_alias


# atomic mass in g/mol
# OBSOLETE
# def _IdealGas(Temp, pressure, atomicmass):
#     A = 2.0 * np.pi * (atomicmass * 0.001 / pc.NA) * pc.kB * Temp / (pc.h**2)
#     return -pc.NkB * Temp * (1.5 * np.log(A) + np.log(pc.kB * Temp / pressure))


def translational_contribution_to_chemical_potential(temperature, mass):
    """Calculate the translational contribution to chemical potential for an ideal gas.

    Args:
        temperature: Temperature in Kelvin
        mass: Molecular mass in g/mol
    """
    return (
        -pc.NkB
        * temperature
        * (3.0 / 2.0)
        * np.log(mass * 0.001 / pc.NA * 2.0 * np.pi * pc.kB * temperature / pc.h**2)
    )


@deprecated_alias("translational_contribution_to_chemical_potential")
def Avalue(temperature, mass):
    pass


# ixx,iyy,izz: moment of inertia in (atomic mass g/mol) x A**2 (Use output of momentofinertia.py)
# c_coeff == coeff comes from moment of inertia


def rotational_contribution_to_chemical_potential(temperature, ixx, iyy, izz):
    """Calculate the rotational contribution to chemical potential for a non-linear molecule.

    Args:
        temperature: Temperature in Kelvin
        ixx, iyy, izz: Principal moments of inertia in (atomic mass g/mol) x A**2
    """
    ixx *= 0.001 / pc.NA * 1e-10**2
    iyy *= 0.001 / pc.NA * 1e-10**2
    izz *= 0.001 / pc.NA * 1e-10**2
    return (
        -pc.NkB
        * temperature
        * (
            (1.0 / 2.0) * np.log(ixx * iyy * izz)
            + (3.0 / 2.0) * np.log(2.0 * np.pi * pc.kB * temperature / pc.h**2)
        )
    )


@deprecated_alias("rotational_contribution_to_chemical_potential")
def Cvalue(temperature, ixx, iyy, izz):
    pass


def linear_molecule_rotational_contribution(temperature, ixx):
    """Calculate the rotational contribution to chemical potential for a linear molecule.

    Args:
        temperature: Temperature in Kelvin
        ixx: Moment of inertia in (atomic mass g/mol) x A**2
    """
    ixx *= 0.001 / pc.NA * 1e-10**2
    return (
        -pc.NkB
        * temperature
        * (np.log(ixx) + np.log(2.0 * np.pi * pc.kB * temperature / pc.h**2))
    )


@deprecated_alias("linear_molecule_rotational_contribution")
def Cvalue2(temperature, ixx):
    pass


# suggested new name: contribution_from_symmetry
def symmetry_correction_to_chemical_potential(temperature, symmetry_number):
    """Calculate the symmetry correction to chemical potential.

    Args:
        temperature: Temperature in Kelvin
        symmetry_number: Number of indistinguishable orientations
    """
    return pc.NkB * temperature * np.log(symmetry_number)


@deprecated_alias("symmetry_correction_to_chemical_potential")
def SymmetryFix(temperature, symmetry):
    pass


# In case dimen == 3, even if you change it,
# both water and guest share the same value and therefore they cancel.
# Correction for empty integration.
def phase_space_integration_correction(temperature, degrees_of_freedom):
    """Calculate the correction for phase space integration.

    Args:
        temperature: Temperature in Kelvin
        degrees_of_freedom: Number of rotational degrees of freedom (1 for linear, 3 for non-linear)
    """
    if degrees_of_freedom == 3:
        return -pc.NkB * temperature * np.log(8 * np.pi**2)
    elif degrees_of_freedom == 1:
        return -pc.NkB * temperature * np.log(4 * np.pi)
    else:
        return 0.0


@deprecated_alias("phase_space_integration_correction")
def IntegrationFixMinus(temperature, dimension):
    pass


# ixx,iyy,izz: moment of inertia in (atomic mass) x A**2 (Use output of momentofinertia.py)
# For use with ideal gas EOS.
# suggested_new_name: corrections_due_to_molecular_shape
def molecular_chemical_potential_corrections(
    temperature, mass, symmetry_number, moment_of_inertia, debug=False
):
    """Calculate all corrections to the chemical potential for a molecule.

    Args:
        temperature: Temperature in Kelvin
        mass: Molecular mass in g/mol
        symmetry_number: Number of indistinguishable orientations
        moment_of_inertia: Tuple of (ixx, iyy, izz) in (atomic mass g/mol) x A**2
        debug: Whether to print debug information
    """
    logger = getLogger()
    ixx, iyy, izz = moment_of_inertia
    value = translational_contribution_to_chemical_potential(temperature, mass)
    if ixx == 0:
        # monatomic
        degrees_of_freedom = 0
        if debug:
            logger.debug(f"{temperature}\tTemperature")
            logger.debug(f"{value}\tTranslational contribution")
        return value
    elif izz == 0:
        # rodlike
        degrees_of_freedom = 1
        if debug:
            logger.debug(f"{temperature}\tTemperature")
            logger.debug(f"{value}\tTranslational contribution")
            logger.debug(
                symmetry_correction_to_chemical_potential(temperature, symmetry_number)
                + linear_molecule_rotational_contribution(temperature, ixx),
                "\tRotational contribution",
            )
        return (
            value
            + symmetry_correction_to_chemical_potential(temperature, symmetry_number)
            + linear_molecule_rotational_contribution(temperature, ixx)
        )
    else:
        # polyatomic
        degrees_of_freedom = 3
        if debug:
            logger.debug(f"{temperature}\tTemperature")
            logger.debug(f"{value}\tTranslational contribution")
            logger.debug(
                symmetry_correction_to_chemical_potential(temperature, symmetry_number)
                + rotational_contribution_to_chemical_potential(
                    temperature, ixx, iyy, izz
                ),
                "\tRotational contribution",
            )
        return (
            value
            + symmetry_correction_to_chemical_potential(temperature, symmetry_number)
            + rotational_contribution_to_chemical_potential(temperature, ixx, iyy, izz)
        )


@deprecated_alias("molecular_chemical_potential_corrections")
def corrections_due_to_molecular_shape(
    temperature, mass, symmetry, moment_of_inertia, debug=False
):
    pass


# #ixx,iyy,izz: moment of inertia in (atomic mass) x A**2 (Use output of momentofinertia.py)
# def StericFix2(T,mass,symm,ixx=0,iyy=0,izz=0,debug=False):
#     value = Avalue(T,mass)
#     if ixx == 0:
#         #monatomic
#         dimen = 0
#         if debug:
#             logger.debug value, "\tAvalue"
#             logger.debug value, "\tAvalue"
#         return value
#     elif izz == 0:
#         #rodlike
#         dimen = 1
#         if debug:
#             logger.debug value, "\tAvalue"
#             logger.debug Cvalue2(T,ixx), "\tCvalue2"
#             logger.debug value, "\tAvalue"
#             logger.debug Cvalue2(T,ixx), "\tCvalue2"
#         return value + Cvalue2(T,ixx)
#     else:
#         #polyatomic
#         dimen = 3
#         if debug:
#             logger.debug value, "\tAvalue"
#             logger.debug Cvalue(T,ixx, iyy, izz), "\tCvalue"
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
#             logger.debug value, "\tAvalue"
#         return value
#     elif izz == 0:
#         #rodlike
#         dimen = 1
#         if debug:
#             logger.debug SymmetryFix(T, symm), "\tSymmetryFix"
#             logger.debug value, "\tAvalue"
#             logger.debug Cvalue2(T,ixx), "\tCvalue2"
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
#             logger.debug SymmetryFix(T, symm), "\tSymmetryFix"
#             logger.debug value, "\tAvalue"
#             logger.debug Cvalue(T,ixx, iyy, izz), "\tCvalue"
#         return value + SymmetryFix(T, symm) + Cvalue(T,ixx, iyy, izz)


# suggested_new_name: chemical_potential_of_ideal_gas
def ideal_gas_chemical_potential(temperature, pressure):
    """Calculate the chemical potential of an ideal gas.

    Args:
        temperature: Temperature in Kelvin
        pressure: Pressure in Pa
    """
    return -pc.NkB * temperature * (np.log(pc.kB * temperature / pressure))


@deprecated_alias("ideal_gas_chemical_potential")
def chempot(temperature, pressure):
    pass


def test():
    import antoine

    basicConfig(level=DEBUG)
    logger = getLogger()

    A = 5.40221
    B = 1838.675
    C = -31.737
    for T in range(273, 293):
        logger.debug(
            [
                T,
                ideal_gas_chemical_potential(T, 101326)
                + translational_contribution_to_chemical_potential(T, 16),
                ideal_gas_chemical_potential(
                    T, antoine.calculate_vapor_pressure(T, A, B, C)
                ),
            ]
        )


if __name__ == "__main__":
    test()
