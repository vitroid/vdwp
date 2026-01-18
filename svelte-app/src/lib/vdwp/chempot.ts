import * as pc from "./physconst";

/**
 * Calculate the translational contribution to chemical potential for an ideal gas.
 * @param temperature Temperature in Kelvin
 * @param mass Molecular mass in g/mol
 */
export function translational_contribution_to_chemical_potential(
  temperature: number,
  mass: number
): number {
  return (
    -pc.NkB *
    temperature *
    (3.0 / 2.0) *
    Math.log(
      ((mass * 0.001) / pc.NA) *
        2.0 *
        Math.PI *
        pc.kB *
        temperature /
        Math.pow(pc.h, 2)
    )
  );
}

/**
 * Calculate the correction for phase space integration.
 * @param temperature Temperature in Kelvin
 * @param degrees_of_freedom Number of rotational degrees of freedom (0 for monatomic, 1 for linear, 3 for non-linear)
 */
export function phase_space_integration_correction(
  temperature: number,
  degrees_of_freedom: number
): number {
  if (degrees_of_freedom === 3) {
    return -pc.NkB * temperature * Math.log(8 * Math.pow(Math.PI, 2));
  } else if (degrees_of_freedom === 1) {
    return -pc.NkB * temperature * Math.log(4 * Math.PI);
  } else {
    return 0.0;
  }
}

/**
 * Calculate the chemical potential of an ideal gas.
 * @param temperature Temperature in Kelvin
 * @param pressure Pressure in Pa
 */
export function ideal_gas_chemical_potential(
  temperature: number,
  pressure: number
): number {
  return -pc.NkB * temperature * Math.log((pc.kB * temperature) / pressure);
}
