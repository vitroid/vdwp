import * as pc from "./physconst";
import * as crystals from "./crystals";
import * as chempot from "./chempot";
import * as vdwp from "./vdwp";
import { fvalue2 } from "../ljd/ljd";
import { inter, type Gas } from "./inter";

export function MultipleClathrate(
  gases: Gas[],
  pressures: number[],
  temperatures: number,
  structures: string[]
): { X: number; Y: number } {
  const beta = 1.0 / (pc.NkB * temperatures);
  const f_list: Record<number, number>[] = [];
  const mu_list: number[] = [];

  for (let i = 0; i < gases.length; i++) {
    const gas = gases[i];
    const pressure = pressures[i];
    if (pressure > 0.0) {
      const mu =
        chempot.ideal_gas_chemical_potential(temperatures, pressure) +
        chempot.phase_space_integration_correction(temperatures, 0);
      mu_list.push(mu);

      const ff: Record<number, number> = {};
      for (const [cageStr, R] of Object.entries(crystals.radii)) {
        const cage = parseInt(cageStr);
        // inter[name] should be used if gas is from gastable
        // Here we assume the passed Gas objects are already from 'inter'
        ff[cage] = fvalue2(
          R,
          crystals.nmemb[cage],
          gas.sig,
          (gas.epsK * 8.314) / 1000,
          beta
        );
      }
      f_list.push(ff);
    }
  }

  const Deltamu = vdwp.calculate_chemical_potential_by_occupation(
    temperatures,
    f_list,
    mu_list,
    structures
  );

  const X = Deltamu["CS1"] - Deltamu["HS1"];
  const Y = Deltamu["CS2"] - Deltamu["HS1"];
  return { X, Y };
}
