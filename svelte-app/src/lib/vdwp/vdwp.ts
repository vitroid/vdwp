import * as pc from "./physconst";
import * as crystals from "./crystals";

/**
 * Calculate chemical potential of water in the hydrate phase.
 */
export function calculate_chemical_potential_by_occupation(
  temperatures: number,
  f_c_list: Record<number, number>[],
  mu_guest_list: number[],
  structures: string[]
): Record<string, number> {
  const Deltamu: Record<string, number> = {};

  for (const structure of structures) {
    const [Nw, cages] = crystals.cage_components[structure];
    let sum = 0;
    const beta = 1.0 / (pc.NkB * temperatures);

    for (const [cageStr, count] of Object.entries(cages)) {
      const cage = parseInt(cageStr);
      const alpha = count / Nw;
      const coeffs = -pc.NkB * temperatures * alpha;
      let sum0 = 1.0;

      for (let i = 0; i < f_c_list.length; i++) {
        const f = f_c_list[i];
        const mu = mu_guest_list[i];
        if (cage in f) {
          sum0 += Math.exp(beta * (mu - f[cage]));
        }
      }
      sum += coeffs * Math.log(sum0);
    }
    Deltamu[structure] = sum;
  }

  return Deltamu;
}
