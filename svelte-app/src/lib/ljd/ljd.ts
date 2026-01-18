/**
 * Lennard-Jones--Devonshire cellpotential function.
 * @param r position of the guest from the center (Angstrom)
 * @param sigma sigma (Angstrom)
 * @param epsilon epsilon (kJ/mol)
 * @param z coordination number
 * @param R free cavity radius (Angstrom)
 * @param a shift (default 0.0)
 */
export function cellpotential(
  r: number,
  sigma: number,
  epsilon: number,
  z: number,
  R: number,
  a: number = 0.0
): number {
  const delta = (N: number): number => {
    return (Math.pow(1 - r / R - a / R, -N) - Math.pow(1 + r / R - a / R, -N)) / N;
  };

  const ans =
    2 *
    z *
    epsilon *
    (Math.pow(sigma, 12) / (Math.pow(R, 11) * r) * (delta(10) + (a * delta(11)) / R) -
      Math.pow(sigma, 6) / (Math.pow(R, 5) * r) * (delta(4) + (a * delta(5)) / R));
  return ans;
}

/**
 * Calculate the free energy of cage occupation using LJD approximation.
 * @param R radius
 * @param z number of water molecules
 * @param sigma sigma
 * @param epsilon epsilon
 * @param beta beta (1/NkB*T)
 */
export function fvalue2(
  R: number,
  z: number,
  sigma: number,
  epsilon: number,
  beta: number
): number {
  const n = 100;
  const dr = R / (n - 1);
  const r = Array.from({ length: n }, (_, i) => i * dr);
  const integrand = new Float64Array(n);

  for (let i = 0; i < n; i++) {
    if (i === 0) {
      // Limit as r -> 0 for cellpotential
      // The cellpotential has a 1/r factor, but delta(N) also has r.
      // For r=0, cellpotential(0) = 2*z*epsilon * [ sigma^12/R^12 - sigma^6/R^6 ]
      // Wait, let's check the limit.
      // (1-x)^-N - (1+x)^-N = [1 + Nx + N(N+1)/2 x^2 + ...] - [1 - Nx + N(N+1)/2 x^2 - ...] = 2Nx + O(x^3)
      // So delta(N) / r = 2/R + O(r^2)
      // cellpotential(0) = 2*z*epsilon * [ sigma^12/R^12 - sigma^6/R^6 ]
      const cp0 = 2 * z * epsilon * (Math.pow(sigma / R, 12) - Math.pow(sigma / R, 6));
      integrand[i] = Math.exp(-beta * cp0) * Math.pow(r[i], 2) * 1e-30;
    } else if (i < n - 1) {
      const cp = cellpotential(r[i], sigma, epsilon, z, R);
      integrand[i] = Math.exp(-beta * cp) * Math.pow(r[i], 2) * 1e-30;
    } else {
      // r = R, potential is infinite
      integrand[i] = 0;
    }
  }

  // Trapezoidal rule
  let integral = 0;
  for (let i = 0; i < n - 1; i++) {
    integral += (integrand[i] + integrand[i + 1]) * dr / 2;
  }

  return -Math.log(4 * Math.PI * integral) / beta;
}
