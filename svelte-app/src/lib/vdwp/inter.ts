export interface Gas {
  name: string;
  sig: number;
  epsK: number;
  TeX: string;
}

export const gastable_raw = `
Methane  3.758 148.6 $\\mathrm{CH}_4$
Ethane   4.520 208.8
C2H4     4.232 205 $\\mathrm{C_2H_4}$
Ne       2.749  35.6
Ar       3.405 119.8
Kr       3.60  171.0
Xe       4.047 231.0
Br2      4.93   540 Q                        # manually set. Becomes sIII at 10 bar.
CO2      4.486 189.0 $\\mathrm{CO}_2$
CS2      4.438 488 $\\mathrm{CS_2}$
cC3H6    4.582733199595731 301.51400454201365 $c\\mathrm{C_3H_6}$ # from CP
`;

export const tip4pice: Gas = {
  name: "tip4pice",
  sig: 3.1668,
  epsK: 106.1,
  TeX: "$\\mathrm{tip4pice}$",
};

export const gases: Record<string, Gas> = {};
export const inter: Record<string, Gas> = {};

const lines = gastable_raw.trim().split("\n");
for (let line of lines) {
  const hashIdx = line.indexOf("#");
  if (hashIdx >= 0) {
    line = line.substring(0, hashIdx);
  }
  const cols = line.trim().split(/\s+/);
  if (cols.length < 3) continue;

  const name = cols[0];
  const sig = parseFloat(cols[1]);
  const epsK = parseFloat(cols[2]);
  let tex = cols.length >= 4 ? cols.slice(3).join(" ") : `$\\mathrm{${name}}$`;

  gases[name] = { name, sig, epsK, TeX: tex };
  inter[name] = {
    name,
    sig: (sig + tip4pice.sig) / 2,
    epsK: Math.sqrt(epsK * tip4pice.epsK),
    TeX: tex,
  };
}
