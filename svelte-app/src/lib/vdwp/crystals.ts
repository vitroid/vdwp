export const names = ["CS1", "CS2", "TS1", "HS1"];
export const aliases: Record<string, string> = {
  CS1: "sI",
  CS2: "sII",
  TS1: "sIII",
  HS1: "sIV",
};

export interface CageComponents {
  nw: number;
  cages: Record<number, number>;
}

export const cage_components: Record<string, [number, Record<number, number>]> = {
  CS2: [136.0, { 16: 8.0, 12: 16.0 }],
  CS1: [46.0, { 12: 2.0, 14: 6.0 }],
  TS1: [172, { 12: 10, 14: 16, 15: 4 }],
  HS1: [80, { 12: 6, 14: 4, 15: 4 }],
};

export const ratios: Record<string, [number, number, number]> = {
  CS1: [1.0, 0.0, 0.0],
  CS2: [0.0, 1.0, 0.0],
  HS1: [0.0, 0.0, 1.0],
  TS1: [23 / 43, 0.0, 20 / 43],
};

const tip4pice_data: Record<string, string> = {
  CS1: "   -59.841454832900  1242  22.723050432819",
  CS2: "   -59.922551474335  1088  23.039907104034",
  TS1: "   -59.775409212738  1376  22.788576341391",
  HS1: "   -59.608521883866  1440  22.826573500133",
};

export const mu_e: Record<string, number> = {};
for (const [structure, s] of Object.entries(tip4pice_data)) {
  const cols = s.trim().split(/\s+/);
  mu_e[structure] = parseFloat(cols[0]);
}

export const radii: Record<number, number> = {
  12: 3.988,
  14: 4.331,
  15: 4.527,
  16: 4.587,
};

export const nmemb: Record<number, number> = {
  12: 20,
  14: 24,
  15: 26,
  16: 28,
};
