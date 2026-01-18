import * as crystals from './crystals';

export function getCoexistenceTraces() {
  const traces: any[] = [];
  const names = ["CS1", "CS2", "HS1"];
  for (let i = 0; i < names.length; i++) {
    for (let j = i + 1; j < names.length; j++) {
      const s1 = names[i];
      const s2 = names[j];
      const A = crystals.ratios[s1][0] - crystals.ratios[s2][0];
      const B = crystals.ratios[s1][1] - crystals.ratios[s2][1];
      const C = crystals.mu_e[s1] - crystals.mu_e[s2];

      const xRange = [-0.4, 0.4];
      let lx: number[] = [];
      let ly: number[] = [];

      if (Math.abs(B) > 1e-9) {
        lx = xRange;
        ly = lx.map((x) => (-C - A * x) / B);
      } else if (Math.abs(A) > 1e-9) {
        const x = -C / A;
        lx = [x, x];
        ly = [-0.1, 0.7];
      }

      if (lx.length > 0) {
        traces.push({
          x: lx,
          y: ly,
          mode: 'lines',
          line: { color: '#ccc', width: 1, dash: 'dash' },
          showlegend: false,
          hoverinfo: 'none',
        });
      }
    }
  }
  return traces;
}

export const baseLayout: any = {
  xaxis: { title: 'Δμc(CS1) - Δμc(HS1) / kJ/mol', range: [-0.4, 0.3], zeroline: false },
  yaxis: {
    title: 'Δμc(CS2) - Δμc(HS1) / kJ/mol',
    range: [-0.05, 0.65],
    scaleanchor: 'x',
    scaleratio: 1,
    zeroline: false,
  },
  width: 500,
  height: 500,
  margin: { t: 50, b: 50, l: 60, r: 30 },
  showlegend: true,
  legend: { x: 0, y: 1 },
};
