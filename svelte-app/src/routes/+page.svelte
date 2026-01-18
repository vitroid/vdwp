<script lang="ts">
    // ... (インポート部分はそのまま) ...
    import { onMount } from 'svelte';
    import * as pc from '$lib/vdwp/physconst';
    import * as crystals from '$lib/vdwp/crystals';
    import * as chempot from '$lib/vdwp/chempot';
    import * as vdwp from '$lib/vdwp/vdwp';
    import { fvalue2 } from '$lib/ljd/ljd';
    import { inter } from '$lib/vdwp/inter';
    import { MultipleClathrate } from '$lib/vdwp/calculations';
  
    let plotSingle: HTMLDivElement;
    let plotMixture: HTMLDivElement;
  
    // Svelte 5 state
    let temperature = $state(273.15);
    let pressureBar = $state(50);
    let isCalculating = $state(false);
    let methanePos = $state({ x: 0, y: 0 }); // 更新確認用
    const stericterm = 0.0;
  
    // ... (getCoexistenceTraces, baseLayout はそのまま) ...
    function getCoexistenceTraces() {
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
  
    const baseLayout: any = {
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
  
    async function updatePlots() {
      if (!plotSingle || !plotMixture) return;
      isCalculating = true;
      
      // 少し遅延を入れてUIに「Calculating」を表示させる（計算が速すぎる場合用）
      await new Promise(resolve => setTimeout(resolve, 0));
  
      const Plotly = (await import('plotly.js-dist-min')).default;
      const beta = 1.0 / (pc.NkB * temperature);
  
      // --- Single Plot ---
      const singleTraces = [...getCoexistenceTraces()];
      for (const [name, gas] of Object.entries(inter)) {
        const f_c: Record<number, number> = {};
        for (const [cageStr, R] of Object.entries(crystals.radii)) {
          f_c[parseInt(cageStr)] =
            fvalue2(R, crystals.nmemb[parseInt(cageStr)], gas.sig, (gas.epsK * 8.314) / 1000, beta) +
            stericterm;
        }
  
        const mu_g =
          chempot.ideal_gas_chemical_potential(temperature, pressureBar * 101326) +
          chempot.phase_space_integration_correction(temperature, 0);
        const Deltamu_c = vdwp.calculate_chemical_potential_by_occupation(
          temperature, [f_c], [mu_g], ['CS1', 'CS2', 'HS1']
        );
  
        const x = Deltamu_c['CS1'] - Deltamu_c['HS1'];
        const y = Deltamu_c['CS2'] - Deltamu_c['HS1'];
        
        if (name === 'Methane') {
          methanePos = { x, y };
        }
  
        singleTraces.push({
          x: [x], y: [y],
          mode: 'markers+text',
          name: name,
          text: [name],
          textposition: 'top center',
          marker: { size: 8 },
        });
      }
      Plotly.react(plotSingle, singleTraces, { ...baseLayout, title: `Guests at ${temperature} K, ${pressureBar} bar` });
  
      // --- Mixture Plot ---
      const mixtureTraces = [...getCoexistenceTraces()];
      const mixturePairs: [string, string, string][] = [
        ['Methane', 'Ethane', 'red'],
        ['Methane', 'C2H4', 'blue'],
        ['Xe', 'Br2', 'green'],
        ['Methane', 'cC3H6', 'orange'],
      ];
  
      const ticks = Array.from({ length: 50 }, (_, i) => i / 49);
      for (const [a, b, color] of mixturePairs) {
        const mx: number[] = [];
        const my: number[] = [];
        for (const frac of ticks) {
          const res = MultipleClathrate(
            [inter[a], inter[b]],
            [pressureBar * (1 - frac) * 101326, pressureBar * frac * 101326],
            temperature,
            ['CS1', 'CS2', 'HS1']
          );
          mx.push(res.X);
          my.push(res.Y);
        }
        mixtureTraces.push({ x: mx, y: my, mode: 'lines', name: `${a}-${b}`, line: { color, width: 2 } });
      }
      Plotly.react(plotMixture, mixtureTraces, { ...baseLayout, title: `Binary Mixtures (Total ${pressureBar} bar)` });
      
      isCalculating = false;
      console.log(`Updated at T=${temperature}, P=${pressureBar}. Methane X: ${methanePos.x.toFixed(4)}`);
    }
  
    $effect(() => {
      // temperature または pressureBar が変わるたびに実行
      updatePlots();
    });
  
    onMount(updatePlots);
</script>
  
<main>
<h1>Hydrate Phase Diagrams {#if isCalculating}<span class="loading">...Calculating...</span>{/if}</h1>

<div class="status-bar">
    Methane Position: X = {methanePos.x.toFixed(5)}, Y = {methanePos.y.toFixed(5)}
</div>

<div class="container">
    <div class="panel">
    <div class="controls">
        <div class="control-group">
        <label for="temp">Temperature: {temperature.toFixed(2)} K</label>
        <input type="range" id="temp" min="200" max="350" step="0.1" bind:value={temperature} />
        <input type="number" bind:value={temperature} step="0.1" />
        </div>
        <div class="control-group">
        <label for="press">Pressure: {pressureBar} bar</label>
        <input type="range" id="press" min="1" max="500" step="1" bind:value={pressureBar} />
        <input type="number" bind:value={pressureBar} step="1" />
        </div>
    </div>
    <div bind:this={plotSingle}></div>
    </div>
    <div class="panel">
    <div bind:this={plotMixture}></div>
    </div>
</div>
</main>

<style>
/* ... (既存のスタイル) ... */
main {
    font-family: sans-serif;
    text-align: center;
}
.loading {
    color: #ff3e00;
    font-size: 0.5em;
    vertical-align: middle;
}
.status-bar {
    background: #333;
    color: #0f0;
    font-family: monospace;
    padding: 5px;
    margin-bottom: 10px;
    font-size: 0.9em;
}
.container {
    display: flex;
    flex-wrap: wrap;
    justify-content: center;
    gap: 20px;
}
.panel {
    background: #f9f9f9;
    border: 1px solid #ddd;
    border-radius: 8px;
    padding: 10px;
    display: flex;
    flex-direction: column;
    align-items: center;
}
.controls {
    width: 100%;
    padding: 10px;
    background: #eee;
    border-radius: 4px;
    margin-bottom: 10px;
    display: flex;
    flex-direction: column;
    gap: 10px;
}
.control-group {
    display: flex;
    align-items: center;
    gap: 10px;
    justify-content: space-between;
}
.control-group label {
    font-weight: bold;
    min-width: 150px;
    text-align: left;
}
.control-group input[type='range'] {
    flex-grow: 1;
}
.control-group input[type='number'] {
    width: 80px;
}
</style>