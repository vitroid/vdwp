<script lang="ts">
  import { onMount } from 'svelte';
  import * as pc from '$lib/vdwp/physconst';
  import * as crystals from '$lib/vdwp/crystals';
  import * as chempot from '$lib/vdwp/chempot';
  import * as vdwp from '$lib/vdwp/vdwp';
  import { fvalue2 } from '$lib/ljd/ljd';
  import { inter } from '$lib/vdwp/inter';
  import { getCoexistenceTraces, baseLayout } from '$lib/vdwp/plot_utils';

  let plotDiv: HTMLDivElement;
  const temperature = 273.15; // Fixed temperature
  let pressureBar = $state(50);
  let isCalculating = $state(false);
  const stericterm = 0.0;

  let gasFreeEnergies: Record<string, Record<number, number>> = {};
  let lastCalcTemperature = 0;

  async function updatePlot() {
    const T = temperature;
    const P = pressureBar;
    if (!plotDiv) return;
    isCalculating = true;

    if (T !== lastCalcTemperature) {
      const beta = 1.0 / (pc.NkB * T);
      const newEnergies: Record<string, Record<number, number>> = {};
      const relevantGases = ['Methane', 'Ethane', 'C2H4', 'Xe', 'Br2', 'cC3H6'];
      for (const name of relevantGases) {
        const gas = inter[name];
        const f_c: Record<number, number> = {};
        for (const [cageStr, R] of Object.entries(crystals.radii)) {
          f_c[parseInt(cageStr)] =
            fvalue2(R, crystals.nmemb[parseInt(cageStr)], gas.sig, (gas.epsK * 8.314) / 1000, beta) +
            stericterm;
        }
        newEnergies[name] = f_c;
      }
      gasFreeEnergies = newEnergies;
      lastCalcTemperature = T;
    }

    const Plotly = (await import('plotly.js-dist-min')).default;
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
      const f_a = gasFreeEnergies[a];
      const f_b = gasFreeEnergies[b];
      
      for (const frac of ticks) {
        const mu_a = chempot.ideal_gas_chemical_potential(T, P * (1 - frac) * 101325) + chempot.phase_space_integration_correction(T, 0);
        const mu_b = chempot.ideal_gas_chemical_potential(T, P * frac * 101325) + chempot.phase_space_integration_correction(T, 0);
        
        const Deltamu_c = vdwp.calculate_chemical_potential_by_occupation(
          T, [f_a, f_b], [mu_a, mu_b], ['CS1', 'CS2', 'HS1']
        );
        mx.push(Deltamu_c['CS1'] - Deltamu_c['HS1']);
        my.push(Deltamu_c['CS2'] - Deltamu_c['HS1']);
      }
      mixtureTraces.push({ x: mx, y: my, mode: 'lines', name: `${a}-${b}`, line: { color, width: 2 } });
    }

    Plotly.react(plotDiv, mixtureTraces, {
      ...baseLayout,
      title: `Binary Mixtures (${T.toFixed(2)} K, Total ${P} bar)`
    });
    isCalculating = false;
  }

  $effect(() => {
    temperature;
    pressureBar;
    updatePlot();
  });

  onMount(updatePlot);
</script>

<div class="panel">
  <div class="controls">
    <h3>Binary Mixture Plot {#if isCalculating}<span class="loading">...</span>{/if}</h3>
    <div class="control-group">
      <label>Total Pressure: {pressureBar} bar</label>
      <input type="range" min="1" max="500" step="1" bind:value={pressureBar} />
      <input type="number" bind:value={pressureBar} step="1" />
    </div>
  </div>
  <div bind:this={plotDiv}></div>
</div>

<style>
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
    gap: 5px;
  }
  .control-group {
    display: flex;
    align-items: center;
    gap: 10px;
    justify-content: space-between;
  }
  .control-group label {
    font-size: 0.9em;
    min-width: 120px;
    text-align: left;
  }
  .loading {
    color: #ff3e00;
  }
  h3 { margin: 0 0 10px 0; font-size: 1.1em; }
</style>
