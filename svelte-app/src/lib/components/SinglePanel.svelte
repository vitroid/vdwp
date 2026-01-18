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
      for (const [name, gas] of Object.entries(inter)) {
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
    const singleTraces = [...getCoexistenceTraces()];

    for (const [name, gas] of Object.entries(inter)) {
      const f_c = gasFreeEnergies[name];
      const mu_g =
        chempot.ideal_gas_chemical_potential(T, P * 101325) +
        chempot.phase_space_integration_correction(T, 0);
      
      const Deltamu_c = vdwp.calculate_chemical_potential_by_occupation(
        T, [f_c], [mu_g], ['CS1', 'CS2', 'HS1']
      );

      singleTraces.push({
        x: [Deltamu_c['CS1'] - Deltamu_c['HS1']],
        y: [Deltamu_c['CS2'] - Deltamu_c['HS1']],
        mode: 'markers+text',
        name: name,
        text: [name],
        textposition: 'top center',
        marker: { size: 8 },
      });
    }

    Plotly.react(plotDiv, singleTraces, {
      ...baseLayout,
      title: `Single Guests (${T.toFixed(2)} K, ${P} bar)`,
      showlegend: false
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
    <h3>Single Component Plot {#if isCalculating}<span class="loading">...</span>{/if}</h3>
    <div class="control-group">
      <label>Pressure: {pressureBar} bar</label>
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
