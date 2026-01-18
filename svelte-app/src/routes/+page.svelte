<script lang="ts">
  import { dndzone } from 'svelte-dnd-action';
  import { flip } from 'svelte/animate';
  import SinglePanel from '$lib/components/SinglePanel.svelte';
  import MixturePanel from '$lib/components/MixturePanel.svelte';

  // Svelte 5 state for the list of panels
  let items = $state([
    { id: 'single', component: SinglePanel },
    { id: 'mixture', component: MixturePanel }
  ]);

  function handleDndConsider(e: any) {
    items = e.detail.items;
  }

  function handleDndFinalize(e: any) {
    items = e.detail.items;
  }
</script>

<main>
  <h1>Hydrate Phase Diagrams (Independent Panels)</h1>
  <p class="hint">ドラッグしてパネルの順序を入れ替えられます</p>
  
  <section 
    class="container" 
    use:dndzone={{items, flipDurationMs: 300, dropTargetStyle: {outline: 'none'}}} 
    onconsider={handleDndConsider} 
    onfinalize={handleDndFinalize}
  >
    {#each items as item (item.id)}
      <div class="panel-wrapper" animate:flip={{duration: 300}}>
        <item.component />
      </div>
    {/each}
  </section>
</main>

<style>
  main {
    font-family: sans-serif;
    text-align: center;
    padding: 20px;
  }
  .hint {
    color: #666;
    font-size: 0.9em;
    margin-bottom: 20px;
  }
  .container {
    display: flex;
    flex-wrap: wrap;
    justify-content: center;
    gap: 20px;
    min-height: 500px;
    padding: 10px;
  }
  .panel-wrapper {
    outline: none;
  }
  /* パネル全体を掴めるようにする */
  :global(.panel) {
    cursor: grab;
  }
  :global(.panel:active) {
    cursor: grabbing;
  }
</style>
