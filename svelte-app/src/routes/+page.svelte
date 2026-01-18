<script lang="ts">
  import SinglePanel from '$lib/components/SinglePanel.svelte';
  import MixturePanel from '$lib/components/MixturePanel.svelte';
</script>

<main>
  <h1>Hydrate Phase Diagrams (Independent Panels)</h1>
  
  <div class="container">
    <SinglePanel />
    <MixturePanel />
  </div>
</main>

<style>
  main {
    font-family: sans-serif;
    text-align: center;
    padding: 20px;
  }
  .container {
    display: flex;
    flex-wrap: wrap;
    justify-content: center;
    gap: 20px;
  }
</style>
