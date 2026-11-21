/* FastTree runs in a dedicated worker so long inferences never freeze the UI
   and can be cancelled (the main thread simply terminates this worker).
   Two builds are available: single precision (FastTree's default) and double
   precision (-DUSE_DOUBLE, like bioconda's binary). Each is loaded on first
   use, so the double build is only downloaded if someone selects it. */
const BUILDS = {
  single: '../assets/fasttree/phylo_fasttree_wasm.js',
  double: '../assets/fasttree/phylo_fasttree_wasm_double.js',
};
const loaded = {};   // precision -> Promise<{mod, wasm}>

function load(precision = 'single') {
  if (!BUILDS[precision]) throw new Error(`unknown precision: ${precision}`);
  loaded[precision] ||= import(BUILDS[precision]).then(async mod => ({ mod, wasm: await mod.default() }));
  return loaded[precision];
}

/* FastTree's exit() paths throw a JS error straight out of WASM. That skips
   the C epilogues, so the shadow-stack pointer is never restored and drifts
   down with every failure until memory is corrupted (~380 errors). Restore it
   ourselves after any throw. */
function guardedBuildTree({ mod, wasm }, fasta, args) {
  const sp = wasm.__wbindgen_add_to_stack_pointer(0);
  try {
    return mod.buildTree(fasta, args);
  } catch (e) {
    wasm.__wbindgen_add_to_stack_pointer(sp - wasm.__wbindgen_add_to_stack_pointer(0));
    throw e;
  }
}

self.onmessage = async ({ data }) => {
  try {
    const build = await load(data.precision);
    if (data.type === 'version') {
      self.postMessage({ id: data.id, ok: true, version: build.mod.fastTreeVersion(),
        precision: build.mod.fastTreePrecision?.() || data.precision || 'single' });
      return;
    }
    const t0 = performance.now();
    const newick = guardedBuildTree(build, data.fasta, data.args);
    self.postMessage({ id: data.id, ok: true, newick, ms: performance.now() - t0,
      precision: build.mod.fastTreePrecision?.() || data.precision || 'single' });
  } catch (e) {
    self.postMessage({ id: data.id, ok: false, error: e?.message || String(e) });
  }
};
