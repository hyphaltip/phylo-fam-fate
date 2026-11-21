/* FastTree runs in a dedicated worker so long inferences never freeze the UI
   and can be cancelled (the main thread simply terminates this worker). */
import init, { buildTree, fastTreeVersion } from '../assets/fasttree/phylo_fasttree_wasm.js';

let wasm = null;
const ready = init().then(exports => { wasm = exports; });

/* FastTree's exit() paths throw a JS error straight out of WASM. That skips
   the C epilogues, so the shadow-stack pointer is never restored and drifts
   down with every failure until memory is corrupted (~380 errors). Restore it
   ourselves after any throw. */
function guardedBuildTree(fasta, args) {
  const sp = wasm.__wbindgen_add_to_stack_pointer(0);
  try {
    return buildTree(fasta, args);
  } catch (e) {
    wasm.__wbindgen_add_to_stack_pointer(sp - wasm.__wbindgen_add_to_stack_pointer(0));
    throw e;
  }
}

self.onmessage = async ({ data }) => {
  try {
    await ready;
    if (data.type === 'version') {
      self.postMessage({ id: data.id, ok: true, version: fastTreeVersion() });
      return;
    }
    const t0 = performance.now();
    const newick = guardedBuildTree(data.fasta, data.args);
    self.postMessage({ id: data.id, ok: true, newick, ms: performance.now() - t0 });
  } catch (e) {
    self.postMessage({ id: data.id, ok: false, error: e?.message || String(e) });
  }
};
