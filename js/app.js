/* Main application: orchestrates the search -> align -> trim -> tree pipeline. */

import { retrieve, parseFasta, toFasta } from './api.js';
import { runAligner, runBlast, runPhmmer, parseBlastHits, parsePhmmerHits } from './tools.js';
import { trimAlignment } from './trim.js';
import { parseNewick, drawTree, treeStats } from './tree.js';
import initFastTree, { buildTree, fastTreeVersion } from '../assets/fasttree/phylo_fasttree_wasm.js';

const $ = (id) => document.getElementById(id);

// Pipeline state
const state = {
  rawFasta: '',       // original input sequences
  searchFasta: '',    // sequences from search hits
  alignedFasta: '',   // aligned sequences
  trimmedFasta: '',   // trimmed alignment
  newick: '',
  records: [],
  metadata: [],
};

const VERSION = '0.1.0';

// ---------------------------------------------------------------------------
// UI helpers
// ---------------------------------------------------------------------------
function status(msg, kind = '') {
  const bar = $('global-status');
  bar.textContent = msg;
  bar.className = `status-bar ${kind}`;
}
function stepStatus(id, text, kind = '') {
  const el = $(`status-${id}`);
  el.textContent = text;
  el.className = `step-status ${kind}`;
}
function recordMeta(step, tool, version, params, inputDesc, outputDesc, ms) {
  state.metadata.push({
    step, tool, version, params,
    input: inputDesc, output: outputDesc,
    time: ms !== null ? `${(ms / 1000).toFixed(2)}s` : '—',
    timestamp: new Date().toISOString(),
  });
  renderMetadata();
}
function renderMetadata() {
  const tbody = document.querySelector('#metadata-table tbody');
  tbody.innerHTML = state.metadata.map(m => `
    <tr>
      <td>${m.step}</td><td>${m.tool}</td><td>${m.version}</td>
      <td title="${escapeAttr(m.params)}">${escapeHtml(trunc(m.params, 60))}</td>
      <td>${escapeHtml(m.input)}</td><td>${escapeHtml(m.output)}</td>
      <td>${m.time}</td>
    </tr>`).join('');
}
function escapeHtml(s) { return String(s).replace(/[&<>]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;' }[c])); }
function escapeAttr(s) { return escapeHtml(s).replace(/"/g, '&quot;'); }
function trunc(s, n) { return String(s).length > n ? String(s).slice(0, n) + '…' : String(s); }

function seqType(fasta) {
  const residues = fasta.split(/\r?\n/)
    .filter(l => l.trim() && !l.trim().startsWith('>'))
    .join('').toUpperCase();
  if (!residues) return 'nt';
  // Residues that are protein-specific in practice (never IUPAC nucleotide codes)
  if (/[EFILPQZXJ]/.test(residues)) return 'aa';
  // Otherwise, nucleotide alignments are overwhelmingly ACGT(U)
  const acgt = (residues.match(/[ACGTU]/g) || []).length;
  return acgt / residues.length < 0.8 ? 'aa' : 'nt';
}

// ---------------------------------------------------------------------------
// Initialization
// ---------------------------------------------------------------------------
try {
  await initFastTree();
  $('tree-run').disabled = false;
  status(`FastTree ${fastTreeVersion()} and Biowasm tools ready`, 'success');
  stepStatus('tree', 'ready');
} catch (e) {
  status(`Failed to load FastTree WASM: ${e}`, 'error');
}

// ---------------------------------------------------------------------------
// STEP 1: Input
// ---------------------------------------------------------------------------
document.querySelectorAll('#input-tabs button').forEach(btn => {
  btn.addEventListener('click', () => {
    document.querySelectorAll('#input-tabs button').forEach(b => b.classList.remove('active'));
    btn.classList.add('active');
    document.querySelectorAll('.tab-panel').forEach(p => p.classList.add('hidden'));
    $(`tab-${btn.dataset.tab}`).classList.remove('hidden');
  });
});

const dropzone = $('dropzone');
const fileInput = $('file-input');
dropzone.addEventListener('click', () => fileInput.click());
dropzone.addEventListener('dragover', e => { e.preventDefault(); dropzone.classList.add('dragover'); });
dropzone.addEventListener('dragleave', () => dropzone.classList.remove('dragover'));
dropzone.addEventListener('drop', e => {
  e.preventDefault();
  dropzone.classList.remove('dragover');
  if (e.dataTransfer.files.length) loadInputFile(e.dataTransfer.files[0]);
});
fileInput.addEventListener('change', () => {
  if (fileInput.files.length) loadInputFile(fileInput.files[0]);
});

function loadInputFile(file) {
  const reader = new FileReader();
  reader.onload = () => {
    setRawFasta(reader.result, `${file.name} (${(file.size / 1024).toFixed(1)} KB)`);
  };
  reader.onerror = () => status('Could not read file', 'error');
  reader.readAsText(file);
}

$('paste-fasta').addEventListener('input', () => setRawFasta($('paste-fasta').value, 'pasted text'));

function setRawFasta(fasta, source) {
  state.rawFasta = fasta.trim();
  state.records = parseFasta(state.rawFasta);
  state.alignedFasta = '';
  state.trimmedFasta = '';
  state.newick = '';
  $('newick-output').textContent = '';
  $('tree-results').classList.add('hidden');

  if (!state.records.length) {
    $('sequence-summary').textContent = '';
    stepStatus('input', 'waiting');
    return;
  }
  const type = seqType(state.rawFasta);
  const lengths = state.records.map(r => r.seq.length);
  const min = Math.min(...lengths), max = Math.max(...lengths);
  $('sequence-summary').innerHTML = `
    <strong>${state.records.length} sequences</strong> from ${escapeHtml(source)} ·
    ${type === 'nt' ? 'nucleotide' : 'protein'} ·
    length ${min === max ? min : `${min}–${max}`} aa/nt ·
    ${state.rawFasta.length.toLocaleString()} characters`;
  $('file-info').textContent = source;
  stepStatus('input', `${state.records.length} seq`, 'done');
  status(`Loaded ${state.records.length} sequences — you can align or build a tree`, 'success');
  stepStatus('align', 'waiting');
  stepStatus('trim', 'waiting');
  stepStatus('tree', 'waiting');
  $('search-enable').checked = false;
  $('search-options').classList.add('hidden');
  $('align-input').value = 'original';
}

// API retrieval
$('api-fetch').addEventListener('click', async () => {
  const source = $('api-source').value;
  const query = $('api-query').value.trim();
  const limit = parseInt($('api-limit').value, 10) || 25;
  if (!query) { status('Enter a query first', 'error'); return; }

  stepStatus('input', 'retrieving…', 'running');
  status(`Retrieving from ${source}…`);
  $('api-fetch').disabled = true;
  try {
    const records = await retrieve(source, query, limit);
    if (!records.length) {
      status('No sequences found', 'error');
      stepStatus('input', 'no hits');
      return;
    }
    const fasta = toFasta(records.slice(0, limit));
    $('api-preview').textContent = fasta.slice(0, 2000) + (fasta.length > 2000 ? '\n…' : '');
    $('api-preview').classList.remove('hidden');
    setRawFasta(fasta, `${source} API query “${trunc(query, 40)}”`);
    status(`Retrieved ${records.length} sequences from ${source}`, 'success');
  } catch (e) {
    status(`API error: ${e.message}`, 'error');
    stepStatus('input', 'error');
  } finally {
    $('api-fetch').disabled = false;
  }
});

// ---------------------------------------------------------------------------
// STEP 2: Search (optional)
// ---------------------------------------------------------------------------
$('search-enable').addEventListener('change', () => {
  $('search-options').classList.toggle('hidden', !$('search-enable').checked);
  if ($('search-enable').checked && !$('search-db').files.length) {
    status('Select a database FASTA file for the local search', 'info');
  }
});

$('search-run').addEventListener('click', async () => {
  const tool = $('search-tool').value;
  const dbFile = $('search-db').files[0];
  if (!state.rawFasta) { status('Load input sequences first', 'error'); return; }
  if (!dbFile) { status('Select a database FASTA file', 'error'); return; }

  const dbFasta = await dbFile.text();
  const maxHits = parseInt($('search-maxhits').value, 10) || 50;
  stepStatus('search', 'running…', 'running');
  status(`Running ${tool} locally…`);
  $('search-run').disabled = true;
  const t0 = performance.now();
  try {
    let hits;
    if (tool === 'phmmer') {
      const out = await runPhmmer(state.rawFasta, dbFasta, maxHits);
      hits = parsePhmmerHits(out);
      $('search-preview').textContent = out.slice(0, 3000);
    } else {
      const out = await runBlast(tool, state.rawFasta, dbFasta, maxHits);
      hits = parseBlastHits(out);
      $('search-preview').textContent = out.slice(0, 3000);
    }
    $('search-preview').classList.remove('hidden');

    // Extract hit sequences from the database FASTA
    const dbRecords = parseFasta(dbFasta);
    const hitRecords = dbRecords.filter(r => hits.includes(r.name) || hits.some(h => h.startsWith(r.name)));
    if (!hitRecords.length) {
      status(`Search ran but no database sequences could be matched`, 'error');
      stepStatus('search', 'no hits');
      return;
    }
    state.searchFasta = toFasta(hitRecords);
    $('align-input').value = 'search';
    $('align-input').querySelector('option[value="search"]').textContent = `Search hits (${hitRecords.length})`;
    stepStatus('search', `${hitRecords.length} hits`, 'done');
    recordMeta('Search', tool, 'biowasm', `maxHits=${maxHits}`,
      `${state.records.length} queries + ${dbRecords.length} db`, `${hitRecords.length} hits`,
      performance.now() - t0);
    status(`Found ${hitRecords.length} hits — align them in step 3`, 'success');
  } catch (e) {
    status(`Search failed: ${e.message}`, 'error');
    stepStatus('search', 'error');
  } finally {
    $('search-run').disabled = false;
  }
});

// ---------------------------------------------------------------------------
// STEP 3: Align
// ---------------------------------------------------------------------------
$('align-run').addEventListener('click', async () => {
  const source = $('align-input').value;
  const fasta = source === 'search' ? state.searchFasta : state.rawFasta;
  if (!fasta) { status(`No ${source === 'search' ? 'search hits' : 'input sequences'} to align`, 'error'); return; }

  const tool = $('align-tool').value;
  const type = seqType(fasta);
  stepStatus('align', 'running…', 'running');
  status(`Aligning with ${tool}…`);
  $('align-run').disabled = true;
  const t0 = performance.now();
  try {
    const aligned = await runAligner(tool, fasta);
    if (!aligned || !parseFasta(aligned).length) throw new Error('aligner produced no output');
    state.alignedFasta = aligned.trim();
    const recs = parseFasta(state.alignedFasta);
    const alen = recs[0].seq.length;
    $('align-preview').textContent = state.alignedFasta.slice(0, 2500) +
      (state.alignedFasta.length > 2500 ? '\n…' : '');
    $('align-preview').classList.remove('hidden');
    stepStatus('align', `${recs.length} × ${alen}`, 'done');
    recordMeta('Align', tool, 'biowasm', `--auto (${type})`,
      `${recs.length} seq`, `${recs.length} × ${alen} alignment`, performance.now() - t0);
    status(`Aligned ${recs.length} sequences (${alen} columns)`, 'success');
    stepStatus('trim', 'ready');
    stepStatus('tree', 'ready');
  } catch (e) {
    status(`Alignment failed: ${e.message}`, 'error');
    stepStatus('align', 'error');
  } finally {
    $('align-run').disabled = false;
  }
});

// ---------------------------------------------------------------------------
// STEP 4: Trim
// ---------------------------------------------------------------------------
function updateRangeLabels() {
  $('trim-gap-val').textContent = parseFloat($('trim-gap').value).toFixed(2);
  $('trim-cons-val').textContent = parseFloat($('trim-cons').value).toFixed(2);
}
$('trim-gap').addEventListener('input', updateRangeLabels);
$('trim-cons').addEventListener('input', updateRangeLabels);
updateRangeLabels();

$('trim-run').addEventListener('click', () => {
  if (!state.alignedFasta) { status('Align first, or load an aligned FASTA', 'error'); return; }
  const opts = {
    gapThreshold: parseFloat($('trim-gap').value),
    minConservation: parseFloat($('trim-cons').value),
    trimEnds: $('trim-ends').checked,
  };
  try {
    const { fasta, stats } = trimAlignment(state.alignedFasta, opts);
    state.trimmedFasta = fasta;
    $('trim-stats').innerHTML = `
      <div class="stat"><span class="stat-label">Columns before</span><span class="stat-value">${stats.before}</span></div>
      <div class="stat"><span class="stat-label">Columns after</span><span class="stat-value">${stats.after}</span></div>
      <div class="stat"><span class="stat-label">Removed</span><span class="stat-value">${stats.removed}</span></div>
      <div class="stat"><span class="stat-label">Sequences</span><span class="stat-value">${stats.sequences}</span></div>`;
    $('trim-stats').classList.remove('hidden');
    $('trim-preview').textContent = fasta.slice(0, 2500) + (fasta.length > 2500 ? '\n…' : '');
    $('trim-preview').classList.remove('hidden');
    stepStatus('trim', `${stats.after} cols`, 'done');
    recordMeta('Trim', 'trimAl-like (JS)', VERSION,
      `gap≤${opts.gapThreshold}, cons≥${opts.minConservation}, ends=${opts.trimEnds}`,
      `${stats.before} cols`, `${stats.after} cols`, null);
    status(`Trimmed to ${stats.after} columns (removed ${stats.removed})`, 'success');
  } catch (e) {
    status(`Trim failed: ${e.message}`, 'error');
    stepStatus('trim', 'error');
  }
});

// ---------------------------------------------------------------------------
// STEP 5: Tree
// ---------------------------------------------------------------------------
$('tree-run').addEventListener('click', () => {
  const fasta = state.trimmedFasta || state.alignedFasta || state.rawFasta;
  if (!fasta) { status('Load or align sequences first', 'error'); return; }

  const model = $('tree-model').value;
  const type = model === 'auto' ? seqType(fasta) : (model === 'nt' ? 'nt' : 'aa');
  const args = [];
  if (type === 'nt') args.push('-nt');
  if ($('tree-gtr').checked && type === 'nt') args.push('-gtr');
  if (model === 'lg') args.push('-lg');
  if (model === 'wag') args.push('-wag');
  if ($('tree-gamma').checked) args.push('-gamma');
  if ($('tree-nosupport').checked) args.push('-nosupport');

  stepStatus('tree', 'running…', 'running');
  status(`Building tree with FastTree (${type === 'nt' ? 'nucleotide' : 'protein'})…`);
  $('tree-run').disabled = true;

  // Let the status paint before the synchronous WASM call
  setTimeout(() => {
    const t0 = performance.now();
    try {
      state.newick = buildTree(fasta, args);
      const tree = parseNewick(state.newick);
      const stats = treeStats(tree);
      $('newick-output').textContent = state.newick;
      $('tree-results').classList.remove('hidden');
      drawTree($('tree-svg'), tree, { showSupport: !$('tree-nosupport').checked });
      stepStatus('tree', `${stats.tips} tips`, 'done');
      recordMeta('Tree', 'FastTree', fastTreeVersion(), args.join(' ') || '(defaults)',
        `${stats.tips} seq alignment`, `tree, ${stats.tips} tips`, performance.now() - t0);
      status(`Tree inferred: ${stats.tips} tips in ${((performance.now() - t0) / 1000).toFixed(2)}s`, 'success');
      $('tree-svg').scrollIntoView({ behavior: 'smooth', block: 'nearest' });
    } catch (e) {
      status(`Tree inference failed: ${e}`, 'error');
      stepStatus('tree', 'error');
    } finally {
      $('tree-run').disabled = false;
    }
  }, 50);
});

$('newick-copy').addEventListener('click', async () => {
  try {
    await navigator.clipboard.writeText(state.newick);
    status('Newick copied to clipboard', 'success');
  } catch (e) {
    status('Clipboard access denied', 'error');
  }
});
$('newick-download').addEventListener('click', () => {
  const blob = new Blob([state.newick + '\n'], { type: 'text/plain' });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = 'phylo-fam-fate-tree.nwk';
  a.click();
  URL.revokeObjectURL(a.href);
});

// ---------------------------------------------------------------------------
// Metadata export
// ---------------------------------------------------------------------------
$('metadata-export').addEventListener('click', () => {
  const data = {
    tool: 'phylo-fam-fate', version: VERSION,
    generated: new Date().toISOString(),
    sequenceCount: state.records.length,
    srcType: seqType(state.rawFasta || 'ACGT'),
    steps: state.metadata,
  };
  const blob = new Blob([JSON.stringify(data, null, 2)], { type: 'application/json' });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = 'phylo-fam-fate-metadata.json';
  a.click();
  URL.revokeObjectURL(a.href);
});

// ---------------------------------------------------------------------------
// Headless test hook: ?autotest=1 runs a full local pipeline
// ---------------------------------------------------------------------------
if (new URLSearchParams(location.search).has('autotest')) {
  const testFasta = '>a\nACGTACGTACGTACGTACGT\n>b\nACGTACGTACGTACGTACGA\n>c\nACGTACGTACGAACGTACGA\n>d\nACGAACGTACGAACGTACGA\n>e\nTCGAACGTACGAACGTACGA\n';
  $('paste-fasta').value = testFasta;
  setRawFasta(testFasta, 'autotest');
  if (new URLSearchParams(location.search).get('autotest') === 'full') {
    const alignTool = new URLSearchParams(location.search).get('align') || 'kalign';
    $('align-tool').value = alignTool;
    setTimeout(() => $('align-run').click(), 100);
    // after alignment finishes, trim and build the tree
    const wait = setInterval(() => {
      if (state.alignedFasta) {
        clearInterval(wait);
        $('trim-run').click();
        setTimeout(() => $('tree-run').click(), 50);
      }
    }, 200);
  } else {
    setTimeout(() => $('tree-run').click(), 100);
  }
}
