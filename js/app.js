/* Main application: orchestrates the search -> align -> trim -> tree pipeline. */

import { retrieve, parseFasta, toFasta } from './api.js';
import { runAligner, runBlast, runPhmmer, parseBlastHits, parsePhmmerHits, TOOL_VERSIONS } from './tools.js';
import { trimAlignment } from './trim.js';
import { parseNewick, drawTree, treeStats, midpointRoot, ladderize, rootOnOutgroup, madRoot, labelEvents, leafNames } from './tree.js';
import { annotateLabels, fetchPfamDomains, accessionOf, taxonAt, colourGroups, domainColour, RANKS, OTHER } from './annotate.js';
import { toNewick, annotationsTsv, itolColourStrip, itolDomains, eventsTsv } from './export.js';
import { renderMsa } from './msa.js';

const $ = (id) => document.getElementById(id);

// Pipeline state
const state = {
  rawFasta: '',       // original input sequences
  dbFasta: '',        // search database held in memory (custom / phmmer)
  dbIndex: null,      // Map name -> seq for dbFasta (built once, reused by phmmer)
  dbSources: [],      // streamed DB URLs [{url, label}] (BLAST presets / proteome set)
  dbFastaKey: '',     // which dbSources the in-memory dbFasta was downloaded from
  proteomeSet: [],    // reference proteomes chosen for the search set
  dbSource: '',       // human-readable DB source
  dbPreset: 'uniprot-sprot',
  proteomes: [],      // reference proteome catalogue
  searchHits: [],     // structured hit records (one per HSP)
  hitGroups: [],      // hits grouped per subject protein
  expandedHits: new Set(),
  selectedHits: new Set(),
  searchFasta: '',    // sequences from selected hits
  alignedFasta: '',   // aligned sequences
  trimmedFasta: '',   // trimmed alignment
  newick: '',
  tree: null,         // parsed Newick (unrooted, as FastTree wrote it)
  descs: new Map(),   // sequence name -> FASTA description line (OS=, GN=…)
  annotations: new Map(), // leaf label -> annotation (annotate.js)
  pfam: new Map(),    // accession -> Pfam domains (InterPro)
  rootClade: null,    // Set of leaf names: user-chosen root branch
  selectedClade: null,// Set of leaf names highlighted in the tree
  shownTree: null,    // tree as last drawn (rooted/ladderized copy)
  shownNodes: [],     // node index -> node for click handling
  ftVersion: '',
  records: [],
  metadata: [],
  timings: [],        // {step, ms} for the runtime strip
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
  const nav = $(`nav-${id}`);
  if (nav) {
    nav.textContent = text;
    nav.parentElement.className = kind;
  }
}
function recordMeta(step, tool, version, params, inputDesc, outputDesc, ms) {
  state.metadata.push({
    step, tool, version, params,
    input: inputDesc, output: outputDesc,
    time: ms !== null ? `${(ms / 1000).toFixed(2)}s` : '—',
    timestamp: new Date().toISOString(),
  });
  if (ms != null) {
    state.timings.push({ step, ms });
    console.log(`[phylo-fam-fate] ${step} (${tool}): ${(ms / 1000).toFixed(2)}s | ${inputDesc} → ${outputDesc}`);
    debugLog(`${step} (${tool}): ${(ms / 1000).toFixed(2)}s | ${inputDesc} → ${outputDesc}`);
    $('timing-strip').innerHTML = state.timings.map(t =>
      `<span class="timing-item"><strong>${escapeHtml(t.step)}</strong> ${(t.ms / 1000).toFixed(2)}s</span>`).join('');
  }
  renderMetadata();
}
function renderMetadata() {
  const tbody = document.querySelector('#metadata-table tbody');
  tbody.innerHTML = state.metadata.map(m => `
    <tr>
      <td>${escapeHtml(m.step)}</td><td>${escapeHtml(m.tool)}</td><td>${escapeHtml(m.version)}</td>
      <td title="${escapeAttr(m.params)}">${escapeHtml(trunc(m.params, 60))}</td>
      <td>${escapeHtml(m.input)}</td><td>${escapeHtml(m.output)}</td>
      <td>${m.time}</td>
    </tr>`).join('');
}
function escapeHtml(s) { return String(s).replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c])); }
const escapeAttr = escapeHtml;
function download(name, text, type = 'text/plain') {
  const a = document.createElement('a');
  a.href = URL.createObjectURL(new Blob([text], { type }));
  a.download = name;
  a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 0);
}
function trunc(s, n) { return String(s).length > n ? String(s).slice(0, n) + '…' : String(s); }

// ---------------------------------------------------------------------------
// Debug log — errors at the bottom of the pipeline are otherwise only shown in
// the status bar at the top of the page, so they look like "nothing happened".
// ---------------------------------------------------------------------------
const DEBUG_LIMIT = 200;
function debugLog(msg, kind = '') {
  const el = $('debug-log');
  if (!el) return;
  const ts = new Date().toISOString().slice(11, 23);
  const lines = el.textContent ? el.textContent.split('\n') : [];
  lines.push(`${ts}${kind ? ` [${kind}]` : ''} ${msg}`);
  while (lines.length > DEBUG_LIMIT) lines.shift();
  el.textContent = lines.join('\n');
  el.scrollTop = el.scrollHeight;
}
let debugErrors = 0;
function revealDebug() {
  $('debug-panel')?.classList.remove('collapsed');
  $('debug-toggle').textContent = 'Hide';
}
function bumpDebugBadge() {
  debugErrors++;
  $('debug-badge').textContent = debugErrors;
  $('debug-badge').classList.remove('hidden');
}
function showStepError(step, err) {
  const body = $(`step-${step}`)?.querySelector('.step-body');
  if (body) {
    let box = $(`error-${step}`);
    if (!box) {
      box = document.createElement('div');
      box.id = `error-${step}`;
      box.className = 'error-box';
      body.appendChild(box);
    }
    box.textContent = err?.stack || String(err);
    box.classList.remove('hidden');
  }
  debugLog(`${step} ERROR: ${err?.stack || err}`, 'ERROR');
  bumpDebugBadge();
}
function clearStepError(step) { $(`error-${step}`)?.classList.add('hidden'); }

// Capture anything the browser or CDN tools report, even if it never reaches
// a step handler. debugLog does not call console, so this cannot recurse.
window.addEventListener('error', e =>
  debugLog(`${e.message} (${e.filename}:${e.lineno})`, 'PAGE ERROR'));
window.addEventListener('unhandledrejection', e =>
  debugLog(`Unhandled rejection: ${e.reason?.stack || e.reason}`, 'PAGE ERROR'));
for (const level of ['warn', 'error']) {
  const orig = console[level].bind(console);
  console[level] = (...args) => {
    debugLog(args.map(a => a?.stack || String(a)).join(' '), level.toUpperCase());
    if (level === 'error') bumpDebugBadge();
    orig(...args);
  };
}
$('debug-toggle').addEventListener('click', () => {
  const collapsed = $('debug-panel').classList.toggle('collapsed');
  $('debug-toggle').textContent = collapsed ? 'Show' : 'Hide';
});
$('debug-clear').addEventListener('click', () => {
  $('debug-log').textContent = '';
  debugErrors = 0;
  $('debug-badge').classList.add('hidden');
});
$('debug-copy').addEventListener('click', async () => {
  try {
    await navigator.clipboard.writeText($('debug-log').textContent);
    status('Debug log copied to clipboard', 'success');
  } catch (e) {
    debugLog(`Clipboard failed: ${e}`, 'WARN');
  }
});
debugLog(`phylo-fam-fate v${VERSION} starting`);

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
// FastTree worker: inference runs off the main thread so the page stays
// responsive and a run can be cancelled (terminate + respawn the worker).
// ---------------------------------------------------------------------------
let ftWorker = null, ftSeq = 0;
const ftPending = new Map();
function spawnFastTree() {
  ftWorker = new Worker(new URL('./fasttree-worker.js', import.meta.url), { type: 'module' });
  ftWorker.onmessage = ({ data }) => {
    const p = ftPending.get(data.id);
    if (!p) return;
    ftPending.delete(data.id);
    data.ok ? p.resolve(data) : p.reject(new Error(data.error));
  };
  ftWorker.onerror = e => {
    const err = new Error(`FastTree worker error: ${e.message || 'failed to load'}`);
    for (const p of ftPending.values()) p.reject(err);
    ftPending.clear();
  };
}
function ftCall(msg) {
  const id = ++ftSeq;
  return new Promise((resolve, reject) => {
    ftPending.set(id, { resolve, reject });
    ftWorker.postMessage({ ...msg, id });
  });
}
function cancelFastTree() {
  ftWorker?.terminate();
  for (const p of ftPending.values()) p.reject(new Error('cancelled'));
  ftPending.clear();
  spawnFastTree();
}

spawnFastTree();
ftCall({ type: 'version' }).then(({ version }) => {
  state.ftVersion = version;
  $('tree-run').disabled = false;
  status(`FastTree ${version} ready — load sequences to begin`, 'success');
  stepStatus('tree', 'ready');
}).catch(e => {
  status(`Failed to load FastTree WASM: ${e.message}`, 'error');
  showStepError('tree', e);
  debugLog('FastTree WASM init failed — the Build tree button stays disabled', 'FATAL');
});

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

let pasteTimer = 0;
$('paste-fasta').addEventListener('input', () => {
  clearTimeout(pasteTimer);
  pasteTimer = setTimeout(() => setRawFasta($('paste-fasta').value, 'pasted text'), 300);
});

async function loadExample(file, label) {
  try {
    const r = await fetch(`assets/examples/${file}`);
    if (!r.ok) throw new Error(`HTTP ${r.status}`);
    setRawFasta(await r.text(), label);
  } catch (e) {
    status(`Could not load example: ${e.message}`, 'error');
  }
}
$('load-example-met12').addEventListener('click', async () => {
  await loadExample('MET12_P46151.fa', 'example query: yeast MET12 (UniProt P46151)');
  if (!state.records.length) return;
  // A single query is meant for searching, so open step 2 ready to go
  $('search-enable').checked = true;
  $('search-options').classList.remove('hidden');
  $('search-tool').value = 'blastp';
  status('Loaded MET12 (P46151). Next: step 2 → Run search against Swiss-Prot, select hits, then “Use selected for alignment”.', 'success');
  $('step-search').scrollIntoView({ behavior: 'smooth', block: 'start' });
});
$('load-example-aa').addEventListener('click', () => loadExample('proteins.fa', 'example proteins (unaligned)'));
$('load-example-nt').addEventListener('click', () => loadExample('nucleotide_aligned.fa', 'example nucleotide alignment'));

function showMsa(id, fasta, type) {
  renderMsa($(id), fasta, { type, maxCols: 5000 });
  $(id).classList.remove('hidden');
}

function setRawFasta(fasta, source) {
  state.rawFasta = fasta.trim();
  state.records = parseFasta(state.rawFasta);
  state.descs = new Map(state.records.filter(r => r.desc).map(r => [r.name, r.desc]));
  state.annotations = new Map(); state.rootClade = null; state.selectedClade = null;
  state.alignedFasta = '';
  state.trimmedFasta = '';
  state.searchFasta = '';
  state.searchHits = []; state.hitGroups = []; state.expandedHits = new Set();
  state.selectedHits = new Set();
  state.newick = '';
  state.tree = null;
  $('newick-output').textContent = '';
  for (const id of ['tree-results', 'align-preview', 'align-actions', 'trim-preview', 'trim-actions',
    'trim-stats', 'hit-controls', 'hit-results', 'search-preview']) $(id).classList.add('hidden');
  $('align-input').querySelector('option[value="search"]').textContent = 'Search hits';
  $('trim-run').disabled = true;
  ['input', 'align', 'trim', 'tree', 'search'].forEach(clearStepError);

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
  const aligned = new Set(lengths).size === 1 && state.records.length > 1;
  status(`Loaded ${state.records.length} sequences — ` +
    (aligned ? 'already aligned: you can trim or build a tree directly' : 'next, align them (step 3)'), 'success');
  stepStatus('align', aligned ? 'optional' : 'ready');
  stepStatus('trim', 'waiting');
  stepStatus('tree', aligned ? 'ready' : 'needs alignment');
  if (aligned) {
    // An aligned input can go straight to trimming
    state.alignedFasta = state.rawFasta;
    $('trim-run').disabled = false;
    showMsa('align-preview', state.alignedFasta, type);
  }
  stepStatus('search', 'optional');
  syncModelControls();
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
const UNIPROT_FTP = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase';
const SPROT_URL = `${UNIPROT_FTP}/complete/uniprot_sprot.fasta.gz`;
const PROTEOME_README = `${UNIPROT_FTP}/reference_proteomes/README`;
const PDB_URL = 'https://files.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz';

/** Decode a FASTA buffer, transparently decompressing gzip (magic bytes 1f 8b). */
async function decodeFastaBuffer(buf, label) {
  const bytes = new Uint8Array(buf);
  const gzipped = bytes.length > 2 && bytes[0] === 0x1f && bytes[1] === 0x8b;
  let text;
  if (gzipped) {
    if (typeof DecompressionStream === 'undefined') {
      throw new Error('gzip requires a browser with DecompressionStream support');
    }
    const stream = new Blob([buf]).stream().pipeThrough(new DecompressionStream('gzip'));
    text = await new Response(stream).text();
  } else {
    text = new TextDecoder().decode(buf);
  }
  const records = parseFasta(text);
  if (!records.length) throw new Error(`${label} does not contain any FASTA records`);
  return { text, records, gzipped, bytes: bytes.length };
}

/** Fetch + decode a remote FASTA (needed for phmmer, which cannot stream a DB). */
async function fetchDb(url, label) {
  const t0 = performance.now();
  status(`Downloading ${label}…`);
  const resp = await fetch(url);
  if (!resp.ok) throw new Error(`HTTP ${resp.status} from ${new URL(url).host}`);
  const decoded = await decodeFastaBuffer(await resp.arrayBuffer(), label);
  const ms = performance.now() - t0;
  console.log(`[phylo-fam-fate] DB download+decode (${label}): ${(ms / 1000).toFixed(2)}s`);
  return decoded;
}

/** Describe a URL-backed DB (one or more streamed files) in the UI. */
function setDbSources(sources, label, note = '') {
  state.dbSources = sources;
  state.dbFasta = '';
  state.dbFastaKey = '';
  state.dbIndex = null; state.dbDescs = null;
  state.dbSource = label;
  $('search-db-info').innerHTML =
    `<strong>${escapeHtml(label)}</strong>${note ? ' · ' + escapeHtml(note) : ''} · streamed into Biowasm at search time (nothing uploaded)`;
  stepStatus('search', 'db ready');
}
const setDbUrl = (url, label, note) => setDbSources([{ url, label }], label, note);
const dbReady = () => state.dbSources.length > 0 || !!state.dbFasta;

/** Describe an in-memory DB in the UI. */
function setSearchDb({ text, records, gzipped, bytes }, source) {
  state.dbIndex = new Map(records.map(r => [r.name, r.seq]));
  state.dbDescs = new Map(records.filter(r => r.desc).map(r => [r.name, r.desc]));
  state.dbFasta = text;
  state.dbSources = [];
  state.dbFastaKey = '';
  state.dbSource = source;
  const size = gzipped ? `${bytes.toLocaleString()} B gzipped → ${text.length.toLocaleString()} B` : `${text.length.toLocaleString()} B`;
  $('search-db-info').innerHTML = `<strong>${escapeHtml(source)}</strong> · ${records.length.toLocaleString()} sequences · ${size} · ready to search`;
  stepStatus('search', `${records.length} db seq`, 'done');
  status(`Database loaded: ${records.length.toLocaleString()} sequences (${size})`, 'success');
}

// --- Database preset selection ------------------------------------------------
function applyDbPreset() {
  const preset = $('search-db-preset').value;
  state.dbPreset = preset;
  $('search-db-custom').classList.toggle('hidden', preset !== 'custom');
  $('proteome-picker').classList.toggle('hidden', preset !== 'uniprot-reference');
  $('pdb-options').classList.toggle('hidden', preset !== 'pdb');

  if (preset === 'uniprot-sprot') {
    setDbUrl(SPROT_URL, 'UniProt Swiss-Prot',
      `${UNIPROT_FTP.split('/')[2]}, 94 MB gzipped, ~565k reviewed proteins`);
  } else if (preset === 'uniprot-reference') {
    loadProteomeCatalog();
    applyProteomeSet();
  } else if (preset === 'pdb') {
    const proxy = $('search-proxy').value.trim();
    const url = proxy ? proxy + encodeURIComponent(PDB_URL) : PDB_URL;
    setDbUrl(url, 'PDB SEQRES',
      proxy ? 'via CORS proxy' : 'direct URL (CORS-blocked — see note)');
  } else {
    state.dbSources = []; state.dbFasta = ''; state.dbFastaKey = '';
    $('search-db-info').textContent = 'Provide a FASTA database as a local file or a CORS-enabled URL.';
    stepStatus('search', 'choose db');
  }
}
$('search-db-preset').addEventListener('change', applyDbPreset);
$('search-proxy').addEventListener('input', () => { if (state.dbPreset === 'pdb') applyDbPreset(); });

// --- Reference proteome catalogue ---------------------------------------------
async function loadProteomeCatalog() {
  if (state.proteomes.length) { renderProteomeOptions(); return; }
  $('proteome-count').textContent = 'Loading catalogue…';
  try {
    const t = await (await fetch(PROTEOME_README)).text();
    state.proteomes = t.split('\n')
      .filter(l => /^UP0/.test(l))
      .map(l => {
        const f = l.split('\t');
        return { id: f[0], taxid: f[1], division: f[3], count: parseInt(f[4], 10) || 0, name: f[7] || '' };
      })
      .filter(p => p.name);
    renderProteomeOptions();
  } catch (e) {
    $('proteome-count').textContent = `Could not load catalogue: ${e.message}`;
  }
}

function renderProteomeOptions() {
  const q = $('proteome-filter').value.trim().toLowerCase();
  const matched = state.proteomes.filter(p =>
    !q || p.name.toLowerCase().includes(q) || p.id.toLowerCase().includes(q) || p.taxid.includes(q));
  const shown = matched.slice(0, 300);
  $('proteome-select').innerHTML = shown.map(p =>
    `<option value="${p.id}">${escapeHtml(p.name)} · ${p.id} · ${p.count.toLocaleString()} proteins</option>`).join('');
  $('proteome-count').textContent = `${matched.length.toLocaleString()} proteomes` +
    (matched.length > 300 ? ' (showing first 300 — refine filter)' : '');
}
$('proteome-filter').addEventListener('input', renderProteomeOptions);

function proteomeUrl(p) {
  const division = p.division.charAt(0).toUpperCase() + p.division.slice(1);
  return `${UNIPROT_FTP}/reference_proteomes/${division}/${p.id}/${p.id}_${p.taxid}.fasta.gz`;
}

// Above this the in-browser makeblastdb step gets slow and memory-hungry.
const SET_WARN_PROTEINS = 1_500_000;

/** Turn the chosen proteomes (+ optional Swiss-Prot) into the search DB. */
function applyProteomeSet() {
  const set = state.proteomeSet;
  const withSprot = $('proteome-with-sprot').checked;
  $('proteome-set').innerHTML = set.map(p => `
    <li title="${escapeAttr(p.name)}">${escapeHtml(trunc(p.name, 40))} <small>${p.id} · ${p.count.toLocaleString()}</small>
      <button type="button" data-remove="${escapeAttr(p.id)}" aria-label="Remove ${escapeAttr(p.name)}">×</button></li>`).join('');
  const proteins = set.reduce((n, p) => n + p.count, 0) + (withSprot ? 570_000 : 0);
  const summary = $('proteome-set-summary');
  summary.classList.toggle('warn', proteins > SET_WARN_PROTEINS);
  if (!set.length && !withSprot) {
    summary.textContent = 'empty — add one or more proteomes';
    state.dbSources = []; state.dbFasta = ''; state.dbFastaKey = '';
    $('search-db-info').textContent = 'Filter the catalogue, pick a proteome and click “Add to search set”. Repeat to search several at once.';
    stepStatus('search', 'choose proteomes');
    return;
  }
  summary.textContent = `${set.length} proteome${set.length === 1 ? '' : 's'}${withSprot ? ' + Swiss-Prot' : ''} · ~${proteins.toLocaleString()} proteins` +
    (proteins > SET_WARN_PROTEINS ? ' — large: indexing in the browser may be slow' : '');
  const sources = set.map(p => ({ url: proteomeUrl(p), label: `${p.name} (${p.id})`, id: p.id }));
  if (withSprot) sources.unshift({ url: SPROT_URL, label: 'UniProt Swiss-Prot', id: 'sprot' });
  const label = sources.length === 1 ? sources[0].label
    : `${sources.length} databases: ${sources.map(x => x.id).join(', ')}`;
  setDbSources(sources, label, `~${proteins.toLocaleString()} proteins`);
}

$('proteome-add').addEventListener('click', () => {
  const id = $('proteome-select').value;
  const p = state.proteomes.find(x => x.id === id);
  if (!p) { status('Select a reference proteome first', 'error'); return; }
  if (state.proteomeSet.some(x => x.id === id)) { status(`${p.name} is already in the search set`, ''); return; }
  state.proteomeSet.push(p);
  applyProteomeSet();
  status(`Added ${p.name} — ${state.proteomeSet.length} proteome${state.proteomeSet.length === 1 ? '' : 's'} in the search set`, 'success');
});
$('proteome-set').addEventListener('click', e => {
  const id = e.target.closest('[data-remove]')?.dataset.remove;
  if (!id) return;
  state.proteomeSet = state.proteomeSet.filter(p => p.id !== id);
  applyProteomeSet();
});
$('proteome-clear').addEventListener('click', () => { state.proteomeSet = []; applyProteomeSet(); });
$('proteome-with-sprot').addEventListener('change', applyProteomeSet);
// Double-click an option to add it straight away
$('proteome-select').addEventListener('dblclick', () => $('proteome-add').click());

// --- Custom DB file / URL ------------------------------------------------------
$('search-db').addEventListener('change', async () => {
  const file = $('search-db').files[0];
  if (!file) return;
  try {
    status(`Reading database ${file.name}…`);
    setSearchDb(await decodeFastaBuffer(await file.arrayBuffer(), file.name), file.name);
  } catch (e) {
    status(`Database error: ${e.message}`, 'error');
    stepStatus('search', 'db error');
  }
});

$('search-db-url-load').addEventListener('click', async () => {
  const url = $('search-db-url').value.trim();
  if (!url) { status('Enter a database URL first', 'error'); return; }
  $('search-db-url-load').disabled = true;
  try {
    setSearchDb(await fetchDb(url, 'URL'), new URL(url).host + '/' + url.split('/').pop());
    status('Database URL loaded — ready to search', 'success');
  } catch (e) {
    const hint = e.message.toLowerCase().includes('failed to fetch')
      ? ' (the server likely lacks CORS headers)' : '';
    status(`URL database error: ${e.message}${hint}`, 'error');
    stepStatus('search', 'db error');
  } finally {
    $('search-db-url-load').disabled = false;
  }
});

$('search-enable').addEventListener('change', () => {
  $('search-options').classList.toggle('hidden', !$('search-enable').checked);
  if ($('search-enable').checked && !dbReady()) {
    status('Choose a database source for the local search', 'info');
  }
});

// --- Live hit filtering + selection --------------------------------------------
// BLAST reports one row per HSP, so one subject protein can appear several
// times (several local alignments, or hits from several queries). Hits are
// grouped by subject: one selectable row per protein, its HSPs nested below.
// Each group keeps both the full-length subject (from blastdbcmd, or the
// in-memory DB for phmmer) and the HSP-derived fragment; the "Full-length"
// checkbox chooses which one is aligned.
// Group keys are stable numeric ids assigned after sorting (not names, which
// can contain characters awkward in HTML attributes).
const hitKey = g => String(g.id);

function passesCutoffs(h) {
  const maxE = parseFloat($('hit-evalue').value);
  const minScore = parseFloat($('hit-score').value) || 0;
  const minPid = parseFloat($('hit-pident').value) || 0;
  const minLen = parseInt($('hit-length').value, 10) || 0;
  return (!isNaN(maxE) && h.evalue <= maxE)
    && h.score >= minScore
    && (h.pident == null || h.pident >= minPid)
    && (h.length == null || h.length >= minLen);
}
/** A protein passes if its best HSP does. */
const groupPasses = g => passesCutoffs(g.best);
const useFull = () => $('hit-full-length').checked;
const groupSeq = g => (useFull() && g.fullSeq) || g.hspSeq;

/** Combine one subject's HSPs into the longest contiguous stretch of subject
 *  residues that contains the best HSP. Overlapping/adjacent HSPs extend it;
 *  a gap between HSPs can't be filled (BLAST only reports aligned residues). */
function mergeHspSeq(hsps, best) {
  const usable = hsps.filter(h => h.seq && h.sstart <= h.send && h.seq.length === h.send - h.sstart + 1);
  if (!usable.includes(best)) return best.seq;
  const res = new Map();
  for (const h of usable) for (let i = 0; i < h.seq.length; i++) {
    const pos = h.sstart + i;
    if (!res.has(pos)) res.set(pos, h.seq[i]);
  }
  let lo = best.sstart, hi = best.send;
  while (res.has(lo - 1)) lo--;
  while (res.has(hi + 1)) hi++;
  let seq = '';
  for (let p = lo; p <= hi; p++) seq += res.get(p);
  return seq;
}

function groupHits(hits) {
  const byName = new Map();
  for (const h of hits) {            // hits are E-value sorted: first = best
    let g = byName.get(h.name);
    if (!g) { g = { name: h.name, best: h, hsps: [] }; byName.set(h.name, g); }
    g.hsps.push(h);
  }
  const groups = [...byName.values()];
  groups.forEach((g, i) => {
    g.id = i;
    g.queries = [...new Set(g.hsps.map(h => h.query))];
    g.hspSeq = mergeHspSeq(g.hsps, g.best);
    g.fullSeq = g.hsps.find(h => h.fullSeq)?.fullSeq || '';
    g.slen = g.best.slen || (g.fullSeq ? g.fullSeq.length : null);
    const spans = g.hsps.filter(h => h.sstart != null).map(h => [Math.min(h.sstart, h.send), Math.max(h.sstart, h.send)]);
    g.range = spans.length ? [Math.min(...spans.map(s => s[0])), Math.max(...spans.map(s => s[1]))] : null;
  });
  return groups;
}

function fmtE(v) { return v < 0.001 ? v.toExponential(1) : v.toFixed(4); }

function renderHits() {
  const maxHits = parseInt($('search-maxhits').value, 10) || 50;
  const shown = state.hitGroups.filter(groupPasses).slice(0, maxHits);
  const html = shown.map(g => {
    const key = hitKey(g);
    const h = g.best;
    const checked = state.selectedHits.has(key) ? 'checked' : '';
    const n = g.hsps.length;
    const open = state.expandedHits.has(key);
    const toggle = n > 1
      ? ` <button type="button" class="hsp-toggle" data-key="${escapeAttr(key)}" aria-expanded="${open}" title="Show the ${n} alignments to this protein">${n} HSPs ${open ? '▾' : '▸'}</button>`
      : '';
    const queries = g.queries.length > 1 ? `${escapeHtml(trunc(h.query, 16))} +${g.queries.length - 1}` : escapeHtml(h.query);
    let row = `<tr class="hit-group">
      <td><input type="checkbox" class="hit-check" data-key="${escapeAttr(key)}" ${checked}></td>
      <td title="${escapeAttr(g.name)}">${escapeHtml(trunc(g.name, 28))}${toggle}</td>
      <td title="${escapeAttr(g.queries.join(', '))}">${queries}</td>
      <td>${fmtE(h.evalue)}</td>
      <td>${h.score.toFixed(1)}</td>
      <td>${h.pident == null ? '—' : h.pident.toFixed(1)}</td>
      <td title="${escapeAttr(`${groupSeq(g).length} residues will be aligned (${useFull() && g.fullSeq ? 'full-length' : 'HSP region'})`)}">${h.length == null ? '—' : h.length}</td>
      <td>${g.range ? `${g.range[0]}–${g.range[1]}` : '—'}${g.slen ? ` <span class="hint">/ ${g.slen}</span>` : ''}${useFull() && !g.fullSeq ? ' <span class="hsp-only" title="Full-length sequence unavailable; the HSP region will be aligned">HSP only</span>' : ''}</td>
    </tr>`;
    if (open) row += g.hsps.map((x, i) => `<tr class="hsp-row${passesCutoffs(x) ? '' : ' filtered-out'}">
      <td></td>
      <td>↳ HSP ${i + 1}</td>
      <td>${escapeHtml(x.query)}</td>
      <td>${fmtE(x.evalue)}</td>
      <td>${x.score.toFixed(1)}</td>
      <td>${x.pident == null ? '—' : x.pident.toFixed(1)}</td>
      <td>${x.length == null ? '—' : x.length}</td>
      <td>${x.sstart != null ? `${x.sstart}–${x.send}` : '—'}</td>
    </tr>`).join('');
    return row;
  }).join('');
  $('hit-table-body').innerHTML = html || '<tr><td colspan="8" class="hint">No hits pass the current cutoffs</td></tr>';
  $('hit-controls').classList.remove('hidden');
  $('hit-results').classList.remove('hidden');
  updateHitSummary();
}

function updateHitSummary() {
  const chosen = state.hitGroups.filter(g => state.selectedHits.has(hitKey(g)) && groupPasses(g)).length;
  const nHsp = state.searchHits.length, nProt = state.hitGroups.length;
  $('hit-summary').textContent = `${chosen} selected · ${nProt} proteins` + (nHsp > nProt ? ` (${nHsp} HSPs)` : '');
  $('hit-use-selected').textContent = `Use ${chosen} selected for alignment →`;
  $('hit-use-selected').disabled = chosen === 0;
}

$('hit-table-body').addEventListener('change', e => {
  const cb = e.target.closest('.hit-check');
  if (!cb) return;
  if (cb.checked) state.selectedHits.add(cb.dataset.key);
  else state.selectedHits.delete(cb.dataset.key);
  updateHitSummary();
});
$('hit-table-body').addEventListener('click', e => {
  const t = e.target.closest('.hsp-toggle');
  if (!t) return;
  const k = t.dataset.key;
  state.expandedHits.has(k) ? state.expandedHits.delete(k) : state.expandedHits.add(k);
  renderHits();
});

$('hit-full-length').addEventListener('change', renderHits);
['hit-evalue', 'hit-score', 'hit-pident', 'hit-length'].forEach(id => {
  $(id).addEventListener('input', renderHits);
});
$('hit-select-all').addEventListener('click', () => {
  state.hitGroups.filter(groupPasses).forEach(g => state.selectedHits.add(hitKey(g)));
  renderHits();
});
$('hit-select-none').addEventListener('click', () => {
  state.selectedHits.clear(); renderHits();
});

$('hit-use-selected').addEventListener('click', () => {
  const picked = state.hitGroups.filter(g => state.selectedHits.has(hitKey(g)) && groupPasses(g) && groupSeq(g));
  // Groups are already one per subject; queries may share a name with a hit
  // (FastTree rejects duplicate names), so skip those.
  const seen = new Set();
  const out = [];
  if ($('hit-include-query').checked) {
    for (const r of state.records) { if (!seen.has(r.name)) { seen.add(r.name); out.push(r); } }
  }
  for (const g of picked) {
    if (seen.has(g.name)) continue;
    seen.add(g.name);
    out.push({ name: g.name, seq: groupSeq(g) });
  }
  const chosen = out;
  if (!chosen.length) { status('No usable sequences in the selection', 'error'); return; }
  state.searchFasta = toFasta(chosen);
  $('align-input').value = 'search';
  $('align-input').querySelector('option[value="search"]').textContent = `Search hits (${chosen.length})`;
  stepStatus('search', `${chosen.length} selected`, 'done');
  status(`${chosen.length} sequences ready for alignment`, 'success');
  const nFull = useFull() ? picked.filter(g => g.fullSeq).length : 0;
  const merged = picked.filter(g => g.hsps.length > 1 && !(useFull() && g.fullSeq)).length;
  recordMeta('Select hits', 'browser', VERSION, 'cutoffs + manual selection; HSPs grouped per subject',
    `${state.hitGroups.length} proteins (${state.searchHits.length} HSPs)`,
    `${chosen.length} sequences (${picked.length} proteins: ${nFull} full-length, ${picked.length - nFull} HSP region${merged ? `, ${merged} with multiple HSPs merged` : ''}${$('hit-include-query').checked ? ', + queries' : ''})`, null);
});

// --- Run the search ------------------------------------------------------------
$('search-run').addEventListener('click', async () => {
  const tool = $('search-tool').value;
  if (!state.rawFasta) { status('Load input sequences first', 'error'); return; }
  if (!dbReady()) { status('Choose a database source first', 'error'); return; }

  const maxHits = parseInt($('search-maxhits').value, 10) || 50;
  const searchE = parseFloat($('search-evalue').value);
  if (!(searchE > 0)) { status('Search E-value must be a positive number, e.g. 1e-3', 'error'); return; }
  stepStatus('search', 'running…', 'running');
  status(`Running ${tool} locally…`);
  $('search-run').disabled = true;
  const t0 = performance.now();
  try {
    let hits;
    if (tool === 'phmmer') {
      // phmmer needs the DB as a readable file; fetch/decode it if URL-backed.
      // phmmer cannot stream: download every source and concatenate them.
      const key = state.dbSources.map(x => x.url).join('\n');
      if (state.dbSources.length && state.dbFastaKey !== key) {
        const texts = [];
        for (const [i, src] of state.dbSources.entries()) {
          status(`Downloading database ${i + 1}/${state.dbSources.length}: ${src.label}…`);
          texts.push((await fetchDb(src.url, src.label)).text);
        }
        state.dbFasta = texts.map(t => t.trimEnd()).join('\n') + '\n';
        state.dbFastaKey = key; // cache for subsequent phmmer runs
        state.dbIndex = null; state.dbDescs = null;
      }
      const dbFasta = state.dbFasta;
      status('Running phmmer locally…');
      const out = await runPhmmer(state.rawFasta, dbFasta, searchE);
      hits = parsePhmmerHits(out);
      // phmmer does not report sequences; recover them from the in-memory DB.
      state.dbIndex ||= new Map(parseFasta(dbFasta).map(r => [r.name, r.seq]));
      hits.forEach(h => { h.seq = h.fullSeq = state.dbIndex.get(h.name) || ''; });
      if (!state.dbDescs) state.dbDescs = new Map(parseFasta(dbFasta).filter(r => r.desc).map(r => [r.name, r.desc]));
      hits.forEach(h => { const d = state.dbDescs.get(h.name); if (d) state.descs.set(h.name, d); });
      $('search-preview').textContent = out.slice(0, 3000);
    } else {
      // BLAST embeds full subject sequences (sseq) in its tabular output.
      const sources = state.dbSources.length
        ? state.dbSources
        : [{ data: state.dbFasta, label: state.dbSource }];
      const { text: out, sequences, full, titles } = await runBlast(tool, state.rawFasta, sources, maxHits, msg => status(msg), searchE);
      debugLog(`BLAST searched ${sources.length} database(s)${sequences != null ? `, ${sequences.toLocaleString()} sequences` : ''}: ` +
        sources.map(x => x.label).join('; '));
      hits = parseBlastHits(out);
      for (const h of hits) {
        h.fullSeq = full.get(h.acc) || full.get(h.name) || '';
        const title = titles.get(h.acc) || titles.get(h.name);
        if (title) state.descs.set(h.name, title);
      }
      const nFull = new Set(hits.filter(h => h.fullSeq).map(h => h.name)).size;
      const nSubj = new Set(hits.map(h => h.name)).size;
      debugLog(`Full-length sequences retrieved for ${nFull}/${nSubj} subjects (blastdbcmd)`);
      $('search-preview').textContent = out.slice(0, 3000);
    }
    $('search-preview').classList.remove('hidden');

    hits.sort((a, b) => a.evalue - b.evalue);
    state.searchHits = hits;
    state.hitGroups = groupHits(hits);
    state.expandedHits = new Set();
    state.selectedHits = new Set(state.hitGroups.filter(groupPasses).slice(0, maxHits).map(hitKey));
    renderHits();
    const nProt = state.hitGroups.length;
    stepStatus('search', `${nProt} proteins`, 'done');
    const ms = performance.now() - t0;
    recordMeta('Search', tool, TOOL_VERSIONS[tool],
      `maxHits=${maxHits}; search evalue≤${searchE}; display evalue≤${$('hit-evalue').value}, score≥${$('hit-score').value}, pid≥${$('hit-pident').value}, len≥${$('hit-length').value}`,
      `${state.records.length} queries + ${state.dbSource}`, `${nProt} proteins (${hits.length} HSPs)`, ms);
    status(`Found ${nProt} proteins${hits.length > nProt ? ` (${hits.length} HSPs)` : ''} — filter and select which to align`, 'success');
  } catch (e) {
    showStepError('search', `Search failed: ${e.stack || e.message || e}`);
    status(`Search failed: ${e.message}`, 'error');
    stepStatus('search', 'error', 'error');
  } finally {
    $('search-run').disabled = false;
  }
});

// Initialize preset UI on load
applyDbPreset();

// ---------------------------------------------------------------------------
// STEP 3: Align
// ---------------------------------------------------------------------------
$('align-download').addEventListener('click', () => download('phylo-fam-fate-aligned.fa', state.alignedFasta + '\n'));
$('trim-download').addEventListener('click', () => download('phylo-fam-fate-trimmed.fa', state.trimmedFasta));
$('align-run').addEventListener('click', async () => {
  const source = $('align-input').value;
  const fasta = source === 'search' ? state.searchFasta : state.rawFasta;
  if (!fasta) { status(`No ${source === 'search' ? 'search hits' : 'input sequences'} to align`, 'error'); return; }

  // Aligners need at least two sequences. The common trap: one query was
  // searched, hits were selected, but "Use selected" was never clicked, so
  // the aligner would get just the single query.
  const nSeq = parseFasta(fasta).length;
  if (nSeq < 2) {
    const pendingHits = source === 'original' && state.selectedHits.size > 0;
    const msg = pendingHits
      ? `Only 1 input sequence. You selected ${state.selectedHits.size} search hits — click ` +
        '"Use selected for alignment →" in step 2 (or switch Input to "Search hits") first.'
      : `Need at least 2 sequences to align; ${source === 'search' ? 'the selected hits give' : 'the input has'} ${nSeq}. ` +
        'Load more sequences, or run a search (step 2) and use the selected hits.';
    showStepError('align', msg);
    status('Alignment needs at least 2 sequences', 'error');
    stepStatus('align', 'need ≥2 seqs', 'error');
    return;
  }

  const tool = $('align-tool').value;
  const type = seqType(fasta);
  if (tool === 'mafft' && nSeq > 200 &&
      !confirm(`MAFFT E-INS-i aligns every pair of sequences, so ${nSeq} sequences may take a long time in the browser. Continue? (Kalign is much faster.)`)) return;
  stepStatus('align', 'running…', 'running');
  status(`Aligning with ${tool}…`);
  $('align-run').disabled = true;
  clearStepError('align');
  const t0 = performance.now();
  try {
    const aligned = await runAligner(tool, fasta);
    if (!aligned || !parseFasta(aligned).length) throw new Error('aligner produced no output');
    state.alignedFasta = aligned.trim();
    // A new alignment invalidates any earlier trim (the tree step prefers the
    // trimmed alignment, so a stale one would silently be used instead).
    state.trimmedFasta = '';
    $('trim-stats').classList.add('hidden');
    $('trim-preview').classList.add('hidden');
    $('trim-actions').classList.add('hidden');
    const recs = parseFasta(state.alignedFasta);
    const alen = recs[0].seq.length;
    showMsa('align-preview', state.alignedFasta, type);
    $('align-actions').classList.remove('hidden');
    $('trim-run').disabled = false;
    stepStatus('align', `${recs.length} × ${alen}`, 'done');
    recordMeta('Align', tool, TOOL_VERSIONS[tool],
      tool === 'muscle' ? `-align (defaults; ${type})`
        : tool === 'mafft' ? `E-INS-i: tbfast + dvtditr, 16 iterations (≈ mafft --genafpair --ep 0 --maxiterate 16; ${type})`
        : `-i -o (defaults; ${type})`,
      `${recs.length} seq`, `${recs.length} × ${alen} alignment`, performance.now() - t0);
    status(`Aligned ${recs.length} sequences (${alen} columns)`, 'success');
    stepStatus('trim', 'ready');
    stepStatus('tree', 'ready');
    $('step-trim').scrollIntoView({ behavior: 'smooth', block: 'nearest' });
  } catch (e) {
    showStepError('align', `Alignment failed: ${e.stack || e.message || e}`);
    status(`Alignment failed: ${e.message}`, 'error');
    stepStatus('align', 'error', 'error');
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
  const t0 = performance.now();
  try {
    const { fasta, stats } = trimAlignment(state.alignedFasta, opts);
    const trimMs = performance.now() - t0;
    state.trimmedFasta = fasta;
    $('trim-stats').innerHTML = `
      <div class="stat"><span class="stat-label">Columns before</span><span class="stat-value">${stats.before}</span></div>
      <div class="stat"><span class="stat-label">Columns after</span><span class="stat-value">${stats.after}</span></div>
      <div class="stat"><span class="stat-label">Removed</span><span class="stat-value">${stats.removed}</span></div>
      <div class="stat"><span class="stat-label">Sequences</span><span class="stat-value">${stats.sequences}</span></div>`;
    $('trim-stats').classList.remove('hidden');
    showMsa('trim-preview', fasta, seqType(fasta));
    $('trim-actions').classList.remove('hidden');
    clearStepError('trim');
    stepStatus('trim', `${stats.after} cols`, 'done');
    recordMeta('Trim', 'trimAl-like (JS)', VERSION,
      `gap≤${opts.gapThreshold}, cons≥${opts.minConservation}, ends=${opts.trimEnds}`,
      `${stats.before} cols`, `${stats.after} cols`, trimMs);
    status(`Trimmed to ${stats.after} columns (removed ${stats.removed})`, 'success');
  } catch (e) {
    showStepError('trim', `Trim failed: ${e.stack || e.message || e}`);
    status(`Trim failed: ${e.message}`, 'error');
    stepStatus('trim', 'error', 'error');
  }
});

// ---------------------------------------------------------------------------
// STEP 5: Tree
// ---------------------------------------------------------------------------
// Rooting is stored as a leaf-name split, which identifies an edge in any
// orientation of the same tree, so it survives re-rendering from the
// unrooted original.
let madCache = null;   // {tree: state.tree, result}
function rootedTree() {
  const t = structuredClone(state.tree);
  const mode = $('tv-root').value;
  $('tv-outgroup').classList.toggle('hidden', mode !== 'outgroup');
  if (mode === 'none') return { t, note: 'unrooted (as inferred)' };
  if (mode === 'branch' && state.rootClade) {
    const r = rootOnOutgroup(t, state.rootClade);
    return { t: r.tree, note: `rooted on chosen branch (${state.rootClade.size} tips)` };
  }
  if (mode === 'outgroup') {
    const taxon = $('tv-outgroup').value;
    const [rank, name] = taxon.split('\t');
    if (!name) return { t: midpointRoot(t), note: 'choose an outgroup taxon (annotate first) — midpoint shown' };
    const set = new Set(leafNames(t).filter(l => taxonAt(state.annotations.get(l), rank) === name));
    const r = rootOnOutgroup(t, set);
    if (r.invalid) return { t: midpointRoot(t), note: `${name} covers ${set.size ? 'every' : 'no'} leaf — midpoint shown` };
    return { t: r.tree, note: r.exact ? `rooted on ${name} (${set.size} tips, monophyletic)`
      : `⚠ ${name} is not monophyletic: rooted on the best split (${r.impurity} tip${r.impurity === 1 ? '' : 's'} on the wrong side)` };
  }
  if (mode === 'mad') {
    if (madCache?.tree !== state.tree) madCache = { tree: state.tree, result: madRoot(state.tree) };
    const r = madCache.result;
    if (!r) return { t: midpointRoot(t), note: 'MAD needs 3–900 tips — midpoint shown' };
    return { t: structuredClone(r.tree), note: `MAD root (deviation ${r.score.toFixed(3)}, ambiguity ${r.ambiguity.toFixed(2)}${r.ambiguity > 0.95 ? ' — nearly tied with another root' : ''})` };
  }
  return { t: midpointRoot(t), note: 'midpoint root' };
}

const ann = leaf => state.annotations.get(leaf.name);
function leafText(leaf) {
  const a = ann(leaf);
  switch ($('tv-label').value) {
    case 'species': return a?.species || leaf.name;
    case 'id-species': return a?.species ? `${leaf.name} ${a.species}` : leaf.name;
    case 'gene': return a?.genes?.[0] ? `${a.genes[0]} (${a.species || leaf.name})` : leaf.name;
    case 'protein': return a?.protein ? `${a.protein} · ${a.species || ''}` : leaf.name;
    default: return leaf.name;
  }
}

/** Domain list for a leaf: InterPro Pfam if fetched, else UniProt features. */
function leafDomains(name) {
  const a = state.annotations.get(name);
  const acc = a?.accession || accessionOf(name);
  const pf = state.pfam.get(acc);
  return pf?.length ? pf : a?.domains || [];
}

function renderTree() {
  if (!state.tree) return;
  const { t, note } = rootedTree();
  if ($('tv-ladderize').checked) ladderize(t);
  const rooted = $('tv-root').value !== 'none';

  // Colour groups at the chosen rank
  const rank = $('tv-color').value;
  let colourOf = null, legend = '';
  if (rank && state.annotations.size) {
    const counts = new Map();
    for (const l of leafNames(t)) {
      const g = taxonAt(state.annotations.get(l), rank) || 'unknown';
      counts.set(g, (counts.get(g) || 0) + 1);
    }
    const { map, sorted } = colourGroups(counts);
    map.set('unknown', null);
    colourOf = leaf => map.get(taxonAt(ann(leaf), rank) || 'unknown') || null;
    legend = `<strong>${escapeHtml(rank)}</strong>` + sorted.map(([g, n]) =>
      `<span class="lg-item"><i style="background:${g === 'unknown' ? '#e2e8f0' : map.get(g) || OTHER}"></i>${escapeHtml(g)} <small>${n}</small></span>`).join('');
  }

  // Duplication / speciation labels (needs a rooted tree and species)
  let events = null;
  if ($('tv-events').checked) {
    if (!rooted) legend += '<span class="lg-note">Duplications need a rooted tree.</span>';
    else if (!state.annotations.size) legend += '<span class="lg-note">Duplications need species: click “Annotate from UniProt”.</span>';
    else {
      events = labelEvents(t, leaf => ann(leaf)?.taxid || ann(leaf)?.species || null);
      legend += `<span class="lg-item"><i class="lg-dup"></i>duplication <small>${events.dup}</small></span>` +
        `<span class="lg-item"><i class="lg-spec"></i>speciation <small>${events.spec}</small></span>`;
    }
  }

  const showDomains = $('tv-domains').checked;
  const domainKeys = new Map();
  const { matches, nodes } = drawTree($('tree-svg'), t, {
    showSupport: $('tv-support').checked && !$('tree-nosupport').checked,
    cladogram: $('tv-cladogram').checked,
    highlight: $('tv-find').value,
    leafText,
    leafColor: colourOf,
    showEvents: !!events,
    selected: state.selectedClade ? n => n.children.length && sameSet(leafNames(n), state.selectedClade) : null,
    domains: showDomains ? leaf => {
      const doms = leafDomains(leaf.name);
      const len = ann(leaf)?.length || (doms.length ? Math.max(...doms.map(d => d.end)) : 0);
      if (!len) return null;
      return { len, doms: doms.map(d => {
        const key = d.id || d.name;
        domainKeys.set(key, d.name);
        return { ...d, color: domainColour(key), label: `${d.name}${d.id ? ` (${d.id})` : ''}` };
      }) };
    } : null,
  });
  if (showDomains && domainKeys.size) {
    legend += '<strong>domains</strong>' + [...domainKeys].map(([k, n]) =>
      `<span class="lg-item"><i style="background:${domainColour(k)}"></i>${escapeHtml(trunc(n, 34))}${k !== n ? ` <small>${escapeHtml(k)}</small>` : ''}</span>`).join('');
  }
  state.shownTree = t;
  state.shownNodes = nodes;
  state.shownEvents = events;
  $('tree-legend').innerHTML = `<span class="lg-note">${escapeHtml(note)}</span>` + legend;
  $('tree-legend').classList.remove('hidden');
  const q = $('tv-find').value.trim();
  $('tv-find-count').textContent = q ? `${matches} match${matches === 1 ? '' : 'es'}` : '';
}
function sameSet(arr, set) { return arr.length === set.size && arr.every(x => set.has(x)); }

['tv-root', 'tv-outgroup', 'tv-ladderize', 'tv-cladogram', 'tv-support', 'tv-label', 'tv-color', 'tv-events']
  .forEach(id => $(id).addEventListener('change', () => { hidePopup(); renderTree(); }));
$('tv-find').addEventListener('input', renderTree);
let resizeTimer = 0;
window.addEventListener('resize', () => { clearTimeout(resizeTimer); resizeTimer = setTimeout(renderTree, 150); });

// --- Annotation --------------------------------------------------------------
async function annotateTree({ quiet = false } = {}) {
  if (!state.tree) return;
  const labels = leafNames(state.tree);
  const btn = $('tv-annotate');
  btn.disabled = true;
  const t0 = performance.now();
  try {
    const say = m => { $('tv-annot-status').textContent = m; };
    say('Annotating…');
    state.annotations = await annotateLabels(labels, state.descs, say);
    const n = state.annotations.size;
    const withLineage = [...state.annotations.values()].filter(a => a.lineage.length).length;
    say(`${n}/${labels.length} leaves annotated${withLineage < n ? ` (${n - withLineage} from FASTA headers, species only)` : ''}`);
    populateOutgroups();
    if (n && !$('tv-color').value) $('tv-color').value = withLineage ? 'phylum' : 'genus';
    recordMeta('Annotate', 'UniProt REST', 'uniprotkb/accessions', 'organism, lineage, names, Pfam, domains, function',
      `${labels.length} leaves`, `${n} annotated`, performance.now() - t0);
    debugLog(`Annotated ${n}/${labels.length} leaves (${withLineage} with lineage)`);
    if ($('tv-domains').checked) await loadDomains();
    renderTree();
  } catch (e) {
    $('tv-annot-status').textContent = `Annotation failed: ${e.message}`;
    if (!quiet) status(`Annotation failed: ${e.message}`, 'error');
    debugLog(`Annotation failed: ${e.stack || e}`, 'error');
  } finally {
    btn.disabled = false;
  }
}
$('tv-annotate').addEventListener('click', () => annotateTree());

async function loadDomains() {
  if (!state.tree) return;
  const accs = leafNames(state.tree).map(l => state.annotations.get(l)?.accession || accessionOf(l)).filter(Boolean);
  if (!accs.length) { $('tv-annot-status').textContent = 'No UniProt accessions: domain track unavailable'; return; }
  if (accs.length > 600) { $('tv-annot-status').textContent = `${accs.length} proteins: using UniProt domain features only (InterPro limited to 600)`; return; }
  const t0 = performance.now();
  const got = await fetchPfamDomains(accs, m => { $('tv-annot-status').textContent = m; });
  for (const [a, d] of got) state.pfam.set(a, d);
  $('tv-annot-status').textContent = `Pfam domains for ${got.size}/${accs.length} proteins (InterPro)`;
  recordMeta('Domains', 'InterPro API', 'entry/pfam', 'Pfam per protein', `${accs.length} proteins`, `${got.size} with data`, performance.now() - t0);
}
$('tv-domains').addEventListener('change', async () => {
  if ($('tv-domains').checked) {
    if (!state.annotations.size) await annotateTree();
    else await loadDomains();
  }
  renderTree();
});

function populateOutgroups() {
  const sel = $('tv-outgroup');
  const prev = sel.value;
  const opts = [];
  const labels = leafNames(state.tree);
  for (const rank of RANKS) {
    const counts = new Map();
    for (const l of labels) {
      const g = taxonAt(state.annotations.get(l), rank);
      if (g) counts.set(g, (counts.get(g) || 0) + 1);
    }
    // An outgroup must leave something on the other side.
    const usable = [...counts].filter(([, c]) => c < labels.length).sort((a, b) => b[1] - a[1]);
    if (usable.length) opts.push(`<optgroup label="${escapeAttr(rank)}">` + usable.map(([g, c]) =>
      `<option value="${escapeAttr(`${rank}\t${g}`)}">${escapeHtml(g)} (${c})</option>`).join('') + '</optgroup>');
  }
  sel.innerHTML = opts.length ? `<option value="">— choose taxon —</option>${opts.join('')}` : '<option value="">annotate first</option>';
  if ([...sel.options].some(o => o.value === prev)) sel.value = prev;
}

// --- Click: leaf details / branch actions -----------------------------------
function hidePopup() { $('tree-popup').classList.add('hidden'); }
function showPopup(html, evt) {
  const pop = $('tree-popup'), box = $('tree-scroll');
  pop.innerHTML = `<button class="popup-close" type="button" aria-label="Close">×</button>${html}`;
  pop.classList.remove('hidden');
  const r = box.getBoundingClientRect();
  let x = evt.clientX - r.left + box.scrollLeft + 12, y = evt.clientY - r.top + box.scrollTop + 12;
  const w = pop.offsetWidth, h = pop.offsetHeight;
  if (x + w > box.scrollLeft + box.clientWidth - 8) x = Math.max(box.scrollLeft + 8, evt.clientX - r.left + box.scrollLeft - w - 12);
  if (y + h > box.scrollTop + box.clientHeight - 8) y = Math.max(box.scrollTop + 8, box.scrollTop + box.clientHeight - h - 8);
  pop.style.left = `${x}px`; pop.style.top = `${y}px`;
}

function domainSvg(len, doms) {
  if (!len) return '';
  const W = 300, sx = W / len;
  const rects = doms.map(d => `<rect x="${((d.start - 1) * sx).toFixed(1)}" y="3" width="${Math.max(2, (d.end - d.start + 1) * sx).toFixed(1)}" height="12" rx="3" fill="${domainColour(d.id || d.name)}"><title>${escapeHtml(`${d.name} ${d.start}–${d.end}`)}</title></rect>`).join('');
  return `<svg class="pop-domains" viewBox="0 0 ${W} 30" width="${W}" height="30"><line x1="0" y1="9" x2="${W}" y2="9" stroke="#94a3b8" stroke-width="2"/>${rects}<text x="0" y="28" font-size="9" fill="#64748b">1</text><text x="${W}" y="28" font-size="9" fill="#64748b" text-anchor="end">${len} aa</text></svg>`;
}

async function leafPopup(leaf, evt) {
  const a = ann(leaf);
  const acc = a?.accession || accessionOf(leaf.name);
  const render = () => {
    const doms = leafDomains(leaf.name);
    const lineage = a?.lineage?.length
      ? a.lineage.filter(x => RANKS.includes(x.rank) && x.rank !== 'species').map(x => `<span title="${escapeAttr(x.rank)}">${escapeHtml(x.name)}</span>`).join(' › ')
      : '';
    showPopup(`
      <div class="pop-title">${escapeHtml(a?.protein || leaf.name)}</div>
      <div class="pop-sub">${escapeHtml(leaf.name)}${a?.genes?.length ? ` · gene <b>${escapeHtml(a.genes.slice(0, 3).join(', '))}</b>` : ''}</div>
      ${a ? `<div class="pop-species"><i>${escapeHtml(a.species || a.organism)}</i>${a.taxid ? ` <small>taxid ${escapeHtml(a.taxid)}</small>` : ''}</div>` : '<p class="hint">No annotation for this leaf. Click “Annotate from UniProt”, or label leaves with UniProt accessions.</p>'}
      ${lineage ? `<div class="pop-lineage">${lineage}</div>` : ''}
      ${doms.length || a?.length ? `<div class="pop-sec">Domains${state.pfam.get(acc)?.length ? ' (Pfam via InterPro)' : a?.domains?.length ? ' (UniProt)' : ''}</div>${domainSvg(a?.length || Math.max(0, ...doms.map(d => d.end)), doms)}
        <div class="pop-domlist">${doms.map(d => `<span><i style="background:${domainColour(d.id || d.name)}"></i>${escapeHtml(d.name)} ${d.start}–${d.end}</span>`).join('') || '<span class="hint">none reported</span>'}</div>` : ''}
      ${a?.pfam?.length && !doms.length ? `<div class="pop-sec">Pfam</div><div>${a.pfam.map(escapeHtml).join(', ')}</div>` : ''}
      ${a?.func ? `<div class="pop-sec">Function</div><div class="pop-func">${escapeHtml(trunc(a.func, 600))}</div>` : ''}
      <div class="pop-foot">${a?.score != null ? `Annotation score ${a.score}/5 · ` : ''}${a?.source ? `${escapeHtml(a.source)} · ` : ''}${acc ? `<a href="https://www.uniprot.org/uniprotkb/${encodeURIComponent(acc)}/entry" target="_blank" rel="noopener">UniProt ${escapeHtml(acc)} ↗</a>` : ''}</div>`, evt);
  };
  render();
  // Fetch this protein's Pfam domains on demand if not already loaded.
  if (acc && !state.pfam.has(acc)) {
    const got = await fetchPfamDomains([acc]);
    if (got.has(acc)) { state.pfam.set(acc, got.get(acc)); if (!$('tree-popup').classList.contains('hidden')) render(); }
  }
}

function branchPopup(node, evt) {
  const names = leafNames(node);
  const species = new Set(names.map(n => state.annotations.get(n)?.species).filter(Boolean));
  const ev = node.event === 'dup' ? 'Duplication node' : node.event === 'spec' ? 'Speciation node' : '';
  showPopup(`
    <div class="pop-title">Clade of ${names.length} tip${names.length === 1 ? '' : 's'}</div>
    <div class="pop-sub">${node.support != null ? `support ${escapeHtml(node.supportLabel || String(node.support))} · ` : ''}branch length ${(+node.length).toPrecision(3)}${ev ? ` · <b>${ev}</b>` : ''}</div>
    ${species.size ? `<div class="pop-species">${species.size} species: <i>${[...species].slice(0, 6).map(escapeHtml).join(', ')}${species.size > 6 ? '…' : ''}</i></div>` : ''}
    <div class="pop-actions">
      <button class="btn small primary" type="button" data-act="root">Root on this branch</button>
      <button class="btn small" type="button" data-act="select">${state.selectedClade && sameSet(names, state.selectedClade) ? 'Clear highlight' : 'Highlight clade'}</button>
      <button class="btn small ghost" type="button" data-act="names">Copy tip names</button>
    </div>`, evt);
  $('tree-popup').dataset.clade = JSON.stringify(names);
}

$('tree-svg').addEventListener('click', e => {
  const el = e.target.closest('[data-n]');
  if (!el) { hidePopup(); return; }
  const node = state.shownNodes[+el.dataset.n];
  if (!node) return;
  if (!node.children.length && !el.classList.contains('t-hitarea')) leafPopup(node, e);
  else branchPopup(node, e);
});
$('tree-popup').addEventListener('click', async e => {
  if (e.target.closest('.popup-close')) { hidePopup(); return; }
  const act = e.target.closest('[data-act]')?.dataset.act;
  if (!act) return;
  const names = JSON.parse($('tree-popup').dataset.clade || '[]');
  if (act === 'root') {
    state.rootClade = new Set(names);
    const opt = $('tv-root').querySelector('option[value="branch"]');
    opt.disabled = false;
    $('tv-root').value = 'branch';
  } else if (act === 'select') {
    state.selectedClade = state.selectedClade && sameSet(names, state.selectedClade) ? null : new Set(names);
  } else if (act === 'names') {
    try { await navigator.clipboard.writeText(names.join('\n')); status(`Copied ${names.length} tip names`, 'success'); }
    catch { status('Clipboard access denied', 'error'); }
    return;
  }
  hidePopup();
  renderTree();
});
document.addEventListener('keydown', e => { if (e.key === 'Escape') hidePopup(); });

function syncModelControls() {
  const model = $('tree-model').value;
  const fasta = state.trimmedFasta || state.alignedFasta || state.rawFasta;
  const type = model === 'auto' ? (fasta ? seqType(fasta) : 'nt') : (model === 'nt' ? 'nt' : 'aa');
  $('tree-gtr').disabled = type !== 'nt';
  $('tree-gtr-wrap').style.opacity = type === 'nt' ? '' : '0.45';
}
$('tree-model').addEventListener('change', syncModelControls);

$('tree-cancel').addEventListener('click', () => {
  cancelFastTree();
  debugLog('Tree inference cancelled by user');
});

$('tree-run').addEventListener('click', async () => {
  const source = state.trimmedFasta ? 'trimmed' : (state.alignedFasta ? 'aligned' : 'raw');
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

  // FastTree requires an alignment. Give a clear, visible message instead of
  // letting the crate's length validator fail at the top of the page.
  const recs = parseFasta(fasta);
  const lengths = [...new Set(recs.map(r => r.seq.length))];
  if (lengths.length > 1) {
    const hint = `Sequences have unequal lengths (${lengths.slice(0, 6).join(', ')}${lengths.length > 6 ? ', …' : ''}). ` +
      'FastTree needs an alignment — run step 3 (Align) first, or load an already-aligned FASTA.';
    showStepError('tree', hint);
    status('Tree needs an alignment — align first', 'error');
    stepStatus('tree', 'needs alignment', 'error');
    return;
  }
  const dupes = recs.map(r => r.name).filter((n, i, a) => a.indexOf(n) !== i);
  if (dupes.length) {
    showStepError('tree', `Duplicate sequence names: ${[...new Set(dupes)].slice(0, 5).join(', ')}. FastTree requires unique names.`);
    stepStatus('tree', 'duplicate names', 'error');
    return;
  }

  clearStepError('tree');
  debugLog(`buildTree: source=${source}, seqs=${recs.length}, length=${lengths[0]}, ` +
    `args=[${args.join(' ') || 'defaults'}], model=${model}, type=${type}, version=${state.ftVersion}, precision=${$('tree-precision').value}`);
  stepStatus('tree', 'running…', 'running');
  status(`Building tree with FastTree (${type === 'nt' ? 'nucleotide' : 'protein'}, ${recs.length} × ${lengths[0]})…`);
  $('tree-run').disabled = true;
  $('tree-cancel').classList.remove('hidden');
  $('tree-elapsed').classList.remove('hidden');
  const t0 = performance.now();
  const tick = setInterval(() => {
    $('tree-elapsed').textContent = `${((performance.now() - t0) / 1000).toFixed(1)}s`;
  }, 100);

  try {
    const precision = $('tree-precision').value;
    const { newick, ms } = await ftCall({ fasta, args, precision });
    const stats = adoptTree(newick);
    stepStatus('tree', `${stats.tips} tips`, 'done');
    recordMeta('Tree', 'FastTree', `${state.ftVersion} (${precision} precision)`, args.join(' ') || '(defaults)',
      `${recs.length} × ${lengths[0]} ${source} alignment`, `tree, ${stats.tips} tips`, ms);
    debugLog(`Tree OK: ${stats.tips} tips, ${newick.length} Newick chars`);
    status(`Tree inferred: ${stats.tips} tips in ${(ms / 1000).toFixed(2)}s (${precision} precision)`, 'success');
    $('tree-results').scrollIntoView({ behavior: 'smooth', block: 'start' });
  } catch (e) {
    if (e.message === 'cancelled') {
      status('Tree inference cancelled', '');
      stepStatus('tree', 'cancelled');
    } else {
      showStepError('tree', `Tree inference failed: ${e.message}`);
      status(`Tree inference failed: ${e.message}`, 'error');
      stepStatus('tree', 'error', 'error');
    }
  } finally {
    clearInterval(tick);
    $('tree-run').disabled = false;
    $('tree-cancel').classList.add('hidden');
  }
});

/** Show a new tree (built here or loaded from a file) in the viewer. */
function adoptTree(newick) {
  const tree = parseNewick(newick.trim());
  const stats = treeStats(tree);
  if (stats.tips < 2) throw new Error('not a Newick tree with at least 2 tips');
  state.newick = newick.trim();
  state.tree = tree;
  // A new tree invalidates branch-based choices from the previous one.
  state.rootClade = null; state.selectedClade = null; madCache = null;
  $('tv-root').querySelector('option[value="branch"]').disabled = true;
  if ($('tv-root').value === 'branch') $('tv-root').value = 'midpoint';
  hidePopup();
  $('newick-output').textContent = state.newick;
  $('tree-results').classList.remove('hidden');
  renderTree();
  // Annotate automatically when leaves look like UniProt records or carry
  // OS= descriptions; keep any annotations already fetched for these labels.
  const labels = stats.names;
  if (labels.some(l => accessionOf(l) || /\bOS=/.test(state.descs.get(l) || ''))) {
    if (!labels.every(l => state.annotations.has(l))) annotateTree({ quiet: true });
  }
  return stats;
}

$('tree-load').addEventListener('change', async e => {
  const file = e.target.files?.[0];
  e.target.value = '';
  if (!file) return;
  try {
    const text = await file.text();
    // Take the first tree if the file holds several (e.g. bootstrap sets)
    const first = text.trim().split(/;\s*/)[0] + ';';
    const stats = adoptTree(first);
    const support = /\)[0-9.]+\/[0-9.]+/.test(first) ? 'SH-aLRT/UFBoot' : /\)[0-9.]+/.test(first) ? 'single values' : 'none';
    stepStatus('tree', `${stats.tips} tips (loaded)`, 'done');
    recordMeta('Tree', 'loaded file', '—', file.name, `${(file.size / 1024).toFixed(1)} kB Newick`, `tree, ${stats.tips} tips; support: ${support}`, 0);
    const known = state.records.length ? stats.names.filter(n => state.records.some(r => r.name === n)).length : 0;
    status(`Loaded ${file.name}: ${stats.tips} tips` + (state.records.length ? ` (${known} match loaded sequences)` : ''), 'success');
    $('tree-results').scrollIntoView({ behavior: 'smooth', block: 'start' });
  } catch (err) {
    showStepError('tree', `Could not read ${file.name} as Newick: ${err.message}`);
    status(`Could not load tree: ${err.message}`, 'error');
  }
});

/* Hand-off: a self-contained shell script with the alignment embedded, which
   runs famfate.sh (IQ-TREE / RAxML-NG) via Apptainer, Docker, a local
   install, or the ML programs directly. */
$('handoff-download').addEventListener('click', () => {
  const fasta = state.trimmedFasta || state.alignedFasta || state.rawFasta;
  if (!fasta) { status('Align sequences first (step 3): the script needs an alignment', 'error'); return; }
  const recs = parseFasta(fasta);
  const lens = new Set(recs.map(r => r.seq.length));
  if (lens.size !== 1) { status('Sequences are not aligned yet: run step 3 first', 'error'); return; }
  const source = state.trimmedFasta ? 'trimmed' : state.alignedFasta ? 'aligned' : 'input';
  const methods = $('handoff-method').value;
  const boot = Math.max(100, parseInt($('handoff-boot').value, 10) || 1000);
  const tag = methods.split(',')[0];
  const body = recs.map(r => `>${r.name}\n${r.seq}`).join('\n');
  const script = handoffScript({ methods, boot, tag, body, n: recs.length, len: [...lens][0], source, type: seqType(fasta) });
  download(`run_famfate_${tag}.sh`, script, 'text/x-shellscript');
  recordMeta('Hand-off', 'famfate.sh', VERSION, `--tree ${methods} --bootstrap ${boot}`, `${recs.length} × ${[...lens][0]} ${source} alignment`, `run_famfate_${tag}.sh`, 0);
  status(`Saved run_famfate_${tag}.sh: run "bash run_famfate_${tag}.sh", then load famfate_${tag}/family.nwk here`, 'success');
});

function handoffScript({ methods, boot, tag, body, n, len, source, type }) {
  const iqOnly = methods === 'iqtree', rxOnly = methods === 'raxml-ng';
  return `#!/usr/bin/env bash
# phylo-fam-fate hand-off: ${methods} on ${n} sequences × ${len} columns (${source} alignment)
# Generated ${new Date().toISOString()} by ${location.origin}${location.pathname}
#
#   bash run_famfate_${tag}.sh [THREADS]
#
# Writes famfate_${tag}/ with family.nwk (the ${tag} tree), family.rooted.nwk,
# SVG, annotations and iTOL files. Load family.nwk in the website's step 5,
# "Load a tree (Newick)", to explore it with rooting, colours, domains and
# duplications. Runs famfate.sh from (first found): a local install, Apptainer/
# Singularity, or Docker (image ghcr.io/hyphaltip/phylo-fam-fate-cli), else
# ${iqOnly ? 'IQ-TREE' : rxOnly ? 'RAxML-NG' : 'IQ-TREE / RAxML-NG'} directly if installed.
set -euo pipefail
THREADS=\${1:-$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)}
IMAGE=\${FAMFATE_IMAGE:-ghcr.io/hyphaltip/phylo-fam-fate-cli:edge}
OUT=famfate_${tag}

cat > family.aln.fa <<'FAMFATE_ALIGNMENT_EOF'
${body}
FAMFATE_ALIGNMENT_EOF

# The browser already trimmed/aligned this, so don't trim again (--gap none).
ARGS=(--aligned family.aln.fa --gap none --tree ${methods} --bootstrap ${boot} -T "$THREADS" -o "$OUT")

if command -v famfate.sh >/dev/null 2>&1; then
  famfate.sh "\${ARGS[@]}"
elif A=$(command -v apptainer || command -v singularity); then
  "$A" exec "docker://$IMAGE" famfate.sh "\${ARGS[@]}"
elif command -v docker >/dev/null 2>&1; then
  docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":/data "$IMAGE" famfate.sh "\${ARGS[@]}"
${rxOnly ? `elif command -v raxml-ng >/dev/null 2>&1; then
  mkdir -p "$OUT"
  raxml-ng --all --msa family.aln.fa --model ${'$'}{RAXML_MODEL:-${type === 'nt' ? 'GTR+G4' : 'LG+G4'}} --prefix "$OUT/family" --threads "$THREADS" --seed 12345 --bs-trees "autoMRE{${boot}}" --bs-metric tbe
  for f in "$OUT/family.raxml.supportTBE" "$OUT/family.raxml.support" "$OUT/family.raxml.bestTree"; do
    [ -s "$f" ] && { cp "$f" "$OUT/family.nwk"; break; }
  done` :
  `elif IQ=$(command -v iqtree3 || command -v iqtree2 || command -v iqtree); then
  mkdir -p "$OUT"
  "$IQ" -s family.aln.fa -m MFP -B ${Math.max(1000, boot)} -alrt 1000 -T "$THREADS" --seed 12345 --prefix "$OUT/family"
  cp "$OUT/family.treefile" "$OUT/family.nwk"`}
else
  echo "Need one of: famfate.sh, apptainer/singularity, docker, ${rxOnly ? 'raxml-ng' : 'iqtree3/iqtree2'}" >&2
  echo "Install e.g.: mamba install -c conda-forge -c bioconda ${rxOnly ? 'raxml-ng' : 'iqtree'}" >&2
  exit 1
fi
echo "Done. Load $OUT/family.nwk in phylo-fam-fate (step 5 → Load a tree)."
`;
}

$('newick-copy').addEventListener('click', async () => {
  try {
    await navigator.clipboard.writeText(state.newick);
    status('Newick copied to clipboard', 'success');
  } catch (e) {
    status('Clipboard access denied', 'error');
  }
});
$('newick-download').addEventListener('click', () => download('phylo-fam-fate-tree.nwk', state.newick + '\n'));
$('newick-shown-download').addEventListener('click', () => {
  if (!state.shownTree) return;
  download('phylo-fam-fate-tree.rooted.nwk', toNewick(state.shownTree) + '\n');
});
$('svg-download').addEventListener('click', () => {
  const svg = $('tree-svg').cloneNode(true);
  svg.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
  download('phylo-fam-fate-tree.svg', new XMLSerializer().serializeToString(svg), 'image/svg+xml');
});

$('annot-download').addEventListener('click', () => {
  if (!state.tree) return;
  download('phylo-fam-fate-annotations.tsv',
    annotationsTsv(leafNames(state.tree), state.annotations, leafDomains), 'text/tab-separated-values');
});
$('itol-color-download').addEventListener('click', () => {
  if (!state.tree) return;
  const rank = $('tv-color').value || 'phylum';
  download(`phylo-fam-fate-itol-${rank}.txt`, itolColourStrip(leafNames(state.tree), state.annotations, rank));
});
$('itol-domain-download').addEventListener('click', () => {
  if (!state.tree) return;
  const txt = itolDomains(leafNames(state.tree), state.annotations, leafDomains);
  if (!txt) { status('No domain data yet: tick “Domains” (or annotate) first', 'error'); return; }
  download('phylo-fam-fate-itol-domains.txt', txt);
});
$('events-download').addEventListener('click', () => {
  if (!state.shownTree || !state.shownEvents) { status('Tick “Duplications” (on a rooted, annotated tree) first', 'error'); return; }
  download('phylo-fam-fate-events.tsv', eventsTsv(state.shownTree, state.annotations), 'text/tab-separated-values');
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
  download('phylo-fam-fate-metadata.json', JSON.stringify(data, null, 2), 'application/json');
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
        const ready = setInterval(() => {
          if (!$('tree-run').disabled) { clearInterval(ready); $('tree-run').click(); }
        }, 100);
      }
    }, 200);
  } else {
    const ready = setInterval(() => {
      if (!$('tree-run').disabled) { clearInterval(ready); $('tree-run').click(); }
    }, 100);
  }
}
