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
  dbFasta: '',        // search database held in memory (custom / phmmer)
  dbUrl: '',          // streamed DB URL (BLAST presets)
  dbSource: '',       // human-readable DB source
  dbPreset: 'uniprot-sprot',
  proteomes: [],      // reference proteome catalogue
  searchHits: [],     // structured hit records
  selectedHits: new Set(),
  searchFasta: '',    // sequences from selected hits
  alignedFasta: '',   // aligned sequences
  trimmedFasta: '',   // trimmed alignment
  newick: '',
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
    $('timing-strip').innerHTML = state.timings.map(t =>
      `<span class="timing-item"><strong>${escapeHtml(t.step)}</strong> ${(t.ms / 1000).toFixed(2)}s</span>`).join('');
  }
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

/** Describe a URL-backed DB in the UI. */
function setDbUrl(url, label, note = '') {
  state.dbUrl = url;
  state.dbFasta = '';
  state.dbSource = label;
  $('search-db-info').innerHTML =
    `<strong>${escapeHtml(label)}</strong>${note ? ' · ' + note : ''} · streamed into Biowasm at search time (nothing uploaded)`;
  stepStatus('search', 'db ready', 'done');
}

/** Describe an in-memory DB in the UI. */
function setSearchDb({ text, records, gzipped, bytes }, source) {
  state.dbFasta = text;
  state.dbUrl = '';
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
    state.dbUrl = ''; state.dbFasta = '';
    $('search-db-info').textContent = 'Choose a reference proteome below, then click “Load proteome”.';
    stepStatus('search', 'choose proteome');
    loadProteomeCatalog();
  } else if (preset === 'pdb') {
    const proxy = $('search-proxy').value.trim();
    const url = proxy ? proxy + encodeURIComponent(PDB_URL) : PDB_URL;
    setDbUrl(url, 'PDB SEQRES',
      proxy ? 'via CORS proxy' : 'direct URL (CORS-blocked — see note)');
  } else {
    state.dbUrl = ''; state.dbFasta = '';
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

$('proteome-load').addEventListener('click', () => {
  const id = $('proteome-select').value;
  const p = state.proteomes.find(x => x.id === id);
  if (!p) { status('Select a reference proteome first', 'error'); return; }
  const division = p.division.charAt(0).toUpperCase() + p.division.slice(1);
  const url = `${UNIPROT_FTP}/reference_proteomes/${division}/${p.id}/${p.id}_${p.taxid}.fasta.gz`;
  setDbUrl(url, `${p.name} (${p.id})`, `${p.count.toLocaleString()} proteins`);
  status(`Selected ${p.name} — ready to search`, 'success');
});

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
  if ($('search-enable').checked && !state.dbUrl && !state.dbFasta) {
    status('Choose a database source for the local search', 'info');
  }
});

// --- Live hit filtering + selection --------------------------------------------
// Hit keys are stable numeric ids assigned after sorting; using an index (rather
// than a query\0name string) avoids embedding control characters in HTML attrs.
const hitKey = h => String(h.id);

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

function fmtE(v) { return v < 0.001 ? v.toExponential(1) : v.toFixed(4); }

function renderHits() {
  const maxHits = parseInt($('search-maxhits').value, 10) || 50;
  const eligible = state.searchHits.filter(passesCutoffs);
  const shown = eligible.slice(0, maxHits);
  const html = shown.map(h => {
    const key = hitKey(h);
    const checked = state.selectedHits.has(key) ? 'checked' : '';
    return `<tr class="${passesCutoffs(h) ? '' : 'filtered-out'}">
      <td><input type="checkbox" class="hit-check" data-key="${escapeAttr(key)}" ${checked}></td>
      <td title="${escapeAttr(h.name)}">${escapeHtml(trunc(h.name, 28))}</td>
      <td>${escapeHtml(h.query)}</td>
      <td>${fmtE(h.evalue)}</td>
      <td>${h.score.toFixed(1)}</td>
      <td>${h.pident == null ? '—' : h.pident.toFixed(1)}</td>
      <td>${h.length == null ? '—' : h.length}</td>
    </tr>`;
  }).join('');
  $('hit-table-body').innerHTML = html || '<tr><td colspan="7" class="hint">No hits pass the current cutoffs</td></tr>';
  document.querySelectorAll('.hit-check').forEach(cb => {
    cb.addEventListener('change', () => {
      if (cb.checked) state.selectedHits.add(cb.dataset.key);
      else state.selectedHits.delete(cb.dataset.key);
      updateHitSummary();
    });
  });
  $('hit-controls').classList.remove('hidden');
  $('hit-results').classList.remove('hidden');
  updateHitSummary();
}

function updateHitSummary() {
  const chosen = state.searchHits.filter(h => state.selectedHits.has(hitKey(h)) && passesCutoffs(h)).length;
  $('hit-summary').textContent = `${chosen} selected · ${state.searchHits.length} total hits`;
  $('hit-use-selected').textContent = `Use ${chosen} selected for alignment →`;
  $('hit-use-selected').disabled = chosen === 0;
}

['hit-evalue', 'hit-score', 'hit-pident', 'hit-length'].forEach(id => {
  $(id).addEventListener('input', renderHits);
});
$('hit-select-all').addEventListener('click', () => {
  state.searchHits.filter(passesCutoffs).forEach(h => state.selectedHits.add(hitKey(h)));
  renderHits();
});
$('hit-select-none').addEventListener('click', () => {
  state.selectedHits.clear(); renderHits();
});

$('hit-use-selected').addEventListener('click', () => {
  const chosen = state.searchHits.filter(h => state.selectedHits.has(hitKey(h)) && passesCutoffs(h) && h.seq);
  if (!chosen.length) { status('No usable sequences in the selection', 'error'); return; }
  state.searchFasta = toFasta(chosen.map(h => ({ name: h.name, seq: h.seq })));
  $('align-input').value = 'search';
  $('align-input').querySelector('option[value="search"]').textContent = `Search hits (${chosen.length})`;
  stepStatus('search', `${chosen.length} selected`, 'done');
  status(`${chosen.length} sequences ready for alignment`, 'success');
  recordMeta('Select hits', 'browser', VERSION, 'cutoffs + manual selection',
    `${state.searchHits.length} hits`, `${chosen.length} sequences`, null);
});

// --- Run the search ------------------------------------------------------------
$('search-run').addEventListener('click', async () => {
  const tool = $('search-tool').value;
  if (!state.rawFasta) { status('Load input sequences first', 'error'); return; }
  if (!state.dbUrl && !state.dbFasta) { status('Choose a database source first', 'error'); return; }

  const maxHits = parseInt($('search-maxhits').value, 10) || 50;
  stepStatus('search', 'running…', 'running');
  status(`Running ${tool} locally…`);
  $('search-run').disabled = true;
  const t0 = performance.now();
  try {
    let hits;
    if (tool === 'phmmer') {
      // phmmer needs the DB as a readable file; fetch/decode it if URL-backed.
      let dbFasta = state.dbFasta;
      if (!dbFasta && state.dbUrl) {
        const decoded = await fetchDb(state.dbUrl, state.dbSource);
        dbFasta = decoded.text;
        state.dbFasta = dbFasta; // cache for subsequent phmmer runs
      }
      const out = await runPhmmer(state.rawFasta, dbFasta);
      hits = parsePhmmerHits(out);
      // phmmer does not report sequences; recover them from the in-memory DB.
      const byName = new Map(parseFasta(dbFasta).map(r => [r.name, r.seq]));
      hits.forEach(h => { h.seq = byName.get(h.name) || ''; });
      $('search-preview').textContent = out.slice(0, 3000);
    } else {
      // BLAST embeds full subject sequences (sseq) in its tabular output.
      const out = await runBlast(tool, state.rawFasta, state.dbFasta, maxHits, state.dbUrl);
      hits = parseBlastHits(out);
      $('search-preview').textContent = out.slice(0, 3000);
    }
    $('search-preview').classList.remove('hidden');

    hits.sort((a, b) => a.evalue - b.evalue);
    hits.forEach((h, i) => { h.id = i; });
    state.searchHits = hits;
    state.selectedHits = new Set(hits.filter(passesCutoffs).slice(0, maxHits).map(hitKey));
    renderHits();
    stepStatus('search', `${hits.length} hits`, 'done');
    const ms = performance.now() - t0;
    recordMeta('Search', tool, 'biowasm',
      `maxHits=${maxHits}; evalue≤${$('hit-evalue').value}, score≥${$('hit-score').value}, pid≥${$('hit-pident').value}, len≥${$('hit-length').value}`,
      `${state.records.length} queries + ${state.dbSource}`, `${hits.length} hits`, ms);
    status(`Found ${hits.length} hits — filter and select which to align`, 'success');
  } catch (e) {
    status(`Search failed: ${e.message}`, 'error');
    stepStatus('search', 'error');
  } finally {
    $('search-run').disabled = false;
  }
});

// Initialize preset UI on load
applyDbPreset();

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
    $('trim-preview').textContent = fasta.slice(0, 2500) + (fasta.length > 2500 ? '\n…' : '');
    $('trim-preview').classList.remove('hidden');
    stepStatus('trim', `${stats.after} cols`, 'done');
    recordMeta('Trim', 'trimAl-like (JS)', VERSION,
      `gap≤${opts.gapThreshold}, cons≥${opts.minConservation}, ends=${opts.trimEnds}`,
      `${stats.before} cols`, `${stats.after} cols`, trimMs);
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
