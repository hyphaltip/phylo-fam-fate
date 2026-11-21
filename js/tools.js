/* Biowasm/Aioli tool wrappers. Aioli is loaded from CDN in index.html. */

let cli = null;
const requestedTools = new Set();
// BLAST databases already built in the current Aioli filesystem, keyed by
// source URL + type. Cleared whenever Aioli is re-created (new filesystem).
const builtDbs = new Map();

function parseSpec(spec) {
  const p = spec.split('/');
  return p.length === 3
    ? { tool: p[0], version: p[2], program: p[1] }
    : { tool: p[0], version: p[1] };
}

/** Normalize Aioli exec() output (plain string, or {stdout, stderr} object). */
function execText(result) {
  if (typeof result === 'string') return result;
  return result?.stdout ?? '';
}

export async function initTools(tools) {
  if (typeof Aioli === 'undefined') {
    throw new Error('Aioli failed to load from CDN — check your network connection');
  }
  const before = requestedTools.size;
  tools.forEach(t => requestedTools.add(t));
  if (cli && requestedTools.size === before) return cli;
  // Aioli's filesystem is shared per instance; re-init with the full tool set
  // whenever a new tool is added so earlier tools remain available.
  // Use Aioli's default output mode: exec() returns a plain string. Some
  // programs (notably BLAST) do not populate stdout when printInterleaved=false.
  cli = await new Aioli([...requestedTools].map(parseSpec));
  builtDbs.clear();
  return cli;
}

/** Exact Biowasm tool versions, recorded in run metadata for provenance. */
export const TOOL_VERSIONS = {
  muscle: '5.1.0', kalign: '3.3.1', blastp: '2.17.0', blastn: '2.17.0', blastdbcmd: '2.17.0', phmmer: '3.4',
};

/** Mount a FASTA string as a virtual file, run an aligner, return aligned FASTA. */
export async function runAligner(tool, fasta) {
  await initTools([`${tool}/${TOOL_VERSIONS[tool]}`]);
  // Unique file names per run: a stale aligned.fa from an earlier run would
  // otherwise be returned silently if this run fails to write its output.
  const tag = `${Date.now().toString(36)}${Math.random().toString(36).slice(2, 6)}`;
  const inputName = `input_${tag}.fa`;
  const outputName = `aligned_${tag}.fa`;
  await cli.mount([{ name: inputName, data: fasta }]);

  let cmd;
  switch (tool) {
    case 'muscle':
      cmd = `muscle -align ${inputName} -output ${outputName}`;
      break;
    case 'kalign':
      cmd = `kalign -i ${inputName} -o ${outputName}`;
      break;
    default:
      throw new Error(`Unknown aligner: ${tool}`);
  }

  const result = await cli.exec(cmd);
  // Biowasm tools run with cwd=/shared/data; output files are relative to it.
  let out = '';
  try { out = await cli.cat(outputName); } catch { /* no output file */ }
  if (typeof out !== 'string' || !out.trim()) {
    // Keep the whole log: the line that says *why* (e.g. "Only 1 sequence
    // was found") comes before the stack of "Function ... failed" lines.
    const log = execText(result).trim();
    throw new Error(`${tool} produced no alignment.\n${log.slice(-2000)}`);
  }
  return out;
}

/** Short stable id for a string (file names in the WASM filesystem). */
function shortHash(str) {
  let h = 2166136261;
  for (let i = 0; i < str.length; i++) { h ^= str.charCodeAt(i); h = Math.imul(h, 16777619); }
  return (h >>> 0).toString(36);
}

/** Mount one FASTA source and build a v4 BLAST database from it (cached). */
async function ensureBlastDb(source, dbType, onProgress) {
  const key = `${dbType}\u0000${source.url || 'inline:' + shortHash(source.data)}`;
  if (builtDbs.has(key)) return builtDbs.get(key);
  const id = shortHash(key);
  const file = source.url ? `db_${id}.fa.gz` : `db_${id}.fa`;
  const name = `bdb_${id}`;
  onProgress?.(`Indexing ${source.label || file}…`);
  await cli.mount([source.url ? { name: file, url: source.url } : { name: file, data: source.data }]);
  // v5 databases use LMDB, which is unavailable in the WASM filesystem;
  // v4 uses traditional .pin/.psq index files and works in-browser.
  // -parse_seqids indexes the IDs so blastdbcmd can return full-length
  // subjects later. It rejects some custom headers (e.g. duplicate IDs), so
  // fall back to an unindexed build rather than failing the search.
  const build = async parse => {
    const log = execText(await cli.exec(
      `makeblastdb -in ${file} -dbtype ${dbType} -out ${name} -blastdb_version 4${parse ? ' -parse_seqids' : ''}`));
    const added = log.match(/added\s+(\d+)\s+sequences/i);
    // The WASM build always prints a harmless
    // "Error: (106.19) ... SetFullVersion() should be used ..." line, so a
    // reported sequence count is the real success signal. Only fail when
    // there is no count *and* a genuine error remains after dropping that line.
    const realErrors = log.split('\n')
      .filter(l => /BLAST options error|Error:/i.test(l) && !/SetFullVersion/.test(l));
    return { log, added, ok: !!added || !realErrors.length };
  };
  let parsed = true;
  let { log, added, ok } = await build(true);
  if (!ok) {
    onProgress?.(`${source.label || file}: IDs could not be indexed; full-length sequences unavailable for it`);
    parsed = false;
    ({ log, added, ok } = await build(false));
  }
  if (!ok) {
    throw new Error(`makeblastdb failed for ${source.label || file}:\n${log.trim().slice(-1500)}`);
  }
  const info = { name, parsed, sequences: added ? parseInt(added[1], 10) : null };
  builtDbs.set(key, info);
  return info;
}

/**
 * Run local BLAST against one or more FASTA sources.
 * @param {{url?: string, data?: string, label?: string}[]} sources - each URL
 *   is streamed into the WASM filesystem (never fully held in JS memory) and
 *   indexed separately; blastp/blastn then searches them together.
 */
export async function runBlast(program, queryFasta, sources, maxHits, onProgress, evalue = 1e-3) {
  // Biowasm exposes each BLAST binary as a separate program (there is no
  // monolithic `blast.js`), so request the chosen program plus makeblastdb.
  await initTools([`blast/${program}/2.17.0`, 'blast/makeblastdb/2.17.0', 'blast/blastdbcmd/2.17.0']);
  const isNucleotide = program === 'blastn';
  const dbType = isNucleotide ? 'nucl' : 'prot';
  const dbs = [];
  for (const [i, src] of sources.entries()) {
    onProgress?.(`Preparing database ${i + 1}/${sources.length}: ${src.label || 'custom'}…`);
    dbs.push(await ensureBlastDb(src, dbType, onProgress));
  }
  await cli.mount([{ name: 'query.fa', data: queryFasta }]);
  onProgress?.(`Running ${program} against ${dbs.length} database${dbs.length > 1 ? 's' : ''}…`);
  // Disable low-complexity masking so tandem-repeat queries aren't filtered
  // out entirely. The search E-value removes noise up front; the UI's live
  // cutoffs can tighten it further.
  const filter = isNucleotide ? ['-task', 'blastn', '-dust', 'no'] : ['-seg', 'no'];
  // Array form: each flag is one argv element. A multi-database -db is a
  // single space-separated argument, which BLAST treats as one combined set
  // (E-values are computed against the total size of all of them).
  const result = await cli.exec(program, [
    '-query', 'query.fa', '-db', dbs.map(d => d.name).join(' '), ...filter,
    '-evalue', String(evalue), '-max_target_seqs', String(maxHits),
    '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq sacc slen',
  ]);
  const text = execText(result);
  if (!text.trim()) throw new Error(`${program} found no hits at E-value ≤ ${evalue} — try a looser search E-value or a different database`);
  const known = dbs.every(d => d.sequences != null);
  onProgress?.('Retrieving full-length subject sequences…');
  const full = await fetchFullSubjects(text, dbs, sources);
  return { text, full, sequences: known ? dbs.reduce((n, d) => n + d.sequences, 0) : null };
}

/** Full-length sequences for every subject in a BLAST tabular result, keyed
 *  by accession (and sseqid). Uses blastdbcmd on the ID-indexed databases;
 *  inline sources that could not be indexed are looked up from their text. */
async function fetchFullSubjects(tabular, dbs, sources) {
  const full = new Map();
  const ids = new Set();
  for (const line of tabular.split('\n')) {
    const f = line.split('\t');
    if (f.length >= 14 && f[13]) ids.add(f[13]);
  }
  if (!ids.size) return full;
  const parsedDbs = dbs.filter(d => d.parsed).map(d => d.name);
  if (parsedDbs.length) {
    await cli.mount([{ name: 'subject_ids.txt', data: [...ids].join('\n') + '\n' }]);
    // %a = accession, %s = sequence. IDs absent from one of the databases
    // are reported on stderr and skipped, which is expected with several.
    const out = execText(await cli.exec('blastdbcmd', [
      '-db', parsedDbs.join(' '), '-entry_batch', 'subject_ids.txt', '-outfmt', '%a\t%s']));
    for (const line of out.split('\n')) {
      const tab = line.indexOf('\t');
      if (tab < 1) continue;
      const acc = line.slice(0, tab).trim();
      const seq = line.slice(tab + 1).trim();
      if (!/^[A-Za-z*-]+$/.test(seq)) continue;
      // Key by accession both with and without a version suffix (X.1 / X).
      for (const k of [acc, acc.replace(/\.\d+$/, '')]) if (!full.has(k)) full.set(k, seq);
    }
  }
  // Inline (uploaded/pasted) sources that failed ID indexing: scan their text.
  for (const [i, d] of dbs.entries()) {
    if (d.parsed || !sources[i]?.data) continue;
    let name = null, seq = [];
    const flush = () => { if (name && !full.has(name)) full.set(name, seq.join('')); };
    for (const line of sources[i].data.split(/\r?\n/)) {
      if (line.startsWith('>')) { flush(); name = line.slice(1).trim().split(/\s+/)[0]; seq = []; }
      else if (name) seq.push(line.trim());
    }
    flush();
  }
  return full;
}

/** Run local phmmer against a user-provided FASTA database.
 *  The WASM build is `phmmer-simd`; it supports -E/-T but not --maxhits/-M. */
export async function runPhmmer(queryFasta, dbFasta, evalue = 1e-3) {
  await initTools(['hmmer/phmmer/3.4']);
  await cli.mount([
    { name: 'query.fa', data: queryFasta },
    { name: 'db.fa', data: dbFasta },
  ]);
  // --noali keeps output compact; the score table and domain coordinates we
  // parse are still present. Truncation happens client-side.
  const result = await cli.exec(`phmmer -E ${evalue} --domE ${evalue} --noali query.fa db.fa`);
  const text = execText(result);
  if (!text.trim()) throw new Error(`phmmer found no hits at E-value ≤ ${evalue} — try a looser search E-value or a different database`);
  return text;
}

/** Parse BLAST outfmt 6 -> array of structured hit records. */
export function parseBlastHits(tabular) {
  const hits = [];
  for (const line of tabular.split('\n')) {
    if (!line.trim() || line.startsWith('#')) continue;
    const f = line.split('\t');
    if (f.length < 12) continue;
    hits.push({
      name: f[1],
      query: f[0],
      pident: parseFloat(f[2]),
      length: parseInt(f[3], 10),
      evalue: parseFloat(f[10]),
      score: parseFloat(f[11]),
      // sseq is the *aligned region* of the subject (HSP only, with '-' for
      // gaps), not the full database sequence. Strip gaps so the aligner
      // gets raw residues.
      seq: (f[12] || '').replace(/[-.]/g, ''),
      sstart: parseInt(f[8], 10),
      send: parseInt(f[9], 10),
      acc: f[13] || null,
      slen: f[14] ? parseInt(f[14], 10) : null,
    });
  }
  return hits;
}

/** Parse phmmer text output -> array of structured hit records.
 *  Reads the "Scores for complete sequences" table and the per-domain
 *  alignment coordinates (for alignment length). Identity is not reported by
 *  phmmer's tabular/text summary, so pident is null for these hits. */
export function parsePhmmerHits(output) {
  const hits = [];
  let query = '';
  let currentTarget = '';
  const lengths = new Map(); // `${query}\0${name}` -> aligned length

  for (const line of output.split('\n')) {
    // Query header: "Query:       query  [L=151]"
    let m = line.match(/^\s*Query:\s+(\S+)/);
    if (m) { query = m[1]; continue; }

    // Domain annotation rows, inside ">> target  description" blocks.
    // e.g. "   1 !  329.1  0.6  2.4e-103  2.4e-103  1  151 []  1  151 []  1  151 [] 1.00"
    if (/^\s*\d+\s+[!?.]\s/.test(line)) {
      // HMMER writes domain-end markers as [] or [. — strip brackets and the
      // stray '.' tokens so column indices stay stable.
      const parts = line.replace(/[\[\]]/g, ' ').trim().split(/\s+/).filter(t => t !== '.');
      // parts: #, mark, score, bias, c-E, i-E, hmmfrom, hmmto, alifrom, alito, ...
      const from = parseInt(parts[8], 10);
      const to = parseInt(parts[9], 10);
      if (!isNaN(from) && !isNaN(to) && to >= from && currentTarget) {
        const key = `${query}\u0000${currentTarget}`;
        lengths.set(key, (lengths.get(key) || 0) + (to - from + 1));
      }
      continue;
    }
    m = line.match(/^>>\s+(\S+)/);
    if (m) { currentTarget = m[1]; continue; }

    // Full-sequence score table row (>= 9 fields):
    // fullE fullScore bias domE domScore domBias exp N name [description...]
    m = line.match(/^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\S+)/);
    if (!m) continue;
    const evalue = parseFloat(m[1]);
    const score = parseFloat(m[2]);
    const name = m[9];
    if (isNaN(evalue) || !name || name === 'Sequence' || !/^[A-Za-z0-9]/.test(name)) continue;
    hits.push({ name, query: query || 'query', evalue, score, pident: null, length: null });
  }

  // Deduplicate by (query, name), keeping the best E-value; attach lengths.
  const best = new Map();
  for (const h of hits) {
    h.length = lengths.get(`${h.query}\u0000${h.name}`) || null;
    const k = `${h.query}\u0000${h.name}`;
    if (!best.has(k) || h.evalue < best.get(k).evalue) best.set(k, h);
  }
  return [...best.values()];
}
