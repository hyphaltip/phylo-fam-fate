/* Leaf annotations: species, taxonomy, protein/gene names, domains, function.
 *
 * Primary source is the UniProt REST batch endpoint (CORS-enabled), queried
 * once per tree with every accession found in the leaf labels. Labels that
 * aren't UniProt accessions fall back to whatever the FASTA description line
 * says (UniProt-style "OS=… OX=… GN=…" headers), which needs no network. */

const UNIPROT_FIELDS = [
  'accession', 'id', 'organism_name', 'organism_id', 'lineage', 'protein_name',
  'gene_names', 'length', 'xref_pfam', 'ft_domain', 'cc_function', 'annotation_score',
].join(',');
const BATCH = 400;                       // keeps the GET URL well under 8 kB

// UniProt accession format (https://www.uniprot.org/help/accession_numbers)
const ACC_RE = /^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:[.-]\d+)?$/;

export const RANKS = ['domain', 'kingdom', 'phylum', 'subphylum', 'class', 'order', 'family', 'genus', 'species'];

/** UniProt accession in a leaf label: "sp|P46151|MET12_YEAST", "P46151", "P46151.2". */
export function accessionOf(label) {
  const parts = String(label).split('|');
  const candidates = parts.length >= 2 && /^(sp|tr)$/i.test(parts[0]) ? [parts[1]] : parts;
  for (const c of candidates) {
    const m = c.trim().match(ACC_RE);
    if (m) return m[1];
  }
  return null;
}

/** Parse "cellular organisms (no rank), Eukaryota (domain), …" → [{name, rank}]. */
export function parseLineage(text) {
  if (!text) return [];
  return text.split(/,\s*(?![^()]*\))/).map(s => {
    const m = s.trim().match(/^(.*?)\s*\(([^()]+)\)$/);
    return m ? { name: m[1], rank: m[2] === 'superkingdom' ? 'domain' : m[2] } : { name: s.trim(), rank: 'no rank' };
  }).filter(x => x.name && x.name !== 'cellular organisms');
}

/** Parse ft_domain: 'DOMAIN 10..291; /note="…"; /evidence="…Pfam:PF02219"; DOMAIN …'. */
export function parseDomains(text) {
  const out = [];
  if (!text) return out;
  const re = /DOMAIN\s+[<>?]?(\d+)\.\.[<>?]?(\d+);\s*(?:\/note="([^"]*)")?(?:;\s*\/evidence="([^"]*)")?/g;
  let m;
  while ((m = re.exec(text))) {
    const src = (m[4] || '').match(/(Pfam|PROSITE|SMART|PROSITE-ProRule):([A-Z0-9_]+)/i);
    out.push({ start: +m[1], end: +m[2], name: m[3] || 'domain', db: src?.[1] || null, id: src?.[2] || null });
  }
  return out;
}

function parseTsv(text) {
  const lines = text.split('\n').filter(Boolean);
  if (!lines.length) return [];
  const head = lines[0].split('\t');
  return lines.slice(1).map(l => {
    const f = l.split('\t');
    return Object.fromEntries(head.map((h, i) => [h, f[i] ?? '']));
  });
}

/** Strain/serotype detail is noise on a tree: "Neurospora crassa (strain …)" → "Neurospora crassa". */
function shortSpecies(name) { return String(name || '').replace(/\s*\(.*$/, '').trim(); }

function fromUniprotRow(r) {
  const lineage = parseLineage(r['Taxonomic lineage']);
  const organism = r['Organism'];
  const species = lineage.find(x => x.rank === 'species')?.name || shortSpecies(organism);
  return {
    source: 'UniProt',
    accession: r['Entry'],
    entryName: r['Entry Name'],
    organism,
    species,
    taxid: r['Organism (ID)'] || null,
    lineage,
    protein: (r['Protein names'] || '').replace(/\s*\(EC [^)]*\)/g, '').replace(/\s*\(.*$/, '') || null,
    proteinFull: r['Protein names'] || null,
    genes: (r['Gene Names'] || '').split(/\s+/).filter(Boolean),
    length: parseInt(r['Length'], 10) || null,
    pfam: (r['Pfam'] || '').split(';').map(s => s.trim()).filter(Boolean),
    domains: parseDomains(r['Domain [FT]']),
    func: (r['Function [CC]'] || '').replace(/^FUNCTION:\s*/, '').replace(/\s*\{ECO:[^}]*\}/g, '').trim() || null,
    score: r['Annotation'] ? parseFloat(r['Annotation']) : null,
  };
}

/** Annotation from a UniProt-style FASTA description (no network). */
export function fromDescription(desc) {
  if (!desc) return null;
  const os = desc.match(/\bOS=(.+?)(?=\s+[A-Z]{2}=|$)/)?.[1];
  if (!os) return null;
  const words = desc.split(/\s+/);
  const protein = desc.split(/\s+OS=/)[0].split(/\s+/).slice(1).join(' ') || null;
  return {
    source: 'FASTA header',
    accession: accessionOf(words[0]),
    organism: os,
    species: shortSpecies(os),
    taxid: desc.match(/\bOX=(\d+)/)?.[1] || null,
    lineage: [],
    protein,
    genes: desc.match(/\bGN=(\S+)/) ? [desc.match(/\bGN=(\S+)/)[1]] : [],
    length: null, pfam: [], domains: [], func: null, score: null,
  };
}

/**
 * Annotate leaf labels.
 * @param {string[]} labels
 * @param {Map<string,string>} descs  label → FASTA description line
 * @param {(msg: string) => void} [onProgress]
 * @returns {Promise<Map<string, object>>} label → annotation
 */
export async function annotateLabels(labels, descs = new Map(), onProgress) {
  const out = new Map();
  const byAcc = new Map();
  for (const l of labels) {
    const acc = accessionOf(l);
    if (acc) (byAcc.get(acc) || byAcc.set(acc, []).get(acc)).push(l);
  }
  const accs = [...byAcc.keys()];
  for (let i = 0; i < accs.length; i += BATCH) {
    const chunk = accs.slice(i, i + BATCH);
    onProgress?.(`Fetching UniProt annotations ${i + 1}–${i + chunk.length} of ${accs.length}…`);
    const url = `https://rest.uniprot.org/uniprotkb/accessions?accessions=${chunk.join(',')}&fields=${UNIPROT_FIELDS}&format=tsv`;
    const r = await fetch(url);
    if (!r.ok) throw new Error(`UniProt returned HTTP ${r.status}`);
    for (const row of parseTsv(await r.text())) {
      const ann = fromUniprotRow(row);
      for (const l of byAcc.get(ann.accession) || []) out.set(l, ann);
    }
  }
  // Fill gaps (non-UniProt labels, obsolete accessions) from description lines.
  for (const l of labels) {
    if (out.has(l)) continue;
    const ann = fromDescription(descs.get(l));
    if (ann) out.set(l, ann);
  }
  return out;
}

/** Taxon of an annotation at a rank ('species' and 'genus' also work for header-only annotations). */
export function taxonAt(ann, rank) {
  if (!ann) return null;
  if (rank === 'species') return ann.species || null;
  const hit = ann.lineage.find(x => x.rank === rank);
  if (hit) return hit.name;
  if (rank === 'genus' && ann.species) return ann.species.split(' ')[0];
  return null;
}

/* Categorical palette (colour-blind-aware, distinguishable on white). */
export const PALETTE = ['#2563eb', '#dc2626', '#16a34a', '#9333ea', '#ea580c', '#0891b2',
  '#ca8a04', '#db2777', '#4d7c0f', '#7c2d12', '#1e40af', '#be185d'];
export const OTHER = '#94a3b8';

/** Assign colours to groups, most frequent first; the long tail shares grey. */
export function colourGroups(counts) {
  const sorted = [...counts.entries()].sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0]));
  const map = new Map();
  sorted.forEach(([g], i) => map.set(g, i < PALETTE.length ? PALETTE[i] : OTHER));
  return { map, sorted };
}

/** Stable colour per domain id/name. */
const DOMAIN_COLOURS = ['#6366f1', '#f59e0b', '#10b981', '#ef4444', '#8b5cf6', '#06b6d4', '#84cc16', '#f97316', '#ec4899', '#14b8a6'];
export function domainColour(key) {
  let h = 0;
  for (const ch of String(key)) h = (h * 31 + ch.charCodeAt(0)) >>> 0;
  return DOMAIN_COLOURS[h % DOMAIN_COLOURS.length];
}

/**
 * Pfam domain coordinates from InterPro, one request per protein (the API has
 * no batch form), run with limited concurrency and cached for the session.
 * Used for the domain track so every protein is drawn from the same source;
 * UniProt's own ft_domain features are the fallback.
 * @returns {Promise<Map<string, {start,end,name,id,db}[]>>} accession → domains
 */
const pfamCache = new Map();
export async function fetchPfamDomains(accessions, onProgress, concurrency = 3) {
  const todo = [...new Set(accessions)].filter(a => !pfamCache.has(a));
  let done = 0;
  const one = async acc => {
    try {
      // InterPro answers 503 ("no free Elasticsearch slot") or 408/429 under
      // load; those are transient, so back off and retry a few times.
      let r;
      for (let attempt = 0; ; attempt++) {
        r = await fetch(`https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/${acc}/?page_size=100`);
        if (![408, 429, 500, 502, 503, 504].includes(r.status) || attempt >= 4) break;
        await new Promise(res => setTimeout(res, 600 * 2 ** attempt + Math.random() * 400));
      }
      if (r.status === 204 || r.status === 404) { pfamCache.set(acc, []); return; }
      if (!r.ok) throw new Error(`HTTP ${r.status}`);
      const doms = [];
      for (const res of (await r.json()).results || []) {
        const m = res.metadata;
        for (const p of res.proteins || []) for (const loc of p.entry_protein_locations || []) for (const f of loc.fragments || []) {
          doms.push({ start: f.start, end: f.end, name: m.name, id: m.accession, db: 'Pfam' });
        }
      }
      doms.sort((a, b) => a.start - b.start);
      pfamCache.set(acc, doms);
    } catch {
      /* leave uncached so a later attempt can retry; caller falls back */
    } finally {
      onProgress?.(`Fetching Pfam domains from InterPro: ${++done}/${todo.length}…`);
    }
  };
  const queue = [...todo];
  await Promise.all(Array.from({ length: Math.min(concurrency, queue.length) }, async () => {
    while (queue.length) await one(queue.shift());
  }));
  return new Map(accessions.filter(a => pfamCache.has(a)).map(a => [a, pfamCache.get(a)]));
}
