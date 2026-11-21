/* Remote sequence retrieval APIs — all CORS-enabled, no proxy needed. */

const NCBI_BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils';
const UNIPROT_BASE = 'https://rest.uniprot.org';
const INTERPRO_BASE = 'https://www.ebi.ac.uk/interpro/api';

const sleep = (ms) => new Promise(r => setTimeout(r, ms));

/** Parse a FASTA string into [{name, seq}]. */
export function parseFasta(text) {
  const records = [];
  let name = '', seq = '';
  for (const line of text.split(/\r?\n/)) {
    const t = line.trim();
    if (!t) continue;
    if (t.startsWith('>')) {
      if (name) records.push({ name, seq });
      name = t.slice(1).trim().split(/\s+/)[0];
      seq = '';
    } else {
      seq += t;
    }
  }
  if (name) records.push({ name, seq });
  return records;
}

export function toFasta(records) {
  return records.map(r => `>${r.name}\n${r.seq.match(/.{1,80}/g).join('\n')}`).join('\n') + '\n';
}

async function fetchJson(url) {
  const r = await fetch(url);
  if (!r.ok) throw new Error(`HTTP ${r.status} from ${new URL(url).host}`);
  return r.json();
}

/** UniProt text search -> FASTA records. */
export async function uniprotSearch(query, limit) {
  const url = `${UNIPROT_BASE}/uniprotkb/search?query=${encodeURIComponent(query)}&format=fasta&size=${limit}`;
  const r = await fetch(url);
  if (!r.ok) throw new Error(`UniProt HTTP ${r.status}`);
  return parseFasta(await r.text());
}

/** NCBI Entrez: esearch -> efetch -> FASTA records. */
export async function ncbiSearch(db, query, limit) {
  const searchUrl = `${NCBI_BASE}/esearch.fcgi?db=${db}&term=${encodeURIComponent(query)}&retmode=json&retmax=${limit}`;
  const data = await fetchJson(searchUrl);
  const ids = data.esearchresult?.idlist || [];
  if (!ids.length) return [];
  await sleep(350); // NCBI rate limit: 3 req/s without API key
  const fetchUrl = `${NCBI_BASE}/efetch.fcgi?db=${db}&id=${ids.join(',')}&rettype=fasta&retmode=text`;
  const r = await fetch(fetchUrl);
  if (!r.ok) throw new Error(`NCBI efetch HTTP ${r.status}`);
  return parseFasta(await r.text());
}

/** InterPro entry -> protein FASTA (via UniProt cross-reference). */
export async function interproSearch(accession, limit) {
  // Look up the entry to confirm and get its name
  const entryUrl = `${INTERPRO_BASE}/entry/${encodeURIComponent(accession)}`;
  const entry = await fetchJson(entryUrl);
  const name = entry?.data?.[0]?.name || accession;
  // UniProt can query by cross-reference: xref:interpro-PF00069
  const query = `xref:interpro-${accession}`;
  const records = await uniprotSearch(query, limit);
  return { records, name };
}

export async function retrieve(source, query, limit) {
  switch (source) {
    case 'uniprot': return uniprotSearch(query, limit);
    case 'ncbi-protein': return ncbiSearch('protein', query, limit);
    case 'ncbi-nucleotide': return ncbiSearch('nuccore', query, limit);
    case 'interpro': {
      const { records } = await interproSearch(query.trim(), limit);
      return records;
    }
    default: throw new Error(`Unknown API source: ${source}`);
  }
}
