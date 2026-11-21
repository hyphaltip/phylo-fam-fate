/* Lightweight alignment trimming in plain JS (trimAl-like heuristics). */
import { parseFasta, toFasta } from './api.js';

const GAP_CHARS = new Set(['-', '.', '?', '~']);
const GAP_CODES = new Uint8Array(128);
for (const c of GAP_CHARS) GAP_CODES[c.charCodeAt(0)] = 1;

/**
 * Trim an alignment.
 * @param {string} fasta - aligned FASTA
 * @param {object} opts - {gapThreshold, minConservation, trimEnds}
 * @returns {{fasta, stats}}
 */
export function trimAlignment(fasta, opts) {
  const records = parseFasta(fasta);
  if (!records.length) return { fasta, stats: { sequences: 0, before: 0, after: 0 } };

  const nSeq = records.length;
  const nPos = records[0].seq.length;
  if (records.some(r => r.seq.length !== nPos)) {
    throw new Error('Alignment is not rectangular — cannot trim');
  }

  // Column statistics: one pass over each sequence with a reusable count
  // table (much faster than per-column objects on large alignments).
  const seqs = records.map(r => r.seq.toUpperCase());
  const counts = new Uint32Array(128);
  const cols = new Array(nPos);
  for (let i = 0; i < nPos; i++) {
    counts.fill(0);
    let gaps = 0, maxCount = 0;
    for (let k = 0; k < nSeq; k++) {
      const code = seqs[k].charCodeAt(i) & 127;
      if (GAP_CODES[code]) { gaps++; continue; }
      const c = ++counts[code];
      if (c > maxCount) maxCount = c;
    }
    const gapFrac = gaps / nSeq;
    const conservation = nSeq > gaps ? maxCount / (nSeq - gaps) : 0;
    cols[i] = { i, gapFrac, conservation };
  }

  // Decide which columns to keep
  const keep = cols.filter(c =>
    c.gapFrac <= opts.gapThreshold && c.conservation >= opts.minConservation
  );

  let start = 0, end = keep.length;
  if (opts.trimEnds && keep.length) {
    // Trim leading/trailing columns where first/last sequence has gaps
    const firstSeq = records[0].seq;
    const lastSeq = records[records.length - 1].seq;
    while (start < keep.length && GAP_CHARS.has(firstSeq[keep[start].i])) start++;
    while (end > start && GAP_CHARS.has(lastSeq[keep[end - 1].i])) end--;
  }
  const selected = keep.slice(start, end);
  if (!selected.length) {
    throw new Error('These thresholds remove every column — raise the gap threshold or lower min conservation');
  }

  const trimmed = records.map(r => ({
    name: r.name,
    seq: selected.map(c => r.seq[c.i]).join(''),
  })).filter(r => r.seq.length > 0);

  return {
    fasta: toFasta(trimmed),
    stats: {
      sequences: nSeq,
      before: nPos,
      after: selected.length,
      removed: nPos - selected.length,
      kept: selected.length,
    },
  };
}
