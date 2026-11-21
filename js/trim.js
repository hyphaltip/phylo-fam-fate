/* Lightweight alignment trimming in plain JS (trimAl-like heuristics). */
import { parseFasta, toFasta } from './api.js';

const GAP_CHARS = new Set(['-', '.', '?', '~']);

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

  // Column statistics
  const cols = [];
  for (let i = 0; i < nPos; i++) {
    let gaps = 0;
    const counts = {};
    for (const r of records) {
      const c = r.seq[i].toUpperCase();
      if (GAP_CHARS.has(c)) gaps++;
      else counts[c] = (counts[c] || 0) + 1;
    }
    const gapFrac = gaps / nSeq;
    const maxCount = Math.max(0, ...Object.values(counts));
    const conservation = nSeq > gaps ? maxCount / (nSeq - gaps) : 0;
    cols.push({ i, gapFrac, conservation });
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
