/* Compact, colored alignment preview (Clustal-like residue classes). */
import { parseFasta } from './api.js';

const AA_CLASS = {};
for (const c of 'AVLIMFWC') AA_CLASS[c] = 'hyd';
for (const c of 'KR') AA_CLASS[c] = 'pos';
for (const c of 'DE') AA_CLASS[c] = 'neg';
for (const c of 'STNQ') AA_CLASS[c] = 'pol';
for (const c of 'G') AA_CLASS[c] = 'gly';
for (const c of 'P') AA_CLASS[c] = 'pro';
for (const c of 'HY') AA_CLASS[c] = 'aro';
const NT_CLASS = { A: 'nA', C: 'nC', G: 'nG', T: 'nT', U: 'nT' };

const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

/**
 * Render an alignment preview into `el`.
 * @param {HTMLElement} el
 * @param {string} fasta
 * @param {{type?: 'nt'|'aa', maxSeqs?: number, maxCols?: number}} opts
 */
export function renderMsa(el, fasta, { type = 'aa', maxSeqs = 40, maxCols = 160 } = {}) {
  const recs = parseFasta(fasta);
  if (!recs.length) { el.innerHTML = ''; return; }
  const nCols = Math.max(...recs.map(r => r.seq.length));
  const shown = recs.slice(0, maxSeqs);
  const cols = Math.min(nCols, maxCols);
  const cls = type === 'nt' ? NT_CLASS : AA_CLASS;

  const rows = shown.map(r => {
    let html = '';
    let run = '', runCls = null;
    const flush = () => {
      if (!run) return;
      html += runCls ? `<span class="r-${runCls}">${run}</span>` : run;
      run = '';
    };
    for (let i = 0; i < cols; i++) {
      const ch = r.seq[i] || ' ';
      const k = cls[ch.toUpperCase()] || null;
      if (k !== runCls) { flush(); runCls = k; }
      run += esc(ch);
    }
    flush();
    return `<div class="msa-row"><span class="msa-name" title="${esc(r.name)}">${esc(r.name)}</span><span class="msa-seq">${html}</span></div>`;
  }).join('');

  const ruler = Array.from({ length: cols }, (_, i) => ((i + 1) % 10 === 0 ? '|' : ' ')).join('');
  const note = [
    `${recs.length} sequences × ${nCols} columns`,
    recs.length > shown.length ? `showing first ${shown.length} sequences` : '',
    nCols > cols ? `columns 1–${cols}` : '',
  ].filter(Boolean).join(' · ');

  el.innerHTML = `<div class="msa-note">${note}</div>
    <div class="msa-scroll"><div class="msa">
      <div class="msa-row msa-ruler"><span class="msa-name"></span><span class="msa-seq">${ruler}</span></div>
      ${rows}
    </div></div>`;
}
