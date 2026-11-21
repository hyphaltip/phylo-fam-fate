/* Export formats shared by the browser UI (app.js) and the command-line
 * tool (cli/famfate-js.mjs), so both write byte-identical files. */
import { leafNames } from './tree.js';
import { RANKS, taxonAt, colourGroups, domainColour, OTHER } from './annotate.js';

export const tsvCell = v => String(v ?? '').replace(/[\t\n\r]+/g, ' ');

/** Serialise a tree (with support values) back to Newick. */
export function toNewick(node) {
  const q = s => /[\s(),:;'\[\]]/.test(s) ? `'${s.replace(/'/g, "''")}'` : s;
  const walk = n => {
    const len = n.length ? `:${+n.length.toPrecision(6)}` : '';
    if (!n.children.length) return q(n.name) + len;
    return `(${n.children.map(walk).join(',')})${n.support != null ? n.support : ''}${len}`;
  };
  return walk(node) + ';';
}

/** Domains for a leaf: InterPro Pfam when present, else UniProt features. */
export function domainsFor(label, annotations, pfam) {
  const a = annotations.get(label);
  const pf = a?.accession ? pfam.get(a.accession) : null;
  return pf?.length ? pf : a?.domains || [];
}

/** One row per leaf: taxonomy, names, Pfam, domain coordinates, function. */
export function annotationsTsv(labels, annotations, domainsOf) {
  const rankCols = RANKS.filter(r => r !== 'species');
  const head = ['label', 'accession', 'species', 'taxid', ...rankCols, 'protein', 'genes', 'length', 'pfam', 'domains', 'function', 'annotation_score', 'source'];
  const rows = labels.map(l => {
    const a = annotations.get(l);
    const doms = domainsOf(l).map(d => `${d.id || d.name}:${d.start}-${d.end}`).join(';');
    return [l, a?.accession, a?.species, a?.taxid, ...rankCols.map(r => taxonAt(a, r)), a?.protein,
      a?.genes?.join(' '), a?.length, a?.pfam?.join(';'), doms, a?.func, a?.score, a?.source].map(tsvCell).join('\t');
  });
  return [head.join('\t'), ...rows].join('\n') + '\n';
}

/* iTOL dataset files (https://itol.embl.de/help.cgi#datasets). Leaf IDs must
   match the Newick labels, so these pair with either Newick export. */
export function itolColourStrip(labels, annotations, rank) {
  const counts = new Map();
  for (const l of labels) { const g = taxonAt(annotations.get(l), rank) || 'unknown'; counts.set(g, (counts.get(g) || 0) + 1); }
  const { map, sorted } = colourGroups(counts);
  const col = g => g === 'unknown' ? '#e2e8f0' : map.get(g) || OTHER;
  const lines = ['DATASET_COLORSTRIP', 'SEPARATOR TAB', `DATASET_LABEL\t${rank}`, 'COLOR\t#2563eb',
    `LEGEND_TITLE\t${rank}`, `LEGEND_SHAPES\t${sorted.map(() => 1).join('\t')}`,
    `LEGEND_COLORS\t${sorted.map(([g]) => col(g)).join('\t')}`, `LEGEND_LABELS\t${sorted.map(([g]) => tsvCell(g)).join('\t')}`,
    'DATA', ...labels.map(l => { const g = taxonAt(annotations.get(l), rank) || 'unknown'; return `${l}\t${col(g)}\t${tsvCell(g)}`; })];
  return lines.join('\n') + '\n';
}

/** @returns {string|null} null when no leaf has a protein length/domains */
export function itolDomains(labels, annotations, domainsOf) {
  const lines = ['DATASET_DOMAINS', 'SEPARATOR TAB', 'DATASET_LABEL\tPfam domains', 'COLOR\t#6366f1', 'DATA'];
  let n = 0;
  for (const l of labels) {
    const doms = domainsOf(l);
    const len = annotations.get(l)?.length || (doms.length ? Math.max(...doms.map(d => d.end)) : 0);
    if (!len) continue;
    n++;
    // SHAPE|START|END|COLOR|LABEL ; RE = rounded rectangle
    lines.push([l, len, ...doms.map(d => `RE|${d.start}|${d.end}|${domainColour(d.id || d.name)}|${tsvCell(d.name).replace(/\|/g, '/')}`)].join('\t'));
  }
  return n ? lines.join('\n') + '\n' : null;
}

/** Duplication/speciation nodes of a tree already labelled by labelEvents(). */
export function eventsTsv(tree, annotations) {
  const rows = ['node\tevent\tsupport\ttips\tspecies\ttip_names'];
  let i = 0;
  const walk = n => {
    if (!n.children.length) return;
    if (n.event) {
      const names = leafNames(n);
      const sp = new Set(names.map(x => annotations.get(x)?.species).filter(Boolean));
      rows.push([`n${++i}`, n.event === 'dup' ? 'duplication' : 'speciation', n.support ?? '', names.length, [...sp].join('; '), names.join(',')].map(tsvCell).join('\t'));
    }
    n.children.forEach(walk);
  };
  walk(tree);
  return rows.join('\n') + '\n';
}
