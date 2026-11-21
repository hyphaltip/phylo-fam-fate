#!/usr/bin/env node
/* famfate-js: runs the website's own JavaScript (js/*.js) from the command
 * line, so CLI results match the browser exactly for the steps that the
 * site implements itself: trimming, rooting (midpoint / MAD / outgroup),
 * UniProt + InterPro annotation, duplication labelling, SVG and iTOL export.
 * It also runs the browser's FastTree WebAssembly build, for parity checks.
 *
 *   famfate-js trim ALN.fa [--gap 0.5] [--cons 0] [--ends]        > trimmed.fa
 *   famfate-js seqtype FASTA                                       → nt | aa
 *   famfate-js post --tree T.nwk [--fasta FAMILY.fa] [--out PREFIX]
 *                   [--root midpoint|mad|none|outgroup:RANK:NAME]
 *                   [--rank phylum] [--label id|species|id-species|gene|protein]
 *                   [--no-annotate] [--no-domains] [--width 1100]
 *   famfate-js fasttree-wasm ALN.fa [FastTree args…]               > tree.nwk
 *   famfate-js compare-trees A.nwk B.nwk [--tol 1e-4]              (exit 1 if different)
 */
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const HERE = path.dirname(fileURLToPath(import.meta.url));
// Inside the container the site lives at /opt/phylo-fam-fate; in a checkout it is ../
const SITE = process.env.FAMFATE_SITE || path.resolve(HERE, '..');
const mod = p => import(path.join(SITE, p));

const { parseFasta } = await mod('js/api.js');
const T = await mod('js/tree.js');

function die(msg, code = 2) { console.error(`famfate-js: ${msg}`); process.exit(code); }
function parseArgs(argv) {
  const pos = [], opt = {};
  for (let i = 0; i < argv.length; i++) {
    const a = argv[i];
    if (a.startsWith('--')) {
      const [k, v] = a.slice(2).split('=', 2);
      if (v !== undefined) opt[k] = v;
      else if (argv[i + 1] !== undefined && !argv[i + 1].startsWith('--')) opt[k] = argv[++i];
      else opt[k] = true;
    } else pos.push(a);
  }
  return { pos, opt };
}
const read = f => { try { return fs.readFileSync(f, 'utf8'); } catch (e) { die(`cannot read ${f}: ${e.message}`); } };

/** Same rule as the browser (app.js seqType). */
function seqType(fasta) {
  const residues = fasta.split(/\r?\n/).filter(l => l.trim() && !l.trim().startsWith('>')).join('').toUpperCase();
  if (!residues) return 'nt';
  if (/[EFILPQZXJ]/.test(residues)) return 'aa';
  const acgt = (residues.match(/[ACGTU]/g) || []).length;
  return acgt / residues.length < 0.8 ? 'aa' : 'nt';
}

const [cmd, ...rest] = process.argv.slice(2);
const { pos, opt } = parseArgs(rest);

switch (cmd) {
  case 'trim': {
    const { trimAlignment } = await mod('js/trim.js');
    if (!pos[0]) die('usage: trim ALN.fa [--gap 0.5] [--cons 0] [--ends]');
    const { fasta, stats } = trimAlignment(read(pos[0]), {
      gapThreshold: parseFloat(opt.gap ?? 0.5), minConservation: parseFloat(opt.cons ?? 0), trimEnds: !!opt.ends,
    });
    process.stdout.write(fasta);
    console.error(`trim: ${stats.sequences} sequences, ${stats.before} → ${stats.after} columns`);
    break;
  }

  case 'seqtype':
    if (!pos[0]) die('usage: seqtype FASTA');
    console.log(seqType(read(pos[0])));
    break;

  case 'fasttree-wasm': {
    if (!pos[0]) die('usage: fasttree-wasm ALN.fa [FastTree args…]');
    const ft = await mod('assets/fasttree/phylo_fasttree_wasm.js');
    ft.initSync(fs.readFileSync(path.join(SITE, 'assets/fasttree/phylo_fasttree_wasm_bg.wasm')));
    const args = process.argv.slice(4);       // pass FastTree flags through verbatim
    process.stdout.write(ft.buildTree(read(pos[0]), args).trim() + '\n');
    break;
  }

  case 'compare-trees': {
    if (pos.length < 2) die('usage: compare-trees A.nwk B.nwk [--tol 1e-4]');
    const tol = parseFloat(opt.tol ?? 1e-4);
    const [a, b] = pos.map(f => T.parseNewick(read(f).trim()));
    const res = compareTrees(a, b, tol);
    console.log(JSON.stringify(res, null, 2));
    process.exit(res.sameTopology && res.maxLengthDiff <= tol ? 0 : 1);
  }

  case 'post': await post(); break;

  default:
    console.error(fs.readFileSync(fileURLToPath(import.meta.url), 'utf8').split('*/')[0].replace(/^#!.*\n\/\*/, '').replace(/^ \* ?/gm, ''));
    process.exit(cmd ? 2 : 0);
}

/** Unrooted comparison: same leaf set, same bipartitions, branch lengths per split. */
function compareTrees(a, b, tol) {
  const splits = root => {
    const leaves = T.leafNames(root).sort();
    const idx = new Map(leaves.map((l, i) => [l, i]));
    const out = new Map();
    const walk = n => {
      if (!n.children.length) {
        const key = canon(new Set([idx.get(n.name)]), leaves.length);
        out.set(key, (out.get(key) || 0) + n.length);
        return [idx.get(n.name)];
      }
      const below = n.children.flatMap(walk);
      if (n !== root) {
        const key = canon(new Set(below), leaves.length);
        // The two edges at a bifurcating root are one unrooted edge: sum them.
        out.set(key, (out.get(key) || 0) + n.length);
      }
      return below;
    };
    walk(root);
    return { leaves, splits: out };
  };
  // Canonical split key: the side not containing leaf 0.
  const canon = (set, n) => {
    const side = set.has(0) ? [...Array(n).keys()].filter(i => !set.has(i)) : [...set];
    return side.sort((x, y) => x - y).join(',');
  };
  const A = splits(a), B = splits(b);
  const sameLeaves = A.leaves.join('\n') === B.leaves.join('\n');
  let maxLengthDiff = 0, missing = 0;
  for (const [k, len] of A.splits) {
    if (!k) continue;
    if (!B.splits.has(k)) { missing++; continue; }
    maxLengthDiff = Math.max(maxLengthDiff, Math.abs(len - B.splits.get(k)));
  }
  const extra = [...B.splits.keys()].filter(k => k && !A.splits.has(k)).length;
  return { tips: A.leaves.length, sameLeaves, sameTopology: sameLeaves && !missing && !extra,
    splitsOnlyInA: missing, splitsOnlyInB: extra, maxLengthDiff };
}

async function post() {
  const A = await mod('js/annotate.js');
  const E = await mod('js/export.js');
  if (!opt.tree) die('post needs --tree T.nwk');
  const tree = T.parseNewick(read(opt.tree).trim());
  const labels = T.leafNames(tree);
  const prefix = opt.out || opt.tree.replace(/\.(nwk|newick|tree|tre)$/i, '');
  const rank = opt.rank || 'phylum';
  const log = m => console.error(`post: ${m}`);
  const descs = new Map();
  if (opt.fasta) for (const r of parseFasta(read(opt.fasta))) if (r.desc) descs.set(r.name, r.desc);

  // 1. Annotation (UniProt batch; FASTA OS= fallback)
  let annotations = new Map();
  if (!opt['no-annotate']) {
    try {
      annotations = await A.annotateLabels(labels, descs, log);
      log(`annotated ${annotations.size}/${labels.length} leaves`);
    } catch (e) {
      log(`annotation failed (${e.message}); continuing without`);
    }
  }
  const pfam = new Map();
  if (!opt['no-annotate'] && !opt['no-domains']) {
    const accs = labels.map(l => annotations.get(l)?.accession || A.accessionOf(l)).filter(Boolean);
    if (accs.length) {
      let done = 0;
      const got = await A.fetchPfamDomains(accs, () => { if (++done % 25 === 0) log(`InterPro ${done}/${accs.length}`); });
      for (const [k, v] of got) pfam.set(k, v);
      log(`Pfam domains for ${got.size}/${accs.length} proteins`);
    }
  }

  // 2. Rooting — same functions as the browser's Root menu
  const how = String(opt.root || 'midpoint');
  let rooted, note;
  if (how === 'none') { rooted = structuredClone(tree); note = 'unrooted (as inferred)'; }
  else if (how === 'mad') {
    const r = T.madRoot(tree);
    if (r) { rooted = r.tree; note = `MAD root (deviation ${r.score.toFixed(4)}, ambiguity ${r.ambiguity.toFixed(3)})`; }
    else { rooted = T.midpointRoot(structuredClone(tree)); note = 'MAD needs 3–900 tips; midpoint used'; }
  } else if (how.startsWith('outgroup:')) {
    const [, r, ...nameParts] = how.split(':');
    const name = nameParts.join(':');
    const set = new Set(labels.filter(l => A.taxonAt(annotations.get(l), r) === name));
    const res = T.rootOnOutgroup(structuredClone(tree), set);
    if (res.invalid) { rooted = T.midpointRoot(structuredClone(tree)); note = `outgroup ${name} (${r}) matched ${set.size} of ${labels.length} tips; midpoint used`; }
    else { rooted = res.tree; note = res.exact ? `rooted on ${name} (${set.size} tips, monophyletic)` : `rooted on best split for non-monophyletic ${name}: ${res.impurity} tips on the wrong side`; }
  } else { rooted = T.midpointRoot(structuredClone(tree)); note = 'midpoint root'; }
  T.ladderize(rooted);
  log(note);

  // 3. Duplications (species overlap), on the rooted tree
  let events = null;
  if (how !== 'none' && annotations.size) {
    events = T.labelEvents(rooted, leaf => annotations.get(leaf.name)?.taxid || annotations.get(leaf.name)?.species || null);
    log(`${events.dup} duplications, ${events.spec} speciations`);
  }

  // 4. Outputs
  const domainsOf = l => E.domainsFor(l, annotations, pfam);
  const written = [];
  const out = (suffix, text) => { if (text == null) return; fs.writeFileSync(prefix + suffix, text); written.push(prefix + suffix); };
  out('.rooted.nwk', E.toNewick(rooted) + '\n');
  if (annotations.size) {
    out('.annotations.tsv', E.annotationsTsv(labels, annotations, domainsOf));
    out(`.itol_${rank}.txt`, E.itolColourStrip(labels, annotations, rank));
    out('.itol_domains.txt', E.itolDomains(labels, annotations, domainsOf));
  }
  if (events) out('.events.tsv', E.eventsTsv(rooted, annotations));

  // SVG drawn by the browser's renderer, with a minimal stand-in for <svg>
  const svg = { attrs: {}, setAttribute(k, v) { this.attrs[k] = v; }, innerHTML: '', parentElement: null };
  const counts = new Map();
  for (const l of labels) { const g = A.taxonAt(annotations.get(l), rank) || 'unknown'; counts.set(g, (counts.get(g) || 0) + 1); }
  const { map } = A.colourGroups(counts);
  const ann = leaf => annotations.get(leaf.name);
  const labelMode = opt.label || (annotations.size ? 'id-species' : 'id');
  const leafText = leaf => {
    const a = ann(leaf);
    switch (labelMode) {
      case 'species': return a?.species || leaf.name;
      case 'id-species': return a?.species ? `${leaf.name} ${a.species}` : leaf.name;
      case 'gene': return a?.genes?.[0] ? `${a.genes[0]} (${a.species || leaf.name})` : leaf.name;
      case 'protein': return a?.protein ? `${a.protein} · ${a.species || ''}` : leaf.name;
      default: return leaf.name;
    }
  };
  T.drawTree(svg, rooted, {
    width: parseInt(opt.width || 1100, 10),
    leafText,
    leafColor: annotations.size ? leaf => map.get(A.taxonAt(ann(leaf), rank) || 'unknown') || null : null,
    showEvents: !!events,
    domains: pfam.size || [...annotations.values()].some(a => a.domains?.length) ? leaf => {
      const doms = domainsOf(leaf.name);
      const len = ann(leaf)?.length || (doms.length ? Math.max(...doms.map(d => d.end)) : 0);
      return len ? { len, doms: doms.map(d => ({ ...d, color: A.domainColour(d.id || d.name), label: `${d.name}${d.id ? ` (${d.id})` : ''}` })) } : null;
    } : null,
  });
  out('.svg', `<svg xmlns="http://www.w3.org/2000/svg" viewBox="${svg.attrs.viewBox}" width="${svg.attrs.width}" height="${svg.attrs.height}" style="background:#fff">${svg.innerHTML}</svg>\n`);

  out('.post.json', JSON.stringify({ tree: opt.tree, tips: labels.length, root: how, rootNote: note,
    annotated: annotations.size, pfamProteins: pfam.size, events, outputs: written.concat(prefix + '.post.json') }, null, 2) + '\n');
  for (const f of written) console.log(f);
}
