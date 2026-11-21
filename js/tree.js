/* Newick parser and SVG phylogenetic tree renderer. */

const NUMERIC = /^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$/;

/** Parse a Newick string into a nested tree object. */
export function parseNewick(s) {
  s = s.trim();
  let pos = 0;

  function readLabel() {
    let label = '';
    if (s[pos] === "'" || s[pos] === '"') {
      const quote = s[pos++];
      while (pos < s.length) {
        if (s[pos] === quote) {
          if (s[pos + 1] === quote) { label += quote; pos += 2; continue; } // escaped ''
          break;
        }
        label += s[pos++];
      }
      pos++; // closing quote
      return label;
    } else {
      while (pos < s.length && !':(),;[]'.includes(s[pos])) label += s[pos++];
    }
    return label.trim();
  }

  function parseNode() {
    const node = { name: '', length: 0, support: null, children: [] };
    if (s[pos] === '[') {
      // skip square-bracket comments
      let depth = 0;
      while (pos < s.length) {
        if (s[pos] === '[') depth++;
        else if (s[pos] === ']') { depth--; if (depth === 0) { pos++; break; } }
        pos++;
      }
    }
    if (s[pos] === '(') {
      pos++;
      while (pos < s.length && s[pos] !== ')') {
        node.children.push(parseNode());
        if (s[pos] === ',') pos++;
      }
      pos++; // ')'
      const label = readLabel();
      if (label) {
        // FastTree writes support as e.g. "0.950"; String(0.95) !== "0.950",
        // so test the lexical form rather than round-tripping the number.
        if (NUMERIC.test(label)) node.support = parseFloat(label);
        else node.name = label;
      }
    } else {
      node.name = readLabel();
    }
    if (s[pos] === ':') {
      pos++;
      let num = '';
      while (pos < s.length && !',();[]'.includes(s[pos])) num += s[pos++];
      node.length = parseFloat(num) || 0;
    }
    return node;
  }

  const root = parseNode();
  return root;
}

/** Compute summary statistics for a tree. */
export function treeStats(root) {
  const leaves = [];
  let nodes = 0, maxDepth = 0, totalLength = 0;
  function walk(node, depth) {
    nodes++;
    maxDepth = Math.max(maxDepth, depth);
    totalLength += node.length;
    if (!node.children.length) leaves.push(node);
    else node.children.forEach(c => walk(c, depth + 1));
  }
  walk(root, 0);
  return {
    tips: leaves.length,
    nodes,
    maxDepth,
    totalLength,
    names: leaves.map(l => l.name),
  };
}

/* ---------------------------------------------------------------------------
 * Rooting. FastTree trees are unrooted (drawn from an arbitrary trifurcation).
 * Every method finds an edge plus a position on it and rebuilds the tree from
 * that point. Support values describe edges (bipartitions), so they travel
 * with their edge when its direction flips.
 * ------------------------------------------------------------------------- */

/** Undirected view of a tree: adjacency lists of {to, e}, e = shared edge. */
function undirected(root) {
  const nodes = collectAll(root);
  const adj = new Map(nodes.map(n => [n, []]));
  for (const n of nodes) {
    for (const c of n.children) {
      const e = { len: c.length, support: c.support, child: c, parent: n };
      adj.get(n).push({ to: c, e });
      adj.get(c).push({ to: n, e });
    }
  }
  return { nodes, adj, leaves: nodes.filter(n => !n.children.length) };
}

/** Rebuild a tree rooted on edge `e` (between e.child and e.parent), with the
 *  new root `toChild` away from the child end. Extra node fields (e.g. `id`)
 *  are copied so callers can keep referring to nodes. */
function rootOnEdge(adj, e, toChild) {
  toChild = Math.max(0, Math.min(e.len, toChild));
  function build(node, from, length, support) {
    const out = { ...node, length, support, children: [] };
    for (const { to, e: edge } of adj.get(node)) {
      if (to === from || edge === e) continue;
      out.children.push(build(to, node, edge.len, edge.support));
    }
    // A former internal node that is now a leaf can't keep a support label;
    // nor can a degree-2 node (old root) which simply passes through.
    if (!out.children.length) out.support = null;
    return out;
  }
  const a = build(e.child, e.parent, toChild, e.support);
  const b = build(e.parent, e.child, e.len - toChild, e.support);
  return { name: '', length: 0, support: null, children: [a, b] };
}

/** Collapse degree-2 internal nodes (the old root after rerooting). */
function suppressUnary(node) {
  node.children = node.children.map(c => {
    while (c.children.length === 1) {
      const only = c.children[0];
      only.length += c.length;
      c = only;
    }
    return suppressUnary(c);
  });
  return node;
}

/** Path-based distances from one node to every other node. */
function distancesFrom(adj, from) {
  const dist = new Map([[from, 0]]), prev = new Map();
  const stack = [from];
  while (stack.length) {
    const u = stack.pop();
    for (const { to, e } of adj.get(u)) {
      if (dist.has(to)) continue;
      dist.set(to, dist.get(u) + e.len);
      prev.set(to, { u, e });
      stack.push(to);
    }
  }
  return { dist, prev };
}

/** Midpoint root: halfway along the longest tip-to-tip path. */
export function midpointRoot(root) {
  const { adj, leaves } = undirected(root);
  if (leaves.length < 3) return root;
  const far = (from) => {
    const { dist, prev } = distancesFrom(adj, from);
    let best = from;
    for (const l of leaves) if (dist.get(l) > dist.get(best)) best = l;
    return { best, dist, prev };
  };
  const A = far(leaves[0]).best;
  const { best: B, dist, prev } = far(A);
  const half = dist.get(B) / 2;
  if (!(half > 0)) return root;
  // Walk back from B toward A until we cross the midpoint.
  let v = B;
  while (prev.has(v) && dist.get(prev.get(v).u) > half) v = prev.get(v).u;
  const { e } = prev.get(v);            // midpoint lies on edge u—v
  const toV = dist.get(v) - half;        // distance from midpoint to v
  return suppressUnary(rootOnEdge(adj, e, e.child === v ? toV : e.len - toV));
}

/** Count, for each node, its leaves and how many are in `set` (by name). */
function cladeCounts(root, set) {
  const all = collectAll(root);
  const m = new Map(), k = new Map();
  for (let i = all.length - 1; i >= 0; i--) {    // children before parents
    const n = all[i];
    if (!n.children.length) { m.set(n, 1); k.set(n, set.has(n.name) ? 1 : 0); continue; }
    let mm = 0, kk = 0;
    for (const c of n.children) { mm += m.get(c); kk += k.get(c); }
    m.set(n, mm); k.set(n, kk);
  }
  return { all, m, k, N: m.get(root) };
}

/**
 * Root on the edge that best separates the leaves named in `outgroup` from
 * the rest. Works on any orientation of the same tree, because an edge is
 * identified by the leaf split it induces.
 * @returns {{tree, exact: boolean, inGroup: number, impurity: number}}
 *   exact = the outgroup is a clade; otherwise `impurity` leaves are on the
 *   wrong side of the chosen edge.
 */
export function rootOnOutgroup(root, outgroup) {
  const set = outgroup instanceof Set ? outgroup : new Set(outgroup);
  const { adj } = undirected(root);
  const { all, m, k, N } = cladeCounts(root, set);
  const S = k.get(root);
  if (!S || S === N) return { tree: root, exact: false, inGroup: S, impurity: 0, invalid: true };
  let best = null;
  for (const c of all) {
    if (c === root) continue;
    // Outgroup either below c (k of m) or above it (S-k of N-m).
    const below = { wrong: (m.get(c) - k.get(c)) + (S - k.get(c)) };
    const above = { wrong: k.get(c) + ((N - m.get(c)) - (S - k.get(c))) };
    const w = Math.min(below.wrong, above.wrong);
    if (!best || w < best.w) best = { c, w };
  }
  const edge = adj.get(best.c).find(x => x.e.child === best.c).e;
  return { tree: suppressUnary(rootOnEdge(adj, edge, edge.len / 2)), exact: best.w === 0, inGroup: S, impurity: best.w };
}

/** Root on the edge above the clade of leaves named in `clade` (exact split). */
export function rootOnClade(root, clade) {
  return rootOnOutgroup(root, clade).tree;
}

/**
 * Minimal Ancestor Deviation rooting (Tria, Landan & Dagan 2017, Nat Ecol
 * Evol 1:0193). For every candidate edge and position, each tip pair's
 * ancestor should sit at the midpoint of their path; the root minimising the
 * squared relative deviation from that clock-like expectation is chosen.
 * O(n³) overall, so it is capped at `maxTips`.
 * @returns {{tree, score: number, ambiguity: number} | null}
 */
export function madRoot(root, maxTips = 900) {
  const { nodes, adj, leaves } = undirected(root);
  const n = leaves.length;
  if (n < 3 || n > maxTips) return null;
  const idx = new Map(nodes.map((x, i) => [x, i]));
  // D[l][node] = path length from leaf l to node
  const D = leaves.map(l => {
    const { dist } = distancesFrom(adj, l);
    const row = new Float64Array(nodes.length);
    for (const [node, d] of dist) row[idx.get(node)] = d;
    return row;
  });
  const leafIdx = leaves.map(l => idx.get(l));
  const dij = (i, j) => D[i][leafIdx[j]];

  // Membership of each leaf below each edge's child (original orientation).
  const below = new Map();
  for (let i = nodes.length - 1; i >= 0; i--) {
    const x = nodes[i];
    if (!x.children.length) { const b = new Uint8Array(n); b[leaves.indexOf(x)] = 1; below.set(x, b); continue; }
    const b = new Uint8Array(n);
    for (const c of x.children) { const cb = below.get(c); for (let t = 0; t < n; t++) b[t] |= cb[t]; }
    below.set(x, b);
  }

  const results = [];
  for (const x of nodes) {
    if (x === root) continue;
    const e = adj.get(x).find(a => a.e.child === x).e;
    const v = idx.get(e.child), u = idx.get(e.parent), L = e.len;
    const inC = below.get(x);
    // Sums for the cross pairs, whose ancestor is the root itself:
    // r = 2(dv_i + t)/d − 1, minimised analytically over t ∈ [0, L].
    let sA = 0, sB = 0, fixedSq = 0, pairs = 0;
    const cross = [];
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
      const d = dij(i, j);
      if (!(d > 0)) continue;
      pairs++;
      if (inC[i] !== inC[j]) {
        const a = inC[i] ? D[i][v] : D[j][v];   // the tip on the child side
        cross.push([a, d]);
        sA += (d - 2 * a) / (d * d);
        sB += 1 / (d * d);
      } else {
        // Same side: the ancestor is where their path is closest to the root
        // end of this edge (v for the child side, u for the rest).
        const end = inC[i] ? v : u;
        const anc = (D[i][end] + d - D[j][end]) / 2;
        const r = 2 * anc / d - 1;
        fixedSq += r * r;
      }
    }
    if (!pairs) continue;
    let t = sB > 0 ? sA / (2 * sB) : L / 2;
    t = Math.max(0, Math.min(L, t));
    let sq = fixedSq;
    for (const [a, d] of cross) { const r = 2 * (a + t) / d - 1; sq += r * r; }
    results.push({ e, t, score: Math.sqrt(sq / pairs) });
  }
  if (!results.length) return null;
  results.sort((a, b) => a.score - b.score);
  const best = results[0];
  const ambiguity = results.length > 1 ? best.score / results[1].score : 0;
  return { tree: suppressUnary(rootOnEdge(adj, best.e, best.t)), score: best.score, ambiguity };
}

/**
 * Label internal nodes of a *rooted* tree as duplications or speciations with
 * the species-overlap rule: a node is a duplication when two of its child
 * clades share a species. Sets node.event = 'dup' | 'spec' (null when the
 * species of a subtree is unknown for every leaf).
 * @param {(leaf) => string|null} speciesOf
 * @returns {{dup: number, spec: number}}
 */
export function labelEvents(root, speciesOf) {
  let dup = 0, spec = 0;
  const walk = node => {
    if (!node.children.length) {
      const s = speciesOf(node);
      return s ? new Set([s]) : new Set();
    }
    const sets = node.children.map(walk);
    const union = new Set();
    let overlap = false, known = 0;
    for (const st of sets) {
      if (st.size) known++;
      for (const x of st) { if (union.has(x)) overlap = true; union.add(x); }
    }
    node.event = overlap ? 'dup' : known >= 2 ? 'spec' : null;
    if (node.event === 'dup') dup++;
    else if (node.event === 'spec') spec++;
    return union;
  };
  walk(root);
  return { dup, spec };
}

/** Names of the leaves below a node. */
export function leafNames(node) {
  return collectAll(node).filter(n => !n.children.length).map(n => n.name);
}

/** Sort children so larger clades come last ("ladderize"), like FigTree. */
export function ladderize(node) {
  if (!node.children.length) return 1;
  const sizes = new Map(node.children.map(c => [c, ladderize(c)]));
  node.children.sort((a, b) => sizes.get(a) - sizes.get(b));
  return [...sizes.values()].reduce((x, y) => x + y, 0);
}

/**
 * Render a tree into an SVG element.
 * @param {SVGSVGElement} svg
 * @param {object} root - parsed tree
 * @param {object} opts
 *   showSupport, cladogram, highlight, labelChars, showScale
 *   leafText(leaf)   → label string (default: leaf.name)
 *   leafColor(leaf)  → CSS colour or null (swatch + branch colouring)
 *   colorBranches    → colour edges whose whole clade shares one colour
 *   domains(leaf)    → {len, doms:[{start,end,color,label}]} | null; enables
 *                      a domain-architecture track to the right of the labels
 *   showEvents       → draw node.event ('dup'|'spec') markers
 *   selected(node)   → highlight that node's clade
 * Every edge, node and leaf carries data-n = index into the returned `nodes`.
 */
export function drawTree(svg, root, opts = {}) {
  opts = { showSupport: true, showScale: true, cladogram: false, highlight: '', labelChars: 40,
    leafText: l => l.name, leafColor: null, colorBranches: true, domains: null, showEvents: false, selected: null, ...opts };
  const stats = treeStats(root);

  function height(node) {
    node._h = node.children.length ? 1 + Math.max(...node.children.map(height)) : 0;
    return node._h;
  }
  height(root);
  function assign(node, x, depth) {
    node._x0 = x;
    node._depth = depth;
    node.children.forEach(c => assign(c, x + c.length, depth + 1));
  }
  assign(root, 0, 0);

  const useLength = !opts.cladogram && stats.totalLength > 1e-9;
  const all = collectAll(root);
  const leaves = all.filter(n => !n.children.length);
  const xOf = n => useLength ? n._x0 : root._h - n._h;
  const maxX = useLength ? Math.max(...all.map(n => n._x0)) : root._h;

  const texts = new Map(leaves.map(l => [l, String(opts.leafText(l) ?? l.name)]));
  const colours = new Map(opts.leafColor ? leaves.map(l => [l, opts.leafColor(l)]) : []);
  const tracks = new Map(opts.domains ? leaves.map(l => [l, opts.domains(l)]) : []);
  const hasTrack = [...tracks.values()].some(Boolean);
  const maxLen = hasTrack ? Math.max(1, ...[...tracks.values()].filter(Boolean).map(t => t.len)) : 1;

  const longest = Math.min(opts.labelChars, Math.max(4, ...[...texts.values()].map(t => t.length)));
  const fontPx = stats.tips > 400 ? 9 : stats.tips > 150 ? 10.5 : 12;
  const width = opts.width || Math.max(640, (svg.parentElement?.clientWidth || 900) - 4);
  const rowH = Math.max(fontPx + 3, Math.min(30, 560 / Math.max(1, stats.tips)));
  const swatch = colours.size ? fontPx * 0.8 + 4 : 0;
  const labelW = 14 + swatch + longest * fontPx * 0.62;
  const trackW = hasTrack ? Math.min(260, Math.max(140, width * 0.22)) : 0;
  const margin = { top: 20, right: Math.min(width * 0.6, labelW + (hasTrack ? trackW + 16 : 0)), bottom: useLength ? 40 : 15, left: 16 };
  const height_ = stats.tips * rowH + margin.top + margin.bottom;
  const treeW = Math.max(60, width - margin.left - margin.right);
  const xScale = treeW / (maxX || 1);
  const trackX = width - trackW - 8;

  let leafIdx = 0;
  function layoutY(node) {
    if (!node.children.length) {
      node._y = margin.top + leafIdx++ * rowH + rowH / 2;
    } else {
      node.children.forEach(layoutY);
      node._y = (node.children[0]._y + node.children[node.children.length - 1]._y) / 2;
    }
    node._px = margin.left + xOf(node) * xScale;
  }
  layoutY(root);

  // Clade colour: a node takes a colour only when every leaf below shares it.
  if (colours.size && opts.colorBranches) {
    for (let i = all.length - 1; i >= 0; i--) {
      const n = all[i];
      if (!n.children.length) { n._col = colours.get(n) || null; continue; }
      const c0 = n.children[0]._col;
      n._col = c0 && n.children.every(c => c._col === c0) ? c0 : null;
    }
  }
  const nIndex = new Map(all.map((n, i) => [n, i]));
  const sel = opts.selected ? all.filter(opts.selected) : [];

  const q = opts.highlight.trim().toLowerCase();
  let matches = 0;
  const parts = [];

  if (opts.showScale && useLength) {
    const niceStep = niceTick(maxX / 6);
    for (let v = 0; v <= maxX + 1e-12; v += niceStep) {
      const x = margin.left + v * xScale;
      parts.push(`<line x1="${x}" y1="${margin.top - 8}" x2="${x}" y2="${height_ - margin.bottom + 5}" class="t-grid"/>`);
      parts.push(`<text x="${x}" y="${height_ - margin.bottom + 19}" class="t-tick" text-anchor="middle">${formatTick(v)}</text>`);
    }
    parts.push(`<text x="${margin.left + treeW / 2}" y="${height_ - 4}" class="t-axis" text-anchor="middle">substitutions / site</text>`);
  }

  // Selected clade background
  for (const n of sel) {
    const ls = collectAll(n).filter(x => !x.children.length);
    const y0 = Math.min(...ls.map(l => l._y)) - rowH / 2, y1 = Math.max(...ls.map(l => l._y)) + rowH / 2;
    parts.push(`<rect x="${(n._px - (n === root ? 0 : (n.length * xScale || 4))).toFixed(1)}" y="${y0.toFixed(1)}" width="${(width - n._px).toFixed(1)}" height="${(y1 - y0).toFixed(1)}" class="t-sel"/>`);
  }

  // Edges, grouped into one path per colour
  const byCol = new Map();
  const addSeg = (col, d) => (byCol.get(col) || byCol.set(col, []).get(col)).push(d);
  const hits = [];
  for (const node of all) {
    for (const c of node.children) {
      addSeg(c._col || '', `M${node._px.toFixed(1)} ${c._y.toFixed(1)}H${c._px.toFixed(1)}`);
      hits.push(`<path d="M${node._px.toFixed(1)} ${c._y.toFixed(1)}H${Math.max(c._px, node._px + 3).toFixed(1)}" class="t-hitarea" data-n="${nIndex.get(c)}"/>`);
    }
    if (node.children.length > 1) {
      const ys = node.children.map(c => c._y);
      addSeg(node._col || '', `M${node._px.toFixed(1)} ${Math.min(...ys).toFixed(1)}V${Math.max(...ys).toFixed(1)}`);
    }
  }
  for (const [col, ds] of byCol) {
    parts.push(`<path d="${ds.join('')}" class="t-edge${col ? ' t-col' : ''}"${col ? ` style="stroke:${col}"` : ''}/>`);
  }
  parts.push(...hits);

  if (opts.showSupport) {
    for (const node of all) {
      if (node.support === null || !node.children.length || node === root) continue;
      const weak = node.support < 0.7 || (node.support > 1 && node.support < 70);
      parts.push(`<text x="${(node._px - 3).toFixed(1)}" y="${(node._y - 4).toFixed(1)}" class="t-support${weak ? ' weak' : ''}" text-anchor="end">${formatSupport(node.support)}</text>`);
    }
  }

  // Duplication (square) / speciation (circle) markers on internal nodes
  if (opts.showEvents) {
    for (const node of all) {
      if (!node.children.length || !node.event) continue;
      const x = node._px.toFixed(1), y = node._y.toFixed(1);
      parts.push(node.event === 'dup'
        ? `<rect x="${(node._px - 4).toFixed(1)}" y="${(node._y - 4).toFixed(1)}" width="8" height="8" class="t-dup" data-n="${nIndex.get(node)}"><title>Duplication</title></rect>`
        : `<circle cx="${x}" cy="${y}" r="3.2" class="t-spec" data-n="${nIndex.get(node)}"><title>Speciation</title></circle>`);
    }
  }

  for (const node of leaves) {
    const text = texts.get(node);
    const hit = q && (node.name.toLowerCase().includes(q) || text.toLowerCase().includes(q));
    if (hit) matches++;
    const y = node._y;
    const labelX = useLength ? node._px : margin.left + treeW;
    if (!useLength && labelX > node._px) {
      parts.push(`<line x1="${node._px}" y1="${y}" x2="${labelX}" y2="${y}" class="t-guide"/>`);
    }
    const col = colours.get(node);
    const tx = labelX + 7 + swatch;
    if (hit) parts.push(`<rect x="${tx - 3}" y="${y - rowH / 2 + 1}" width="${Math.min(longest, text.length) * fontPx * 0.62 + 6}" height="${rowH - 2}" rx="3" class="t-hit"/>`);
    parts.push(`<circle cx="${labelX}" cy="${y}" r="2.4" class="t-tip"${col ? ` style="fill:${col}"` : ''}/>`);
    if (swatch) {
      const sz = fontPx * 0.8;
      parts.push(`<rect x="${labelX + 6}" y="${(y - sz / 2).toFixed(1)}" width="${sz.toFixed(1)}" height="${sz.toFixed(1)}" rx="2" style="fill:${col || '#e2e8f0'}"/>`);
    }
    parts.push(`<text x="${tx}" y="${y + fontPx * 0.35}" class="t-label" font-size="${fontPx}" data-n="${nIndex.get(node)}"><title>${escapeXml(node.name === text ? text : `${text}\n${node.name}`)}</title>${escapeXml(truncate(text, opts.labelChars))}</text>`);

    // Domain architecture, drawn to a common scale (longest protein = trackW)
    const tr = tracks.get(node);
    if (tr) {
      const sx = trackW / maxLen;
      const h = Math.max(4, Math.min(12, rowH * 0.55));
      parts.push(`<line x1="${trackX}" y1="${y}" x2="${(trackX + tr.len * sx).toFixed(1)}" y2="${y}" class="t-prot"/>`);
      for (const d of tr.doms) {
        parts.push(`<rect x="${(trackX + (d.start - 1) * sx).toFixed(1)}" y="${(y - h / 2).toFixed(1)}" width="${Math.max(1.5, (d.end - d.start + 1) * sx).toFixed(1)}" height="${h.toFixed(1)}" rx="2.5" style="fill:${d.color}" class="t-dom" data-n="${nIndex.get(node)}"><title>${escapeXml(`${d.label} ${d.start}–${d.end}`)}</title></rect>`);
      }
    }
  }
  if (hasTrack) {
    parts.push(`<text x="${trackX}" y="${margin.top - 6}" class="t-tick">domains (scale: ${maxLen} aa)</text>`);
  }

  svg.setAttribute('viewBox', `0 0 ${width} ${height_}`);
  svg.setAttribute('width', width);
  svg.setAttribute('height', height_);
  svg.innerHTML = `<style>${TREE_CSS}</style>` + parts.join('');
  return { width, height: height_, stats, matches, nodes: all };
}

/* Inlined so a downloaded SVG looks the same as on screen. */
const TREE_CSS = `
.t-edge{fill:none;stroke:#334155;stroke-width:1.4;stroke-linecap:square}
.t-edge.t-col{stroke-width:2}
.t-hitarea{fill:none;stroke:transparent;stroke-width:9;cursor:pointer}
.t-hitarea:hover{stroke:rgba(37,99,235,.18)}
.t-sel{fill:#dbeafe;opacity:.55}
.t-grid{stroke:#eef2f7;stroke-width:1}
.t-guide{stroke:#cbd5e1;stroke-width:0.8;stroke-dasharray:2 3}
.t-tick{font:10px system-ui,sans-serif;fill:#94a3b8}
.t-axis{font:11px system-ui,sans-serif;fill:#64748b}
.t-support{font:600 9px system-ui,sans-serif;fill:#7c3aed}
.t-support.weak{fill:#94a3b8;font-weight:400}
.t-tip{fill:#2563eb}
.t-label{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;fill:#1e293b;cursor:pointer}
.t-label:hover{text-decoration:underline}
.t-hit{fill:#fde68a}
.t-dup{fill:#dc2626;stroke:#fff;stroke-width:1;cursor:pointer}
.t-spec{fill:#2563eb;stroke:#fff;stroke-width:1;cursor:pointer}
.t-prot{stroke:#94a3b8;stroke-width:1.5}
.t-dom{stroke:rgba(0,0,0,.25);stroke-width:.6;cursor:pointer}
`;

function formatSupport(v) {
  return v <= 1 ? v.toFixed(2).replace(/^0/, '') : String(Math.round(v));
}

function collectAll(root) {
  const out = [], stack = [root];
  while (stack.length) {
    const n = stack.pop();
    out.push(n);
    for (let i = n.children.length - 1; i >= 0; i--) stack.push(n.children[i]);
  }
  return out;
}
function niceTick(x) {
  const mag = Math.pow(10, Math.floor(Math.log10(x || 1)));
  const norm = x / mag;
  return (norm < 1.5 ? 1 : norm < 3 ? 2 : norm < 7 ? 5 : 10) * mag;
}
function formatTick(v) {
  if (v === 0) return '0';
  return v < 0.01 ? v.toExponential(1) : String(+v.toFixed(4));
}
function truncate(s, n) { return s.length > n ? s.slice(0, n - 1) + '…' : s; }
function escapeXml(s) {
  return String(s).replace(/[<>&"']/g, c => ({ '<': '&lt;', '>': '&gt;', '&': '&amp;', '"': '&quot;', "'": '&apos;' }[c]));
}
