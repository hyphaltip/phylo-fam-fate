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

/**
 * Midpoint-root a tree: place the root halfway along the longest tip-to-tip
 * path. FastTree trees are unrooted (drawn from an arbitrary trifurcation),
 * so this usually gives a far more readable picture. Support values describe
 * edges, so they travel with their edge when the direction flips.
 */
export function midpointRoot(root) {
  // Undirected adjacency: each edge {a, b, len, support}
  const nodes = collectAll(root);
  const adj = new Map(nodes.map(n => [n, []]));
  for (const n of nodes) {
    for (const c of n.children) {
      const e = { len: c.length, support: c.support };
      adj.get(n).push({ to: c, e });
      adj.get(c).push({ to: n, e });
    }
  }
  const leaves = nodes.filter(n => !n.children.length);
  if (leaves.length < 3) return root;

  function farthest(from) {
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
    let best = from;
    for (const l of leaves) if (dist.get(l) > dist.get(best)) best = l;
    return { best, dist, prev };
  }
  const A = farthest(leaves[0]).best;
  const { best: B, dist, prev } = farthest(A);
  const half = dist.get(B) / 2;
  if (!(half > 0)) return root;

  // Walk back from B toward A until we cross the midpoint.
  let v = B;
  while (prev.has(v) && dist.get(prev.get(v).u) > half) v = prev.get(v).u;
  const { u, e } = prev.get(v);          // midpoint lies on edge u—v
  const toV = dist.get(v) - half;         // distance from midpoint to v
  const toU = e.len - toV;

  function build(node, parent, length, support) {
    const out = { name: node.name, length, support, children: [] };
    for (const { to, e: edge } of adj.get(node)) {
      if (to === parent || edge === e) continue;
      out.children.push(build(to, node, edge.len, edge.support));
    }
    // A node that lost its only other neighbour keeps no stale support label
    if (!out.children.length) out.support = null;
    return out;
  }
  return {
    name: '', length: 0, support: null,
    children: [build(v, u, toV, e.support), build(u, v, toU, e.support)],
  };
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
 * @param {object} opts - {width, showSupport, cladogram, highlight, labelChars}
 */
export function drawTree(svg, root, opts = {}) {
  opts = { showSupport: true, showScale: true, cladogram: false, highlight: '', labelChars: 40, ...opts };
  const stats = treeStats(root);

  // Leaf depth in edges, used for cladogram layout (tips aligned at the right)
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
  const xOf = n => useLength ? n._x0 : root._h - n._h;
  const maxX = useLength ? Math.max(...all.map(n => n._x0)) : root._h;

  const longest = Math.min(opts.labelChars, Math.max(4, ...stats.names.map(n => n.length)));
  const fontPx = stats.tips > 400 ? 9 : stats.tips > 150 ? 10.5 : 12;
  const width = opts.width || Math.max(640, (svg.parentElement?.clientWidth || 900) - 4);
  const rowH = Math.max(fontPx + 3, Math.min(30, 560 / Math.max(1, stats.tips)));
  const margin = { top: 20, right: Math.min(width * 0.45, 14 + longest * fontPx * 0.62), bottom: useLength ? 40 : 15, left: 16 };
  const height_ = stats.tips * rowH + margin.top + margin.bottom;
  const treeW = Math.max(60, width - margin.left - margin.right);
  const xScale = treeW / (maxX || 1);

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

  const edgeParts = [];
  for (const node of all) {
    for (const c of node.children) {
      edgeParts.push(`M${node._px.toFixed(1)} ${c._y.toFixed(1)}H${c._px.toFixed(1)}`);
    }
    if (node.children.length > 1) {
      const ys = node.children.map(c => c._y);
      edgeParts.push(`M${node._px.toFixed(1)} ${Math.min(...ys).toFixed(1)}V${Math.max(...ys).toFixed(1)}`);
    }
  }
  parts.push(`<path d="${edgeParts.join('')}" class="t-edge"/>`);

  if (opts.showSupport) {
    for (const node of all) {
      if (node.support === null || !node.children.length || node === root) continue;
      const weak = node.support < 0.7 || (node.support > 1 && node.support < 70);
      parts.push(`<text x="${(node._px - 3).toFixed(1)}" y="${(node._y - 4).toFixed(1)}" class="t-support${weak ? ' weak' : ''}" text-anchor="end">${formatSupport(node.support)}</text>`);
    }
  }

  for (const node of all) {
    if (node.children.length) continue;
    const hit = q && node.name.toLowerCase().includes(q);
    if (hit) matches++;
    const y = node._y;
    const labelX = useLength ? node._px : margin.left + treeW;
    if (!useLength && labelX > node._px) {
      parts.push(`<line x1="${node._px}" y1="${y}" x2="${labelX}" y2="${y}" class="t-guide"/>`);
    }
    if (hit) parts.push(`<rect x="${labelX + 4}" y="${y - rowH / 2 + 1}" width="${Math.min(longest, node.name.length) * fontPx * 0.62 + 6}" height="${rowH - 2}" rx="3" class="t-hit"/>`);
    parts.push(`<circle cx="${labelX}" cy="${y}" r="2.4" class="t-tip"/>`);
    parts.push(`<text x="${labelX + 7}" y="${y + fontPx * 0.35}" class="t-label" font-size="${fontPx}"><title>${escapeXml(node.name)}</title>${escapeXml(truncate(node.name, opts.labelChars))}</text>`);
  }

  svg.setAttribute('viewBox', `0 0 ${width} ${height_}`);
  svg.setAttribute('width', width);
  svg.setAttribute('height', height_);
  svg.innerHTML = `<style>${TREE_CSS}</style>` + parts.join('');
  return { width, height: height_, stats, matches };
}

/* Inlined so a downloaded SVG looks the same as on screen. */
const TREE_CSS = `
.t-edge{fill:none;stroke:#334155;stroke-width:1.4;stroke-linecap:square}
.t-grid{stroke:#eef2f7;stroke-width:1}
.t-guide{stroke:#cbd5e1;stroke-width:0.8;stroke-dasharray:2 3}
.t-tick{font:10px system-ui,sans-serif;fill:#94a3b8}
.t-axis{font:11px system-ui,sans-serif;fill:#64748b}
.t-support{font:600 9px system-ui,sans-serif;fill:#7c3aed}
.t-support.weak{fill:#94a3b8;font-weight:400}
.t-tip{fill:#2563eb}
.t-label{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;fill:#1e293b}
.t-hit{fill:#fde68a}
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
