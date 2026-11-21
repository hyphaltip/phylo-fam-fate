/* Newick parser and SVG phylogenetic tree renderer. */

/** Parse a Newick string into a nested tree object. */
export function parseNewick(s) {
  s = s.trim();
  let pos = 0;

  function readLabel() {
    let label = '';
    if (s[pos] === "'" || s[pos] === '"') {
      const quote = s[pos++];
      while (pos < s.length && s[pos] !== quote) {
        if (s[pos] === "''") label += "'", pos += 2;
        else label += s[pos++];
      }
      pos++; // closing quote
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
        const num = parseFloat(label);
        if (!isNaN(num) && label === String(num)) node.support = num;
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
 * Render a tree into an SVG element.
 * @param {SVGSVGElement} svg
 * @param {object} root - parsed tree
 * @param {object} opts - {width, showSupport, showLength, colorBy}
 */
export function drawTree(svg, root, opts = {}) {
  opts = { showSupport: true, showScale: true, ...opts };
  const stats = treeStats(root);

  // Compute cumulative x positions (branch lengths) and depth (topology)
  function assign(node, x, depth) {
    node._x0 = x;
    node._depth = depth;
    const childX = x + node.length;
    node.children.forEach(c => assign(c, childX, depth + 1));
  }
  assign(root, 0, 0);

  // For a cladogram (no meaningful lengths), use depth as x
  const useLength = stats.totalLength > 0.001;
  const maxX = useLength
    ? Math.max(...collectAll(root).map(n => n._x0 + n.length))
    : stats.maxDepth;

  const width = opts.width || Math.max(700, svg.parentElement?.clientWidth - 4 || 900);
  const rowH = Math.max(20, Math.min(38, 520 / Math.max(1, stats.tips)));
  const height = stats.tips * rowH + 70;
  const margin = { top: 25, right: 210, bottom: 35, left: 20 };
  const treeW = width - margin.left - margin.right - (opts.showScale ? 30 : 0);
  const xScale = treeW / (maxX || 1);

  // Leaf y positions, internal nodes at mean of children
  let leafIdx = 0;
  function layoutY(node) {
    if (!node.children.length) {
      node._y = margin.top + leafIdx++ * rowH + rowH / 2;
    } else {
      node.children.forEach(layoutY);
      node._y = (node.children[0]._y + node.children[node.children.length - 1]._y) / 2;
    }
    node._px = margin.left + (useLength ? node._x0 : node._depth) * xScale;
  }
  layoutY(root);

  const parts = [];

  // Grid lines / scale axis
  if (opts.showScale && useLength) {
    const niceStep = niceTick(maxX / 6);
    for (let v = 0; v <= maxX; v += niceStep) {
      const x = margin.left + v * xScale;
      parts.push(`<line x1="${x}" y1="${margin.top - 8}" x2="${x}" y2="${height - margin.bottom + 5}" stroke="#eef2f7" stroke-width="1"/>`);
      parts.push(`<text x="${x}" y="${height - margin.bottom + 20}" font-size="10" fill="#94a3b8" text-anchor="middle">${formatTick(v)}</text>`);
    }
    parts.push(`<text x="${margin.left + treeW / 2}" y="${height - 5}" font-size="11" fill="#64748b" text-anchor="middle">substitutions / site</text>`);
  }

  // Edges
  function edges(node) {
    for (const c of node.children) {
      const cx = margin.left + (useLength ? c._x0 : c._depth) * xScale;
      // elbow connector
      parts.push(`<path d="M ${node._px} ${node._y} L ${node._px} ${c._y} L ${cx} ${c._y}" fill="none" stroke="#475569" stroke-width="1.4" stroke-linejoin="round"/>`);
      edges(c);
    }
  }
  edges(root);

  // Support values on internal branches
  if (opts.showSupport) {
    function supports(node) {
      if (node.support !== null && node.children.length) {
        parts.push(`<text x="${node._px - 3}" y="${node._y - 4}" font-size="9" fill="#7c3aed" text-anchor="end" font-weight="600">${node.support}</text>`);
      }
      node.children.forEach(supports);
    }
    supports(root);
  }

  // Leaf labels with hover tooltip
  function labels(node) {
    if (!node.children.length) {
      const y = node._y;
      parts.push(`<circle cx="${node._px}" cy="${y}" r="2.5" fill="#2563eb"/>`);
      parts.push(`<text x="${node._px + 7}" y="${y + 4}" font-size="12" font-family="ui-monospace,monospace" fill="#1e293b"><title>${escapeXml(node.name)}</title>${escapeXml(truncate(node.name, 32))}</text>`);
    }
    node.children.forEach(labels);
  }
  labels(root);

  svg.setAttribute('viewBox', `0 0 ${width} ${height}`);
  svg.setAttribute('preserveAspectRatio', 'xMinYMin meet');
  svg.innerHTML = parts.join('');
  return { width, height, stats };
}

function collectAll(node) {
  const out = [node];
  node.children.forEach(c => out.push(...collectAll(c)));
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
