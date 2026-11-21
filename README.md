# phylo-fam-fate

A browser-based phylogenetics workbench. Upload sequences, search/align/trim them,
and infer a tree — all locally in the browser, with no server-side compute.

## Pipeline

```
          ┌──────────────────────────────────────────────┐
          │  Input: upload / paste / UniProt / NCBI /   │
          │         InterPro API retrieval              │
          └──────────────────┬───────────────────────────┘
                             ↓
          ┌──────────────────────────────────────────────┐
          │  Search (optional): local BLAST / phmmer    │
          │  (Biowasm WASM; DB from file, URL, or .gz)  │
          └──────────────────┬───────────────────────────┘
                             ↓
          ┌──────────────────────────────────────────────┐
          │  Align: Kalign / MUSCLE (Biowasm WASM)      │
          └──────────────────┬───────────────────────────┘
                             ↓
          ┌──────────────────────────────────────────────┐
          │  Trim: gap / conservation thresholds (JS)   │
          └──────────────────┬───────────────────────────┘
                             ↓
          ┌──────────────────────────────────────────────┐
          │  Tree: FastTree 2.1.11 WASM                 │
          │  (phylo-fasttree-wasm, direct-linked)        │
          └──────────────────────────────────────────────┘
```

Each step records **metadata**: tool, version, parameters, input/output sizes,
timing, and provenance.

## Data & limits

### Local (in-browser) data
- All compute runs in WebAssembly with a ~2 GB practical memory ceiling
- Alignments up to a few hundred MB are feasible; beyond that, use subsampling
- Biowasm's `CLI.mount()` streams remote files — large references need not be
  fully loaded into memory
- GitHub Pages: 1 GB site/repo limit, 100 MB max individual file

### Search database sources
Choose a **Database source** from the dropdown in step 2:

1. **UniProt Swiss-Prot** (default) — the reviewed UniProtKB FASTA
   (~94 MB gzipped, ~565k proteins) is streamed directly from
   `ftp.uniprot.org`. Biowasm builds a local BLAST index in the browser; the
   file is never fully loaded into JS memory.
2. **UniProt reference proteomes** — loads the
   [`reference_proteomes/README`](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README)
   catalogue (~37,600 proteomes) and shows a searchable dropdown. Filter by
   species name, proteome ID (`UP…`), or taxon ID, then click **Add to search
   set** (or double-click). Repeat to build a **search set** of several
   proteomes, optionally plus Swiss-Prot. Each proteome's
   `<ID>_<TAXID>.fasta.gz` is streamed into Biowasm and indexed as its own
   BLAST database; BLAST then searches them together (E-values reflect the
   combined size). Built indexes are cached for the session, so adding one
   more proteome only indexes the new one. Sets above ~1.5M proteins get a
   warning because in-browser indexing becomes slow.
3. **PDB SEQRES** — `pdb_seqres.txt.gz` from `files.wwpdb.org`. That server does
   **not** send CORS headers, so the browser cannot fetch it directly. Either
   download the file and use *Custom file*, or provide a CORS-proxy prefix
   (e.g. your own Cloudflare Worker).
4. **Custom file / URL** — pick a FASTA (`.fa`, `.fasta`, `.faa`, …) from your
   device, or paste an `https://…` URL and click **Load URL**. The server must
   send CORS headers (`Access-Control-Allow-Origin`). Compressed inputs work
   too: gzip is detected by magic bytes (`1f 8b`), not filename, and decompressed
   with the browser's native `DecompressionStream('gzip')`.

phmmer cannot stream a URL-backed DB, so a URL/custom DB is downloaded and
decompressed into browser memory on first phmmer run and then cached for later
runs. Nothing is ever uploaded to a server. The DB source, sequence count, and
sizes are recorded in the metadata table.

### Hit cutoffs and selection
The **Search E-value** (default `1e-3`) is passed to BLAST (`-evalue`) and
phmmer (`-E`/`--domE`), so weak hits are never reported. The live cutoffs below
can only tighten it; to see weaker hits, loosen it and search again.

After a search, every hit is shown in a table with **E-value**, **score**,
**% identity** (BLAST only — phmmer does not report identity), and **alignment
length**. Live cutoffs (`E-value ≤`, `Score ≥`, `% identity ≥`, `Align length ≥`)
filter the table as you type. Use **Select all** / **Select none** or individual
checkboxes, then **Use selected for alignment** to pass only those sequences to
step 3. phmmer hits are recovered from the in-memory DB by subject ID.

Hits are grouped **one row per subject protein**. When BLAST reports several
HSPs (local alignments) to the same protein, or several queries hit it, the row
shows the best HSP and an **“N HSPs ▸”** toggle that expands the individual
alignments with their subject coordinates. A protein passes the cutoffs if its
best HSP does. The **Subject range / length** column shows which part of the
protein the HSPs cover.

**Full-length sequences** (on by default) aligns each selected protein's whole
sequence. Databases are built with `makeblastdb -parse_seqids`, and after the
search `blastdbcmd -entry_batch` pulls every subject by accession from the same
in-browser databases, with no extra download. Untick it to align only the
regions BLAST matched (useful for single-domain trees): the best HSP, extended
by any overlapping or adjacent HSPs. A database whose headers can't be
ID-indexed (e.g. duplicate IDs) is rebuilt without `-parse_seqids`. Its hits
fall back to the HSP region and are flagged **HSP only**. phmmer hits are
always full-length.

### Runtime timing
Each step's elapsed time is printed to the browser console
(`[phylo-fam-fate] …`) and shown in the **timing strip** above the pipeline.
The same numbers appear in the run metadata table and in the exported JSON.

### Tree inference and viewer
FastTree runs in a Web Worker (`js/fasttree-worker.js`), so the page stays
responsive during long inferences; a live timer and **Cancel** button appear
while it runs. The tree view supports ladderizing, a cladogram layout, toggling
support values and finding tips by ID or species. Aligned and trimmed
alignments are shown in a colored viewer and can be downloaded as FASTA.

**Rooting** (`js/tree.js`). FastTree trees are unrooted, and every option
reroots the original tree, so switching between them is lossless. Support
values move with their bipartition.
- *Midpoint*: halfway along the longest tip-to-tip path (default).
- *MAD*, minimal ancestor deviation
  ([Tria et al. 2017](https://doi.org/10.1038/s41559-017-0193)): picks the
  root that makes the tree most clock-like across all tip pairs. It copes
  better than midpoint with uneven rates. The legend reports the deviation and
  ambiguity (a ratio near 1 means another root is almost as good). O(n³), so it
  is limited to 900 tips.
- *Outgroup taxon*: any taxon from the annotations, at any rank. If it isn't
  monophyletic, the tree is rooted on the split with the fewest tips on the
  wrong side, and the legend warns how many.
- *Chosen branch*: click any branch and choose **Root on this branch**.

**Annotation** (`js/annotate.js`). After a tree is built, leaves with UniProt
accessions (`sp|P46151|…`, `tr|Q7S830|…` or bare accessions) are annotated in
one batched call to UniProt's `uniprotkb/accessions` endpoint. It returns
organism, taxid, ranked lineage, protein and gene names, length, Pfam, domain
features, function and annotation score. Other labels fall back to
`OS=`/`OX=`/`GN=` in their FASTA description line, which needs no network and
gives species but no lineage. **Annotate from UniProt** reruns it.
- *Label*: ID, species, ID + species, gene or protein name.
- *Colour by*: domain of life, kingdom, phylum, subphylum, class, order,
  family, genus or species. Each leaf gets a swatch, and a branch is coloured
  when its whole clade shares one taxon, so a group that isn't a single clade
  stands out. The legend gives counts. The 12 most frequent groups get
  distinct colours; the rest are grey.
- *Click a leaf* for a card with protein and gene, species and lineage, a
  domain diagram, function text and a UniProt link. *Click a branch* to root
  on it, highlight its clade, or copy its tip names.
- *Domains*: Pfam domain architecture to a common scale beside each leaf.
  Coordinates come from the InterPro API, one request per protein with
  retries, up to 600 proteins; UniProt's own domain features are the fallback.
  Fusions, losses and truncated copies stand out next to the tree.
- *Duplications*: on a rooted, annotated tree, each internal node is labelled
  with the species-overlap rule: a duplication (red square) if its child
  clades share a species, otherwise a speciation (blue dot). The legend counts
  each. Results depend on the rooting, so compare midpoint, MAD and outgroup
  roots.

**Export** menu: SVG as shown; Newick as inferred (unrooted) or as shown
(rooted, ladderized); annotations TSV with one row per leaf (ranks, names,
Pfam, domain coordinates, function); iTOL colour-strip and protein-domain
datasets that match the Newick labels; and a TSV of duplication and
speciation nodes with their tips and species.

### Troubleshooting / debug log
A **Debug log** panel is pinned (collapsed) to the bottom-right of the page; a
red badge counts errors. It captures
console warnings/errors, uncaught exceptions, FastTree arguments and input
preview, and each step's timing. Use **Copy log** to share it. Step failures
also appear inline, directly underneath the relevant step (so an error at the
tree step is visible without scrolling back to the top status bar).

FastTree requires an **alignment**: if you build a tree from unaligned FASTA
with sequences of different lengths, step 5 shows an inline message asking you
to run the align step first.

### Remote APIs (CORS-enabled, callable directly from the browser)
| API | Use | Rate limit |
|---|---|---|
| [UniProt REST](https://rest.uniprot.org) | Protein sequence search & FASTA retrieval | ~10 req/s |
| [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/) | Nucleotide/protein search & retrieval | 3 req/s (no key) |
| [InterPro](https://www.ebi.ac.uk/interpro/api/) | Pfam/InterPro family search; Pfam domain coordinates for the tree | fair use (503s retried) |
| [UniProt `uniprotkb/accessions`](https://rest.uniprot.org/docs/) | Batched leaf annotation (≤400 accessions per request) | ~10 req/s |

NCBI BLAST URL API and EBI HMMER web service do **not** send CORS headers, so
remote BLAST/HMMER requires a proxy (planned).

## Run locally

```sh
python3 -m http.server 8000
# open http://localhost:8000/
```

No build step — the site uses ES modules and CDN-hosted tools.

## Deploy

GitHub Actions (`.github/workflows/deploy.yml`) deploys `main` to GitHub Pages.

## License

The site code is MIT. Bundled FastTree WASM is GPL-2.0-or-later (see
`assets/fasttree/LICENSE`). Biowasm tools carry their own licenses.
