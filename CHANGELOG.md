# Changelog

## v0.1.0 (2026-11-21)

First release. The browser workbench (https://hyphaltip.github.io/phylo-fam-fate/)
and its command-line twin, `famfate.sh`, share code and default settings.

### Browser
- **Input:** upload or paste FASTA, or retrieve from UniProt, NCBI or InterPro.
  Examples include a single-query MET12 (P46151) workflow.
- **Search:** in-browser BLASTP/BLASTN (BLAST+ 2.17.0) and phmmer (HMMER 3.4).
  - Databases: Swiss-Prot, one or more UniProt reference proteomes searched
    together, PDB, or a custom file or URL. Indexes are cached for the
    session.
  - Search E-value 1e-3 by default, with live cutoffs in the hit table.
  - Hits are grouped one row per protein, with HSPs nested underneath.
  - Full-length subjects are retrieved with `blastdbcmd`.
- **Alignment:** Kalign 3.3.1, MUSCLE 5.1, MAFFT 7.520 E-INS-i, with a
  coloured alignment viewer. **Trimming** by gap fraction and conservation.
- **Tree:** FastTree 2.1.11 compiled to WebAssembly, in a worker that can be
  cancelled. **Single or double precision**; both run at the same speed in
  WASM.
- **Tree viewer:**
  - Rooting: midpoint, MAD, outgroup taxon, or a clicked branch.
  - Ladderize, cladogram layout, and finding tips by ID or species.
  - Loading external Newick trees (IQ-TREE, RAxML-NG, …), keeping paired
    `SH-aLRT/UFBoot` support labels.
- **Annotation:**
  - Batched UniProt lookup: species, ranked lineage, protein/gene, Pfam,
    function. Offline fallback to `OS=` in FASTA headers.
  - Colour by any rank, with clade-aware branch colouring.
  - Leaf and branch popups.
  - Pfam domain track from InterPro.
  - Duplication/speciation labelling by species overlap.
- **Export:** SVG, Newick (unrooted or rooted), annotation TSV, iTOL
  colour-strip and domain datasets, duplication-node TSV, and run-metadata
  JSON.
- **Hand-off to IQ-TREE / RAxML-NG:** downloads a self-contained script, with
  the alignment embedded, that runs on your own machine or cluster.
- Debug log panel and per-step timing.

### Command line (`cli/famfate.sh`) and container
- One-command pipeline with the browser's defaults: query → search
  (Swiss-Prot and/or proteomes) → full-length hits → align → trim → tree →
  root → annotate → duplications → exports, plus `run.json` provenance.
- Tree methods: FastTree (`--fasttree browser`, byte-identical to the site,
  or `--fasttree double`), IQ-TREE 3.1.4 (ModelFinder + UFBoot/SH-aLRT) and
  RAxML-NG 2.0.4 (bootstrap, TBE). `--tree iqtree,raxml-ng,fasttree` compares
  them in `trees_compare.tsv`.
- `famfate-js` runs the site's own JavaScript for trimming, rooting (incl.
  MAD), annotation, SVG and iTOL export.
- Container `ghcr.io/hyphaltip/phylo-fam-fate-cli:0.1.0` (linux/amd64), with
  tools pinned to the browser's versions. Kalign 3.3.1 is built from source
  without CPU-specific AVX code, because bioconda's kalign3 binary crashed
  with "illegal instruction" on some CPUs. It works with Docker and with
  Apptainer/Singularity (`apptainer pull docker://…`).
- CI runs the following in the container on every change:
  - the offline pipelines;
  - FastTree native-vs-WASM parity;
  - a live MET12 search;
  - ML trees;
  - an Apptainer run.

### Known limitations
- The container is x86-64 only, because BLAST+ 2.17.0 has no linux-arm64
  build on bioconda.
- MAFFT in the browser runs MAFFT's internal programs directly. The CLI's
  `mafft --genafpair --ep 0 --maxiterate 16` is equivalent, not
  byte-identical.
- Annotation needs UniProt accessions in leaf labels, or `OS=` in FASTA
  headers. NCBI taxonomy lookup is not implemented yet.
- Browser memory limits mean very large searches and alignments belong on the
  CLI (guidance in the README).
