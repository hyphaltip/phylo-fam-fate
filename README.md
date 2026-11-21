# phylo-fam-fate

A browser-based workbench for gene-family phylogenetics. You can start from a
single protein, find its homologs in Swiss-Prot or any set of UniProt reference
proteomes, then align, trim and infer a tree. The tree can be annotated with
taxonomy, Pfam domains and function, rooted several ways, and marked with
gene duplications. Everything runs **in your browser** (WebAssembly), with no
server, account or upload.

**Live site:** <https://hyphaltip.github.io/phylo-fam-fate/> · **Current release:**
[v0.1.0](https://github.com/hyphaltip/phylo-fam-fate/releases/tag/v0.1.0)
([changelog](CHANGELOG.md); container `ghcr.io/hyphaltip/phylo-fam-fate-cli:0.1.0`)

Every step uses a standard tool (BLAST+, HMMER, Kalign, MUSCLE, FastTree) at a
recorded version and parameters. This README gives the **equivalent
command-line run** for each step, so a browser analysis can be reproduced,
scaled up or scripted on a workstation or cluster.

---

## Contents

1. [Two ways to run it](#two-ways-to-run-it)
2. [Quick start: MET12 in the browser](#quick-start-met12-in-the-browser)
3. [The pipeline, step by step](#the-pipeline-step-by-step) (browser use and CLI equivalent)
   - [Step 1 · Input sequences](#step-1--input-sequences)
   - [Step 2 · Homology search](#step-2--homology-search-optional)
   - [Step 3 · Multiple sequence alignment](#step-3--multiple-sequence-alignment)
   - [Step 4 · Alignment trimming](#step-4--alignment-trimming)
   - [Step 5 · Tree inference](#step-5--tree-inference)
   - [Step 6 · Rooting](#step-6--rooting)
   - [Step 7 · Annotation: taxonomy, domains, function](#step-7--annotation-taxonomy-domains-function)
   - [Step 8 · Duplications vs speciations](#step-8--duplications-vs-speciations)
   - [Step 9 · Export and provenance](#step-9--export-and-provenance)
4. [Command-line tool: `famfate.sh`](#command-line-tool-famfatesh) (container, Apptainer/Singularity, conda)
5. [Browser vs command line: choosing](#browser-vs-command-line-choosing)
6. [Data sources and limits](#data-sources-and-limits)
7. [Troubleshooting](#troubleshooting)
8. [How it works (architecture)](#how-it-works-architecture)
9. [Run locally / deploy / license](#run-locally)

---

## Two ways to run it

| | **Browser** (this site) | **Command line** (`famfate.sh`) |
|---|---|---|
| Install | Nothing; open the page | Container (Docker / Apptainer) or conda environment |
| Compute | Your machine, in WebAssembly (single-threaded) | Native, multi-threaded |
| Data | Streamed from UniProt into the browser; never uploaded | Downloaded to disk |
| Best for | Exploring one family; teaching; quick looks; sharing a method | Large families (>1–2k sequences), many proteomes, batch runs, publication runs |
| Output | SVG, Newick, FASTA, TSV, iTOL files, run-metadata JSON | Same formats |

Tool versions used in the browser and pinned in the CLI container:

| Step | Tool | Browser | CLI container |
|---|---|---|---|
| Search | BLAST+ `blastp` / `blastn`, `makeblastdb`, `blastdbcmd` | 2.17.0 | 2.17.0 |
| Search | HMMER `phmmer` | 3.4 | 3.4 |
| Align | Kalign | 3.3.1 | 3.3.1, built from source without AVX autodetection (`cli/build-kalign-portable.sh`): bioconda's binary can crash with SIGILL on CPUs lacking the build host's SIMD set |
| Align | MUSCLE | 5.1.0 | 5.1.0 |
| Align | MAFFT (E-INS-i) | 7.520 | 7.520 |
| ML tree (CLI) | IQ-TREE / RAxML-NG | — (hand-off script) | 3.1.4 / 2.0.4 |
| Trim | site's JS (`js/trim.js`) | — | same code, run by Node |
| Tree | FastTree, via [phylo-fasttree-wasm](https://github.com/hyphaltip/phylo-fasttree-wasm) | 2.1.11, single (default) or double precision | 2.1.11 `FastTree-browser` (default: byte-identical to the site's single precision) or bioconda's `FastTree`/`FastTreeMP` (`--fasttree double`: matches the site's double precision) |
| Root / annotate / duplications / export | site's JS (`js/tree.js`, `js/annotate.js`, `js/export.js`) | — | same code, run by Node |

The command-line tool is described in
[Command-line tool: `famfate.sh`](#command-line-tool-famfatesh). The
**CLI equivalent** blocks under each step below show the underlying commands,
for anyone who wants to run the steps by hand or adapt them.

---

## Quick start: MET12 in the browser

This walks through yeast MET12 (methylenetetrahydrofolate reductase,
UniProt [P46151](https://www.uniprot.org/uniprotkb/P46151)) from query to
annotated tree in about five minutes.

1. **Step 1:** click **Load MET12 query (search example)**. This loads a
   single query sequence; search is switched on with BLASTP selected.
2. **Step 2:** keep **UniProt Swiss-Prot** as the database, or switch to
   **UniProt reference proteomes** and add a few, e.g. *Neurospora crassa*,
   *Homo sapiens*, *Arabidopsis thaliana*, *Escherichia coli*. Keep the
   **Search E-value** at `1e-3`. Click **Run search**. The first run streams
   and indexes the database, which takes from tens of seconds to a few
   minutes depending on its size; later searches in the same session reuse
   the index.
3. Review the hit table: one row per protein, with **N HSPs ▸** to expand
   multiple local alignments. Adjust the live cutoffs, keep **Full-length
   sequences** and **Include query** ticked, then click **Use N selected for
   alignment →**.
4. **Step 3:** choose **Kalign** (fast) or **MUSCLE** (more accurate on
   divergent sets). Input = **Search hits**. Click **Align**.
5. **Step 4 (optional):** **Trim** with gap threshold 0.5 to drop columns that
   are mostly gaps.
6. **Step 5:** model **Protein (LG)**, tick **Gamma rates**, click **Build
   tree**.
7. The tree is annotated automatically from UniProt. Set **Colour by →
   Phylum**, **Label → ID + species**, tick **Domains** and **Duplications**,
   and try **Root → MAD** or **Root → Outgroup taxon… → Bacteria**.
8. Click any leaf for its details card. Use **Export ▾** to save the SVG,
   rooted Newick, annotation TSV or iTOL datasets.

---

## The pipeline, step by step

```
 Input ─► Search (optional) ─► Align ─► Trim (optional) ─► Tree ─► Root ─► Annotate ─► Duplications ─► Export
 file /   BLAST / phmmer        Kalign /  gap & conserv.    Fast-   mid/MAD/  UniProt +    species         SVG, Newick,
 paste /  vs Swiss-Prot or      MUSCLE    thresholds        Tree    outgroup/ InterPro     overlap         TSV, iTOL,
 API      proteome set                                              branch                                 metadata
```

The pipeline stepper at the top shows each step's state. Every step writes a
row to the **Run metadata** table (tool, version, parameters, input/output and
time), which can be exported as JSON.

### Step 1 · Input sequences

**Browser.** Three ways to load sequences:

- **Upload FASTA:** drag a file onto the drop zone (`.fa`, `.fasta`, `.faa`,
  `.fna`, …).
- **Paste FASTA:** paste into the text box.
- **Retrieve from API:**
  - *UniProt*, with full query syntax, e.g. `gene:MET12 AND reviewed:true` or
    `family:"MTHFR family" AND taxonomy_id:4751`.
  - *NCBI Protein / Nucleotide*: an Entrez query.
  - *InterPro*: an entry accession such as `PF02219`, which fetches member
    proteins via UniProt cross-references.

The *Limit* field caps the number of records. Example buttons load MET12 (a
single query, for the search workflow), a small protein set, or a ready-made
nucleotide alignment.

Sequence names are the first word of each FASTA header. The rest of the
header is kept: UniProt-style `OS=`, `OX=`, `GN=` fields feed the annotation
step even when offline. Sequence type (protein or nucleotide) is detected
automatically.

Two paths through the pipeline:
- **Family from a query:** load one or a few proteins, then use step 2 to
  collect homologs.
- **Family you already have:** load the full set and go straight to step 3, or
  to step 4/5 if it is already aligned.

**CLI equivalent**

```sh
# A single query
curl -s https://rest.uniprot.org/uniprotkb/P46151.fasta > query.fa

# A family by UniProt query (the same syntax as the browser box)
curl -s -G https://rest.uniprot.org/uniprotkb/search \
     --data-urlencode 'query=family:"MTHFR family" AND reviewed:true' \
     --data-urlencode 'format=fasta' --data-urlencode 'size=500' > family.fa

# Members of a Pfam/InterPro entry
curl -s -G https://rest.uniprot.org/uniprotkb/search \
     --data-urlencode 'query=xref:interpro-PF02219 AND reviewed:true' \
     --data-urlencode 'format=fasta' --data-urlencode 'size=500' > pf02219.fa

# NCBI Entrez
ids=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term=MTHFR%5BGene%5D&retmode=json&retmax=50" \
      | jq -r '.esearchresult.idlist | join(",")')
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=$ids&rettype=fasta&retmode=text" > ncbi.fa
```

(For more than 500 UniProt records, follow the `Link: …; rel="next"` header,
or use the `/uniprotkb/stream` endpoint.)

### Step 2 · Homology search (optional)

**Browser.**

1. Pick the program:
   - **BLASTP**: protein.
   - **BLASTN**: nucleotide.
   - **phmmer**: protein, HMMER. More sensitive for distant homologs; slower;
     needs the database in memory.
2. Pick the **Database source**:

   | Source | What it is | Notes |
   |---|---|---|
   | **UniProt Swiss-Prot** (default) | Reviewed UniProtKB (~570k proteins, ~94 MB gz) | Streamed from `ftp.uniprot.org`; never fully held in JS memory |
   | **UniProt reference proteomes** | Searchable catalogue of ~37,600 proteomes | Filter by species, `UP…` ID or taxid; **Add to search set** (or double-click), repeat to search many at once; optionally **+ Swiss-Prot** |
   | **PDB SEQRES** | Sequences of PDB structures | The server sends no CORS headers: download it and use *Custom*, or give a CORS-proxy prefix |
   | **Custom file / URL** | Any FASTA, plain or gzip | URLs must allow CORS; gzip is detected from the file's magic bytes |

   With several proteomes, each is indexed as its own BLAST database and they
   are searched **together** in one run. E-values then reflect the combined
   database size, as a single CLI search against all of them would. Built
   indexes are cached for the session, so adding one more proteome indexes
   only that one. Sets above ~1.5M proteins get a warning because indexing in
   the browser becomes slow.
3. Set **Max hits** (BLAST `-max_target_seqs`) and **Search E-value**
   (default `1e-3`; BLAST `-evalue`, phmmer `-E`/`--domE`). Low-complexity
   filtering is off (`-seg no`, or `-dust no` for blastn) so repeat-rich
   proteins aren't masked away.
4. **Run search.** The hit table has:
   - **One row per subject protein.** Several HSPs to one protein, or hits
     from several queries, are grouped. **N HSPs ▸** expands them, with
     subject coordinates.
   - **Subject range / length**: the part of the protein the HSPs cover, out
     of its full length.
   - **Live cutoffs** (E-value, default `1e-5`; score; % identity; alignment
     length) that filter as you type. They can only tighten the search
     E-value; to see weaker hits, loosen it and search again.
   - **Select all / none** or per-row checkboxes.
5. Choose what goes into the alignment:
   - **Full-length sequences** (default on): each protein's whole sequence,
     pulled from the same in-browser BLAST database with `blastdbcmd`.
     Proteins whose database couldn't be ID-indexed are flagged **HSP only**.
   - Unticked: only the matched region (best HSP, extended by any overlapping
     or adjacent HSPs). Useful for single-domain trees.
   - **Include query sequence(s)**: add your queries to the set.
6. Click **Use N selected for alignment →**.

**CLI equivalent**

```sh
UP=https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase

# Database: Swiss-Prot...
curl -s $UP/complete/uniprot_sprot.fasta.gz | gunzip > sprot.fa
# ...and/or reference proteomes: <Division>/<UPID>/<UPID>_<TAXID>.fasta.gz
# (look up ID, taxid and division in $UP/reference_proteomes/README)
curl -s $UP/reference_proteomes/Eukaryota/UP000001805/UP000001805_367110.fasta.gz | gunzip > ncrassa.fa
curl -s $UP/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz   | gunzip > human.fa

# Index each database, with IDs indexed so full-length subjects can be fetched
for db in sprot ncrassa human; do
  makeblastdb -in $db.fa -dbtype prot -out $db -parse_seqids
done

# Search all of them together (one combined database size, as in the browser)
blastp -query query.fa -db "sprot ncrassa human" -seg no \
       -evalue 1e-3 -max_target_seqs 50 -num_threads 8 \
       -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sacc slen' \
       > hits.tsv

# One entry per subject protein (the browser's grouping), with the live E-value cutoff
awk -F'\t' '$11 <= 1e-5' hits.tsv | cut -f13 | sort -u > subject_ids.txt

# Full-length subject sequences (the "Full-length sequences" option)
blastdbcmd -db "sprot ncrassa human" -entry_batch subject_ids.txt > hits.fa
cat query.fa hits.fa | seqkit rmdup -n > family.fa     # "Include query", deduplicated

# --- phmmer alternative (more sensitive for distant homologs)
cat sprot.fa ncrassa.fa human.fa > searchset.fa
phmmer -E 1e-3 --domE 1e-3 --noali --tblout phmmer.tbl query.fa searchset.fa > /dev/null
grep -v '^#' phmmer.tbl | awk '{print $1}' | sort -u > phmmer_ids.txt
esl-sfetch --index searchset.fa && esl-sfetch -f searchset.fa phmmer_ids.txt > hits.fa
```

### Step 3 · Multiple sequence alignment

**Browser.**
- **Tool**:
  - **Kalign**: very fast; a good default.
  - **MUSCLE 5**: slower, usually more accurate for divergent families.
  - **MAFFT E-INS-i**: MAFFT's most accurate mode for proteins with several
    domains or long unalignable regions, which suits domain-shuffled gene
    families. It aligns every pair, so it's slow: best up to ~200 sequences
    (the page asks before running more).
    - Biowasm ships MAFFT's core programs but not the `mafft` wrapper, so the
      page runs them directly with E-INS-i parameters: `tbfast`, then
      `dvtditr` refinement (16 iterations).
- **Input**: **Search hits** (from step 2) or **Original input** (step 1).
- **Align** shows a coloured alignment viewer: residue classes, a ruler,
  every column scrollable, the first 40 sequences. **Download aligned FASTA**
  saves it.

The aligner needs at least 2 sequences. If you searched but forgot **Use
selected**, the error message says so.

**CLI equivalent**

```sh
kalign -i family.fa -o family.aln.fa                  # browser: Kalign 3.3.1
muscle -align family.fa -output family.aln.fa         # browser: MUSCLE 5.1 (CLI: add -threads 8)
mafft --genafpair --ep 0 --maxiterate 16 family.fa > family.aln.fa   # browser: MAFFT E-INS-i, 16 iterations (close, not byte-identical)
# Large sets (CLI only): muscle -super5 family.fa -output family.aln.fa, or mafft --auto
```

### Step 4 · Alignment trimming

**Browser.** Columns are removed when:
- **Gap threshold**: the fraction of gaps is above this value (0.5 = drop
  columns that are more than half gaps).
- **Min conservation**: the most common residue makes up less than this
  fraction of the non-gap residues (0 = off).
- **Trim leading/trailing gaps**: also drop end columns where the first or
  last sequence has a gap. This is a crude heuristic; prefer the gap
  threshold.

The stats panel reports columns kept and removed. Trimming that would remove
every column is refused with a message. The tree step uses the most processed
data available: trimmed, then aligned, then input.

**CLI equivalent**

```sh
# trimAl's -gt is the fraction of sequences WITHOUT a gap, so gap threshold 0.5 ⇔ -gt 0.5
trimal -in family.aln.fa -out family.trim.fa -gt 0.5
# Automated alternatives widely used in published studies: trimal -automated1, or ClipKIT
```

An exact re-implementation of the browser's rule (gap fraction ≤ G and
majority-residue fraction ≥ C):

```python
# trim.py: python3 trim.py family.aln.fa 0.5 0.0 > family.trim.fa
import sys
names, seqs = [], []
for line in open(sys.argv[1]):
    line = line.strip()
    if line.startswith('>'): names.append(line[1:].split()[0]); seqs.append('')
    elif line: seqs[-1] += line
G, C = float(sys.argv[2]), float(sys.argv[3])
keep = []
for i in range(len(seqs[0])):
    col = [s[i].upper() for s in seqs]
    res = [c for c in col if c not in '-.?~']
    gap = 1 - len(res) / len(col)
    cons = max(res.count(c) for c in set(res)) / len(res) if res else 0
    if gap <= G and cons >= C: keep.append(i)
for n, s in zip(names, seqs):
    print(f'>{n}\n' + ''.join(s[i] for i in keep))
```

### Step 5 · Tree inference

**Browser.**
- **Model**:
  - *Auto*: detects protein or nucleotide.
  - *Nucleotide*: JC by default; tick **GTR** for GTR.
  - *Protein (JTT)*, *Protein (LG)* or *Protein (WAG)*.
- **Gamma rates**: `-gamma`, a discrete-gamma-rescaled likelihood. Use it when
  you want branch lengths and likelihoods comparable across runs.
- **Precision**: *Single* (FastTree's default) or *Double* (`-DUSE_DOUBLE`,
  the build bioconda ships).
  - Double allows much shorter minimum branch lengths, so choose it for very
    similar sequences (e.g. strains of one species, recent paralogs).
  - Nucleotide trees with near-zero branches can differ between the two;
    protein trees usually differ only in the 4th decimal place.
  - In WebAssembly both run at the same speed (single precision is faster
    natively only because of SSE code, which the WASM build doesn't use).
    Double needs more memory, and its 307 KB module loads only when selected.
- **No support values**: `-nosupport`. Otherwise you get FastTree's SH-like
  local support (0–1, from 1000 resamples).
- **Build tree** runs FastTree in a background worker. A timer and a
  **Cancel** button appear; the page stays usable.

The viewer offers ladderize, cladogram layout, support values, and **Find
tip / species**. Weak supports are greyed: below 0.7 (or 70) for single
values, and below SH-aLRT 80 / UFBoot 95 for IQ-TREE's paired values. Its
rooting, colouring, domain and duplication controls are described in steps
6–8.

**Maximum-likelihood trees (IQ-TREE, RAxML-NG).** These are too heavy to run
in a browser tab, so the tree step hands off to your own machine:
1. **Run IQ-TREE / RAxML-NG on your machine… → Download run script** saves
   `run_famfate_<method>.sh`. It contains the alignment the tree step uses
   (trimmed → aligned → input) and chooses a runner automatically:
   - a local `famfate.sh`;
   - Apptainer/Singularity or Docker with the `phylo-fam-fate-cli` image;
   - otherwise IQ-TREE or RAxML-NG directly, if installed.
2. Run it (`bash run_famfate_iqtree.sh 8`) on a laptop or cluster.
3. **Load a tree (Newick)…** with the resulting `famfate_<method>/family.nwk`.

The method choices are:
- *IQ-TREE 3*: ModelFinder plus ultrafast bootstrap and SH-aLRT;
- *RAxML-NG*: ML plus bootstrap, with TBE support;
- *all three*, with IQ-TREE as the main tree.

**Load a tree** accepts any Newick: IQ-TREE `.treefile`/`.contree`,
RAxML-NG `.support`, or trees from other programs. A loaded tree gets the
same rooting, annotation, domain and duplication views. Paired support
labels such as `95.3/98` are kept through rerooting and export. Leaf names
should match the sequence names, so annotation and descriptions line up.

**CLI equivalent**

```sh
FastTree -lg -gamma family.trim.fa > family.nwk               # protein, LG+CAT with Gamma20 rescaling
FastTree -wag family.trim.fa > family.nwk                     # protein, WAG
FastTree family.trim.fa > family.nwk                          # protein, JTT (FastTree default)
FastTree -nt -gtr -gamma family.trim.fa > family.nwk          # nucleotide GTR
# Multithreaded on the CLI: FastTreeMP (same options)

# Maximum likelihood (what famfate.sh --tree iqtree / raxml-ng runs)
iqtree3 -s family.trim.fa -m MFP -B 1000 -alrt 1000 -T AUTO --seed 12345 --prefix iqtree/family
raxml-ng --all --msa family.trim.fa --model LG+G4 --bs-trees 'autoMRE{1000}' --bs-metric tbe --seed 12345 --prefix raxml/family
```

To reproduce the browser's trees on the CLI, match its precision:
- *Single* (the default): use a single-precision, strict-floating-point
  build (`cli/build-fasttree-browser.sh`, or `FastTree-browser` in the
  container). Trees are byte-identical.
- *Double*: use bioconda's `fasttree`. Same topology, numbers within 1e-5.

The exact arguments and precision are recorded in the metadata table.

### Step 6 · Rooting

FastTree trees are **unrooted**. The **Root** menu reroots the original tree.
Switching methods is lossless, and support values stay with their branch.

| Option | What it does | When to use |
|---|---|---|
| **Midpoint** (default) | Root halfway along the longest tip-to-tip path | Quick default; assumes roughly equal rates |
| **MAD** | Minimal Ancestor Deviation ([Tria et al. 2017](https://doi.org/10.1038/s41559-017-0193)): the root that makes the tree most clock-like across all tip pairs | Better than midpoint with uneven rates; no outgroup needed. The legend shows the score and *ambiguity* (near 1 = another root is almost as good). Up to 900 tips |
| **Outgroup taxon…** | Root on the branch separating a chosen taxon (any rank from the annotations, e.g. *Bacteria*, *Microsporidia*) | You know which lineage branches off first. If it isn't monophyletic, the tree is rooted on the best split and the legend warns how many tips are on the wrong side |
| **Chosen branch** | Click any branch → **Root on this branch** | Manual rooting on a known duplication or a clade |
| **As inferred (unrooted)** | FastTree's arbitrary trifurcation | To see the raw output |

**CLI equivalent**

```sh
# Midpoint (ETE3)
python3 -c "
from ete3 import Tree
t = Tree('family.nwk'); t.set_outgroup(t.get_midpoint_outgroup()); t.write(outfile='family.mid.nwk')"

# Outgroup: one tip, or several (rooted on their common ancestor)
nw_reroot family.nwk 'sp|P0AEZ1|METF_ECOLI' > family.og.nwk          # newick_utils
python3 -c "
from ete3 import Tree
t = Tree('family.nwk'); og = [l for l in t if 'ECOLI' in l.name or 'BACSU' in l.name]
t.set_outgroup(t.get_common_ancestor(og) if len(og) > 1 else og[0]); t.write(outfile='family.og.nwk')"

# MAD: the authors' program (Tria, Landan & Dagan 2017), e.g.
python3 mad.py family.nwk          # writes the MAD-rooted tree next to the input
```

### Step 7 · Annotation: taxonomy, domains, function

**Browser.** After a tree is built, leaves are annotated automatically when
their names look like UniProt records: `sp|P46151|MET12_YEAST`,
`tr|Q7S830|Q7S830_NEUCR`, or a bare accession such as `Q7S830`.
**Annotate from UniProt** reruns it.

- **Source**: one batched request per 400 leaves to UniProt's
  `uniprotkb/accessions` endpoint returns organism, taxid, ranked lineage,
  protein and gene names, length, Pfam, domain features, function and
  annotation score. Leaves that aren't UniProt fall back to `OS=`/`OX=`/`GN=`
  in their FASTA header, which needs no network and gives species only, no
  lineage.
- **Label**: ID, species, ID + species, gene (species), or protein name.
- **Colour by**: domain of life, kingdom, phylum, subphylum, class, order,
  family, genus or species.
  - Each leaf gets a colour swatch.
  - A **branch is coloured only when its whole clade shares that taxon**, so a
    taxon that isn't a single clade looks patchy. That pattern is a first hint
    of paralogy or horizontal transfer.
  - The legend lists groups with counts. The 12 most frequent get distinct
    colours; the rest are grey.
- **Click a leaf** for a card with protein and gene names, species and
  lineage, a to-scale **domain diagram**, function text, annotation score and
  a UniProt link.
- **Click a branch** for the clade's size, support and species, with
  **Root on this branch**, **Highlight clade** and **Copy tip names**.
- **Domains** draws the Pfam domain architecture for every leaf beside the
  tree, on a common scale.
  - Coordinates come from InterPro, one request per protein, with retries
    when InterPro is busy; up to 600 proteins. UniProt's own domain features
    are the fallback.
  - Domain fusions, losses and truncated copies stand out against the tree.
  - The track always draws the full UniProt protein, even if you aligned only
    the HSP region.

**CLI equivalent**

```sh
# Leaf labels → accessions (handles sp|ACC|NAME and tr|ACC|NAME labels)
grep -oE '(sp|tr)\|[A-Z0-9]+\|' family.nwk | cut -d'|' -f2 | sort -u > accs.txt

# Batched annotation (≤400 accessions per request), the same fields as the browser
curl -s -G https://rest.uniprot.org/uniprotkb/accessions \
  --data-urlencode "accessions=$(paste -sd, accs.txt)" \
  --data-urlencode 'fields=accession,id,organism_name,organism_id,lineage,protein_name,gene_names,length,xref_pfam,ft_domain,cc_function,annotation_score' \
  --data-urlencode 'format=tsv' > annotations.tsv

# Pfam domain coordinates per protein (InterPro)
while read acc; do
  curl -s "https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/$acc/?page_size=100" \
  | jq -r --arg a "$acc" '.results[]? | .metadata as $m | .proteins[].entry_protein_locations[].fragments[]
        | [$a, $m.accession, $m.name, .start, .end] | @tsv'
  sleep 0.2
done < accs.txt > pfam_domains.tsv
# Fully local alternative: hmmscan --domtblout doms.tbl Pfam-A.hmm family.fa
```

### Step 8 · Duplications vs speciations

**Browser.** Tick **Duplications** on a rooted, annotated tree. Each internal
node is labelled with the **species-overlap** rule:
- a **duplication** (red square) when two of its child clades share a species;
- otherwise a **speciation** (blue dot).

The legend counts each, and clicking a node shows the event in its branch
card. Results depend on the root, so compare Midpoint, MAD and a sensible
outgroup. Species are identified by NCBI taxid when available, otherwise by
species name.

**CLI equivalent** (ETE3 implements the same species-overlap algorithm):

```python
from ete3 import PhyloTree
import csv
sp = {r['Entry']: r['Organism (ID)'] for r in csv.DictReader(open('annotations.tsv'), delimiter='\t')}
acc = lambda name: name.split('|')[1] if '|' in name else name
t = PhyloTree('family.mid.nwk', sp_naming_function=lambda name: sp.get(acc(name), name))
events = t.get_descendant_evol_events()          # species-overlap
print('duplications', sum(e.etype == 'D' for e in events),
      'speciations',  sum(e.etype == 'S' for e in events))
```

For reconciliation against a species tree (losses, transfers), use
GeneRax or NOTUNG on the CLI.

### Step 9 · Export and provenance

**Browser.** In the tree **Export ▾** menu:

| Export | Contents |
|---|---|
| **SVG (as shown)** | The tree exactly as displayed (colours, domains, markers); styles inlined |
| **Newick (FastTree, unrooted)** | FastTree's raw output |
| **Newick (rooted, as shown)** | After the current rooting and ladderizing, with support values |
| **Copy Newick** | Raw Newick to the clipboard |
| **Annotations (TSV)** | One row per leaf: accession, species, taxid, domain…genus ranks, protein, genes, length, Pfam, domain coordinates, function, score, source |
| **iTOL colour strip** | `DATASET_COLORSTRIP` for the current *Colour by* rank, with legend |
| **iTOL protein domains** | `DATASET_DOMAINS` with Pfam domains to scale |
| **Duplication / speciation nodes (TSV)** | Each labelled node: event, support, number of tips, species, tip names |

To use iTOL, upload the Newick at <https://itol.embl.de>, then drag the
dataset files onto the tree; leaf IDs match the Newick labels.

The alignment, trimmed alignment and tree also have their own download
buttons. **Run metadata → Export JSON** records every step's tool, version,
parameters, input/output sizes and timing, for a methods section or to rerun
the analysis with the CLI commands above.

---

## Command-line tool: `famfate.sh`

`cli/famfate.sh` runs the whole pipeline in one command, with the browser's
defaults:

1. Search (BLAST+ or phmmer).
2. Full-length hits via `blastdbcmd`.
3. Align with Kalign or MUSCLE.
4. Trim.
5. FastTree.
6. Root.
7. Annotate from UniProt and InterPro.
8. Label duplications.
9. Export.

Steps 4 and 6–9 run the **website's own JavaScript** under Node
(`cli/famfate-js.mjs` imports `js/trim.js`, `js/tree.js`, `js/annotate.js`
and `js/export.js`). They produce exactly what the browser would, including
MAD rooting and the same SVG renderer. The search, alignment and tree steps
use the native tools at the versions in the table above.

**Two FastTree builds.** FastTree can be compiled in single or double
precision, and the two give different trees on alignments with very short
branches. In CI, on the bundled 24-sequence nucleotide example, the
double-precision build gave a different topology (one split) and branch
lengths up to 0.09 apart with GTR+Γ; protein trees agreed to within 3e-4.
- **`--fasttree browser`** (default) uses `FastTree-browser`: upstream
  FastTree 2.1.11 compiled the way the website's WASM build is (single
  precision, no SSE, strict floating point `-ffp-contract=off`), by
  `cli/build-fasttree-browser.sh`. Its trees are **byte-identical to the
  website's**, which CI checks.
- **`--fasttree double`** uses bioconda's `FastTree`/`FastTreeMP`, the
  double-precision build. It is multithreaded and matches the website's
  **Precision → Double** option (CI checks that they give the same topology,
  with numbers within 1e-5; the small differences come from different math
  libraries).

### Get it

**Container** (recommended; linux/amd64), built and tested by CI
(`.github/workflows/cli-container.yml`) and published to GitHub Container
Registry:

| Tag | Built from |
|---|---|
| `:edge` | the latest `main` |
| `:sha-<commit>` | every published build |
| `:X.Y.Z`, `:X.Y`, `:latest` | release tags `vX.Y.Z` |

```sh
# Docker: mount the working directory at /data and run as yourself so outputs are yours
docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":/data \
    ghcr.io/hyphaltip/phylo-fam-fate-cli:edge famfate.sh -q P46151

# Apptainer / Singularity (HPC): pull once into a .sif, then exec
apptainer pull famfate.sif docker://ghcr.io/hyphaltip/phylo-fam-fate-cli:edge
apptainer exec famfate.sif famfate.sh -q P46151 -T 8
```

Under Apptainer the current directory and `$HOME` are bound automatically.
Downloaded and indexed databases are cached in `--dbdir` (default
`./famfate_db`, or `$FAMFATE_DB`). On a cluster, point that at shared scratch
so every job reuses one Swiss-Prot or proteome index:

```sh
export FAMFATE_DB=/scratch/$USER/famfate_db
apptainer exec famfate.sif famfate.sh -q P46151 -p UP000001805,UP000005640 -T $SLURM_CPUS_PER_TASK
```

**Conda** (no container; linux-64 or macOS):

```sh
git clone https://github.com/hyphaltip/phylo-fam-fate && cd phylo-fam-fate
mamba env create -f environment.yml && mamba activate famfate
./cli/build-fasttree-browser.sh        # compiles FastTree-browser into $CONDA_PREFIX/bin (needs gcc)
./cli/famfate.sh -h
```

**Build the image yourself:**
`docker build --platform linux/amd64 -t phylo-fam-fate-cli .`. BLAST+ 2.17.0
has no linux-arm64 build on bioconda, so ARM machines either run the amd64
image under emulation or use conda on macOS.

### Usage

```text
famfate.sh (-q ACC | -i QUERY.fa | -f FAMILY.fa | -a ALIGNED.fa) [options]

Input (one of):
  -q, --query ACC[,ACC…]   UniProt accession(s) of the query protein(s)
  -i, --query-fasta FILE   query FASTA
  -f, --family FILE        unaligned family FASTA: skip the search
  -a, --aligned FILE       aligned FASTA: skip search and alignment
Search:      -p UPID,UPID…  --no-sprot  --db FILE[,FILE]  -t blastp|blastn|phmmer
             -e 1e-3 (search E)  -E 1e-5 (keep-hit E)  -n 50  --hsp-only  --no-query
Align:       -A kalign|muscle|mafft  -g 0.5|none  -c 0
Tree:        --tree fasttree|iqtree|raxml-ng[,…]  (first = main tree; others compared)
             FastTree: -m lg|wag|jtt|nt|gtr  --no-gamma  --fasttree browser|double  --fasttree-args "…"
             IQ-TREE / RAxML-NG: --ml-model MFP|<model>  --raxml-model M  --bootstrap 1000  --seed 12345
                                 --iqtree-args "…"  --raxml-args "…"
Post:        -r midpoint|mad|none|outgroup:RANK:NAME  --rank phylum  --no-annotate  --no-domains
General:     -o OUTDIR  -d DBDIR  -T THREADS  -h  -V
```

The defaults match the browser:
- Search: Swiss-Prot, BLASTP, search E-value 1e-3, hits kept at E ≤ 1e-5,
  50 hits, full-length subjects, query included.
- Alignment and tree: Kalign; trimming at gap fraction 0.5; FastTree LG with
  Gamma (`-lg -gamma`).
- Post-processing: midpoint root; colours and iTOL at phylum.

**ML trees.** `--tree iqtree` runs IQ-TREE 3:
- ModelFinder (`-m MFP`) chooses the substitution model;
- support is ultrafast bootstrap (`-B 1000`) plus SH-aLRT (`-alrt 1000`),
  written as `SH-aLRT/UFBoot` (e.g. `95.3/98`).

`--tree raxml-ng` runs `raxml-ng --all`:
- ML search, then bootstrap until converged (`autoMRE{N}`), with transfer
  bootstrap (TBE) support;
- RAxML-NG has no ModelFinder. Its model is `--raxml-model` if given, else
  `--ml-model` if not MFP, else IQ-TREE's chosen model when IQ-TREE ran in
  the same command and RAxML-NG accepts it, else LG+G4 (protein) / GTR+G4
  (nucleotide).

Several methods, e.g. `--tree iqtree,raxml-ng,fasttree`:
- the first is the **main tree**, which gets rooted, annotated and labelled
  with duplications;
- every tree is saved as `family.<method>.nwk`;
- `trees_compare.tsv` lists, for each other method against the main tree:
  same topology?, splits found in only one of the two, the Robinson–Foulds
  distance, and the largest branch-length difference on shared splits.

Duplication calls on splits that the methods disagree on deserve less trust
than those on splits all three recover.

### Examples

```sh
# MET12 against Swiss-Prot (the browser's quick-start)
famfate.sh -q P46151

# Against a set of reference proteomes instead, MUSCLE, MAD rooting
famfate.sh -q P46151 -p UP000001805,UP000005640,UP000006548,UP000000625 --no-sprot \
           -A muscle -r mad -o met12_proteomes

# More sensitive search; root on bacteria; colour by class
famfate.sh -q P46151 -t phmmer -r outgroup:domain:Bacteria --rank class

# You already have the family: align → tree → annotate
famfate.sh -f my_family.fa -m wag -o myfam

# You already have an alignment (e.g. from MAFFT): tree + annotation only, offline
famfate.sh -a my_family.aln.fa --no-annotate -r midpoint

# Domain-only tree: matched regions, no trimming
famfate.sh -q P46151 --hsp-only -g none

# Publication-style: MAFFT E-INS-i, IQ-TREE main tree, compared with RAxML-NG and FastTree
famfate.sh -q P46151 -p UP000001805,UP000005640,UP000000625 -A mafft \
           --tree iqtree,raxml-ng,fasttree -r outgroup:domain:Bacteria -T 16

# RAxML-NG only, with a chosen model
famfate.sh -a family.aln.fa --tree raxml-ng --raxml-model LG+I+G4 --bootstrap 500
```

### Outputs (in `-o OUTDIR`)

| File | Contents |
|---|---|
| `query.fa` | Query sequence(s) |
| `hits.tsv` | BLAST tabular output, one row per HSP (`sseq`, `sacc`, `slen` columns added) |
| `best_hsp.tsv`, `subject_ids.txt` | One row / ID per subject protein passing `-E` (the browser's grouped hit table) |
| `hits.fa`, `hits.summary` | Subject sequences: full-length from `blastdbcmd`, or HSP region where unavailable |
| `phmmer.tbl` | phmmer per-sequence table (with `-t phmmer`) |
| `family.fa` | Sequences that were aligned (query + hits, deduplicated by name) |
| `family.aln.fa`, `family.trim.fa` | Alignment and trimmed alignment |
| `family.nwk` | Main tree, unrooted (first `--tree` method) |
| `family.<method>.nwk` | Each tree built (`fasttree`, `iqtree`, `raxml-ng`) |
| `iqtree/`, `raxml/` | Full IQ-TREE / RAxML-NG output (`.iqtree` report with the chosen model, `.log`, bootstrap trees) |
| `trees_compare.tsv` | With several `--tree` methods: topology agreement and RF distance vs the main tree |
| `family.rooted.nwk` | Rooted, ladderized tree (the browser's "Newick (rooted, as shown)") |
| `family.svg` | Tree drawn by the browser's renderer: taxon colours, domain track, duplication markers |
| `family.annotations.tsv` | One row per leaf: ranks, protein/gene, length, Pfam, domain coordinates, function |
| `family.itol_<rank>.txt`, `family.itol_domains.txt` | iTOL colour-strip and protein-domain datasets |
| `family.events.tsv` | Duplication / speciation nodes, with support, tips and species |
| `family.post.json` | Rooting method and note (e.g. MAD score, outgroup monophyly), counts |
| `run.json`, `run.log` | Provenance (command, each step's tool, version, parameters, seconds) and full log |

To view or refine the result interactively, open the site, load
`family.trim.fa` in step 1 and build the tree in step 5.

### `famfate-js`: the website's code as standalone commands

The container also puts `famfate-js` on the PATH:

```sh
famfate-js trim ALN.fa --gap 0.5 --cons 0 > trimmed.fa         # browser trimming rule
famfate-js post --tree T.nwk --fasta FAMILY.fa --root mad --out T # root + annotate + duplications + SVG/TSV/iTOL
famfate-js fasttree-wasm ALN.fa -lg -gamma > tree.nwk             # the browser's FastTree WASM build (add --double for double precision)
famfate-js compare-trees A.nwk B.nwk --tol 1e-4                  # same unrooted topology? max branch-length difference
famfate-js seqtype FASTA                                         # nt | aa (the browser's detection rule)
```

### What CI checks

On every push touching the CLI, JS or container (and on pull requests), the
workflow builds the image and runs `cli/tests/ci.sh` inside it:

1. **Offline pipeline:** family → Kalign → trim → FastTree, and MUSCLE/WAG
   without trimming.
2. **FastTree parity**, on nucleotide (JC, GTR+Γ) and protein (JTT, LG+Γ,
   WAG) alignments:
   - `FastTree-browser` and the website's single-precision WASM must give
     **byte-identical** Newick;
   - bioconda's double-precision `FastTree` and the website's double WASM
     must give the same topology, with numbers within 1e-5;
   - `famfate.sh`'s own tree must equal the website's.
3. **Live search:** MET12 against the *S. cerevisiae* and *E. coli*
   proteomes. It checks that:
   - all hits come back full-length;
   - MET13 and MetF are recovered;
   - UniProt taxonomy is filled in;
   - the MET12/MET13 duplication is labelled;
   - iTOL files are written;
   - MAD rooting is applied.
4. **phmmer** with outgroup rooting.
5. **ML trees:** IQ-TREE + RAxML-NG + FastTree on one alignment. It checks
   that all three trees are written, that IQ-TREE's is the main tree, that
   the `SH-aLRT/UFBoot` labels survive rooting, that `trees_compare.tsv` is
   filled in, and that the chosen model is recorded. Also MAFFT E-INS-i.
6. **Option errors.**
7. **Apptainer:** the image is converted to a SIF and run with
   `apptainer exec`.

The image is pushed only if all of these pass.

---

## Browser vs command line: choosing

- **Use the browser** for exploring a family: which proteomes have it, how
  many copies, which domain architectures. Also for teaching, and for
  checking a hypothesis before a larger run. Comfortable sizes:
  - ≤ 500 sequences to align;
  - ≤ ~1.5M proteins in the search set;
  - trees up to a few thousand tips (MAD ≤ 900 tips; domain track ≤ 600
    proteins).
- **Use the CLI** for large searches (all reference proteomes, nr), large
  alignments (MUSCLE `-super5`, MAFFT), ML with bootstrap (IQ-TREE,
  RAxML-NG), reconciliation against a species tree (GeneRax, NOTUNG), or
  anything to be rerun automatically.
- **Mix them.** Run `famfate.sh` (or just the heavy search and alignment) on
  a cluster, then load `family.trim.fa` into the browser for the interactive
  tree, annotation and duplication view.
  Or explore in the browser and copy the parameters from the exported
  metadata JSON into the CLI script.

---

## Data sources and limits

### In the browser
- All compute runs in WebAssembly, with a practical memory ceiling of about
  2 GB per tab.
- Databases are streamed into Biowasm's virtual filesystem (`CLI.mount`), not
  copied into JavaScript memory. BLAST uses v4 databases because v5 needs
  LMDB, which is unavailable in WASM.
- phmmer can't stream, so its database is downloaded, decompressed in memory
  and cached for later runs.
- Nothing is uploaded. Only public APIs are contacted: UniProt, InterPro,
  NCBI and the UniProt FTP.

### Remote APIs (CORS-enabled, called directly from the browser)
| API | Use | Rate limit |
|---|---|---|
| [UniProt REST](https://rest.uniprot.org) | Sequence search and FASTA retrieval | ~10 req/s |
| [UniProt `uniprotkb/accessions`](https://rest.uniprot.org/docs/) | Batched leaf annotation (≤ 400 accessions per request) | ~10 req/s |
| [UniProt FTP](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/) | Swiss-Prot and reference proteomes | — |
| [InterPro](https://www.ebi.ac.uk/interpro/api/) | Family retrieval; Pfam domain coordinates | Fair use; busy responses (503) are retried |
| [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/) | Nucleotide/protein retrieval | 3 req/s without a key |

NCBI BLAST's URL API and EBI's HMMER web service send no CORS headers, so
remote BLAST/HMMER would need a proxy; the in-browser search is used instead.

---

## Troubleshooting

A **Debug log** panel (bottom-right; a red badge counts errors) captures
console output, uncaught errors, tool logs, FastTree arguments and step
timings. **Copy log** copies it for a bug report. Step failures also appear
inline under the step. After a new release, hard-refresh (Ctrl/Cmd+Shift+R)
to drop cached scripts.

| Symptom | Cause / fix |
|---|---|
| `makeblastdb … SetFullVersion()` line in the log | A harmless warning from the WASM build; the database built if the log says `added N sequences` |
| Aligner: "Only 1 sequence" / need ≥ 2 | You searched but didn't click **Use N selected for alignment →**, or Input is set to *Original input* |
| A hit marked **HSP only** | Its database's headers couldn't be ID-indexed, so only the matched region is available |
| FastTree: "needs an alignment" | Sequences have different lengths; run step 3 first |
| No colours / outgroup list says "annotate first" | Leaves aren't UniProt accessions and headers have no `OS=`; use UniProt-derived sequences, then **Annotate from UniProt** |
| Domain track sparse | InterPro was busy or the protein has no Pfam match; toggle **Domains** off and on to retry |
| "found no hits at E-value ≤ …" | Loosen **Search E-value** or add databases |
| Tree step stops responding after many failed runs | Reload the page; FastTree's error path is guarded, but WASM memory can still be exhausted after hundreds of failures |

---

## How it works (architecture)

```
index.html ── js/app.js              UI, state, pipeline wiring
               ├─ js/api.js          FASTA parse/format; UniProt / NCBI / InterPro retrieval
               ├─ js/tools.js        Biowasm/Aioli wrappers: BLAST (+makeblastdb, blastdbcmd), phmmer, Kalign, MUSCLE
               ├─ js/trim.js         column trimming (typed-array column stats)
               ├─ js/msa.js          coloured alignment viewer
               ├─ js/fasttree-worker.js  Web Worker running FastTree WASM (stack-pointer guard on errors)
               ├─ js/tree.js         Newick parser, SVG renderer, midpoint/MAD/outgroup rooting, species-overlap events
               ├─ js/annotate.js     UniProt batch annotation, lineage/domain parsing, InterPro Pfam, palettes
               └─ js/export.js       Newick / annotation TSV / iTOL / events writers (shared with the CLI)
cli/famfate.sh       command-line pipeline (BLAST+/HMMER/Kalign/MUSCLE/FastTree + famfate-js)
cli/famfate-js.mjs   runs js/*.js under Node: trim, post (root/annotate/events/SVG), fasttree-wasm, compare-trees
cli/tests/ci.sh      container tests run by CI
Dockerfile, environment.yml   pinned CLI toolchain (container / conda)
assets/fasttree/   FastTree 2.1.11 compiled to wasm32 (from phylo-fasttree-wasm), GPL-2.0
assets/examples/   example inputs (MET12 query, protein set, nucleotide alignment)
```

- Biowasm tools are loaded on demand from the Biowasm CDN into one shared
  virtual filesystem.
- FastTree is a separate, directly linked WASM module (no Emscripten), built
  from [phylo-fasttree-wasm](https://github.com/hyphaltip/phylo-fasttree-wasm).
- Each rooting stores the chosen branch as a split of leaf names rather than
  a pointer to a node, so it survives re-rendering and works on any
  orientation of the same tree.

---

## Run locally

```sh
git clone https://github.com/hyphaltip/phylo-fam-fate
cd phylo-fam-fate
python3 -m http.server 8000      # any static server; ES modules need http://, not file://
# open http://localhost:8000/
```

There is no build step: the site uses ES modules and CDN-hosted tools.
`?autotest=1` runs a tiny end-to-end pipeline on load, as a smoke test.

## Releasing

1. Bump the version in `js/version.js` and `cli/famfate.sh` (`VERSION=`); CI
   fails if they differ. The site shows it in the header and footer.
2. Add a section to `CHANGELOG.md`.
3. Tag and push: `git tag -a vX.Y.Z -m "vX.Y.Z" && git push origin vX.Y.Z`.
   CI publishes the container as `:X.Y.Z`, `:X.Y` and `:latest`.
4. Create the GitHub release from the changelog section.

## Deploy

GitHub Actions (`.github/workflows/deploy.yml`) deploys `main` to GitHub
Pages.

## License

The site code is MIT. The bundled FastTree WASM is GPL-2.0-or-later (see
`assets/fasttree/LICENSE`). Biowasm tools carry their own licenses.

If you use this in published work, cite the underlying tools: BLAST+
(Camacho et al. 2009), HMMER (Eddy 2011), Kalign 3 (Lassmann 2020),
MUSCLE 5 (Edgar 2022), FastTree 2 (Price et al. 2010), MAD rooting (Tria et
al. 2017), UniProt and InterPro (their current Nucleic Acids Research
database-issue papers).
