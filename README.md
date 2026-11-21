# phylo-fam-fate

A browser-based workbench for gene-family phylogenetics. You can start from a
single protein, find its homologs in Swiss-Prot or any set of UniProt reference
proteomes, then align, trim and infer a tree. The tree can be annotated with
taxonomy, Pfam domains and function, rooted several ways, and marked with
gene duplications. Everything runs **in your browser** (WebAssembly), with no
server, account or upload.

**Live site:** <https://hyphaltip.github.io/phylo-fam-fate/>

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
4. [Complete command-line pipeline](#complete-command-line-pipeline)
5. [Browser vs command line: choosing](#browser-vs-command-line-choosing)
6. [Data sources and limits](#data-sources-and-limits)
7. [Troubleshooting](#troubleshooting)
8. [How it works (architecture)](#how-it-works-architecture)
9. [Run locally / deploy / license](#run-locally)

---

## Two ways to run it

| | **Browser** (this site) | **Command line** (equivalent tools) |
|---|---|---|
| Install | Nothing; open the page | conda/mamba environment (below) |
| Compute | Your machine, in WebAssembly (single-threaded) | Native, multi-threaded |
| Data | Streamed from UniProt into the browser; never uploaded | Downloaded to disk |
| Best for | Exploring one family; teaching; quick looks; sharing a method | Large families (>1–2k sequences), many proteomes, batch runs, publication runs |
| Output | SVG, Newick, FASTA, TSV, iTOL files, run-metadata JSON | Same formats |

Tool versions used in the browser (use the same on the CLI to match results):

| Step | Tool | Browser version |
|---|---|---|
| Search | BLAST+ `blastp` / `blastn`, `makeblastdb`, `blastdbcmd` | 2.17.0 |
| Search | HMMER `phmmer` | 3.4 |
| Align | Kalign | 3.3.1 |
| Align | MUSCLE | 5.1.0 |
| Trim | built-in (JS); `trimal` is the CLI analogue | — |
| Tree | FastTree, via [phylo-fasttree-wasm](https://github.com/hyphaltip/phylo-fasttree-wasm) | 2.1.11 |

CLI environment, one time:

```sh
mamba create -n famfate -c conda-forge -c bioconda \
    blast hmmer kalign3 muscle=5.1 trimal fasttree \
    seqkit ete3 newick_utils curl jq
mamba activate famfate
```

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
- **Tool**: **Kalign** (very fast; a good default) or **MUSCLE 5** (slower,
  usually more accurate for divergent families).
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
# Large sets (CLI only): muscle -super5 family.fa -output family.aln.fa
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
- **No support values**: `-nosupport`. Otherwise you get FastTree's SH-like
  local support (0–1, from 1000 resamples).
- **Build tree** runs FastTree in a background worker. A timer and a
  **Cancel** button appear; the page stays usable.

The viewer offers ladderize, cladogram layout, support values (weak ones
under 0.7 in grey), and **Find tip / species**. Its rooting, colouring,
domain and duplication controls are described in steps 6–8.

**CLI equivalent**

```sh
FastTree -lg -gamma family.trim.fa > family.nwk               # protein, LG+CAT with Gamma20 rescaling
FastTree -wag family.trim.fa > family.nwk                     # protein, WAG
FastTree family.trim.fa > family.nwk                          # protein, JTT (FastTree default)
FastTree -nt -gtr -gamma family.trim.fa > family.nwk          # nucleotide GTR
# Multithreaded on the CLI: FastTreeMP (same options)
# For publication, consider ML with ultrafast bootstrap: iqtree2 -s family.trim.fa -m MFP -B 1000 -T AUTO
```

The browser compiles FastTree 2.1.11's C source to WebAssembly unchanged, so
a CLI FastTree 2.1.11 run with the same options and input should give the
same tree. The exact arguments are recorded in the metadata table.

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

## Complete command-line pipeline

A script equivalent to the default browser run (BLASTP vs Swiss-Prot plus
proteomes, full-length subjects, Kalign, trimming at 0.5 gaps, FastTree LG+Γ,
annotation):

```sh
#!/usr/bin/env bash
# famfate.sh  QUERY_ACC  "Division/UPID/UPID_TAXID ..."
set -euo pipefail
Q=${1:-P46151}                     # UniProt accession of the query
PROTEOMES=${2:-"Eukaryota/UP000001805/UP000001805_367110 Eukaryota/UP000005640/UP000005640_9606"}
EVALUE=1e-3; MAXHITS=50; THREADS=8
UP=https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase

curl -s https://rest.uniprot.org/uniprotkb/$Q.fasta > query.fa

dbs=(sprot)
[ -s sprot.pin ] || { curl -s $UP/complete/uniprot_sprot.fasta.gz | gunzip > sprot.fa
                      makeblastdb -in sprot.fa -dbtype prot -out sprot -parse_seqids; }
for p in $PROTEOMES; do
  n=$(basename $p)
  [ -s $n.pin ] || { curl -s $UP/reference_proteomes/$p.fasta.gz | gunzip > $n.fa
                     makeblastdb -in $n.fa -dbtype prot -out $n -parse_seqids; }
  dbs+=($n)
done

blastp -query query.fa -db "${dbs[*]}" -seg no -evalue $EVALUE -max_target_seqs $MAXHITS \
       -num_threads $THREADS \
       -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sacc slen' \
       > hits.tsv
cut -f13 hits.tsv | sort -u > ids.txt
blastdbcmd -db "${dbs[*]}" -entry_batch ids.txt > hits.fa
cat query.fa hits.fa | seqkit rmdup -n > family.fa

kalign -i family.fa -o family.aln.fa
trimal -in family.aln.fa -out family.trim.fa -gt 0.5
FastTree -lg -gamma family.trim.fa > family.nwk

curl -s -G https://rest.uniprot.org/uniprotkb/accessions \
  --data-urlencode "accessions=$( (echo $Q; cat ids.txt) | sort -u | paste -sd,)" \
  --data-urlencode 'fields=accession,organism_name,organism_id,lineage,protein_name,gene_names,length,xref_pfam' \
  --data-urlencode 'format=tsv' > annotations.tsv

echo "Done: family.nwk ($(grep -c '>' family.fa) sequences), annotations.tsv"
```

The outputs (`family.nwk`, `annotations.tsv`) can be opened in iTOL, or
analysed with ETE3 as in steps 6–8. Alternatively, load `family.trim.fa` into
the browser (step 1, then step 5) to use its interactive viewer (rooting,
colours, domains, duplications) on a CLI-made alignment.

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
- **Mix them.** Run the heavy search and alignment on the CLI, then load the
  alignment into the browser for the tree, annotation and duplication view.
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
               └─ js/annotate.js     UniProt batch annotation, lineage/domain parsing, InterPro Pfam, palettes
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
