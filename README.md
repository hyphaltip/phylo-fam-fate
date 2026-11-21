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
2. **UniProt reference proteome** — loads the
   [`reference_proteomes/README`](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README)
   catalogue (~37,600 proteomes) and shows a searchable dropdown. Filter by
   species name, proteome ID (`UP…`), or taxon ID, then click **Load proteome**.
   The selected proteome's `<ID>_<TAXID>.fasta.gz` is streamed for BLAST.
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
After a search, every hit is shown in a table with **E-value**, **score**,
**% identity** (BLAST only — phmmer does not report identity), and **alignment
length**. Live cutoffs (`E-value ≤`, `Score ≥`, `% identity ≥`, `Align length ≥`)
filter the table as you type. Use **Select all** / **Select none** or individual
checkboxes, then **Use selected for alignment** to pass only those sequences to
step 3. BLAST subject sequences are taken from the `sseq` column of the tabular
output, and phmmer hits are recovered from the in-memory DB by subject ID.

### Runtime timing
Each step's elapsed time is printed to the browser console
(`[phylo-fam-fate] …`) and shown in the **timing strip** above the pipeline.
The same numbers appear in the run metadata table and in the exported JSON.

### Remote APIs (CORS-enabled, callable directly from the browser)
| API | Use | Rate limit |
|---|---|---|
| [UniProt REST](https://rest.uniprot.org) | Protein sequence search & FASTA retrieval | ~10 req/s |
| [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/) | Nucleotide/protein search & retrieval | 3 req/s (no key) |
| [InterPro](https://www.ebi.ac.uk/interpro/api/) | Pfam/InterPro family search | fair use |

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
