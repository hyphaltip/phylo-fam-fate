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
          │  (Biowasm WASM, user-provided database)     │
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
