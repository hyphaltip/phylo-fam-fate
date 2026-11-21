#!/usr/bin/env bash
# famfate.sh — command-line version of the phylo-fam-fate browser pipeline.
#
# Same tools and defaults as https://hyphaltip.github.io/phylo-fam-fate/ :
#   search (BLAST+/phmmer) → full-length hits → align (Kalign/MUSCLE) → trim
#   → FastTree → root (midpoint/MAD/outgroup) → UniProt/InterPro annotation
#   → duplication labels → SVG, Newick, TSV and iTOL exports.
# Trimming, rooting, annotation and export run the website's own JavaScript
# (cli/famfate-js.mjs), so they match the browser exactly.
#
# Run with -h for options. Every step is logged to OUTDIR/run.log and
# OUTDIR/run.json records tools, versions and parameters.
set -euo pipefail

VERSION=0.1.0
HERE=$(cd "$(dirname "$(readlink -f "$0")")" && pwd)
FAMFATE_JS=${FAMFATE_JS:-"node $HERE/famfate-js.mjs"}
UP=https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase

usage() {
  cat <<EOF
famfate.sh $VERSION — gene-family phylogeny pipeline (CLI twin of the phylo-fam-fate website)

Input (one of):
  -q, --query ACC[,ACC…]   UniProt accession(s) of the query protein(s)        e.g. P46151
  -i, --query-fasta FILE   query FASTA (one or more sequences)
  -f, --family FILE        unaligned family FASTA: skip the search
  -a, --aligned FILE       aligned FASTA: skip search and alignment

Search (default: blastp vs Swiss-Prot):
  -p, --proteomes LIST     reference proteomes, comma-separated UniProt IDs    e.g. UP000001805,UP000005640
      --no-sprot           don't include Swiss-Prot (needs --proteomes or --db)
      --db FILE            extra FASTA database(s), comma-separated (.fa or .fa.gz)
  -t, --tool blastp|blastn|phmmer   search program                              [blastp]
  -e, --evalue X           search E-value (BLAST -evalue, phmmer -E/--domE)    [1e-3]
  -E, --hit-evalue X       E-value a hit must reach to be kept (browser's live cutoff) [1e-5]
  -n, --max-hits N         max subject sequences (BLAST -max_target_seqs)      [50]
      --hsp-only           use only the matched region, not full-length subjects
      --no-query           don't add the query sequence(s) to the family

Align / trim / tree:
  -A, --aligner kalign|muscle                                                  [kalign]
  -g, --gap X              trim columns with gap fraction > X; 'none' = no trim [0.5]
  -c, --cons X             trim columns with majority-residue fraction < X     [0]
  -m, --model jtt|lg|wag|nt|gtr   FastTree model (nt/gtr for nucleotides)      [lg]
      --no-gamma           omit -gamma (on by default)
      --fasttree browser|double  FastTree build: 'browser' = single precision, byte-identical
                           to the website's default (FastTree-browser); 'double' = bioconda's
                           double-precision FastTree/FastTreeMP (multithreaded), matching the
                           website's "Double" precision option              [browser]
      --fasttree-args "…"  extra FastTree arguments

Post-processing:
  -r, --root midpoint|mad|none|outgroup:RANK:NAME                              [midpoint]
                           e.g. outgroup:domain:Bacteria, outgroup:phylum:Microsporidia
      --rank RANK          colour/iTOL rank (domain…species)                   [phylum]
      --no-annotate        skip UniProt/InterPro lookups (offline)
      --no-domains         skip InterPro Pfam domains

General:
  -o, --outdir DIR         output directory                                    [famfate_out]
  -d, --dbdir DIR          cache for downloaded/indexed databases              [\$FAMFATE_DB or ./famfate_db]
  -T, --threads N          threads for BLAST/phmmer/MUSCLE/FastTreeMP          [nproc]
  -h, --help / -V, --version
EOF
}

# ---------------------------------------------------------------- options
ORIG_ARGS=("$@")
QUERY= QFASTA= FAMILY= ALIGNED= PROTEOMES= EXTRA_DB= SPROT=1
TOOL=blastp EVALUE=1e-3 HIT_EVALUE=1e-5 MAXHITS=50 HSP_ONLY=0 ADD_QUERY=1
ALIGNER=kalign GAP=0.5 CONS=0 MODEL=lg GAMMA=1 FT_EXTRA= FT_BUILD=browser
ROOT=midpoint RANK=phylum ANNOTATE=1 DOMAINS=1
OUTDIR=famfate_out DBDIR=${FAMFATE_DB:-famfate_db} THREADS=$(nproc 2>/dev/null || echo 2)

while [ $# -gt 0 ]; do
  case $1 in
    -q|--query) QUERY=$2; shift ;;
    -i|--query-fasta) QFASTA=$2; shift ;;
    -f|--family) FAMILY=$2; shift ;;
    -a|--aligned) ALIGNED=$2; shift ;;
    -p|--proteomes) PROTEOMES=$2; shift ;;
    --no-sprot) SPROT=0 ;;
    --db) EXTRA_DB=$2; shift ;;
    -t|--tool) TOOL=$2; shift ;;
    -e|--evalue) EVALUE=$2; shift ;;
    -E|--hit-evalue) HIT_EVALUE=$2; shift ;;
    -n|--max-hits) MAXHITS=$2; shift ;;
    --hsp-only) HSP_ONLY=1 ;;
    --no-query) ADD_QUERY=0 ;;
    -A|--aligner) ALIGNER=$2; shift ;;
    -g|--gap) GAP=$2; shift ;;
    -c|--cons) CONS=$2; shift ;;
    -m|--model) MODEL=$2; shift ;;
    --no-gamma) GAMMA=0 ;;
    --fasttree-args) FT_EXTRA=$2; shift ;;
    --fasttree) FT_BUILD=$2; shift ;;
    -r|--root) ROOT=$2; shift ;;
    --rank) RANK=$2; shift ;;
    --no-annotate) ANNOTATE=0 ;;
    --no-domains) DOMAINS=0 ;;
    -o|--outdir) OUTDIR=$2; shift ;;
    -d|--dbdir) DBDIR=$2; shift ;;
    -T|--threads) THREADS=$2; shift ;;
    -h|--help) usage; exit 0 ;;
    -V|--version) echo "famfate.sh $VERSION"; exit 0 ;;
    *) echo "famfate.sh: unknown option '$1' (see -h)" >&2; exit 2 ;;
  esac
  shift
done

n_inputs=0
for v in "$QUERY" "$QFASTA" "$FAMILY" "$ALIGNED"; do [ -n "$v" ] && n_inputs=$((n_inputs + 1)); done
[ "$n_inputs" -eq 1 ] || { echo "famfate.sh: give exactly one of --query, --query-fasta, --family, --aligned (see -h)" >&2; exit 2; }
case $TOOL in blastp|blastn|phmmer) ;; *) echo "famfate.sh: --tool must be blastp, blastn or phmmer" >&2; exit 2 ;; esac
case $ALIGNER in kalign|muscle) ;; *) echo "famfate.sh: --aligner must be kalign or muscle" >&2; exit 2 ;; esac

# Input paths are resolved before we cd into OUTDIR
abspath() { [ -z "$1" ] && return; [ -e "$1" ] || { echo "famfate.sh: no such file: $1" >&2; exit 2; }; echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"; }
QFASTA=$(abspath "$QFASTA"); FAMILY=$(abspath "$FAMILY"); ALIGNED=$(abspath "$ALIGNED")
if [ -n "$EXTRA_DB" ]; then
  _dbs=; for f in ${EXTRA_DB//,/ }; do _dbs+="${_dbs:+,}$(abspath "$f")"; done; EXTRA_DB=$_dbs
fi
mkdir -p "$OUTDIR" "$DBDIR"
OUTDIR=$(cd "$OUTDIR" && pwd); DBDIR=$(cd "$DBDIR" && pwd)
LOG=$OUTDIR/run.log
: > "$LOG"
STEPS=()
say()  { printf '[famfate %s] %s\n' "$(date +%H:%M:%S)" "$*" | tee -a "$LOG" >&2; }
die()  { say "ERROR: $*"; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "'$1' not found on PATH (use the container, or: mamba env create -f environment.yml)"; }
# record STEP TOOL VERSION PARAMS — appended to run.json at the end
record() { STEPS+=("$(printf '{"step":"%s","tool":"%s","version":"%s","params":"%s","seconds":%s}' "$1" "$2" "$3" "${4//\"/\\\"}" "$5")"); }
# timed CMD… : run CMD with output to the log; print its runtime (s). On failure,
# show the command's exit status and last lines of output, then stop.
timed()  {
  local t0=$SECONDS out rc=0
  out=$(mktemp)
  "$@" >"$out" 2>&1 || rc=$?
  cat "$out" >> "$LOG"
  if [ "$rc" -ne 0 ]; then
    say "command failed (exit $rc): $*"
    tail -n 25 "$out" | sed 's/^/    | /' | tee -a "$LOG" >&2
    rm -f "$out"; die "stopping (full log: $LOG)"
  fi
  rm -f "$out"; echo $((SECONDS - t0))
}
nseq()   { grep -c '^>' "$1" || true; }

need node
cd "$OUTDIR"

# ---------------------------------------------------------------- 1. input
if [ -n "$QUERY" ]; then
  : > query.fa
  for acc in ${QUERY//,/ }; do
    curl -fsS "https://rest.uniprot.org/uniprotkb/$acc.fasta" >> query.fa || die "could not fetch UniProt $acc"
  done
  say "query: $(nseq query.fa) sequence(s) from UniProt ($QUERY)"
elif [ -n "$QFASTA" ]; then
  cp "$QFASTA" query.fa; say "query: $(nseq query.fa) sequence(s) from $QFASTA"
fi

# ---------------------------------------------------------------- 2. search
fetch_db() {   # fetch_db NAME URL → $DBDIR/NAME.fa (cached)
  local name=$1 url=$2
  if [ ! -s "$DBDIR/$name.fa" ]; then
    say "downloading $name"
    curl -fsS "$url" | gunzip -c > "$DBDIR/$name.fa.tmp" || die "download failed: $url"
    mv "$DBDIR/$name.fa.tmp" "$DBDIR/$name.fa"
  fi
}
proteome_url() {   # UPID → FTP URL, via the reference_proteomes README
  [ -s "$DBDIR/reference_proteomes.README" ] || curl -fsS "$UP/reference_proteomes/README" -o "$DBDIR/reference_proteomes.README"
  awk -F'\t' -v id="$1" '$1==id { d=toupper(substr($4,1,1)) substr($4,2); print d "/" $1 "/" $1 "_" $2 ".fasta.gz"; exit }' \
      "$DBDIR/reference_proteomes.README"
}

if [ -s query.fa ]; then
  DBS=()
  if [ "$SPROT" = 1 ]; then fetch_db sprot "$UP/complete/uniprot_sprot.fasta.gz"; DBS+=(sprot); fi
  for p in ${PROTEOMES//,/ }; do
    rel=$(proteome_url "$p"); [ -n "$rel" ] || die "proteome $p not in the UniProt reference_proteomes README"
    fetch_db "$p" "$UP/reference_proteomes/$rel"; DBS+=("$p")
  done
  for f in ${EXTRA_DB//,/ }; do
    n=$(basename "$f"); n=${n%%.*}
    [ -s "$DBDIR/$n.fa" ] || { case $f in *.gz) gunzip -c "$f" ;; *) cat "$f" ;; esac > "$DBDIR/$n.fa"; }
    DBS+=("$n")
  done
  [ ${#DBS[@]} -gt 0 ] || die "no database: drop --no-sprot or give --proteomes/--db"
  say "search set: ${DBS[*]}"

  if [ "$TOOL" = phmmer ]; then
    need phmmer; need esl-sfetch
    files=(); for d in "${DBS[@]}"; do files+=("$DBDIR/$d.fa"); done
    cat "${files[@]}" > searchset.fa
    s=$(timed phmmer --cpu "$THREADS" -E "$EVALUE" --domE "$EVALUE" --noali --tblout phmmer.tbl -o /dev/null query.fa searchset.fa)
    record search phmmer "$(phmmer -h | grep -m1 -o 'HMMER [0-9.]*')" "-E $EVALUE --domE $EVALUE" "$s"
    # Columns: target query … full-seq E-value (5th). Keep hits passing --hit-evalue, best first.
    grep -v '^#' phmmer.tbl | awk -v e="$HIT_EVALUE" '$5+0 <= e+0 {print $1}' | awk '!seen[$0]++' | head -n "$MAXHITS" > subject_ids.txt
    [ -s subject_ids.txt ] || die "phmmer found no hits at E ≤ $HIT_EVALUE"
    esl-sfetch --index searchset.fa >>"$LOG" 2>&1
    esl-sfetch -f searchset.fa subject_ids.txt > hits.fa
    [ "$HSP_ONLY" = 1 ] && say "note: --hsp-only applies to BLAST; phmmer hits are full-length"
  else
    need "$TOOL"; need makeblastdb; need blastdbcmd
    dbtype=prot; [ "$TOOL" = blastn ] && dbtype=nucl
    for d in "${DBS[@]}"; do
      if [ ! -s "$DBDIR/$d.${dbtype:0:1}in" ] && [ ! -s "$DBDIR/$d.00.${dbtype:0:1}in" ]; then
        say "indexing $d"
        makeblastdb -in "$DBDIR/$d.fa" -dbtype $dbtype -out "$DBDIR/$d" -parse_seqids >>"$LOG" 2>&1 \
          || { say "  -parse_seqids rejected $d; indexing without (those hits will be HSP-only)";
               makeblastdb -in "$DBDIR/$d.fa" -dbtype $dbtype -out "$DBDIR/$d" >>"$LOG" 2>&1 || die "makeblastdb failed for $d"; }
      fi
    done
    dbarg=$(printf "$DBDIR/%s " "${DBS[@]}"); dbarg=${dbarg% }
    filter=(-seg no); [ "$TOOL" = blastn ] && filter=(-task blastn -dust no)
    s=$(timed "$TOOL" -query query.fa -db "$dbarg" "${filter[@]}" -evalue "$EVALUE" -max_target_seqs "$MAXHITS" \
        -num_threads "$THREADS" -out hits.tsv \
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq sacc slen')
    record search "$TOOL" "$("$TOOL" -version | head -n1)" "-evalue $EVALUE -max_target_seqs $MAXHITS ${filter[*]}" "$s"
    # One entry per subject (best HSP first), passing the hit E-value cutoff
    sort -t$'\t' -k11,11g -k12,12gr hits.tsv | awk -F'\t' -v e="$HIT_EVALUE" '$11+0 <= e+0 && !seen[$14]++' > best_hsp.tsv
    [ -s best_hsp.tsv ] || die "$TOOL found no hits at E ≤ $HIT_EVALUE"
    cut -f14 best_hsp.tsv > subject_ids.txt
    say "$(wc -l < hits.tsv) HSPs → $(wc -l < subject_ids.txt) subject proteins at E ≤ $HIT_EVALUE"
    # Name each hit by its sseqid (as the browser does). Full-length sequence and
    # description (OS=…, used for offline species) come from blastdbcmd by
    # accession; subjects it can't return fall back to the best HSP region.
    : > full.tsv
    if [ "$HSP_ONLY" = 0 ]; then
      blastdbcmd -db "$dbarg" -entry_batch subject_ids.txt -outfmt $'%a\t%s\t%t' > full.tsv 2>>"$LOG" || true
    fi
    awk -F'\t' '
      NR == FNR { a = $1; sub(/\.[0-9]+$/, "", a); if (!(a in seq)) { seq[a] = $2; title[a] = $3 }; next }
      { acc = $14; sub(/\.[0-9]+$/, "", acc)
        if (acc in seq) { print ">" $2 (title[acc] != "" ? " " title[acc] : "") "\n" seq[acc]; full++ }
        else { s = $13; gsub(/[-.]/, "", s); print ">" $2 "\n" s; hsp++ } }
      END { printf "%d full-length, %d HSP-region\n", full, hsp > "/dev/stderr" }' full.tsv best_hsp.tsv > hits.fa 2> hits.summary
    say "  hits: $(cat hits.summary)"
  fi
  if [ "$ADD_QUERY" = 1 ]; then cat query.fa hits.fa > family.raw.fa; else cp hits.fa family.raw.fa; fi
  # Drop duplicate names (a query can also be a hit), keeping the first
  awk '/^>/ { split(substr($0,2), w, " "); keep = !seen[w[1]]++ } keep' family.raw.fa > family.fa
  say "family: $(nseq family.fa) sequences"
elif [ -n "$FAMILY" ]; then
  cp "$FAMILY" family.fa; say "family: $(nseq family.fa) sequences from $FAMILY"
fi

# ---------------------------------------------------------------- 3. align
if [ -n "$ALIGNED" ]; then
  cp "$ALIGNED" family.aln.fa; say "alignment: $(nseq family.aln.fa) sequences from $ALIGNED"
else
  [ "$(nseq family.fa)" -ge 2 ] || die "need at least 2 sequences to align (got $(nseq family.fa))"
  need "$ALIGNER"
  if [ "$ALIGNER" = kalign ]; then
    s=$(timed kalign -i family.fa -o family.aln.fa)
    record align kalign "$(kalign --version 2>&1 | grep -m1 -o '[0-9][0-9.]*')" "-i -o (defaults)" "$s"
  else
    s=$(timed muscle -align family.fa -output family.aln.fa -threads "$THREADS")
    record align muscle "$(muscle -version 2>&1 | grep -m1 -o '[0-9][0-9.]*')" "-align (defaults)" "$s"
  fi
  say "aligned with $ALIGNER (${s}s)"
fi

# ---------------------------------------------------------------- 4. trim
TYPE=$($FAMFATE_JS seqtype family.aln.fa)
if [ "$GAP" = none ]; then
  cp family.aln.fa family.trim.fa
else
  t0=$SECONDS
  $FAMFATE_JS trim family.aln.fa --gap "$GAP" --cons "$CONS" > family.trim.fa 2>>"$LOG" || die "trimming failed (see $LOG)"
  record trim famfate-js "$VERSION" "gap<=$GAP cons>=$CONS" $((SECONDS - t0))
  say "$(tail -n1 "$LOG")"
fi

# ---------------------------------------------------------------- 5. tree
case $FT_BUILD in
  browser)
    FT=FastTree-browser
    command -v "$FT" >/dev/null || die "FastTree-browser not found: run cli/build-fasttree-browser.sh, or use --fasttree double" ;;
  double)
    FT=FastTree; [ "$THREADS" -gt 1 ] && command -v FastTreeMP >/dev/null && FT=FastTreeMP
    need "$FT" ;;
  *) die "--fasttree must be 'browser' or 'double'" ;;
esac
ftargs=()
case $MODEL in
  nt)  ftargs=(-nt) ;;
  gtr) ftargs=(-nt -gtr) ;;
  lg)  ftargs=(-lg) ;;
  wag) ftargs=(-wag) ;;
  jtt) ftargs=() ;;
  *) die "--model must be jtt, lg, wag, nt or gtr" ;;
esac
if [ "$TYPE" = nt ] && [[ ! " ${ftargs[*]} " =~ " -nt " ]]; then ftargs=(-nt); say "nucleotide alignment detected: using -nt"; fi
[ "$GAMMA" = 1 ] && ftargs+=(-gamma)
# shellcheck disable=SC2206
[ -n "$FT_EXTRA" ] && ftargs+=($FT_EXTRA)
t0=$SECONDS
OMP_NUM_THREADS=$THREADS "$FT" "${ftargs[@]}" family.trim.fa > family.nwk 2>>"$LOG" || die "FastTree failed (see $LOG)"
record tree "$FT" "$("$FT" -expert 2>&1 | head -n1 | sed 's/^Detailed usage for FastTree //; s/:$//')" "${ftargs[*]}" $((SECONDS - t0))
say "tree: $(grep -o ',' family.nwk | wc -l | awk '{print $1+1}') tips ($FT ${ftargs[*]})"

# ---------------------------------------------------------------- 6–9. root, annotate, duplications, export
postargs=(post --tree family.nwk --fasta family.fa --root "$ROOT" --rank "$RANK" --out family)
[ -s family.fa ] || postargs=(post --tree family.nwk --fasta family.aln.fa --root "$ROOT" --rank "$RANK" --out family)
[ "$ANNOTATE" = 0 ] && postargs+=(--no-annotate)
[ "$DOMAINS" = 0 ] && postargs+=(--no-domains)
t0=$SECONDS
$FAMFATE_JS "${postargs[@]}" > post_outputs.txt 2> >(tee -a "$LOG" >&2) || die "post-processing failed (see $LOG)"
record post famfate-js "$VERSION" "root=$ROOT rank=$RANK annotate=$ANNOTATE domains=$DOMAINS" $((SECONDS - t0))

# ---------------------------------------------------------------- provenance
{
  printf '{\n  "tool": "famfate.sh", "version": "%s", "generated": "%s",\n' "$VERSION" "$(date -u +%FT%TZ)"
  printf '  "command": "%s",\n  "steps": [\n    ' "$(printf '%q ' "$0" "${ORIG_ARGS[@]:-}" | sed 's/"/\\"/g')"
  sep=; for st in "${STEPS[@]}"; do printf '%s%s' "$sep" "$st"; sep=$',\n    '; done
  printf '\n  ]\n}\n'
} > run.json

say "done. Outputs in $OUTDIR:"
for f in family.fa family.aln.fa family.trim.fa family.nwk $(cat post_outputs.txt | xargs -n1 basename 2>/dev/null) run.json run.log; do
  [ -e "$f" ] && printf '  %s\n' "$f" >&2
done
