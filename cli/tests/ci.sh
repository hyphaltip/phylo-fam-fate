#!/usr/bin/env bash
# Container smoke + parity tests for famfate.sh (run by .github/workflows/cli-container.yml):
#   docker run --rm phylo-fam-fate-cli:test bash /opt/phylo-fam-fate/cli/tests/ci.sh
# Needs network for the search/annotation tests (UniProt FTP + REST, InterPro).
set -euo pipefail
SITE=${FAMFATE_SITE:-/opt/phylo-fam-fate}
EX=$SITE/assets/examples
W=$(mktemp -d)
cd "$W"
export FAMFATE_DB=$W/db
fail=0
ok()   { printf '  \033[32mPASS\033[0m %s\n' "$*"; }
bad()  { printf '  \033[31mFAIL\033[0m %s\n' "$*"; fail=1; }
check(){ local d=$1; shift; if "$@"; then ok "$d"; else bad "$d"; fi; }
tips() { grep -o '[,(][^(),:;]*:' "$1" | wc -l; }
section(){ printf '\n== %s\n' "$*"; }

section "tool versions"
blastp -version | head -n1
phmmer -h | grep -m1 HMMER
kalign --version 2>&1 | head -n1 || true
muscle -version 2>&1 | head -n1
FastTree -expert 2>&1 | head -n1
node --version
famfate.sh --version

section "1. family → kalign → trim → FastTree (offline)"
famfate.sh --family "$EX/proteins.fa" --no-annotate -o t1 -T 2 >/dev/null 2>t1.err || { cat t1.err; bad "pipeline ran"; }
check "aligned FASTA written"   test -s t1/family.aln.fa
check "tree has 16 tips"        test "$(tips t1/family.nwk)" -eq 16
check "rooted tree + SVG"       test -s t1/family.rooted.nwk -a -s t1/family.svg
check "run.json records steps"  grep -q '"step":"tree"' t1/run.json

section "2. MUSCLE path, no trimming, WAG"
famfate.sh --family "$EX/proteins.fa" --aligner muscle --gap none --model wag --no-gamma --no-annotate -o t2 -T 2 >/dev/null 2>t2.err || { cat t2.err; bad "pipeline ran"; }
check "MUSCLE alignment"        test -s t2/family.aln.fa
check "tree has 16 tips"        test "$(tips t2/family.nwk)" -eq 16

section "3. FastTree parity: native vs the browser's WASM build"
parity() {   # parity NAME ALN ARGS…
  local name=$1 aln=$2; shift 2
  FastTree "$@" "$aln" > "$name.native.nwk" 2>/dev/null
  famfate-js fasttree-wasm "$aln" "$@" > "$name.wasm.nwk"
  if cmp -s "$name.native.nwk" "$name.wasm.nwk"; then ok "$name: identical Newick"; return; fi
  # Not byte-identical (floating-point build flags can differ): require the same
  # topology and branch lengths within 1e-3.
  if famfate-js compare-trees "$name.native.nwk" "$name.wasm.nwk" --tol 1e-3 > "$name.cmp.json"; then
    ok "$name: same topology, lengths within 1e-3 ($(grep maxLengthDiff "$name.cmp.json" | tr -d ' ,'))"
  else
    bad "$name: trees differ"; cat "$name.cmp.json"
  fi
}
parity nt_default "$EX/nucleotide_aligned.fa" -nt
parity nt_gtr_gamma "$EX/nucleotide_aligned.fa" -nt -gtr -gamma
parity aa_default t1/family.aln.fa
parity aa_lg_gamma t1/family.aln.fa -lg -gamma

section "4. MET12 vs two small proteomes (BLAST, full-length hits, annotation, MAD)"
# S. cerevisiae (UP000002311) + E. coli K-12 (UP000000625): small, quick to index
famfate.sh -q P46151 -p UP000002311,UP000000625 --no-sprot -E 1e-5 --root mad -o t4 -T 2 2>t4.err >/dev/null || { tail -n 30 t4.err; bad "pipeline ran"; }
check "BLAST hits found"                   test -s t4/hits.tsv
check "all hits full-length (blastdbcmd)"  grep -q ' 0 HSP-region' t4/hits.summary
check "MET13 (P53128) recovered"           grep -q 'P53128' t4/family.fa
check "E. coli MetF (P0AEZ1) recovered"    grep -q 'P0AEZ1' t4/family.fa
check "annotations TSV"                    test "$(wc -l < t4/family.annotations.tsv)" -ge 4
check "taxonomy present (Ascomycota)"      grep -q Ascomycota t4/family.annotations.tsv
check "MET12/MET13 duplication labelled"   grep -q duplication t4/family.events.tsv
check "iTOL files"                         test -s t4/family.itol_phylum.txt -a -s t4/family.itol_domains.txt
check "MAD rooting used"                   grep -q '"root": "mad"' t4/family.post.json

section "5. phmmer path + outgroup rooting"
famfate.sh -q P46151 -p UP000002311,UP000000625 --no-sprot -t phmmer -E 1e-5 \
  --root outgroup:domain:Bacteria --no-domains -o t5 -T 2 2>t5.err >/dev/null || { tail -n 30 t5.err; bad "pipeline ran"; }
check "phmmer hits"                        test -s t5/phmmer.tbl
check "outgroup rooting reported"          grep -q 'Bacteria' t5/family.post.json

section "6. option errors"
check "rejects two inputs"   bash -c '! famfate.sh -q P46151 --family x.fa 2>/dev/null'
check "rejects bad tool"     bash -c '! famfate.sh -q P46151 -t diamond 2>/dev/null'

echo
if [ "$fail" = 0 ]; then echo "All container tests passed"; else echo "Some container tests FAILED"; fi
exit $fail
