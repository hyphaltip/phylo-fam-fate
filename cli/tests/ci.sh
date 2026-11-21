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
mafft --version 2>&1 | head -n1
iqtree3 --version 2>&1 | head -n1
raxml-ng --version 2>&1 | grep -m1 RAxML
FastTree -expert 2>&1 | head -n1
FastTree-browser -expert 2>&1 | head -n1
node --version
famfate.sh --version

section "1. family → kalign → trim → FastTree (offline)"
famfate.sh --family "$EX/proteins.fa" --no-annotate -o t1 -T 2 >/dev/null 2>t1.err || { cat t1.err; bad "pipeline ran"; }
check "aligned FASTA written"   test -s t1/family.aln.fa
check "browser FastTree used"   grep -q '"tool":"FastTree-browser"' t1/run.json
check "tree has 16 tips"        test "$(tips t1/family.nwk)" -eq 16
check "rooted tree + SVG"       test -s t1/family.rooted.nwk -a -s t1/family.svg
check "run.json records steps"  grep -q '"step":"tree"' t1/run.json

section "2. MUSCLE path, no trimming, WAG"
famfate.sh --family "$EX/proteins.fa" --aligner muscle --gap none --model wag --no-gamma --no-annotate -o t2 -T 2 >/dev/null 2>t2.err || { cat t2.err; bad "pipeline ran"; }
check "MUSCLE alignment"        test -s t2/family.aln.fa
check "tree has 16 tips"        test "$(tips t2/family.nwk)" -eq 16

section "2b. double-precision FastTree option"
famfate.sh --family "$EX/proteins.fa" --fasttree double --no-annotate -o t2b -T 2 >/dev/null 2>t2b.err || { cat t2b.err; bad "pipeline ran"; }
check "bioconda FastTree used"  grep -qE '"tool":"FastTree(MP)?"' t2b/run.json

section "3. FastTree parity: native builds vs the website's WASM builds"
# Single precision: FastTree-browser is compiled like the website's default WASM
# build (strict FP), so its Newick must be byte-identical.
# Double precision: bioconda's FastTree vs the website's double WASM build.
# Same topology required; numbers may differ slightly (bioconda vs Rust libm).
parity() {   # parity NAME ALN ARGS…
  local name=$1 aln=$2; shift 2
  FastTree-browser "$@" "$aln" > "$name.native.nwk" 2>/dev/null || true
  famfate-js fasttree-wasm "$aln" "$@" > "$name.wasm.nwk" || true
  check "$name [$*] single: identical Newick" cmp -s "$name.native.nwk" "$name.wasm.nwk"
  FastTree "$@" "$aln" > "$name.double.nwk" 2>/dev/null || true
  famfate-js fasttree-wasm "$aln" --double "$@" > "$name.wasm-double.nwk" || true
  if cmp -s "$name.double.nwk" "$name.wasm-double.nwk"; then ok "$name [$*] double: identical Newick"
  else check "$name [$*] double: same topology, lengths within 1e-5" \
         famfate-js compare-trees "$name.double.nwk" "$name.wasm-double.nwk" --tol 1e-5
  fi
}
parity nt_jc "$EX/nucleotide_aligned.fa" -nt
parity nt_gtr_gamma "$EX/nucleotide_aligned.fa" -nt -gtr -gamma
AA=$SITE/cli/tests/data/aa_aligned.fa
parity aa_jtt "$AA"
parity aa_lg_gamma "$AA" -lg -gamma
parity aa_wag "$AA" -wag
# The pipeline itself (default --fasttree browser) must reproduce the website's tree
famfate.sh --aligned "$AA" --no-annotate -o t3 >/dev/null 2>t3.err || { cat t3.err; bad "pipeline ran"; }
famfate-js fasttree-wasm t3/family.trim.fa -lg -gamma > t3.wasm.nwk
check "famfate.sh tree == website tree" cmp -s t3/family.nwk t3.wasm.nwk

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

section "6. ML trees: IQ-TREE + RAxML-NG + FastTree, compared"
famfate.sh --aligned "$SITE/cli/tests/data/aa_aligned.fa" --tree iqtree,raxml-ng,fasttree --no-annotate -o t6 -T 2 2>t6.err >/dev/null || { tail -n 40 t6.err; bad "pipeline ran"; }
check "IQ-TREE tree"                       test -s t6/family.iqtree.nwk
check "RAxML-NG tree"                      test -s t6/family.raxml-ng.nwk
check "FastTree tree"                      test -s t6/family.fasttree.nwk
check "main tree is IQ-TREE's"             cmp -s t6/family.nwk t6/family.iqtree.nwk
check "SH-aLRT/UFBoot labels kept after rooting" grep -qE '\)[0-9.]+/[0-9.]+:' t6/family.rooted.nwk
check "trees_compare.tsv has 2 rows"       test "$(tail -n +2 t6/trees_compare.tsv | wc -l)" -eq 2
check "IQ-TREE model recorded"             grep -q 'model ' t6/run.json
famfate.sh --family "$EX/proteins.fa" --aligner mafft --no-annotate -o t6b -T 2 >/dev/null 2>t6b.err || { cat t6b.err; bad "pipeline ran"; }
check "MAFFT alignment, 16 sequences"      test "$(grep -c '>' t6b/family.aln.fa)" -eq 16

section "7. option errors"
check "rejects two inputs"   bash -c '! famfate.sh -q P46151 --family x.fa 2>/dev/null'
check "rejects bad tool"     bash -c '! famfate.sh -q P46151 -t diamond 2>/dev/null'
check "rejects bad tree"     bash -c '! famfate.sh -a /opt/phylo-fam-fate/cli/tests/data/aa_aligned.fa --tree phyml --no-annotate -o bad 2>/dev/null'

echo
if [ "$fail" = 0 ]; then echo "All container tests passed"; else echo "Some container tests FAILED"; fi
exit $fail
