#!/usr/bin/env bash
# Build "FastTree-browser": FastTree 2.1.11 compiled exactly like the website's
# WebAssembly build: single precision (no -DUSE_DOUBLE), no SSE, and strict
# floating point (-ffp-contract=off, so no fused multiply-adds). It produces
# trees byte-identical to the browser; CI checks this.
#
# Bioconda's FastTree is the double-precision variant (-DUSE_DOUBLE), which
# changes minimum branch lengths and can change nucleotide topologies; use it
# with `famfate.sh --fasttree double`.
#
#   cli/build-fasttree-browser.sh [PREFIX]      # installs PREFIX/bin/FastTree-browser (default: $CONDA_PREFIX or /usr/local)
set -euo pipefail
PREFIX=${1:-${CONDA_PREFIX:-/usr/local}}
URL=https://raw.githubusercontent.com/morgannprice/fasttree/main/old/FastTree-2.1.11.c
SHA256=9026ae550307374be92913d3098f8d44187d30bea07902b9dcbfb123eaa2050f
CC=${CC:-gcc}
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
curl -fsSL "$URL" -o "$tmp/FastTree.c"
echo "$SHA256  $tmp/FastTree.c" | sha256sum -c - >/dev/null || { echo "checksum mismatch for $URL" >&2; exit 1; }
"$CC" -O3 -DNO_SSE -DNDEBUG -ffp-contract=off -w -o "$tmp/FastTree-browser" "$tmp/FastTree.c" -lm
mkdir -p "$PREFIX/bin"
install -m 0755 "$tmp/FastTree-browser" "$PREFIX/bin/FastTree-browser"
echo "installed $PREFIX/bin/FastTree-browser"
