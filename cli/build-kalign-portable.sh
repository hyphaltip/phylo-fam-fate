#!/usr/bin/env bash
# Build Kalign 3.3.1 (the version the website runs) as a portable binary.
#
# Kalign's configure (AX_EXT) enables every SIMD extension the *build*
# machine's CPU has. A binary built on an AVX-512 machine then dies with
# SIGILL ("illegal instruction", exit 132) on CPUs without it, as bioconda's
# kalign3 3.3.2 did on some CI runners. This build turns off AVX/AVX2/AVX-512
# detection, so only SSE up to 4.2 is used (any x86-64 CPU from ~2009 on), and
# runs anywhere. On ARM it is a plain build.
#
#   cli/build-kalign-portable.sh [PREFIX]   # installs PREFIX/bin/kalign (default: $CONDA_PREFIX or /usr/local)
# Needs: gcc, make, autoconf, automake, libtool, curl.
set -euo pipefail
PREFIX=${1:-${CONDA_PREFIX:-/usr/local}}
VERSION=3.3.1
URL=https://github.com/TimoLassmann/kalign/archive/refs/tags/v$VERSION.tar.gz
SHA256=7f10acf9a3fa15deabbc0304e7c14efa25cea39108318c9f02b47257de2d7390
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
curl -fsSL "$URL" -o "$tmp/kalign.tar.gz"
echo "$SHA256  $tmp/kalign.tar.gz" | sha256sum -c - >/dev/null || { echo "checksum mismatch for $URL" >&2; exit 1; }
tar xzf "$tmp/kalign.tar.gz" -C "$tmp"
cd "$tmp/kalign-$VERSION"
./autogen.sh >/dev/null 2>&1 || autoreconf -fi
# Pretend the OS lacks AVX support: disables -mavx/-mavx2/-mavx512* and HAVE_AVX*.
./configure --prefix="$tmp/install" ax_cv_have_avx_os_support_ext=no ax_cv_have_avx512_os_support_ext=no >/dev/null
grep -E '^SIMD_FLAGS' Makefile | head -n1
make -j"$(nproc 2>/dev/null || echo 2)" >/dev/null
make install >/dev/null
mkdir -p "$PREFIX/bin"
install -m 0755 "$tmp/install/bin/kalign" "$PREFIX/bin/kalign"
"$PREFIX/bin/kalign" --version 2>&1 | grep -m1 -i 'kalign' || true
echo "installed $PREFIX/bin/kalign ($VERSION, portable SIMD)"
