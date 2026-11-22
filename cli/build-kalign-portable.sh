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
# Pretend the OS lacks AVX support (disables -mavx/-mavx2/-mavx512*/FMA and
# HAVE_AVX*), and skip SHA/AES and AMD-only SSE4a, which not every x86-64 CPU
# has. The result is the same instruction set whichever runner builds it:
# MMX + SSE … SSE4.2.
# --without-openmp: no libgomp runtime dependency (the browser's Kalign is
# single-threaded too), so the binary runs in any image.
./configure --prefix="$tmp/install" --without-openmp \
  ax_cv_have_avx_os_support_ext=no ax_cv_have_avx512_os_support_ext=no \
  ax_cv_have_sha_cpu_ext=no ax_cv_have_aes_cpu_ext=no ax_cv_have_sse4a_cpu_ext=no >/dev/null
flags=$(grep -E '^SIMD_FLAGS' Makefile | head -n1); echo "$flags"
# Guard: nothing beyond SSE4.2 may slip in, whatever CPU this build runs on
if echo "$flags" | grep -qE -- '-m(avx|fma|sha|aes|sse4a|xop|bmi|f16c)'; then
  echo "non-portable SIMD flags enabled: $flags" >&2; exit 1
fi
make -j"$(nproc 2>/dev/null || echo 2)" >/dev/null
make install >/dev/null
mkdir -p "$PREFIX/bin"
install -m 0755 "$tmp/install/bin/kalign" "$PREFIX/bin/kalign"
# Fail loudly if the binary can't run (missing library, bad instructions)
# (judge by the output: some versions exit non-zero even for --version)
ver=$("$PREFIX/bin/kalign" --version 2>&1 || true)
echo "$ver" | grep -q "$VERSION" || { echo "built kalign does not run or reports an unexpected version: $ver" >&2; exit 1; }
if command -v ldd >/dev/null; then echo "kalign links: $(ldd "$PREFIX/bin/kalign" | awk '{print $1}' | tr '\n' ' ')"; fi
echo "installed $PREFIX/bin/kalign ($VERSION, portable SIMD, no OpenMP)"
