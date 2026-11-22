# phylo-fam-fate command-line container: famfate.sh + pinned BLAST+/HMMER/
# Kalign/MUSCLE/FastTree + the website's JS (trim, rooting, annotation, export).
# Two FastTree 2.1.11 builds: FastTree-browser (single precision, identical to
# the website; the default) and bioconda's FastTree/FastTreeMP (double precision).
#   docker run --rm -u $(id -u):$(id -g) -v "$PWD":/data ghcr.io/hyphaltip/phylo-fam-fate-cli famfate.sh -q P46151
#   apptainer exec docker://ghcr.io/hyphaltip/phylo-fam-fate-cli famfate.sh -q P46151
# Tools compiled from source for portability / browser parity:
#  - FastTree-browser: single precision, strict FP, as the website (cli/build-fasttree-browser.sh)
#  - Kalign 3.3.1 without AVX autodetection (cli/build-kalign-portable.sh): bioconda's
#    kalign3 binary can die with SIGILL on CPUs lacking the build host's SIMD set.
FROM debian:bookworm-slim AS build
RUN apt-get update && apt-get install -y --no-install-recommends \
      gcc libc6-dev make autoconf automake libtool curl ca-certificates \
 && rm -rf /var/lib/apt/lists/*
COPY cli/build-fasttree-browser.sh cli/build-kalign-portable.sh /tmp/
RUN bash /tmp/build-fasttree-browser.sh /out && bash /tmp/build-kalign-portable.sh /out

FROM mambaorg/micromamba:1.5.10

LABEL org.opencontainers.image.title="phylo-fam-fate-cli" \
      org.opencontainers.image.description="Gene-family phylogeny pipeline: CLI twin of the phylo-fam-fate website" \
      org.opencontainers.image.source="https://github.com/hyphaltip/phylo-fam-fate" \
      org.opencontainers.image.licenses="MIT AND GPL-2.0-or-later"

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
 && micromamba clean --all --yes \
 && rm /tmp/environment.yml

# Only what the CLI needs from the site: JS modules, FastTree WASM (for the
# parity check), examples and the CLI itself.
COPY js/                 /opt/phylo-fam-fate/js/
COPY assets/fasttree/    /opt/phylo-fam-fate/assets/fasttree/
COPY assets/examples/    /opt/phylo-fam-fate/assets/examples/
COPY cli/                /opt/phylo-fam-fate/cli/
COPY --from=build /out/bin/FastTree-browser /usr/local/bin/FastTree-browser
# Replace conda's kalign (3.3.2, host-tuned SIMD) with the portable 3.3.1 build
COPY --from=build /out/bin/kalign /opt/conda/bin/kalign
USER root
# The source-built tools must run in *this* image (catches missing libraries)
RUN (/opt/conda/bin/kalign --version 2>&1 || true) | grep 3.3.1 \
 && (/usr/local/bin/FastTree-browser -expert 2>&1 || true) | grep -m1 'FastTree 2.1.11'
RUN chmod -R a+rX /opt/phylo-fam-fate \
 && chmod a+x /opt/phylo-fam-fate/cli/famfate.sh /opt/phylo-fam-fate/cli/famfate-js.mjs \
 && ln -s /opt/phylo-fam-fate/cli/famfate-js.mjs /usr/local/bin/famfate-js \
 && ln -s /opt/phylo-fam-fate/cli/famfate.sh /usr/local/bin/famfate.sh \
 && mkdir -p /data && chmod 1777 /data
USER $MAMBA_USER

# PATH is set explicitly (not only via the micromamba entrypoint) so that
# `apptainer exec` / `singularity exec`, which skip ENTRYPOINT, find the tools.
# FAMFATE_DB is deliberately unset: the default ./famfate_db is created in the
# working directory: /data under Docker (mounted), your own cwd under Apptainer
# (where the image itself, including /data, is read-only).
ENV PATH=/opt/conda/bin:/usr/local/bin:$PATH \
    FAMFATE_SITE=/opt/phylo-fam-fate \
    HOME=/tmp
WORKDIR /data
CMD ["famfate.sh", "--help"]
