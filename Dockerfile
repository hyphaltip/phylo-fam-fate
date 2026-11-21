# phylo-fam-fate command-line container: famfate.sh + pinned BLAST+/HMMER/
# Kalign/MUSCLE/FastTree + the website's JS (trim, rooting, annotation, export).
# Two FastTree 2.1.11 builds: FastTree-browser (single precision, identical to
# the website; the default) and bioconda's FastTree/FastTreeMP (double precision).
#   docker run --rm -u $(id -u):$(id -g) -v "$PWD":/data ghcr.io/hyphaltip/phylo-fam-fate-cli famfate.sh -q P46151
#   apptainer exec docker://ghcr.io/hyphaltip/phylo-fam-fate-cli famfate.sh -q P46151
# FastTree as the browser runs it: single precision, strict FP (see cli/build-fasttree-browser.sh)
FROM debian:bookworm-slim AS fasttree
RUN apt-get update && apt-get install -y --no-install-recommends gcc libc6-dev curl ca-certificates \
 && rm -rf /var/lib/apt/lists/*
COPY cli/build-fasttree-browser.sh /tmp/
RUN bash /tmp/build-fasttree-browser.sh /out

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
COPY --from=fasttree /out/bin/FastTree-browser /usr/local/bin/FastTree-browser
USER root
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
