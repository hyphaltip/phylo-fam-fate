# phylo-fam-fate command-line container: famfate.sh + pinned BLAST+/HMMER/
# Kalign/MUSCLE/FastTree + the website's JS (trim, rooting, annotation, export).
#   docker run --rm -u $(id -u):$(id -g) -v "$PWD":/data ghcr.io/hyphaltip/phylo-fam-fate-cli famfate.sh -q P46151
#   apptainer exec docker://ghcr.io/hyphaltip/phylo-fam-fate-cli famfate.sh -q P46151
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
USER root
RUN chmod -R a+rX /opt/phylo-fam-fate \
 && chmod a+x /opt/phylo-fam-fate/cli/famfate.sh /opt/phylo-fam-fate/cli/famfate-js.mjs \
 && ln -s /opt/phylo-fam-fate/cli/famfate-js.mjs /usr/local/bin/famfate-js \
 && ln -s /opt/phylo-fam-fate/cli/famfate.sh /usr/local/bin/famfate.sh \
 && mkdir -p /data && chmod 1777 /data
USER $MAMBA_USER

# PATH is set explicitly (not only via the micromamba entrypoint) so that
# `apptainer exec` / `singularity exec`, which skip ENTRYPOINT, find the tools.
ENV PATH=/opt/conda/bin:/usr/local/bin:$PATH \
    FAMFATE_SITE=/opt/phylo-fam-fate \
    FAMFATE_DB=/data/famfate_db \
    HOME=/tmp
WORKDIR /data
CMD ["famfate.sh", "--help"]
