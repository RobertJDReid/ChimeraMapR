# ChimeraMapR command-line interface in a Linux container.
#
# Build (from the repo root):
#   docker build -t chimeramapr .
# Run (mount the folder that holds your inputs at /data):
#   docker run --rm -v "$(pwd):/data" chimeramapr \
#     -o /data/results/ /data/reads.csv.gz /data/snps.vcf.gz /data/genome.fasta.fai
#   (PowerShell: use -v "${PWD}:/data")
FROM rocker/r-ver:4.5.1

# igraph's Linux binary links against GLPK and libxml2. Without them
# loadNamespace("igraph") fails and the events/chain step cannot run.
RUN apt-get update \
 && apt-get install -y --no-install-recommends libglpk40 libxml2 \
 && rm -rf /var/lib/apt/lists/*

RUN install2.r --error --skipinstalled \
      data.table ggplot2 pracma Matrix Rcpp igraph optparse patchwork R.utils

COPY . /app
WORKDIR /app

# Compile src/loh_hmm.cpp once at build time so every run starts warm,
# and confirm igraph loads.
RUN rm -rf /app/src/.rcpp_cache \
 && Rscript -e "source('/app/chimera_functions.R'); invisible(loadNamespace('igraph')); cat('ChimeraMapR build check OK\n')"

ENTRYPOINT ["Rscript", "/app/chimera_cli.R"]
CMD ["--help"]
