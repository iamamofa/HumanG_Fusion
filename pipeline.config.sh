#!/usr/bin/env bash
# pipeline.config.sh - set variables here

# Input
SAMPLES_DIR="/path/to/fastq"        # contains sample_R1.fastq.gz sample_R2.fastq.gz
SAMPLES_FILE="/path/to/samples.txt" # optional: tab-delim sample_id    R1    R2

# Output
OUTDIR="/path/to/outdir"
LOGDIR="${OUTDIR}/logs"
WORKDIR="${OUTDIR}/work"

# References (HUMAN ONLY)
GENOME_FASTA="/refs/GRCh38/GRCh38.primary_assembly.genome.fa"
GTF="/refs/GRCh38/gencode.v41.annotation.gtf"
STAR_INDEX="/refs/STAR_index_GRCh38"  # must be prebuilt with --sjdbGTFfile ${GTF}
STAR_FUSION_CTAT_LIB="/refs/ctat_resource_lib"  # resource lib for STAR-Fusion

# Tools (either module names, paths, or commands inside containers)
FASTQC="fastqc"
CUTADAPT="cutadapt"
STAR="STAR"
STAR_FUSION="STAR-Fusion"   # ensure in PATH or use full path
ARRIBA="arriba"
FUSIONCATCHER="fusioncatcher" # optional

# Options
THREADS=8
TRIM_ARGS="--cores ${THREADS} -q 20 -m 30"   # cutadapt example
STAR_THREADS=${THREADS}
ARRIBA_THREADS=${THREADS}

# Filters
MIN_SPLIT_READS=3
MIN_SPANNING_READS=2

# Other
FORCE=false  # set to true to overwrite outputs
