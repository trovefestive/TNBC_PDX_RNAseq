#!/usr/bin/env bash
# =============================================================================
# trackB_salmon/config.sh — settings for Track B (reprocess FASTQs with Salmon).
#   Loads the main project config first (sample sheet, FASTQ folder, SLURM
#   account), then adds Track B's own paths. Paths below are relative to this
#   trackB_salmon/ folder.
# =============================================================================
source ../config.sh                 # main project settings (+ config.local.sh)

# Shared inputs from the main project (one level up)
export FASTQ_DIR="../fastq"
export RUNINFO="../metadata/GSE183187_runinfo.tsv"

# Track B's own folders
export META_DIR="metadata"          # Track B tool versions, QC tables
export LOG_DIR="logs"

export RNASEQ_ENV="tnbc-rnaseq"     # salmon, fastqc, multiqc, cutadapt (scripts/setup_env_rnaseq.slurm)

# --- Reference files (step 03 downloads, step 04 builds the Salmon index) ----
# Human and mouse transcripts go in ONE index, so mouse (stroma) reads are
# assigned to mouse genes instead of being forced onto human genes.
export REF_DIR="data/ref"
export GENCODE_HUMAN="50"           # human GENCODE release (GRCh38)
export GENCODE_MOUSE="M39"          # mouse GENCODE release (GRCm39)
export SALMON_INDEX="${REF_DIR}/salmon_hs${GENCODE_HUMAN}_mm${GENCODE_MOUSE}_ercc_k31"
