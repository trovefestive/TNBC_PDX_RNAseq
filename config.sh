#!/usr/bin/env bash
# =============================================================================
# config.sh — settings shared by every script in this project.
#   Every script runs `source config.sh` first, so change things HERE, not in
#   each script. All paths are relative to the project folder (where this file is).
# =============================================================================

# --- Which GEO dataset we are downloading -----------------------------------
export GSE="GSE183187"            # GEO series
export BIOPROJECT="PRJNA759581"   # the BioProject/SRA study GEO links to (SRA study SRP335205)

# --- Where things go ---------------------------------------------------------
export META_DIR="metadata"        # sample sheet + run list
export FASTQ_DIR="fastq"          # downloaded *.fastq.gz files
export CHECK_DIR="fastq/.checks"  # small "this file is done" marker files written by step 01
export LOG_DIR="logs"             # SLURM log files

export RUNINFO="${META_DIR}/${GSE}_runinfo.tsv"     # full sample sheet (written by step 00)
export RUN_LIST="${META_DIR}/runs_to_download.txt"  # one SRR per line; step 01 downloads these

# --- How to download ---------------------------------------------------------
# "sra" = NCBI SRA via prefetch + fasterq-dump + pigz. NCBI keeps the data on
#         Amazon cloud storage; from our cluster this ran at ~55 MB/s (Oct 2026).
# "ena" = ready-made .fastq.gz from ENA over HTTPS, md5-checked. Simpler, but
#         it ran at only ~1.4 MB/s from here (2.3 TB would take days).
# Whichever you pick, step 01 falls back to the other one if it fails.
export DOWNLOAD_METHOD="sra"

# --- Software ----------------------------------------------------------------
export CONDA_ENV="tnbc-download"  # conda env with sra-tools + pigz (made by scripts/setup_env.slurm)
export R_ENV="tnbc-r"             # conda env with R, DESeq2, apeglm, fgsea, ... (scripts/setup_env_r.slurm)

# --- Site-specific settings (SLURM account/partition) -----------------------
# These live in config.local.sh, which is NOT tracked by git, so the scripts
# never contain an allocation name. Copy config.local.sh.example to start.
# (found next to this file, so it also works from a subfolder such as trackB_salmon/)
_CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "${_CONFIG_DIR}/config.local.sh" ]]; then
    source "${_CONFIG_DIR}/config.local.sh"
fi
