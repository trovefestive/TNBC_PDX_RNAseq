#!/usr/bin/env bash
# =============================================================================
# submit.sh — a small wrapper around `sbatch`.
#   usage:  scripts/submit.sh [extra sbatch options] scripts/NN_step.slurm
#   example: scripts/submit.sh --array=1-68%6 scripts/01_download_fastq.slurm
#
# Why a wrapper? It loads config.sh (+ config.local.sh), which sets
# SBATCH_ACCOUNT and SBATCH_PARTITION. sbatch picks those up from the
# environment, so no script has to hard-code your allocation name.
# =============================================================================
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."      # always run from the project folder
source ./config.sh
: "${SBATCH_ACCOUNT:?Please create config.local.sh (see config.local.sh.example)}"
mkdir -p "${LOG_DIR}"
sbatch "$@"
