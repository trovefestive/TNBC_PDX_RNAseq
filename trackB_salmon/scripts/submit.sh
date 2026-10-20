#!/usr/bin/env bash
# =============================================================================
# submit.sh — a small wrapper around `sbatch`.
#   usage:  scripts/submit.sh [extra sbatch options] scripts/NN_step.slurm
#   example: (run from trackB_salmon/)  scripts/submit.sh scripts/04_build_salmon_index.slurm
#
# Why a wrapper? It loads config.sh (+ config.local.sh), which sets
# SBATCH_ACCOUNT and SBATCH_PARTITION. sbatch picks those up from the
# environment, so no script has to hard-code your allocation name.
# =============================================================================
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."      # always run from the trackB_salmon folder
source ./config.sh
: "${SBATCH_ACCOUNT:?Please create config.local.sh (see config.local.sh.example)}"
mkdir -p "${LOG_DIR}"
sbatch "$@"
