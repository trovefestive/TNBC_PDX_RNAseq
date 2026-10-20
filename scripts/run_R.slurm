#!/usr/bin/env bash
# =============================================================================
# run_R.slurm — run one R script as a SLURM job, inside the tnbc-r conda env.
#
#   scripts/submit.sh scripts/run_R.slurm scripts/12_deseq2_trackA.R
#
# Why: R analyses can take more than a few minutes or several GB of memory,
# which is not allowed on the shared login node. Anything after the script
# name is passed on to R.
# The log is logs/R_<jobid>.log.
# =============================================================================
#SBATCH --job-name=run_R
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/R_%j.log

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
source ./config.sh
SCRIPT="${1:?usage: scripts/submit.sh scripts/run_R.slurm scripts/<script>.R}"
shift

module load miniforge
set +u; source activate "${R_ENV}"; set -u
export OPENBLAS_NUM_THREADS=1          # avoid BLAS thread over-subscription

echo "[$(date '+%F %T')] Running ${SCRIPT} on $(hostname) (R $(Rscript -e 'cat(as.character(getRversion()))'))"
Rscript "${SCRIPT}" "$@"
echo "[$(date '+%F %T')] Finished ${SCRIPT}"
