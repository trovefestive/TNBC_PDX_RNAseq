#!/usr/bin/env bash
# =============================================================================
# 07_multiqc.slurm — one combined QC web page for FastQC + Salmon.
#
#   cd trackB_salmon && scripts/submit.sh scripts/07_multiqc.slurm
#
# Output: qc/multiqc_report.html (open it in a web browser).
# =============================================================================
#SBATCH --job-name=multiqc
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/07_multiqc_%j.log

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
source ./config.sh
module load miniforge
set +u; source activate "${RNASEQ_ENV}"; set -u
multiqc --force --outdir qc --filename multiqc_report.html "${FASTQC_DIR}" "${QUANT_DIR}"
echo "[$(date +%T)] Wrote qc/multiqc_report.html"
