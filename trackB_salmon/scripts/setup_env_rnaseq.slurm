#!/usr/bin/env bash
# =============================================================================
# setup_env_rnaseq.slurm — conda environment for read QC and quantification (run once).
#
#   cd trackB_salmon && scripts/submit.sh scripts/setup_env_rnaseq.slurm
#
# Installs (bioconda / conda-forge):
#   salmon   : builds the transcript index and quantifies reads (steps 04 and 06)
#   fastqc   : per-file read quality report
#   multiqc  : combines all FastQC / Salmon reports into one web page
#   cutadapt : adapter trimming, only if FastQC shows adapters
# Versions are appended to metadata/software_versions.txt.
# Safe to rerun: skips the build if the environment already exists.
# =============================================================================
#SBATCH --job-name=setup_rnaseq_env
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/setup_env_rnaseq_%j.log

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
source ./config.sh
module load miniforge

if conda env list | awk '{ print $1 }' | grep -qx "${RNASEQ_ENV}"; then
    echo "[$(date +%T)] conda env ${RNASEQ_ENV} already exists; skipping the build."
else
    echo "[$(date +%T)] Creating conda env ${RNASEQ_ENV} ..."
    conda create -y -n "${RNASEQ_ENV}" -c conda-forge -c bioconda salmon fastqc multiqc cutadapt
fi

set +u; source activate "${RNASEQ_ENV}"; set -u
echo "[$(date +%T)] Recording versions"
{
    echo
    echo "# conda env ${RNASEQ_ENV} (scripts/setup_env_rnaseq.slurm, $(date +%F))"
    echo "salmon:   $(salmon --version 2>&1 | awk '{ print $2 }')"
    echo "fastqc:   $(fastqc --version 2>&1 | awk '{ print $2 }')"
    echo "multiqc:  $(multiqc --version 2>&1 | awk '{ print $NF }')"
    echo "cutadapt: $(cutadapt --version 2>&1)"
} | tee -a "${META_DIR}/software_versions.txt"
echo "[$(date +%T)] Done."
