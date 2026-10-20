#!/usr/bin/env bash
# =============================================================================
# setup_env.slurm — make the conda environment with the download tools (run once).
#
#   scripts/submit.sh scripts/setup_env.slurm
#
# Installs (from the bioconda / conda-forge channels):
#   sra-tools : prefetch + fasterq-dump, the fallback downloader in step 01
#   pigz      : parallel gzip, to compress / decompress FASTQ quickly
# and writes the installed versions to metadata/software_versions.txt.
# Building an environment takes a few minutes, so it runs as a job, not on the
# login node. Safe to rerun: it skips the build if the environment exists.
# =============================================================================
#SBATCH --job-name=setup_env
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/setup_env_%j.log

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
source ./config.sh
module load miniforge

if conda env list | awk '{ print $1 }' | grep -qx "${CONDA_ENV}"; then
    echo "[$(date +%T)] conda env ${CONDA_ENV} already exists; skipping the build."
else
    echo "[$(date +%T)] Creating conda env ${CONDA_ENV} ..."
    conda create -y -n "${CONDA_ENV}" -c conda-forge -c bioconda sra-tools pigz
fi

set +u                       # conda's activate script uses unset variables
source activate "${CONDA_ENV}"
set -u

echo "[$(date +%T)] Recording versions in ${META_DIR}/software_versions.txt"
{
    echo "# Software versions (written by scripts/setup_env.slurm on $(date +%F))"
    echo "conda env: ${CONDA_ENV}"
    echo "sra-tools (prefetch):     $(prefetch --version | grep -Eo '[0-9]+\.[0-9]+\.[0-9]+' | head -1)"
    echo "sra-tools (fasterq-dump): $(fasterq-dump --version | grep -Eo '[0-9]+\.[0-9]+\.[0-9]+' | head -1)"
    echo "pigz:                     $(pigz --version 2>&1 | awk '{ print $2 }')"
    echo "wget (system):            $(wget --version | head -1 | awk '{ print $3 }')"
    echo "curl (system):            $(curl --version | head -1 | awk '{ print $2 }')"
    echo "md5sum (system):          $(md5sum --version | head -1 | awk '{ print $NF }')"
} > "${META_DIR}/software_versions.txt"
cat "${META_DIR}/software_versions.txt"
echo "[$(date +%T)] Done."
