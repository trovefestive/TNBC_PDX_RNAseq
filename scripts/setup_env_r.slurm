#!/usr/bin/env bash
# =============================================================================
# setup_env_r.slurm — conda environment with R + Bioconductor for the analysis (run once).
#
#   scripts/submit.sh scripts/setup_env_r.slurm
#
# Installs R and these packages:
#   DESeq2     : differential expression
#   apeglm     : shrinks noisy fold changes (better ranking, nicer plots)
#   tximport   : reads Salmon output (Track B)
#   fgsea      : fast gene set enrichment analysis (GSEA)
#   ggplot2, ggrepel, pheatmap, patchwork : plots
#   data.table, matrixStats, R.utils, bit64 : fast tables (incl. .gz) and row statistics
# Versions are appended to metadata/software_versions.txt.
# Solving an R environment can take 10-30 minutes, hence a job.
# =============================================================================
#SBATCH --job-name=setup_r_env
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/setup_env_r_%j.log

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
source ./config.sh
module load miniforge

if conda env list | awk '{ print $1 }' | grep -qx "${R_ENV}"; then
    echo "[$(date +%T)] conda env ${R_ENV} already exists; skipping the build."
else
    echo "[$(date +%T)] Creating conda env ${R_ENV} ..."
    conda create -y -n "${R_ENV}" -c conda-forge -c bioconda \
        r-base=4.4 bioconductor-deseq2 bioconductor-apeglm bioconductor-tximport \
        bioconductor-fgsea r-ggplot2 r-ggrepel r-pheatmap r-patchwork \
        r-data.table r-matrixstats r-r.utils r-bit64
fi

set +u; source activate "${R_ENV}"; set -u
export OPENBLAS_NUM_THREADS=1
echo "[$(date +%T)] Recording versions"
{
    echo
    echo "# conda env ${R_ENV} (scripts/setup_env_r.slurm, $(date +%F))"
    Rscript -e 'cat("R:", as.character(getRversion()), "\n");
                for (p in c("DESeq2","apeglm","tximport","fgsea","ggplot2","ggrepel","pheatmap","patchwork","data.table","matrixStats"))
                    cat(p, ": ", as.character(packageVersion(p)), "\n", sep = "")'
} | tee -a "${META_DIR}/software_versions.txt"
echo "[$(date +%T)] Done."
