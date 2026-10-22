#!/usr/bin/env bash
# =============================================================================
# 05_fastqc.slurm — read quality report for each PDX model (SLURM array).
#
#   cd trackB_salmon && scripts/submit.sh --array=1-53%10 scripts/05_fastqc.slurm
#
# Task N runs FastQC on the FASTQ files of the model on data row N of
# ../metadata/pdx_samples.tsv. Look especially at "Adapter Content": if adapters
# are above ~5%, trim with cutadapt before step 06.
# Safe to rerun: files that already have a report are skipped.
# =============================================================================
#SBATCH --job-name=fastqc
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=logs/05_fastqc_%A_%a.log

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
source ./config.sh
mkdir -p "${FASTQC_DIR}"

# Row N of the sample sheet (row 1 is the header, so take line N + 1)
ROW=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${SAMPLES}")
MODEL=$(cut -f1 <<< "${ROW}"); SRR=$(cut -f3 <<< "${ROW}")
echo "[$(date +%T)] Task ${SLURM_ARRAY_TASK_ID}: ${MODEL} (${SRR})"

module load miniforge
set +u; source activate "${RNASEQ_ENV}"; set -u

shopt -s nullglob
for fq in "${FASTQ_DIR}/${SRR}"_1.fastq.gz "${FASTQ_DIR}/${SRR}"_2.fastq.gz "${FASTQ_DIR}/${SRR}".fastq.gz; do
    [[ -f "${fq}" ]] || continue
    report="${FASTQC_DIR}/$(basename "${fq}" .fastq.gz)_fastqc.html"
    if [[ -s "${report}" ]]; then echo "  have ${report}"; continue; fi
    echo "[$(date +%T)]   fastqc $(basename "${fq}") ..."
    fastqc --threads "${SLURM_CPUS_PER_TASK}" --outdir "${FASTQC_DIR}" --quiet "${fq}"
done
echo "[$(date +%T)] Done."
