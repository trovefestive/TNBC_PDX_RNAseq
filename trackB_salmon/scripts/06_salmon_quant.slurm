#!/usr/bin/env bash
# =============================================================================
# 06_salmon_quant.slurm — quantify one PDX model with Salmon (SLURM array).
#
#   cd trackB_salmon && scripts/submit.sh --array=1-53%6 scripts/06_salmon_quant.slurm
#
# Task N quantifies the model on data row N of ../metadata/pdx_samples.tsv
# against the combined human + mouse + ERCC index from step 04. Every read
# goes to whichever species' transcript fits it best, so mouse stroma reads
# end up on mouse genes and can be dropped later.
#
# Options:
#   -l A               detect the library type automatically (expected: IU, unstranded)
#   --validateMappings selective alignment (more accurate than plain k-mer matching)
#   --seqBias --gcBias correct for sequence-specific and GC-content biases
#
# Special case: SRR15687248 (BCM-3807) is mostly single-end. We give Salmon its
# unpaired file plus mate 1 as single-end reads, i.e. one read per spot.
#
# Output: quant/<model>/quant.sf (+ logs, lib_format_counts.json).
# Safe to rerun: models with a finished quant.sf are skipped.
# =============================================================================
#SBATCH --job-name=salmon_quant
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=logs/06_salmon_quant_%A_%a.log

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
source ./config.sh

ROW=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${SAMPLES}")
MODEL=$(cut -f1 <<< "${ROW}"); SRR=$(cut -f3 <<< "${ROW}")
OUT="${QUANT_DIR}/${MODEL}"
echo "[$(date '+%F %T')] Task ${SLURM_ARRAY_TASK_ID}: ${MODEL} (${SRR}) on $(hostname)"

[[ -f "${SALMON_INDEX}.done" ]] || { echo "ERROR: index not finished (run step 04)"; exit 1; }
if [[ -s "${OUT}/quant.sf" ]]; then echo "${OUT}/quant.sf exists; skipping."; exit 0; fi

module load miniforge
set +u; source activate "${RNASEQ_ENV}"; set -u

R1="${FASTQ_DIR}/${SRR}_1.fastq.gz"; R2="${FASTQ_DIR}/${SRR}_2.fastq.gz"; U="${FASTQ_DIR}/${SRR}.fastq.gz"
if [[ -f "${U}" ]]; then
    echo "  mostly single-end run: using ${U} + ${R1} as single-end reads"
    READS=(-r "${U}" "${R1}"); BIAS=(--seqBias)          # GC-bias model needs pairs
else
    READS=(-1 "${R1}" -2 "${R2}"); BIAS=(--seqBias --gcBias)
fi

rm -rf "${OUT}.tmp"
salmon quant -i "${SALMON_INDEX}" -l A "${READS[@]}" \
    --validateMappings "${BIAS[@]}" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -o "${OUT}.tmp"
mv "${OUT}.tmp" "${OUT}"                 # only a finished run gets the real name

echo "  library type: $(grep -o '"expected_format": "[A-Z]*"' "${OUT}/lib_format_counts.json")"
echo "  mapping rate: $(grep -o '"percent_mapped": [0-9.]*' "${OUT}/aux_info/meta_info.json")"
echo "[$(date '+%F %T')] Done: ${OUT}/quant.sf"
