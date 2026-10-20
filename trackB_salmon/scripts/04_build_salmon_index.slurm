#!/usr/bin/env bash
# =============================================================================
# 04_build_salmon_index.slurm — build ONE Salmon index for human + mouse + ERCC.
#
#   cd trackB_salmon && scripts/submit.sh scripts/04_build_salmon_index.slurm
#   (needs step 03's files and the tnbc-rnaseq conda env)
#
# WHAT IT DOES
#   1. Makes the "gentrome": all transcripts (human, mouse, ERCC) followed by
#      both genomes. The genomes act as decoys (see 03_download_refs.slurm).
#      Human and mouse chromosomes have the same names (chr1, chr2, ...), so
#      we rename them hs_chr1 ... and mm_chr1 ... to keep them apart.
#   2. Writes decoys.txt, the list of decoy sequence names.
#   3. Writes tx2gene: one line per transcript, giving its gene and species.
#      tximport needs this later to add transcripts up into genes, and we use
#      the species column to keep human genes only.
#   4. Runs `salmon index` (k = 31, the usual choice for reads of ~100 bp).
#
# RESOURCES: indexing two genomes as decoys needs a lot of memory. We ask for
# 16 CPUs and 128 GB; it usually takes a few hours.
# Safe to rerun: skipped if the index is already finished (marker file .done).
# =============================================================================
#SBATCH --job-name=salmon_index
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=logs/04_salmon_index_%j.log

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
source ./config.sh
THREADS="${SLURM_CPUS_PER_TASK:-16}"
log() { echo "[$(date '+%F %T')] $*"; }

if [[ -f "${SALMON_INDEX}.done" ]]; then
    log "Index ${SALMON_INDEX} already built (found ${SALMON_INDEX}.done). Nothing to do."
    exit 0
fi

HS_TX="${REF_DIR}/gencode.v${GENCODE_HUMAN}.transcripts.fa.gz"
MM_TX="${REF_DIR}/gencode.v${GENCODE_MOUSE}.transcripts.fa.gz"
HS_GENOME="${REF_DIR}/GRCh38.primary_assembly.genome.fa.gz"
MM_GENOME="${REF_DIR}/GRCm39.primary_assembly.genome.fa.gz"
ERCC="${REF_DIR}/ERCC92.fa"
for f in "${HS_TX}" "${MM_TX}" "${HS_GENOME}" "${MM_GENOME}" "${ERCC}"; do
    [[ -s "${f}" ]] || { log "ERROR: missing ${f}. Run 03_download_refs.slurm first."; exit 1; }
done

GENTROME="${REF_DIR}/gentrome_hs${GENCODE_HUMAN}_mm${GENCODE_MOUSE}_ercc.fa.gz"
DECOYS="${REF_DIR}/decoys_hs_mm.txt"
TX2GENE="${REF_DIR}/tx2gene_hs${GENCODE_HUMAN}_mm${GENCODE_MOUSE}_ercc.tsv"

# =============================================================================
# Step 1: gentrome = transcripts first, then the (renamed) genomes
# =============================================================================
# A genome header looks like ">chr1 1"; we keep the first word and add a prefix.
log "Step 1: writing ${GENTROME} ..."
{
    zcat "${HS_TX}"
    zcat "${MM_TX}"
    tr -d '\r' < "${ERCC}"          # remove Windows line endings, if any
    zcat "${HS_GENOME}" | sed -E 's/^>([^ ]+).*/>hs_\1/'
    zcat "${MM_GENOME}" | sed -E 's/^>([^ ]+).*/>mm_\1/'
} | pigz -p "${THREADS}" > "${GENTROME}"

# =============================================================================
# Step 2: decoy names = every genome sequence name (with its prefix)
# =============================================================================
log "Step 2: writing ${DECOYS} ..."
{
    zcat "${HS_GENOME}" | grep '^>' | sed -E 's/^>([^ ]+).*/hs_\1/'
    zcat "${MM_GENOME}" | grep '^>' | sed -E 's/^>([^ ]+).*/mm_\1/'
} > "${DECOYS}"
log "  $(wc -l < "${DECOYS}") decoy sequences"

# =============================================================================
# Step 3: tx2gene from the FASTA headers
# =============================================================================
# A GENCODE transcript header is "|"-separated:
#   >ENST00000456328.2|ENSG00000290825.1|-|-|DDX11L2-202|DDX11L2|1657|lncRNA|
#     1 transcript ID   2 gene ID              5 tx name 6 gene name 8 biotype
# ERCC spike-ins are their own "gene".
log "Step 3: writing ${TX2GENE} ..."
{
    printf 'tx_id\tgene_id\tgene_name\tgene_type\tspecies\n'
    zcat "${HS_TX}" | awk -F'|' -v OFS='\t' '/^>/ { sub(/^>/, "", $1); print $1, $2, $6, $8, "human" }'
    zcat "${MM_TX}" | awk -F'|' -v OFS='\t' '/^>/ { sub(/^>/, "", $1); print $1, $2, $6, $8, "mouse" }'
    tr -d '\r' < "${ERCC}" | awk -v OFS='\t' '/^>/ { id = substr($1, 2); print id, id, id, "spike_in", "ERCC" }'
} > "${TX2GENE}"
log "  transcripts per species:"
awk -F'\t' 'NR > 1 { n[$5]++ } END { for (s in n) print "    " s ": " n[s] }' "${TX2GENE}"

# =============================================================================
# Step 4: salmon index
# =============================================================================
# --gencode : cut transcript names at the first "|" so they become just ENST...
# -k 31     : k-mer length; 31 suits reads of 75 bp or longer
module load miniforge
set +u; source activate "${RNASEQ_ENV}"; set -u
log "Step 4: salmon $(salmon --version | awk '{ print $2 }') index with ${THREADS} threads ..."
rm -rf "${SALMON_INDEX}"
salmon index \
    --transcripts "${GENTROME}" \
    --decoys "${DECOYS}" \
    --index "${SALMON_INDEX}" \
    --kmerLen 31 \
    --gencode \
    --threads "${THREADS}"

# The index is finished: mark it, and delete the big temporary gentrome file.
date '+%F %T' > "${SALMON_INDEX}.done"
rm -f "${GENTROME}"
log "Index size: $(du -sh "${SALMON_INDEX}" | cut -f1)"
if [[ -f "${SALMON_INDEX}/duplicate_clusters.tsv" ]]; then
    log "Transcripts with identical sequences that Salmon merged: $(($(wc -l < "${SALMON_INDEX}/duplicate_clusters.tsv") - 1))"
fi
log "Done. Index: ${SALMON_INDEX}"
