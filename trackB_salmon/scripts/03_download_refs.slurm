#!/usr/bin/env bash
# =============================================================================
# 03_download_refs.slurm — download the reference files for the Salmon index.
#
#   cd trackB_salmon && scripts/submit.sh scripts/03_download_refs.slurm
#
# WHAT IT DOWNLOADS (into data/ref/)
#   Human, GENCODE release 50 (GRCh38):
#     transcripts FASTA, primary-assembly genome FASTA, annotation GTF
#   Mouse, GENCODE release M39 (GRCm39):
#     transcripts FASTA, primary-assembly genome FASTA, annotation GTF
#   ERCC92: the 92 synthetic spike-in RNAs added to these libraries
#
# WHY BOTH SPECIES?
#   A PDX is a human tumour growing in a mouse. Its RNA-seq reads are partly
#   human (tumour) and partly mouse (stroma). If we gave Salmon only human
#   transcripts, mouse reads from conserved genes would be counted as human.
#   With both in one index, each read goes to the species it matches best.
#
# WHY THE GENOMES?
#   Step 04 uses them as "decoys": reads that look more like a genome region
#   than like any transcript (e.g. intronic reads) are thrown away instead of
#   being forced onto a similar transcript.
#
# Release numbers come from config.sh. Safe to rerun: existing files are
# kept, and each download is checked with `gzip -t`.
# Total about 2.5 GB to download.
# =============================================================================
#SBATCH --job-name=dl_refs
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=logs/03_download_refs_%j.log

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
source ./config.sh
mkdir -p "${REF_DIR}"

log() { echo "[$(date '+%F %T')] $*"; }

HS="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_HUMAN}"
MM="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${GENCODE_MOUSE}"

# get URL  -> downloads into data/ref/ unless the file is already there and complete
get() {
    local url="$1" dest="${REF_DIR}/$(basename "$1")"
    if [[ -s "${dest}" ]] && gzip -t "${dest}" 2>/dev/null; then
        log "have $(basename "${dest}")"
        return 0
    fi
    log "downloading $(basename "${dest}") ..."
    wget -q -c --tries=10 --waitretry=30 -O "${dest}.part" "${url}"
    gzip -t "${dest}.part"            # fails (and stops the script) if the file is broken
    mv "${dest}.part" "${dest}"
}

# =============================================================================
# Step 1: GENCODE human and mouse
# =============================================================================
get "${HS}/gencode.v${GENCODE_HUMAN}.transcripts.fa.gz"
get "${HS}/GRCh38.primary_assembly.genome.fa.gz"
get "${HS}/gencode.v${GENCODE_HUMAN}.primary_assembly.annotation.gtf.gz"
get "${MM}/gencode.v${GENCODE_MOUSE}.transcripts.fa.gz"
get "${MM}/GRCm39.primary_assembly.genome.fa.gz"
get "${MM}/gencode.v${GENCODE_MOUSE}.primary_assembly.annotation.gtf.gz"

# GENCODE publishes md5 checksums for every file: check ours against them.
log "checking md5 sums against GENCODE's MD5SUMS files ..."
for pair in "${HS}|gencode.v${GENCODE_HUMAN}.transcripts.fa.gz GRCh38.primary_assembly.genome.fa.gz gencode.v${GENCODE_HUMAN}.primary_assembly.annotation.gtf.gz" \
            "${MM}|gencode.v${GENCODE_MOUSE}.transcripts.fa.gz GRCm39.primary_assembly.genome.fa.gz gencode.v${GENCODE_MOUSE}.primary_assembly.annotation.gtf.gz"; do
    base="${pair%%|*}"; files="${pair#*|}"
    sums=$(curl -sS --fail "${base}/MD5SUMS")
    for f in ${files}; do
        want=$(awk -v f="${f}" '$2 == f { print $1 }' <<< "${sums}")
        got=$(md5sum "${REF_DIR}/${f}" | cut -d' ' -f1)
        if [[ "${want}" == "${got}" ]]; then log "  md5 OK  ${f}"
        else log "  md5 MISMATCH ${f} (expected ${want}, got ${got}); delete it and rerun"; exit 1; fi
    done
done

# =============================================================================
# Step 2: ERCC spike-in sequences (Thermo Fisher)
# =============================================================================
if [[ ! -s "${REF_DIR}/ERCC92.fa" ]]; then
    log "downloading ERCC92 ..."
    wget -q -O "${REF_DIR}/ERCC92.zip" \
        "https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip"
    # unzip may not be installed everywhere; Python can unzip too
    python3 -c "import zipfile,sys; zipfile.ZipFile(sys.argv[1]).extractall(sys.argv[2])" \
        "${REF_DIR}/ERCC92.zip" "${REF_DIR}"
    rm -f "${REF_DIR}/ERCC92.zip"
fi
log "ERCC92.fa has $(grep -c '^>' "${REF_DIR}/ERCC92.fa") sequences (expected 92)"

# =============================================================================
# Summary
# =============================================================================
log "Files in ${REF_DIR}:"
ls -lh "${REF_DIR}"
log "Done. Next: scripts/submit.sh scripts/04_build_salmon_index.slurm"
