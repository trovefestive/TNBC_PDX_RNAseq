#!/usr/bin/env bash
# =============================================================================
# 01_download_fastq.slurm — download the FASTQ files of ONE sequencing run.
#
# HOW IT IS USED
#   This is a SLURM *array* job: SLURM starts many copies of this script, and
#   each copy gets its own number in $SLURM_ARRAY_TASK_ID (1, 2, 3, ...).
#   Copy number N downloads the run on line N of metadata/runs_to_download.txt.
#
#     scripts/submit.sh --array=1-68%6 scripts/01_download_fastq.slurm
#                               ^^^^ ^^
#                 tasks 1 to 68 ─┘    └─ at most 6 running at the same time
#
# WHAT ONE TASK DOES
#   1. Looks up its run (SRR) in the sample sheet: download links, md5 checksums,
#      file sizes and number of spots.
#   2. Downloads the run with the method set in config.sh (DOWNLOAD_METHOD):
#      "sra": SRA-toolkit `prefetch` fetches the compact .sra file from NCBI
#             (prefetch checks its md5 itself), `vdb-validate` double-checks it,
#             `fasterq-dump` turns it into FASTQ, `pigz` compresses it.
#      "ena": downloads ready-made .fastq.gz from ENA over HTTPS and checks
#             each file against ENA's md5 (wrong md5 = broken: delete, retry).
#      If the chosen method fails, it tries the other one.
#   3. Counts the reads in each file (must equal the number of spots) and
#      writes a small "done" marker file:   fastq/.checks/<SRR>.done
#
# SAFE TO RERUN
#   A run with a .done marker is skipped, and an ENA download that stopped
#   halfway is resumed (`wget -c`). If a task fails, just resubmit its number.
# =============================================================================
#SBATCH --job-name=dl_fastq
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/01_download_%A_%a.log
# (%A = the array job ID, %a = the task number, so every task gets its own log.)

set -euo pipefail        # stop on any error, unset variable or failed pipe
cd "${SLURM_SUBMIT_DIR}" # the project folder (submit.sh submits from there)
source ./config.sh

THREADS="${SLURM_CPUS_PER_TASK:-8}"
TASK="${SLURM_ARRAY_TASK_ID:?This script must be submitted as an array job (--array=...)}"
mkdir -p "${FASTQ_DIR}" "${CHECK_DIR}"

log() { echo "[$(date '+%F %T')] $*"; }

# =============================================================================
# Step 1: which run am I? Look it up in the sample sheet.
# =============================================================================
SRR=$(sed -n "${TASK}p" "${RUN_LIST}")
if [[ -z "${SRR}" ]]; then
    log "No line ${TASK} in ${RUN_LIST}; nothing to do."
    exit 0
fi

# Find this run's row (column 3 = srr) and pull out the columns we need.
ROW=$(awk -F'\t' -v r="${SRR}" '$3 == r' "${RUNINFO}")
if [[ -z "${ROW}" ]]; then
    log "ERROR: ${SRR} is not in ${RUNINFO}. Rerun scripts/00_get_run_list.sh?"
    exit 1
fi
GSM=$(cut -f1  <<< "${ROW}")
TITLE=$(cut -f2  <<< "${ROW}")
SPOTS=$(cut -f9  <<< "${ROW}")
URLS=$(cut -f15 <<< "${ROW}")    # several values are separated by ";"
MD5S=$(cut -f16 <<< "${ROW}")
BYTES=$(cut -f17 <<< "${ROW}")

log "Task ${TASK}: ${SRR} (${GSM}, ${TITLE}), ${SPOTS} spots expected"
log "Running on $(hostname) with ${THREADS} CPUs"

if [[ -f "${CHECK_DIR}/${SRR}.done" ]]; then
    log "${SRR} already finished (found ${CHECK_DIR}/${SRR}.done). Skipping."
    exit 0
fi

# Turn "a;b" into bash arrays: URL_ARR[0]=a, URL_ARR[1]=b
IFS=';' read -r -a URL_ARR  <<< "${URLS}"
IFS=';' read -r -a MD5_ARR  <<< "${MD5S}"
IFS=';' read -r -a SIZE_ARR <<< "${BYTES}"

# =============================================================================
# Step 2a: method "ena" — ready-made .fastq.gz from ENA, md5-checked
# =============================================================================
# download_one URL EXPECTED_MD5 DESTINATION  -> returns 0 if the file is good
download_one() {
    local url="https://$1" md5="$2" dest="$3"
    local attempt got

    # Already downloaded by an earlier run of this task? Check it, don't refetch.
    if [[ -f "${dest}" ]]; then
        log "  ${dest} exists; checking its md5 ..."
        got=$(md5sum "${dest}" | cut -d' ' -f1)
        if [[ "${got}" == "${md5}" ]]; then log "  md5 OK"; return 0; fi
        log "  md5 mismatch (${got}); downloading again"
        rm -f "${dest}"
    fi

    for attempt in 1 2 3; do
        log "  downloading $(basename "${dest}") (attempt ${attempt}) ..."
        # -c = continue a partial file; the ".part" name means "not finished yet"
        if wget -c --tries=10 --waitretry=30 --timeout=120 --progress=dot:giga \
                -O "${dest}.part" "${url}"; then
            got=$(md5sum "${dest}.part" | cut -d' ' -f1)
            if [[ "${got}" == "${md5}" ]]; then
                mv "${dest}.part" "${dest}"     # finished and correct: give it its real name
                log "  md5 OK (${md5})"
                return 0
            fi
            log "  md5 mismatch: expected ${md5}, got ${got}; deleting and retrying"
            rm -f "${dest}.part"
        else
            log "  wget failed; will retry"
        fi
        sleep 60
    done
    return 1
}

get_from_ena() {
    local i dest
    if [[ -z "${URLS}" ]]; then log "ENA has no FASTQ links for ${SRR}."; return 1; fi
    for i in "${!URL_ARR[@]}"; do
        dest="${FASTQ_DIR}/$(basename "${URL_ARR[$i]}")"
        log "File $((i + 1)) of ${#URL_ARR[@]}: ${dest} ($(awk -v b="${SIZE_ARR[$i]}" 'BEGIN{printf "%.1f GB", b/1e9}'))"
        download_one "${URL_ARR[$i]}" "${MD5_ARR[$i]}" "${dest}" || return 1
    done
}

# =============================================================================
# Step 2b: method "sra" — NCBI SRA toolkit, then compress with pigz
# =============================================================================
get_from_sra() {
    local TMP="${FASTQ_DIR}/tmp_${SRR}"     # working space; deleted at the end

    # Finished .fastq.gz files from an earlier attempt? They are only moved into
    # fastq/ after pigz finished, so they are complete: skip to counting reads.
    local f all_there=true
    for f in "${URL_ARR[@]}"; do [[ -s "${FASTQ_DIR}/$(basename "${f}")" ]] || all_there=false; done
    if [[ -n "${URLS}" && "${all_there}" == true ]]; then
        log "FASTQ files for ${SRR} already in ${FASTQ_DIR}; not downloading again."
        return 0
    fi
    rm -rf "${TMP}"; mkdir -p "${TMP}"

    module load miniforge
    set +u                                 # conda's activate script uses unset variables
    source activate "${CONDA_ENV}"         # provides prefetch, vdb-validate, fasterq-dump, pigz
    set -u

    log "prefetch ${SRR}: download the compact .sra file (prefetch checks its md5) ..."
    prefetch --max-size u --output-directory "${TMP}" "${SRR}" || return 1

    log "vdb-validate: check the .sra file is complete and not corrupted ..."
    vdb-validate "${TMP}/${SRR}" || return 1

    log "fasterq-dump: .sra -> FASTQ (--split-3 writes _1 and _2 for paired reads) ..."
    fasterq-dump --split-3 --threads "${THREADS}" --temp "${TMP}" \
        --outdir "${TMP}" "${TMP}/${SRR}" || return 1

    log "pigz: compress the FASTQ files with ${THREADS} threads ..."
    local fq
    for fq in "${TMP}/${SRR}"*.fastq; do
        pigz -p "${THREADS}" "${fq}" || return 1
        mv "${fq}.gz" "${FASTQ_DIR}/"
    done
    rm -rf "${TMP}"
}

# =============================================================================
# Step 2c: run the chosen method; if it fails, try the other one
# =============================================================================
clean_partial() { rm -rf "${FASTQ_DIR}/tmp_${SRR}" "${FASTQ_DIR}/${SRR}"_*.fastq.gz.part; }

if [[ "${DOWNLOAD_METHOD}" == "ena" ]]; then FIRST="ena"; SECOND="sra"
else                                          FIRST="sra"; SECOND="ena"; fi

METHOD="${FIRST}"
log "Downloading with method '${FIRST}' ..."
if ! "get_from_${FIRST}"; then
    log "Method '${FIRST}' failed; trying '${SECOND}' instead."
    clean_partial
    METHOD="${SECOND}"
    if ! "get_from_${SECOND}"; then
        log "ERROR: both methods failed for ${SRR}. See the messages above."
        clean_partial
        exit 1
    fi
fi

# =============================================================================
# Step 3: count reads (a FASTQ record is 4 lines) and write the "done" marker
# =============================================================================
# For paired-end data, _1 and _2 must both contain one read per spot.
log "Counting reads ..."
COUNTS=""
shopt -s nullglob        # a pattern that matches nothing expands to nothing
for fq in "${FASTQ_DIR}/${SRR}"_*.fastq.gz "${FASTQ_DIR}/${SRR}".fastq.gz; do
    [[ -f "${fq}" ]] || continue     # SRR.fastq.gz (unpaired reads) usually doesn't exist
    n=$(( $(pigz -dc -p 2 "${fq}" | wc -l) / 4 ))
    log "  $(basename "${fq}"): ${n} reads"
    [[ "${n}" -eq "${SPOTS}" ]] || log "  WARNING: expected ${SPOTS} (the number of spots in SRA)"
    COUNTS="${COUNTS}$(basename "${fq}")=${n};"
done

# The marker is tab-separated: SRR, method, finished at, expected spots, counts
printf '%s\t%s\t%s\t%s\t%s\n' "${SRR}" "${METHOD}" "$(date '+%F %T')" "${SPOTS}" "${COUNTS%;}" \
    > "${CHECK_DIR}/${SRR}.done"
log "Done: ${SRR} via ${METHOD}. Marker written to ${CHECK_DIR}/${SRR}.done"
