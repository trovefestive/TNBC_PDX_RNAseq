#!/usr/bin/env bash
# =============================================================================
# 02_check_downloads.sh — did every run download completely and correctly?
#
# For every run in metadata/runs_to_download.txt it checks:
#   1. step 01 finished it      (the marker fastq/.checks/<SRR>.done exists)
#   2. every expected file exists, and its size matches ENA's file size
#   3. integrity was checked    (ena: step 01 matched ENA's md5 before writing the
#                                marker, and --md5 recomputes it; sra: prefetch
#                                and vdb-validate checked the .sra file)
#   4. read counts: _1 and _2 are equal, and pairs (+ unpaired reads, if the
#      submitter deposited any) add up to the number of spots SRA reports
# Results go to metadata/download_check.tsv; the last line says PASS or FAIL.
#
# TWO WAYS TO RUN IT
#   Quick check (seconds; login node is fine; reads only small files):
#       bash scripts/02_check_downloads.sh
#   Full check that recomputes every md5 (hours for ~2 TB, so use SLURM):
#       scripts/submit.sh scripts/02_check_downloads.sh --md5
#   The #SBATCH lines below are just comments to bash, but sbatch reads them.
# =============================================================================
#SBATCH --job-name=check_fastq
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=logs/02_check_%j.log

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")/..}"
source ./config.sh

RECHECK_MD5=false
if [[ "${1:-}" == "--md5" ]]; then
    RECHECK_MD5=true
    if [[ -z "${SLURM_JOB_ID:-}" ]]; then
        echo "--md5 reads every file (~2 TB). Please submit it with sbatch:"
        echo "    scripts/submit.sh scripts/02_check_downloads.sh --md5"
        exit 1
    fi
fi

OUT="${META_DIR}/download_check.tsv"
echo -e "srr\tgsm\tstatus\tmethod\tfiles_ok\tmd5\tspots_expected\tread_counts\tproblem" > "${OUT}"
echo "[$(date +%T)] Checking $(wc -l < "${RUN_LIST}") runs (recompute md5: ${RECHECK_MD5}) ..."

# md5 recompute helper: runs md5sum on several files in parallel (one per CPU)
if [[ "${RECHECK_MD5}" == true ]]; then
    echo "[$(date +%T)] Recomputing md5 of all fastq.gz files with ${SLURM_CPUS_PER_TASK:-8} parallel workers ..."
    ls "${FASTQ_DIR}"/*.fastq.gz | xargs -P "${SLURM_CPUS_PER_TASK:-8}" -n 1 md5sum \
        > "${META_DIR}/fastq_md5_recomputed.txt"
fi

n_pass=0; n_fail=0
while read -r SRR; do
    ROW=$(awk -F'\t' -v r="${SRR}" '$3 == r' "${RUNINFO}")
    GSM=$(cut -f1 <<< "${ROW}");  SPOTS=$(cut -f9 <<< "${ROW}")
    IFS=';' read -r -a URL_ARR  <<< "$(cut -f15 <<< "${ROW}")"
    IFS=';' read -r -a MD5_ARR  <<< "$(cut -f16 <<< "${ROW}")"
    IFS=';' read -r -a SIZE_ARR <<< "$(cut -f17 <<< "${ROW}")"
    problem=""; method="NA"; counts="NA"; md5="not_checked"; files_ok=0

    # --- 1. marker from step 01 -------------------------------------------
    MARKER="${CHECK_DIR}/${SRR}.done"
    if [[ -f "${MARKER}" ]]; then
        method=$(cut -f2 "${MARKER}")
        counts=$(cut -f5 "${MARKER}")
        # ena: step 01 checked each file's md5; sra: prefetch + vdb-validate checked the .sra
        if [[ "${method}" == "sra" ]]; then md5="sra_validated"; else md5="ok_in_step01"; fi
    else
        problem+="not finished (no .done marker); "
    fi

    # --- 2. files exist with the right size; 3. md5 -------------------------
    for i in "${!URL_ARR[@]}"; do
        f="${FASTQ_DIR}/$(basename "${URL_ARR[$i]}")"
        if [[ ! -s "${f}" ]]; then
            problem+="missing $(basename "${f}"); "; continue
        fi
        # files made from the .sra by fasterq-dump + pigz differ in size from ENA's: skip size check
        if [[ "${method}" != "sra" && $(stat -c %s "${f}") -ne "${SIZE_ARR[$i]}" ]]; then
            problem+="size differs for $(basename "${f}"); "; continue
        fi
        if [[ "${RECHECK_MD5}" == true && "${method}" != "sra" ]]; then
            got=$(awk -v f="${f}" '$2 == f { print $1 }' "${META_DIR}/fastq_md5_recomputed.txt")
            if [[ "${got}" != "${MD5_ARR[$i]}" ]]; then
                problem+="md5 differs for $(basename "${f}"); "; md5="FAIL"; continue
            fi
            md5="ok_recomputed"
        fi
        files_ok=$((files_ok + 1))
    done

    # --- 4. read counts must add up to the number of spots ------------------
    # Normal paired run: _1 and _2 each hold one read per spot.
    # Some runs were deposited partly unpaired (SRR.fastq.gz = spots without a
    # mate). Then: _1 must equal _2, and pairs + unpaired reads = spots.
    if [[ "${counts}" != "NA" ]]; then
        n1=0; n2=0; nsingle=0
        for pair in ${counts//;/ }; do               # e.g. SRR123_1.fastq.gz=2000
            case "${pair%%=*}" in
                *_1.fastq.gz) n1="${pair#*=}" ;;
                *_2.fastq.gz) n2="${pair#*=}" ;;
                *)            nsingle="${pair#*=}" ;;
            esac
        done
        if [[ "${n1}" -ne "${n2}" ]]; then
            problem+="_1 (${n1}) and _2 (${n2}) differ; "
        elif [[ $((n1 + nsingle)) -ne "${SPOTS}" ]]; then
            problem+="pairs + unpaired = $((n1 + nsingle)), expected ${SPOTS}; "
        fi
    fi

    if [[ -z "${problem}" ]]; then status="PASS"; n_pass=$((n_pass + 1))
    else                           status="FAIL"; n_fail=$((n_fail + 1)); fi
    echo -e "${SRR}\t${GSM}\t${status}\t${method}\t${files_ok}/${#URL_ARR[@]}\t${md5}\t${SPOTS}\t${counts}\t${problem%; }" >> "${OUT}"
done < "${RUN_LIST}"

echo
echo "Results written to ${OUT}"
if [[ "${n_fail}" -gt 0 ]]; then
    echo "Runs that FAILED:"
    awk -F'\t' '$3 == "FAIL" { print "  " $1 ": " $9 }' "${OUT}"
    echo "To redo one, find its line number with:  grep -n SRRxxxxx ${RUN_LIST}"
    echo "then:  scripts/submit.sh --array=<line number> scripts/01_download_fastq.slurm"
fi
du -sh "${FASTQ_DIR}" 2>/dev/null | awk '{ print "Space used by " $2 ": " $1 }'
echo "PASS: ${n_pass}   FAIL: ${n_fail}   => $([[ ${n_fail} -eq 0 ]] && echo PASS || echo FAIL)"
