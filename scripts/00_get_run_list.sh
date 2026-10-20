#!/usr/bin/env bash
# =============================================================================
# 00_get_run_list.sh — build the sample sheet for GEO series GSE183187.
#
# WHAT IT DOES
#   1. Asks ENA (the European copy of SRA) for every sequencing run (SRR) in the
#      BioProject linked to this GEO series, with read counts, file links and
#      md5 checksums.
#   2. Asks GEO for each sample's description (PDX vs patient tumour, subtype).
#   3. Joins the two into one table:     metadata/GSE183187_runinfo.tsv
#      and writes the list of runs:      metadata/runs_to_download.txt
#
# WHY
#   GEO describes the *samples* (GSM...), but the reads live in SRA as *runs*
#   (SRR...). We need the GSM -> SRR mapping, the download links and the md5
#   checksums before we download anything.
#
# WHERE TO RUN
#   On the login node is fine: it only downloads a few small text files.
#     bash scripts/00_get_run_list.sh
#
# Safe to rerun: it overwrites the two metadata files with fresh copies.
# =============================================================================
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."     # run from the project folder
source ./config.sh
mkdir -p "${META_DIR}"

# =============================================================================
# Step 1: run table from ENA
# =============================================================================
# The ENA "filereport" API returns one row per run as a tab-separated table.
# read_count = number of spots (for paired-end data, one spot = one read pair).
echo "[$(date +%T)] Step 1: asking ENA for the runs in ${BIOPROJECT} ..."
ENA_FIELDS="run_accession,sample_alias,sample_title,library_layout,library_strategy,instrument_model,scientific_name,read_count,base_count,fastq_ftp,fastq_md5,fastq_bytes"
curl -sS --fail --retry 3 \
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${BIOPROJECT}&result=read_run&fields=${ENA_FIELDS}&format=tsv" \
    > "${META_DIR}/ena_filereport.tsv"
echo "    ENA returned $(($(wc -l < "${META_DIR}/ena_filereport.tsv") - 1)) runs"

# =============================================================================
# Step 2: sample descriptions from GEO
# =============================================================================
# GEO's "brief text" view lists every sample (GSM) with lines such as
#   !Sample_characteristics_ch1 = tumor subtype: human TNBC tumor
# GEO sends Windows line endings (\r\n); `tr -d '\r'` removes the \r.
echo "[$(date +%T)] Step 2: asking GEO for the sample descriptions of ${GSE} ..."
curl -sS --fail --retry 3 \
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${GSE}&targ=gsm&form=text&view=brief" \
    | tr -d '\r' > "${META_DIR}/geo_samples.txt"

# Turn the GEO text into a small table: GSM, sample type, tumour subtype, host, organisms.
awk '
    function flush() { if (gsm != "") print gsm "\t" type "\t" subtype "\t" host "\t" org }
    /^\^SAMPLE = /                 { flush(); gsm = $3; type = "patient_tumor"; subtype = "NA"; host = "NA"; org = ""; next }
    /^!Sample_organism_ch1 = /     { sub(/^!Sample_organism_ch1 = /, ""); org = (org == "" ? $0 : org "+" $0); next }
    /characteristics_ch1 = cell type: .*graft/ { type = "PDX"; next }
    /characteristics_ch1 = tumor subtype: /     { sub(/.*tumor subtype: human /, ""); sub(/ tumor$/, ""); subtype = $0; next }
    /characteristics_ch1 = strain: /           { sub(/.*strain: /, ""); host = $0; next }
    END { flush() }
' "${META_DIR}/geo_samples.txt" > "${META_DIR}/geo_samples.tsv"
echo "    GEO returned $(wc -l < "${META_DIR}/geo_samples.tsv") samples"

# =============================================================================
# Step 3: join ENA runs with GEO samples (key = GSM, which ENA calls sample_alias)
# =============================================================================
echo "[$(date +%T)] Step 3: writing ${RUNINFO} ..."
awk -F'\t' -v OFS='\t' '
    # first file: GEO table, remembered by GSM
    FNR == NR { type[$1] = $2; subtype[$1] = $3; host[$1] = $4; org[$1] = $5; next }
    # second file: ENA table. Column numbers follow ENA_FIELDS above.
    FNR == 1  { print "gsm", "sample_title", "srr", "sample_type", "tumor_subtype", "mouse_host",
                      "geo_organism", "layout", "spots", "reads", "bases", "instrument",
                      "n_files", "fastq_gb", "fastq_url", "fastq_md5", "fastq_bytes"; next }
    {
        gsm = $2
        n = split($12, size, ";"); gb = 0
        for (i = 1; i <= n; i++) gb += size[i] / 1e9
        reads = ($4 == "PAIRED") ? 2 * $8 : $8          # paired-end: 2 reads per spot
        print gsm, $3, $1, type[gsm], subtype[gsm], host[gsm], org[gsm], $4, $8, reads, $9, $6,
              n, sprintf("%.1f", gb), $10, $11, $12
    }
' "${META_DIR}/geo_samples.tsv" "${META_DIR}/ena_filereport.tsv" \
    | (read -r header; echo "${header}"; sort -t$'\t' -k1,1) > "${RUNINFO}"

# The run list for step 01: every run. Edit this file if you only want some.
cut -f3 "${RUNINFO}" | tail -n +2 > "${RUN_LIST}"

# =============================================================================
# Step 4: sanity checks + summary
# =============================================================================
n_runs=$(wc -l < "${RUN_LIST}")
n_missing_url=$(awk -F'\t' 'NR > 1 && $15 == ""' "${RUNINFO}" | wc -l)
n_no_geo=$(awk -F'\t' 'NR > 1 && $4 == ""' "${RUNINFO}" | wc -l)
echo
echo "Summary"
echo "  runs:                    ${n_runs}"
echo "  runs without ENA fastq:  ${n_missing_url}   (these would need the SRA-toolkit fallback)"
echo "  runs not matched to GEO: ${n_no_geo}"
echo "  by sample type / subtype / layout:"
awk -F'\t' 'NR > 1 { print "    " $4 "\t" $5 "\t" $8 }' "${RUNINFO}" | sort | uniq -c
awk -F'\t' 'NR > 1 { gb += $14 } END { printf "  total download size:     %.0f GB (%.2f TB)\n", gb, gb / 1000 }' "${RUNINFO}"
echo
echo "Next: check ${RUN_LIST}, then submit step 01 (see README.md)."
