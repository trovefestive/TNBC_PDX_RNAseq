#!/usr/bin/env bash
# =============================================================================
# 10_get_author_counts.sh — Track A, step 1: download the authors' gene counts
#                           and the gene set collections.
#
# WHAT IT DOWNLOADS (all small, < 50 MB in total)
#   data/geo_counts/
#     GSE183187_BCM_PDX_RNAseq_RSEM_ExpectedCounts.txt.gz
#         The authors' count matrix: genes x PDX models. Mouse reads were
#         removed with Xenome, human reads counted with RSEM (hg38).
#     GSM5552651_BCM-3807_human.genes.results.txt.gz
#         BCM-3807 is missing from the matrix, but its own RSEM file is on GEO.
#     GSM5552631_BCM-0002_human.genes.results.txt.gz
#         A model that IS in the matrix, to check the two formats agree.
#   data/genesets/  MSigDB gene sets (human gene symbols):
#     h.all (Hallmark; positive control), c5.go.bp (GO biological process),
#     c2.cp.reactome, c2.cp.wikipathways
#
# WHERE TO RUN: login node is fine (small text files).
#     bash scripts/10_get_author_counts.sh
# Safe to rerun: files that already exist are not downloaded again.
# =============================================================================
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
source ./config.sh

GEO_DIR="data/geo_counts"
GS_DIR="data/genesets"
MSIGDB="2026.2.Hs"                 # MSigDB release (pinned for reproducibility)
mkdir -p "${GEO_DIR}" "${GS_DIR}"

log() { echo "[$(date +%T)] $*"; }

# get URL FOLDER -> download unless already there
get() {
    local dest="$2/$(basename "$1")"
    if [[ -s "${dest}" ]]; then log "have ${dest}"; return 0; fi
    log "downloading ${dest} ..."
    curl -sS --fail --retry 3 -o "${dest}.part" "$1"
    mv "${dest}.part" "${dest}"
}

# =============================================================================
# Step 1: the authors' counts from GEO
# =============================================================================
GEO="https://ftp.ncbi.nlm.nih.gov/geo"
get "${GEO}/series/GSE183nnn/GSE183187/suppl/GSE183187_BCM_PDX_RNAseq_RSEM_ExpectedCounts.txt.gz" "${GEO_DIR}"
get "${GEO}/samples/GSM5552nnn/GSM5552651/suppl/GSM5552651_BCM-3807_human.genes.results.txt.gz" "${GEO_DIR}"
get "${GEO}/samples/GSM5552nnn/GSM5552631/suppl/GSM5552631_BCM-0002_human.genes.results.txt.gz" "${GEO_DIR}"

# =============================================================================
# Step 2: MSigDB gene sets (.gmt = one gene set per line: name, link, genes...)
# =============================================================================
MS="https://data.broadinstitute.org/gsea-msigdb/msigdb/release/${MSIGDB}"
for c in h.all c5.go.bp c2.cp.reactome c2.cp.wikipathways; do
    get "${MS}/${c}.v${MSIGDB}.symbols.gmt" "${GS_DIR}"
done
for f in "${GS_DIR}"/*.gmt; do log "  $(basename "${f}"): $(wc -l < "${f}") gene sets"; done

log "Done. Next: Rscript scripts/11_prepare_author_counts.R"
