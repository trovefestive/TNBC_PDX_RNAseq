#!/usr/bin/env Rscript
# =============================================================================
# 11_prepare_author_counts.R — Track A, step 2: sample sheet + count matrix.
#
# WHAT IT DOES
#   1. Sample sheet for the 53 PDX models: GSM, SRR, model, subtype, group
#      (TNBC / nonTNBC), sequencing instrument  -> metadata/pdx_samples.tsv
#   2. Reads the authors' RSEM count matrix (93 models) and keeps our 53:
#      52 come from the matrix, BCM-3807 from its own RSEM file.
#   3. Checks that a per-sample RSEM file gives the same numbers as the matrix
#      (BCM-0002 is in both), so adding BCM-3807 from its own file is fair.
#   4. Writes results/counts/authors_rsem_expected_counts.tsv.gz
#      (genes x 53 models; RSEM "expected counts", which are not whole numbers).
#
# RUN (small; login node is fine):
#   module load miniforge && source activate tnbc-r
#   Rscript scripts/11_prepare_author_counts.R
# =============================================================================
suppressPackageStartupMessages(library(data.table))
msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

geo_dir <- "data/geo_counts"
dir.create("results/counts", recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Step 1: sample sheet
# =============================================================================
msg("Step 1: sample sheet")
ri <- fread("metadata/GSE183187_runinfo.tsv")
pdx <- ri[sample_type == "PDX"]
samples <- data.table(
    model      = sub(".*model ", "", pdx$sample_title),     # "BCM-0002"
    gsm        = pdx$gsm,
    srr        = pdx$srr,
    subtype    = pdx$tumor_subtype,                          # TNBC, ER+, HER2+, ER+/HER2+
    group      = ifelse(pdx$tumor_subtype == "TNBC", "TNBC", "nonTNBC"),
    instrument = gsub("Illumina | ", "", pdx$instrument),   # HiSeq2000, HiSeq2500, NovaSeq6000
    layout     = pdx$layout,
    spots      = pdx$spots
)
setorder(samples, group, model)
stopifnot(nrow(samples) == 53, !anyDuplicated(samples$model))
print(samples[, .N, by = .(instrument, group)][order(instrument, group)])

# =============================================================================
# Step 2: the authors' matrix
# =============================================================================
msg("Step 2: authors' RSEM matrix")
mat <- fread(file.path(geo_dir, "GSE183187_BCM_PDX_RNAseq_RSEM_ExpectedCounts.txt.gz"))
setnames(mat, 1, "gene")
msg("  matrix: ", nrow(mat), " genes x ", ncol(mat) - 1, " models")

# Undo Excel's date conversion of gene names (same rules as script 14):
#   "5-Sep" -> SEPTIN5, "3-Mar" -> MARCHF3, "1-Dec" -> DELEC1.
# "1-Mar" and "2-Mar" are dropped: Excel turned BOTH MARCH1 and MARC1 (now MTARC1)
# into "1-Mar" (same for 2), so those rows mix two different genes.
# The per-sample RSEM files were never opened in Excel and use the old symbols
# (MARCH5, SEPT5, MARC1); fix_symbols() maps those to the same current names so
# that per-sample files and the matrix match by gene.
fix_symbols <- function(g) {
    g <- sub("^([0-9]+)-Sep$", "SEPTIN\\1", g)
    g <- sub("^([0-9]+)-Mar$", "MARCHF\\1", g)
    g <- sub("^SEPT([0-9]+)$",  "SEPTIN\\1", g)
    g <- sub("^MARCH([0-9]+)$", "MARCHF\\1", g)
    g <- sub("^MARC([12])$",    "MTARC\\1",  g)
    g[g == "1-Dec"] <- "DELEC1"
    g
}
ambiguous <- c("1-Mar", "2-Mar")
msg("  dropping ambiguous Excel-mangled rows: ", paste(mat[gene %in% ambiguous, gene], collapse = ", "))
mat <- mat[!gene %in% ambiguous]
mangled <- grepl("^[0-9]+-(Sep|Mar|Dec)$", mat$gene)
msg("  renaming ", sum(mangled), " Excel-mangled symbols (SEPTIN*, MARCHF*, DELEC1)")
mat[, gene := fix_symbols(gene)]
stopifnot(!any(grepl("^[0-9]+-[A-Z][a-z]{2}$", mat$gene)))

# Any symbol still duplicated would be summed (none expected after the fix above).
dups <- mat[duplicated(gene), unique(gene)]
msg("  duplicated gene symbols (summed): ", if (length(dups)) paste(dups, collapse = ", ") else "none")
mat <- mat[, lapply(.SD, sum), by = gene]

samples[, in_author_matrix := model %in% names(mat)]
msg("  our models in the matrix: ", sum(samples$in_author_matrix), " / 53; missing: ",
    paste(samples[in_author_matrix == FALSE, model], collapse = ", "))

# =============================================================================
# Step 3: per-sample RSEM files -- same numbers as the matrix?
# =============================================================================
read_rsem <- function(f) {
    x <- fread(f, select = c("gene_id", "expected_count"))
    x[, gene_id := fix_symbols(gene_id)]
    x[, .(expected_count = sum(expected_count)), by = .(gene = gene_id)]
}
msg("Step 3: checking a per-sample file (BCM-0002) against the matrix")
chk <- merge(read_rsem(file.path(geo_dir, "GSM5552631_BCM-0002_human.genes.results.txt.gz")),
             mat[, .(gene, matrix = `BCM-0002`)], by = "gene", all = TRUE)
msg("  genes: file ", sum(!is.na(chk$expected_count)), ", matrix ", sum(!is.na(chk$matrix)),
    ", in both ", sum(complete.cases(chk)))
both <- chk[complete.cases(chk)]
max_diff <- max(abs(both$expected_count - both$matrix))
msg("  largest absolute difference: ", signif(max_diff, 3),
    "; correlation: ", signif(cor(both$expected_count, both$matrix), 6))
if (max_diff > 0.01) stop("Per-sample file and matrix disagree: do not merge them.")

# =============================================================================
# Step 4: add BCM-3807 and write the 53-model matrix
# =============================================================================
msg("Step 4: adding BCM-3807 from its own RSEM file")
b3807 <- read_rsem(file.path(geo_dir, "GSM5552651_BCM-3807_human.genes.results.txt.gz"))
setnames(b3807, "expected_count", "BCM-3807")
msg("  BCM-3807 genes: ", nrow(b3807), "; not in matrix: ", sum(!b3807$gene %in% mat$gene),
    "; matrix genes missing from it: ", sum(!mat$gene %in% b3807$gene))
counts <- merge(mat[, c("gene", intersect(samples$model, names(mat))), with = FALSE],
                b3807, by = "gene", all.x = TRUE)
counts[is.na(`BCM-3807`), `BCM-3807` := 0]    # gene absent from its file = 0 counts
setcolorder(counts, c("gene", samples$model))
stopifnot(identical(names(counts)[-1], samples$model))

fwrite(samples, "metadata/pdx_samples.tsv", sep = "\t")
fwrite(counts, "results/counts/authors_rsem_expected_counts.tsv.gz", sep = "\t")
msg("Wrote metadata/pdx_samples.tsv (", nrow(samples), " models) and ",
    "results/counts/authors_rsem_expected_counts.tsv.gz (", nrow(counts), " genes x ",
    ncol(counts) - 1, " models)")

# Library sizes, a first sanity check
libs <- data.table(model = samples$model, group = samples$group, instrument = samples$instrument,
                   million_counts = round(colSums(counts[, -1]) / 1e6, 1))
print(libs[order(million_counts)][c(1:5, (.N - 2):.N)])
msg("Median counts per model: ", median(libs$million_counts), " M")
