#!/usr/bin/env Rscript
# =============================================================================
# 14_make_teaching_files.R — build the input files for the undergraduate
#                            DESeq2 exercise in teaching_deseq2/.
#
# Uses the authors' RSEM counts prepared in step 11 (Track A):
#   - 52 PDX models (BCM-3807 left out: only 1.4 M counts in the authors' output)
#   - counts rounded to whole numbers (DESeq2 needs integers)
#   - genes with zero counts in every model removed (smaller file)
# Writes:
#   teaching_deseq2/data/counts_TNBC_PDX.csv.gz   genes (rows) x models (columns)
#   teaching_deseq2/data/sample_info.csv          one row per model
#
# RUN (small; login node is fine):  Rscript scripts/14_make_teaching_files.R
# =============================================================================
suppressPackageStartupMessages(library(data.table))
dir.create("teaching_deseq2/data", recursive = TRUE, showWarnings = FALSE)

samples <- fread("metadata/pdx_samples.tsv")[model != "BCM-3807"]
counts  <- fread("results/counts/authors_rsem_expected_counts.tsv.gz")

# Undo Excel's date conversion of gene names in the authors' matrix:
#   "5-Sep" -> SEPTIN5, "3-Mar" -> MARCHF3, "1-Dec" -> DELEC1.
# "1-Mar" and "2-Mar" are dropped: Excel turned BOTH MARCH1 and MARC1 (now MTARC1)
# into "1-Mar" (same for 2), so those rows mix two different genes.
counts <- counts[!gene %in% c("1-Mar", "2-Mar")]
counts[grepl("^[0-9]+-Sep$", gene), gene := paste0("SEPTIN", sub("-Sep", "", gene))]
counts[grepl("^[0-9]+-Mar$", gene), gene := paste0("MARCHF", sub("-Mar", "", gene))]
counts[gene == "1-Dec", gene := "DELEC1"]
stopifnot(!anyDuplicated(counts$gene))

m <- round(as.matrix(counts[, samples$model, with = FALSE]))
rownames(m) <- counts$gene
m <- m[rowSums(m) > 0, ]

# Sample sheet: friendly column names, sorted so TNBC and non-TNBC are grouped
info <- samples[, .(sample = model, group = group, subtype = subtype, instrument = instrument, geo_accession = gsm)]
setorder(info, group, sample)
m <- m[, info$sample]
stopifnot(identical(colnames(m), info$sample))

fwrite(data.table(gene = rownames(m), m), "teaching_deseq2/data/counts_TNBC_PDX.csv.gz")
fwrite(info, "teaching_deseq2/data/sample_info.csv")
message("Wrote ", nrow(m), " genes x ", ncol(m), " samples; groups: ",
        paste(names(table(info$group)), table(info$group), sep = " = ", collapse = ", "))
