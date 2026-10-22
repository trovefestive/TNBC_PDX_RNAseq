#!/usr/bin/env Rscript
# =============================================================================
# 08_import_salmon.R — Track B: read the Salmon results, measure how much of
#                      each PDX is human vs mouse, keep human genes only.
#
# WHAT IT DOES
#   1. Per model QC from Salmon's own files: library type, mapping rate.
#   2. Species split: share of mapped reads that went to human, mouse and ERCC
#      transcripts. The mouse share is a measure of stroma content.
#   3. tximport: transcripts -> genes (all species), then keep HUMAN genes only.
#      Mouse and ERCC genes are dropped before any normalization.
#   4. Flags the human genes that have a transcript identical to a mouse one
#      (Salmon merged those, so mouse reads can land on them).
#
# OUTPUT (results/):
#   qc/salmon_qc.tsv                    one row per model
#   qc/mouse_fraction_by_group.(pdf|png)
#   counts/txi_human.rds                tximport object, human genes only (for DESeq2)
#   counts/human_gene_counts.tsv.gz     gene x model estimated counts
#   counts/human_gene_annotation.tsv    gene_id, gene_name, gene_type, mouse_identical_tx
#
# RUN (from trackB_salmon/):  scripts/submit.sh scripts/run_R.slurm scripts/08_import_salmon.R
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(tximport); library(ggplot2) })
msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

for (d in c("results/qc", "results/counts")) dir.create(d, recursive = TRUE, showWarnings = FALSE)
group_cols <- c(TNBC = "#D55E00", nonTNBC = "#0072B2")
index   <- Sys.getenv("SALMON_INDEX", "data/ref/salmon_hs50_mmM39_ercc_k31")
tx2gene <- fread(Sys.glob("data/ref/tx2gene_hs*_mm*_ercc.tsv"))

# =============================================================================
# Step 1: which models have finished Salmon output?
# =============================================================================
samples <- fread("../metadata/pdx_samples.tsv")
samples[, quant := file.path("quant", model, "quant.sf")]
missing <- samples[!file.exists(quant), model]
if (length(missing)) stop("Salmon output missing for: ", paste(missing, collapse = ", "))
msg("Step 1: ", nrow(samples), " models with quant.sf")

# Small helper: pull one number/word out of a Salmon JSON file without a JSON package
json_field <- function(file, field) {
    x <- paste(readLines(file, warn = FALSE), collapse = " ")
    sub(paste0('.*"', field, '":\\s*"?([^",}]+)"?.*'), "\\1", x)
}
samples[, lib_type     := sapply(model, function(m) json_field(file.path("quant", m, "lib_format_counts.json"), "expected_format"))]
samples[, reads_total  := as.numeric(sapply(model, function(m) json_field(file.path("quant", m, "aux_info", "meta_info.json"), "num_processed")))]
samples[, pct_mapped   := as.numeric(sapply(model, function(m) json_field(file.path("quant", m, "aux_info", "meta_info.json"), "percent_mapped")))]
samples[, pct_decoy    := 100 * as.numeric(sapply(model, function(m) json_field(file.path("quant", m, "aux_info", "meta_info.json"), "num_decoy_fragments"))) / reads_total]
print(samples[, .N, by = lib_type])

# =============================================================================
# Step 2: species split (share of assigned reads per species)
# =============================================================================
msg("Step 2: human / mouse / ERCC read shares")
species_of <- setNames(tx2gene$species, tx2gene$tx_id)
split <- rbindlist(lapply(samples$model, function(m) {
    q <- fread(file.path("quant", m, "quant.sf"), select = c("Name", "NumReads"))
    s <- tapply(q$NumReads, species_of[q$Name], sum)
    data.table(model = m, reads_human = s[["human"]], reads_mouse = s[["mouse"]], reads_ERCC = s[["ERCC"]])
}))
split[, assigned := reads_human + reads_mouse + reads_ERCC]
split[, `:=`(pct_human = 100 * reads_human / assigned,
             pct_mouse = 100 * reads_mouse / assigned,
             pct_ERCC  = 100 * reads_ERCC  / assigned)]
qc <- merge(samples[, .(model, srr, subtype, group, instrument, lib_type, reads_total, pct_mapped, pct_decoy)],
            split, by = "model")
fwrite(qc[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)], "results/qc/salmon_qc.tsv", sep = "\t")
print(qc[, .(median_pct_mapped = median(pct_mapped), median_pct_human = median(pct_human),
             median_pct_mouse = median(pct_mouse), median_pct_ERCC = median(pct_ERCC)), by = group])

# Does the mouse (stroma) share differ between groups? Linear model with instrument.
fit <- summary(lm(pct_mouse ~ instrument + group, data = qc))$coefficients
msg("  mouse share, TNBC vs nonTNBC (adjusted for instrument): ",
    round(fit["groupTNBC", "Estimate"], 2), " percentage points, p = ", signif(fit["groupTNBC", "Pr(>|t|)"], 3))

p <- ggplot(qc, aes(group, pct_mouse, colour = group)) +
    geom_boxplot(outlier.shape = NA, colour = "grey40") + geom_jitter(aes(shape = instrument), width = 0.15, size = 2) +
    scale_colour_manual(values = group_cols, guide = "none") + theme_classic(base_size = 11) +
    labs(x = NULL, y = "% of assigned reads that are mouse", title = "Mouse (stroma) read share per PDX model")
for (ext in c("pdf", "png")) ggsave(paste0("results/qc/mouse_fraction_by_group.", ext), p, width = 5, height = 4.5, dpi = 300)

# =============================================================================
# Step 3: tximport, then keep human genes only
# =============================================================================
msg("Step 3: tximport (transcripts -> genes)")
files <- setNames(samples$quant, samples$model)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene[, .(tx_id, gene_id)], countsFromAbundance = "no")

genes <- unique(tx2gene[, .(gene_id, gene_name, gene_type, species)])
human_ids <- intersect(rownames(txi$counts), genes[species == "human", gene_id])
txi_h <- txi
for (slot in c("abundance", "counts", "length")) txi_h[[slot]] <- txi[[slot]][human_ids, ]
msg("  genes: ", nrow(txi$counts), " total -> ", length(human_ids), " human")

# Human genes with a transcript identical to a mouse transcript (merged by Salmon)
dup <- fread(file.path(index, "duplicate_clusters.tsv"))
cross <- dup[grepl("^ENST", RetainedRef) & grepl("^ENSMUST", DuplicateRef)]
flag_genes <- unique(tx2gene[tx_id %in% cross$RetainedRef, gene_id])
ann <- genes[gene_id %in% human_ids, .(gene_id, gene_name, gene_type)]
ann[, mouse_identical_tx := gene_id %in% flag_genes]
msg("  human genes with a mouse-identical transcript (flagged): ", sum(ann$mouse_identical_tx))

saveRDS(txi_h, "results/counts/txi_human.rds")
fwrite(ann, "results/counts/human_gene_annotation.tsv", sep = "\t")
fwrite(data.table(gene_id = human_ids, round(txi_h$counts, 2)), "results/counts/human_gene_counts.tsv.gz", sep = "\t")
msg("Done.")
