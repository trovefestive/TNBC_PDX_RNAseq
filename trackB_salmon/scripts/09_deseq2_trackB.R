#!/usr/bin/env Rscript
# =============================================================================
# 09_deseq2_trackB.R — Track B: TNBC vs non-TNBC with DESeq2 on Salmon counts
#                      (human genes only), then comparison with Track A.
#
# SIGN CONVENTION: non-TNBC is the reference; log2FoldChange = TNBC - nonTNBC.
#
# Same model as Track A (~ instrument + group, apeglm shrinkage, >=10 counts in
# >=3 models, BH FDR < 0.05). All 53 models are used: unlike the authors'
# output, Salmon quantified all reads of BCM-3807.
#
# OUTPUT (results/):
#   de/DE_TNBC_vs_nonTNBC_instrument_adjusted.tsv   (gene_id + gene_name)
#   de/sensitivity_summary.tsv
#   de/concordance_with_trackA.tsv, plots/trackA_vs_trackB_log2FC.(pdf|png)
#   plots/volcano_TNBC_vs_nonTNBC.(pdf|png), plots/pca_subtype_instrument.(pdf|png)
#
# RUN (from trackB_salmon/):  scripts/submit.sh scripts/run_R.slurm scripts/09_deseq2_trackB.R
# =============================================================================
suppressPackageStartupMessages({
    library(data.table); library(DESeq2); library(ggplot2); library(ggrepel); library(matrixStats)
})
msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)
set.seed(1)
for (d in c("results/de", "results/plots")) dir.create(d, recursive = TRUE, showWarnings = FALSE)
group_cols   <- c(TNBC = "#D55E00", nonTNBC = "#0072B2")
subtype_cols <- c(TNBC = "#D55E00", "ER+" = "#0072B2", "HER2+" = "#009E73", "ER+/HER2+" = "#CC79A7")
theme_set(theme_classic(base_size = 11))
save_plot <- function(p, name, w = 6, h = 5) for (ext in c("pdf", "png"))
    ggsave(file.path("results/plots", paste0(name, ".", ext)), p, width = w, height = h, dpi = 300)

# =============================================================================
# Step 1: load tximport (human genes) and sample sheet
# =============================================================================
msg("Step 1: loading Salmon gene counts (human only)")
txi <- readRDS("results/counts/txi_human.rds")
ann <- fread("results/counts/human_gene_annotation.tsv")
samples <- fread("../metadata/pdx_samples.tsv")
coldata <- as.data.frame(samples); rownames(coldata) <- coldata$model
coldata <- coldata[colnames(txi$counts), ]
coldata$group      <- factor(coldata$group, levels = c("nonTNBC", "TNBC"))
coldata$instrument <- factor(coldata$instrument, levels = c("HiSeq2000", "HiSeq2500", "NovaSeq6000"))
coldata$subtype    <- factor(coldata$subtype, levels = names(subtype_cols))
libs <- colSums(txi$counts) / 1e6
msg("  human counts per model (M): min ", round(min(libs), 1), " (", names(which.min(libs)), "), median ",
    round(median(libs), 1), ", max ", round(max(libs), 1))
msg("  BCM-3807: ", round(libs[["BCM-3807"]], 1), " M human counts")

dds <- DESeqDataSetFromTximport(txi, coldata, design = ~ instrument + group)
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
msg("  genes after filter: ", nrow(dds))

# =============================================================================
# Step 2: PCA
# =============================================================================
vsd <- vst(dds, blind = TRUE)
top <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 2000)
pca <- prcomp(t(assay(vsd)[top, ])); pve <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pc <- data.table(coldata, pca$x[, 1:4])
save_plot(ggplot(pc, aes(PC1, PC2, colour = subtype, shape = instrument)) + geom_point(size = 2.5) +
              scale_colour_manual(values = subtype_cols) +
              labs(x = paste0("PC1 (", pve[1], "%)"), y = paste0("PC2 (", pve[2], "%)"),
                   title = "Track B (Salmon, human genes): VST, top 2000 genes"),
          "pca_subtype_instrument", w = 6.5, h = 4.5)

# =============================================================================
# Step 3: DESeq2
# =============================================================================
msg("Step 3: DESeq2 ~ instrument + group")
dds <- DESeq(dds, quiet = TRUE)
coef_name <- "group_TNBC_vs_nonTNBC"
res <- results(dds, name = coef_name, alpha = 0.05)
shr <- lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE)
de <- data.table(gene_id = rownames(res), baseMean = res$baseMean, log2FC = res$log2FoldChange,
                 lfcSE = res$lfcSE, stat = res$stat, pvalue = res$pvalue, padj = res$padj,
                 log2FC_shrunk = shr$log2FoldChange)
de <- merge(ann, de, by = "gene_id")
stopifnot(!anyDuplicated(de$gene_id))          # exactly one row per gene
setorder(de, padj, pvalue, na.last = TRUE)
de[, sig := !is.na(padj) & padj < 0.05]

# Flag "sporadic" genes: detected (count > 0) in fewer than half the models of
# BOTH groups. They can get huge raw fold changes (e.g. zero in every non-TNBC
# model) from a handful of samples; apeglm shrinks them to ~0. Flagged, not removed.
nz <- counts(dds) > 0
det_tnbc <- rowMeans(nz[, coldata$group == "TNBC"]); det_non <- rowMeans(nz[, coldata$group == "nonTNBC"])
de[, `:=`(detected_frac_TNBC = round(det_tnbc[gene_id], 2), detected_frac_nonTNBC = round(det_non[gene_id], 2))]
de[, sporadic := pmax(detected_frac_TNBC, detected_frac_nonTNBC) < 0.5]

esr1 <- de[gene_name == "ESR1", log2FC]
msg("  sign check: ESR1 log2FC = ", round(esr1, 2))
if (length(esr1) != 1 || esr1 >= 0) stop("ESR1 is not lower in TNBC: check the group coding!")

f <- "results/de/DE_TNBC_vs_nonTNBC_instrument_adjusted.tsv"
writeLines(c("# Sign convention: log2FC = log2(TNBC / nonTNBC); reference level = nonTNBC; positive = higher in TNBC",
             "# Design: ~ instrument + group; DESeq2 Wald test; BH padj; log2FC_shrunk = apeglm",
             "# Counts: Salmon (human GENCODE v50 + mouse vM39 + ERCC index, genome decoys), human genes only; mouse_identical_tx = gene has a transcript identical to a mouse one; sporadic = detected in < 50% of models of both groups (raw log2FC unreliable; use log2FC_shrunk)"), f)
fwrite(de, f, sep = "\t", append = TRUE, col.names = TRUE)
msg("  significant: ", sum(de$sig), " (", de[sig & log2FC > 0, .N], " up, ", de[sig & log2FC < 0, .N], " down in TNBC)")
msg("  significant and not sporadic: ", de[sig & !sporadic, .N], " (", de[sig & !sporadic & log2FC > 0, .N], " up, ",
    de[sig & !sporadic & log2FC < 0, .N], " down); sporadic among significant: ", de[sig & sporadic, .N])
print(de[1:15, .(gene_name, log2FC = round(log2FC, 2), padj = signif(padj, 3))])

# =============================================================================
# Step 4: volcano plot
# =============================================================================
vol <- de[!is.na(padj)]
vol[, direction := fifelse(sig & log2FC > 0, "higher in TNBC", fifelse(sig & log2FC < 0, "lower in TNBC", "n.s."))]
lab <- head(vol[sig == TRUE, gene_name], 20)
save_plot(ggplot(vol, aes(log2FC_shrunk, -log10(padj), colour = direction)) +
              geom_point(size = 0.7, alpha = 0.6) +
              scale_colour_manual(values = c("higher in TNBC" = "#D55E00", "lower in TNBC" = "#0072B2", "n.s." = "grey75")) +
              geom_hline(yintercept = -log10(0.05), linetype = 2) +
              geom_text_repel(data = vol[gene_name %in% lab], aes(label = gene_name), size = 2.8, colour = "black", max.overlaps = 30) +
              labs(x = "log2 fold change, TNBC vs non-TNBC (apeglm)", y = "-log10 adjusted p", colour = NULL,
                   title = paste0("Track B (Salmon): TNBC (n = ", sum(coldata$group == "TNBC"), ") vs non-TNBC (n = ",
                                  sum(coldata$group == "nonTNBC"), ")")),
          "volcano_TNBC_vs_nonTNBC", w = 7.5, h = 6)

# =============================================================================
# Step 5: sensitivity — drop HiSeq2500; drop BCM-3807; same 52 models as Track A
# =============================================================================
msg("Step 5: sensitivity analyses")
refit <- function(keep_models, design) {
    sub <- txi                                   # same tximport object, fewer models
    for (slot in c("abundance", "counts", "length")) sub[[slot]] <- txi[[slot]][, keep_models]
    d <- DESeqDataSetFromTximport(sub, droplevels(coldata[keep_models, ]), design = design)
    d <- DESeq(d[rowSums(counts(d) >= 10) >= 3, ], quiet = TRUE)
    r <- results(d, name = coef_name, alpha = 0.05)
    data.table(gene_id = rownames(r), log2FC = r$log2FoldChange, padj = r$padj)
}
compare <- function(s, label) {
    m <- merge(de[, .(gene_id, main = log2FC, main_sig = sig)], s, by = "gene_id")
    ms <- m[main_sig == TRUE & !is.na(log2FC)]
    data.table(analysis = label, n_sig = sum(s$padj < 0.05, na.rm = TRUE),
               main_sig_reproduced = sum(ms$padj < 0.05, na.rm = TRUE),
               main_sig_same_sign_pct = round(100 * mean(sign(ms$log2FC) == sign(ms$main)), 1),
               spearman_all = round(cor(m$main, m$log2FC, method = "spearman", use = "complete"), 3))
}
all_m <- rownames(coldata)
sens <- rbindlist(list(
    compare(refit(all_m[coldata$instrument != "HiSeq2500"], ~ instrument + group), "drop HiSeq2500"),
    compare(refit(setdiff(all_m, "BCM-3807"), ~ instrument + group), "drop BCM-3807 (= Track A's 52 models)")))
fwrite(sens, "results/de/sensitivity_summary.tsv", sep = "\t")
print(sens)

# =============================================================================
# Step 6: concordance with Track A (authors' Xenome + RSEM counts)
# =============================================================================
# Track A uses RefSeq gene symbols, Track B GENCODE IDs: match on gene name
# (genes whose name appears once in each table only).
msg("Step 6: concordance with Track A")
ta <- fread("../results/trackA/de/DE_TNBC_vs_nonTNBC_instrument_adjusted.tsv", skip = "gene\t")
tb <- de[!is.na(stat) & !gene_name %in% gene_name[duplicated(gene_name)]]
m <- merge(ta[, .(gene_name = gene, A_log2FC = log2FC, A_padj = padj, A_sig = sig)],
           tb[, .(gene_name, B_log2FC = log2FC, B_padj = padj, B_sig = sig, mouse_identical_tx)], by = "gene_name")
conc <- data.table(
    genes_matched = nrow(m),
    sig_A = sum(m$A_sig), sig_B = sum(m$B_sig), sig_both = sum(m$A_sig & m$B_sig),
    pct_A_sig_also_B = round(100 * sum(m$A_sig & m$B_sig) / sum(m$A_sig), 1),
    same_sign_among_both_pct = round(100 * mean(sign(m[A_sig & B_sig, A_log2FC]) == sign(m[A_sig & B_sig, B_log2FC])), 1),
    pearson_log2FC = round(cor(m$A_log2FC, m$B_log2FC, use = "complete"), 3),
    spearman_log2FC = round(cor(m$A_log2FC, m$B_log2FC, method = "spearman", use = "complete"), 3))
fwrite(conc, "results/de/concordance_with_trackA.tsv", sep = "\t")
fwrite(m, "results/de/trackA_vs_trackB_genes.tsv", sep = "\t")
print(conc)

m[, class := fifelse(A_sig & B_sig, "both", fifelse(A_sig, "Track A only", fifelse(B_sig, "Track B only", "neither")))]
save_plot(ggplot(m[order(class != "neither")], aes(A_log2FC, B_log2FC, colour = class)) +
              geom_point(size = 0.6, alpha = 0.6) + geom_abline(linetype = 2) +
              scale_colour_manual(values = c(both = "#CC79A7", "Track A only" = "#E69F00", "Track B only" = "#56B4E9", neither = "grey80")) +
              labs(x = "Track A log2FC (authors' Xenome + RSEM)", y = "Track B log2FC (Salmon human+mouse)",
                   colour = "padj < 0.05", title = paste0("TNBC vs non-TNBC: r = ", conc$pearson_log2FC)),
          "trackA_vs_trackB_log2FC", w = 6.5, h = 5.5)

saveRDS(dds, "results/de/dds_trackB.rds")
msg("Done.")
