#!/usr/bin/env Rscript
# =============================================================================
# 12_deseq2_trackA.R — Track A, step 3: QC and differential expression,
#                      TNBC vs non-TNBC PDX models (authors' RSEM counts).
#
# SIGN CONVENTION: non-TNBC is the reference level, so
#   log2FoldChange = log2(TNBC / non-TNBC): positive = higher in TNBC.
#
# WHAT IT DOES
#   1. Loads the 53-model count matrix and sample sheet (step 11).
#   2. Normalized matrices: upper-quartile + log2 (as in the paper), and VST.
#   3. QC: PCA, how much each PC is explained by instrument vs group, and a
#      subtype marker panel (ESR1, PGR, ERBB2, KRT5, ...) to check the labels.
#   4. DESeq2 with design ~ instrument + group; apeglm-shrunk fold changes.
#   5. Volcano, MA plot, heatmap of the top 50 genes.
#   6. Sensitivity analyses + exploratory subtype contrasts.
#
# RUN:  scripts/submit.sh scripts/run_R.slurm scripts/12_deseq2_trackA.R
# =============================================================================
suppressPackageStartupMessages({
    library(data.table); library(DESeq2); library(ggplot2); library(ggrepel)
    library(pheatmap); library(patchwork); library(matrixStats)
})
msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)
set.seed(1)

out <- "results/trackA"
for (d in c("qc", "de", "plots")) dir.create(file.path(out, d), recursive = TRUE, showWarnings = FALSE)

# Colours used in every plot of this project
group_cols   <- c(TNBC = "#D55E00", nonTNBC = "#0072B2")
subtype_cols <- c(TNBC = "#D55E00", "ER+" = "#0072B2", "HER2+" = "#009E73", "ER+/HER2+" = "#CC79A7")
theme_set(theme_classic(base_size = 11) + theme(legend.position = "right"))

# save_plot: every figure as PDF and PNG
save_plot <- function(p, name, w = 6, h = 5) {
    ggsave(file.path(out, "plots", paste0(name, ".pdf")), p, width = w, height = h)
    ggsave(file.path(out, "plots", paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}

# =============================================================================
# Step 1: load data
# =============================================================================
msg("Step 1: loading counts and sample sheet")
samples <- fread("metadata/pdx_samples.tsv")
counts  <- fread("results/counts/authors_rsem_expected_counts.tsv.gz")
cm <- as.matrix(counts[, -1]); rownames(cm) <- counts$gene
stopifnot(identical(colnames(cm), samples$model))

# RSEM expected counts are fractional; DESeq2 needs whole numbers.
cm <- round(cm)
storage.mode(cm) <- "integer"

coldata <- as.data.frame(samples)
rownames(coldata) <- coldata$model
coldata$group      <- factor(coldata$group, levels = c("nonTNBC", "TNBC"))   # reference first
coldata$instrument <- factor(coldata$instrument, levels = c("HiSeq2000", "HiSeq2500", "NovaSeq6000"))
coldata$subtype    <- factor(coldata$subtype, levels = names(subtype_cols))

# BCM-3807 has only ~1.4 M counts in the authors' RSEM output (median ~42 M):
# their RSEM apparently used only its 7.7 M paired spots (the run is mostly
# single-end). It is left out of the main model and tested in sensitivity 4.
exclude_main <- "BCM-3807"
cm_all <- cm; coldata_all <- coldata
cm <- cm[, !colnames(cm) %in% exclude_main]
coldata <- droplevels(coldata[colnames(cm), ])
msg("  main analysis: ", sum(coldata$group == "TNBC"), " TNBC vs ", sum(coldata$group == "nonTNBC"),
    " non-TNBC (excluded: ", paste(exclude_main, collapse = ", "), ")")

# Pre-filter (lab default): at least 10 counts in at least 3 models
keep <- rowSums(cm >= 10) >= 3
msg("  genes: ", nrow(cm), " -> ", sum(keep), " after keeping >=10 counts in >=3 models")
cm <- cm[keep, ]

# =============================================================================
# Step 2: normalized expression matrices
# =============================================================================
# Upper-quartile (UQ): divide each model by the 75th percentile of its non-zero
# gene counts, scale to 1000 (the TCGA RSEM convention), then log2(x + 1).
msg("Step 2: upper-quartile log2 and VST matrices")
uq <- apply(cm, 2, function(x) quantile(x[x > 0], 0.75))
log2uq <- log2(sweep(cm, 2, uq, "/") * 1000 + 1)
fwrite(data.table(gene = rownames(log2uq), round(log2uq, 4)),
       file.path(out, "qc", "log2_UQ_expression.tsv.gz"), sep = "\t")

dds <- DESeqDataSetFromMatrix(cm, coldata, design = ~ instrument + group)
vsd <- vst(dds, blind = TRUE)
fwrite(data.table(gene = rownames(vsd), round(assay(vsd), 4)),
       file.path(out, "qc", "vst_expression.tsv.gz"), sep = "\t")

# =============================================================================
# Step 3: QC — PCA and subtype markers
# =============================================================================
msg("Step 3: PCA")
top <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 2000)   # 2000 most variable genes
pca <- prcomp(t(assay(vsd)[top, ]))
pve <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pc  <- data.table(coldata, pca$x[, 1:10])

# How much of each PC does instrument or group explain? (R-squared of a linear model)
r2 <- rbindlist(lapply(1:10, function(i) data.table(
    PC = paste0("PC", i), pct_variance = pve[i],
    R2_instrument = summary(lm(pc[[paste0("PC", i)]] ~ pc$instrument))$r.squared,
    R2_group      = summary(lm(pc[[paste0("PC", i)]] ~ pc$group))$r.squared)))
fwrite(r2[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)],
       file.path(out, "qc", "pca_variance_explained.tsv"), sep = "\t")
print(r2[1:5])

p1 <- ggplot(pc, aes(PC1, PC2, colour = subtype, shape = instrument)) +
    geom_point(size = 2.5) + scale_colour_manual(values = subtype_cols) +
    labs(x = paste0("PC1 (", pve[1], "%)"), y = paste0("PC2 (", pve[2], "%)"),
         title = "VST, 2000 most variable genes")
p2 <- ggplot(pc, aes(PC3, PC4, colour = subtype, shape = instrument)) +
    geom_point(size = 2.5) + scale_colour_manual(values = subtype_cols) +
    labs(x = paste0("PC3 (", pve[3], "%)"), y = paste0("PC4 (", pve[4], "%)"))
save_plot(p1 + p2 + plot_layout(guides = "collect"), "pca_subtype_instrument", w = 11, h = 4.5)

# Marker panel: do the labels agree with the expression of ER / HER2 / basal genes?
msg("Step 3b: subtype marker panel")
markers <- c("ESR1", "PGR", "GATA3", "FOXA1", "ERBB2", "GRB7", "KRT5", "KRT14", "FOXC1")
markers <- intersect(markers, rownames(log2uq))
mk <- melt(data.table(coldata[, c("model", "subtype", "instrument")], t(log2uq[markers, ])),
           id.vars = c("model", "subtype", "instrument"), variable.name = "gene", value.name = "log2UQ")
pm <- ggplot(mk, aes(subtype, log2UQ, colour = subtype)) +
    geom_boxplot(outlier.shape = NA, colour = "grey40") +
    geom_jitter(width = 0.2, size = 1.4) +
    facet_wrap(~ gene, scales = "free_y", nrow = 3) +
    scale_colour_manual(values = subtype_cols) +
    labs(x = NULL, y = "log2 UQ expression") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
save_plot(pm, "subtype_markers", w = 9, h = 8)

# Flag models whose markers contradict the label:
#   ER+ (incl. ER+/HER2+) with ESR1 at or below the TNBC 90th percentile,
#   HER2+ (incl. ER+/HER2+) with ERBB2 at or below the TNBC 90th percentile,
#   TNBC with ESR1 or ERBB2 above the median of the positive groups.
wide <- dcast(mk, model + subtype + instrument ~ gene, value.var = "log2UQ")
tn <- wide[subtype == "TNBC"]
er_pos  <- wide$subtype %in% c("ER+", "ER+/HER2+")
her_pos <- wide$subtype %in% c("HER2+", "ER+/HER2+")
wide[, flag := ""]
wide[er_pos  & ESR1  <= quantile(tn$ESR1, 0.9),  flag := paste(flag, "ER+ label but ESR1 low;")]
wide[her_pos & ERBB2 <= quantile(tn$ERBB2, 0.9), flag := paste(flag, "HER2+ label but ERBB2 low;")]
wide[subtype == "TNBC" & ESR1  > median(wide[er_pos,  ESR1]),  flag := paste(flag, "TNBC but ESR1 high;")]
wide[subtype == "TNBC" & ERBB2 > median(wide[her_pos, ERBB2]), flag := paste(flag, "TNBC but ERBB2 high;")]
wide[, flag := trimws(flag)]
fwrite(wide[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)],
       file.path(out, "qc", "subtype_marker_check.tsv"), sep = "\t")
msg("  models with a marker/label conflict: ", sum(wide$flag != ""))
if (any(wide$flag != "")) print(wide[flag != "", .(model, subtype, ESR1 = round(ESR1, 1), ERBB2 = round(ERBB2, 1), flag)])

# =============================================================================
# Step 4: DESeq2, TNBC vs non-TNBC, adjusting for instrument
# =============================================================================
# Note: all 9 HiSeq2500 models are TNBC. With instrument in the model they
# help estimate gene variability but carry no information about TNBC vs
# non-TNBC; that comparison comes from HiSeq2000 (22 v 7) and NovaSeq (8 v 7).
msg("Step 4: DESeq2  ~ instrument + group")
dds <- DESeq(dds, quiet = TRUE)
coef_name <- "group_TNBC_vs_nonTNBC"
stopifnot(coef_name %in% resultsNames(dds))
res <- results(dds, name = coef_name, alpha = 0.05)
shr <- lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE)

# Hard check of the sign convention: ESR1 must be LOWER in TNBC (log2FC < 0).
esr1_lfc <- res["ESR1", "log2FoldChange"]
msg("  sign check: ESR1 log2FC (TNBC - nonTNBC) = ", round(esr1_lfc, 2))
if (!is.finite(esr1_lfc) || esr1_lfc >= 0) stop("ESR1 is not lower in TNBC: check the group coding!")

grp_mean <- function(g) rowMeans(log2uq[, coldata$group == g, drop = FALSE])
de <- data.table(gene = rownames(res), baseMean = res$baseMean,
                 log2FC = res$log2FoldChange, lfcSE = res$lfcSE, stat = res$stat,
                 pvalue = res$pvalue, padj = res$padj, log2FC_shrunk = shr$log2FoldChange,
                 mean_log2UQ_TNBC = grp_mean("TNBC"), mean_log2UQ_nonTNBC = grp_mean("nonTNBC"))
setorder(de, padj, pvalue, na.last = TRUE)   # genes without padj go last
de[, sig := !is.na(padj) & padj < 0.05]
de[, direction := fifelse(sig & log2FC > 0, "up_in_TNBC", fifelse(sig & log2FC < 0, "down_in_TNBC", "ns"))]

write_de <- function(dt, file, design_txt) {
    f <- file.path(out, "de", file)
    writeLines(c("# Sign convention: log2FC = log2(TNBC / nonTNBC); reference level = nonTNBC; positive = higher in TNBC",
                 paste0("# Design: ", design_txt, "; DESeq2 Wald test; BH-adjusted padj; log2FC_shrunk = apeglm"),
                 "# Counts: authors' RSEM expected counts (Xenome-filtered human reads, GEO GSE183187), rounded"), f)
    fwrite(dt, f, sep = "\t", append = TRUE, col.names = TRUE)
}
write_de(de, "DE_TNBC_vs_nonTNBC_instrument_adjusted.tsv", "~ instrument + group")
saveRDS(dds, file.path(out, "de", "dds_instrument_adjusted.rds"))

n_up <- sum(de$direction == "up_in_TNBC"); n_dn <- sum(de$direction == "down_in_TNBC")
msg("  significant (padj < 0.05): ", n_up + n_dn, "  (", n_up, " up, ", n_dn, " down in TNBC); ",
    "with |log2FC| > 1: ", de[sig & abs(log2FC) > 1, .N])
print(de[1:15, .(gene, log2FC = round(log2FC, 2), padj = signif(padj, 3))])

# =============================================================================
# Step 5: plots — volcano, MA, heatmap of the top 50 genes
# =============================================================================
msg("Step 5: plots")
vol <- copy(de[!is.na(padj)])
vol[, nlp := -log10(padj)]
lab <- unique(c(vol[sig == TRUE][order(padj)][1:15, gene],
                vol[sig == TRUE][order(-abs(log2FC_shrunk))][1:15, gene]))
pv <- ggplot(vol, aes(log2FC_shrunk, nlp)) +
    geom_point(aes(colour = direction), size = 0.8, alpha = 0.7) +
    scale_colour_manual(values = c(up_in_TNBC = group_cols[["TNBC"]], down_in_TNBC = group_cols[["nonTNBC"]], ns = "grey75"),
                        labels = c(up_in_TNBC = paste0("higher in TNBC (", n_up, ")"),
                                   down_in_TNBC = paste0("lower in TNBC (", n_dn, ")"), ns = "n.s.")) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey50") +
    geom_vline(xintercept = c(-1, 1), linetype = 2, colour = "grey50") +
    geom_text_repel(data = vol[gene %in% lab], aes(label = gene), size = 2.8, max.overlaps = 40) +
    labs(x = "log2 fold change, TNBC vs non-TNBC (apeglm)", y = "-log10 adjusted p", colour = NULL,
         title = paste0("TNBC (n = ", sum(coldata$group == "TNBC"), ") vs non-TNBC (n = ",
                        sum(coldata$group == "nonTNBC"), ") PDX"), subtitle = "DESeq2 ~ instrument + group")
save_plot(pv, "volcano_TNBC_vs_nonTNBC", w = 7.5, h = 6)

pma <- ggplot(vol, aes(log10(baseMean), log2FC_shrunk, colour = direction)) +
    geom_point(size = 0.6, alpha = 0.7) +
    scale_colour_manual(values = c(up_in_TNBC = group_cols[["TNBC"]], down_in_TNBC = group_cols[["nonTNBC"]], ns = "grey75")) +
    geom_hline(yintercept = 0) + labs(x = "log10 mean normalized count", y = "log2FC (apeglm)", colour = NULL)
save_plot(pma, "MA_TNBC_vs_nonTNBC", w = 6.5, h = 4.5)

# Heatmap: VST corrected for instrument (for display only), rows scaled
top50 <- head(de[sig == TRUE, gene], 50)
if (length(top50) >= 2) {
    vsd_adj <- assay(vsd)
    vsd_adj <- t(apply(vsd_adj, 1, function(y) residuals(lm(y ~ coldata$instrument)) + mean(y)))
    ann <- data.frame(subtype = coldata$subtype, instrument = coldata$instrument, row.names = coldata$model)
    ord <- order(coldata$group, coldata$subtype)
    hm <- pheatmap(vsd_adj[top50, ord], scale = "row", cluster_cols = FALSE, show_colnames = FALSE,
                   annotation_col = ann, annotation_colors = list(subtype = subtype_cols),
                   color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
                   breaks = seq(-2.5, 2.5, length.out = 102), fontsize_row = 7,
                   main = "Top 50 DE genes (VST, instrument-corrected, row z-score)", silent = TRUE)
    for (ext in c("pdf", "png")) {
        if (ext == "pdf") pdf(file.path(out, "plots", "heatmap_top50.pdf"), 8, 9)
        else png(file.path(out, "plots", "heatmap_top50.png"), 8, 9, units = "in", res = 300)
        grid::grid.draw(hm$gtable); dev.off()
    }
}

# =============================================================================
# Step 6: sensitivity analyses — do the main results hold up?
# =============================================================================
# run_de: refit DESeq2 on a subset of models with a given design; returns the
# Wald statistic, log2FC and padj for the TNBC vs non-TNBC (or other) coefficient.
run_de <- function(keep_models, design, coef = coef_name, ref = "nonTNBC", grp_col = "group",
                   mat = cm, cdata = coldata) {
    cd <- droplevels(cdata[keep_models, ])
    cd$grp <- relevel(factor(cd[[grp_col]]), ref = ref)
    d <- DESeqDataSetFromMatrix(mat[, keep_models], cd, design = design)
    d <- d[rowSums(counts(d) >= 10) >= 3, ]
    d <- DESeq(d, quiet = TRUE)
    r <- results(d, name = coef, alpha = 0.05)
    data.table(gene = rownames(r), log2FC = r$log2FoldChange, stat = r$stat, padj = r$padj)
}
compare_to_main <- function(s, label, n_tnbc, n_other) {
    m <- merge(de[, .(gene, main_lfc = log2FC, main_sig = sig)], s, by = "gene")
    ms <- m[main_sig == TRUE & !is.na(log2FC)]
    data.table(analysis = label, n_TNBC = n_tnbc, n_comparator = n_other,
               n_sig = sum(s$padj < 0.05, na.rm = TRUE),
               overlap_with_main_sig = sum(ms$padj < 0.05, na.rm = TRUE),
               main_sig_same_sign_pct = round(100 * mean(sign(ms$log2FC) == sign(ms$main_lfc)), 1),
               spearman_log2FC_all_genes = round(cor(m$main_lfc, m$log2FC, method = "spearman", use = "complete"), 3))
}

msg("Step 6: sensitivity analyses")
all_m  <- coldata$model
grp_n  <- function(ms) c(sum(coldata_all[ms, "group"] == "TNBC"), sum(coldata_all[ms, "group"] == "nonTNBC"))
sens <- list()

ms <- all_m[coldata$instrument != "HiSeq2500"]; n <- grp_n(ms)
s1 <- run_de(ms, ~ instrument + grp, coef = "grp_TNBC_vs_nonTNBC")
sens[[1]] <- compare_to_main(s1, "1. drop HiSeq2500 (~ instrument + group)", n[1], n[2])

n <- grp_n(all_m)
s2 <- run_de(all_m, ~ grp, coef = "grp_TNBC_vs_nonTNBC")
sens[[2]] <- compare_to_main(s2, "2. no instrument covariate (~ group)", n[1], n[2])

for (inst in c("HiSeq2000", "NovaSeq6000")) {
    ms <- all_m[coldata$instrument == inst]; n <- grp_n(ms)
    s3 <- run_de(ms, ~ grp, coef = "grp_TNBC_vs_nonTNBC")
    sens[[length(sens) + 1]] <- compare_to_main(s3, paste0("3. ", inst, " only (~ group)"), n[1], n[2])
}

ms <- coldata_all$model; n <- grp_n(ms)
s4 <- run_de(ms, ~ instrument + grp, coef = "grp_TNBC_vs_nonTNBC",
             mat = cm_all[rownames(cm), ], cdata = coldata_all)
sens[[length(sens) + 1]] <- compare_to_main(s4, "4. include BCM-3807 (1.4 M counts)", n[1], n[2])

sens <- rbindlist(sens)
fwrite(sens, file.path(out, "de", "sensitivity_summary.tsv"), sep = "\t")
print(sens)

# =============================================================================
# Step 7: exploratory subtype contrasts (small comparator groups)
# =============================================================================
# The 2 ER+/HER2+ models are left out of both, so each comparator is one subtype.
msg("Step 7: exploratory contrasts TNBC vs ER+ and TNBC vs HER2+")
for (cmp in c("ER+", "HER2+")) {
    ms <- all_m[coldata$subtype %in% c("TNBC", cmp)]
    cd_inst <- droplevels(coldata[ms, "instrument"])
    design <- if (nlevels(cd_inst) > 1) ~ instrument + grp else ~ grp
    tag <- gsub("[^A-Za-z0-9]", "", cmp)                          # "ER", "HER2"
    # "ER+" is not a valid R name, so code the comparator as "other"
    coldata$cmp_grp <- ifelse(coldata$subtype == "TNBC", "TNBC", "other")
    s <- run_de(ms, design, coef = "grp_TNBC_vs_other", ref = "other", grp_col = "cmp_grp")
    setorder(s, padj, na.last = TRUE)
    write_de(s, paste0("DE_exploratory_TNBC_vs_", tag, ".tsv"),
             paste0(deparse(design), "; TNBC (n = ", sum(coldata[ms, "subtype"] == "TNBC"),
                    ") vs ", cmp, " (n = ", sum(coldata[ms, "subtype"] == cmp), "); EXPLORATORY; log2FC = TNBC - ", cmp))
    msg("  TNBC vs ", cmp, ": ", sum(s$padj < 0.05, na.rm = TRUE), " genes padj < 0.05")
}

writeLines(capture.output(sessionInfo()), file.path(out, "de", "sessionInfo_12.txt"))
msg("Done. Results in ", out, "/{qc,de,plots}")
