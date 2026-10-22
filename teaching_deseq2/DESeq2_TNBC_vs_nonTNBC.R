# =============================================================================
#  Differential gene expression with DESeq2
#  TNBC vs non-TNBC breast cancer PDX models (GEO GSE183187)
# =============================================================================
#
#  QUESTION: which genes are expressed differently in triple-negative breast
#  cancer (TNBC) compared with other breast cancers (ER+, HER2+)?
#
#  DATA (in the data/ folder):
#    counts_TNBC_PDX.csv.gz  how many RNA-seq reads came from each gene (rows)
#                            in each tumour model (columns). Raw counts, NOT
#                            normalized: DESeq2 does the normalization itself.
#    sample_info.csv         one row per tumour model: its group (TNBC or
#                            nonTNBC), subtype and the sequencing machine used.
#
#  The models are patient-derived xenografts (PDX): human tumours grown in mice.
#  The authors already removed the mouse reads, so these counts are human only.
#
#  HOW TO USE THIS SCRIPT
#    Open the RStudio project (TNBC_DESeq2.Rproj) so file paths work, then run
#    the script ONE SECTION AT A TIME (select lines + Ctrl/Cmd + Enter). Read
#    the comments and look at each result before moving on.
# =============================================================================


# =============================================================================
# 1. Packages
# =============================================================================
# Run the install lines ONCE (remove the "#" first). Bioconductor packages
# (DESeq2, apeglm) are installed with BiocManager, not install.packages().
#
# install.packages(c("BiocManager", "ggplot2", "ggrepel", "pheatmap"))
# BiocManager::install(c("DESeq2", "apeglm"))

library(DESeq2)     # the statistics
library(ggplot2)    # plots
library(ggrepel)    # gene labels on plots that don't overlap
library(pheatmap)   # heatmaps

dir.create("results", showWarnings = FALSE)   # all output goes here


# =============================================================================
# 2. Read the data
# =============================================================================
counts <- read.csv("data/counts_TNBC_PDX.csv.gz", row.names = 1, check.names = FALSE)
info   <- read.csv("data/sample_info.csv", row.names = 1)

dim(counts)          # genes x samples
counts[1:5, 1:4]     # a small corner of the table
head(info)
table(info$group)                    # how many models per group?
table(info$subtype, info$group)      # which subtypes are "nonTNBC"?
table(info$instrument, info$group)   # which machine sequenced which group?

# QUESTION: Look at the last table. Is every machine used for both groups?


# =============================================================================
# 3. Check that samples match  (the most common source of silent mistakes!)
# =============================================================================
# The columns of the count table must be the same samples, in the same order,
# as the rows of the sample table. If not, every result would be wrong.
all(colnames(counts) == rownames(info))     # must print TRUE

# If it printed FALSE, this line puts the count columns in the right order:
counts <- counts[, rownames(info)]
stopifnot(all(colnames(counts) == rownames(info)))


# =============================================================================
# 4. Tell R which group is the reference
# =============================================================================
# Fold changes are always "group B compared with reference group A". We make
# nonTNBC the reference, so:
#     log2FoldChange > 0  ->  gene is HIGHER in TNBC
#     log2FoldChange < 0  ->  gene is LOWER in TNBC
info$group      <- factor(info$group, levels = c("nonTNBC", "TNBC"))
info$instrument <- factor(info$instrument)
levels(info$group)   # the first level is the reference


# =============================================================================
# 5. Build the DESeq2 object
# =============================================================================
# The "design" tells DESeq2 what to compare:
#   ~ instrument + group
# means: compare TNBC with nonTNBC (group) while correcting for differences
# between sequencing machines (instrument). The variable you want to test goes
# LAST. Correcting for the machine matters because samples run on different
# machines/years can differ for technical reasons (a "batch effect").
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = info,
                              design    = ~ instrument + group)
dds

# Remove genes with very few reads: keep genes with at least 10 reads in at
# least 3 samples. Fewer genes = faster, and fewer tests to correct for.
keep <- rowSums(counts(dds) >= 10) >= 3
table(keep)
dds <- dds[keep, ]


# =============================================================================
# 6. Quality check: PCA (do samples group by biology or by machine?)
# =============================================================================
# vst() puts counts on a log-like scale so that samples can be compared.
# PCA then squeezes thousands of genes into 2 dimensions (PC1, PC2).
vsd <- vst(dds, blind = TRUE)
pca <- plotPCA(vsd, intgroup = c("subtype", "instrument"), returnData = TRUE)
pct <- round(100 * attr(pca, "percentVar"))

ggplot(pca, aes(PC1, PC2, colour = subtype, shape = instrument)) +
    geom_point(size = 3) +
    labs(x = paste0("PC1: ", pct[1], "% of variance"),
         y = paste0("PC2: ", pct[2], "% of variance"),
         title = "PCA of 52 breast cancer PDX models") +
    theme_classic()
ggsave("results/PCA.png", width = 6, height = 4.5)

# QUESTION: Do TNBC models separate from the others along PC1 or PC2?
#           Do you see clusters by machine (shape)?


# =============================================================================
# 7. Run DESeq2
# =============================================================================
# One command does all the work:
#   1. size factors  : corrects for different sequencing depth per sample
#   2. dispersions   : how variable each gene is between samples of a group
#   3. statistics    : a test for every gene (negative binomial model)
# It takes a minute or two.
dds <- DESeq(dds)
resultsNames(dds)    # the comparisons DESeq2 can report


# =============================================================================
# 8. Get the results: TNBC vs nonTNBC
# =============================================================================
res <- results(dds, contrast = c("group", "TNBC", "nonTNBC"), alpha = 0.05)
summary(res)

# Columns you will use:
#   baseMean        average normalized count (how highly expressed the gene is)
#   log2FoldChange  1 = 2x higher in TNBC, -1 = 2x lower, 2 = 4x higher ...
#   pvalue          probability of a difference this big by chance alone
#   padj            p-value corrected for testing ~20,000 genes (use THIS one)
#
# We call a gene "significant" if padj < 0.05.
res <- res[order(res$padj), ]     # most significant first
head(res, 20)
sum(res$padj < 0.05, na.rm = TRUE)                                # how many?
sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)       # higher in TNBC
sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)       # lower in TNBC

# SANITY CHECK: TNBC tumours lack the estrogen receptor, so ESR1 must be LOWER
# in TNBC (negative log2FoldChange). ERBB2 (= HER2) should be lower too.
res[c("ESR1", "PGR", "ERBB2"), c("baseMean", "log2FoldChange", "padj")]

# NOTE: 9 of the 38 TNBC models come from one machine (HiSeq2500) that ran no
# nonTNBC models. With instrument in the design, those 9 samples cannot
# tell DESeq2 anything about TNBC vs nonTNBC - the comparison comes from the
# other two machines. Good to know when you interpret the results.


# =============================================================================
# 9. Shrink fold changes (better for plots and ranking)
# =============================================================================
# Genes with few reads can show huge but unreliable fold changes. "Shrinkage"
# pulls those noisy estimates towards 0; well-measured genes barely change.
res_shrunk <- lfcShrink(dds, coef = "group_TNBC_vs_nonTNBC", type = "apeglm")
res_shrunk <- res_shrunk[rownames(res), ]       # same order as res


# =============================================================================
# 10. Plots
# =============================================================================
# --- MA plot: fold change vs expression level -------------------------------
plotMA(res_shrunk, ylim = c(-6, 6), main = "TNBC vs nonTNBC (blue = padj < 0.05)")

# --- Volcano plot: fold change vs significance -------------------------------
vol <- data.frame(gene  = rownames(res),
                  lfc   = res_shrunk$log2FoldChange,
                  padj  = res$padj)
vol <- vol[!is.na(vol$padj), ]
vol$direction <- "not significant"
vol$direction[vol$padj < 0.05 & vol$lfc > 0] <- "higher in TNBC"
vol$direction[vol$padj < 0.05 & vol$lfc < 0] <- "lower in TNBC"
top_genes <- head(vol$gene, 20)                 # vol is sorted by padj

ggplot(vol, aes(lfc, -log10(padj), colour = direction)) +
    geom_point(size = 0.8, alpha = 0.6) +
    scale_colour_manual(values = c("higher in TNBC" = "#D55E00",
                                   "lower in TNBC" = "#0072B2",
                                   "not significant" = "grey75")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_text_repel(data = subset(vol, gene %in% top_genes), aes(label = gene),
                    size = 3, colour = "black", max.overlaps = 30) +
    labs(x = "log2 fold change (TNBC vs nonTNBC)", y = "-log10 adjusted p-value",
         colour = NULL, title = "Volcano plot") +
    theme_classic()
ggsave("results/volcano.png", width = 7, height = 5.5)

# --- One gene at a time: normalized counts per sample ------------------------
plot_gene <- function(g) {
    d <- plotCounts(dds, gene = g, intgroup = c("group", "subtype"), returnData = TRUE)
    ggplot(d, aes(group, count, colour = subtype)) +
        geom_jitter(width = 0.15, size = 2) +
        scale_y_log10() +
        labs(title = g, x = NULL, y = "normalized count (log scale)") +
        theme_classic()
}
plot_gene("ESR1")
plot_gene("ERBB2")
plot_gene(rownames(res)[1])       # the most significant gene

# --- Heatmap of the 30 most significant genes --------------------------------
top30 <- head(rownames(res), 30)
annotation <- info[, c("group", "subtype")]
pheatmap(assay(vsd)[top30, order(info$group)],
         scale = "row",                 # z-score per gene: colours show relative expression
         cluster_cols = FALSE,
         show_colnames = FALSE,
         annotation_col = annotation,
         main = "Top 30 genes, TNBC vs nonTNBC",
         filename = "results/heatmap_top30.png", width = 7, height = 7)


# =============================================================================
# 11. Save the results table
# =============================================================================
out <- data.frame(gene = rownames(res), as.data.frame(res),
                  log2FoldChange_shrunk = res_shrunk$log2FoldChange)
write.csv(out, "results/DESeq2_TNBC_vs_nonTNBC.csv", row.names = FALSE)
# Open this file in Excel to browse it. Careful: Excel turns some gene names
# into dates (e.g. SEPTIN5 is safe, but old names like "SEPT5" become 5-Sep)!


# =============================================================================
# 12. Exercises
# =============================================================================
# 1. How many genes are significant (padj < 0.05) AND change more than 2-fold
#    (|log2FoldChange| > 1)?
# 2. Look up the top 5 genes higher in TNBC and the top 5 lower in TNBC
#    (e.g. on genecards.org). Do they make biological sense?
# 3. Re-run sections 5-8 with design = ~ group (no machine correction).
#    Do you get more or fewer significant genes? Why might that be?
# 4. Compare TNBC with ER+ models only: keep only those samples,
#        keep_s <- info$subtype %in% c("TNBC", "ER+")
#    and use counts[, keep_s] and info[keep_s, ]. Which genes now come on top?
# 5. Many of the top genes lower in TNBC sit next to ERBB2 on chromosome 17
#    (STARD3, PGAP3, MIEN1, GRB7). Why would that be? (Hint: HER2 amplification.)

sessionInfo()   # records package versions: always include in your report
