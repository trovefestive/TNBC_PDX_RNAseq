#!/usr/bin/env Rscript
# =============================================================================
# 13_gsea_trackA.R — Track A, step 4: gene set enrichment (GSEA), TNBC vs non-TNBC.
#
# SIGN CONVENTION: genes are ranked by the DESeq2 Wald statistic of
#   TNBC vs non-TNBC (reference = non-TNBC), so
#   NES > 0 = gene set higher in TNBC, NES < 0 = lower in TNBC.
#
# Collections (MSigDB 2026.2, human symbols; downloaded by step 10):
#   GO biological process, Reactome, WikiPathways (the paper's collections)
#   Hallmark (positive control: estrogen response must be LOWER in TNBC)
#
# RUN:  scripts/submit.sh scripts/run_R.slurm scripts/13_gsea_trackA.R
# =============================================================================
suppressPackageStartupMessages({
    library(data.table); library(fgsea); library(ggplot2); library(patchwork)
})
msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)
set.seed(1)

out <- "results/trackA"
dir.create(file.path(out, "gsea"), recursive = TRUE, showWarnings = FALSE)
group_cols <- c(TNBC = "#D55E00", nonTNBC = "#0072B2")
theme_set(theme_classic(base_size = 10))
save_plot <- function(p, name, w = 7, h = 5) {
    ggsave(file.path(out, "plots", paste0(name, ".pdf")), p, width = w, height = h)
    ggsave(file.path(out, "plots", paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}

# =============================================================================
# Step 1: ranked gene list from DESeq2
# =============================================================================
msg("Step 1: ranking genes by the DESeq2 Wald statistic")
de <- fread(file.path(out, "de", "DE_TNBC_vs_nonTNBC_instrument_adjusted.tsv"), skip = "gene\t")
de <- de[!is.na(stat)]
ranks <- setNames(de$stat, de$gene)
ranks <- ranks[!duplicated(names(ranks))]
msg("  ", length(ranks), " ranked genes")

# =============================================================================
# Step 2: fgsea for each collection
# =============================================================================
collections <- c(GO_BP        = "c5.go.bp",
                 Reactome     = "c2.cp.reactome",
                 WikiPathways = "c2.cp.wikipathways",
                 Hallmark     = "h.all")
all_res <- list()
for (nm in names(collections)) {
    gmt <- Sys.glob(file.path("data/genesets", paste0(collections[[nm]], ".v*.symbols.gmt")))
    stopifnot(length(gmt) == 1)
    sets <- gmtPathways(gmt)
    msg("Step 2: ", nm, " (", length(sets), " sets)")
    r <- fgseaMultilevel(sets, ranks, minSize = 15, maxSize = 500, eps = 0)
    r[, collection := nm]
    r[, leadingEdge := vapply(leadingEdge, function(g) paste(head(g, 25), collapse = ","), "")]
    r <- r[order(padj, -abs(NES), na.last = TRUE)]
    f <- file.path(out, "gsea", paste0("GSEA_", nm, "_TNBC_vs_nonTNBC.tsv"))
    writeLines("# NES > 0 = higher in TNBC; NES < 0 = lower in TNBC (ranked by DESeq2 Wald stat, ~ instrument + group); leadingEdge = first 25 genes", f)
    fwrite(r, f, sep = "\t", append = TRUE, col.names = TRUE)
    msg("  padj < 0.05: ", r[padj < 0.05, .N], " (", r[padj < 0.05 & NES > 0, .N], " up, ",
        r[padj < 0.05 & NES < 0, .N], " down in TNBC)")
    all_res[[nm]] <- r
}
res <- rbindlist(all_res)

# =============================================================================
# Step 3: positive control — estrogen response must be lower in TNBC
# =============================================================================
pc <- res[pathway %in% c("HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_ESTROGEN_RESPONSE_LATE"),
          .(pathway, NES = round(NES, 2), padj = signif(padj, 3))]
print(pc)
if (any(pc$NES >= 0)) warning("Estrogen response is not lower in TNBC: check the ranking sign!")

# =============================================================================
# Step 4: plots — top 10 up and top 10 down per collection, enrichment curves
# =============================================================================
msg("Step 4: plots")
bar_plot <- function(r, nm) {
    s <- r[padj < 0.05]
    top <- rbind(head(s[NES > 0][order(-NES)], 10), head(s[NES < 0][order(NES)], 10))
    if (nrow(top) == 0) return(NULL)
    top[, label := substr(gsub("_", " ", sub("^[A-Z]+_", "", pathway)), 1, 60)]
    top[, label := factor(label, levels = label[order(NES)])]
    ggplot(top, aes(NES, label, fill = NES > 0)) + geom_col() +
        scale_fill_manual(values = c(`TRUE` = group_cols[["TNBC"]], `FALSE` = group_cols[["nonTNBC"]]),
                          labels = c(`TRUE` = "higher in TNBC", `FALSE` = "lower in TNBC"), name = NULL) +
        labs(x = "NES", y = NULL, title = paste0(nm, " (padj < 0.05; top 10 each way)"))
}
for (nm in names(all_res)) {
    p <- bar_plot(all_res[[nm]], nm)
    if (!is.null(p)) save_plot(p, paste0("gsea_", nm), w = 8.5, h = 6)
}

sets_h <- gmtPathways(Sys.glob("data/genesets/h.all.v*.symbols.gmt"))
curves <- lapply(c("HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_E2F_TARGETS"), function(pw)
    plotEnrichment(sets_h[[pw]], ranks) + labs(title = pw, x = "rank (higher in TNBC  <-  ->  lower in TNBC)"))
save_plot(wrap_plots(curves, ncol = 2), "gsea_enrichment_curves_hallmark", w = 10, h = 3.8)

# One summary table: significant sets from all collections
fwrite(res[padj < 0.05][order(collection, padj)], file.path(out, "gsea", "GSEA_all_collections_padj0.05.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "gsea", "sessionInfo_13.txt"))
msg("Done.")
