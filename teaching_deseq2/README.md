# Hands-on: differential gene expression with DESeq2 in RStudio

**Question:** which genes differ between **triple-negative breast cancer (TNBC)** and other breast cancers (ER+,
HER2+)?

**Data:** RNA-seq from 52 patient-derived xenograft (PDX) models (GEO
[GSE183187](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183187)): 38 TNBC and 14 non-TNBC. Each model is a
patient's tumour grown in a mouse. The authors removed the mouse reads and counted human reads per gene with RSEM.

## What's in this folder

```
teaching_deseq2/
├── TNBC_DESeq2.Rproj            double-click (or File > Open Project) to open in RStudio
├── DESeq2_TNBC_vs_nonTNBC.R     the analysis, step by step, with explanations
├── data/
│   ├── counts_TNBC_PDX.csv.gz   raw read counts: 25,143 genes (rows) x 52 samples (columns)
│   └── sample_info.csv          one row per sample: group, subtype, sequencing machine
└── results/                     created by the script (plots + results table)
```

## Before you start

1. **Open the project.** In RStudio choose *File → Open Project…* and pick `TNBC_DESeq2.Rproj`. Working this way means
   paths such as `data/sample_info.csv` resolve without `setwd()`.
2. **Install the packages once** (about 5–10 minutes):

   ```r
   install.packages(c("BiocManager", "ggplot2", "ggrepel", "pheatmap"))
   BiocManager::install(c("DESeq2", "apeglm"))
   ```

   If RStudio asks *"Update all/some/none?"*, answer `n` (none) to save time.
3. **Open `DESeq2_TNBC_vs_nonTNBC.R` and run it section by section.** Select lines and press Ctrl+Enter (Mac:
   Cmd+Enter). Read each comment and look at each output before moving on.

## The workflow

| Section | Step | What you learn |
|---|---|---|
| 2–3 | Read counts and sample sheet; check that they match | The two input tables of every RNA-seq analysis; why matching sample order matters |
| 4 | Set the reference group (non-TNBC) | What the sign of a fold change means |
| 5 | `DESeqDataSetFromMatrix(..., design = ~ instrument + group)`; filter low counts | Design formulas, and correcting for batch (sequencing machine) |
| 6 | VST + PCA | Quality control: do samples group by biology or by batch? |
| 7 | `DESeq()` | Size factors, dispersion, negative binomial test |
| 8 | `results()` | log2 fold change, p-value vs adjusted p-value (FDR) |
| 9 | `lfcShrink()` | Why noisy fold changes are shrunk |
| 10 | MA plot, volcano, single-gene plots, heatmap | Showing DE results |
| 11 | Save the table | Exporting results (and Excel's gene-name pitfall) |
| 12 | Exercises | Interpretation and changing the design |

**Sanity check:** TNBC tumours lack the estrogen receptor, so **ESR1 must come out lower in TNBC** (negative log2 fold
change). If it doesn't, the groups are coded the wrong way round.

**Expected numbers** (with design `~ instrument + group`): about 1,600 genes at padj < 0.05, roughly half higher and
half lower in TNBC. ESR1 log2FC is about −4.2. The ERBB2 (HER2) amplicon genes are near the top of the list.

## Key concepts in one line each

- **Raw counts in, not TPM/FPKM.** DESeq2 normalizes for sequencing depth itself.
- **Reference level.** log2FC = TNBC vs non-TNBC: +1 means 2× higher in TNBC, −1 means 2× lower.
- **padj (FDR), not pvalue.** With about 20,000 tests, about 1,000 genes would reach p < 0.05 by chance alone.
- **Batch.** 9 of the 38 TNBC models were sequenced on a machine that ran no non-TNBC models. Putting `instrument` in
  the design keeps machine differences from being mistaken for TNBC biology.
- **PDX caveat.** Stroma and immune cells in these tumours come from the mouse and were removed, so immune or stromal
  genes can't be studied here.

## Running on the cluster instead of a laptop

The same folder works in RStudio Server (e.g. Open OnDemand) on the cluster. Open the folder there, then follow the
steps above.
