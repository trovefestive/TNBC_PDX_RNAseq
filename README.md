# TNBC PDX RNA-seq: transcriptomic comparison of triple-negative and non-TNBC breast cancer PDX models

A reanalysis of public total RNA-seq from 53 breast cancer patient-derived xenograft (PDX) models from the Baylor
College of Medicine collection ([GEO GSE183187](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183187)). The
primary question is which human tumour genes and pathways differ between **triple-negative (TNBC, n = 39)** and
**non-TNBC (ER+, HER2+, ER+/HER2+; n = 14)** models.

The pipeline runs on a SLURM HPC cluster. Every step is a numbered, commented script so that it can be rerun and
taught.

## Dataset

| | |
|---|---|
| GEO / BioProject / SRA | GSE183187 / PRJNA759581 / SRP335205 |
| Samples | 53 PDX models (human tumour in the mammary fat pad of SCID/Bg mice) and 15 patient tumours |
| Library | total RNA, NuGEN Ovation RNA-Seq v2 with ERCC spike-ins; paired-end, about 100 bp |
| Sequencers | Illumina HiSeq 2000, HiSeq 2500, NovaSeq 6000 |
| Raw data | 68 runs, 2.3 TB `.fastq.gz` (downloaded and verified; not stored in this repository) |
| Original studies | Petrosyan *et al.* 2022, *NPJ Breast Cancer* ([doi:10.1038/s41523-022-00476-0](https://doi.org/10.1038/s41523-022-00476-0)); Lei *et al.* 2024, *bioRxiv* ([doi:10.1101/2024.12.09.627518](https://doi.org/10.1101/2024.12.09.627518)) |

PDX models by subtype and sequencer:

| Sequencer | TNBC | non-TNBC |
|---|---|---|
| HiSeq 2000 | 22 | 7 |
| HiSeq 2500 | 9 | 0 |
| NovaSeq 6000 | 8 | 7 |

Each PDX library mixes reads from the human tumour and the mouse stroma. Human signal must therefore be separated
from mouse signal before expression analysis.

## Analysis strategy

Two independent quantification tracks feed the same downstream analysis:

| Track | Input | Human/mouse separation | Status |
|---|---|---|---|
| **A** (`scripts/1*`) | Authors' RSEM expected counts deposited on GEO | Xenome (authors) | QC, DE and GSEA done |
| **B** (`trackB_salmon/`) | Our FASTQ download | One Salmon index with human GENCODE v50 + mouse GENCODE vM39 + ERCC, with both genomes as decoys | Quantification, DE, GSEA and concordance done |

The downstream steps are the same for both tracks:

1. **Sample QC.** PCA on VST expression, and how much of each PC is explained by sequencer vs subtype. A subtype
   marker panel (ESR1, PGR, GATA3, FOXA1, ERBB2, GRB7, KRT5, KRT14, FOXC1) checks the labels.
2. **Differential expression.** DESeq2 with `~ instrument + group`, using non-TNBC as the reference, and
   apeglm-shrunk fold changes. Filter: ≥10 counts in ≥3 models. Significance: BH FDR < 0.05.
3. **Robustness checks:**
   - drop the HiSeq 2500 batch;
   - no sequencer covariate;
   - fits within a single sequencer;
   - include or exclude a low-depth model.

   Exploratory contrasts: TNBC vs ER+ and TNBC vs HER2+.
4. **GSEA.** fgsea on the DESeq2 Wald statistic against MSigDB 2026.2 GO Biological Process, Reactome and
   WikiPathways. Hallmark sets serve as a positive control.
5. **Track A vs Track B concordance.** Fold-change and NES correlation, plus overlap of significant genes and
   pathways.

**Sign convention:** every log2 fold change and NES is **TNBC − non-TNBC**, so positive means higher in TNBC. Each
result file states this in its header. As a guard, the DE script stops unless ESR1 comes out lower in TNBC.

## Results so far

### Track A: authors' counts

| | |
|---|---|
| Models analysed | 38 TNBC vs 14 non-TNBC. BCM-3807 is excluded because it has only 1.4 M counts in the authors' output; it is tested in a sensitivity analysis. |
| Genes tested | 21,428 |
| Positive control | ESR1 log2FC = −4.19 (lower in TNBC) |
| DE genes (FDR < 0.05) | **1,572** (744 higher, 828 lower in TNBC); 1,244 with \|log2FC\| > 1 |
| Label check | 1 of 53 models (BCM-7441, labelled ER+) has low ESR1 and high ERBB2 expression |

Robustness of the main result. Columns 4 and 5 count the 1,572 main-model DE genes that the alternative analysis
reproduces (FDR < 0.05) or matches in direction.

| Analysis | TNBC vs non-TNBC | DE genes | Main DE genes reproduced | Same direction | Spearman ρ (all genes) |
|---|---|---|---|---|---|
| Drop HiSeq 2500 batch | 29 v 14 | 1,773 | 1,475 | 100% | 0.98 |
| No sequencer covariate | 38 v 14 | 3,247 | 1,283 | 99.3% | 0.81 |
| HiSeq 2000 only | 21 v 7 | 324 | 267 | 96.7% | 0.78 |
| NovaSeq 6000 only | 8 v 7 | 1,204 | 584 | 97.7% | 0.74 |
| Include BCM-3807 | 39 v 14 | 1,497 | 1,440 | 100% | 0.99 |

The top DE genes lower in TNBC are the **17q12 ERBB2 amplicon** (ERBB2, MIEN1, STARD3, PGAP3, ORMDL3, MED1), which
reflects the HER2+ models in the comparator.

**GSEA** (fgsea, ranked by the DESeq2 Wald statistic; BH padj < 0.05):

| Collection | Gene sets significant | Higher in TNBC | Lower in TNBC |
|---|---|---|---|
| GO Biological Process | 243 | 173 | 70 |
| Reactome | 101 | 80 | 21 |
| WikiPathways | 49 | 43 | 6 |
| Hallmark (control) | 22 | 16 | 6 |

- **Positive control passes:** estrogen response is strongly lower in TNBC. Hallmark early: NES −2.46,
  padj 3e−14. Late: NES −2.33, padj 6e−12.
- **Higher in TNBC:**
  - proliferation (Hallmark E2F targets, G2M checkpoint, MYC targets V1/V2; WikiPathways retinoblastoma gene in
    cancer);
  - epithelial–mesenchymal transition (Hallmark NES 2.56);
  - collagen formation and fibril assembly (Reactome NES about 2.5). The stroma in these tumours is mouse, so these
    human-gene signals come from the tumour cells themselves.
- **Lower in TNBC:**
  - luminal mammary epithelial lineage (Reactome);
  - 17q12 copy-number region (WikiPathways NES −3.03), i.e. the ERBB2 amplicon;
  - interferon-α response;
  - insulin-secretion/secretory GO terms.

Full tables are in [`results/trackA/gsea/`](results/trackA/gsea/).

Figures are in [`results/trackA/plots/`](results/trackA/plots/) and tables in [`results/trackA/`](results/trackA/).

### Track B: Salmon reprocessing of all 53 models

| | |
|---|---|
| Quantification | 53 of 53 models. Library type is unstranded (IU); BCM-3807 runs as single-end. Median 70% of reads mapped. |
| Species split | Mouse (stroma) share is a median 8% of mapped reads in TNBC and 10% in non-TNBC. The difference is not significant after adjusting for sequencer (p = 0.65), so it is not used as a covariate. |
| ERCC spike-ins | Essentially absent (median < 0.1% of reads), so they are not usable for QC. |
| Models analysed | 39 TNBC vs 14 non-TNBC. BCM-3807 is included: Salmon assigned it 125 M human counts. |
| Genes tested | 46,750 GENCODE human genes |
| Positive control | ESR1 log2FC = −4.3 |
| DE genes (FDR < 0.05) | **2,139** (982 higher, 1,157 lower in TNBC). 2,055 remain after removing 84 "sporadic" genes (detected in < 50% of models in both groups; their raw fold changes are unreliable, and apeglm shrinks them to ~0). |
| Robustness | Dropping the HiSeq 2500 batch or BCM-3807 gives ρ = 0.98–0.99, and 100% of main DE genes keep the same direction. |
| Top genes | The 17q12 ERBB2 amplicon again (STARD3, ERBB2, ORMDL3, MIEN1, PGAP3, MED1, CDK12), plus ANKRD30A and MUC6 lower in TNBC |

**Agreement with Track A.** The two tracks used different read assignment (Xenome vs a joint Salmon index) and
different annotation (RefSeq vs GENCODE). 17,842 genes are matched by name between the two tables.

| | |
|---|---|
| log2FC correlation | Pearson 0.89, Spearman 0.91 |
| Significant-gene overlap | 72.5% of Track A DE genes are also significant in Track B. All 988 genes significant in both have the same direction. |
| GSEA NES correlation | GO BP 0.86, Reactome 0.88, WikiPathways 0.87, Hallmark 0.97. All gene sets significant in both tracks have the same direction. |

Track B reproduces the Track A pathway picture:
- **Higher in TNBC:** EMT (NES 2.58), E2F targets (2.54), G2M checkpoint (2.49), MYC targets.
- **Lower in TNBC:** estrogen response early (−2.56) and late (−2.48).

Track B tables and figures are in [`trackB_salmon/results/`](trackB_salmon/results/).

### Caveats

- **Sequencer is partly confounded with subtype.** The HiSeq 2500 batch contains only TNBC, so the subtype effect is
  estimated within HiSeq 2000 and NovaSeq models.
- **The non-TNBC group is mixed.** It combines ER+ and HER2+ biology, so many DE genes reflect the ER or HER2
  programmes. That is expected, not a discovery.
- **No human stroma or immune cells.** In PDX tissue those cells come from the mouse, so stromal and immune
  pathways cannot be compared using human genes.

## Repository layout

```
├── README.md                     project overview (this file)
├── docs/download_guide.md        step-by-step guide to downloading the data (written for new students)
├── config.sh                     shared settings (paths, environments, download method)
├── config.local.sh.example       template for site-specific SLURM settings (the real file is not tracked)
├── scripts/
│   ├── submit.sh                 sbatch wrapper that adds account/partition from config.local.sh
│   ├── setup_env*.slurm          conda environments (download tools; R/Bioconductor)
│   ├── 00_get_run_list.sh        GEO + ENA metadata → run sheet
│   ├── 01_download_fastq.slurm   SLURM array: SRA (prefetch/fasterq-dump/pigz) or ENA, with fallback
│   ├── 02_check_downloads.sh     verifies files, checksums and read counts
│   ├── 10_get_author_counts.sh   Track A: authors' RSEM counts + MSigDB gene sets
│   ├── 11_prepare_author_counts.R  sample sheet and 53-model count matrix
│   ├── 12_deseq2_trackA.R        QC, DESeq2, sensitivity analyses, plots
│   ├── 13_gsea_trackA.R          fgsea (GO BP, Reactome, WikiPathways, Hallmark)
│   └── run_R.slurm               runs any R script as a SLURM job
├── teaching_deseq2/              self-contained RStudio exercise: DESeq2 TNBC vs non-TNBC for students
├── metadata/                     run sheet, PDX sample sheet, download checks, software versions
├── results/trackA/               DE tables, GSEA tables, QC tables, figures
└── trackB_salmon/                Track B: Salmon reprocessing (own config, scripts, README and results)
```

Raw FASTQs, reference files, intermediate matrices, logs and site settings are excluded by `.gitignore`.

## Reproducing the analysis

```bash
cp config.local.sh.example config.local.sh
scripts/submit.sh scripts/setup_env.slurm
scripts/submit.sh scripts/setup_env_r.slurm
bash scripts/00_get_run_list.sh
scripts/submit.sh --array=1-68%6 scripts/01_download_fastq.slurm
bash scripts/02_check_downloads.sh
bash scripts/10_get_author_counts.sh
Rscript scripts/11_prepare_author_counts.R
scripts/submit.sh scripts/run_R.slurm scripts/12_deseq2_trackA.R
scripts/submit.sh scripts/run_R.slurm scripts/13_gsea_trackA.R
```

1. **Setup.** Copy the settings template, then set your SLURM account and partition in `config.local.sh`. The two
   `setup_env` jobs build the conda environments.
2. **Download the raw data** (optional for Track A): `00_get_run_list.sh` writes the run sheet, the array job
   downloads 68 runs (2.3 TB), and `02_check_downloads.sh` verifies them.
3. **Track A.** `10_get_author_counts.sh` and `11_prepare_author_counts.R` are small and run on the login node. The
   last two commands submit DESeq2 and GSEA as jobs.

[`docs/download_guide.md`](docs/download_guide.md) explains each download step, how to monitor jobs, and how to
recover from failures. Tool versions are recorded in
[`metadata/software_versions.txt`](metadata/software_versions.txt).

## Teaching material

[`teaching_deseq2/`](teaching_deseq2/) is a stand-alone RStudio project for undergraduates. It contains the count
matrix (52 models × 25,143 genes) and a sample sheet, both built by `scripts/14_make_teaching_files.R`. A commented
DESeq2 script runs TNBC vs non-TNBC from reading the data through PCA, DESeq2, shrinkage, volcano plot and heatmap,
and ends with exercises. It reproduces the Track A result: 1,572 DE genes, ESR1 log2FC −4.19.

## Data and citation

All data come from GEO GSE183187. Please cite the original studies above when using these data.
