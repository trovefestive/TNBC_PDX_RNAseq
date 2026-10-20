# Step-by-step guide: downloading GEO GSE183187

This guide walks through downloading the raw RNA-seq reads (FASTQ files) of GEO series
[GSE183187](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183187):
*Total RNA-seq from Patient-Derived Xenograft models of Breast Cancer developed at
Baylor College of Medicine*. The scripts are kept short and commented so that a
new student can read each one top to bottom.

## The dataset

| | |
|---|---|
| GEO series | GSE183187 |
| BioProject / SRA study | PRJNA759581 / SRP335205 |
| Samples = runs | 68 (one SRR run per GSM sample) |
| Layout | paired-end (2 files per run: `_1` and `_2`), except SRR15687248, which is mostly single-end (see Notes) |
| Library | total RNA, NuGEN Ovation RNA-Seq v2 + ERCC spike-ins; Illumina HiSeq 2000/2500, NovaSeq 6000 |
| Total download | **about 2.3 TB** of `.fastq.gz` |

| Sample type | Runs | Size |
|---|---|---|
| PDX, TNBC | 39 | 1.08 TB |
| PDX, ER+ | 6 | |
| PDX, HER2+ | 6 | |
| PDX, ER+/HER2+ | 2 | |
| Patient tumour (no mouse) | 15 | 0.85 TB |

**Human and mouse reads are mixed.** Each PDX run is human tumour grown in the
mammary fat pad of a SCID/Bg mouse. Its reads come from both the human tumour and
the mouse stroma, all in the same FASTQ file. GEO lists both organisms for these
samples, while SRA labels every run "Homo sapiens". The authors split the reads
with Xenome before aligning to hg38 and mm10. We have to do the same kind of split
(for example Xenome, XenofilteR or Disambiguate) before any analysis. The 15
patient-tumour runs are human only.

The full sample sheet is `metadata/GSE183187_runinfo.tsv`, one row per run:

| Column | Meaning |
|---|---|
| `gsm` | GEO sample ID |
| `sample_title` | GEO sample title, e.g. "Patient Derived Xenograft model BCM-0002" |
| `srr` | SRA run ID: the name of the FASTQ files |
| `sample_type` | `PDX` or `patient_tumor` |
| `tumor_subtype` | TNBC, ER+, HER2+, ER+/HER2+ (NA for patient tissue) |
| `mouse_host`, `geo_organism` | host strain; organisms listed in GEO |
| `layout` | PAIRED or SINGLE |
| `spots` | read pairs (one spot = one read pair) |
| `reads` | individual reads (2 × spots for paired-end data) |
| `bases`, `instrument` | as reported by SRA |
| `n_files`, `fastq_gb` | number and total size of the FASTQ files |
| `fastq_url`, `fastq_md5`, `fastq_bytes` | ENA download links, checksums and sizes (`;`-separated) |

## Folder layout

```
TNBC_PDX_RNAseq/
├── README.md               <- project overview
├── docs/download_guide.md  <- you are here
├── config.sh               <- settings used by every script (paths, env name, download method)
├── config.local.sh         <- YOUR SLURM account/partition (not in git; see .example)
├── scripts/
│   ├── submit.sh               wrapper around sbatch (adds account/partition)
│   ├── setup_env.slurm         one-time: conda env with sra-tools + pigz
│   ├── 00_get_run_list.sh      step 0: build the sample sheet (login node OK)
│   ├── 01_download_fastq.slurm step 1: download, one array task per run (SLURM)
│   └── 02_check_downloads.sh   step 2: check that everything is there and correct
├── metadata/               <- sample sheet, run list, check results, software versions
├── fastq/                  <- the downloaded *.fastq.gz (+ fastq/.checks/ markers)
└── logs/                   <- one SLURM log file per job / array task
```

## Cluster rules (read first)

- **The login node is shared.** Only run tiny things there: `ls`, `head`, and
  `00_get_run_list.sh`, which fetches a few kB of metadata. Downloads, md5 checks
  and read counting go through `sbatch`.
- Always run commands **from the project folder** (the top folder of this repository).
- Submit with `scripts/submit.sh`, not plain `sbatch`. The wrapper loads
  `config.local.sh`, which holds the SLURM account and partition, so the scripts
  never contain an allocation name.

## Step by step

### Once: site settings and software

```bash
cp config.local.sh.example config.local.sh
nano config.local.sh
scripts/submit.sh scripts/setup_env.slurm
```

In `config.local.sh`, put your allocation in `SBATCH_ACCOUNT` and your partition
in `SBATCH_PARTITION`. The `setup_env.slurm` job builds a conda env, records the
tool versions in `metadata/software_versions.txt`, and takes a few minutes.

### Step 0: get the list of runs (login node, a few seconds)

```bash
bash scripts/00_get_run_list.sh
```

This writes `metadata/GSE183187_runinfo.tsv` (the sample sheet) and
`metadata/runs_to_download.txt` (one SRR per line, all 68 runs). It prints a
summary that includes the total size.

**You only want some runs?** Make a shorter list, then use its length in the
`--array` range of step 1. For example, TNBC PDX only (39 runs):

```bash
awk -F'\t' '$5 == "TNBC" {print $3}' metadata/GSE183187_runinfo.tsv > metadata/runs_to_download.txt
```

### Step 1: download (SLURM array job)

```bash
wc -l metadata/runs_to_download.txt
scripts/submit.sh --array=1-68%6 scripts/01_download_fastq.slurm
```

The first command shows how many runs are in the list (68).

- `--array=1-68` starts 68 tasks. Task N downloads line N of the run list.
- `%6` means at most 6 tasks run at the same time. This is polite to ENA and
  to the shared file system.
- How each task downloads is set by `DOWNLOAD_METHOD` in `config.sh`:
  - `sra` (default): `prefetch` fetches the compact `.sra` file from NCBI, which
    stores it on Amazon cloud storage. `prefetch` checks its md5 and
    `vdb-validate` checks it again. `fasterq-dump` writes `_1`/`_2` FASTQ, and
    `pigz` compresses them.
  - `ena`: downloads ENA's ready-made `.fastq.gz` over HTTPS and checks each
    file against ENA's md5. On a mismatch it deletes the file and retries.
  - If the chosen method fails, the task tries the other one.
  - Why `sra` is the default: in a speed test on a compute node (Oct 2026), NCBI
    gave about 55 MB/s and ENA only about 1.4 MB/s. At ENA's speed 2.3 TB takes
    days.
- Every task then counts the reads (one per spot in each of `_1` and `_2`) and
  writes `fastq/.checks/<SRR>.done`.
- While it works, a task needs temporary space in `fastq/tmp_<SRR>/`: up to
  about 3× the run's size for the uncompressed FASTQ. The space is freed when the
  task finishes.
- One run is 11–67 GB of `.fastq.gz`. Expect roughly 1–3 h per task with the `sra` method (download + convert + compress + count).

To test on a single small run first, submit only its line number:

```bash
scripts/submit.sh --array=4 scripts/01_download_fastq.slurm
```

### Watching your jobs

These commands are run one at a time:

| Command | What it shows |
|---|---|
| `squeue -u $USER` | jobs that are waiting (PD) or running (R) |
| `sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS,ExitCode` | how each task ended (COMPLETED / FAILED / TIMEOUT) |
| `tail -f logs/01_download_<jobid>_<task>.log` | follow one task's log live (Ctrl-C to stop) |
| `ls fastq/.checks/ \| wc -l` | how many runs are finished |
| `scancel <jobid>` | cancel a whole job |
| `scancel <jobid>_<task>` | cancel one task |
| `seff <jobid>_<task>` | CPU / memory efficiency of a finished task |

Logs are in `logs/`. Step 1 writes one per task, `01_download_<jobid>_<task>.log`,
and the log is the first place to look when something fails.

### Step 2: check the downloads

```bash
bash scripts/02_check_downloads.sh
scripts/submit.sh scripts/02_check_downloads.sh --md5
```

The first command is the quick check: seconds, login node OK. The second is
optional. It recomputes every md5, which takes hours, so it runs as a job.

For every run, step 2 checks four things:

- the `.done` marker exists;
- both files exist with ENA's exact byte size;
- the md5 was verified;
- each file holds one read per spot.

It writes `metadata/download_check.tsv` and ends with `PASS` or `FAIL`. Failed
runs are listed together with the command to redo them.

### Something failed?

Just resubmit. The scripts are **safe to rerun**: finished runs are skipped
(they have a `.done` marker) and half-downloaded files are resumed (`wget -c`).
To redo one run:

```bash
grep -n SRR15687233 metadata/runs_to_download.txt
scripts/submit.sh --array=4 scripts/01_download_fastq.slurm
```

The first command finds the run's line number (here 4). Then submit only that
task number. To redo a run that already has a `.done` marker, delete
`fastq/.checks/<SRR>.done` and its files first.

## Notes

- Scratch space is purged after a period of inactivity. Copy anything you need
  to keep, such as the metadata or processed results, to project storage.
- `.gitignore` is a whitelist. Only the scripts, README and metadata tables can
  be committed. FASTQs, logs and `config.local.sh` never are.
- Software versions are in `metadata/software_versions.txt`.
- **SRR15687248 (GSM5552651, BCM-3807) is mostly single-end.** It was deposited
  that way, and ENA has the same 3 files: 7.7 M read pairs (`_1`/`_2`) plus
  173.6 M unpaired 101 bp reads (`SRR15687248.fastq.gz`). Treat it as single-end
  (use the unpaired file, or all 3 files) or leave it out. All other 67 runs are
  normal paired-end.
- Download status (2026-10-20): all 68 runs downloaded (2.3 TB). `02_check_downloads.sh`
  gives PASS for all 68: files present, `.sra` validated, read counts equal SRA's.
- FASTQ files made by `fasterq-dump` name reads `SRRxxxx.1`, `SRRxxxx.2`, ...
  instead of the original Illumina read names. For RNA-seq analysis this does
  not matter.
