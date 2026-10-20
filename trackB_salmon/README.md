# Track B: reprocess the PDX FASTQs with Salmon (human + mouse + ERCC)

Track B is an independent replication of Track A, which uses the authors' Xenome/RSEM counts. It quantifies the
downloaded FASTQs (`../fastq/`) against **one** Salmon index that holds human GENCODE v50, mouse GENCODE vM39 and the
ERCC spike-in transcripts. Mouse (stroma) reads therefore go to mouse genes instead of being counted as human. Both
genomes serve as decoys.

Run every command from this folder (`cd trackB_salmon`). `config.sh` here loads the main project's `../config.sh`
first, which provides the SLURM account and the FASTQ folder.

| Step | Command | What it does |
|---|---|---|
| once | `scripts/submit.sh scripts/setup_env_rnaseq.slurm` | conda env `tnbc-rnaseq`: salmon, fastqc, multiqc, cutadapt |
| 03 | `scripts/submit.sh scripts/03_download_refs.slurm` | GENCODE human v50 and mouse vM39 transcripts, genomes and GTFs (md5-checked), plus ERCC92 |
| 04 | `scripts/submit.sh scripts/04_build_salmon_index.slurm` | gentrome + decoys + tx2gene → Salmon index (16 CPU / 128 GB / a few h) |
| 05+ | (to write) | FastQC/MultiQC, `salmon quant` array, tximport → human genes, DESeq2 as in Track A |

Run 03 and the env job at the same time, then 04 once both have finished:

```bash
a=$(scripts/submit.sh --parsable scripts/setup_env_rnaseq.slurm)
b=$(scripts/submit.sh --parsable scripts/03_download_refs.slurm)
scripts/submit.sh --dependency=afterok:$a:$b scripts/04_build_salmon_index.slurm
```

Logs go to `logs/`, references to `data/ref/` (not tracked by git), and tool versions to
`metadata/software_versions.txt`.

Status: scripts written, not run yet.
