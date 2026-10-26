# Running the DESeq2 exercise on the cluster (interactive session)

This guide runs the same analysis as the RStudio version (`DESeq2_TNBC_vs_nonTNBC.R`), but on the university's
high-performance computing (HPC) cluster, from a terminal. You will learn the three things every cluster user needs:

1. **Where to work:** the login node vs a compute node.
2. **How to get a compute node:** an *interactive job*.
3. **How to set up software:** modules and a conda environment.

The R analysis itself is identical. Read the comments in the script exactly as you would in RStudio.

> Replace the words in `<angle brackets>` with your own values. Your instructor will give you the
> **allocation** (the account that pays for compute time) and the **login address**.

---

## 0. The one rule: don't compute on the login node

When you log in, you land on a **login node**, a shared machine where hundreds of people edit files and submit
work. Running R analyses there slows it down for everyone, and the system may kill your program.

Real work happens on **compute nodes**. You ask for one, use it, and give it back:

```
your laptop ──ssh──> login node ──ijob──> compute node   (R runs here)
                     (files, editing,     (CPUs and memory reserved for you)
                      asking for nodes)
```

## 1. Log in

From a terminal (macOS/Linux: *Terminal*; Windows: *PowerShell* or the cluster's web portal *Shell*):

```bash
ssh <computing_id>@<cluster-login-address>
```

## 2. Get the exercise files (once)

Work in your **scratch** folder: it is large and fast. Note that scratch is cleaned out after a period of inactivity,
so copy anything you want to keep.

```bash
cd /scratch/$USER
git clone https://github.com/trovefestive/TNBC_PDX_RNAseq.git
cd TNBC_PDX_RNAseq/teaching_deseq2
ls
```

You should see `DESeq2_TNBC_vs_nonTNBC.R`, `run_section.R`, `environment.yml`, this guide and the `data/` folder.

## 3. Start an interactive job

```bash
ijob -A <allocation> -p interactive -c 2 --mem=16G -t 02:00:00
```

| Option | Meaning |
|---|---|
| `-A <allocation>` | the account charged for the compute time |
| `-p interactive` | the *partition* (queue) for interactive work |
| `-c 2` | 2 CPU cores |
| `--mem=16G` | 16 GB of memory (DESeq2 on 52 samples needs a few GB) |
| `-t 02:00:00` | the time limit (2 h). When it runs out, the session ends, so save your work. |

`ijob` waits until a node is free (usually seconds to minutes), then logs you into it. **Check where you are:**

```bash
hostname          # a compute-node name now, not the login node
squeue -u $USER   # your running job, with its JOBID and time left
```

If your cluster has no `ijob` command, `salloc` takes the same options:
`salloc -A <allocation> -p interactive -c 2 --mem=16G -t 02:00:00`.

## 4. The software: modules and a conda environment

The cluster offers software as **modules**, which you switch on with `module load`. The R packages for this exercise
(DESeq2, apeglm, ggplot2, ggrepel, pheatmap) live in a **conda environment**: a folder holding exactly these
software versions.

**Option A, in class (recommended):** your instructor has built the environment once for everyone. You will be given
its path; nothing to install.

**Option B, on your own:** build your own copy from `environment.yml`. This takes about **20–25 minutes** and
**1.6 GB** of your home folder, and it must run **inside the interactive job**. You only need to do it once:

```bash
module load miniforge
conda env create -f environment.yml      # creates an environment called tnbc-deseq2
```

> **Instructors:** build one shared copy before class, inside a job, in a folder students can read:
> `conda env create -p <shared_folder>/tnbc-deseq2 -f environment.yml`, then
> `chmod -R a+rX <shared_folder>/tnbc-deseq2`. Students then activate it by its path (step 5).

## 5. Start R

Every session (each new `ijob`) starts from a clean shell, so run these lines each time:

```bash
module load miniforge
source activate <shared_folder>/tnbc-deseq2   # Option A: the path from your instructor
# source activate tnbc-deseq2                 # Option B: your own environment
export OPENBLAS_NUM_THREADS=1   # keeps the maths library from grabbing more cores than you asked for
```

Then make sure you are in the exercise folder and start R:

```bash
cd /scratch/$USER/TNBC_PDX_RNAseq/teaching_deseq2
R
```

## 6. Run the analysis one section at a time

There is no RStudio here, so you can't select lines and press Ctrl+Enter. Instead, a small helper runs one numbered
section of the script and prints every line with its output, just as RStudio would:

```r
source("run_section.R")
sections()     # lists the 12 sections
section(1)     # packages
section(2)     # read the data
section(3)     # ...and so on up to section(11)
```

After each section, **stop and read**: the comments explain what happened and ask questions. You can also type your
own commands between sections, because everything stays in your R workspace. For example:

```r
dim(counts)
head(res)
res["GATA3", ]
```

To read the script itself, open a second terminal (log in again; you can stay on the login node just to read) and
run `less DESeq2_TNBC_vs_nonTNBC.R`. Press `q` to quit `less`.

Run the sections **in order**. Each one uses objects made by the ones before it. The whole analysis takes about
2 minutes of computer time; most of that is section 7 (`DESeq()`) and section 9 (shrinkage), each about a minute.
The `>` prompt comes back when a section has finished.

## 7. Look at your plots

A terminal can't show pictures. The script saves its plots as files in `results/`:
`PCA.png`, `volcano.png` and `heatmap_top30.png`. Plots that are only drawn on screen in RStudio (the MA plot and
the single-gene plots) go to `Rplots.pdf` here instead.

To view them, either:
- open the cluster's **web portal → Files**, browse to the folder and click the image; or
- copy them to your laptop. Run this **on your laptop**, not on the cluster:

  ```bash
  scp "<computing_id>@<cluster-login-address>:/scratch/<computing_id>/TNBC_PDX_RNAseq/teaching_deseq2/results/*.png" .
  ```

## 8. Check your result

You should get the same numbers as the RStudio version:
- about **1,572** genes with padj < 0.05 (744 higher, 828 lower in TNBC);
- **ESR1 log2FC ≈ −4.19**, lower in TNBC. This is the sanity check: if ESR1 comes out *positive*, the groups are
  coded the wrong way round.

## 9. Finish: give the node back

```r
q()        # leave R; answer n (no need to save the workspace)
```

```bash
exit       # leave the compute node and end the job
squeue -u $USER   # should now be empty
```

Ending the job when you are done matters. The allocation is charged for the time you *hold* the node, not just the
time you use it.

---

## Bonus: the same analysis as a batch job

An interactive job is for exploring. Once a script works, you normally submit it as a **batch job**: it runs on its
own, even after you log out. Save this as `run_deseq2.slurm` in the exercise folder:

```bash
#!/usr/bin/env bash
#SBATCH -A <allocation>
#SBATCH -p standard
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH -o deseq2_%j.log
set -euo pipefail
module load miniforge
source activate tnbc-deseq2
export OPENBLAS_NUM_THREADS=1
Rscript DESeq2_TNBC_vs_nonTNBC.R
```

```bash
sbatch run_deseq2.slurm          # prints "Submitted batch job <JOBID>"
squeue -u $USER                  # PD = waiting, R = running
cat deseq2_<JOBID>.log           # everything R printed
sacct -j <JOBID> --format=JobID,State,Elapsed,MaxRSS   # did it finish, and how much memory did it use?
```

**Question:** what are the advantages of each way of working? When would you choose an interactive job over a
batch job?

## Troubleshooting

| Problem | Fix |
|---|---|
| `conda: command not found` or `source: activate not found` | Run `module load miniforge` first. |
| `there is no package called 'DESeq2'` | The environment isn't active: run the `source activate ...` line from step 5. |
| `object 'counts' not found` (or `info`, `dds`, `res` ...) | You skipped a section, or R was restarted. Run the sections in order from `section(1)`. |
| `OpenBLAS error: Memory allocation still failed` | You forgot `export OPENBLAS_NUM_THREADS=1`, or you are on the login node (check with `hostname`). |
| `cannot open file 'data/...': No such file` | You started R in the wrong folder. In R: `getwd()`, then `setwd("/scratch/<computing_id>/TNBC_PDX_RNAseq/teaching_deseq2")`. |
| `ijob` waits a long time | The interactive partition is busy. Wait, or ask for less (`-c 1 --mem=8G`). |
| The session closed by itself | The time limit (`-t`) ran out, or your connection dropped. Start a new `ijob` and re-run the sections from `section(1)`; it takes about 2 minutes. |
| `Disk quota exceeded` while creating the environment | Conda stores packages in your home folder. Ask your instructor; `conda clean --all` frees space. |

## Words you met

| Term | Meaning |
|---|---|
| login node | the shared entry machine; for files and editing only |
| compute node | a machine whose CPUs/memory are reserved for your job |
| SLURM | the scheduler that hands out compute nodes |
| allocation / account | the budget your compute time is charged to |
| partition | a queue of nodes with its own rules (`interactive`, `standard`, ...) |
| interactive job (`ijob`, `salloc`) | a compute node you type on directly |
| batch job (`sbatch`) | a script that runs on a compute node without you |
| module | a ready-installed program you switch on with `module load` |
| conda environment | your own set of software versions, switched on with `source activate` |
