# CellCoal 2.4.1 — what changed

**The simulated data of all existing runs are identical to 2.4.0.**

## New
- `--amp-cell-var` now also works in reads mode: every cell draws its own mean amplification error (the concentration of `-A` is kept) and the molecules of the BAM/FASTQ follow it (test 2g40). Default data are unchanged.
- `--blocks`/`--fasta`/`--germline-vcf` now also work in reads mode (`--bam`, `--fastq`): one BAM and FASTQ pair per cell along the whole chromosome, written block by block; `ref.fa`, `.fai` and `truth.vcf` in chromosome coordinates. In a chromosome in blocks the per-cell coverage multipliers and amplification means (`-J`, `--amp-cell-var`) are now drawn once and shared by all the blocks (before, every block redrew them). Default data are unchanged (reads without blocks are bit-identical to the previous commit). See manual 8.29, test 2g41.
- `--qual-sd S`, `--qual-slope B` (reads mode): variable, calibrated base qualities (error of a base = -E x lognormal x position factor; Q written = -10 log10 p). `--strand-bias P` (reads mode): a fraction P of the amplification errors is seen only on one strand. `tools/benchmark.py pileup` gets `--min-bq` and `--both-strands`. Reads mode now stops with an explanation when dropout is not in segments (`-D` without `--ado-block` left almost empty BAM files). Default data are unchanged. Manual 8.30, test 2g42.
