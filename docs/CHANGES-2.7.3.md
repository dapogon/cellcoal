# CellCoal 2.7.3 — what changed

- Indels in the reads mode: `--indel-max` up to 100; `--indel-tail` and `--indel-tail-exp` (power-law lengths, long tail); `--indel-repeat` (indels enriched in homopolymers and short tandem repeats, as slippages of one unit); `--germline-indel-rate` (germline indels of the healthy genome, in all the cells and in the bulk, flag `GERMLINE` in `indels_truth.vcf`). Tests 2g64.
- Reads mode: `--mappability FILE` (BED with a score from 0 to 1: the MAPQ of the BAM records follows the mappability) and `--map-drop F` (fraction of the pairs with MAPQ 0 left out of the BAM; the FASTQ files keep every read). Manual Fig. 16, tests 2g65.
- Faster reads mode: the BGZF writer reuses the zlib stream of the thread for every block and computes the CRC with zlib; the BAM, BCF and FASTQ.gz files are identical byte by byte and the reads mode is about 23% faster (10 cells of 1 Mb at 30X: 10.5 s to 8.2 s).
- Names of folders and files: up to 1023 characters (was 119).

## Code

- `cellcoal.c` is split into `validate.c` (checks of the options), `replicate_output.c` (output of every replicate and the summary) and `runlog.c` (run record); the BAM writer is in `bam.c`. Same results for the same seed and options.
