# CellCoal 2.7.2 — what changed

## Fixes

- Pileup (`-9`): every read is written in the bases column, also the sequencing and amplification errors to bases that are in no genotype of the site. Before, the column could be shorter than the depth and these errors were undercounted. The counts now match the `RC` field of the VCF.
- Windows: the build did not find zlib, so the BAM, BCF and FASTQ.gz files were written without compression and `--compress-level` had no effect. Now they are compressed as in Linux and macOS.
- `-C` above 1000000 is an error (it overflowed an integer).

## Tests and CI

- The tests run on macOS and Windows (`wc`, `stat`, UTF-8, CRLF, paths); the CI passes in all the platforms, with sanitizers, clang and without OpenMP.
