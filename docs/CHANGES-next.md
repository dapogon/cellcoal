# Unreleased (after 2.2.3)

- Internal: `main` is split further (`RunReplicate`), the statistics accumulated over replicates are in a `RunStats` structure (about 30 fewer globals), and `RemoveDir` no longer needs `nftw`.
- Windows (MSYS2/MinGW-w64): the program builds and runs; the random number generator and the index arithmetic use 64-bit integers explicitly (`long` is 32 bits on Windows). New CI job `windows` (build, smoke tests, whole suite, all checks pass). See section 6.10j of the manual.

The output on Linux and macOS is identical to 2.2.3.

## Truth metrics and better sweeps

- `--summary` adds `called_fraction`, `genotype_error`, `snv_precision` and `snv_recall` against the true genotypes (`NA` without reads); stored and streamed reads give the same values.
- `--sweep-jobs N` runs N combinations at the same time; `--sweep-resume` continues an interrupted sweep (the summaries of finished runs are kept in `<dir>/runs`); `sweep.tsv` handles `NA` metrics.
- `cellcoal_plots.py sweep` draws the table of a sweep.
- Configuration files (`--config`) are read robustly on Windows (CRLF in text mode).
