# CellCoal 2.2.4 — what changed

**The output of all existing runs is bit-identical to 2.2.3.**

## Truth metrics and better sweeps

- `--summary` adds `called_fraction`, `genotype_error`, `snv_precision` and `snv_recall`, computed against the true genotypes (`NA` without reads); stored and streamed reads give the same values.
- `--sweep-jobs N` runs N combinations at the same time (same results as running them one by one); `--sweep-resume` continues an interrupted sweep; `sweep.tsv` handles `NA` metrics; values can be ranges (`LOW..HIGH..STEP`).
- `cellcoal_plots.py sweep` draws the table of a sweep.

## Speed

- Vectorized genotype likelihoods and a lighter dense memory layout: runs without reads are about 35% faster; runs with reads are unchanged.
- The CI has a speed check (`tests/bench.py`) against the latest release.

## Windows

- The program builds and runs on Windows with MSYS2/MinGW-w64, and the CI runs the whole test suite there (see section 6.10j of the manual). The random number generator and the index arithmetic use 64-bit integers explicitly because `long` is 32 bits on Windows; configuration files are read correctly with CRLF line endings.

## Internal

- `main` is split further (`RunReplicate`), the statistics accumulated over replicates are in a `RunStats` structure (about 30 fewer globals), and `RemoveDir` no longer needs `nftw`.
