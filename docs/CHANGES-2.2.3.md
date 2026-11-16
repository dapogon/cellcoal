# CellCoal 2.2.3 — what changed

**The output of all existing runs is bit-identical to 2.2.2.**

## New options

- `--cov-doublet-corr RHO` and `--ado-doublet-corr RHO`: correlation of being a doublet (`-B`) with the coverage (`-J`) and the ADO rate (`-P`) of the cells. The coverage, ADO rate and doublet status of every cell are written to `cell_params`.
- `--summary`: `summary.tsv`, one row of metrics per replicate.
- `--sweep "C=2,5,10 E=0.001,0.01"`: runs a grid of parameter values with a common seed and writes `sweep.tsv` (mean and standard error per combination) and `sweep_replicates.tsv`; `--sweep-keep` keeps the run folders. Also available as `sweep:` in `--config` files. See section 6.10i of the manual.

## Fix

- `--ado-doublet-corr` without `--cov-ado-corr` crashed (segmentation fault) because the ADO rates of the cells were not kept.

## Code and checks

- `main` is split into functions (default options, option combinations, summary statistics, streamed reads, per-replicate output); no change in the output.
- The CI has a `speed` job (`tests/bench.py`): it builds the latest release and the new code on the same runner and fails if a workload is more than 25% slower.
