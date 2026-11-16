# Unreleased (after 2.2.2)

- `--cov-doublet-corr RHO` and `--ado-doublet-corr RHO`: correlation of being a doublet (`-B`) with the coverage (`-J`) and the ADO rate (`-P`) of the cells; the per-cell values are written to `cell_params`.
- `--summary`: `summary.tsv` with one row of metrics per replicate.
- `--sweep "C=2,5,10 E=0.001,0.01"`: runs a grid of parameter values with a common seed and writes `sweep.tsv` (mean and standard error of the metrics per combination) and `sweep_replicates.tsv`; `--sweep-keep` keeps the run folders. See section 6.10i of the manual.

The default output is identical to 2.2.2.
