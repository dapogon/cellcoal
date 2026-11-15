# CellCoal 2.2.0 — what changed

**Without the new options the output of all existing runs is bit-identical to 2.1.9.**

## Fast mode (`--fast`)

`--fast` (or `fast: true` in a `--config` file) simulates the reads about 10-15% faster: Poisson numbers by inversion/PTRS, binomial numbers of error reads by geometric jumps, and polynomial approximations of the base-10 logarithm and power in the genotype likelihoods (relative error < 1e-12). The distributions are the same, but the random numbers of the reads are consumed differently, so a seed gives different reads. The tree, mutations and true genotypes do not change, and results do not depend on the number of threads. `tests/fast_stats.py` (run by `make test`) compares depth, read counts, error reads and likelihoods with and without `--fast` with Kolmogorov-Smirnov and z tests, in three scenarios.

## ADO rate of every cell (`ado_cells`)

With `-P` the ADO rate drawn for every cell is now written to `ado_cells` (`ado_cells_dir/ado_cells.NNNN` with `-v`), one `[#replicate]` section per replicate with the columns `cell` and `ado_rate` (the outgroup included).

## Correlation between coverage and ADO (`--cov-ado-corr RHO`)

Optional correlation (-1 to 1) between the coverage multiplier of the cells (`-J`) and their ADO rate (`-P`), as a Gaussian copula on the ranks; both distributions are unchanged. Negative: cells with more dropout have less coverage. Also `cov_ado_correlation` in a configuration file.

## Clearer `--config` errors

Errors give the file, line and column and show the line with a caret; a misspelled key suggests the closest option (`unknown option: cell (did you mean 'cells'?)`).

## Plots, CI and packaging

`cellcoal_plots.py all Results` writes every plot and a self-contained `report.html`. CI installs matplotlib and runs the plot tests; a packaging workflow tests the Homebrew formula and the conda recipe.
