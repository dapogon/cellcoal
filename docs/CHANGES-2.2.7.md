# CellCoal 2.2.7 — what changed

**Without the new option, the simulated data are identical to 2.2.6.**

- New option `--amp-cell-var V` (`amplification_error_cell_variance` in configuration files): the mean amplification error (`-A`) varies among cells, drawn from a beta distribution with the mean of `-A` and variance `V`; within a cell, sites and haplotypes keep the concentration of `-A`. The mean of every cell is written to the new column `amplification_error` of `cell_params`. It needs `-A` and `V < mean (1 - mean)`. Manual section 6.10c2; test 2g10b. It was contributed as pull request #12.
- Fix: `run.json` did not list the long options (`--fast`, `--bcf`, `--summary`, `--cov-ado-corr`, `--cov-doublet-corr`, `--ado-doublet-corr`, `--amp-cell-var`), so `--config Results/run.json` did not repeat a run that used them. They are now written with the names of the configuration files.
- The README, the manual (section 1.1), `CITATION.cff` and the header of the source state that CellCoal 2 was developed together with Claude, an AI model made by Anthropic, working as a coding agent in collaboration with David Posada.
