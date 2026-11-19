# Unreleased (after 2.2.6)
- Documentation: the README, the manual, `CITATION.cff` and the source header state that CellCoal 2 was developed together with Claude (Anthropic).
- New option `--amp-cell-var V` (`amplification_error_cell_variance` in configuration files): the mean amplification error (`-A`) varies among cells (beta distribution with the mean of `-A` and variance `V`); the mean of every cell is written to `cell_params`. Without the option the output is unchanged. Manual section 6.10c2; test 2g10b.
