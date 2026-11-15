# CellCoal 2.1.9 — what changed

**Output of all existing runs is bit-identical to 2.1.8, with the same speed.**

## Source split into modules

`cellcoal.c` (11,600 lines) is divided into `tree.c`, `mutation.c`, `observation.c`, `reads.c`, `reads_output.c`, `output.c`, `rng.c`, `options.c`, `config.c` and `globals.c`; the main program stays in `cellcoal.c`. The shared declarations are in `cellcoal.h` (global variables are now `extern` and defined once in `globals.c`). Compiled output, tests (including ASan/UBSan) and timings are the same as before the split.

## Configuration files (--config)

`cellcoal --config file.json|yaml` reads the options from a file. Keys are option letters or long names (`cells`, `coverage`, `amplification_error`, ...), values are numbers, strings, booleans or lists; options after `--config` override the file. `Results/run.json` has a new `options` object, including the seed, that is accepted by `--config`, so `cellcoal --config Results/run.json` repeats a run. An example is in `examples/config.example.yaml`.

## Summary plots

`tools/cellcoal_plots.py` (also installed as `cellcoal-plots`): `depth` (depth per cell, per site and per genotype, with the dispersion statistics), `missing` (genotypes without reads) and `heatmap` (cells x SNVs genotype heatmap).

## Tests

`make test` checks that `--config` with run.json and with YAML reproduces a command line, that options after `--config` override it, that unknown options are errors, and that the plots are written.
