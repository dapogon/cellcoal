# CellCoal 2.6.0 — what changed

**Output of all existing runs is bit-identical to 2.5.0.**

## Run record

Every run writes `Results/run.json`: program version, command line, seed, number of replicates, cells and sites, threads, the `CELLCOAL_*` and `OMP_NUM_THREADS` environment variables that were set, the genome storage (dense or sparse) and read mode (stored or streamed) used, compiler, OpenMP, and elapsed and CPU seconds. With the version, the command line and the seed, it is what is needed to repeat a run.

## Elapsed time

The time printed at the end of a run is now the elapsed (wall clock) time, followed by the CPU time; before, it was the CPU time summed over all the threads, which overstated the duration of parallel runs. It is also kept as a double instead of a float.

## Build

The object files now depend on the header files, so changing a header rebuilds the program (before, an old object could be linked with a changed header or version string).

## Manual

The manual now has a single source in Markdown, `docs/manual/cellcoal-manual.md` (figures in `docs/manual/media`). `python3 tools/build_manual.py` builds the self-contained HTML (`docs/cellcoal.manual.v2.0.html`) with pandoc and the PDF (with Chromium through playwright, if available). Equations are now rendered (the previous PDF showed their LaTeX source), and the section on what is new is reorganized by version.
