# Unreleased (after 2.6.0)
- A mistyped long option now says which one it is and suggests the closest (`Unknown option: --seg-lenn (did you mean --seg-len?)`) instead of `Incorrect parameter: -`.
- Fix: `--check` together with `--sweep` ran the grid and failed; it now checks the options and does not run it.
- `tests/manual_check.py` (also in `tests/run_tests.sh`) runs the commands of the manual and of the README and checks that every long option of `--help` is described in the manual.
- `cnv_segments` is written 90 times faster with sparse genomes (5,000 cells x 400,000 sites: 12.4 s to 0.14 s): only the variable sites and the ends of the extra-copy segments are visited, not every site of every cell. The file is identical. Losses and gains with 10,000 cells x 1,000,000 sites take about 0.7 s and 260 MB, as without them.
