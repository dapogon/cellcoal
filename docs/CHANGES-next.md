# Unreleased (after 2.6.0)
- A mistyped long option now says which one it is and suggests the closest (`Unknown option: --seg-lenn (did you mean --seg-len?)`) instead of `Incorrect parameter: -`.
- Fix: `--check` together with `--sweep` ran the grid and failed; it now checks the options and does not run it.
- `tests/manual_check.py` (also in `tests/run_tests.sh`) runs the commands of the manual and of the README and checks that every long option of `--help` is described in the manual.
