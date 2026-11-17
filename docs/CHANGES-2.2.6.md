# CellCoal 2.2.6 — what changed

**The simulated data of all existing runs are bit-identical to 2.2.5** (only two words of the log changed).

- Manual: how to open the macOS program (Gatekeeper warning for unsigned programs).
- The template `parameters` file has 100 cells and 1000 sites (it had 10 and 10).
- The macOS package has `libomp.dylib` writable, so `xattr -dr com.apple.quarantine .` works on the whole folder.
- `CITATION.cff` (the *Cite this repository* button of GitHub) and issue templates (bug report, feature request) that ask for the version, system, command and `run.json`.
- `cellcoal --version` prints the version and `--help` the usage (both exit at once).
- `--check`: reads and validates all the options (command line, `--config` or `parameters`) and shows what a run would do (replicates, cells, sites, expected mutations and tree depth, outputs, results folder) without simulating anything. Invalid options give a non-zero exit code.
- Log text: "No read counts were produced" and "numAltAlleles (in read counts)" (CellCoal does not simulate individual reads).
- `LICENSE` file added.
- Binary for Intel Macs (`macos-x86_64`) in the releases.
