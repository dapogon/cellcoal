# CellCoal 2.8.0 — what changed

**Output of all existing runs is bit-identical to 2.7.0.**

## Coverage dispersion between sites (-O)

`-O alpha` multiplies the expected depth of every cell at a site by a gamma variate (mean 1, shape alpha) drawn once per site, so all the cells share the site's depth bias (mappability, capture). It combines with `-V` (dispersion within a site). With `-O3` the squared coefficient of variation of the bulk depth across sites is 0.33 (1/alpha), as expected. Stored and streamed reads give the same results.

## ADO at invariant sites (sparse mode)

The dropout of alleles at sites that do not vary is now counted when the site information is printed (`-y3`); the fraction of dropped alleles matches the requested rate (0.166 vs 0.163 expected for `-D0.3`).

## Packaging

- `make install [PREFIX=...]` installs `cellcoal`, `cellcoal-bulk-vaf`, the manual and the examples.
- The Release workflow attaches tarballs for Linux x86_64 (static, OpenMP) and macOS arm64 (libomp bundled), built by `.github/workflows/binaries.yml`.
- After each release, the workflow writes `Formula/cellcoal.rb` (Homebrew tap: `brew tap dapogon/cellcoal https://github.com/dapogon/cellcoal`) and `conda/meta.yaml` from the templates in `packaging/`.
