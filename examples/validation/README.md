# Validation with CellPhy

Checks that a maximum-likelihood tree inferred by [CellPhy](https://github.com/amkozlov/cellphy) from CellCoal's simulated VCF recovers the true genealogy.

```
./run_cellphy_validation.sh [clean|noisy] [number of replicates]
```

- Needs a compiled CellCoal (`make` in the repository root) and a copy of the CellPhy repository (it includes the `raxml-ng-cellphy` binaries for Linux and macOS): `git clone https://github.com/amkozlov/cellphy ~/cellphy`.
- Environment variables: `CELLCOAL` (binary, default the one in `../../bin`), `CELLPHY` (default `~/cellphy`), `THREADS` (default 2), `WORKDIR`.
- It simulates 20 cells and 3000 sites per replicate, infers the tree with `GT10+FO`, removes the outgroup and compares the unrooted topologies with `rf_distance.py` (Robinson-Foulds distance; no dependencies besides Python 3).

Scenarios: `clean` = `-C30 -V5 -E0.001`; `noisy` = `-C15 -V5 -E0.01 -D0.1 -A0.05 0.01 0`.

Results obtained with CellCoal 2.0.1 (normalized RF, 0 = identical topology; a binary tree with 20 cells has 17 internal splits):

| Scenario | Replicates | SNVs per replicate | Mean normalized RF |
|---|---|---|---|
| clean | 5 | 315–594 | 0.047 |
| noisy | 10 | about 1,830–1,985 | 0.171 |

Topology only; branch lengths are not compared. Small sample sizes: use more replicates for a serious benchmark.
