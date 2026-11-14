# CellCoal 2.1.0 — what changed

Scaling release: 10,000 cells x 1,000,000 sites runs in seconds on a laptop (genotype outputs). **Output in dense mode is bit-identical to 2.0.3.**

## Sparse genome storage

When the dense genomes would need more than 2048 MB (`CELLCOAL_DENSE_MB` changes the limit), genomes are stored sparsely: only the sites that vary keep a per-node column; all other sites are read from the healthy root. `CELLCOAL_SPARSE=1` forces it, `CELLCOAL_SPARSE=0` disables it. Without errors, sparse and dense outputs are identical for the same seed.

## Observation layer for ADO and genotyping errors

In sparse mode allelic dropout is decided lazily by hashing (cell, site), and genotyping errors are drawn as sparse events (geometric skipping), so they are not stored in the genomes. The distributions match the dense mode; the random stream differs. ADO at invariant sites is not counted.

## Faster output

SNV genotypes are printed in blocks of cells; dense ADO no longer allocates a cells x sites rate matrix.

## Limits

- Reads, VCF, CATG and PILEUP are not supported in sparse mode (automatic mode keeps dense storage for them; `CELLCOAL_SPARSE=1` with `-C` is an error).
- Full-sequence outputs (-4, -5, -Y) stop with an error on data sets above the dense limit.

Reference: `-n1 -s10000 -l1000000 -u1e-8 -e10000 -g0 -b1 -D0.1 -G1e-6 1e-13 -1`: about 7 s, 360 MB.
