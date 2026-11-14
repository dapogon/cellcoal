# CellCoal 2.0.1 — what changed

Patch release on top of 2.0.0. Output is identical to 2.0.0 for the same seed and options.

## Memory

The per-cell, per-site blocks (`cell[i].site`, 48 bytes per cell and site) were always allocated and filled, even when no option used them. They are now allocated only when needed: read counts (NGS), SNV/true/full haplotypes, full genotypes, CATG, PILEUP or ancestral cells. SNV genotypes without reads are printed from the genotype matrix and do not need them.

| Run (SNV genotypes, genotyping error and ADO, 1 thread) | 2.0.0 | 2.0.1 |
|---|---|---|
| 500 cells × 20,000 sites: peak memory | 695 MB | 236 MB |
| 500 cells × 20,000 sites: time | 3.3 s | 2.7 s |

## Validation

- 11 option combinations (with and without reads, with haplotypes, full genotypes, ancestors, true haplotypes, CATG, PILEUP, doublets, ADO, genotyping error, rate switches): identical output and exit codes in 2.0.0 and 2.0.1.
- AddressSanitizer and UBSan: no errors in these runs and in the example parameter files. Two examples (`parameters.JC.usertree`, `parameters.signature.noreads`) still report one-off start-up allocations that are not freed at exit (also in 2.0.0); they are harmless and do not affect the results.
