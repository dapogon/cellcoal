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
- AddressSanitizer, UBSan and LeakSanitizer: no reports in these runs and in the five example parameter files.
- Example parameter files with a fixed seed: all result files identical to the previous build.

## Memory freed at exit

Two start-up allocations were never freed (harmless, but reported by LeakSanitizer): the node and cell names read from a user tree (`-T`), and the genetic signature tables (`-t` signatures). They are now freed. The two tree-node names that were overwritten with string literals (`healthyroot`, `healthy`) are now copied into their buffers.
