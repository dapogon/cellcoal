# CellCoal 2.1.8 — what changed

**Output is bit-identical to 2.1.7 except for data sets with doublets (`-B`).**

## Coverage between cells (-J)

`-J alpha` multiplies the expected depth of every cell, at all its sites, by a gamma variate (mean 1, shape alpha) drawn once per cell and replicate. The expected reads of an allele are `coverage x site factor (-O) x cell factor (-J)`, and `-V` adds independent Negative Binomial noise on top. With `-J4` the squared coefficient of variation of the mean depth between cells is about 0.25. Stored and streamed reads give the same results.

## Doublets use the stream of each site

The reads of the donor of a doublet are now simulated site by site from the same random stream as the cell reads, in stored and in streamed mode. Both modes therefore give identical files (VCF, PILEUP, CATG, genotypes, bulk) with doublets, for the same seed and any number of threads. The choice of doublet cells and donors now happens before the per-site streams start, so data sets with doublets are different from those of 2.1.7 (same distribution). The unused sequential doublet code was removed.

## Tests

Stored = streamed with doublets and with `-J`; between-cell depth dispersion with `-J4`; the doublet test uses 300 cells for a tighter check of the doublet fraction.
