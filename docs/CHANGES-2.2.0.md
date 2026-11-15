# CellCoal 2.2.0 — what changed

Reads (VCF and PILEUP) for data sets whose per-cell, per-site data would not fit in memory. **Output of all existing runs is bit-identical to 2.1.0.**

## Streamed reads

When `(cells+1) x sites x 48 bytes` exceeds 2048 MB (`CELLCOAL_STREAM_MB` changes the limit), the reads are simulated and printed in blocks of sites and then discarded, so memory stays small (about 220 MB for 500 cells x 100,000 sites instead of 2.7 GB, and 1.5x faster). `CELLCOAL_STREAM=1` forces it, `CELLCOAL_STREAM=0` disables it. Streamed reads imply sparse genome storage. Reads, genotype likelihoods and VCF/PILEUP lines are the same as without streaming for the same seed (this is tested), for any number of threads. Doublets are simulated from the stream of each site, so with doublets (`-B`) the random numbers differ from the non-streamed mode (same distribution).

## Limits

- Only the VCF and PILEUP outputs are supported while streaming (automatic mode keeps the old behaviour when other outputs such as `-1 -2 -3 -4 -5 -Y`, CATG or ancestors are requested; forcing `CELLCOAL_STREAM=1` with them is an error).
- Sparse genome storage now also works with reads when forced (`CELLCOAL_SPARSE=1`).
- The cost is still proportional to cells x sites (about 4 microseconds per cell-site and core), and the VCF contains every apparent SNV, which with sequencing error is almost every site: 500 cells x 100,000 sites wrote an 18 GB VCF. For 10,000 cells x 1,000,000 sites expect many hours of CPU and terabytes of output.
