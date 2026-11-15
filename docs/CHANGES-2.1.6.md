# CellCoal 2.1.6 — what changed

**Output of all existing runs is bit-identical to 2.1.5.**

## CATG and ML genotypes with streamed reads

When reads are streamed (very large data sets, or `CELLCOAL_STREAM=1`), these outputs are now supported:

- `-8` CATG: the site lines are written block by block to a temporary file and copied behind the header once the number of SNVs is known.
- `-1` / `-2` SNV genotypes and haplotypes from reads: the maximum-likelihood alleles are kept in a compact store (one byte per cell and SNV).

The files equal those of the stored mode for the same seed and any number of threads (tested with 1 and 3 threads); with doublets the random numbers differ, as before. Example: 1,000 cells x 30,000 sites with `-1 -8`: 191 MB instead of 1.7 GB.

Still unavailable while streaming: `-3`, `-4`, `-5`, `-Y` (forcing streaming with them is an error).

## Tests

`make test` checks stored vs streamed CATG, SNV genotypes and haplotypes with 1 and 3 threads.
