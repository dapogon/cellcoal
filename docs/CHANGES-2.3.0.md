# CellCoal 2.3.0 — what changed

**Output of all existing runs is bit-identical to 2.2.0.**

## New option -0: pooled (bulk) read counts and VAF

`-0` writes the file `bulk` (or `bulk_dir/bulk.NNNN` with `-v`), one line per site, with the reads of the four nucleotides pooled over the tumour cells (as in a bulk sample): site, reference, readsA, readsC, readsG, readsT, depth, alt_reads, bulk_vaf, true_alt_alleles, true_vaf. Only sites with alternate reads or a true alternate allele are listed. It needs simulated reads (`-C`), works with and without streamed reads (same output) and its size depends on the number of sites only, not on the number of cells, so the VAF of 5,000 or 10,000 cells can be studied without writing the huge VCF.

`python3 tools/bulk_vaf.py --bulk Results/bulk --min-reads 10 --out vaf.png` plots the true and the read-based VAF distributions.
