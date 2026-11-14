# CellCoal 2.0.3 — what changed

Performance release on top of 2.0.2. **Output is bit-identical to 2.0.2** for the same parameters and seed.

## Faster allele counting

`CountTrueVariants` and `CountAllelesInObservedGenotypes` walked the genomes site by site, jumping from one cell's genome to the next at every access, which is very unfriendly to the CPU cache. They now traverse each cell's genome contiguously (in blocks of sites).

## One byte per site

Genome matrices were stored as 4-byte integers although every value (nucleotide, deletion, ADO, missing) fits in one byte. They are now stored as `signed char`, dividing the memory used by the genomes by four.

## Effect

With 1000 cells and 50,000 sites (genotypes only): about 1.5 s to 0.6 s and 775 MB to 203 MB. Runs that simulate reads are not noticeably affected.

## Validation

The five example parameter files (fixed seed) and five further configurations (reads, signatures, deletions, cnLOH, binary alphabet, user genome) give identical output to 2.0.2. AddressSanitizer and UBSan report nothing, and `make test` passes on Linux and macOS.
