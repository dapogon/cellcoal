# CellCoal 2.0.2 — what changed

Patch release on top of 2.0.1.

## Bug fix: germline SNPs (-c)

In 2.0.0 and 2.0.1, the germline SNPs requested with `-c` were written only to the genome of the root (the ancestor of the healthy cell and the tumor), after the genomes of all the other nodes had already been copied from it. As a result, the sampled cells and the healthy cell never inherited the germline variants. (Providing the germline genome directly with `-U` was not affected.) Thanks to the user who reported it.

The SNPs are now passed down to every node of the genealogy, for both the binary and the DNA alphabets. With `-c` > 0, simulated data change; with `-c` = 0 (the default) output is identical to 2.0.1 for the same seed.

## Validation

With `-c0.1 -u0` (no somatic mutations), every cell and the healthy cell now carry the same heterozygous germline sites (maternal and paternal genomes differ at the SNP sites); before the fix all genomes were identical.
