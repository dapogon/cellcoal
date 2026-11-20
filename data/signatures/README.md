# Mutational signature profiles (96 trinucleotide contexts, COSMIC format)

Used with `--signature-file FILE` (the columns of the file are the signatures 1, 2, ... of `-S`).

- `COSMIC_v3.6_SBS1_SBS5_GRCh37.txt`: SBS1 (column 1) and SBS5 (column 2) from COSMIC v3.6 (GRCh37), taken from the
  reference signatures of SigProfilerAssignment (`COSMIC_v3.6_SBS_GRCh37.txt`).
- `PTA_artifact_SCAN2.txt`: the SNV artifact signature of primary template-directed amplification (PTA) used by SCAN2
  (Luquette et al. 2022, Nat Genet), from `data/snv.artifact.signature.v3.rda` of the `r-scan2` package
  (github.com/parklab/r-scan2). 67.7% of its mass is C>T; it resembles SBS30 / SBS102.

Each column sums to 1. Rows are `5'[ref>alt]3'` with the pyrimidine as reference, as in COSMIC.
