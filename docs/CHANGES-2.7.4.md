# CellCoal 2.7.4 — what changed

- Targeted sequencing in the reads mode: `--targets FILE` (BED of target regions; scWES, gene panels or a few chromosome fragments; it needs a real chromosome, `--fasta`) and `--off-target F` (depth off target relative to the targets, 0.05 by default). `-C` is the depth of the targets. Manual Fig. 17, tests 2g66.
- Fix: `--off-target` accepts 0.
