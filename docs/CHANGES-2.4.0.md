# CellCoal 2.4.0 — what changed

**Output of all existing runs is bit-identical to 2.3.0.**

## New option -Z: no VCF

With simulated reads the VCF is always written, and it is by far the largest output (about 360 bytes per cell and site). `-Z` skips it when another read output is requested (`-8` CATG, `-9` PILEUP or `-0` bulk); `-Z` alone is an error. Together with `-0` it makes the VAF distribution of thousands of cells cheap: 5,000 cells x 20,000 sites with reads takes 19 s, 244 MB and writes a 1 MB file (`-n1 -s5000 -l20000 -u2e-7 -e10000 -g1e-4 -b1 -C2 -V1 -E0.001 -0 -Z`).

## tools/bulk_vaf.py

New `--min-vaf` (default 0.01, hides the sequencing-error floor of non-mutated sites, which otherwise dominates the plot) and `--logy`.
