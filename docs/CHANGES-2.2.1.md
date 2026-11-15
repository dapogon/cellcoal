# CellCoal 2.2.1 — what changed

**Without the new options the output of all existing runs is bit-identical to 2.2.0.**

## BCF output (`--bcf`)

The genotype likelihoods are written as a BCF 2.2 file (`vcf.bcf`, `vcf_dir/vcf.NNNN.bcf` with `-v`) instead of the VCF, about the size of `vcf.gz`. The program also writes the CSI index (`vcf.bcf.csi`), so `bcftools view -r` works at once. Verified with an independent decoder (`tests/bcf_check.py`) and with bcftools in CI.

## Streamed reads with `-3 -4 -5 -Y`

Large data sets can now be streamed with the full ML genotypes (`-3`) and haplotypes (`-4`), the ancestors (`-5`) and the true haplotypes (`-Y`): the ML alleles go to a temporary file and the files are written row by row. Identical files to the stored mode, much less memory (200 cells x 60,000 sites with `-3 -4 -Y`: 106 MB instead of 825 MB).

## Faster

A per-thread exact cache of genotype likelihoods makes the default mode about 30% faster at low coverage without `-A` (bit-identical output; `CELLCOAL_NO_GL_CACHE=1` disables it). `--fast` is now 25-30% faster than the default (vectorizable batch logarithms and powers, shorter polynomials).
