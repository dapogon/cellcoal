# CellCoal 2.5.0 — what changed

**Output of all existing runs is bit-identical to 2.4.0.**

## Faster bulk runs

When only the bulk file is requested (`-0 -Z`), the genotype likelihoods, which take about half of the time of the read simulation and are not needed, are skipped, and reads are streamed automatically: 5,000 cells x 20,000 sites went from 19 s / 244 MB to 12.8 s / 102 MB. The bulk file is identical to the one produced with the likelihoods (tested), because the likelihoods use no random numbers. The bulk file now also includes the reads of doublets.

## VCF: compression and filtering

- `-z` writes the VCF compressed with gzip (`vcf.gz`, `vcf_dir/vcf.NNNN.gz`); gzip must be installed. The content is identical to the plain VCF.
- `-M f` lists in the VCF only the apparent SNVs whose pooled VAF in reads (alternate reads / reads of all the cells) is at least `f` (for example `-M0.02`), which removes the sites that only carry sequencing errors.

## tools/bulk_vaf.py

With `--bulk`, new `--depth D` samples the pooled reads at a realistic bulk depth (Poisson, binomial sampling) and `--purity p` dilutes them with healthy cells, so the VAF distribution looks like a real bulk experiment (for example `--depth 100 --purity 0.8`).
