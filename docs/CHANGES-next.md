# Unreleased (after 2.6.2)
- `--germline-pop FILE`: germline variation of a population. A VCF of sites with the allele frequency (`INFO/AF`, as in gnomAD or 1000 Genomes) gives the genotype of one individual at every biallelic SNV of the chromosome, drawn with Hardy-Weinberg proportions (heterozygous alleles on a random haplotype; a new individual in every replicate; sites independent, no linkage disequilibrium). Used with `--fasta`; not with `-c` or `--germline-vcf`. The sites are flagged `GERMLINE`. Without the option nothing changes.

