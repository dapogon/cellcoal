# Unreleased (after 2.2.4)

- `--summary` adds `het_error`, `hom_error`, `het_to_hom` (the footprint of allelic dropout), `doublet_error`, `singlet_error` and `ado_simulated` (the allelic dropout simulated in the genotypes, to compare with `het_to_hom`).
- Windows: OpenMP is enabled automatically, and the release workflow builds a Windows zip (program plus the compiler DLLs).
- Windows: PowerShell splits arguments such as `-E0.001` into `-E0` and `.001`, so runs launched from PowerShell stopped with a parameter error; the program now joins them again (test 2g21). Error messages of failed allocations printed `size_t` values with `%ld` (wrong on Windows). The Windows CI job now fails on any compiler warning and runs the plot tests, and the Windows zip is tested from plain PowerShell without MSYS2.
- Speed: the VCF formatting and the read simulation were profiled again; no change was worth its risk, the output is identical.
- Manual: new chapter 3, a quickstart with installation (prebuilt programs, Homebrew, conda, source, Windows), a first simulation, worked examples and a tutorial on `--sweep`; every command in it is run.
- Fixes in the plotting tools: `cellcoal_plots.py heatmap` and `all` failed with the SNV genotype file of a single replicate, `bulk_vaf.py --vcf` failed with the VCF of a single replicate, and the legend of the heatmap overlapped the axis label. New tests (2g20).
- Terminology: CellCoal simulates read counts (the number of reads of each nucleotide in each cell and site), not individual reads (no sequences, FASTQ or BAM). The manual, the README, the example parameter files, the help (`-C`), the error messages and the tools now say so.
