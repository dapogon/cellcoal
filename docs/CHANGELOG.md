# CellCoal changelog

Every release has its own file with all the details (`docs/CHANGES-X.Y.Z.md`); this page lists them, newest first. The README has a short summary of each series.

## 2.7.2
Fixes: the pileup (`-9`) wrote only the reads of bases present in some genotype of the site, so the bases column could be shorter than the depth; now every read is written and the counts match the `RC` field of the VCF. Windows builds now find zlib, so BAM, BCF and FASTQ.gz are compressed and `--compress-level` has an effect. `-C` above 1000000 is an error. Tests and CI: portable on macOS and Windows, sanitizers pass.

## 2.7.1
Bug-fix and documentation release: `-d` (single-nucleotide deletions of isolated sites) is refused with a chromosome in blocks, a reference or reads, and the error points to `--indel-rate`; manual with section 5.4.9 (how ADO, gain and coverage fit together, MDA and PTA), new figures (main flow, tumor regions, blocks and regional rates, indels), every figure numbered and with a caption, page numbers in the PDF; MDA and PTA examples.

## 2.7.0
Tumor regions (`--demes`: a structured coalescent with colonizations, without migration), small indels in the reads (`--indel-rate`, with CIGAR, `indels_truth.vcf` and left-alignment), sequencing indel errors in homopolymers (`--homopolymer-error`), population germline variation (`--germline-pop`) and a cleaner code (no `-Wshadow` warnings). Details in [CHANGES-2.7.0.md](CHANGES-2.7.0.md).

## 2.6.2
`--damage` (library damage with orientation bias) and a reads mode that works in tiles of the genome: the memory does not depend on the length of the genome (4 Mb, 4 cells, 30X: 0.96 GB instead of 2.4 GB; a chromosome in blocks, 87 MB). The reads of a given seed change (same statistics).

## 2.6.1
Doublets (`-B`) in reads mode, `--region-rates` (mutation rate that changes along the chromosome, needs `--blocks` or `--fasta`), `--seq-context` (sequencing errors that depend on the sequence), faster `cnv_segments` with sparse genomes, suggestions for mistyped options, `tests/manual_check.py`; the data are identical to 2.6.0.

## 2.6.0
Reads mode parallel by cells (the reads of a given seed change; same statistics), `--compress-level`, valid zero and negative values accepted by the long options, `--check` shows copy number, `--help` fixes. See [CHANGES-2.6.0.md](CHANGES-2.6.0.md).

## 2.5.3
Documentation: the manual is reorganized from top to bottom (copy number with the mutation models, starting from a real chromosome, one chapter per topic). Output identical to 2.5.2. See [CHANGES-2.5.3.md](CHANGES-2.5.3.md).

## 2.5.2
Documentation: a rewritten README and an updated manual. Output identical to 2.5.1. See [CHANGES-2.5.2.md](CHANGES-2.5.2.md).

## 2.5.1
Fix: the Windows build works again (2.4.0 to 2.5.0 have no Windows zip). Output identical to 2.5.0. See [CHANGES-2.5.1.md](CHANGES-2.5.1.md).

## 2.5.0
Copy number: segment losses and gains interleaved with the SNVs, deletions and cnLOH, branch by branch, over copies with known phase; `--phase`, haplotypes per copy and `cnv_segments`; every substitution model. See [CHANGES-2.5.0.md](CHANGES-2.5.0.md).

## 2.4.2
One command-line parser, grouped options (`--ado`, `--depth`, `--output`, ...), a reorganized `--help`, and about 2,300 lines less code. Output identical to 2.4.1. See [CHANGES-2.4.2.md](CHANGES-2.4.2.md).

## 2.4.1
`--amp-cell-var` and chromosome blocks in reads mode, variable base qualities and strand bias in the reads. See [CHANGES-2.4.1.md](CHANGES-2.4.1.md).

## 2.4.0
Allelic dropout in contiguous segments (`--ado-block`), signature files (`--signature-file`), signature-shaped amplification errors (`--amp-signature`), autocorrelated coverage (`--gain-sd`), the read-level simulator (`--bam`, `--fastq`, `truth.vcf`) and a chromosome in blocks inside the program (`--blocks`, `--fasta`, `--region`, `--germline-vcf`) and explicit modes (genotypes, read counts, reads) with option checks. See [CHANGES-2.4.0.md](CHANGES-2.4.0.md).

## 2.3.0
Tumor-only runs (`--no-outgroup`), calling filters that keep singletons (`--min-alt-reads`, `--min-alt-vaf`), `--coverage` (depth of the non-variable sites), `tools/chrom.py` for long chromosomes, and fixes for user trees (`-T`). See [CHANGES-2.3.0.md](CHANGES-2.3.0.md).

## 2.2.7

`--amp-cell-var`: the amplification error varies among cells. Fix: `run.json` now records the long options, so `--config Results/run.json` repeats any run. The README, the manual and `CITATION.cff` now state that CellCoal 2 was developed together with Claude (Anthropic). Without the new option the output is identical to 2.2.6. See [docs/CHANGES-2.2.7.md](CHANGES-2.2.7.md).

## 2.2.6

`cellcoal --version`; `--check` validates the options and shows what a run would do without simulating; a binary for Intel Macs; `LICENSE` (GPL-2), `CITATION.cff` and issue templates; the template `parameters` file has 100 cells and 1000 sites; the macOS package can be unquarantined with `xattr -dr`. Output identical to 2.2.5. See [docs/CHANGES-2.2.6.md](CHANGES-2.2.6.md).

## 2.2.5

A quickstart chapter in the manual (installation, a first simulation, worked examples, a `--sweep` tutorial); `--summary` reports the error by genotype type and doublets, and the simulated allelic dropout; the release has a Windows zip, and CellCoal builds with OpenMP on Windows; the plotting tools work with the files of a single replicate. The manual and README now say clearly that CellCoal simulates read counts, not individual reads. Output identical to 2.2.4.

## 2.2.4

`--summary` also reports metrics against the true genotypes (called fraction, genotype error, SNV precision and recall); `--sweep` runs combinations in parallel (`--sweep-jobs N`), resumes (`--sweep-resume`) and can be plotted (`cellcoal_plots.py sweep`). Runs without read counts are about 35% faster, and the program builds and passes the whole test suite on Windows (MSYS2). Output identical to 2.2.3.

## 2.2.3

Three options for experiments with many conditions: `--cov-doublet-corr RHO` and `--ado-doublet-corr RHO` correlate being a doublet (`-B`) with the coverage (`-J`) and the ADO rate (`-P`) of the cells (the per-cell values go to `cell_params`); `--summary` writes `summary.tsv` with one row of metrics per replicate; and `--sweep "C=2,5,10 E=0.001,0.01"` runs a grid of parameter values with a common seed and collects the mean and standard error of the metrics in `sweep.tsv`. Fixed a crash with `--ado-doublet-corr` alone. Internally `main` is split into functions and the CI now checks speed against the latest release. Output of existing runs is bit-identical to 2.2.2. See section 6.10i of the manual and [docs/CHANGES-2.2.3.md](CHANGES-2.2.3.md).

## 2.2.2

Fixes found by the new CI checks: a 2-byte heap overflow in the BCF writer (2.2.0 and 2.2.1; files were valid but memory was corrupted at the end of a record), hangs with enormous gamma shapes (`-a1e308`, `-J`, `-O`, `-V`, `-i`) and two integer overflows with absurd values. The CI now also builds with clang and without OpenMP, runs the tests under AddressSanitizer/UBSan, fuzzes the options and `--config`, reports coverage, and fails when `make test` fails (before, `make test | tee` hid the exit status). `tests/theory.py` validates the model against theory (see section 6.10h of the manual). See [docs/CHANGES-2.2.2.md](CHANGES-2.2.2.md).

## 2.2.1

`--bcf` writes the likelihoods as a BCF file (`vcf.bcf`) instead of the VCF; it reads with bcftools/htslib, and has about the size of `vcf.gz` (4-5 times smaller than the VCF). The program also writes its CSI index (`vcf.bcf.csi`), so region queries work at once. **Streamed read counts now work with `-3 -4 -5 -Y`** (full ML genotypes and haplotypes, ancestors, true haplotypes): 200 cells x 60,000 sites with `-3 -4 -Y` use 106 MB instead of 825 MB, with identical files. A per-thread **cache of genotype likelihoods** (cells and sites with the same read counts and amplification errors share them) makes the read counts about 30% faster at low coverage without `-A`, with bit-identical output. `--fast` is now about 25-30% faster than the default: the genotype likelihoods are computed in `double` with vectorizable batch logarithms and powers (`tests/fast_stats.py` compares them with the exact code). Without `--fast` nothing changes. See [docs/CHANGES-2.2.1.md](CHANGES-2.2.1.md).

## 2.2.0

`--fast` (or `fast: true` in a `--config` file) simulates the read counts about 10-15% faster with the same distributions but different random numbers (so a seed gives different read counts); the default output is unchanged. `tests/fast_stats.py` checks that depth, read counts and likelihoods have the same distributions. `--cov-ado-corr RHO` correlates the coverage (`-J`) and the ADO rate (`-P`) of the cells. With `-P`, the ADO rate drawn for each cell is written to the new file `ado_cells`. `--config` errors now give the line and column and suggest the correct key for a misspelled one. Also in this release: `plots all` (every plot plus a `report.html`), CI and packaging tests (Homebrew, conda). Without the new options the output is bit-identical to 2.1.9. See [docs/CHANGES-2.2.0.md](CHANGES-2.2.0.md).

## 2.1.9

The source is split into modules (`src/tree.c`, `mutation.c`, `observation.c`, `reads.c`, `reads_output.c`, `output.c`, `rng.c`, `options.c`, `config.c`, `globals.c`; the main program is `cellcoal.c`); output and speed are unchanged. **Configuration files**: `bin/cellcoal-2.1.9 --config examples/config.example.yaml` reads the options from a JSON or YAML file (option letters or long names such as `cells`, `coverage`), and `Results/run.json` now has an `options` object that can be passed to `--config` to repeat a run. **Summary plots**: `python3 tools/cellcoal_plots.py depth|missing|heatmap ...` (and `all Results` for every plot plus a `report.html`, from the next release) (depth per cell, site and genotype; genotypes without read counts; genotype heatmap). See [docs/CHANGES-2.1.9.md](CHANGES-2.1.9.md).

## 2.1.8

New option `-J alpha`: coverage varies between cells (shared by all the sites of a cell; e.g. `-J4`), complementing `-O` (between sites) and `-V` (within a genotype). The read counts of doublets now come from the per-site random stream in stored and streamed mode, so both give identical data; data sets with doublets differ from 2.1.7, the rest are identical. See [docs/CHANGES-2.1.8.md](CHANGES-2.1.8.md).

## 2.1.7

New option `-O alpha`: coverage varies between sites, shared by all the cells (e.g. `-O3`). ADO of invariant sites is counted in sparse mode. `make install`, prebuilt Linux/macOS tarballs in each release, a Homebrew formula (`brew tap dapogon/cellcoal https://github.com/dapogon/cellcoal && brew install cellcoal`) and a conda recipe (`conda/`). Existing output is identical to 2.1.6. See [docs/CHANGES-2.1.7.md](CHANGES-2.1.7.md).

## 2.1.6

Streamed read counts now also write the CATG file (`-8`) and the maximum-likelihood SNV genotypes and haplotypes (`-1`, `-2`) obtained from read counts, so large data sets with these outputs run in small memory (1,000 cells x 30,000 sites: 191 MB instead of 1.7 GB). The files are identical to the stored mode. See [docs/CHANGES-2.1.6.md](CHANGES-2.1.6.md).

## 2.1.5

`Results/run.json` records how every run was made (version, command line, seed, threads, environment, modes, times); the end-of-run time is now elapsed time; the Makefile rebuilds when a header changes; and the manual has a single Markdown source (`docs/manual/cellcoal-manual.md`, built with `python3 tools/build_manual.py`). Existing output is identical to 2.1.4. See [docs/CHANGES-2.1.5.md](CHANGES-2.1.5.md).

## 2.1.4

Faster bulk runs (`-0 -Z` skips the genotype likelihoods), `-z` (gzip-compressed VCF), `-M` (minimum pooled VAF of the sites in the VCF) and `--depth`/`--purity` in `tools/bulk_vaf.py`. Existing output is identical to 2.1.3. See [docs/CHANGES-2.1.4.md](CHANGES-2.1.4.md).

## 2.1.3

New option `-Z`: do not write the VCF (needs `-8`, `-9` or `-0`); with `-0` it gives the bulk VAF of thousands of cells in seconds. `tools/bulk_vaf.py` gets `--min-vaf` and `--logy`. Existing output is identical to 2.1.2. See [docs/CHANGES-2.1.3.md](CHANGES-2.1.3.md).

## 2.1.2

New option `-0`: pooled (bulk) read counts and VAF per site, without the per-cell VCF. Existing output is identical to 2.1.1. See [docs/CHANGES-2.1.2.md](CHANGES-2.1.2.md).

## 2.1.1

Streamed read counts: for very large data sets, read counts, genotype likelihoods and the VCF/PILEUP files are produced in blocks of sites with small memory. Existing output is identical to 2.1.0. See [docs/CHANGES-2.1.1.md](CHANGES-2.1.1.md).

## 2.1.0

Sparse genome storage and a lazy ADO / genotyping-error layer: 10,000 cells x 1,000,000 sites in about 7 s and 360 MB (genotype outputs; read counts/VCF not yet). Dense output is identical to 2.0.3. See [docs/CHANGES-2.1.0.md](CHANGES-2.1.0.md).

## 2.0.3

Performance: faster allele counting and genomes stored with one byte per site (about 2.5x faster and 4x less memory for large genotype-only runs). Output is identical to 2.0.2. See [docs/CHANGES-2.0.3.md](CHANGES-2.0.3.md).

## 2.0.2

Fix: germline SNPs (`-c`) are now inherited by all cells (before, they were written only to the root genome). See [docs/CHANGES-2.0.2.md](CHANGES-2.0.2.md).

## 2.0.1

Runs that do not need per-cell site data (for example, only SNV genotypes) use about 66% less memory, with identical output. See [docs/CHANGES-2.0.1.md](CHANGES-2.0.1.md).

## 2.0.0

Bug fixes (genotyping error, allelic dropout variation, memory safety), 3–4× faster with about 7× less memory, and parallel read-count generation (`OMP_NUM_THREADS`). Output for a given seed differs from versions ≤ 1.3.4 (but is identical for any number of threads). VCF likelihoods that are exactly zero are written as `.` instead of `-inf`, and `PL`/`PLN` are now −10·log10 (non-negative), as in the VCF specification. See [docs/CHANGES-2.0.0.md](CHANGES-2.0.0.md).
