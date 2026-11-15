## CellCoal 2
CellCoal: coalescent simulation of single-cell NGS genotypes 

&copy;2020 David Posada, University of Vigo, Spain <http://darwin.uvigo.es>

## About CellCoal
CellCoal simulates the somatic evolution of single-cells. CellCoal generates a coalescent genealogy for a sample of somatic cells –no recombination– obtained from a growing population, together with a another cell as outgroup, introduces mutations along this genealogy, and produces single-cell diploid genotypes (single-nuclotide variants or SNVs). CellCoal implements multiple mutations models (0/1, DNA, infinite and finite site models, deletion, copy-neutral LOH, 30 cancer signatures) and is able to generate read counts and genotype likelihoods considering allelic dropout, sequencing and amplification error, plus doublet cells. Rate variation among lineages can be simulated with a gamma distribution or with a few *rate switches* (`-N`): each switch picks a branch at random, with probability proportional to its length, and multiplies the length of that branch and of all its descendants by a given factor, creating a clade that evolves faster or slower (for example, `-N2 4 0.2`).

## Getting started

1. **Download**: get the program from the GitHub repository at <https://github.com/dapogon/cellcoal>. Under the Code tab you will see a section call *release*. Click on it download the source code in .zip (or tar.gz) format. Then unzip the *cellcoal-x.y.z.zip* file (for example, cellcoal-2.1.9.zip). You should see now a folder called *cellcoal-x.y.z*. Move into the folder to compile and run the program.

2. **Compile**: type `make`. The program should compile without much problem in Linux/MacOSX. A *Makefile* is provided. The executable file will be located in the *bin* folder (e.g., `bin/cellcoal-2.1.9`). OpenMP (multithreading) is detected automatically: `-fopenmp` if the compiler supports it, otherwise Homebrew's `libomp` on macOS (`brew install libomp`), otherwise the program is built single-threaded with identical results. Use `OMP_NUM_THREADS=n` to choose the number of threads, or `make OMPFLAG=` to build without OpenMP.

3. **Run**: `bin/cellcoal-x.y.z`. In this case, the program arguments will be read from a file called *parameters*. Alternatively, arguments can be entered directly in the command line, or from a specific parameter file with a different name.

4. **Inspect results**: Once the run is finished, different output files will be written by default inside a folder called *Results* (or the folder given with `-o`).

## Bulk-like VAF distribution

Easiest: option `-0` makes CellCoal write the reads pooled over the tumour cells, one line per site (VAF from the simulated reads and the true VAF; size independent of the number of cells), and `python3 tools/bulk_vaf.py --bulk Results/bulk --out vaf.png` plots it (add `-Z` to skip the VCF when you only need the bulk file; add `--depth 100 --purity 0.8` to mimic a real bulk experiment). The details below cover the other inputs.

`tools/bulk_vaf.py` (Python with numpy and matplotlib) draws the variant allele frequency (VAF) distribution of the sample as if it were sequenced in bulk.

- **From the simulated reads** (VCF made with `-C`, `-V`, ...): the reads (RC field) of all the tumour cells are pooled, and VAF = alternate reads / total reads at each site. The healthy cell is excluded.
  `python3 tools/bulk_vaf.py --vcf Results/vcf_dir/vcf.0001 --min-reads 10 --out vaf.png --tsv vaf.tsv`
- **From the genotypes** (`-1`): VAF is the fraction of alternate alleles among the tumour cells; a bulk experiment is mimicked by sampling reads at a given depth and tumour purity.
  `python3 tools/bulk_vaf.py Results/snv_gen --depth 100 --purity 1.0 --out vaf.png --tsv vaf.tsv`

With simulated reads, sequencing and amplification errors give almost every site a few alternate reads, so the VAF of non-mutated sites sits at the error rate near zero; for a clean plot of the true SNVs use `-E0`, or look at VAF above the error level. Pooling thousands of cells at a few reads each gives a very high pooled depth, so the VAF is very precise. The VCF grows with cells x sites (about 360 bytes per cell and site): for thousands of cells keep the number of sites low (a few thousand).

## Unreleased (next 2.1.x)

`--fast` (or `fast: true` in a `--config` file) simulates the reads about 10-15% faster with the same distributions but different random numbers (so a seed gives different reads); the default output is unchanged. `tests/fast_stats.py` checks that depth, read counts and likelihoods have the same distributions. `--cov-ado-corr RHO` correlates the coverage (`-J`) and the ADO rate (`-P`) of the cells. With `-P`, the ADO rate drawn for each cell is written to the new file `ado_cells`. `--config` errors now give the line and column and suggest the correct key for a misspelled one.

## What's new in 2.1.9

The source is split into modules (`src/tree.c`, `mutation.c`, `observation.c`, `reads.c`, `reads_output.c`, `output.c`, `rng.c`, `options.c`, `config.c`, `globals.c`; the main program is `cellcoal.c`); output and speed are unchanged. **Configuration files**: `bin/cellcoal-2.1.9 --config examples/config.example.yaml` reads the options from a JSON or YAML file (option letters or long names such as `cells`, `coverage`), and `Results/run.json` now has an `options` object that can be passed to `--config` to repeat a run. **Summary plots**: `python3 tools/cellcoal_plots.py depth|missing|heatmap ...` (and `all Results` for every plot plus a `report.html`, from the next release) (depth per cell, site and genotype; genotypes without reads; genotype heatmap). See [docs/CHANGES-2.1.9.md](docs/CHANGES-2.1.9.md).

## What's new in 2.1.8

New option `-J alpha`: coverage varies between cells (shared by all the sites of a cell; e.g. `-J4`), complementing `-O` (between sites) and `-V` (within a genotype). The reads of doublets now come from the per-site random stream in stored and streamed mode, so both give identical data; data sets with doublets differ from 2.1.7, the rest are identical. See [docs/CHANGES-2.1.8.md](docs/CHANGES-2.1.8.md).

## What's new in 2.1.7

New option `-O alpha`: coverage varies between sites, shared by all the cells (e.g. `-O3`). ADO of invariant sites is counted in sparse mode. `make install`, prebuilt Linux/macOS tarballs in each release, a Homebrew formula (`brew tap dapogon/cellcoal https://github.com/dapogon/cellcoal && brew install cellcoal`) and a conda recipe (`conda/`). Existing output is identical to 2.1.6. See [docs/CHANGES-2.1.7.md](docs/CHANGES-2.1.7.md).

## What's new in 2.1.6

Streamed reads now also write the CATG file (`-8`) and the maximum-likelihood SNV genotypes and haplotypes (`-1`, `-2`) obtained from reads, so large data sets with these outputs run in small memory (1,000 cells x 30,000 sites: 191 MB instead of 1.7 GB). The files are identical to the stored mode. See [docs/CHANGES-2.1.6.md](docs/CHANGES-2.1.6.md).

## What's new in 2.1.5

`Results/run.json` records how every run was made (version, command line, seed, threads, environment, modes, times); the end-of-run time is now elapsed time; the Makefile rebuilds when a header changes; and the manual has a single Markdown source (`docs/manual/cellcoal-manual.md`, built with `python3 tools/build_manual.py`). Existing output is identical to 2.1.4. See [docs/CHANGES-2.1.5.md](docs/CHANGES-2.1.5.md).

## What's new in 2.1.4

Faster bulk runs (`-0 -Z` skips the genotype likelihoods), `-z` (gzip-compressed VCF), `-M` (minimum pooled VAF of the sites in the VCF) and `--depth`/`--purity` in `tools/bulk_vaf.py`. Existing output is identical to 2.1.3. See [docs/CHANGES-2.1.4.md](docs/CHANGES-2.1.4.md).

## What's new in 2.1.3

New option `-Z`: do not write the VCF (needs `-8`, `-9` or `-0`); with `-0` it gives the bulk VAF of thousands of cells in seconds. `tools/bulk_vaf.py` gets `--min-vaf` and `--logy`. Existing output is identical to 2.1.2. See [docs/CHANGES-2.1.3.md](docs/CHANGES-2.1.3.md).

## What's new in 2.1.2

New option `-0`: pooled (bulk) read counts and VAF per site, without the per-cell VCF. Existing output is identical to 2.1.1. See [docs/CHANGES-2.1.2.md](docs/CHANGES-2.1.2.md).

## What's new in 2.1.1

Streamed reads: for very large data sets, reads, genotype likelihoods and the VCF/PILEUP files are produced in blocks of sites with small memory. Existing output is identical to 2.1.0. See [docs/CHANGES-2.1.1.md](docs/CHANGES-2.1.1.md).

## What's new in 2.1.0

Sparse genome storage and a lazy ADO / genotyping-error layer: 10,000 cells x 1,000,000 sites in about 7 s and 360 MB (genotype outputs; reads/VCF not yet). Dense output is identical to 2.0.3. See [docs/CHANGES-2.1.0.md](docs/CHANGES-2.1.0.md).

## What's new in 2.0.3

Performance: faster allele counting and genomes stored with one byte per site (about 2.5x faster and 4x less memory for large genotype-only runs). Output is identical to 2.0.2. See [docs/CHANGES-2.0.3.md](docs/CHANGES-2.0.3.md).

## What's new in 2.0.2

Fix: germline SNPs (`-c`) are now inherited by all cells (before, they were written only to the root genome). See [docs/CHANGES-2.0.2.md](docs/CHANGES-2.0.2.md).

## What's new in 2.0.1

Runs that do not need per-cell site data (for example, only SNV genotypes) use about 66% less memory, with identical output. See [docs/CHANGES-2.0.1.md](docs/CHANGES-2.0.1.md).

## What's new in 2.0.0

Bug fixes (genotyping error, allelic dropout variation, memory safety), 3–4× faster with about 7× less memory, and parallel read-count generation (`OMP_NUM_THREADS`). Output for a given seed differs from versions ≤ 1.3.4 (but is identical for any number of threads). VCF likelihoods that are exactly zero are written as `.` instead of `-inf`, and `PL`/`PLN` are now −10·log10 (non-negative), as in the VCF specification. See [docs/CHANGES-2.0.0.md](docs/CHANGES-2.0.0.md).

## Citation

- If you use CellCoal, please cite: Posada D. (2020) CellCoal: coalescent simulation of single-cell sequencing samples. *Molecular Biology and Evolution* 37: 1535-1542. <a href="https://academic.oup.com/mbe/article/37/5/1535/5728646" target="_blank">https://academic.oup.com/mbe/article/37/5/1535/5728646</a>.

## Tests

`make test` runs a small regression suite (`tests/run_tests.sh`): the example parameter files, identical results with 1 and 2 threads, reproducibility with the same seed, and inheritance of germline SNPs (`-c`), and that the observed rates of allelic dropout (`-D`), genotyping error (`-G`) the fixed number of mutations (`-j`), allelic imbalance (`-I`), doublets (`-B`), ADO per cell (`-P`) and the spectrum of mutational signature 1 match the expected values. It also checks that sparse genome storage gives the same output as dense storage and the same ADO and genotyping-error rates, that streamed reads give the same VCF, PILEUP and bulk file as stored reads, that `-0`, `-Z`, `-z` and `-M` behave, that the pooled (bulk) VAF of the reads agrees with the true VAF (ratio of means and correlation), that doublets multiply the pooled depth by the expected factor with stored and streamed reads, and the content of `run.json`. It does not depend on exact floating-point output, so it can be run on any platform.

## Validation
A script that simulates with CellCoal, infers the tree with CellPhy and compares it with the true tree is in [examples/validation](examples/validation).

## Manual
A detailed manual is available in [html](https://dapogon.github.io/cellcoal/) or in [pdf](docs/cellcoal.manual.v2.0.pdf).

