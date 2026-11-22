## CellCoal 2
CellCoal: coalescent simulation of single-cell NGS genotypes 

&copy;2020 David Posada, University of Vigo, Spain <http://darwin.uvigo.es>

## About CellCoal
CellCoal simulates the somatic evolution of single-cells. CellCoal generates a coalescent genealogy for a sample of somatic cells –no recombination– obtained from a growing population, together with a another cell as outgroup, introduces mutations along this genealogy, and produces single-cell diploid genotypes (single-nuclotide variants or SNVs). CellCoal implements multiple mutations models (0/1, DNA, infinite and finite site models, deletion, copy-neutral LOH, 30 cancer signatures) and is able to generate read counts and genotype likelihoods considering allelic dropout, sequencing and amplification error, plus doublet cells. **CellCoal does not simulate individual reads (there are no sequences, no FASTQ and no BAM files): it simulates the number of reads of each nucleotide in each cell and site (read counts).** Rate variation among lineages can be simulated with a gamma distribution or with a few *rate switches* (`-N`): each switch picks a branch at random, with probability proportional to its length, and multiplies the length of that branch and of all its descendants by a given factor, creating a clade that evolves faster or slower (for example, `-N2 4 0.2`).

## Developed with Claude

**CellCoal 2 was developed together with Claude.** The CellCoal 2 series (versions 2.0.0 onwards) was developed by Claude, an AI model made by Anthropic, working as a coding agent (Claude Code) in collaboration with David Posada. Claude took part in all the work on the code: finding and fixing bugs, making the program faster and lighter in memory, splitting the source into modules, new features (sparse and streamed genomes, BCF output, `--summary`, `--sweep`, `--config`, `--check`), the test suite, the continuous integration, the packages for Linux, macOS and Windows, the releases and the documentation, including this text. David Posada defined the goals, made the scientific and design decisions, tested the program and approved every release, and is the author of CellCoal and responsible for it. The 2020 article describes version 1, which was written before this work.

## Getting started

New to CellCoal? Chapter 3 of the manual ([docs/cellcoal.manual.v2.0.html](docs/cellcoal.manual.v2.0.html)) is a quickstart: installation on Linux, macOS and Windows, a first simulation, worked examples (read counts and VCF, bulk VAF, doublets, summaries, parameter and configuration files) and a tutorial on `--sweep`.

1. **Download**: get the program from the GitHub repository at <https://github.com/dapogon/cellcoal>. Under the Code tab you will see a section call *release*. Click on it download the source code in .zip (or tar.gz) format. Then unzip the *cellcoal-x.y.z.zip* file (for example, cellcoal-2.2.3.zip). You should see now a folder called *cellcoal-x.y.z*. Move into the folder to compile and run the program.

2. **Compile**: type `make`. The program should compile without much problem in Linux/MacOSX. A *Makefile* is provided. The executable file will be located in the *bin* folder (e.g., `bin/cellcoal-2.2.3`). OpenMP (multithreading) is detected automatically: `-fopenmp` if the compiler supports it, otherwise Homebrew's `libomp` on macOS (`brew install libomp`), otherwise the program is built single-threaded with identical results. Use `OMP_NUM_THREADS=n` to choose the number of threads, or `make OMPFLAG=` to build without OpenMP.

3. **Run**: `bin/cellcoal-x.y.z`. In this case, the program arguments will be read from a file called *parameters*. Alternatively, arguments can be entered directly in the command line, or from a specific parameter file with a different name.

4. **Inspect results**: Once the run is finished, different output files will be written by default inside a folder called *Results* (or the folder given with `-o`).

## Bulk-like VAF distribution

Easiest: option `-0` makes CellCoal write the read counts pooled over the tumour cells, one line per site (VAF from the simulated read counts and the true VAF; size independent of the number of cells), and `python3 tools/bulk_vaf.py --bulk Results/bulk --out vaf.png` plots it (add `-Z` to skip the VCF when you only need the bulk file; add `--depth 100 --purity 0.8` to mimic a real bulk experiment). The details below cover the other inputs.

`tools/bulk_vaf.py` (Python with numpy and matplotlib) draws the variant allele frequency (VAF) distribution of the sample as if it were sequenced in bulk.

- **From the simulated read counts** (VCF made with `-C`, `-V`, ...): the read counts (RC field) of all the tumour cells are pooled, and VAF = alternate reads / total reads at each site. The healthy cell is excluded.
  `python3 tools/bulk_vaf.py --vcf Results/vcf_dir/vcf.0001 --min-reads 10 --out vaf.png --tsv vaf.tsv`
- **From the genotypes** (`-1`): VAF is the fraction of alternate alleles among the tumour cells; a bulk experiment is mimicked by sampling reads at a given depth and tumour purity.
  `python3 tools/bulk_vaf.py Results/snv_gen --depth 100 --purity 1.0 --out vaf.png --tsv vaf.tsv`

With simulated read counts, sequencing and amplification errors give almost every site a few alternate reads, so the VAF of non-mutated sites sits at the error rate near zero; for a clean plot of the true SNVs use `-E0`, or look at VAF above the error level. Pooling thousands of cells at a few reads each gives a very high pooled depth, so the VAF is very precise. The VCF grows with cells x sites (about 360 bytes per cell and site): for thousands of cells keep the number of sites low (a few thousand).

## Versions

A short summary of each series; the details of every release are in [docs/CHANGELOG.md](docs/CHANGELOG.md) and in `docs/CHANGES-X.Y.Z.md`.

- **2.4: realistic amplification, Illumina reads and a chromosome in blocks.** Allelic dropout in contiguous segments, signature files, signature-shaped amplification errors, autocorrelated coverage; paired-end reads (BAM, FASTQ) with PCR duplicates, base qualities and a truth VCF; a chromosome in blocks inside the program (`--blocks`, `--fasta`, `--region`); explicit modes (genotypes, read counts, reads); grouped options (`--ado`, `--depth`, `--output`, ...) and a cleaner code base (2.4.1 and 2.4.2: about 2,300 lines less, same output).
- **2.3: tumor-only runs and long chromosomes.** `--no-outgroup`, calling filters that keep singletons (`--min-alt-reads`, `--min-alt-vaf`), the depth of the non-variable sites (`--coverage`) and `tools/chrom.py` to simulate a long chromosome by blocks.
- **2.2: speed, BCF, summaries and sweeps.** Faster read counts (`--fast`), BCF output, streaming with every output, `--summary` with metrics against the true genotypes, `--sweep` parameter grids (parallel and resumable), validation against theory and a fuzzer-hardened code, `--check`, `--amp-cell-var`, a quickstart chapter, and binaries for Linux, macOS and Windows.
- **2.1: very large data sets.** Sparse genome storage (10,000 cells x 1,000,000 sites in about 7 s), streamed read counts with small memory, bulk VAF (`-0`), coverage varying between sites and cells (`-O`, `-J`), configuration files, a modular source, `run.json` run records, Homebrew and conda packages.
- **2.0: a faster, safer CellCoal.** Bug fixes (genotyping error, allelic dropout, memory safety), 3-4 times faster with about 7 times less memory, and parallel read counts. The output for a given seed differs from versions 1.3.4 and earlier.

## Citation

- If you use CellCoal, please cite: Posada D. (2020) CellCoal: coalescent simulation of single-cell sequencing samples. *Molecular Biology and Evolution* 37: 1535-1542. <a href="https://academic.oup.com/mbe/article/37/5/1535/5728646" target="_blank">https://academic.oup.com/mbe/article/37/5/1535/5728646</a>.

## Tests

The CI also builds with clang and without OpenMP, runs the suite under AddressSanitizer/UBSan, fuzzes the options and the `--config` parser (`tests/fuzz.py`), reports code coverage, and `tests/theory.py` (part of `make test`) validates the model against theory: coalescent waiting times, growth, number and site frequency spectrum of mutations, depth, sequencing error and dropout. See section 6.10h of the manual.

`make test` runs a small regression suite (`tests/run_tests.sh`): the example parameter files, identical results with 1 and 2 threads, reproducibility with the same seed, and inheritance of germline SNPs (`-c`), and that the observed rates of allelic dropout (`-D`), genotyping error (`-G`) the fixed number of mutations (`-j`), allelic imbalance (`-I`), doublets (`-B`), ADO per cell (`-P`) and the spectrum of mutational signature 1 match the expected values. It also checks that sparse genome storage gives the same output as dense storage and the same ADO and genotyping-error rates, that streamed read counts give the same VCF, PILEUP and bulk file as stored read counts, that `-0`, `-Z`, `-z` and `-M` behave, that the pooled (bulk) VAF of the read counts agrees with the true VAF (ratio of means and correlation), that doublets multiply the pooled depth by the expected factor with stored and streamed read counts, and the content of `run.json`. It does not depend on exact floating-point output, so it can be run on any platform.

## Validation
A script that simulates with CellCoal, infers the tree with CellPhy and compares it with the true tree is in [examples/validation](examples/validation).

## Manual
A detailed manual is available in [html](https://dapogon.github.io/cellcoal/) or in [pdf](docs/cellcoal.manual.v2.0.pdf). To start, read the [quickstart](https://dapogon.github.io/cellcoal/cellcoal.manual.v2.0.html#toc-quickstart) (installation and examples).

## License

CellCoal is free software, distributed under the GNU General Public License, version 2 or (at your option) any later version; see [LICENSE](LICENSE).
