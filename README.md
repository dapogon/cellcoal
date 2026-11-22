# CellCoal 2

**Coalescent simulation of single-cell sequencing data: from the genealogy of a tumour to the reads of every cell.**

&copy;2020 David Posada, University of Vigo, Spain <http://darwin.uvigo.es>

## What CellCoal does

CellCoal simulates the somatic evolution of a sample of cells and the single-cell sequencing experiment that observes them. It generates a coalescent genealogy of the cells (no recombination; a growing or birth-death population, or a tree you provide), together with a healthy outgroup cell, and evolves their diploid genomes along it with:

- **point mutations**: binary (0/1) or DNA data, infinite- and finite-sites models (JC, HKY, GTR), rate variation among sites and lineages, germline SNPs, and the 30 COSMIC mutational signatures or your own (`--signature-file`);
- **copy-number events**: point deletions and copy-neutral LOH, and **losses and gains of chromosome segments** (up to 4 copies per site by default) that happen on the branches of the tree, in random order together with the SNVs, so that the order of events matters. The phase of every copy is known and can be written out.

It then simulates the experiment, at three levels of detail (you choose with the options, `--check` tells you which one a run is):

| Mode | What it simulates | Main outputs |
|---|---|---|
| **Genotypes** | the true genotypes, with allelic dropout (ADO) and genotyping error | genotype and haplotype files, trees |
| **Read counts** (`-C`) | the number of reads of each nucleotide in every cell and site, with sequencing and amplification error, allelic imbalance, doublets, and the genotype likelihoods that follow | VCF/BCF, CATG, PILEUP, bulk VAF |
| **Reads** (`-C --bam --fastq`) | paired-end Illumina-like reads of the amplified genome of every cell, as molecules with PCR duplicates, amplification errors and gain, base qualities and strand bias | BAM + index, FASTQ, `ref.fa`, `truth.vcf` |

Other things it can do:

- **Allelic dropout in blocks of sites**, correlated with the coverage and with doublets; coverage that varies between cells, sites and along the genome (`--gain-sd`); amplification errors with a trinucleotide signature (`--amp-signature`).
- **A real or a long chromosome**: simulate a chromosome in blocks on one tree (`--blocks`), take the reference from a FASTA file and a region (`--fasta`, `--region`) and put in your own germline variants (`--germline-vcf`), writing one VCF with real chromosome coordinates, or one BAM per cell along the whole chromosome.
- **Tumour only** (`--no-outgroup`) and the calling filters of a real analysis (`--min-alt-reads`, `--min-alt-vaf`).
- **Very large data sets**: 10,000 cells x 1,000,000 sites in seconds and a few hundred MB (sparse genomes), reads streamed in blocks, parallel with OpenMP and identical results for any number of threads.
- **Experiments**: `--summary` compares what a caller would see with the truth (called fraction, genotype error, precision and recall of the SNVs), and `--sweep` runs a grid of parameter values in parallel and resumes, with plots.
- **Reproducibility**: seeds, `Results/run.json` (version, command line, environment) and configuration files (YAML/JSON) that repeat any run.

CellCoal is validated against theory (coalescent times, site frequency spectrum, depth, errors) and with a test suite that runs on Linux, macOS and Windows; see [Tests](#tests) and [Validation](#validation).

## Getting started

The manual ([html](https://dapogon.github.io/cellcoal/) or [pdf](docs/cellcoal.manual.v2.0.pdf)) has a quickstart (chapter 3) with installation on Linux, macOS and Windows, a first simulation, worked examples and a tutorial on `--sweep`.

**Install.** Download the binary for your system from the [releases](https://github.com/dapogon/cellcoal/releases) page (Linux, macOS and Windows), or install it with Homebrew (`brew tap dapogon/cellcoal https://github.com/dapogon/cellcoal && brew install cellcoal`) or conda (recipe in `conda/`). To compile from the source (the zip or tar.gz of a release, or `git clone`), type `make`: the program is written in C, has no dependencies and appears in `bin/` as `cellcoal-x.y.z`. OpenMP is detected automatically (on macOS `brew install libomp`; `make OMPFLAG=` builds without it). `OMP_NUM_THREADS=n` chooses the threads.

**Run.** Options go on the command line, in a parameter file (`bin/cellcoal-x.y.z` with no arguments reads the file called *parameters*) or in a configuration file (`--config`, see [examples](examples)). Results go to the folder *Results* (or the one given with `-o`). A few examples:

```
# genotypes of 20 cells with allelic dropout, SNV genotypes and trees
cellcoal -n1 -s20 -l1000 -u1e-6 -e10000 -D0.1 -1 -6

# read counts at 10x, a VCF with genotype likelihoods, and a summary against the truth
cellcoal -n10 -s20 -l1000 -u1e-6 -e10000 -C10 -E0.001 --summary

# paired-end reads (BAM and FASTQ) of 20 amplified cells, with a truth VCF
cellcoal -n1 -s20 -l30000 -u1e-7 -e10000 -C10 --ado mean=0.1,block=3000 --reads bam,fastq

# a chromosome of 50 blocks of 200,000 sites on one tree, tumour only; or a region of a real reference
cellcoal -s50 -l200000 -u1e-8 -e10000 -C10 --tree no-outgroup --chromosome blocks=50
cellcoal -s50 -l200000 -u1e-8 -e10000 -C10 --chromosome fasta=hg38.fa.gz,region=chr1:10000001-20000000

# copy-number losses and gains of segments, with the phase in the genotypes
cellcoal -n1 -s20 -l5000 -u1e-6 -e1000 --mutation loss-rate=1e-3,gain-rate=1e-3,segment-length=500 -1 --phase
```

`cellcoal --help` lists the everyday options and the groups of parameters (`--tree`, `--mutation`, `--ado`, `--errors`, `--depth`, `--reads`, `--chromosome`, `--output`, ...); `--help-all` lists every parameter with its configuration key.

## Developed with Claude

**CellCoal 2 was developed together with Claude.** The CellCoal 2 series (versions 2.0.0 onwards) was developed by Claude, an AI model made by Anthropic, working as a coding agent (Claude Code) in collaboration with David Posada. Claude took part in all the work on the code: finding and fixing bugs, making the program faster and lighter in memory, splitting the source into modules, new features (sparse and streamed genomes, BCF output, `--summary`, `--sweep`, `--config`, `--check`, the read-level simulator, a chromosome in blocks and the copy-number model), the test suite, the continuous integration, the packages for Linux, macOS and Windows, the releases and the documentation, including this text. David Posada defined the goals, made the scientific and design decisions, tested the program and approved every release, and is the author of CellCoal and responsible for it. The 2020 article describes version 1, which was written before this work.

## Bulk-like VAF distribution (tools/bulk_vaf.py)

Easiest: option `-0` makes CellCoal write the read counts pooled over the tumour cells, one line per site (VAF from the simulated read counts and the true VAF; size independent of the number of cells), and `python3 tools/bulk_vaf.py --bulk Results/bulk --out vaf.png` plots it (add `-Z` to skip the VCF when you only need the bulk file; add `--depth 100 --purity 0.8` to mimic a real bulk experiment). The details below cover the other inputs.

`tools/bulk_vaf.py` (Python with numpy and matplotlib) draws the variant allele frequency (VAF) distribution of the sample as if it were sequenced in bulk.

- **From the simulated read counts** (VCF made with `-C`, `-V`, ...): the read counts (RC field) of all the tumour cells are pooled, and VAF = alternate reads / total reads at each site. The healthy cell is excluded.
  `python3 tools/bulk_vaf.py --vcf Results/vcf_dir/vcf.0001 --min-reads 10 --out vaf.png --tsv vaf.tsv`
- **From the genotypes** (`-1`): VAF is the fraction of alternate alleles among the tumour cells; a bulk experiment is mimicked by sampling reads at a given depth and tumour purity.
  `python3 tools/bulk_vaf.py Results/snv_gen --depth 100 --purity 1.0 --out vaf.png --tsv vaf.tsv`

With simulated read counts, sequencing and amplification errors give almost every site a few alternate reads, so the VAF of non-mutated sites sits at the error rate near zero; for a clean plot of the true SNVs use `-E0`, or look at VAF above the error level. Pooling thousands of cells at a few reads each gives a very high pooled depth, so the VAF is very precise. The VCF grows with cells x sites (about 360 bytes per cell and site): for thousands of cells keep the number of sites low (a few thousand).

## Versions

A short summary of each series; the details of every release are in [docs/CHANGELOG.md](docs/CHANGELOG.md) and in `docs/CHANGES-X.Y.Z.md`.

- **2.5, copy number on the tree.** Losses and gains of chromosome segments are simulated together with the SNVs, point deletions and cnLOH, branch by branch and in random order, over copies with known phase (up to 4 copies per site by default) and with every substitution model. Genotypes with several copies, the phase (`--phase`), one haplotype row per copy and a truth file of the copies per cell and segment (`cnv_segments`).
- **2.4, reads, amplification and whole chromosomes.** Paired-end Illumina-like reads (BAM, FASTQ) with PCR duplicates, base qualities and a truth VCF; allelic dropout in blocks of sites; signature files and signature-shaped amplification errors; coverage that varies along the genome; a chromosome in blocks, from a real FASTA reference and with your germline variants (`--blocks`, `--fasta`, `--region`, `--germline-vcf`); explicit modes (genotypes, read counts, reads); grouped options (`--ado`, `--depth`, `--output`, ...) and a cleaner code base (2.4.1 and 2.4.2: about 2,300 lines less, same output).
- **2.3, tumour-only runs and long chromosomes.** `--no-outgroup`, calling filters that keep singletons (`--min-alt-reads`, `--min-alt-vaf`), the depth of the non-variable sites (`--coverage`) and `tools/chrom.py`.
- **2.2, experiments at scale.** Faster read counts (`--fast`), BCF output, streaming with every output, `--summary` with metrics against the true genotypes, `--sweep` parameter grids (parallel and resumable), validation against theory, `--check`, and binaries for Linux, macOS and Windows.
- **2.1, very large data sets.** Sparse genome storage (10,000 cells x 1,000,000 sites in about 7 s), streamed read counts with small memory, bulk VAF (`-0`), coverage varying between sites and cells (`-O`, `-J`), configuration files, a modular source, `run.json` run records, Homebrew and conda packages.
- **2.0, a faster, safer CellCoal.** Bug fixes (genotyping error, allelic dropout, memory safety), 3-4 times faster with about 7 times less memory, and parallel read counts. The output for a given seed differs from versions 1.3.4 and earlier.

## Citation

- If you use CellCoal, please cite: Posada D. (2020) CellCoal: coalescent simulation of single-cell sequencing samples. *Molecular Biology and Evolution* 37: 1535-1542. <a href="https://academic.oup.com/mbe/article/37/5/1535/5728646" target="_blank">https://academic.oup.com/mbe/article/37/5/1535/5728646</a>.

## Tests

`make test` runs the regression suite (`tests/run_tests.sh`): the example files, reproducibility and identical results for 1 and 2 threads, the observed rates of ADO, genotyping error, doublets and allelic imbalance, the mutation models and signatures (also with copy-number events), the equivalence of dense, sparse and streamed genomes, the read counts, the reads and the copy-number outputs. It does not depend on exact floating-point output, so it runs on any platform. The CI also builds with clang and without OpenMP, runs the suite under AddressSanitizer/UBSan, fuzzes the options and the `--config` parser (`tests/fuzz.py`) and reports code coverage; `tests/theory.py` validates the model against theory (coalescent waiting times, growth, number and spectrum of mutations, depth, sequencing error and dropout). See section 6.10h of the manual.

## Validation
A script that simulates with CellCoal, infers the tree with CellPhy and compares it with the true tree is in [examples/validation](examples/validation).

## Manual
A detailed manual is available in [html](https://dapogon.github.io/cellcoal/) or in [pdf](docs/cellcoal.manual.v2.0.pdf). To start, read the [quickstart](https://dapogon.github.io/cellcoal/cellcoal.manual.v2.0.html#toc-quickstart) (installation and examples).

## License

CellCoal is free software, distributed under the GNU General Public License, version 2 or (at your option) any later version; see [LICENSE](LICENSE).
