# CellCoal 2.4.2 — what changed

**The simulated data of all existing runs are identical to 2.4.1.**

## Command line and help
- `--help` lists every option by sections and ends with the configuration-file key of each one. New examples `examples/config.genotypes.yaml`, `config.readcounts.yaml`, `config.reads.yaml` (manual 6.10, test 2g44). Help typos fixed (`-h`, `-k` examples) and a NUL byte in the defaults listing.
- Related options can now be given as one grouped option: `--tree`, `--mutation`, `--substitution`, `--ado`, `--errors`, `--amplification`, `--depth`, `--cell-corr`, `--call-filter`, `--reads`, `--chromosome`, `--output`. Syntax: `--GROUP key=value,key,...` (':' separates several numbers), e.g. `--ado cells=0.1:0.001,sites=0.2:0.01 --depth mean=20,cells=0.3 --output snv-genotypes,vcf-gz`.
- The single-letter options keep working and give identical results. The same group names are accepted as keys in configuration files.
- `--help` is organized by group and includes a table of equivalences with the single-letter options; error messages use the group key names.

## Internal
- Internal cleanup (about 2,300 lines less in `src/`): one parser for the command line and the parameter and configuration files (a parameter file is turned into the equivalent command line), `main()` and `RunReplicate` split into setup and replicate phases, `Checked()` for allocations, one table of result files and one row-label helper in `output.c`, `GenotypeLikelihoods` split into one function per model, `SiteReadCounts`, `SimulateGenome`, `VCFCellFields` and the BCF writer split into helpers, dead debug code removed, one function draws the trinucleotide class of a signature mutation.
- The amplification substitution rows are computed once instead of at every use (a negligible speed-up).
