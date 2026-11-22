# Unreleased (after 2.5.3)

- Fix: the long options of the table rejected valid values. `--qual-slope` (documented range -6 to 6) rejected negative values and 0, and `--qual-sd 0`, `--strand-bias 0`, `--dup-mean 0`, `--seg-loss-rate 0`, `--seg-gain-rate 0` were rejected too, which broke `--sweep` over a range that includes 0 and configuration files that write the default explicitly. Values out of range are still rejected.
- `--help`: the long options were listed twice in the table of keys of the configuration file.
- `tests/fuzz.py` also gives random values to the long options (copy number, qualities, dropout in segments, correlations).
- `--help`: the table of configuration keys showed `--sweep_keep` and `--no_outgroup`; the options are `--sweep-keep` and `--no-outgroup`.
- `--check` shows the copy-number settings when losses or gains are on.
- New option `--compress-level N` (1 to 9, default 6, which gives the same files as before): the zlib level of the BGZF files (BAM, BCF and FASTQ.gz). The compression is more than half of the time of the reads mode; level 1 makes it about 40% faster with files about 20% larger. The content is the same.
- **The reads mode (`--bam`, `--fastq`) is parallel by cells.** Every cell and the healthy bulk sample have their own random stream (taken from the main stream, like the sites of the read counts), so the files are the same for any number of threads and the time falls with the number of cores (8 cells, 1 Mb, 30X: 17 s with 1 thread, 8.9 s with 2). **The reads of a given seed differ from those of 2.5.x** (statistically the same: depth, duplicates and error rate); the genotypes, trees, VCF and read counts are unchanged.
