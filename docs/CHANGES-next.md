# Unreleased (after 2.5.3)

- Fix: the long options of the table rejected valid values. `--qual-slope` (documented range -6 to 6) rejected negative values and 0, and `--qual-sd 0`, `--strand-bias 0`, `--dup-mean 0`, `--seg-loss-rate 0`, `--seg-gain-rate 0` were rejected too, which broke `--sweep` over a range that includes 0 and configuration files that write the default explicitly. Values out of range are still rejected.
- `--help`: the long options were listed twice in the table of keys of the configuration file.
- `tests/fuzz.py` also gives random values to the long options (copy number, qualities, dropout in segments, correlations).
- `--help`: the table of configuration keys showed `--sweep_keep` and `--no_outgroup`; the options are `--sweep-keep` and `--no-outgroup`.
- `--check` shows the copy-number settings when losses or gains are on.
- New option `--compress-level N` (1 to 9, default 6, which gives the same files as before): the zlib level of the BGZF files (BAM, BCF and FASTQ.gz). The compression is more than half of the time of the reads mode; level 1 makes it about 40% faster with files about 20% larger. The content is the same.
