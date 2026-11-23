# Unreleased (after 2.6.1)
- `--damage SPEC` (reads mode): damage of the library before the PCR, as `FROM>TO:RATE` items (`G>T:0.002,C>T:0.001`). The change is in one strand of the molecule, so it is seen in only one orientation of the reads (orientation bias of oxidation and deamination) and the PCR duplicates inherit it; it also affects the healthy bulk. Without the option nothing changes.

