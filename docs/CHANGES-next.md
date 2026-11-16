# Unreleased (after 2.2.4)

- `--summary` adds `het_error`, `hom_error`, `het_to_hom` (the footprint of allelic dropout), `doublet_error` and `singlet_error`.
- Windows: OpenMP is enabled automatically, and the release workflow builds a Windows zip (program plus the compiler DLLs).
- Speed: the VCF formatting and the read simulation were profiled again; no change was worth its risk, the output is identical.
