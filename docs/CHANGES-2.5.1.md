# CellCoal 2.5.1 — what changed

**The simulated data of all existing runs are identical to 2.5.0.**

## Fixes
- The Windows build works again: the read-level simulator (`--bam`, `--fastq`) called `mkdir` with a POSIX signature, which MSYS2 rejects (it now uses the portable `MAKE_DIR`). Versions 2.4.0 to 2.5.0 therefore have Linux and macOS binaries but no Windows zip; 2.5.1 has all three.
