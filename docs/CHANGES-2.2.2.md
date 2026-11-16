# CellCoal 2.2.2 — what changed

**The output of all existing runs is bit-identical to 2.2.1.**

## Fixes (found by the new fuzzer and sanitizers)

- **BCF writer** (`--bcf`, 2.2.0 and 2.2.1): a 2-byte heap overflow at the end of a record. The files were valid, but memory was corrupted.
- **Hangs** with enormous gamma shapes (`-a1e308`, `-J`, `-O`, `-V`, `-i`): the gamma generator now uses its normal limit above a shape of 1e7.
- Two **integer overflows** with absurd values (mean of a Poisson above 2e9; printing of the demographic periods).

## Checks

CI with clang and without OpenMP, the whole suite under AddressSanitizer and UBSan, a fuzzer of the options and of `--config` (`tests/fuzz.py`), line coverage (81%), and `tests/theory.py`, which validates the model against theory (coalescent, growth, mutations, site frequency spectrum, depth, sequencing error, dropout). A failing `make test` now fails the CI job (before, `tee` hid the status). See section 6.10h of the manual.
