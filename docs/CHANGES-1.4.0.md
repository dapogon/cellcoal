# CellCoal 1.4.0 — what changed

Release 1.4.0 fixes several bugs, makes the program much faster and lighter on memory, and parallelizes read-count generation. **Output for a given random seed is different from versions ≤ 1.3.4** (see "Reproducibility").

## 1. Bug fixes that change the simulated data

| Area | Problem | Fix |
|---|---|---|
| Genotyping error (`-G`) | The rows of the error-weight table (default `{0,1,1,1}`, sum 3) were not normalized, but the sampler expects probabilities summing to 1, so the first non-zero entry was always chosen. | Rows are normalized before sampling; the erroneous base now follows the intended probabilities. |
| Allelic dropout (`-P`, `-Q`) | ADO variation between cells and between sites was not simulated correctly. | ADO rate is drawn once per cell and once per site and reused (per-cell and per-site arrays). |
| VCF `PL` / `PLN` | Values were printed with `%d` from a `double`; infinite likelihoods gave undefined output. | Printed correctly; non-finite values are written as `.`. |

## 2. Robustness fixes (no change to valid output)

- Heap overflow when reading a user genome, use-after-free when reading a user tree, and memory leaks in the signature simulation (found with AddressSanitizer / UBSan).
- Missing option arguments (e.g. `-n` as last argument) are now reported instead of reading past `argv`.
- Bounded name/file reads, safe `realloc` (no lost block on failure), checks after `calloc`/`fgets`/`fscanf`.
- About 78 compiler warnings removed; the build is warning-free.

## 3. Speed and memory

Benchmark (300 cells × 3000 sites, genotype likelihoods, ADO, doublets, 1 replicate):

| Version | Wall time | Peak memory |
|---|---|---|
| 1.3.4 | 17–19 s | ≈ 1 GB |
| 1.4.0, 1 thread | ≈ 4.9 s | ≈ 134 MB |
| 1.4.0, 2 threads | ≈ 3.4 s | ≈ 136 MB |

What was done:

- Compact genotype-likelihood storage (11 distinct values instead of the full 4×4 tables).
- A budgeted cache of per-site genotype likelihoods. The budget defaults to 256 MB and can be changed with the environment variable `CELLCOAL_GL_CACHE_MB` (`0` = always recompute). It only affects speed/memory, never the output.
- Exact memoization of `log10` calls (bit-identical results).
- Faster VCF text formatting; VCF printing is parallel.
- Read-count generation runs in parallel (OpenMP).
- Compiler flags: `-O3 -ffp-contract=off -fno-math-errno -fopenmp`. `-ffp-contract=off` keeps floating-point results identical across compilers/CPUs; do **not** add `-ffast-math` or `-march=native` if you need identical output.

### Threads

OpenMP is enabled by default. Use `OMP_NUM_THREADS=n` to choose the number of threads. To build without OpenMP: `make OMPFLAG=`.

## 4. Reproducibility

Read counts are now generated with one independent random stream per site (SplitMix64, seeded from the main generator and the site index). Consequences:

- For the same seed and the same options, output is identical regardless of the number of threads and across repeated runs.
- Output for a given seed differs from versions ≤ 1.3.4 (the sampling scheme changed). Distributions are statistically equivalent: mean depth, fraction of zero-depth sites and major-allele fraction were compared over 8 seeds against the previous build and agree within noise.
- Everything else (tree, mutations, true genotypes) is still drawn from the original generator.
- As before, the default seed is time-based; pass `-#N` for reproducible runs. If options are read with `-F file`, the seed must be inside the file.

## 5. Validation performed

- Bit-identical output to the original program for all code paths except the deliberate behaviour changes listed in section 1.
- AddressSanitizer + UBSan clean on the examples and extra configurations.
- ThreadSanitizer: only fork/join boundary reports from the (uninstrumented) OpenMP runtime; no races between worker threads.
- Same output with 1, 2, 3, 4, 5 and 7 threads.
