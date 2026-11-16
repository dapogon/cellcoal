# Unreleased (after 2.2.3)

- Internal: `main` is split further (`RunReplicate`), the statistics accumulated over replicates are in a `RunStats` structure (about 30 fewer globals), and `RemoveDir` no longer needs `nftw`.
- Windows (MSYS2/MinGW-w64): the program builds and runs; the random number generator and the index arithmetic use 64-bit integers explicitly (`long` is 32 bits on Windows). New CI job `windows` (build, smoke tests, whole suite; 6 of 80 checks still fail there and are being investigated). See section 6.10j of the manual.

The output on Linux and macOS is identical to 2.2.3.
