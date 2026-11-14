#!/bin/bash
# CellCoal regression tests (portable: they do not depend on exact floating-point output).
# Usage: tests/run_tests.sh [path/to/cellcoal-binary]   (or: make test)
#  1. every example parameter file runs without error;
#  2. with the same seed, results are identical for 1 and 2 threads (read counts, VCF);
#  3. germline SNPs (-c) are inherited by all cells and the healthy cell;
#  4. the same seed gives the same output twice.
BIN=${1:-$(ls -t "$(dirname "$0")"/../bin/cellcoal-* | head -1)}
BIN=$(cd "$(dirname "$BIN")" && pwd)/$(basename "$BIN")
ROOT=$(cd "$(dirname "$0")/.." && pwd)
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
fail=0
ok()  { echo "  ok   $1"; }
bad() { echo "  FAIL $1"; fail=1; }
# print the results of a run without lines that change between runs (dates, times, command line)
norm() { (cd "$1" && find . -type f | sort | xargs grep -v -i -e "fileDate" -e "Time processing" -e "date" -e "Command" -e "time" -e "CellCoal" -e "[0-9][0-9]:[0-9][0-9]:[0-9][0-9]" 2>/dev/null); }

echo "Binary: $BIN"
echo "1. Example parameter files"
for p in "$ROOT"/examples/parameters.*; do
  n=$(basename "$p"); d="$TMP/ex_$n"; mkdir -p "$d"
  cp "$ROOT"/examples/usergenome.fas "$ROOT"/examples/usertree.0001 "$d"/ 2>/dev/null
  # fix the seed so that the example is deterministic
  sed 's/^\[seed\] \[#[0-9]*\]/[seed] #1234/' "$p" > "$d/par"
  (cd "$d" && "$BIN" -Fpar >out.txt 2>&1) && ok "$n" || bad "$n"
done

echo "2. Same results with 1 and 2 threads"
A="-n2 -s30 -l500 -u1e-6 -e10000 -g1e-4 -b1 -m2 -6 -v -C30 -V5 -E0.001 -D0.1 -A0.05 0.01 0 -#77"
for t in 1 2; do mkdir -p "$TMP/t$t"; (cd "$TMP/t$t" && OMP_NUM_THREADS=$t "$BIN" $A >out.txt 2>&1); done
if [ -d "$TMP/t1/Results" ] && [ "$(norm "$TMP/t1/Results" | md5sum)" = "$(norm "$TMP/t2/Results" | md5sum)" ]; then ok "1 vs 2 threads"; else bad "1 vs 2 threads"; fi

echo "3. Germline SNPs (-c) are inherited by every cell"
for b in 0 1; do
  d="$TMP/germ$b"; mkdir -p "$d"
  (cd "$d" && "$BIN" -n1 -s4 -l200 -c0.1 -u0 -e10000 -b$b -4 -v -#3 >out.txt 2>&1)
  f=$(ls "$d"/Results/full_haplotypes_dir/* 2>/dev/null | head -1)
  res=$(python3 - "$f" <<'PY'
import sys
rows=[l.split() for l in open(sys.argv[1]).read().split('\n')[1:] if l.strip()]
seq={r[0]:r[1] for r in rows}
cells=sorted({k[:-1] for k in seq})
# every cell (and the healthy cell) must carry the same maternal/paternal differences
diff=lambda c:[i for i,(a,b) in enumerate(zip(seq[c+'m'],seq[c+'p'])) if a!=b]
ref=diff(cells[0])
print('OK' if ref and all(diff(c)==ref for c in cells) else 'BAD')
PY
)
  [ "$res" = "OK" ] && ok "alphabet -b$b" || bad "alphabet -b$b"
done

echo "4. Reproducibility"
for r in 1 2; do mkdir -p "$TMP/r$r"; (cd "$TMP/r$r" && "$BIN" $A >out.txt 2>&1); done
[ "$(norm "$TMP/r1/Results" | md5sum)" = "$(norm "$TMP/r2/Results" | md5sum)" ] && ok "same seed, same output" || bad "same seed, same output"

[ $fail -eq 0 ] && echo "All tests passed" || echo "Some tests FAILED"
exit $fail
