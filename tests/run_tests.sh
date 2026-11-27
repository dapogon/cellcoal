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
# portability (macOS): wc pads its numbers with spaces and stat has no -c
wc() { command wc "$@" | sed 's/^ *//'; }
stat() { if [ "$1" = "-c" ] && [ "$2" = "%s" ] && ! command stat -c %s / >/dev/null 2>&1; then shift 2; command stat -f %z "$@"; else command stat "$@"; fi; }
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
fail=0
ok()  { echo "  ok   $1"; }
bad() { echo "  FAIL $1"; fail=1; }
# print the results of a run without lines that change between runs (dates, times, command line)
close() { python3 -c "import sys; sys.exit(0 if abs($1-$2)<=$3 else 1)"; }
norm() { (cd "$1" && find . -type f ! -name run.json | sort | xargs grep -v -i -e "fileDate" -e "Time processing" -e "date" -e "Command" -e "time" -e "CellCoal" -e "[0-9][0-9]:[0-9][0-9]:[0-9][0-9]" 2>/dev/null); }

echo "Binary: $BIN"
echo "1. Example parameter files"
for p in "$ROOT"/examples/parameters.*; do
  n=$(basename "$p"); d="$TMP/ex_$n"; mkdir -p "$d"
  cp "$ROOT"/examples/usergenome.fas "$ROOT"/examples/usertree.0001 "$d"/ 2>/dev/null
  # fix the seed so that the example is deterministic
  sed 's/\[seed\] \[#[0-9]*\]/[seed] #1234/' "$p" > "$d/par"
  (cd "$d" && "$BIN" -Fpar >out.txt 2>&1) && ok "$n" || bad "$n"
done

echo "2. Same results with 1 and 2 threads"
A="-n2 -s30 -l500 -u1e-6 -e10000 -g1e-4 -b1 -m2 -6 -v -C30 -V5 -E0.001 -D0.1 -A0.05 0.01 0 -#77"
for t in 1 2; do mkdir -p "$TMP/t$t"; (cd "$TMP/t$t" && OMP_NUM_THREADS=$t "$BIN" $A >out.txt 2>&1); done
if [ -d "$TMP/t1/Results" ] && [ "$(norm "$TMP/t1/Results" | md5sum)" = "$(norm "$TMP/t2/Results" | md5sum)" ]; then ok "1 vs 2 threads"; else bad "1 vs 2 threads"; fi

echo "2b. Sparse genome storage gives the same results as dense storage (read-free outputs)"
B="-n2 -s40 -l3000 -u1e-6 -e10000 -g1e-4 -b1 -m2 -p0.5 -H1e-7 -d1e-7 -c1e-3 -1 -2 -3 -4 -Y -6 -#12"
mkdir -p "$TMP/dn" "$TMP/sp"
(cd "$TMP/dn" && CELLCOAL_SPARSE=0 "$BIN" $B >out.txt 2>&1); (cd "$TMP/sp" && CELLCOAL_SPARSE=1 "$BIN" $B >out.txt 2>&1)
if [ -d "$TMP/dn/Results" ] && [ "$(norm "$TMP/dn/Results" | md5sum)" = "$(norm "$TMP/sp/Results" | md5sum)" ]; then ok "dense vs sparse"; else bad "dense vs sparse"; fi
echo "2c. A large data set runs in sparse mode (2000 cells x 1,000,000 sites)"
mkdir -p "$TMP/big"
if (cd "$TMP/big" && "$BIN" -n1 -s2000 -l1000000 -u1e-8 -e10000 -g0 -b1 -D0.1 -G1e-6 1e-13 -1 -#5 >out.txt 2>&1) && [ -n "$(find "$TMP/big/Results/snv_gen" -type f -size +1k 2>/dev/null)" ]; then ok "2000 x 1e6 sites"; else bad "2000 x 1e6 sites"; fi

echo "2d. Streamed reads (blocks of sites) give the same VCF and PILEUP as storing all the reads"
C="-n2 -s30 -l3000 -u1e-6 -e10000 -g1e-4 -b1 -m0 -v -C20 -V5 -E0.001 -A0.05 0.01 0 -9 -P0.1 0.01 -#79"
mkdir -p "$TMP/rs0" "$TMP/rs1"
(cd "$TMP/rs0" && CELLCOAL_SPARSE=1 CELLCOAL_STREAM=0 "$BIN" $C >out.txt 2>&1); (cd "$TMP/rs1" && CELLCOAL_STREAM=1 "$BIN" $C >out.txt 2>&1)
if [ -d "$TMP/rs0/Results/vcf_dir" ] && [ "$(norm "$TMP/rs0/Results" | md5sum)" = "$(norm "$TMP/rs1/Results" | md5sum)" ]; then ok "stored vs streamed reads"; else bad "stored vs streamed reads"; fi

echo "2e. Pooled (bulk) read counts and VAF (-0), also streamed"
C="-n2 -s60 -l3000 -u1e-6 -e10000 -g1e-4 -b1 -C5 -V1 -E0 -0 -#7"
mkdir -p "$TMP/bk0" "$TMP/bk1"
(cd "$TMP/bk0" && CELLCOAL_STREAM=0 "$BIN" $C >out.txt 2>&1); (cd "$TMP/bk1" && CELLCOAL_STREAM=1 "$BIN" $C >out.txt 2>&1)
if [ -s "$TMP/bk0/Results/bulk" ] && cmp -s "$TMP/bk0/Results/bulk" "$TMP/bk1/Results/bulk"; then ok "bulk file, stored vs streamed"; else bad "bulk file, stored vs streamed"; fi

echo "2f. -Z: no VCF, only the bulk file"
mkdir -p "$TMP/nz"
(cd "$TMP/nz" && "$BIN" $C -Z >out.txt 2>&1)
if [ -s "$TMP/nz/Results/bulk" ] && [ ! -e "$TMP/nz/Results/vcf" ] && cmp -s "$TMP/nz/Results/bulk" "$TMP/bk0/Results/bulk"; then ok "-Z"; else bad "-Z"; fi

echo "2g. Bulk-only runs skip the likelihoods; -z compresses the VCF; -M filters it"
D="-n2 -s60 -l3000 -u1e-6 -e10000 -g1e-4 -b1 -C5 -V1 -E0.001 -B0.3 0.01 -A0.05 0.01 0 -#7"
mkdir -p "$TMP/ro" "$TMP/st" "$TMP/gz" "$TMP/pl" "$TMP/fl" "$TMP/fs"
(cd "$TMP/ro" && "$BIN" $D -0 -Z >out.txt 2>&1); (cd "$TMP/st" && CELLCOAL_STREAM=1 "$BIN" $D -0 -9 >out.txt 2>&1)
if [ -s "$TMP/ro/Results/bulk" ] && cmp -s "$TMP/ro/Results/bulk" "$TMP/st/Results/bulk"; then ok "bulk only (no likelihoods) = full reads"; else bad "bulk only (no likelihoods) = full reads"; fi
(cd "$TMP/pl" && "$BIN" $D >out.txt 2>&1); (cd "$TMP/gz" && "$BIN" $D -z >out.txt 2>&1)
if [ "$(gzip -dc "$TMP/gz/Results/vcf.gz" | grep -v -i filedate | md5sum)" = "$(grep -v -i filedate "$TMP/pl/Results/vcf" | md5sum)" ]; then ok "-z"; else bad "-z"; fi
(cd "$TMP/fl" && "$BIN" $D -M0.02 >out.txt 2>&1); (cd "$TMP/fs" && CELLCOAL_STREAM=1 "$BIN" $D -M0.02 >out.txt 2>&1)
nf=$(grep -a -vc '^#' "$TMP/fl/Results/vcf"); np=$(grep -a -vc '^#' "$TMP/pl/Results/vcf")
if [ "$nf" -gt 0 ] && [ "$nf" -lt "$np" ]; then ok "-M reduces the VCF ($np -> $nf sites)"; else bad "-M reduces the VCF"; fi

echo "2g2. Streamed CATG and ML genotype/haplotype files (-1 -2 -8) equal the stored ones, any thread count"
C="-n2 -s30 -l3000 -u1e-6 -e10000 -g1e-4 -b1 -m0 -v -C20 -V5 -E0.001 -A0.05 0.01 0 -B0.3 0.001 -1 -2 -8 -9 -#79"
mkdir -p "$TMP/ml0" "$TMP/ml1" "$TMP/ml2"
(cd "$TMP/ml0" && CELLCOAL_STREAM=0 "$BIN" $C >out.txt 2>&1); (cd "$TMP/ml1" && CELLCOAL_STREAM=1 OMP_NUM_THREADS=1 "$BIN" $C >out.txt 2>&1); (cd "$TMP/ml2" && CELLCOAL_STREAM=1 OMP_NUM_THREADS=3 "$BIN" $C >out.txt 2>&1)
if [ -n "$(find "$TMP/ml0/Results/catg_dir" -type f -size +1k)" ]; then
  if [ "$(norm "$TMP/ml0/Results" | md5sum)" = "$(norm "$TMP/ml1/Results" | md5sum)" ] && [ "$(norm "$TMP/ml1/Results" | md5sum)" = "$(norm "$TMP/ml2/Results" | md5sum)" ]; then ok "CATG, SNV genotypes and haplotypes: stored = streamed (1 and 3 threads)"; else bad "CATG/ML genotypes stored vs streamed"; fi
else bad "CATG file not produced"; fi

echo "2g3. Coverage dispersion between sites (-O) and ADO counted at every site in sparse mode"
Q="-n1 -s100 -l2000 -u2e-6 -e10000 -g1e-4 -b1 -C10 -V5 -E0.02 -0 -Z"
mkdir -p "$TMP/o0" "$TMP/o1" "$TMP/o2"
(cd "$TMP/o0" && "$BIN" $Q -#4 >/dev/null 2>&1); (cd "$TMP/o1" && CELLCOAL_STREAM=0 "$BIN" $Q -O3 -#4 >/dev/null 2>&1); (cd "$TMP/o2" && CELLCOAL_STREAM=1 "$BIN" $Q -O3 -#4 >/dev/null 2>&1)
cv() { python3 -c "
import statistics as st,sys
d=[float(l.split('\t')[6]) for l in open(sys.argv[1]) if l[:1].isdigit() and l.count('\t')>8]
m=st.mean(d); print(round(st.pvariance(d)/m**2,4))" "$1"; }
c0=$(cv "$TMP/o0/Results/bulk"); c1=$(cv "$TMP/o1/Results/bulk")
close "$c1" 0.3333 0.06 && ok "-O3: squared coefficient of variation of the site depth $c1 (expected 0.333; $c0 without -O)" || bad "-O3 depth dispersion ($c1)"
if cmp -s "$TMP/o1/Results/bulk" "$TMP/o2/Results/bulk"; then ok "-O, stored vs streamed bulk file"; else bad "-O, stored vs streamed bulk file"; fi
mkdir -p "$TMP/ad"
(cd "$TMP/ad" && CELLCOAL_SPARSE=1 "$BIN" -n1 -s20 -l1000 -u1e-6 -e10000 -g1e-4 -b1 -y3 -D0.3 -#4 >o.txt 2>&1)
fr=$(python3 -c "
import re,sys
t=open(sys.argv[1],errors='ignore').read()
v=[int(x) for x in re.findall(r'countDropped = (\d+)',t)]
print(round(sum(v)/(len(v)*42.0),4))" "$TMP/ad/o.txt")
close "$fr" 0.1634 0.02 && ok "sparse ADO counted at invariant sites: allele dropout fraction $fr (expected 0.1634)" || bad "sparse ADO count ($fr)"

echo "2g4. Coverage heterogeneity: -V is noise independent in every cell (averages out in the bulk), -O is shared by all the cells (does not)"
cvstat() { python3 -c "
import statistics as st,sys
dp=[]
for l in open(sys.argv[1]+'/vcf',errors='ignore'):
    if l.startswith('1\t'):
        f=l.rstrip('\n').split('\t'); dp+=[int(x.split(':')[1]) for x in f[9:-1]]
b=[float(l.split('\t')[6]) for l in open(sys.argv[1]+'/bulk') if l[:1].isdigit() and l.count('\t')>8]
cv=lambda d: st.pvariance(d)/st.mean(d)**2
print(round(cv(dp),4), round(cv(b),4))" "$1"; }
Q="-n1 -s100 -l1500 -u1e-6 -e10000 -g1e-4 -b1 -C10 -E0.05 -0"
for k in V O VO; do mkdir -p "$TMP/cv$k"; done
(cd "$TMP/cvV" && "$BIN" $Q -V4 -#4 >/dev/null 2>&1); (cd "$TMP/cvO" && "$BIN" $Q -O3 -#4 >/dev/null 2>&1); (cd "$TMP/cvVO" && "$BIN" $Q -V4 -O3 -#4 >/dev/null 2>&1)
read cV bV <<< "$(cvstat "$TMP/cvV/Results")"; read cO bO <<< "$(cvstat "$TMP/cvO/Results")"; read cVO bVO <<< "$(cvstat "$TMP/cvVO/Results")"
close "$cV" 0.225 0.03 && ok "-V4: depth CV^2 per cell $cV (expected 0.225)" || bad "-V4 per-cell depth CV^2 ($cV)"
close "$bV" 0.0 0.01 && ok "-V4: bulk depth CV^2 $bV (independent noise averages out)" || bad "-V4 bulk CV^2 ($bV)"
close "$bO" 0.333 0.05 && ok "-O3: bulk depth CV^2 $bO (expected 0.333, shared by all the cells)" || bad "-O3 bulk CV^2 ($bO)"
close "$bVO" 0.333 0.05 && ok "-V4 -O3: bulk depth CV^2 $bVO (only -O remains; per cell $cVO > $cV)" || bad "-V4 -O3 bulk CV^2 ($bVO)"
python3 -c "import sys; sys.exit(0 if $cVO > $cO > $cV else 1)" && ok "per-cell depth CV^2 grows when the sources are combined ($cV, $cO, $cVO)" || bad "combined per-cell CV^2"

echo "2g5. Coverage between cells (-J) and doublets: same data stored and streamed"
Q="-n1 -s100 -l1500 -u1e-6 -e10000 -g1e-4 -b1 -C10 -E0.05 -0 -J4"
mkdir -p "$TMP/j0" "$TMP/j1" "$TMP/jd0" "$TMP/jd1"
(cd "$TMP/j0" && CELLCOAL_STREAM=0 "$BIN" $Q -9 -#4 >/dev/null 2>&1); (cd "$TMP/j1" && CELLCOAL_STREAM=1 "$BIN" $Q -9 -#4 >/dev/null 2>&1)
if [ "$(norm "$TMP/j0/Results" | md5sum)" = "$(norm "$TMP/j1/Results" | md5sum)" ] && [ -s "$TMP/j0/Results/bulk" ]; then ok "-J: stored = streamed"; else bad "-J: stored vs streamed"; fi
(cd "$TMP/jd0" && CELLCOAL_STREAM=0 "$BIN" $Q -B0.3 0.01 -1 -8 -9 -#4 >/dev/null 2>&1); (cd "$TMP/jd1" && CELLCOAL_STREAM=1 "$BIN" $Q -B0.3 0.01 -1 -8 -9 -#4 >/dev/null 2>&1)
if [ "$(norm "$TMP/jd0/Results" | md5sum)" = "$(norm "$TMP/jd1/Results" | md5sum)" ] && [ -s "$TMP/jd0/Results/bulk" ]; then ok "doublets: stored = streamed (same random numbers)"; else bad "doublets: stored vs streamed"; fi
cellcv=$(python3 -c "
import statistics as st,sys
n=[0]*101; k=0
for l in open(sys.argv[1],errors='ignore'):
    if l.startswith('1\t'):
        k+=1
        f=l.rstrip('\n').split('\t')
        for c,x in enumerate(f[9:-1]): n[c]+=int(x.split(':')[1])
m=[x/k for x in n[:100]]
print(round(st.pvariance(m)/st.mean(m)**2,4))" "$TMP/j0/Results/vcf")
close "$cellcv" 0.25 0.08 && ok "-J4: CV^2 of the mean depth between cells $cellcv (expected 0.25)" || bad "-J4 between-cell depth ($cellcv)"

echo "2g6. Configuration files (--config): JSON, YAML and the options of run.json reproduce a command line"
Q="-n2 -s20 -l500 -u1e-6 -e10000 -b1 -C10 -E0.001 -A0.05 0.01 0 -1 -9 -v -#7"
mkdir -p "$TMP/cf0" "$TMP/cf1" "$TMP/cf2" "$TMP/cf3"
(cd "$TMP/cf0" && "$BIN" $Q >/dev/null 2>&1)
(cd "$TMP/cf1" && "$BIN" --config "$TMP/cf0/Results/run.json" >/dev/null 2>&1)
cat > "$TMP/c.yaml" <<'YML'
# same run as a YAML file
replicates: 2
cells: 20
sites: 500
mutation_rate: 1e-6
population_size: 10000
alphabet: 1
coverage: 10
sequencing_error: 0.001
amplification_error: [0.05, 0.01, 0]
print_snv_genotypes: true
print_pileup: true
replicate_folders: true
seed: 7
YML
(cd "$TMP/cf2" && "$BIN" --config "$TMP/c.yaml" >/dev/null 2>&1)
if [ -d "$TMP/cf0/Results/pileup_dir" ] && [ "$(norm "$TMP/cf0/Results" | md5sum)" = "$(norm "$TMP/cf1/Results" | md5sum)" ]; then ok "--config run.json reproduces the run"; else bad "--config run.json"; fi
if [ "$(norm "$TMP/cf0/Results" | md5sum)" = "$(norm "$TMP/cf2/Results" | md5sum)" ]; then ok "--config YAML with long names = command line"; else bad "--config YAML"; fi
(cd "$TMP/cf3" && "$BIN" --config "$TMP/c.yaml" -n1 >/dev/null 2>&1); n=$(ls "$TMP/cf3/Results/snv_genotypes_dir" 2>/dev/null | wc -l | tr -d " ")
[ "$n" = "1" ] && ok "options after --config override the file (-n1)" || bad "override after --config ($n)"
echo '{"foo": 1}' > "$TMP/bad.json"
if "$BIN" --config "$TMP/bad.json" 2>&1 | grep -q 'unknown option: foo'; then ok "unknown option in a config file is an error"; else bad "unknown option in config"; fi

printf '{\n  "cells": 20,\n  "cell": 5\n}\n' > "$TMP/typo.json"
msg=$("$BIN" --config "$TMP/typo.json" 2>&1)
if echo "$msg" | grep -q 'line 3, column 3: unknown option: cell' && echo "$msg" | grep -q "did you mean 'cells'" && echo "$msg" | grep -q '"cell": 5'; then ok "config typo: line, column, source line and 'did you mean'"; else bad "config typo message"; fi
printf '# c\ncells: 20\n  mutation_rte: 1e-6\n' > "$TMP/typo.yaml"
if "$BIN" --config "$TMP/typo.yaml" 2>&1 | grep -q "line 3, column 3: unknown option: mutation_rte.*did you mean 'mutation_rate'"; then ok "config YAML typo is located and corrected"; else bad "config YAML typo message"; fi
printf '{\n  "cells": 20\n  "sites": 100\n}\n' > "$TMP/syn.json"
if "$BIN" --config "$TMP/syn.json" 2>&1 | grep -q "line 3, column 3: ',' or '}' was expected"; then ok "config JSON syntax error gives line and column"; else bad "config JSON syntax message"; fi

echo "2g7. Summary plots (tools/cellcoal_plots.py)"
if python3 -c "import matplotlib, numpy" 2>/dev/null; then
  mkdir -p "$TMP/pl2"
  (cd "$TMP/pl2" && "$BIN" -n1 -s30 -l2000 -u3e-6 -e10000 -b1 -C8 -V5 -E0.001 -D0.1 -1 -v -#3 >/dev/null 2>&1)
  ok_all=1
  for c in "depth Results/vcf_dir/vcf.0001" "missing Results/vcf_dir/vcf.0001" "heatmap Results/snv_genotypes_dir/snv_gen.0001"; do
    set -- $c
    (cd "$TMP/pl2" && python3 "$ROOT/tools/cellcoal_plots.py" $1 $2 --out "$1.png" >/dev/null 2>&1) || ok_all=0
    [ -s "$TMP/pl2/$1.png" ] || ok_all=0
  done
  [ $ok_all = 1 ] && ok "depth, missing and heatmap plots are written" || bad "summary plots"
  mkdir -p "$TMP/pl3"
  (cd "$TMP/pl3" && "$BIN" -n1 -s30 -l2000 -u1e-6 -e10000 -b1 -C8 -V5 -E0.001 -D0.1 -1 -0 -v -#3 >/dev/null 2>&1 && python3 "$ROOT/tools/cellcoal_plots.py" all Results --outdir plots >plots.txt 2>&1)
  if [ -s "$TMP/pl3/plots/report.html" ] && grep -q 'base64' "$TMP/pl3/plots/report.html" && [ -s "$TMP/pl3/plots/vaf.png" ] && [ -s "$TMP/pl3/plots/heatmap.png" ]; then ok "'all' writes the plots and a self-contained report.html"; else bad "cellcoal_plots.py all"; fi
else echo "  skip (matplotlib not installed)"; fi

echo "2g8. Fast mode (--fast)"
FQ="-n1 -s40 -l1000 -u1e-6 -e10000 -b1 -C20 -V5 -E0.005 -D0.1 -A0.05 0.01 0 -B0.1 0.005 -1 -2 -8 -#31"
mkdir -p "$TMP/fa1" "$TMP/fa2" "$TMP/fa3" "$TMP/fa4"
(cd "$TMP/fa1" && OMP_NUM_THREADS=1 "$BIN" --fast $FQ >/dev/null 2>&1)
(cd "$TMP/fa2" && OMP_NUM_THREADS=3 "$BIN" --fast $FQ >/dev/null 2>&1)
(cd "$TMP/fa3" && "$BIN" $FQ >/dev/null 2>&1)
printf 'fast: true\n' > "$TMP/fast.yaml"
(cd "$TMP/fa4" && OMP_NUM_THREADS=1 "$BIN" --config "$TMP/fast.yaml" $FQ >/dev/null 2>&1)
if [ -s "$TMP/fa1/Results/vcf" ] && [ "$(norm "$TMP/fa1/Results" | md5sum)" = "$(norm "$TMP/fa2/Results" | md5sum)" ]; then ok "--fast: same results with 1 and 3 threads"; else bad "--fast threads"; fi
[ "$(norm "$TMP/fa1/Results" | md5sum)" != "$(norm "$TMP/fa3/Results" | md5sum)" ] && ok "--fast changes the random numbers of the reads" || bad "--fast has no effect"
mkdir -p "$TMP/fa5" "$TMP/fa6"
(cd "$TMP/fa5" && "$BIN" --fast -n2 -s30 -l800 -u1e-6 -e10000 -b1 -C10 -Y -#31 >/dev/null 2>&1)
(cd "$TMP/fa6" && "$BIN" -n2 -s30 -l800 -u1e-6 -e10000 -b1 -C10 -Y -#31 >/dev/null 2>&1)
[ -s "$TMP/fa5/Results/true_hap" ] && cmp -s "$TMP/fa5/Results/true_hap" "$TMP/fa6/Results/true_hap" && ok "--fast keeps tree, mutations and true genotypes" || bad "--fast changed the genotypes"
[ "$(norm "$TMP/fa1/Results" | md5sum)" = "$(norm "$TMP/fa4/Results" | md5sum)" ] && ok "fast: true in a config file = --fast" || bad "fast in config"
if python3 -c "import numpy" 2>/dev/null; then
  if python3 "$ROOT/tests/fast_stats.py" "$BIN" "$TMP/fastcmp" > "$TMP/fast_stats.txt" 2>&1; then ok "--fast: DP, read counts and likelihoods have the same distributions (KS and z tests)"; else bad "--fast distributions"; cat "$TMP/fast_stats.txt" | head -30; fi
else echo "  skip (numpy not installed)"; fi

echo "2g9. ADO rate of every cell (-P) is written to ado_cells"
mkdir -p "$TMP/ac1" "$TMP/ac2"
(cd "$TMP/ac1" && "$BIN" -n2 -s300 -l50 -u1e-6 -e10000 -b1 -C5 -P0.2 0.01 -#4 >/dev/null 2>&1)
(cd "$TMP/ac2" && CELLCOAL_SPARSE=1 "$BIN" -n1 -s300 -l50 -u1e-6 -e10000 -b1 -C5 -P0.2 0.01 -#4 >/dev/null 2>&1)
for d in ac1 ac2; do
  r=$(python3 -c "
import sys
rows=[l.split() for l in open(sys.argv[1]) if l[0] != '[' and not l.startswith('cell\tado')][:301]
v=[float(x[1]) for x in rows]
print(len(v), round(sum(v)/len(v),2), all(0<x<1 for x in v))" "$TMP/$d/Results/ado_cells" 2>/dev/null)
  [ "$r" = "301 0.2 True" ] && ok "ado_cells ($d): 300 cells + outgroup, mean 0.2" || bad "ado_cells $d ($r)"
done
[ "$(grep -c '^\[#' "$TMP/ac1/Results/ado_cells")" = 2 ] && ok "ado_cells has one section per replicate" || bad "ado_cells replicates"

echo "2g10. Correlation between the coverage and the ADO rate of the cells (--cov-ado-corr)"
cat > "$TMP/covado.py" <<'PY'
import sys, numpy as np
d = sys.argv[1]
ado = [float(l.split('\t')[1]) for l in open(d + '/Results/ado_cells') if l[0] != '[' and not l.startswith('cell\tado')]
dp = None
for l in open(d + '/Results/vcf'):
    if l.startswith('#CHROM'): dp = np.zeros(len(l.split('\t')) - 9); continue
    if dp is None or l[0] in '#[' or '\t' not in l: continue
    p = l.split('\t'); i = p[8].split(':').index('DP')
    for k, c in enumerate(p[9:]):
        x = c.split(':')
        if len(x) > i and x[i] != '.': dp[k] += int(x[i])
ado = np.array(ado[:len(dp)]); r = lambda a: np.argsort(np.argsort(a))
print(round(float(np.corrcoef(r(ado), r(dp))[0, 1]), 3))
PY
if python3 -c "import numpy" 2>/dev/null; then
  CQ="-n1 -s300 -l3000 -u1e-6 -e10000 -b1 -C20 -J4 -P0.2 0.02 -#8"
  mkdir -p "$TMP/cv0" "$TMP/cv1" "$TMP/cv2"
  (cd "$TMP/cv0" && "$BIN" $CQ >/dev/null 2>&1)
  (cd "$TMP/cv1" && "$BIN" --cov-ado-corr -0.8 $CQ >/dev/null 2>&1)
  printf 'cov_ado_correlation: 0.8\n' > "$TMP/covado.yaml"
  (cd "$TMP/cv2" && "$BIN" --config "$TMP/covado.yaml" $CQ >/dev/null 2>&1)
  c0=$(python3 "$TMP/covado.py" "$TMP/cv0"); c1=$(python3 "$TMP/covado.py" "$TMP/cv1"); c2=$(python3 "$TMP/covado.py" "$TMP/cv2")
  python3 -c "import sys; sys.exit(0 if abs($c0) < 0.25 else 1)" && ok "no correlation without the option ($c0)" || bad "cov-ado baseline ($c0)"
  python3 -c "import sys; sys.exit(0 if $c1 < -0.6 else 1)" && ok "--cov-ado-corr -0.8: more dropout, less coverage ($c1)" || bad "cov-ado negative ($c1)"
  python3 -c "import sys; sys.exit(0 if $c2 > 0.5 else 1)" && ok "cov_ado_correlation: 0.8 in a config file ($c2)" || bad "cov-ado positive ($c2)"
  cmp -s "$TMP/cv0/Results/ado_cells" "$TMP/cv1/Results/ado_cells" && ok "the ADO rates themselves do not change" || bad "cov-ado changed the ADO rates"
else echo "  skip (numpy not installed)"; fi
if "$BIN" --cov-ado-corr 0.5 -s5 2>&1 | grep -q 'needs coverage variation among cells (--depth cells) and ADO variation'; then ok "--cov-ado-corr needs -J and -P"; else bad "cov-ado validation"; fi

echo "2g10c. Pileup has every read (also errors to bases that are in no genotype)"
mkdir -p "$TMP/pu"
(cd "$TMP/pu" && "$BIN" -n1 -s20 -l2000 -j200 -q0 -D0.2 -C10 -E0.02 -A0.05 0.002 0 -9 -#5 >/dev/null 2>&1)
if (cd "$TMP/pu" && python3 -c "
import sys
b = 0; n = 0; vcf = {}
for l in open('Results/vcf'):
    p = l.rstrip('\n').split('\t')
    if l[0] == '#' or len(p) < 10: continue
    k = p[8].split(':').index('RC')
    vcf[p[1]] = [list(map(int, c.split(':')[k].split(','))) for c in p[9:]]
for l in open('Results/pileup'):
    p = l.rstrip('\n').split('\t')
    if len(p) < 6: continue
    for i in range((len(p) - 3) // 3):
        d = int(p[3 + 3 * i]); s = p[4 + 3 * i]; n += 1
        rc = vcf[p[1]][i]; ref = 'ACGT'.index(p[2])
        cnt = [s.count(x) for x in 'ACGT']; cnt[ref] = s.count('.')
        if len(s) != d or cnt != rc: b += 1
sys.exit(0 if n > 0 and b == 0 else 1)
"); then ok "pileup bases = depth and = VCF read counts"; else bad "pileup reads"; fi

echo "2g10b. Amplification error that varies among cells (--amp-cell-var)"
cat > "$TMP/ampcell.py" <<'PY'
import sys, numpy as np
d = sys.argv[1]
amp = [float(l.rstrip('\n').split('\t')[-1]) for l in open(d + '/Results/cell_params') if l[0] != '[' and not l.startswith('cell\t')]
nr = dp = None
for l in open(d + '/Results/pileup'):
    p = l.rstrip('\n').split('\t')
    if len(p) < 6: continue
    if nr is None: n = (len(p) - 3) // 3; nr = np.zeros(n); dp = np.zeros(n)
    for k in range(len(nr)):
        dp[k] += int(p[3 + 3 * k]); nr[k] += sum(ch in 'ACGTacgt' for ch in p[4 + 3 * k])
m = min(len(amp), len(nr)); amp = np.array(amp[:m]); x = (nr / np.maximum(dp, 1))[:m]
print(round(float(np.std(amp)), 4), round(float(np.corrcoef(amp, x)[0, 1]), 3))
PY
if python3 -c "import numpy" 2>/dev/null; then
  AQ="-n1 -s60 -l3000 -u1e-6 -e10000 -b1 -C10 -E0.001 -A0.05 0.002 0 -9 -#4"
  mkdir -p "$TMP/am0" "$TMP/am1" "$TMP/am2"
  (cd "$TMP/am0" && "$BIN" $AQ >/dev/null 2>&1)
  (cd "$TMP/am1" && "$BIN" --amp-cell-var 0.0007 $AQ >/dev/null 2>&1)
  printf 'amplification_error_cell_variance: 0.0007\n' > "$TMP/ampcell.yaml"
  (cd "$TMP/am2" && "$BIN" --config "$TMP/ampcell.yaml" $AQ >/dev/null 2>&1)
  read s1 c1 <<< "$(python3 "$TMP/ampcell.py" "$TMP/am1")"
  python3 -c "import sys; sys.exit(0 if 0.015 < $s1 < 0.04 and $c1 > 0.9 else 1)" && ok "--amp-cell-var: per-cell means in cell_params (sd $s1) drive the error reads (r = $c1)" || bad "amp-cell-var ($s1, $c1)"
  cmp -s "$TMP/am1/Results/pileup" "$TMP/am2/Results/pileup" && ok "amplification_error_cell_variance in a config file" || bad "amp-cell-var config"
  [ ! -e "$TMP/am0/Results/cell_params" ] && ok "no cell_params without the option" || bad "cell_params without --amp-cell-var"
else echo "  skip (numpy not installed)"; fi
if "$BIN" --amp-cell-var 0.001 -s5 2>&1 | grep -q 'needs an amplification error (--amplification error)'; then ok "--amp-cell-var needs -A"; else bad "amp-cell-var validation"; fi

echo "2g11. BCF output (--bcf)"
BQ="-n2 -s30 -l1500 -u3e-6 -e10000 -b1 -C8 -V5 -E0.01 -A0.05 0.01 0 -D0.1 -B0.1 0.005 -J3 -O3 -v -#9"
mkdir -p "$TMP/bf1" "$TMP/bf2" "$TMP/bf3"
(cd "$TMP/bf1" && "$BIN" --bcf $BQ >/dev/null 2>&1)
(cd "$TMP/bf2" && "$BIN" $BQ >/dev/null 2>&1)
(cd "$TMP/bf3" && CELLCOAL_STREAM=1 CELLCOAL_STREAM_MB=1 "$BIN" --bcf $BQ >/dev/null 2>&1)
(cd "$TMP/bf3" && mkdir -p txt && cd txt && CELLCOAL_STREAM=1 CELLCOAL_STREAM_MB=1 "$BIN" $BQ >/dev/null 2>&1)
if [ -s "$TMP/bf1/Results/vcf_dir/vcf.0001.bcf" ] && [ ! -e "$TMP/bf1/Results/vcf_dir/vcf.0001" ]; then ok "--bcf writes vcf.NNNN.bcf (no VCF text left)"; else bad "--bcf files"; fi
for i in 0001 0002; do
  if python3 "$ROOT/tests/bcf_check.py" "$TMP/bf1/Results/vcf_dir/vcf.$i.bcf" "$TMP/bf2/Results/vcf_dir/vcf.$i" >"$TMP/bcfchk.txt" 2>&1; then ok "BCF replicate $i equals the VCF ($(cat "$TMP/bcfchk.txt"))"; else bad "BCF replicate $i"; head -5 "$TMP/bcfchk.txt"; fi
done
if python3 "$ROOT/tests/bcf_check.py" "$TMP/bf3/Results/vcf_dir/vcf.0001.bcf" "$TMP/bf3/txt/Results/vcf_dir/vcf.0001" >"$TMP/bcfchk.txt" 2>&1; then ok "streamed BCF equals the streamed VCF"; else bad "streamed BCF"; head -5 "$TMP/bcfchk.txt"; fi
if command -v bcftools >/dev/null 2>&1; then
  f="$TMP/bf1/Results/vcf_dir/vcf.0001.bcf"; v="$TMP/bf2/Results/vcf_dir/vcf.0001"
  n1=$(bcftools view -H "$f" 2>/dev/null | wc -l | tr -d " "); n2=$(grep -v '^#' "$v" | grep -c "$(printf '\t')")
  bcftools query -f '%POS\t%REF\t%ALT[\t%DP]\n' "$f" 2>/dev/null | md5sum > "$TMP/q1.txt"
  python3 - "$v" <<'PY' | md5sum > "$TMP/q2.txt"
import sys
for l in open(sys.argv[1]):
    if l[0] == '#' or '\t' not in l: continue
    p = l.rstrip('\n').split('\t'); i = p[8].split(':').index('DP')
    print('\t'.join([p[1], p[3], p[4]] + [c.split(':')[i] for c in p[9:]]))
PY
  if [ "$n1" = "$n2" ] && [ "$n1" -gt 100 ] && cmp -s "$TMP/q1.txt" "$TMP/q2.txt"; then ok "bcftools reads the BCF: $n1 records, positions, alleles and DP equal the VCF"; else bad "bcftools reading ($n1 vs $n2)"; bcftools view -H "$f" 2>&1 | head -3; fi
  bcftools index -f "$f" >/dev/null 2>&1 && ok "bcftools can (re)index the BCF" || bad "bcftools index"
else echo "  skip (bcftools not installed)"; fi
echo "2g11b. The BCF index (.bcf.csi) is written by the program and is correct (several bins)"
BX="-n1 -s12 -l70000 -u1e-6 -e10000 -b1 -C8 -V5 -E0.01 -v -#9"
mkdir -p "$TMP/bx1" "$TMP/bx2"
(cd "$TMP/bx1" && "$BIN" --bcf $BX >/dev/null 2>&1); (cd "$TMP/bx2" && "$BIN" $BX >/dev/null 2>&1)
if [ -s "$TMP/bx1/Results/vcf_dir/vcf.0001.bcf.csi" ] && python3 "$ROOT/tests/bcf_check.py" "$TMP/bx1/Results/vcf_dir/vcf.0001.bcf" "$TMP/bx2/Results/vcf_dir/vcf.0001" >"$TMP/bcfchk.txt" 2>&1 && grep -q "CSI ok" "$TMP/bcfchk.txt"; then ok "CSI index checked independently ($(grep 'CSI ok' "$TMP/bcfchk.txt"))"; else bad "CSI index"; head -5 "$TMP/bcfchk.txt"; fi
if command -v bcftools >/dev/null 2>&1; then
  f="$TMP/bx1/Results/vcf_dir/vcf.0001.bcf"; v="$TMP/bx2/Results/vcf_dir/vcf.0001"
  nidx=$(bcftools index -n "$f" 2>/dev/null); ntxt=$(grep -v '^#' "$v" | grep -c "$(printf '\t')")
  r1=$(bcftools view -H -r 1:20000-50000 "$f" 2>/dev/null | cut -f2 | md5sum)
  r2=$(grep -v '^#' "$v" | grep "$(printf '\t')" | awk -F'\t' '$2>=20000 && $2<=50000 {print $2}' | md5sum)
  if [ "$nidx" = "$ntxt" ] && [ "$ntxt" -gt 50 ] && [ "$r1" = "$r2" ]; then ok "bcftools uses our .csi: $nidx records, region query equals the VCF"; else bad "bcftools with our .csi (count $nidx vs $ntxt)"; bcftools view -H -r 1:20000-50000 "$f" 2>&1 | head -3; fi
fi
if "$BIN" --bcf -z -s5 2>&1 | grep -q 'cannot be combined' && "$BIN" --bcf -n3 -s5 2>&1 | grep -q 'needs --output replicate-folders'; then ok "--bcf option checks"; else bad "--bcf option checks"; fi

echo "2g12. Cache of genotype likelihoods (same counts and amplification errors) does not change the results"
GQ="-n1 -s60 -l2000 -u2e-6 -e10000 -b1 -C6 -E0.01 -D0.1 -B0.1 0.005 -1 -8 -#12"
mkdir -p "$TMP/gc1" "$TMP/gc2"
(cd "$TMP/gc1" && "$BIN" $GQ >/dev/null 2>&1)
(cd "$TMP/gc2" && CELLCOAL_NO_GL_CACHE=1 "$BIN" $GQ >/dev/null 2>&1)
if [ -s "$TMP/gc1/Results/vcf" ] && [ "$(norm "$TMP/gc1/Results" | md5sum)" = "$(norm "$TMP/gc2/Results" | md5sum)" ]; then ok "VCF and CATG identical with and without the cache"; else bad "genotype likelihood cache changes the results"; fi

echo "2g13. Streamed reads also work with the full genotype/haplotype files (-3 -4) and -5 -Y: same files as stored reads, any thread count"
C="-n2 -s30 -l3000 -u1e-6 -e10000 -g1e-4 -b1 -m0 -v -C20 -V5 -E0.001 -A0.05 0.01 0 -B0.3 0.001 -P0.1 0.01 -1 -3 -4 -5 -Y -8 -9 -#79"
mkdir -p "$TMP/fu0" "$TMP/fu1" "$TMP/fu2"
(cd "$TMP/fu0" && CELLCOAL_SPARSE=1 CELLCOAL_STREAM=0 "$BIN" $C >out.txt 2>&1); (cd "$TMP/fu1" && CELLCOAL_STREAM=1 OMP_NUM_THREADS=1 "$BIN" $C >out.txt 2>&1); (cd "$TMP/fu2" && CELLCOAL_STREAM=1 OMP_NUM_THREADS=3 "$BIN" $C >out.txt 2>&1)
if [ -n "$(find "$TMP/fu0/Results/full_genotypes_dir" -type f -size +1k)" ] && [ -n "$(find "$TMP/fu1/Results/full_haplotypes_dir" -type f -size +1k)" ] && [ -n "$(find "$TMP/fu1/Results/true_haplotypes_dir" -type f -size +1k)" ]; then
  if [ "$(norm "$TMP/fu0/Results" | md5sum)" = "$(norm "$TMP/fu1/Results" | md5sum)" ] && [ "$(norm "$TMP/fu1/Results" | md5sum)" = "$(norm "$TMP/fu2/Results" | md5sum)" ]; then ok "-3 -4 -5 -Y: stored = streamed (1 and 3 threads)"; else bad "full genotypes/haplotypes stored vs streamed"; fi
else bad "full genotype/haplotype files not produced"; fi

echo "2g14. Biological validation: coalescent, mutations, site frequency spectrum, reads and dropout agree with theory (tests/theory.py)"
if python3 "$ROOT/tests/theory.py" "$BIN" >"$TMP/theory.txt" 2>&1; then ok "$(tail -1 "$TMP/theory.txt")"; else bad "theory"; grep FAIL "$TMP/theory.txt" | head -8; fi

echo "2g15. Enormous gamma shapes (found by tests/fuzz.py) terminate; the fuzzer finds no crash on a short run"
mkdir -p "$TMP/gh"
if (cd "$TMP/gh" && python3 - "$BIN" <<'PY'
import subprocess, sys
for o in ("-a1e308", "-J1e308", "-O1e308", "-V1e308", "-i1e308"):
    subprocess.run([sys.argv[1], "-n1", "-s8", "-l200", "-u1e-6", "-e1000", "-b1", "-C5", "-V1", "-E0.001", o], capture_output=True, timeout=60, check=True)
PY
); then ok "huge gamma shapes (-a -J -O -V -i 1e308) finish"; else bad "huge gamma shape hangs or fails"; fi
if python3 "$ROOT/tests/fuzz.py" "$BIN" --runs 60 --seed 11 >"$TMP/fuzz.txt" 2>&1; then ok "$(tail -1 "$TMP/fuzz.txt")"; else bad "fuzz"; head -6 "$TMP/fuzz.txt"; fi

echo "2g16. Correlation of being a doublet with the coverage (--cov-doublet-corr) and the ADO rate (--ado-doublet-corr) of the cells"
for cfg in "0 0" "0.8 0" "-0.8 0" "0 0.8" "0.8 0.8"; do
  set -- $cfg; d="$TMP/dc_${1}_${2}"; mkdir -p "$d"
  (cd "$d" && "$BIN" -n30 -s100 -l200 -u1e-6 -e10000 -b1 -C8 -V5 -E0.001 -J4 -P0.2 0.02 -B0.3 0.02 --cov-ado-corr 0.5 --cov-doublet-corr $1 --ado-doublet-corr $2 -#9 >/dev/null 2>&1)
done
if python3 - "$TMP" <<'PY'
import math, statistics as st, sys
def load(path):
    reps = []; cur = None
    for l in open(path):
        l = l.rstrip("\n")
        if l.startswith("[#"): cur = []; reps.append(cur); continue
        f = l.split("\t")
        if f[0] != "cell": cur.append((float(f[1]), float(f[2]), int(f[3])))
    return reps
def corr(x, y):
    mx, my = st.mean(x), st.mean(y)
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / math.sqrt(sum((a - mx) ** 2 for a in x) * sum((b - my) ** 2 for b in y))
def stats(cfg):
    reps = load("%s/dc_%s/Results/cell_params" % (sys.argv[1], cfg))
    cd = [corr([c[0] for c in r], [c[2] for c in r]) for r in reps]
    ad = [corr([c[1] for c in r], [c[2] for c in r]) for r in reps]
    return st.mean(cd), st.mean(ad), sum(c[2] for c in reps[0])
z = stats("0_0"); cp = stats("0.8_0"); cn = stats("-0.8_0"); ap = stats("0_0.8"); both = stats("0.8_0.8")
assert abs(z[0]) < 0.1 and abs(z[1]) < 0.1, ("no correlation requested", z)
assert cp[0] > 0.45 and cn[0] < -0.35, ("coverage - doublet", cp, cn)
assert ap[1] > 0.45 and both[0] > 0.45 and both[1] > 0.45, ("ADO - doublet", ap, both)
assert cp[0] > z[0] + 0.4 and ap[1] > z[1] + 0.4
assert len({z[2], cp[2], cn[2], ap[2], both[2]}) == 1, ("number of doublets of replicate 1 must not change", z[2], cp[2], cn[2], ap[2], both[2])
print("r(cov,doublet) %.2f/%.2f/%.2f, r(ado,doublet) %.2f/%.2f, doublets in replicate 1: %d" % (z[0], cp[0], cn[0], z[1], ap[1], z[2]))
PY
then ok "doublets follow coverage and ADO (point-biserial r below the latent 0.8 because 30% are doublets), same number of doublets"; else bad "doublet correlations"; fi
if "$BIN" -n1 -s20 -l50 -C5 -B0.3 0.02 --cov-doublet-corr 0.5 >/dev/null 2>&1; then bad "--cov-doublet-corr without -J accepted"; else ok "--cov-doublet-corr without -J is rejected"; fi
if "$BIN" -n1 -s20 -l50 -C5 -J4 -P0.2 0.02 -B0.3 0.02 --cov-ado-corr 0.5 --cov-doublet-corr -0.8 --ado-doublet-corr 0.8 >/dev/null 2>&1; then bad "incompatible correlations accepted"; else ok "incompatible correlations (not a valid matrix) are rejected"; fi
printf '{"s":20,"l":100,"C":5,"J":4,"P":[0.2,0.02],"B":[0.3,0.02],"cov_doublet_correlation":0.6,"ado_doublet_correlation":0.3}' > "$TMP/dcc.json"
mkdir -p "$TMP/dcc" && (cd "$TMP/dcc" && "$BIN" --config "$TMP/dcc.json" -n2 -#3 >/dev/null 2>&1) && [ -s "$TMP/dcc/Results/cell_params" ] && ok "correlations from a --config file" || bad "correlations from --config"

echo "2h. Run record (run.json)"
mkdir -p "$TMP/rj"
(cd "$TMP/rj" && "$BIN" -n2 -s20 -l300 -u1e-6 -e10000 -C5 -V1 -#99 >out.txt 2>&1)
if python3 -c "import json,sys; d=json.load(open(sys.argv[1])); assert d['seed']==99 and d['replicates']==2 and d['cells']==20 and d['sites']==300 and '-#99' in d['command_line'] and d['wall_seconds']>=0" "$TMP/rj/Results/run.json" 2>/dev/null; then ok "run.json"; else bad "run.json"; fi

echo "2i. Statistical checks of the pooled reads (bulk file)"
cat > "$TMP/bulkstat.py" <<'PY'
import sys
# usage: bulkstat.py <bulk file> -> total depth, mean bulk VAF, mean true VAF, correlation of the two
rows = [l.rstrip("\n").split("\t") for l in open(sys.argv[1])]
rows = [r for r in rows if len(r) == 11 and r[0] != "site"]
dep = sum(int(r[6]) for r in rows)
b = [float(r[8]) for r in rows]; t = [float(r[10]) for r in rows]
n = len(rows); mb = sum(b)/n; mt = sum(t)/n
cov = sum((x-mb)*(y-mt) for x, y in zip(b, t)); vb = sum((x-mb)**2 for x in b); vt = sum((y-mt)**2 for y in t)
print(dep, round(mb, 5), round(mt, 5), round(cov/((vb*vt)**0.5), 4))
PY
Q="-n40 -s100 -l1500 -u2e-7 -e10000 -g0 -b1 -C5 -V1 -0 -Z -#11"
mkdir -p "$TMP/vf" "$TMP/d0" "$TMP/d1" "$TMP/d2"
(cd "$TMP/vf" && "$BIN" $Q -E0 >/dev/null 2>&1); r=($(python3 "$TMP/bulkstat.py" "$TMP/vf/Results/bulk"))
close $(python3 -c "print(${r[1]}/${r[2]})") 1.0 0.05 && ok "bulk VAF vs true VAF: mean ${r[1]} vs ${r[2]}" || bad "bulk VAF vs true VAF: mean ${r[1]} vs ${r[2]}"
close ${r[3]} 1.0 0.03 && ok "bulk VAF vs true VAF: correlation ${r[3]}" || bad "bulk VAF vs true VAF: correlation ${r[3]}"
(cd "$TMP/d0" && CELLCOAL_STREAM=0 "$BIN" $Q -E0.01 >/dev/null 2>&1); (cd "$TMP/d1" && CELLCOAL_STREAM=0 "$BIN" $Q -E0.01 -B0.3 0.001 >/dev/null 2>&1); (cd "$TMP/d2" && CELLCOAL_STREAM=1 "$BIN" $Q -E0.01 -B0.3 0.001 >/dev/null 2>&1)
r0=($(python3 "$TMP/bulkstat.py" "$TMP/d0/Results/bulk")); r1=($(python3 "$TMP/bulkstat.py" "$TMP/d1/Results/bulk")); r2=($(python3 "$TMP/bulkstat.py" "$TMP/d2/Results/bulk"))
close $(python3 -c "print(${r1[0]}/${r0[0]})") 1.3 0.07 && ok "doublets (-B0.3), stored reads: depth x$(python3 -c "print(round(${r1[0]}/${r0[0]},3))") (expected 1.3)" || bad "doublets, stored reads"
close $(python3 -c "print(${r2[0]}/${r0[0]})") 1.3 0.07 && ok "doublets (-B0.3), streamed reads: depth x$(python3 -c "print(round(${r2[0]}/${r0[0]},3))") (expected 1.3)" || bad "doublets, streamed reads"

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

echo "5. Option rates match the requested values (statistical, loose tolerances)"
cat > "$TMP/rates.py" <<'PY'
import sys
# usage: rates.py <Results dir> -> prints the fraction of true heterozygous genotypes changed by
# the observation process, the fraction of all genotypes changed, and the fraction of het with both alleles lost
d=sys.argv[1]
sg=open(d+'/snv_gen').read().split('\n')
k=int(sg[2].split()[1]); sites=list(map(int,sg[3].split()))
obs={l.split()[0]:l.split()[1:] for l in sg[4:] if l.strip()}
seq={l.split()[0]:l.split()[1] for l in open(d+'/true_hap').read().split('\n')[3:] if l.strip()}
het=chg=tot=chgall=miss=0
for c,g in obs.items():
    if c=='outgcell': continue
    for si,s in enumerate(sites):
        t=seq[c+'m'][s-1]+seq[c+'p'][s-1]; tot+=1; chgall+=(g[si]!=t)
        if t[0]!=t[1]:
            het+=1; chg+=(g[si]!=t); miss+=(g[si]=='??')
print(round(chg/het,4), round(chgall/tot,4), round(miss/het,4))
PY
M="-n1 -s30 -l3000 -u3e-6 -e10000 -g0 -b1 -1 -Y -#5"
mkdir -p "$TMP/ado" "$TMP/gen" "$TMP/fj"
(cd "$TMP/ado" && "$BIN" $M -D0.3 >/dev/null 2>&1); r=($(python3 "$TMP/rates.py" "$TMP/ado/Results"))
close ${r[0]} 0.30 0.03 && ok "fixed ADO (-D0.3): ${r[0]} changed het genotypes (expected 0.30)" || bad "fixed ADO: ${r[0]} (expected 0.30)"
close ${r[2]} 0.027 0.012 && ok "ADO both alleles: ${r[2]} (expected 0.027)" || bad "ADO both alleles: ${r[2]} (expected 0.027)"
(cd "$TMP/gen" && "$BIN" $M -G0.1 0.001 >/dev/null 2>&1); r=($(python3 "$TMP/rates.py" "$TMP/gen/Results"))
close ${r[1]} 0.10 0.01 && ok "genotyping error (-G0.1): ${r[1]} (expected 0.10)" || bad "genotyping error: ${r[1]} (expected 0.10)"

echo "2j. ADO and genotyping errors in sparse mode match the requested rates"
mkdir -p "$TMP/sado" "$TMP/sgen"
(cd "$TMP/sado" && CELLCOAL_SPARSE=1 "$BIN" $M -D0.3 >/dev/null 2>&1); r=($(python3 "$TMP/rates.py" "$TMP/sado/Results"))
close ${r[0]} 0.30 0.03 && ok "sparse fixed ADO (-D0.3): ${r[0]} (expected 0.30)" || bad "sparse fixed ADO: ${r[0]} (expected 0.30)"
(cd "$TMP/sgen" && CELLCOAL_SPARSE=1 "$BIN" $M -G0.1 0.001 >/dev/null 2>&1); r=($(python3 "$TMP/rates.py" "$TMP/sgen/Results"))
close ${r[1]} 0.10 0.01 && ok "sparse genotyping error (-G0.1): ${r[1]} (expected 0.10)" || bad "sparse genotyping error: ${r[1]} (expected 0.10)"

(cd "$TMP/fj" && "$BIN" -n20 -s20 -l1000 -u1e-6 -e10000 -g0 -b1 -j50 -y1 -#5 >out.txt 2>&1)
m=$(strings "$TMP/fj/Results/log" | grep "Mean number of mutation events" | awk '{print $NF}')
[ "$m" = "50.00" ] && ok "fixed number of mutations (-j50): $m" || bad "fixed number of mutations: $m (expected 50.00)"

echo "6. Reads and signatures"
mkdir -p "$TMP/imb" "$TMP/dbl" "$TMP/sig" "$TMP/pq"
cat > "$TMP/reads.py" <<'PY'
import sys, statistics as st, collections
mode, f = sys.argv[1], sys.argv[2]
idx={'A':0,'C':1,'G':2,'T':3}; comp={'A':'T','C':'G','G':'C','T':'A'}
names=None; fr=[]; has=collections.Counter(); tot=collections.Counter(); c=collections.Counter()
for l in open(f):
    if l.startswith('#CHROM'): names=l.rstrip('\n').split('\t')[9:]; continue
    if l.startswith('#'): continue
    x=l.rstrip('\n').split('\t'); fmt=x[8].split(':')
    if mode=='sig':
        r,a=x[3],x[4].split(',')[0]
        if r in 'AG': r,a=comp[r],comp[a]
        c[r+'>'+a]+=1; continue
    for n,cell in zip(names,x[9:]):
        v=dict(zip(fmt,cell.split(':')))
        if mode=='imb':
            tg=v['TG']
            if len(tg)==3 and tg[0]!=tg[2]:
                rc=list(map(int,v['RC'].split(','))); dp=sum(rc)
                if dp>0: fr.append(rc[idx[tg[0]]]/dp)
        else:
            tot[n]+=1; has[n]+=(v['DG'] not in ('.','','./.'))
if mode=='imb': print(round(st.mean(fr),4))
elif mode=='dbl': print(round(sum(1 for n in names if tot[n] and has[n]/tot[n]>0.5)/len(names),4))
else:
    n=sum(c.values()); print(round(c['C>T']/n,4), round(c['T>C']/n,4))
PY
(cd "$TMP/imb" && "$BIN" -n1 -s20 -l5000 -u1e-6 -e10000 -g0 -b1 -C60 -V1000 -E0 -I0.7 0.001 -v -#5 >/dev/null 2>&1)
r=$(python3 "$TMP/reads.py" imb "$TMP/imb/Results/vcf_dir/vcf.0001")
close $r 0.70 0.02 && ok "allelic imbalance (-I0.7): maternal read fraction $r" || bad "allelic imbalance: $r (expected 0.70)"
(cd "$TMP/dbl" && "$BIN" -n1 -s300 -l1500 -u1e-6 -e10000 -g0 -b1 -C30 -V1000 -E0 -B0.2 0.001 -v -#5 >/dev/null 2>&1)
r=$(python3 "$TMP/reads.py" dbl "$TMP/dbl/Results/vcf_dir/vcf.0001")
close $r 0.20 0.05 && ok "doublets (-B0.2): fraction of doublet cells $r" || bad "doublets: $r (expected 0.20)"
(cd "$TMP/sig" && "$BIN" -n1 -s20 -l500000 -u1e-8 -e10000 -g0 -b1 -S1 1 1.0 -C20 -V1000 -E0 -v -#5 >/dev/null 2>&1)
r=($(python3 "$TMP/reads.py" sig "$TMP/sig/Results/vcf_dir/vcf.0001"))
close ${r[0]} 0.675 0.05 && ok "signature 1: C>T fraction ${r[0]} (expected 0.675)" || bad "signature 1: C>T ${r[0]} (expected 0.675)"
(cd "$TMP/pq" && "$BIN" $M -P0.2 0.01 >/dev/null 2>&1); r=($(python3 "$TMP/rates.py" "$TMP/pq/Results"))
close ${r[0]} 0.20 0.06 && ok "ADO per cell (-P0.2 0.01): mean ${r[0]} (expected 0.20)" || bad "ADO per cell: ${r[0]} (expected 0.20)"

echo "2g17. --summary and --sweep"
mkdir -p "$TMP/sw" && rm -rf "$TMP/sw"/* && cd "$TMP/sw"
"$BIN" -n3 -s10 -l500 -u1e-6 -e10000 -b1 -C8 -E0.001 --summary -#3 >/dev/null 2>&1
[ "$(wc -l < Results/summary.tsv)" -eq 4 ] && head -1 Results/summary.tsv | grep -q "snv_true" && ok "--summary writes one row per replicate" || bad "--summary table"
"$BIN" -n4 -s10 -l500 -u1e-6 -e10000 -b1 -E0.001 -P0.2 0.02 --sweep "C=2,8 J=1,5" -oSW -#7 >/dev/null 2>&1
[ "$(wc -l < SW/sweep.tsv)" -eq 5 ] && [ ! -d SW/run_0001 ] && ok "--sweep: 2x2 grid, one row per combination, run folders removed" || bad "--sweep grid"
python3 - <<'PY' && ok "--sweep: common seed (same genealogy in every run) and replicates table" || bad "--sweep content"
rows=[l.rstrip("\n").split("\t") for l in open("SW/sweep.tsv")]
h=rows[0]; t=h.index("tmrca")
# C does not change the genealogy: runs with the same J (rows 1,3 and 2,4) have the same mean TMRCA
assert rows[1][t]==rows[3][t] and rows[2][t]==rows[4][t], rows
assert len(open("SW/sweep_replicates.tsv").readlines())==1+4*4
PY
"$BIN" -s8 -l200 -n2 -u1e-6 -e10000 -b1 --sweep "C=2..6..2" -oSR --sweep-keep >/dev/null 2>&1
[ "$(wc -l < SR/sweep.tsv)" -eq 4 ] && [ -d SR/run_0003 ] && ok "--sweep with a range and --sweep-keep" || bad "--sweep range"
printf 'sweep: "C=2,4"\ncells: 6\nsites: 200\nreplicates: 2\nresults_folder: SC\nseed: 5\n' > c.yaml
"$BIN" --config c.yaml >/dev/null 2>&1; [ "$(wc -l < SC/sweep.tsv)" -eq 3 ] && ok "--sweep from a --config file" || bad "--sweep config"
"$BIN" -s5 -l100 --sweep "zzz=1" >/dev/null 2>&1; r1=$?
"$BIN" -s5 -l100 --sweep "n=1,2" >/dev/null 2>&1; r2=$?
"$BIN" -s5 -l100 --sweep >/dev/null 2>&1; r3=$?
[ $r1 -ne 0 ] && [ $r2 -ne 0 ] && [ $r3 -ne 0 ] && ok "--sweep rejects unknown options, -n and a missing specification" || bad "--sweep errors ($r1 $r2 $r3)"
cd "$ROOT"

echo "2g18. --summary against the true genotypes"
mkdir -p "$TMP/tr" && rm -rf "$TMP/tr"/* && cd "$TMP/tr"
TA="-n3 -s30 -l3000 -u1e-6 -e10000 -b1 -C10 -V5 -E0.002 -P0.2 0.02 -B0.1 0.01 --summary -#4"
CELLCOAL_SPARSE=1 CELLCOAL_STREAM=0 "$BIN" $TA -oST >/dev/null 2>&1
CELLCOAL_SPARSE=1 CELLCOAL_STREAM=1 "$BIN" $TA -oSS >/dev/null 2>&1
cmp -s ST/summary.tsv SS/summary.tsv && ok "truth metrics are the same with stored and streamed reads" || bad "truth metrics stored vs streamed"
python3 - <<'PY' && ok "truth metrics are in range and recall is high" || bad "truth metrics values"
rows=[l.rstrip("\n").split("\t") for l in open("ST/summary.tsv")]
h=rows[0]
for r in rows[1:]:
    d=dict(zip(h,r))
    assert 0.9 < float(d["called_fraction"]) <= 1, d
    assert 0 < float(d["genotype_error"]) < 0.1, d
    assert 0 < float(d["snv_precision"]) <= 1 and 0.5 < float(d["snv_recall"]) <= 1, d
PY
"$BIN" -n2 -s10 -l500 -u1e-6 -e10000 -b1 --summary -oSN >/dev/null 2>&1
[ "$(sed -n 2p SN/summary.tsv | cut -f9-18)" = "$(printf 'NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA')" ] && ok "truth metrics are NA without reads" || bad "truth metrics without reads"
"$BIN" -n3 -s40 -l3000 -u1e-6 -e10000 -b1 -C8 -V5 -E0.001 -P0.4 0.05 -B0.3 0.01 --summary -oSD -#6 >/dev/null 2>&1
"$BIN" -n3 -s40 -l3000 -u1e-6 -e10000 -b1 -C8 -V5 -E0.001 --summary -oS0 -#6 >/dev/null 2>&1
python3 - <<'PY' && ok "truth metrics by genotype: allelic dropout makes heterozygous cells look homozygous; doublets raise the error" || bad "truth metrics by genotype and doublets"
def load(f):
    rows=[l.rstrip("\n").split("\t") for l in open(f)]
    return [dict(zip(rows[0],r)) for r in rows[1:]]
for d in load("SD/summary.tsv"):
    het,hom,h2h=float(d["het_error"]),float(d["hom_error"]),float(d["het_to_hom"])
    assert het > 5*hom and 0.1 < het < 0.8, d
    assert h2h <= het + 1e-9 and h2h > 0.8*het, d
    assert float(d["doublet_error"]) > float(d["singlet_error"]), d
    assert 0.15 < float(d["ado_simulated"]) < 0.55, d
for d in load("S0/summary.tsv"):
    assert float(d["het_to_hom"]) < 0.2 and d["doublet_error"] == "NA" and float(d["ado_simulated"]) == 0, d
PY
"$BIN" -n6 -s20 -l2000 -u1e-6 -e10000 -b1 -C10 -V5 -P0.1 0.01 --sweep "E=0.0005,0.02" -oSE -#5 >/dev/null 2>&1
python3 - <<'PY' && ok "--sweep: genotype error grows with the sequencing error" || bad "--sweep genotype_error"
rows=[l.rstrip("\n").split("\t") for l in open("SE/sweep.tsv")]
h=rows[0]; g=h.index("genotype_error")
assert float(rows[2][g]) > 3*float(rows[1][g]), rows
PY
"$BIN" -n2 -s8 -l300 -u1e-6 -e10000 -b1 --sweep "s=8,10" -oSX >/dev/null 2>&1
grep -q "genotype_error" SX/sweep.tsv && sed -n 2p SX/sweep.tsv | grep -q "NA" && ok "--sweep: NA columns without reads" || bad "--sweep NA"
cd "$ROOT"

echo "2g19. --sweep-jobs, --sweep-resume and the sweep plot"
mkdir -p "$TMP/sj" && rm -rf "$TMP/sj"/* && cd "$TMP/sj"
JA="-n4 -s15 -l1500 -u1e-6 -e10000 -b1 -C10 -V5 -P0.1 0.01 -#5"
"$BIN" $JA --sweep "E=0.0005,0.005,0.02" -oS1 >/dev/null 2>&1
"$BIN" $JA --sweep "E=0.0005,0.005,0.02" --sweep-jobs 2 -oS2 >/dev/null 2>&1
cmp -s S1/sweep.tsv S2/sweep.tsv && cmp -s S1/sweep_replicates.tsv S2/sweep_replicates.tsv && ok "--sweep-jobs gives the same tables as running one at a time" || bad "--sweep-jobs"
rm S2/runs/summary_0003.tsv S2/sweep.tsv
"$BIN" $JA --sweep "E=0.0005,0.005,0.02" --sweep-resume -oS2 > resume.log 2>&1
cmp -s S1/sweep.tsv S2/sweep.tsv && grep -q "run 1/3 already done" resume.log && ok "--sweep-resume runs only what is missing" || bad "--sweep-resume"
"$BIN" $JA --sweep "E=0.0005,0.01" --sweep-resume -oS2 >/dev/null 2>&1 && bad "--sweep-resume accepts another sweep" || ok "--sweep-resume rejects a folder with another sweep"
"$BIN" $JA --sweep "E=0.0005,0.02" --sweep-jobs 0 -oS3 >/dev/null 2>&1 && bad "--sweep-jobs 0 accepted" || ok "--sweep-jobs needs a number from 1 to 256"
if python3 -c "import matplotlib" 2>/dev/null; then
python3 "$ROOT/tools/cellcoal_plots.py" sweep S1 --out sw.png >/dev/null 2>&1 && [ -s sw.png ] && ok "cellcoal_plots.py sweep draws the table" || bad "sweep plot"
else echo "  skip (matplotlib not installed)"; fi
cd "$ROOT"

echo "2g20. Plots from the files of a single run (the examples of the quickstart)"
mkdir -p "$TMP/qs" && rm -rf "$TMP/qs"/* && cd "$TMP/qs"
if python3 -c "import matplotlib, numpy" 2>/dev/null; then
"$BIN" -n1 -s12 -l2000 -u1e-6 -e10000 -b1 -C10 -V5 -E0.001 -P0.2 0.02 -1 -0 -#6 -oQ >/dev/null 2>&1
python3 "$ROOT/tools/cellcoal_plots.py" heatmap Q/snv_gen --out h.png >/dev/null 2>&1 && [ -s h.png ] && ok "heatmap from the SNV genotype file of one replicate" || bad "heatmap of a single file"
python3 "$ROOT/tools/bulk_vaf.py" --vcf Q/vcf --out v.png >/dev/null 2>&1 && [ -s v.png ] && ok "bulk_vaf.py reads the VCF of one replicate" || bad "bulk_vaf of a single VCF"
python3 "$ROOT/tools/cellcoal_plots.py" all Q --outdir QP >/dev/null 2>&1 && [ -s QP/report.html ] && ok "cellcoal_plots.py all on a results folder" || bad "plots all"
else echo "  skip (matplotlib not installed)"; fi
cd "$ROOT"

echo "2g21. Decimal values split by PowerShell (-E0 .001) are joined again"
rm -rf "$TMP/psh" && mkdir -p "$TMP/psh/a" "$TMP/psh/b"
(cd "$TMP/psh/a" && "$BIN" -n2 -s10 -l500 -u1e-6 -e10000 -b1 -C8 -E0.001 -P0.2 0.02 -1 -8 -#3 >/dev/null 2>&1)
(cd "$TMP/psh/b" && "$BIN" -n2 -s10 -l500 -u1e-6 -e10000 -b1 -C8 -E0 .001 -P0 .2 0.02 -1 -8 -#3 >/dev/null 2>&1)
if [ -s "$TMP/psh/a/Results/catg" ] && [ "$(norm "$TMP/psh/a/Results" | md5sum)" = "$(norm "$TMP/psh/b/Results" | md5sum)" ]; then ok "-E0 .001 -P0 .2 = -E0.001 -P0.2"; else bad "PowerShell-split decimals"; fi
cd "$ROOT"

echo "2g21. --version and --check"
"$BIN" --version | grep -q "^CellCoal version" && ok "--version" || bad "--version"
mkdir -p "$TMP/chk" && rm -rf "$TMP/chk"/* && cd "$TMP/chk"
"$BIN" --check -n3 -s20 -l1000 -C30 | grep -q "Nothing was simulated" && [ ! -d Results ] && ok "--check validates without simulating" || bad "--check"
"$BIN" --check -n3 -s0 >/dev/null 2>&1 && bad "--check accepts invalid options" || ok "--check rejects invalid options"
cd "$ROOT"

echo "2g22. run.json records the long options, so --config run.json repeats the run"
mkdir -p "$TMP/rj" && rm -rf "$TMP/rj"/* && cd "$TMP/rj"
"$BIN" -n2 -s10 -l2000 -u1e-6 -e10000 -C10 -J2 -P0.2 0.02 -A0.055 0.003 0 --amp-cell-var 0.0007 --fast --summary --cov-ado-corr -0.5 -#3 -oA >/dev/null 2>&1
grep -q '"cov_ado_correlation": -0.5' A/run.json && grep -q '"fast": true' A/run.json && grep -q '"amplification_error_cell_variance"' A/run.json && ok "run.json lists the long options" || bad "run.json long options"
"$BIN" --config A/run.json -oB >/dev/null 2>&1
cmp -s A/summary.tsv B/summary.tsv && cmp -s A/cell_params B/cell_params && ok "--config run.json repeats the run" || bad "--config run.json"
cd "$ROOT"

echo "2g23. user tree: -s ignored (cells = leaves), branch lengths keep precision"
mkdir -p "$TMP/ut" && rm -rf "$TMP/ut"/* && cd "$TMP/ut"
printf '((c1:1e-9,c2:2.5e-8):1e-6,(c3:0.0000003,(c4:4e-10,c5:5e-10):1e-7):1e-5);\n' > tiny.tre
"$BIN" -n1 -s20 -l500 -u1e-6 -e10000 -Ttiny.tre -6 -#4 -oA >o.txt 2>&1 && ok "-T with -s different from the number of leaves" || bad "-T with -s != leaves"
grep -q "c4:4e-10" A/trees && grep -q "c1:1e-09" A/trees && ok "tiny branch lengths are written without rounding to 0" || bad "tree branch precision"
printf '(%s);\n' "$(for i in $(seq 1 100); do printf 'c%d:0.001,' $i; done | sed 's/,$//')" > star.tre
"$BIN" -n1 -s20 -l500 -u1e-6 -e10000 -Tstar.tre -#4 -oB >o2.txt 2>&1 && ok "-T with 100 leaves and default -s" || bad "-T with 100 leaves"
"$BIN" --check -n1 -s20 -l500 -Tstar.tre 2>&1 | grep -q -e "-s20 is ignored: the user tree has 100 leaves" && ok "--check reports that -s is ignored with -T" || bad "--check: note on -s ignored"
echo "2g24. --no-outgroup: tumor only"
mkdir -p "$TMP/no" && rm -rf "$TMP/no"/* && cd "$TMP/no"
"$BIN" -n1 -s8 -l3000 -u1e-6 -e10000 -g1e-4 -b1 -C10 -E0.001 -6 -1 -8 -#3 -oA >/dev/null 2>&1
"$BIN" -n1 -s8 -l3000 -u1e-6 -e10000 -g1e-4 -b1 -C10 -E0.001 -6 -1 -8 --no-outgroup -#3 -oB >/dev/null 2>&1
grep -q outgcell A/vcf && ! grep -q outgcell B/vcf && ! grep -q outgcell B/catg && ! grep -q outgcell B/trees && ! grep -q outgcell B/snv_gen && ok "no outgroup in the VCF, catg, trees and genotypes" || bad "--no-outgroup outputs"
[ "$(grep -a -m1 '^#CHROM' B/vcf | tr '\t' '\n' | wc -l)" = "17" ] && ok "VCF has the 8 cells only" || bad "VCF columns"
cp B/trees B.tre
"$BIN" -n1 -l3000 -u1e-6 -e10000 -g1e-4 -b1 -C10 -E0.001 -TB.tre --no-outgroup -#3 -oC >/dev/null 2>&1 && [ "$(grep -a -m1 '^#CHROM' C/vcf | tr '\t' '\n' | wc -l)" = "17" ] && ok "the tree written without outgroup is read back with -T" || bad "re-reading the tumor-only tree"
"$BIN" --no-outgroup -n1 -s8 -l3000 -u1e-6 -e10000 -C10 -E0.001 -#3 -oD >/dev/null 2>&1; grep -q '"no_outgroup": true' D/run.json && ok "run.json records no_outgroup" || bad "run.json no_outgroup"
printf '((((c1:0.0001,c2:0.0001):0.0001,c3:0.0002):0.0001,(c4:0.0001,c5:0.0001):0.0002):0.0001,c6:0.0004);\n((((c1:0.001,c2:0.001):0.001,c3:0.002):0.001,(c4:0.001,c5:0.001):0.002):0.001,c6:0.004);\n((((c1:0.01,c2:0.01):0.01,c3:0.02):0.01,(c4:0.01,c5:0.01):0.02):0.01,c6:0.04);\n' > three.tre
rm -rf TT; "$BIN" -n3 -s6 -l20000 -Tthree.tre -1 -#4 -oTT > tt.out 2>&1
snv=$(awk '/^\[#/{getline; print $2}' TT/snv_gen | tr '\n' ' ')
set -- $snv; [ "$#" = "3" ] && [ "$1" -lt "$2" ] && [ "$2" -lt "$3" ] && ok "-T with a file of several trees: replicate i uses tree i (SNVs $snv)" || bad "several trees in the user tree file (SNVs $snv)"
"$BIN" -n5 -s6 -l20000 -Tthree.tre -#4 -oTU > tt2.out 2>&1
grep -q "has 3 trees and there are 5 replicates" tt2.out && ok "-T: a file with several trees needs one tree per replicate" || bad "too few trees accepted"
printf '((c1:0.01,c2:0.01):0.01,c3:0.02);\n((c1:0.01,c2:0.01):0.01,(c3:0.01,c4:0.01):0.02);\n' > uneven.tre
"$BIN" -n2 -s3 -l100 -Tuneven.tre -#4 -oTV > tt3.out 2>&1
grep -q "all the trees need the same number of cells" tt3.out && ok "-T: the trees of a file have the same number of cells" || bad "trees with different leaves accepted"
cd "$ROOT"

echo "2g25. --min-alt-reads and the singleton metrics"
mkdir -p "$TMP/mar" && rm -rf "$TMP/mar"/* && cd "$TMP/mar"
C="-n1 -s30 -l20000 -u2e-7 -e10000 -g1e-4 -b1 -C10 -D0.2 -E0.001 --summary -#5"
"$BIN" $C -oA >/dev/null 2>&1
"$BIN" $C --min-alt-reads 3 -oB >/dev/null 2>&1
nA=$(grep -a -vc '^#' A/vcf); nB=$(grep -a -vc '^#' B/vcf)
[ "$nB" -lt "$nA" ] && [ "$nB" -gt 0 ] && ok "--min-alt-reads 3 keeps fewer sites ($nA -> $nB)" || bad "--min-alt-reads sites ($nA -> $nB)"
grep -q '"min_alt_reads": 3' B/run.json && ok "run.json records min_alt_reads" || bad "run.json min_alt_reads"
head -1 A/summary.tsv | grep -q "singleton_sites.singleton_recall" && ok "summary has the singleton columns" || bad "singleton columns"
"$BIN" $C --min-alt-reads 0 -oC >/dev/null 2>&1; cmp -s <(grep -a -v -i -e date -e time -e Command A/vcf) <(grep -a -v -i -e date -e time -e Command C/vcf) && ok "--min-alt-reads 0 = default" || bad "--min-alt-reads 0"
cd "$ROOT"

echo "2g26. tools/chrom.py: chromosome in blocks on one tree"
if command -v python3 >/dev/null 2>&1; then
mkdir -p "$TMP/chr" && rm -rf "$TMP/chr"/* && cd "$TMP/chr"
python3 "$ROOT/tools/chrom.py" --exe "$BIN" --blocks 3 --block-size 5000 --out ch --jobs 2 --seed 4 -- -s15 -u1e-7 -e10000 -g1e-4 -b1 -C10 -D0.2 -E0.001 --min-alt-reads 2 >chrom.log 2>&1
if [ -s ch/chrom.vcf ] && [ -s ch/chrom_summary.tsv ]; then
  nv=$(grep -a -vc '^#' ch/chrom.vcf)
  tot=$(python3 - <<'P'
r=[l.rstrip("\n").split("\t") for l in open("ch/chrom_summary.tsv")]
h=r[0]; t=r[-1]
print(sum(int(t[i]) for i,c in enumerate(h) if c.startswith("inv_")) + int(t[h.index("snv_observed")]))
P
)
  [ "$tot" = "15000" ] && ok "invariant counts + variable sites = chromosome length" || bad "chrom.py totals ($tot)"
  last=$(grep -a -v '^#' ch/chrom.vcf | awk '{print $2}' | sort -n | tail -1)
  [ "$last" -gt 10000 ] && [ "$last" -le 15000 ] && ok "VCF positions run along the whole chromosome" || bad "chrom.py positions ($last)"
else bad "chrom.py did not write its outputs"; fi
cd "$ROOT"
else echo "  skip (python3 not installed)"; fi

echo "2g27. tools/chrom.py with a reference FASTA (several lines, lower case, N)"
if command -v python3 >/dev/null 2>&1; then
mkdir -p "$TMP/chf" && rm -rf "$TMP/chf"/* && cd "$TMP/chf"
python3 - <<'P'
import random
random.seed(5)
s=''.join(random.choice('ACGT') for _ in range(12000))
s=s[:4000]+'N'*1000+s[5000:]
s=s[:8000]+s[8000:9000].lower()+s[9000:]
with open('ref.fa','w') as f:
    f.write('>chrT description\n')
    for i in range(0,len(s),60): f.write(s[i:i+60]+'\n')
open('seq.txt','w').write(s.upper())
P
python3 "$ROOT/tools/chrom.py" --exe "$BIN" --fasta ref.fa --region chrT:1001-11000 --block-size 4000 --out ch --jobs 2 --seed 4 -- -s15 -u1e-6 -e10000 -g1e-4 -b1 -C10 -D0.2 -E0.001 --min-alt-reads 2 >chrom.log 2>&1
if [ -s ch/chrom.vcf ]; then
  python3 - <<'P' && ok "FASTA region: N left out, positions and counts consistent" || bad "chrom.py FASTA checks"
s=open('seq.txt').read()
r=[l.rstrip("\n").split("\t") for l in open("ch/chrom_summary.tsv")]
h=r[0]; t=r[-1]
tot=sum(int(t[i]) for i,c in enumerate(h) if c.startswith("inv_"))+int(t[h.index("snv_observed")])
assert tot==9000, tot   # 10000 sites in the region minus the 1000 N
pos=[int(l.split("\t")[1]) for l in open("ch/chrom.vcf") if not l.startswith("#")]
assert pos and all(1001<=p<=11000 and not (4001<=p<=5000) for p in pos)
assert all(l.split("\t")[0]=="chrT" for l in open("ch/chrom.vcf") if not l.startswith("#"))
P
  nb=$(ls ch/chrom_coverage_bed/*.bed 2>/dev/null | wc -l)
  [ "$nb" = "15" ] && awk -F'\t' '$1!="chrT"||$3<=$2||NF!=5{exit 1}' ch/chrom_coverage_bed/cell0001.bed && ok "one BED per cell with the coverage of the sites not in the VCF" || bad "per-cell BED files ($nb)"
else bad "chrom.py did not write its outputs (FASTA)"; fi
cd "$ROOT"
else echo "  skip (python3 not installed)"; fi

echo "2g28. --coverage: depth per cell of the sites that are not in the VCF"
mkdir -p "$TMP/cov" && rm -rf "$TMP/cov"/* && cd "$TMP/cov"
"$BIN" -n1 -s12 -l20000 -u2e-7 -e10000 -g1e-4 -b1 -C10 -D0.2 -E0.001 --summary --min-alt-reads 2 --coverage 5000 -#5 -oA >/dev/null 2>&1
if [ -s A/coverage.tsv ]; then
  python3 - <<'P' && ok "coverage.tsv: windows, sites and mean depth are consistent with the VCF" || bad "coverage.tsv content"
r=[l.rstrip("\n").split("\t") for l in open("A/coverage.tsv")]
assert r[0][:4]==["replicate","start","end","n_sites"] and len(r[0])==4+12 and len(r)==5, r[0]
s=[l.rstrip("\n").split("\t") for l in open("A/summary.tsv")]
d=dict(zip(s[0],s[1]))
assert sum(int(x[3]) for x in r[1:])+int(d["snv_observed"])==20000
m=[float(v) for x in r[1:] for v in x[4:]]
assert 5 < sum(m)/len(m) < 11, sum(m)/len(m)
P
else bad "--coverage did not write coverage.tsv"; fi
"$BIN" -n1 -s5 -l1000 -u1e-6 -e10000 --coverage 100 -#5 -oB >B.txt 2>&1; grep -q "needs read counts" B.txt && ok "--coverage needs -C" || bad "--coverage without -C"
cd "$ROOT"

echo "2g29. --min-alt-vaf"
mkdir -p "$TMP/vaf" && rm -rf "$TMP/vaf"/* && cd "$TMP/vaf"
C="-n1 -s30 -l20000 -u2e-7 -e10000 -g1e-4 -b1 -C10 -D0.2 -E0.001 --summary -#5"
"$BIN" $C --min-alt-reads 2 -oA >/dev/null 2>&1
"$BIN" $C --min-alt-reads 2 --min-alt-vaf 0.3 -oB >/dev/null 2>&1
nA=$(grep -a -vc '^#' A/vcf); nB=$(grep -a -vc '^#' B/vcf)
[ "$nB" -lt "$nA" ] && [ "$nB" -gt 0 ] && ok "--min-alt-vaf 0.3 keeps fewer sites ($nA -> $nB)" || bad "--min-alt-vaf sites ($nA -> $nB)"
grep -q '"min_alt_vaf": 0.3' B/run.json && ok "run.json records min_alt_vaf" || bad "run.json min_alt_vaf"
"$BIN" $C --min-alt-vaf 1.5 -oC >C.txt 2>&1; grep -qi "between 0 and 1" C.txt && ok "--min-alt-vaf range checked" || bad "--min-alt-vaf range"
cd "$ROOT"

echo "2g30. --ado-block: allelic dropout in contiguous segments"
mkdir -p "$TMP/adob" && rm -rf "$TMP/adob"/* && cd "$TMP/adob"
C="-n1 -s20 -l100000 -u1e-7 -e1000 -b1 -D0.25 --ado-block 2000 --ado-seed 7"
"$BIN" $C -oA -3 -#3 >/dev/null 2>&1
"$BIN" -n1 -s20 -l50000 -u1e-7 -e1000 -b1 -D0.25 --ado-block 2000 --ado-seed 7 --ado-origin 50000 -oB -3 -#9 >/dev/null 2>&1
CELLCOAL_SPARSE=1 "$BIN" $C -oS -3 -#3 >/dev/null 2>&1
CELLCOAL_SPARSE=0 "$BIN" $C -oD -3 -#3 >/dev/null 2>&1
python3 - <<'P' && ok "--ado-block: lost fraction, contiguous segments, continuity with --ado-origin" || bad "--ado-block content"
def rd(f):
    d={}
    for l in open(f).read().split("\n")[3:]:
        if l.strip():
            t=l.split(); d[t[0]]="".join(t[1:])
    return d
a=rd("A/full_gen"); b=rd("B/full_gen")
fr=0; runs=0; n=0
for c in a:
    m=a[c][0::2]; fr+=m.count("?")/len(m); n+=1
    prev=False
    for ch in m:
        q=(ch=="?")
        if q and not prev: runs+=1
        prev=q
fr/=n
assert 0.07 < fr < 0.20, fr             # expected 1-sqrt(0.75) = 0.134
assert runs < 12*n, runs                 # about 7 segments per allele (independent sites would give ~10000)
for c in a:
    x=a[c][2*50000:]; y=b[c]
    assert [ch=="?" for ch in x[0::2]]==[ch=="?" for ch in y[0::2]]
    assert [ch=="?" for ch in x[1::2]]==[ch=="?" for ch in y[1::2]]
P
cmp -s S/full_gen D/full_gen && ok "--ado-block: sparse and dense genomes give the same data" || bad "--ado-block sparse vs dense"
"$BIN" -n1 -s5 -l1000 -u1e-6 -e1000 -b1 -C10 -P0.2 0.01 -Q0.1 0.01 --ado-block 100 -oE >E.txt 2>&1; grep -aq "cannot be combined" E.txt && ok "--ado-block cannot be combined with -Q" || bad "--ado-block with -Q"
cd "$ROOT"

echo "2g31. --signature-file: COSMIC v3.6 SBS1/SBS5 profiles"
mkdir -p "$TMP/sigf" && rm -rf "$TMP/sigf"/* && cd "$TMP/sigf"
SF="$ROOT/data/signatures/COSMIC_v3.6_SBS1_SBS5_GRCh37.txt"
python3 - <<'P'
import random
random.seed(3)
s="".join(random.choice("ACGT") for _ in range(60000))
open("g.fa","w").write(">maternal\n%s\n>paternal\n%s\n" % (s, s))
open("root.txt","w").write(s)
P
"$BIN" -n1 -s20 -l60000 -u5e-8 -e1000 -S1 1 1 --signature-file "$SF" -Ug.fa -oA -Y -#4 >/dev/null 2>&1
python3 - <<'P' && ok "SBS1 (v3.6) gives mostly C>T" || bad "SBS1 from --signature-file is not mostly C>T"
comp={"A":"T","C":"G","G":"C","T":"A"}
root=open("root.txt").read()
n=ct=0
for l in open("A/true_hap").read().split("\n")[3:]:
    t=l.split()
    if len(t)<2 or t[0].startswith("outg"): continue
    for a,b in zip(root,t[1]):
        if a!=b:
            n+=1
            if (a=="C" and b=="T") or (a=="G" and b=="A"): ct+=1
assert n>50 and ct/n>0.9, (n, ct)
P
"$BIN" -n1 -s10 -l1000 -u1e-6 -e1000 -S1 3 1 --signature-file "$SF" -oB >B.txt 2>&1; grep -aq "does not exist" B.txt && ok "signature number checked against the file" || bad "signature number not checked"
"$BIN" -n1 -s10 -l1000 -u1e-6 -e1000 -S1 1 1 --signature-file nofile.txt -oC >C.txt 2>&1; grep -aq "Can't open the signature file" C.txt && ok "missing signature file reported" || bad "missing signature file"
cd "$ROOT"

echo "2g32. --amp-signature: amplification errors with a 96-context profile (PTA artifact signature of SCAN2)"
mkdir -p "$TMP/ampsig" && rm -rf "$TMP/ampsig"/* && cd "$TMP/ampsig"
PS="$ROOT/data/signatures/PTA_artifact_SCAN2.txt"
python3 - <<'P'
import random
random.seed(5)
s="".join(random.choice("ACGT") for _ in range(100000))
open("g.fa","w").write(">maternal\n%s\n>paternal\n%s\n" % (s, s))
open("root.txt","w").write(s)
P
"$BIN" -n1 -s10 -l100000 -u1e-9 -e1000 -b1 -C30 -A0.02 0.0002 0 -Ug.fa -0 -#4 --amp-signature "$PS" -oA >/dev/null 2>&1
python3 - <<'P' && ok "amplification errors follow the PTA signature (about 68% C>T)" || bad "--amp-signature spectrum"
comp={"A":"T","C":"G","G":"C","T":"A"}
root=open("root.txt").read()
ct=n=0
for l in open("A/bulk"):
    t=l.split()
    if len(t)<8 or not t[0].isdigit(): continue
    ref=root[int(t[0])-1]
    for a,c in zip("ACGT",map(int,t[2:6])):
        if a!=ref and c:
            n+=c
            if (ref=="C" and a=="T") or (ref=="G" and a=="A"): ct+=c
assert n>1000 and 0.55 < ct/n < 0.80, (n, ct/n)
P
"$BIN" -n1 -s10 -l1000 -u1e-6 -e1000 -b1 -C10 --amp-signature "$PS" -oB >B.txt 2>&1; grep -aq "needs read counts (-C) and an amplification error (--amplification error)" B.txt && ok "--amp-signature needs -A" || bad "--amp-signature without -A"
cd "$ROOT"

echo "2g33. --gain-sd / --gain-len: autocorrelated coverage from the amplification gain"
mkdir -p "$TMP/gain" && rm -rf "$TMP/gain"/* && cd "$TMP/gain"
B="-n1 -s10 -u1e-9 -e1000 -b1 -C30 --coverage 5000 --ado-seed 9"
"$BIN" $B -l200000 --gain-sd 0.6 --gain-len 40000 -oA -#3 >/dev/null 2>&1
"$BIN" $B -l100000 --gain-sd 0.6 --gain-len 40000 --ado-origin 100000 -oB -#8 >/dev/null 2>&1
"$BIN" $B -l200000 -oC -#3 >/dev/null 2>&1
python3 - <<'P' && ok "gain: mean depth kept, windows autocorrelated, continuous across --ado-origin" || bad "--gain-sd coverage"
import math, statistics as st
def rd(f):
    r=[l.rstrip("\n").split("\t") for l in open(f)]
    return [[float(x) for x in row[4:] if x!="NA"] for row in r[1:]]
def ac(w):
    out=[]
    for c in range(len(w[0])):
        s=[x[c] for x in w]; mu=st.mean(s)
        out.append(sum((s[i]-mu)*(s[i+1]-mu) for i in range(len(s)-1))/sum((v-mu)**2 for v in s))
    return st.mean(out)
a=rd("A/coverage.tsv"); b=rd("B/coverage.tsv"); c=rd("C/coverage.tsv")
m=st.mean(v for x in a for v in x)
assert 25 < m < 36, m
assert ac(a) > 0.6 and abs(ac(c)) < 0.3, (ac(a), ac(c))
h=len(a)//2; x=[v for r in a[h:] for v in r]; y=[v for r in b for v in r]
mx,my=st.mean(x),st.mean(y)
r=sum((p-mx)*(q-my) for p,q in zip(x,y))/math.sqrt(sum((p-mx)**2 for p in x)*sum((q-my)**2 for q in y))
assert r > 0.9, r
P
"$BIN" -n1 -s5 -l1000 -u1e-6 -e1000 -b1 -C10 --gain-sd 0.5 -oD >D.txt 2>&1; grep -aq "go together" D.txt && ok "--gain-sd needs --gain-len" || bad "--gain-sd without --gain-len"
cd "$ROOT"

echo "2g34. --bam / --fastq: paired-end reads of the amplified cells (molecules, PCR duplicates)"
mkdir -p "$TMP/rd" && rm -rf "$TMP/rd"/* && cd "$TMP/rd"
"$BIN" -n1 -s3 -l60000 -u1e-7 -e1000 -b1 -C10 -E0.001 -A0.01 0.0002 0 -D0.2 --ado-block 3000 --ado-seed 3 --gain-sd 0.4 --gain-len 5000 --bam --fastq --bulk-depth 20 --dup-mean 0.3 -oA -#5 >/dev/null 2>&1
ROOT="$ROOT" python3 - <<'P' && ok "reads: sorted paired BAM+BAI, depth, duplicates, errors, FASTQ agrees with the BAM" || bad "read-level simulator"
import sys, gzip, glob, subprocess
sys.path.insert(0, __import__("os").environ["ROOT"] + "/tools")
import bamcheck as b
name, ref = b.read_fasta("A/reads_1/ref.fa")
assert len(ref) == 60000
assert open("A/reads_1/ref.fa.fai").read().split() == ["chr1", "60000", "6", "60", "61"]
_, refs, recs = b.parse_bam("A/reads_1/cell0001.bam")
assert refs[0][1] == 60000 and all(recs[i]["pos"] <= recs[i+1]["pos"] for i in range(len(recs)-1))
dup = sum(1 for r in recs if r["flag"] & 1024) / len(recs)
assert 0.18 < dup < 0.28, dup
depth = sum(len(r["bases"]) for r in recs) / 60000.0
assert 4 < depth < 16, depth
_, _, nrm = b.parse_bam("A/reads_1/normal.bam")
mm = sum(1 for r in nrm for i, c in enumerate(r["bases"]) if ref[r["pos"]+i] != c) / sum(len(r["bases"]) for r in nrm)
assert 0.0006 < mm < 0.0015, mm
mmc = sum(1 for r in recs for i, c in enumerate(r["bases"]) if ref[r["pos"]+i] != c) / sum(len(r["bases"]) for r in recs)
assert 0.005 < mmc < 0.02, mmc
n = sum(1 for _ in gzip.open("A/reads_1/cell0001.R1.fastq.gz", "rt")) // 4
assert n == len(recs) // 2, (n, len(recs))
P
"$BIN" -n1 -s5 -l1000 -u1e-6 -e1000 -b1 --bam -oB >B.txt 2>&1; grep -aq "need the mean coverage" B.txt && ok "--bam needs -C" || bad "--bam without -C"
cd "$ROOT"

echo "2g35. truth.vcf and tools/benchmark.py"
cd "$TMP/rd"
python3 "$ROOT/tools/benchmark.py" eval A/reads_1/truth.vcf A/reads_1/truth.vcf > E.txt 2>&1
grep -q "precision 1.000 recall 1.000" E.txt && ok "benchmark eval: the truth against itself gives 1/1" || bad "benchmark eval"
python3 "$ROOT/tools/benchmark.py" pileup A/reads_1 --min-alt 2 --min-vaf 0.2 > P.txt 2>&1
[ $(grep -c "cell-level" P.txt) -eq 2 ] && ok "benchmark pileup: raw and deduplicated" || bad "benchmark pileup"
cd "$ROOT"

echo "2g36. tools/chrom.py with germline SNPs of the user: REF is the reference base"
mkdir -p "$TMP/gl" && rm -rf "$TMP/gl"/* && cd "$TMP/gl"
python3 - <<'P'
import random
r = random.Random(5); s = "".join(r.choice("ACGT") for _ in range(40000))
open("g.fa", "w").write(">chrT\n" + "\n".join(s[i:i+60] for i in range(0, len(s), 60)) + "\n")
v = ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS"]
for p in sorted(r.sample(range(1, 40001), 60)):
    alt = r.choice([x for x in "ACGT" if x != s[p-1]])
    v.append("chrT\t%d\t.\t%s\t%s\t.\t.\t.\tGT\t%s" % (p, s[p-1], alt, r.choice(["0|1", "1|0", "1|1", "0/1"])))
open("g.vcf", "w").write("\n".join(v) + "\n")
P
python3 "$ROOT/tools/chrom.py" --exe "$BIN" --fasta g.fa --region chrT:1-40000 --germline-vcf g.vcf --block-size 20000 --out o --seed 3 -- -s6 -u1e-6 -e1000 -b1 -C10 > chrom.log 2>&1
python3 - <<'P' && ok "chrom.py: REF = reference base, germline sites flagged GERMLINE" || bad "chrom.py germline REF"
seq = "".join(l.strip() for l in open("g.fa") if not l.startswith(">"))
germ = {int(l.split("\t")[1]) for l in open("g.vcf") if not l.startswith("#")}
n = g = 0
for l in open("o/chrom.vcf"):
    if l.startswith("#"):
        continue
    f = l.split("\t")
    n += 1
    assert f[3] == seq[int(f[1]) - 1], (f[1], f[3])
    if int(f[1]) in germ:
        assert "GERMLINE" in f[7] and "SOMATIC" not in f[7], f[7]
        g += 1
assert n > 0 and g > 20, (n, g)
P
cd "$ROOT"

echo "2g37. --blocks / --fasta / --region / --germline-vcf in the binary"
mkdir -p "$TMP/nb" && rm -rf "$TMP/nb"/* && cd "$TMP/nb"
"$BIN" -s8 -l6000 --blocks 4 -u1e-6 -e1000 -b1 -C10 -E0.001 --min-alt-reads 2 --summary --coverage 3000 -oA -#4 >/dev/null 2>&1
python3 - <<'P' && ok "blocks: one plain VCF along the chromosome, one tree, TOTAL row, coverage in chromosome coordinates" || bad "native --blocks"
rows = [l.rstrip("\n").split("\t") for l in open("A/vcf")]
assert rows[0][0].startswith("##fileformat"), rows[0][:1]
body = [r for r in rows if not r[0].startswith("#")]
assert body and all(len(r) > 9 for r in body), "blank or short lines"
pos = [int(r[1]) for r in body]
assert pos == sorted(pos) and max(pos) > 18000 and max(pos) <= 24000, (min(pos), max(pos))
s = [l.rstrip("\n").split("\t") for l in open("A/summary.tsv")]
h = s[0]; assert s[-1][0] == "TOTAL" and len(s) == 6
t = float(s[1][h.index("tmrca")])
assert all(abs(float(r[h.index("tmrca")]) - t) < 1e-9 for r in s[1:-1]), "the tree changed between blocks"
inv = sum(int(s[-1][i]) for i, c in enumerate(h) if c.startswith("inv_")) + int(s[-1][h.index("snv_observed")])
assert inv == 24000, inv
c = [l.split("\t") for l in open("A/coverage.tsv")][1:]
assert int(c[-1][2]) == 24000 and int(c[0][1]) == 1
P
"$BIN" -s6 -l5000 --fasta "$TMP/gl/g.fa" --region chrT:1-12000 -u1e-6 -e1000 -b1 -C10 --blocks 2 -oB -#4 >/dev/null 2>&1
python3 - <<'P' && ok "--blocks limits the number of blocks of --fasta" || bad "--blocks with --fasta"
pos = [int(l.split("\t")[1]) for l in open("B/vcf") if not l.startswith("#") and l.strip()]
assert pos and max(pos) <= 10000
P
cd "$ROOT"

echo "2g38. built-in signatures = COSMIC v2 table of the repository"
mkdir -p "$TMP/v2" && rm -rf "$TMP/v2"/* && cd "$TMP/v2"
"$BIN" -n1 -s6 -l30000 -u1e-6 -e1000 -b1 -S1 1 1 -oA -#3 -3 >/dev/null 2>&1
"$BIN" -n1 -s6 -l30000 -u1e-6 -e1000 -b1 -S1 1 1 --signature-file "$ROOT/data/signatures/COSMIC_v2_SBS_GRCh37.txt" -oB -#3 -3 >/dev/null 2>&1
[ -s A/full_gen ] && cmp -s A/full_gen B/full_gen && ok "the built-in signatures give the same data as the COSMIC v2 file" || bad "built-in signatures differ from the COSMIC v2 file"
cd "$ROOT"

echo "2g39. the three modes: options that a mode does not use"
mkdir -p "$TMP/md" && rm -rf "$TMP/md"/* && cd "$TMP/md"
"$BIN" -n1 -s4 -l1000 -C10 --bam --check > a.txt 2>&1; grep -aq "mode:.*reads (paired-end" a.txt && ok "--check shows the reads mode" || bad "--check mode reads"
"$BIN" -n1 -s4 -l1000 -C10 --check > b.txt 2>&1; grep -aq "mode:.*read counts" b.txt && ok "--check shows the read counts mode" || bad "--check mode counts"
"$BIN" -n1 -s4 -l1000 --check > c.txt 2>&1; grep -aq "mode:.*genotypes" c.txt && ok "--check shows the genotypes mode" || bad "--check mode genotypes"
"$BIN" -n1 -s4 -l1000 -C10 --bam -V5 --check > d.txt 2>&1; grep -aq "reads mode" d.txt && grep -aq "is an option of the read counts" d.txt && ok "reads mode rejects -V" || bad "reads mode -V"
"$BIN" -n1 -s4 -l1000 -C10 --dup-mean 0.2 --check > e.txt 2>&1; grep -aq "needs --bam" e.txt && ok "counts mode rejects --dup-mean" || bad "counts mode --dup-mean"
"$BIN" -n1 -s4 -l1000 -E0.001 --check > f.txt 2>&1; grep -aq "ignored without -C" f.txt && ok "genotypes mode says that -E is ignored" || bad "genotypes mode -E note"
cd "$ROOT"

echo "2g40. --amp-cell-var in reads mode"
mkdir -p "$TMP/acv" && rm -rf "$TMP/acv"/* && cd "$TMP/acv"
"$BIN" -n1 -s8 -l60000 -u1e-9 -e1000 -b1 -C10 -E0.0001 -A0.01 0.0002 0 --amp-cell-var 0.00003 --bam -oA -#5 >/dev/null 2>&1
"$BIN" -n1 -s8 -l60000 -u1e-9 -e1000 -b1 -C10 -E0.0001 -A0.01 0.0002 0 --bam -oB -#5 >/dev/null 2>&1
ROOT="$ROOT" python3 - <<'P' && ok "reads: the amplification error differs among cells with --amp-cell-var" || bad "--amp-cell-var in reads mode"
import sys, statistics as st
sys.path.insert(0, __import__("os").environ["ROOT"] + "/tools")
import bamcheck as b
def rates(d):
    name, ref = b.read_fasta(d + "/reads_1/ref.fa")
    out = []
    for i in range(1, 9):
        _, _, recs = b.parse_bam("%s/reads_1/cell%04d.bam" % (d, i))
        mm = sum(1 for r in recs if not r["flag"] & 1024 for k, c in enumerate(r["bases"]) if ref[r["pos"]+k] != c)
        n = sum(len(r["bases"]) for r in recs if not r["flag"] & 1024)
        out.append(mm / n)
    return out
a, c = rates("A"), rates("B")
cv = lambda x: st.pstdev(x) / st.mean(x)
assert cv(a) > 1.5 * cv(c) and 0.003 < st.mean(a) < 0.03, (cv(a), cv(c), st.mean(a))
P
cd "$ROOT"

echo "2g41. --blocks / --fasta with reads: one BAM per cell along the chromosome"
mkdir -p "$TMP/br" && rm -rf "$TMP/br"/* && cd "$TMP/br"
python3 - <<'P'
import random
r = random.Random(9); s = [r.choice("ACGT") for _ in range(30000)]
for a, b in ((7000, 7400), (21000, 21050)):
    for i in range(a, b): s[i] = "N"
s = "".join(s)
open("g.fa", "w").write(">chrR\n" + "\n".join(s[i:i+70] for i in range(0, len(s), 70)) + "\n")
v = ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS"]
for p in sorted(r.sample([q for q in range(5001, 30001) if s[q-1] != "N"], 40)):
    alt = r.choice([x for x in "ACGT" if x != s[p-1]])
    v.append("chrR\t%d\t.\t%s\t%s\t.\t.\t.\tGT\t0|1" % (p, s[p-1], alt))
open("g.vcf", "w").write("\n".join(v) + "\n")
P
"$BIN" -n1 -s5 -l10000 --fasta g.fa --region chrR:5001-30000 --germline-vcf g.vcf -u1e-6 -e1000 -b1 -C12 -E0.001 --blocks 2 --bam --fastq -oA -#6 >log.txt 2>&1
ROOT="$ROOT" python3 - <<'P' && ok "blocks + reads: BAM in chromosome coordinates, ref.fa and truth.vcf consistent" || bad "blocks + reads"
import os, sys
sys.path.insert(0, os.environ["ROOT"] + "/tools")
import bamcheck as b
d = "A/reads_1"
seq = "".join(l.strip() for l in open("g.fa") if not l.startswith(">"))
name, ref = b.read_fasta(d + "/ref.fa")
assert name == "chrR" and len(ref) == 30000, (name, len(ref))
assert ref[5000:25000].replace("N", "") == seq[5000:25000].replace("N", "") and set(ref[:5000]) == {"N"} and set(ref[25000:]) == {"N"}
assert ref[7000:7400] == "N" * 400 and ref[:100] == "N" * 100
n = 0
for c in range(1, 6):
    text, refs, recs = b.parse_bam("%s/cell%04d.bam" % (d, c))
    assert refs[0] == ("chrR", 30000), refs
    pos = [r["pos"] for r in recs]
    assert pos == sorted(pos) and min(pos) >= 5000 and len(recs) > 800
    for r in recs:
        for k, ch in enumerate(r["bases"]):
            assert ref[r["pos"] + k] != "N", "read over a masked site"
    n += len(recs)
v = [l.split("\t") for l in open(d + "/truth.vcf") if not l.startswith("#")]
p = [int(x[1]) for x in v]
assert p == sorted(p) and min(p) > 5000 and max(p) > 15000 and max(p) <= 25000
assert all(x[3] == seq[int(x[1]) - 1] for x in v)
assert sum("GERMLINE" in x[7] for x in v) > 20
assert os.path.exists(d + "/cell0001.R1.fastq.gz") and os.path.exists(d + "/cell0001.bam.bai")
P
cd "$ROOT"

echo "2g42. variable base qualities and strand-specific amplification errors"
mkdir -p "$TMP/qs" && rm -rf "$TMP/qs"/* && cd "$TMP/qs"
"$BIN" -n1 -s2 -l60000 -u1e-9 -e1000 -b1 -C20 -E0.004 --qual-sd 0.8 --qual-slope 2 --bam --fastq -oQ -#5 >/dev/null 2>&1
"$BIN" -n1 -s2 -l60000 -u1e-9 -e1000 -b1 -C20 -E0.0002 -A0.01 0.0002 0 --dup-mean 0.001 --bam -oS0 -#5 >/dev/null 2>&1
"$BIN" -n1 -s2 -l60000 -u1e-9 -e1000 -b1 -C20 -E0.0002 -A0.01 0.0002 0 --dup-mean 0.001 --strand-bias 1 --bam -oS1 -#5 >/dev/null 2>&1
ROOT="$ROOT" python3 - <<'P' && ok "qualities: the error rate follows the quality, rises along the read; strand-specific errors are seen on one strand" || bad "qualities / strand bias"
import os, sys, collections
sys.path.insert(0, os.environ["ROOT"] + "/tools")
import bamcheck as b
_, ref = b.read_fasta("Q/reads_1/ref.fa")
_, _, recs = b.parse_bam("Q/reads_1/cell0001.bam")
n = collections.Counter(); m = collections.Counter(); cyc = collections.Counter(); cm = collections.Counter()
for r in recs:
    rev = bool(r["flag"] & 16)
    for i, c in enumerate(r["bases"]):
        q = r["qual"][i] // 4; k = ((len(r["bases"]) - 1 - i) if rev else i) // 50
        n[q] += 1; cyc[k] += 1
        if c != ref[r["pos"] + i]:
            m[q] += 1; cm[k] += 1
tot = sum(m.values()) / sum(n.values())
assert 0.003 < tot < 0.005, tot
for q in n:
    if n[q] > 20000:
        exp = 10 ** (-(4 * q + 1.5) / 10)
        assert 0.6 * exp < m[q] / n[q] < 1.6 * exp, (q, m[q] / n[q], exp)
assert cm[2] / cyc[2] > 1.5 * cm[0] / cyc[0], "the error does not rise along the read"
assert len({x for r in recs[:50] for x in r["qual"]}) > 5
assert os.path.getsize("Q/reads_1/cell0001.R1.fastq.gz") > 1000
def one_strand(d):
    _, ref = b.read_fasta(d + "/reads_1/ref.fa")
    _, _, recs = b.parse_bam(d + "/reads_1/cell0001.bam")
    st = collections.defaultdict(lambda: [0, 0])
    for r in recs:
        for i, c in enumerate(r["bases"]):
            if c != ref[r["pos"] + i]:
                st[(r["pos"] + i, c)][1 if r["flag"] & 16 else 0] += 1
    multi = [v for v in st.values() if sum(v) >= 3]
    return len(multi), sum(1 for v in multi if min(v) == 0) / len(multi)
n0, f0 = one_strand("S0"); n1, f1 = one_strand("S1")
assert n0 > 100 and n1 > 20 and f0 < 0.4 and f1 > 0.8, (n0, f0, n1, f1)
P
"$BIN" -n1 -s4 -l5000 -C10 -D0.1 --bam -oE >e.txt 2>&1; grep -aq "has to come in segments" e.txt && ok "reads mode asks for --ado-block when there is dropout" || bad "reads mode ADO without --ado-block"
cd "$ROOT"

echo "2g43. a parameter file and the same command line give the same result"
mkdir -p "$TMP/pf" && rm -rf "$TMP/pf"/* && cd "$TMP/pf"
printf '[a comment] n2 s6 l300 e1000 b1 u1e-6 h2 1000 100 40000 200 30000 20000 m2 p0.5 w2\n f0.3 0.2 0.2 0.3 t1.7 [x] N2 4 0.2 S2 11 0.3 22 0.7 G0.1 0.01 1 3 oA #7\n' > par
"$BIN" -Fpar >/dev/null 2>&1
"$BIN" -n2 -s6 -l300 -e1000 -b1 -u1e-6 -h2 1000 100 40000 200 30000 20000 -m2 -p0.5 -w2 -f0.3 0.2 0.2 0.3 -t1.7 -N2 4 0.2 -S2 11 0.3 22 0.7 -G0.1 0.01 -1 -3 -oB -#7 >/dev/null 2>&1
[ -s A/snv_gen ] && [ "$(cat A/snv_gen A/full_gen | md5sum)" = "$(cat B/snv_gen B/full_gen | md5sum)" ] && ok "parameter file = command line" || bad "parameter file = command line"
cd "$ROOT"

echo "2g44. --help lists every option and the configuration keys; example configuration files"
h=$("$BIN" --help 2>&1)
echo "$h" | grep -q "^-A:" && echo "$h" | grep -q "^-D:" && echo "$h" | grep -q "amplification_error" && ok "--help lists every option and key" || bad "--help"
"$BIN" -? >/dev/null 2>&1 && ok "-? exits with 0" || bad "-? exit status"
for f in genotypes readcounts reads; do
  mkdir -p "$TMP/ex_$f" && rm -rf "$TMP/ex_$f"/* && cd "$TMP/ex_$f"
  "$BIN" --config "$ROOT/examples/config.$f.yaml" -n1 -l1500 >/dev/null 2>&1 && [ -s Results/run.json ] && ok "examples/config.$f.yaml runs" || bad "examples/config.$f.yaml"
  cd "$ROOT"
done

echo "2g45. Grouped options (--ado, --depth, --output, ...) are equivalent to the single-letter options"
mkdir -p "$TMP/grp" && cd "$TMP/grp" && rm -rf ./*
"$BIN" -n2 -s7 -l800 -u1e-6 -e1000 -b1 -C20 -P0.1 0.001 -J0.3 -Q0.2 0.01 -1 -3 -#5 -oL >/dev/null 2>&1
"$BIN" -n2 -s7 -l800 -u1e-6 -e1000 -b1 --depth mean=20,cells=0.3 --ado cells=0.1:0.001,sites=0.2:0.01 --output snv-genotypes,full-genotypes,dir=G -#5 >/dev/null 2>&1
[ -s L/snv_gen ] && [ "$(cat L/snv_gen L/full_gen | md5sum)" = "$(cat G/snv_gen G/full_gen | md5sum)" ] && ok "group options = legacy letters" || bad "group options = legacy letters"
"$BIN" -n1 -s7 -l800 --ado nokey=1 >E.txt 2>&1; grep -aqi "unknown" E.txt && ok "unknown group key is rejected" || bad "unknown group key"
printf 'cells: 2\nsites: 800\nseed: 7\nmutation: rate=1e-6\n' > c.yaml
"$BIN" --config c.yaml -n2 -l800 -s7 -u1e-6 -oC >/dev/null 2>&1 && ok "config accepts group keys" || bad "config group key"
cd "$ROOT"

echo "2g65. Mappability (--mappability, --map-drop): MAPQ of the reads, pairs left out of the BAM"
cd "$TMP"
cat > mapq.py <<'PY'
import gzip, struct, sys, collections
f = sys.argv[1]
data = gzip.open(f, 'rb').read(); o = 4; lt = struct.unpack_from('<i', data, o)[0]; o += 4 + lt
nref = struct.unpack_from('<i', data, o)[0]; o += 4
for _ in range(nref):
    ln = struct.unpack_from('<i', data, o)[0]; o += 4 + ln + 4
c = collections.defaultdict(list)
while o < len(data):
    bs = struct.unpack_from('<i', data, o)[0]; ref, pos, lrn, mq, bn, nc, fl, ls = struct.unpack_from('<iiBBHHHi', data, o + 4)
    reg = 'low' if 10300 <= pos < 13700 else ('mid' if 30300 <= pos < 35700 else ('normal' if pos < 9700 or 14300 <= pos < 29700 or pos >= 36300 else 'edge'))
    c[reg].append(mq); o += 4 + bs
m = lambda r: sum(c[r]) / max(1, len(c[r]))
what = sys.argv[2]
if what == 'mapq':   # low 0.25 -> MAPQ 0; mid 0.75 -> MAPQ 30; elsewhere 60
    sys.exit(0 if len(c['low']) > 100 and m('low') == 0 and abs(m('mid') - 30) < 1 and abs(m('normal') - 60) < 0.01 else 1)
if what == 'drop':   # no read of the low region is left, the other regions are as before
    sys.exit(0 if len(c['low']) == 0 and len(c['mid']) > 100 and len(c['normal']) > 1000 else 1)
PY
printf "1\t10000\t14000\t0.25\n1\t30000\t36000\t0.75\n" > map.bed
if command -v python3 >/dev/null; then
  A="-n1 -s4 -#5 -l50000 -u1e-8 -e1000 -C20 -E0 --bam --fastq --mappability map.bed"
  "$BIN" $A -oMP0 >/dev/null 2>&1
  "$BIN" $A --map-drop 1 -oMP1 >/dev/null 2>&1
  python3 mapq.py MP0/reads_1/cell0001.bam mapq && ok "--mappability: MAPQ 0 where the score is 0.25, 30 at 0.75 and 60 elsewhere" || bad "mappability MAPQ"
  python3 mapq.py MP1/reads_1/cell0001.bam drop && ok "--map-drop 1: the ambiguous reads (MAPQ 0) are left out of the BAM" || bad "map-drop"
  [ "$(zcat MP0/reads_1/cell0001.R1.fastq.gz | md5sum)" = "$(zcat MP1/reads_1/cell0001.R1.fastq.gz | md5sum)" ] && ok "--map-drop does not change the FASTQ files (reads before the alignment)" || bad "map-drop changes the FASTQ"
else echo "  skip (python3 not installed)"; fi
"$BIN" -n1 -s3 -l1000 -u1e-7 -e1000 -C5 --bam --map-drop 0.5 -oME >/dev/null 2>&1 && bad "--map-drop without --mappability accepted" || ok "--map-drop needs --mappability"
printf "1\t10\t20\t1.5\n" > badmap.bed
"$BIN" -n1 -s3 -l1000 -u1e-7 -e1000 -C5 --bam --mappability badmap.bed -oME >/dev/null 2>&1 && bad "score above 1 accepted" || ok "a mappability score above 1 is rejected"
rm -rf MP0 MP1 ME
cd "$ROOT"

echo "2g66. Targeted sequencing (--targets, --off-target): depth on and off target"
cd "$TMP"
cat > tgchk.py <<'PY'
import gzip, struct, sys
data = gzip.open(sys.argv[1], 'rb').read(); o = 4; lt = struct.unpack_from('<i', data, o)[0]; o += 4 + lt
nref = struct.unpack_from('<i', data, o)[0]; o += 4
for _ in range(nref):
    ln = struct.unpack_from('<i', data, o)[0]; o += 4 + ln + 4
on = off = 0
while o < len(data):
    bs = struct.unpack_from('<i', data, o)[0]; ref, pos = struct.unpack_from('<ii', data, o + 4)
    if 10500 <= pos < 13500 or 30500 <= pos < 35500: on += 1
    elif pos < 8000 or 15000 <= pos < 29000 or pos >= 38000: off += 1
    o += 4 + bs
lon = 3000 + 5000; loff = 8000 + 14000 + 12000
r = (off / loff) / (on / lon) if on else 9
what = sys.argv[2]
if what == 'default': sys.exit(0 if on > 1000 and 0.03 < r < 0.08 else 1)
if what == 'zero': sys.exit(0 if on > 1000 and off == 0 else 1)
PY
python3 - <<'PY'
import random
r = random.Random(5); s = "".join(r.choice("ACGT") for _ in range(50000))
open("t.fa", "w").write(">chrT\n" + "\n".join(s[i:i+70] for i in range(0, len(s), 70)) + "\n")
PY
printf "chrT\t10000\t14000\nchrT\t30000\t36000\nchr2\t0\t50000\n" > targets.bed	# the lines of another chromosome are ignored
if command -v python3 >/dev/null; then
  A="-n1 -s4 -#5 -l50000 --fasta t.fa --region chrT:1-50000 -u1e-8 -e1000 -C30 -E0 --bam --targets targets.bed"
  "$BIN" $A -oTG0 >/dev/null 2>&1
  "$BIN" $A --off-target 0 -oTG1 >/dev/null 2>&1
  python3 tgchk.py TG0/reads_1/cell0001.bam default && ok "--targets: the depth off target is 5% of the depth on target (default --off-target); other chromosomes of the BED are ignored" || bad "targets default off-target"
  python3 tgchk.py TG1/reads_1/cell0001.bam zero && ok "--off-target 0: no reads away from the targets" || bad "targets off-target 0"
else echo "  skip (python3 not installed)"; fi
"$BIN" -n1 -s3 -l1000 -u1e-7 -e1000 -C5 --bam --targets targets.bed -oTE >/dev/null 2>&1 && bad "--targets without --fasta accepted" || ok "--targets needs a real chromosome (--fasta)"
"$BIN" -n1 -s3 -l1000 -u1e-7 -e1000 -C5 --bam --off-target 0.1 -oTE >/dev/null 2>&1 && bad "--off-target without --targets accepted" || ok "--off-target needs --targets"
printf "chrT\t20\t10\n" > badtg.bed
"$BIN" -n1 -s3 -l1000 --fasta t.fa --region chrT:1-5000 -u1e-7 -e1000 -C5 --bam --targets badtg.bed -oTE >/dev/null 2>&1 && bad "start after end accepted" || ok "a target with the start after the end is rejected"
rm -rf TG0 TG1 TE
cd "$ROOT"

echo "2g63. Copy number in the reads: the depth of the BAM follows the true copies of every segment"
cd "$TMP"
cat > cnchk.py <<'PY'
import gzip, struct, sys, collections
def depth(path, L):
    data = gzip.open(path, 'rb').read()
    o = 4; lt = struct.unpack_from('<i', data, o)[0]; o += 4 + lt
    nref = struct.unpack_from('<i', data, o)[0]; o += 4
    for _ in range(nref):
        ln = struct.unpack_from('<i', data, o)[0]; o += 4 + ln + 4
    dep = [0] * (L + 2)
    while o < len(data):
        bs = struct.unpack_from('<i', data, o)[0]
        ref, pos, lrn, mapq, bn, ncig, flag, lseq = struct.unpack_from('<iiBBHHHi', data, o + 4)
        co = o + 4 + 32 + lrn; p = pos
        for k in range(ncig):
            v = struct.unpack_from('<I', data, co + 4 * k)[0]; n = v >> 4; op = v & 15
            if op in (0, 7, 8):
                for q in range(p, min(p + n, L)): dep[q] += 1
                p += n
            elif op in (2, 3): p += n
        o += 4 + bs
    return dep
d, L, C = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
seg = collections.defaultdict(list)
for l in open(d + '/cnv_segments'):
    p = l.rstrip('\n').split('\t')
    if len(p) == 6 and p[1] != 'start': seg[p[0]].append((int(p[1]), int(p[2]), int(p[5])))
acc = collections.defaultdict(list)
for cell, ss in seg.items():
    dep = depth('%s/reads_1/%s.bam' % (d, cell), L)
    for a, b, cn in ss:
        if b - a >= 1500: acc[cn].append(sum(dep[a + 200:b - 200]) / (b - a - 400))
bad = 0; seen = 0
for cn, v in acc.items():
    m = sum(v) / len(v)
    if cn == 0 and m > 0.01: bad += 1
    if cn > 0 and abs(m - cn * C / 2) > 0.12 * cn * C / 2 + 1: bad += 1
    if len(v) >= 3: seen += 1
print('copy classes', sorted(acc), 'bad', bad)
sys.exit(0 if bad == 0 and seen >= 3 else 1)
PY
if command -v python3 >/dev/null; then
  "$BIN" -n1 -s12 -#3 -l60000 -u1e-8 -e1000 -C20 -E0.001 --seg-loss-rate 3e-4 --seg-gain-rate 3e-4 --seg-len 5000 --bam -ocn63 >/dev/null 2>&1
  python3 cnchk.py cn63 60000 20 && ok "BAM depth is proportional to the true number of copies (1 to 4 copies)" || bad "copy number in reads"
else echo "  skip (python3 not installed)"; fi
cd "$ROOT"

echo "2g46. Copy-number losses of segments (--seg-loss-rate, --seg-len)"
mkdir -p "$TMP/loss" && cd "$TMP/loss" && rm -rf ./*
"$BIN" -n1 -s20 -l5000 -u1e-6 -e1000 -b1 -j200 --summary -#5 >/dev/null 2>&1; a=$(cut -f6 Results/summary.tsv | tail -1)
"$BIN" -n1 -s20 -l5000 -u1e-6 -e1000 -b1 -j200 --seg-loss-rate 2e-3 --seg-len 500 --summary -#5 >/dev/null 2>&1; b=$(cut -f6 Results/summary.tsv | tail -1); lost=$(tail -1 Results/summary.tsv | awk -F'\t' '{print $NF}')
[ "$a" = 200 ] && [ "$b" -lt "$a" ] && [ $((a - b)) -eq "$lost" ] && ok "segment losses remove variable sites and the summary counts them" || bad "segment losses ($a $b $lost)"
"$BIN" -n1 -s20 -l5000 -u1e-6 -e1000 -b1 -j200 --mutation loss-rate=2e-3,segment-length=500 --summary -#5 >/dev/null 2>&1; c=$(cut -f6 Results/summary.tsv | tail -1)
[ "$b" = "$c" ] && ok "--mutation loss-rate,loss-len = --seg-loss-rate/--seg-len" || bad "group loss-rate"
"$BIN" -n1 -s20 -l5000 -u1e-6 -e1000 -b1 -C20 --seg-loss-rate 2e-3 >/dev/null 2>&1 && ok "segment losses with read counts" || bad "segment losses with -C"
cd "$ROOT"

echo "2g47. Copy-number gains of segments (--seg-gain-rate)"
mkdir -p "$TMP/gain" && cd "$TMP/gain" && rm -rf ./*
cat > meandepth.py <<'PY'
import sys
t = n = 0
for l in open(sys.argv[1] + '/pileup'):
    p = l.rstrip('\n').split('\t')
    if len(p) < 6: continue
    for k in range((len(p) - 3) // 3):
        t += int(p[3 + 3 * k]); n += 1
print(t / n if n else 0)
PY
"$BIN" -n1 -s8 -l4000 -u1e-5 -e1000 -b1 -C100 -9 -#5 -oA >/dev/null 2>&1
"$BIN" -n1 -s8 -l4000 -u1e-5 -e1000 -b1 -C100 -9 -#5 --seg-gain-rate 1e-3 --seg-len 1500 -oB >/dev/null 2>&1
if command -v python3 >/dev/null; then
  da=$(python3 meandepth.py A); db=$(python3 meandepth.py B)
  python3 -c "import sys; sys.exit(0 if float('$db') > 1.1*float('$da') else 1)" && ok "gained segments raise the read depth" || bad "gain depth ($da $db)"
else echo "  skip (python3 not installed)"; fi
grep -q "^gain" B/cnv_events && ok "cnv_events lists the gains" || bad "cnv_events"
"$BIN" -n1 -s8 -l4000 -u1e-5 -e1000 -b1 -C100 -9 -#5 --mutation gain-rate=1e-3,segment-length=1500 -oG >/dev/null 2>&1
cmp -s B/pileup G/pileup && ok "group gain-rate,segment-length = --seg-gain-rate/--seg-len" || bad "group gain-rate"
"$BIN" -n1 -s8 -l4000 -u1e-5 -e1000 -b1 -d0 --seg-gain-rate 1e-3 -1 --phase >/dev/null 2>&1 && grep -q "|" Results/snv_gen && ok "--seg-gain-rate works without -C (genotypes with copies)" || bad "gain without -C"
cd "$ROOT"

echo "2g48. A gain keeps the alleles of its moment: SNVs before it are in both copies (VAF 2/3), SNVs after it in one (VAF 1/3)"
mkdir -p "$TMP/gainvaf" && cd "$TMP/gainvaf" && rm -rf ./*
cat > vaf.py <<'PY'
import collections
ev = [l.rstrip('\n').split('\t') for l in open('Results/cnv_events') if l[:4] == 'gain']
cs = [(int(s), int(e), set(c.split(','))) for t, h, s, e, n, c in ev]
idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
names = None; hin = collections.Counter(); hout = collections.Counter()
for l in open('Results/vcf'):
    if l[0] in '[1' and '\t' not in l: continue
    p = l.rstrip('\n').split('\t')
    if l.startswith('#CHROM'): names = p[9:]; continue
    if names is None: continue
    pos = int(p[1]); fmt = p[8].split(':'); alt = idx[p[4]]
    for k, c in enumerate(names):
        rc = [int(x) for x in dict(zip(fmt, p[9 + k].split(':')))['RC'].split(',')]
        if sum(rc) < 300: continue
        v = rc[alt] / sum(rc)
        if not 0.1 < v < 0.9: continue
        n = [1 for s, e, cc in cs if s <= pos <= e and c in cc]
        if len(n) == 1: hin[round(v * 6) / 6] += 1
        elif not n: hout[round(v * 6) / 6] += 1
a, b, o = hin.get(1/3, 0), hin.get(2/3, 0), hout.get(0.5, 0)
print(1 if a > 100 and b > 100 and hin.get(0.5, 0) == 0 and o > 100 else 0)
PY
if command -v python3 >/dev/null; then
  "$BIN" -n1 -s10 -l20000 -u3e-5 -e1000 -b1 -C1000 --seg-gain-rate 4e-4 --seg-gain-len 5000 -#7 >/dev/null 2>&1
  [ "$(python3 vaf.py)" = 1 ] && ok "gained segments show VAF 1/3 and 2/3, the others 1/2" || bad "gain VAF modes"
else echo "  skip (python3 not installed)"; fi
"$BIN" -n1 -s8 -l4000 -u1e-5 -e1000 -b1 -j50 --seg-loss-rate 2e-3 --seg-loss-len 100 --summary -#5 >/dev/null 2>&1; l1=$(cut -f6 Results/summary.tsv | tail -1)
"$BIN" -n1 -s8 -l4000 -u1e-5 -e1000 -b1 -j50 --seg-loss-rate 2e-3 --seg-loss-len 2000 --summary -#5 >/dev/null 2>&1; l2=$(cut -f6 Results/summary.tsv | tail -1)
[ "$l2" -lt "$l1" ] && ok "--seg-loss-len sets the length of the losses (longer losses lose more SNVs)" || bad "--seg-loss-len ($l1 $l2)"
cd "$ROOT"

echo "2g49. Copy-number genotypes: phase is shown, the number of copies is limited by --seg-max-copies"
mkdir -p "$TMP/cng" && cd "$TMP/cng" && rm -rf ./*
cat > maxc.py <<'PY'
import sys
lim = int(sys.argv[1]); seen = 0; bad = 0; bars = 0
for l in open('Results/snv_gen'):
    p = l.split()
    if len(p) < 3 or '|' not in p[1]: continue
    for g in p[1:]:
        m, q = g.split('|')
        bars += 1
        n = len(m.replace('-', '')) + len(q.replace('-', ''))
        seen = max(seen, n)
        if n > lim: bad += 1
print(1 if bad == 0 and bars > 0 and seen == lim else 0)
PY
if command -v python3 >/dev/null; then
  "$BIN" -n1 -s10 -l2000 -u1e-6 -e1000 -d0 -#7 --seg-loss-rate 1e-3 --seg-gain-rate 4e-3 --seg-len 300 -1 --phase >/dev/null 2>&1
  [ "$(python3 maxc.py 4)" = 1 ] && ok "genotypes show maternal|paternal copies, up to 4 by default" || bad "copy genotypes (max 4)"
  "$BIN" -n1 -s10 -l2000 -u1e-6 -e1000 -d0 -#7 --seg-loss-rate 1e-3 --seg-gain-rate 4e-3 --seg-len 300 --seg-max-copies 3 -1 --phase >/dev/null 2>&1
  [ "$(python3 maxc.py 3)" = 1 ] && ok "--seg-max-copies limits the copies of a site" || bad "--seg-max-copies"
else echo "  skip (python3 not installed)"; fi
"$BIN" -n1 -s10 -l300 -u2e-5 -e3000 -d0 -#7 --seg-loss-rate 1e-3 >E.txt 2>&1; grep -aq "ISM" E.txt && ok "copy-number events keep the ISM check" || bad "ISM check with copy number"
"$BIN" -n1 -s10 -l2000 -u1e-6 -e1000 -d0 -#7 --seg-loss-rate 1e-3 --seg-gain-rate 2e-3 -1 -2 --phase >A.txt 2>&1; cp Results/snv_gen a.gen
cat > truth.py <<'PY'
# cnv_segments agrees with the genotypes and with the haplotype rows
import collections
segs = collections.defaultdict(list)
for l in open('Results/cnv_segments'):
    p = l.split()
    if len(p) == 6 and p[0] != 'cell': segs[p[0]].append((int(p[1]), int(p[2]), int(p[3]), int(p[4])))
L = open('Results/snv_gen').read().split('\n')
pos = [int(x) for x in L[3].split()]
ok = 1
for l in L[4:]:
    p = l.split()
    if len(p) < 2 or '|' not in p[1]: continue
    for site, g in zip(pos, p[1:]):
        m, q = g.split('|')
        for s, e, cm, cp in segs[p[0]]:
            if s <= site <= e:
                if (len(m.replace('-', '')) != cm) or (len(q.replace('-', '')) != cp): ok = 0
H = open('Results/snv_hap').read().split('\n')
rows = [x for x in H[4:] if x.strip()]
if int(H[2].split()[0]) != len(rows): ok = 0
print(ok)
PY
[ "$(python3 truth.py)" = 1 ] && ok "cnv_segments, genotypes and haplotype rows (one per copy) agree" || bad "cnv_segments / haplotype rows"
"$BIN" -n1 -s10 -l2000 -u1e-6 -e1000 -d0 -#7 --seg-loss-rate 1e-3 --seg-gain-rate 2e-3 -1 -2 --phase >A.txt 2>&1
cmp -s Results/snv_gen a.gen && ok "copy-number runs are reproducible with the same seed" || bad "copy-number reproducibility"
cd "$ROOT"

echo "2g50. Copy-number events work with every substitution model (finite sites, genetic signatures, mixtures)"
mkdir -p "$TMP/cnm" && cd "$TMP/cnm" && rm -rf ./*
cat > models.py <<'PY'
# the number of mutations and of variable sites with a negligible copy-number rate (the new engine) must agree with the usual simulation
import subprocess, os, shutil, statistics as st, sys
B = sys.argv[1]
def run(tag, extra, cn, seed):
    d = tag; shutil.rmtree(d, ignore_errors=True); os.makedirs(d)
    a = [B, '-n300', '-s6', '-l400', '-e1000', '-d0', '-#%d' % seed, '--summary'] + extra
    if cn: a += ['--seg-loss-rate', '1e-13']
    subprocess.run(a, cwd=d, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    rows = [l.split('\t') for l in open(d + '/Results/summary.tsv').read().split('\n')[1:] if l and l[0].isdigit()]
    return [int(r[2]) for r in rows], [int(r[5]) for r in rows]
cases = {'JC': ['-u2e-6', '-m2', '-p1'], 'HKY': ['-u2e-6', '-m2', '-p1', '-t1.7', '-f0.3', '0.2', '0.2', '0.3'], 'Mk2': ['-u2e-6', '-b0', '-m1', '-p1'],
         'mixed': ['-u2e-6', '-m2', '-p0.5', '-w4'], 'gamma': ['-u2e-6', '-m2', '-p1', '-a0.5'], 'signatures': ['-u2e-7', '-S2', '11', '0.3', '22', '0.7']}
bad = 0
for k, e in cases.items():
    a = run(k + 'L', e, False, 11); b = run(k + 'C', e, True, 12)
    for i in (0, 1):
        se = ((st.pvariance(a[i]) + st.pvariance(b[i])) / len(a[i])) ** 0.5
        if len(a[i]) < 250 or len(b[i]) < 250 or abs(st.mean(a[i]) - st.mean(b[i])) > 4 * se: bad += 1
print(bad)
PY
if command -v python3 >/dev/null; then
  [ "$(python3 models.py "$BIN")" = 0 ] && ok "mutations and variable sites agree with the usual simulation (JC, HKY, Mk2, mixed, gamma, signatures)" || bad "copy-number models"
else echo "  skip (python3 not installed)"; fi
"$BIN" -n1 -s8 -l500 -u1e-6 -e1000 -d0 -#7 -m2 -p0.5 -w3 --seg-loss-rate 1e-3 --seg-gain-rate 2e-3 --seg-len 100 -1 --phase >E.txt 2>&1 && grep -q "|" Results/snv_gen && ok "finite sites with losses and gains run and show the phase" || bad "finite sites run"
cd "$ROOT"

echo "2g51. Long options accept the values that their range allows (0, negative slope)"
cd "$TMP"
okopt=1
for o in "--qual-slope -1" "--qual-slope 0" "--qual-sd 0" "--strand-bias 0" "--dup-mean 0" "--seg-loss-rate 0" "--seg-gain-rate 0"; do
  "$BIN" -n1 -s4 -l1000 -u1e-7 -e10000 -C5 -E0.001 -D0.1 --ado-block 300 --bam $o -#1 -oZeroOpts >Z.txt 2>&1 || { okopt=0; echo "    rejected: $o"; }
  rm -rf ZeroOpts
done
[ $okopt = 1 ] && ok "0 and a negative slope are accepted" || bad "long options reject valid values"
for o in "--qual-slope 9" "--strand-bias -1" "--read-len 0"; do
  "$BIN" -n1 -s4 -l1000 -u1e-7 -e10000 -C5 -E0.001 -D0.1 --ado-block 300 --bam $o -#1 -oZeroOpts >Z.txt 2>&1 && { okopt=0; echo "    accepted: $o"; }
  rm -rf ZeroOpts
done
[ $okopt = 1 ] && ok "values out of range are still rejected" || bad "out of range values accepted"
cd "$ROOT"

echo "2g52. --compress-level changes the size of the compressed files, not their content"
cd "$TMP"
A="-n1 -s3 -l20000 -u1e-7 -e10000 -C10 -E0.001 -D0.1 --ado-block 2000 --fastq --bam -#3"
"$BIN" $A -oCL6 >/dev/null 2>&1; "$BIN" $A --compress-level 1 -oCL1 >/dev/null 2>&1; "$BIN" $A --compress-level 9 -oCL9 >/dev/null 2>&1
h() { gzip -dc < "$1" | md5sum; }
if [ "$(h CL6/reads_1/cell0001.R1.fastq.gz)" = "$(h CL1/reads_1/cell0001.R1.fastq.gz)" ] && [ "$(h CL6/reads_1/cell0001.R1.fastq.gz)" = "$(h CL9/reads_1/cell0001.R1.fastq.gz)" ] \
   && [ "$(h CL6/reads_1/cell0001.bam)" = "$(h CL1/reads_1/cell0001.bam)" ] && [ "$(h CL6/reads_1/cell0001.bam)" = "$(h CL9/reads_1/cell0001.bam)" ] \
   && [ "$(stat -c %s CL1/reads_1/cell0001.bam)" -gt "$(stat -c %s CL9/reads_1/cell0001.bam)" ]; then ok "levels 1, 6 and 9 give the same BAM and FASTQ content (level 1 is larger)"; else bad "--compress-level"; echo "    fastq md5: $(h CL6/reads_1/cell0001.R1.fastq.gz) $(h CL1/reads_1/cell0001.R1.fastq.gz) $(h CL9/reads_1/cell0001.R1.fastq.gz)"; echo "    bam md5: $(h CL6/reads_1/cell0001.bam) $(h CL1/reads_1/cell0001.bam) $(h CL9/reads_1/cell0001.bam)"; echo "    sizes: $(stat -c %s CL6/reads_1/cell0001.bam) $(stat -c %s CL1/reads_1/cell0001.bam) $(stat -c %s CL9/reads_1/cell0001.bam)"; fi
"$BIN" $A --compress-level 10 -oCL0 >/dev/null 2>&1 && bad "--compress-level 10 accepted" || ok "--compress-level out of range is rejected"
rm -rf CL6 CL1 CL9 CL0
cd "$ROOT"

echo "2g53. The reads mode is parallel by cells and gives the same files for any number of threads"
cd "$TMP"
A="-n1 -s6 -l20000 -u1e-7 -e10000 -C10 -E0.001 -D0.1 --ado-block 2000 --gain-sd 0.4 --gain-len 3000 -A0.01 0.0002 0 --qual-sd 0.5 --strand-bias 0.5 --fastq --bam --bulk-depth 5 -#4"
OMP_NUM_THREADS=1 "$BIN" $A -oTH1 >/dev/null 2>&1; OMP_NUM_THREADS=3 "$BIN" $A -oTH3 >/dev/null 2>&1
same=1
for f in cell0001.bam cell0006.bam normal.bam cell0002.R1.fastq.gz cell0005.R2.fastq.gz truth.vcf; do cmp -s TH1/reads_1/$f TH3/reads_1/$f || { same=0; echo "    differs: $f"; }; done
[ $same = 1 ] && ok "1 and 3 threads give identical BAM, FASTQ and truth files" || bad "reads mode depends on the number of threads"
rm -rf TH1 TH3
cd "$ROOT"

echo "2g55. Doublets in reads mode: the BAM of a doublet has the reads of its donor, any number of threads gives the same files"
cd "$TMP"
A="-n1 -s9 -l20000 -u1e-7 -e10000 -C10 -E0.001 -D0.1 --ado-block 2000 -B0.4 0.01 --bam --fastq -#6"
OMP_NUM_THREADS=1 "$BIN" $A -oDB1 >/dev/null 2>&1; OMP_NUM_THREADS=3 "$BIN" $A -oDB3 >/dev/null 2>&1
same=1
for f in cell0001.bam cell0004.bam cell0005.bam cell0003.R1.fastq.gz doublets.tsv; do cmp -s DB1/reads_1/$f DB3/reads_1/$f || { same=0; echo "    differs: $f"; }; done
[ $same = 1 ] && ok "1 and 3 threads give identical files with doublets" || bad "doublets in reads mode depend on the number of threads"
nd=$(tail -n +2 DB1/reads_1/doublets.tsv 2>/dev/null | wc -l)
[ "$nd" -gt 0 ] && ok "doublets.tsv lists $nd doublets" || bad "doublets.tsv has no doublets"
if command -v python3 >/dev/null && [ "$nd" -gt 0 ]; then
  python3 - DB1/reads_1 "$ROOT" <<'PYEOF' && ok "doublets have about twice the depth of the other cells, and sorted BAM" || bad "depth or order of the BAM of a doublet"
import sys, subprocess, re, glob
dr, root = sys.argv[1:3]
dbl = set(l.split("\t")[0] for l in open(dr + "/doublets.tsv").read().splitlines()[1:])
ds, ss = [], []
for f in sorted(glob.glob(dr + "/cell*.bam")):
    o = subprocess.run(["python3", root + "/tools/bamcheck.py", dr + "/ref.fa", f], capture_output=True, text=True).stdout.strip().splitlines()[-1]
    if "sorted True" not in o: sys.exit(1)
    (ds if f.split("/")[-1][:-4] in dbl else ss).append(float(re.search(r"depth ([0-9.]+)X", o).group(1)))
m1, m0 = sum(ds) / len(ds), sum(ss) / len(ss)
sys.exit(0 if 1.5 < m1 / m0 < 2.6 else 1)
PYEOF
fi
rm -rf DB1 DB3
cd "$ROOT"

echo "2g56. --region-rates: the density of SNVs follows the multipliers of the regions, only with a chromosome with coordinates"
cd "$TMP"
printf "# test\nchr1 0 5000 4\nchr1 5000 10000 0.25\n" > rr.bed
for m in rr plain; do
  : > rr_$m.txt
  for sd in 1 2 3 4; do
    X=""; [ $m = rr ] && X="--region-rates rr.bed"
    rm -rf RR; "$BIN" -n1 -s10 -#$sd --blocks 4 -l5000 -u2e-7 -e10000 -1 $X -oRR >/dev/null 2>&1
    awk 'NF==2 && $1==11 {print $2}' RR/snv_gen | tr '\n' ' ' >> rr_$m.txt; echo >> rr_$m.txt
  done
done
python3 - <<'P' && ok "blocks with multipliers 4, 0.25, 1, 1 give 2.56, 0.16, 0.64, 0.64 times the SNVs of the plain run (the total does not change)" || bad "regional mutation rates"
def tot(f):
    t = [0, 0, 0, 0]
    for l in open(f):
        v = l.split()
        for i in range(4): t[i] += int(v[i])
    return t
a, b = tot("rr_rr.txt"), tot("rr_plain.txt")
exp = [2.56, 0.16, 0.64, 0.64]
for i in range(4):
    r = a[i] / b[i]
    assert abs(r - exp[i]) < 0.25 * exp[i] + 0.05, (i, r, exp[i])
assert abs(sum(a) / sum(b) - 1) < 0.12, (sum(a), sum(b))
P
"$BIN" -n1 -s10 -l5000 -u2e-7 -e10000 -1 --region-rates rr.bed -oRR2 > rr2.out 2>&1
grep -q "needs a chromosome with coordinates" rr2.out && ok "without --blocks or --fasta the option is an error that explains why" || bad "--region-rates without coordinates was accepted"
printf "chr1 100 50 2\n" > bad.bed
"$BIN" -n1 -s10 -#1 --blocks 2 -l5000 -u2e-7 -e10000 -1 --region-rates bad.bed -oRR3 > rr3.out 2>&1
grep -q "line 1" rr3.out && ok "a wrong line of the file is reported with its number" || bad "bad region file not reported"
rm -rf RR RR2 RR3
cd "$ROOT"

echo "2g57. --seq-context: errors concentrate in the contexts of the file, -E stays the mean error, any number of threads gives the same files"
cd "$TMP"
printf "# a C after GG: 8 times the error\nGGC 8\n" > ctx.txt
A="-n1 -s6 -#3 -l30000 -u1e-7 -e10000 -C20 -E0.004 --bam"
OMP_NUM_THREADS=1 "$BIN" $A --seq-context ctx.txt -oSC1 >/dev/null 2>&1; OMP_NUM_THREADS=3 "$BIN" $A --seq-context ctx.txt -oSC3 >/dev/null 2>&1; "$BIN" $A -oSC0 >/dev/null 2>&1
cmp -s SC1/reads_1/cell0002.bam SC3/reads_1/cell0002.bam && ok "1 and 3 threads give identical BAM files" || bad "--seq-context depends on the number of threads"
if command -v python3 >/dev/null; then
python3 - "$ROOT" <<'P' && ok "GGC has several times the error of the other contexts (both reads of the pair) and the mean error is the same as without the option" || bad "context-dependent sequencing errors"
import sys, glob, collections
sys.path.insert(0, sys.argv[1] + "/tools")
import bamcheck as b
comp = {'A':'T','C':'G','G':'C','T':'A'}
def stats(d):
    name, ref = b.read_fasta(d + "/reads_1/ref.fa")
    tot = collections.Counter(); err = collections.Counter()
    for f in sorted(glob.glob(d + "/reads_1/cell*.bam")):
        t, refs, recs = b.parse_bam(f)
        for r in recs:
            if r['flag'] & 1024: continue
            rev = bool(r['flag'] & 16); p = r['pos']
            for i in range(len(r['bases'])):
                if p + i - 2 < 0: continue
                g = ref[p+i-2:p+i+1] if not rev else ref[p+i:p+i+3]
                if len(g) < 3 or 'N' in g: continue
                c = g if not rev else ''.join(comp[x] for x in reversed(g))
                k = 'GGC' if c == 'GGC' else 'x'
                tot[k] += 1
                if r['bases'][i] != ref[p+i]: err[k] += 1
    return err['GGC'] / tot['GGC'], err['x'] / tot['x'], (err['GGC'] + err['x']) / (tot['GGC'] + tot['x'])
g1, o1, m1 = stats("SC1")
g0, o0, m0 = stats("SC0")
assert g1 > 3.5 * o1, (g1, o1)
assert abs(m1 / m0 - 1) < 0.15, (m1, m0)
assert abs(g0 / o0 - 1) < 0.3, (g0, o0)
P
else echo "  skip (python3 not installed)"; fi
"$BIN" -n1 -s6 -#3 -l30000 -u1e-7 -e10000 -C20 --bam --seq-context ctx.txt -oSC4 > sc4.out 2>&1
grep -q "needs a sequencing error" sc4.out && ok "--seq-context without -E is an error" || bad "--seq-context accepted without -E"
printf "GGX 2\n" > badctx.txt
"$BIN" -n1 -s6 -#3 -l30000 -u1e-7 -e10000 -C20 -E0.004 --bam --seq-context badctx.txt -oSC5 > sc5.out 2>&1
grep -q "line 1" sc5.out && ok "a wrong line of the context file is reported with its number" || bad "bad context file not reported"
rm -rf SC0 SC1 SC3 SC4 SC5
cd "$ROOT"

echo "2g58. --damage: G>T in the forward reads and C>A in the reverse reads (orientation bias), same files for any number of threads"
cd "$TMP"
A="-n1 -s6 -#3 -l30000 -u1e-7 -e10000 -C20 -E0.0005 --bam"
OMP_NUM_THREADS=1 "$BIN" $A --damage "G>T:0.03" -oDM1 >/dev/null 2>&1; OMP_NUM_THREADS=3 "$BIN" $A --damage "G>T:0.03" -oDM3 >/dev/null 2>&1
cmp -s DM1/reads_1/cell0002.bam DM3/reads_1/cell0002.bam && ok "1 and 3 threads give identical BAM files" || bad "--damage depends on the number of threads"
if command -v python3 >/dev/null; then
python3 - "$ROOT" <<'P' && ok "G>T is about 0.03 in the forward reads and C>A in the reverse reads, and the other substitutions stay at the baseline" || bad "orientation bias of --damage"
import sys, glob, collections
sys.path.insert(0, sys.argv[1] + "/tools")
import bamcheck as b
name, ref = b.read_fasta("DM1/reads_1/ref.fa")
c = collections.Counter(); t = collections.Counter()
for f in sorted(glob.glob("DM1/reads_1/cell*.bam")):
    _, _, recs = b.parse_bam(f)
    for r in recs:
        o = 'R' if r['flag'] & 16 else 'F'
        for i, x in enumerate(r['bases']):
            g = ref[r['pos'] + i]
            if g == 'N': continue
            t[o, g] += 1
            if x != g: c[o, g, x] += 1
r = lambda o, g, x: c[o, g, x] / t[o, g]
assert 0.022 < r('F', 'G', 'T') < 0.038, r('F', 'G', 'T')
assert 0.022 < r('R', 'C', 'A') < 0.038, r('R', 'C', 'A')
assert r('R', 'G', 'T') < 0.003 and r('F', 'C', 'A') < 0.003 and r('F', 'A', 'C') < 0.003
P
else echo "  skip (python3 not installed)"; fi
"$BIN" $A --damage "G>G:0.1" -oDM4 > dm4.out 2>&1
grep -q "different ones" dm4.out && ok "a wrong --damage item is reported" || bad "bad --damage item accepted"
rm -rf DM1 DM3 DM4
cd "$ROOT"

echo "2g59. Reads in tiles: the BAM files do not depend on the size of the tile, the reads of a FASTQ are the same"
cd "$TMP"
A="-n1 -s6 -#4 -l60000 -u1e-7 -e10000 -C15 -E0.002 -A0.01 0.0001 0 --damage G>T:0.01 -B0.3 0.01 --bam --fastq"
CELLCOAL_READ_TILE=1000 "$BIN" $A -oTL1 >/dev/null 2>&1; CELLCOAL_READ_TILE=262144 "$BIN" $A -oTL2 >/dev/null 2>&1; CELLCOAL_READ_TILE=1000000000 "$BIN" $A -oTL3 >/dev/null 2>&1
same=1
for f in cell0001.bam cell0003.bam cell0006.bam cell0003.bam.bai; do cmp -s TL1/reads_1/$f TL2/reads_1/$f && cmp -s TL1/reads_1/$f TL3/reads_1/$f || { same=0; echo "    differs: $f"; }; done
[ $same = 1 ] && ok "the BAM files (with doublets, amplification errors and damage) are the same for tiles of 1000, 262144 and 1e9 bases" || bad "the BAM files depend on the size of the tile"
a=$(zcat TL1/reads_1/cell0002.R1.fastq.gz | paste - - - - | sort | md5sum); b=$(zcat TL2/reads_1/cell0002.R1.fastq.gz | paste - - - - | sort | md5sum)
[ "$a" = "$b" ] && ok "the FASTQ files have the same reads (in another order)" || bad "the reads of the FASTQ depend on the size of the tile"
rm -rf TL1 TL2 TL3
cd "$ROOT"

echo "2g60. --germline-pop: genotypes of an individual sampled from allele frequencies (Hardy-Weinberg)"
cd "$TMP"
python3 - <<'P'
import random
random.seed(7)
seq = ''.join(random.choice('ACGT') for _ in range(24000))
open('gp.fa', 'w').write('>chrT\n' + '\n'.join(seq[i:i+60] for i in range(0, len(seq), 60)) + '\n')
out = ['##fileformat=VCFv4.2', '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO']
for p in range(5, len(seq) - 5, 12):
    r = seq[p-1]; a = random.choice([x for x in 'ACGT' if x != r])
    out.append("chrT\t%d\t.\t%s\t%s\t.\tPASS\tDP=10;AF=%s;AC=3" % (p, r, a, '0.3' if (p // 12) % 2 == 0 else '0.02'))
out.append("chrT\t100\t.\tA\tC,G\t.\tPASS\tAF=0.2,0.1")
open('gp.vcf', 'w').write('\n'.join(out) + '\n')
P
: > gp_gt.txt
for sd in 1 2 3 4 5 6 7 8; do
  rm -rf GP; "$BIN" -n1 -s4 -#$sd -l6000 --fasta gp.fa --blocks 4 -u1e-9 -e10000 -C10 -E0.001 --germline-pop gp.vcf --bam -oGP >/dev/null 2>&1
  grep -v "^#" GP/reads_1/truth.vcf | grep GERMLINE | awk '{split($10,g,":"); print int((($2-1)/12)%2), g[1]}' >> gp_gt.txt
done
python3 - <<'P' && ok "the genotypes follow Hardy-Weinberg for AF 0.3 and 0.02, the heterozygous alleles are in both haplotypes" || bad "germline-pop genotypes"
import collections
c = collections.Counter(tuple(l.split()) for l in open('gp_gt.txt'))
tot = 8 * 1000
hi = sum(v for (g, x), v in c.items() if g == '0'); lo = sum(v for (g, x), v in c.items() if g == '1')
assert abs(hi / tot - 0.51) < 0.04, hi / tot            # 1 - 0.7^2
assert abs(lo / tot - 0.0396) < 0.012, lo / tot          # 1 - 0.98^2
het = c[('0', '0|1')] + c[('0', '1|0')]
assert abs(het / hi - 0.824) < 0.05, het / hi
assert 0.4 < c[('0', '0|1')] / het < 0.6
assert abs(c[('0', '1|1')] / hi - 0.176) < 0.05
P
"$BIN" -n1 -s4 -#1 -l6000 --fasta gp.fa --blocks 4 -c1e-3 --germline-pop gp.vcf -C10 --bam -oGPX > gpx.out 2>&1
grep -q "cannot be combined" gpx.out && ok "--germline-pop with -c is an error" || bad "--germline-pop with -c accepted"
"$BIN" -n1 -s4 -#1 -l6000 --blocks 4 --germline-pop gp.vcf -C10 --bam -oGPY > gpy.out 2>&1
grep -q "need --fasta" gpy.out && ok "--germline-pop without --fasta is an error" || bad "--germline-pop without --fasta accepted"
rm -rf GP GPX GPY
cd "$ROOT"

echo "2g61. Structured coalescent (--demes): tumor regions founded by colonization, without migration"
cd "$TMP"
rm -rf DT; "$BIN" -n2000 -#3 -l10 -u1e-9 -e10000 -7 --demes n=3:3,size=1:0.5,time=1,founders=1 -oDT >/dev/null 2>&1
python3 - <<'P' && ok "primary tumor and a monoclonal metastasis (size 0.5, colonized at 1): mean coalescence times 1.0 (primary), 0.43 (metastasis) and 2.0 (between)" || bad "structured coalescent times"
import re, itertools
txt = open('DT/times').read().split('-------------------- Nodes -------------------')[1:]
deme = [0, 0, 0, 1, 1, 1]
s = {}; c = {}
for blk in txt:
    node = {}
    for line in blk.splitlines():
        m = re.match(r'\s*(tip|internal|ingroupMRCA)\s+\d+\s+(\d+)\s+\(\s*(-?\d+)\s+(-?\d+)\s+(-?\d+)\)\s+\|\s+([-\d.]+)', line)
        if m: node[int(m.group(2))] = (int(m.group(5)), float(m.group(6)))
    def path(i):
        p = [i]
        while node[p[-1]][0] in node and node[p[-1]][0] != -1:
            p.append(node[p[-1]][0])
        return p
    for i, j in itertools.combinations(range(1, 7), 2):
        pj = set(path(j))
        t = next(node[x][1] for x in path(i) if x in pj)
        k = (deme[i-1], deme[j-1])
        s[k] = s.get(k, 0) + t / 10000; c[k] = c.get(k, 0) + 1
assert abs(s[(0, 0)] / c[(0, 0)] - 1.0) < 0.1, s[(0, 0)] / c[(0, 0)]
assert abs(s[(1, 1)] / c[(1, 1)] - 0.432) < 0.06, s[(1, 1)] / c[(1, 1)]
assert abs(s[(0, 1)] / c[(0, 1)] - 2.0) < 0.15, s[(0, 1)] / c[(0, 1)]
P
[ "$(cut -f2 DT/demes.tsv | tr -d '\r' | tail -n +2 | tr '\n' ' ')" = "1 1 1 2 2 2 " ] && ok "demes.tsv lists the deme of every cell" || bad "demes.tsv"
"$BIN" -n1 -#1 -l100 -u1e-6 --deme-n 3:3 > dm1.out 2>&1
grep -q "deme-time is needed" dm1.out && ok "--deme-n without colonization times is an error" || bad "missing --deme-time accepted"
"$BIN" -n1 -#1 -l100 -u1e-6 --deme-n 3:3:3 --deme-time 1:2 --deme-source 1:3 > dm2.out 2>&1
grep -q "source has to be older\|has to be another deme" dm2.out && ok "a source that does not exist yet is an error" || bad "bad --deme-source accepted"
"$BIN" -n1 -#1 -l100 -u1e-6 -g1e-4 --deme-n 3:3 --deme-time 1 > dm3.out 2>&1
grep -q "cannot be combined" dm3.out && ok "--demes with growth is an error" || bad "--demes with growth accepted"
rm -rf DT
cd "$ROOT"

echo "2g67. Clones of different fitness (--clones): growth rate by clone, sample at random from the tumor"
cd "$TMP"
python3 "$ROOT/tests/clones_check.py" "$BIN" 3000 > clones_check.out 2>&1 && ok "--clones: the times to coalescence of pairs of cells agree with a forward simulation of the whole tumor, and the clone frequencies are those of the sizes" || { cat clones_check.out; bad "--clones genealogy"; }
FWDX=$(ls "$(dirname "$BIN")"/cellforward-* 2>/dev/null | head -1)
[ -n "$FWDX" ] && { python3 "$ROOT/tests/clones_check.py" "$BIN" 3000 yule "$FWDX" > clones_check2.out 2>&1 && ok "--clones: a pure birth tumor with random sizes agrees with the tumors simulated forward by cellforward (frequencies and times of pairs of cells)" || { cat clones_check2.out; bad "--clones against cellforward"; }; }
rm -rf CL; "$BIN" -n3 -s12 -l50 -u1e-6 --clones n=3,cells=10000,birth=1e-3,time=5000:2000,parent=1:2,fitness=1.5:2 -oCL > cl1.out 2>&1
[ "$(head -1 CL/clones.tsv)" = "$(printf 'replicate\tcell\tclone')" ] && [ "$(tail -n +2 CL/clones.tsv | wc -l | tr -d ' ')" = "36" ] && ok "clones.tsv lists replicate, cell and clone (n=3 replicates x 12 cells)" || bad "clones.tsv"
grep -q "Clone 2: born 5000 generations ago from clone 1" CL/log* && grep -q "Clone 3: .*fitness 2 " CL/log* && ok "the log lists the clones with their size and frequency" || bad "clones in the log"
"$BIN" -n1 -s5 -l50 -u1e-6 --clones n=3,cells=10000,birth=1e-3,time=500 > cl2.out 2>&1
grep -q "needs 2 times" cl2.out && ok "--clones: the lists need n-1 values" || bad "short --clone-time accepted"
"$BIN" -n1 -s5 -l50 -u1e-6 --clones n=2,cells=10000,birth=1e-3,time=500 > cl3.out 2>&1
grep -q "clone-fitness is needed" cl3.out && ok "--clones without fitness is an error" || bad "missing --clone-fitness accepted"
"$BIN" -n1 -s5 -l50 -u1e-6 --clones n=2,cells=10000,birth=1e-3,time=500,fitness=2,founders=3000 > cl4.out 2>&1
grep -q "fewer than its" cl4.out && ok "--clones: a parent with too few cells to found a clone is an error" || bad "too many founders accepted"
"$BIN" -n1 -s5 -l50 -u1e-6 --clones n=2,cells=10000,birth=1e-3,time=500,fitness=2 --deme-n 3:3 --deme-time 1 > cl5.out 2>&1
grep -q "cannot be combined with --demes" cl5.out && ok "--clones with --demes is an error" || bad "--clones with --demes accepted"
"$BIN" -n1 -s5 -l50 -u1e-6 --clones n=2,cells=10000,time=500,fitness=2 > cl6.out 2>&1
grep -q "needs the birth rate" cl6.out && ok "--clones needs the birth rate of the clone 1 (birth=)" || bad "--clones without birth accepted"
"$BIN" -n1 -s5 -l50 -u1e-6 -e10000 --clones n=2,cells=10000,birth=1e-3,time=500,fitness=2 > cl7.out 2>&1
grep -q "ERROR: -e is not used with --clones" cl7.out && ok "--clones with -e is an error (cells= is the tumor size)" || bad "--clones with -e accepted"
"$BIN" -n1 -s5 -l50 -u1e-6 -g1e-3 --clones n=2,cells=10000,birth=1e-3,time=500,fitness=2 > cl8.out 2>&1
grep -q "not as -g, -K or -L" cl8.out && ok "--clones with -g, -K or -L is an error" || bad "--clones with -g accepted"
rm -rf CQ; "$BIN" -n300 -s30 -l10 -u1e-7 --clones n=3,cells=10000,birth=1e-3,death=0,freq=0.2:0.05,fitness=1.5:2,parent=1:2,sizes=fixed -oCQ > cq.out 2>&1
grep -q "Clone 2: born 5067.27 generations ago" CQ/log* && grep -q "Clone 3: .*500 cells.*frequency 0.0500" CQ/log* && ok "--clones freq= with fitness=: the time of birth is derived from the frequency (ln(2000)/1.5e-3 = 5067.27)" || bad "--clone-freq with fitness"
f2=$(awk -F'\t' 'NR>1 && $3==2{n++} END{print n/300/30}' CQ/clones.tsv)
awk -v f="$f2" 'BEGIN{exit !(f>0.18 && f<0.22)}' && ok "--clones freq=: the frequency of the clone in the sample is the one asked ($f2 of 0.2)" || bad "--clone-freq sample frequency ($f2)"
"$BIN" -n1 -s5 -l50 -u1e-6 --clones n=3,cells=10000,birth=1e-3,freq=0.2:0.05,time=5000:2000,parent=1:2 2>&1 | grep -q "Clone 3: born 2000 generations ago from clone 2, 1 founder(s), fitness 3.107" && ok "--clones freq= with time=: the fitness is derived from the frequency" || bad "--clone-freq with time"
"$BIN" -n1 -s5 -l50 -u1e-6 --clones n=2,cells=10000,birth=1e-3,freq=0.2 > cq2.out 2>&1
grep -q "give --clone-time or --clone-fitness" cq2.out && ok "--clones freq= needs time= or fitness=" || bad "--clone-freq alone accepted"
"$BIN" -n1 -s5 -l50 -u1e-6 --clones n=2,cells=10000,birth=1e-3,freq=0.2,time=500,fitness=2 > cq3.out 2>&1
grep -q "give --clone-time or --clone-fitness" cq3.out && ok "--clones freq= with time= and fitness= is an error (one is derived)" || bad "--clone-freq with both accepted"
"$BIN" -n1 -s5 -l50 -u1e-6 --clones n=3,cells=10000,birth=1e-3,freq=0.6:0.5,time=500:300 > cq4.out 2>&1
grep -q "nothing is left for the clone 1" cq4.out && ok "--clones freq=: the frequencies of the clones cannot add up to 1" || bad "--clone-freq sum >= 1 accepted"
rm -rf CQ
"$BIN" -n1 -s5 -l50 -u1e-6 --clones n=2,cells=10000,birth=1e-3,time=500,fitness=2,sizes=big > cl9.out 2>&1
grep -q "sizes= is random" cl9.out && ok "--clones sizes= has to be random or fixed" || bad "--clones bad sizes accepted"
rm -rf CL; "$BIN" -n200 -s30 -l10 -u1e-7 --clones n=2,cells=10000,birth=1e-3,time=5000,fitness=1.5,sizes=random -oCL > cl10.out 2>&1
rm -rf CF; "$BIN" -n200 -s30 -l10 -u1e-7 --clones n=2,cells=10000,birth=1e-3,time=5000,fitness=1.5,sizes=fixed -oCF > cl11.out 2>&1
vr=$(awk 'NR>1 && $3==2{c[$1]++} END{for(r in c){s+=c[r]; q+=c[r]*c[r]} n=200; m=s/n; print (q/n-m*m)}' CL/clones.tsv)
vf=$(awk 'NR>1 && $3==2{c[$1]++} END{for(r in c){s+=c[r]; q+=c[r]*c[r]} n=200; m=s/n; print (q/n-m*m)}' CF/clones.tsv)
awk -v a="$vr" -v b="$vf" 'BEGIN{exit !(a>1.3*b)}' && ok "--clones: with random sizes the number of cells of a clone varies more among replicates than with fixed sizes" || bad "random sizes do not add variance (random $vr, fixed $vf)"
rm -rf CL
cd "$ROOT"

echo "2g68. cellforward: the tumor simulated forward in time, trees for CellCoal"
FWD=$(ls "$(dirname "$BIN")"/cellforward-* 2>/dev/null | head -1)
if [ -z "$FWD" ]; then bad "cellforward is not built (make)"; else
cd "$TMP"
python3 "$ROOT/tests/forward_check.py" "$FWD" "$BIN" 2000 > fwd_check.out 2>&1 && ok "cellforward: neutral birth-death tumors agree with the coalescent of CellCoal (time to the MRCA), and the drivers have the expected number" || { cat fwd_check.out; bad "cellforward check"; }
python3 "$ROOT/tests/forward_theory.py" "$FWD" > fwd_theory.out 2>&1 && ok "cellforward agrees with exact results of the birth-death process (depth of a Yule tree, growth of a clone of fitness F)" || { cat fwd_theory.out; bad "cellforward theory check"; }
mkdir -p lowc && cd lowc && "$BIN" -n1 -s4 -l20000 -e1000 -u1e-7 -C0.5 -#1 -olc > LC.txt 2>&1; rc=$?; "$BIN" -n1 -s4 -l2000 -e1000 -u1e-7 -C0 -#1 -olc0 > LC0.txt 2>&1; rc0=$?; cd ..
[ $rc -eq 0 ] && [ $rc0 -ne 0 ] && ok "-C below 1 (low-pass) is accepted, and -C0 is still an error" || bad "-C below 1 (rc $rc, rc0 $rc0)"
rm -rf FW; "$FWD" --birth 0.5 --death 0.1 --size 5000 --driver-rate 1e-4 --driver-effect 0.1 -s 12 -n 7 --mutation-rate 1e-5 --seed 3 -o FW > fw.out 2>&1
[ "$(wc -l < FW/trees.nwk | tr -d ' ')" = "7" ] && [ "$(head -1 FW/trees.nwk | tr -cd 'c' | wc -c | tr -d ' ')" -ge "12" ] && [ "$(tail -n +2 FW/clones.tsv | wc -l | tr -d ' ')" = "84" ] && ok "cellforward writes one tree per replicate, and the clone of every sampled cell (7 replicates x 12 cells)" || bad "cellforward outputs"
rm -rf FW2; "$FWD" --birth 0.5 --death 0.1 --size 5000 --driver-rate 1e-4 --driver-effect 0.1 -s 12 -n 7 --mutation-rate 1e-5 --seed 3 -o FW2 > /dev/null 2>&1
cmp -s FW/trees.nwk FW2/trees.nwk && ok "cellforward: the same seed gives the same trees" || bad "cellforward is not reproducible"
rm -rf FC; "$BIN" -n7 -s12 -l5000 -u1e-6 -TFW/trees.nwk -1 -#4 -oFC > fc.out 2>&1
[ "$(grep -c '^\[#' FC/snv_gen)" = "7" ] && ok "CellCoal reads the trees of cellforward, one per replicate (-T)" || bad "CellCoal with the trees of cellforward"
rm -rf FT; "$FWD" --birth 0.3 --time 30 -s 10 -n 3 --seed 2 -o FT > ft.out 2>&1
awk 'NR>1 && $3!=30{bad=1} END{exit bad}' FT/summary.tsv && ok "cellforward --time stops at the time asked" || bad "cellforward --time"
"$FWD" --size 100 --time 10 --birth 0.5 > fe1.out 2>&1; grep -q "exactly one of --size" fe1.out && ok "cellforward: --size and --time together are an error" || bad "cellforward --size with --time"
"$FWD" --size 100 > fe2.out 2>&1; grep -q "\-\-birth is needed" fe2.out && ok "cellforward: --birth is needed" || bad "cellforward without --birth"
"$FWD" --birth 0.5 --size 10 -s 20 > fe3.out 2>&1; grep -q "sample (-s) is larger" fe3.out && ok "cellforward: a sample larger than the tumor is an error" || bad "cellforward sample > tumor"
"$FWD" --birth 0.5 --size 100 --driver-rate 0.01 > fe4.out 2>&1; grep -q "driver-effect is needed" fe4.out && ok "cellforward: drivers need an effect" || bad "cellforward drivers without effect"
rm -rf FP; "$FWD" --birth 0.2 --death 0.1 --time 80 --clone-time 30:50 --clone-fitness 2:1.5 --clone-parent 0:1 -s 20 -n 20 --seed 1 -o FP > fp.out 2>&1
awk -F'\t' 'NR>1 && $2==1{n++; if ($5!=30 || $6!=1.5 || $3!=0) bad=1} NR>1 && $2==2{m++; if ($5!=50 || $6!=1.25 || $3!=1) bad=1} END{exit (bad || n==0)}' FP/drivers.tsv && ok "cellforward: clones started by the user (--clone-time, --clone-fitness, --clone-parent) have their time, parent and birth rate death + F (birth - death)" || bad "cellforward clones of the user"
"$FWD" --birth 0.2 --time 80 --clone-time 30:50 --clone-fitness 2 > fe5.out 2>&1; grep -q "same number of values" fe5.out && ok "cellforward: the lists of the clones need the same number of values" || bad "cellforward clone lists"
"$FWD" --birth 0.2 --time 40 --clone-time 50 --clone-fitness 2 > fe6.out 2>&1; grep -q "after the end of the tumor" fe6.out && ok "cellforward: a clone that starts after the end is an error" || bad "cellforward clone after the end"
rm -rf FD; "$FWD" --birth 0.5 --death 0.1 --size 2000 --mutation-rate 0 --division-mutation-rate 1 -s 6 -n 5 --seed 1 -o FD > fd.out 2>&1
tr ':,();' '\n\n\n\n\n' < FD/trees.nwk | grep -E '^[0-9.e+-]+$' | awk '$1!=int($1){bad=1} $1<1{bad=1} END{exit (bad || NR<10)}' && ok "cellforward --division-mutation-rate: the branch lengths are whole numbers of divisions" || bad "cellforward branch lengths in divisions"
for w in growth birth; do rm -rf FG; "$FWD" --birth 0.5 --death 0.1 --size 3000 -s 6 -n 40 --driver-rate 0.02 --driver-dist fixed --driver-effect 0.5 --driver-acts-on $w --seed 2 -o FG > /dev/null 2>&1; eval "r_$w=\$(awk -F'\t' 'NR>1 && \$4==1{print \$6}' FG/drivers.tsv | sort -u | head -1)"; done
[ "$r_growth" = "1.4" ] && [ "$r_birth" = "1.5" ] && ok "cellforward --driver-acts-on: a driver of effect 0.5 gives the relative birth rate 1.4 (net growth x 1.5) or 1.5 (birth x 1.5)" || bad "cellforward driver action (growth $r_growth, birth $r_birth)"
rm -rf FM; "$FWD" --birth 0.5 --death 0.1 --size 3000 -s 6 -n 30 --driver-rate 0.02 --driver-effect 0.3 --max-clones 3 --seed 4 -o FM > /dev/null 2>&1; mx=$(awk -F'\t' 'NR>1{print $5}' FM/summary.tsv | sort -n | tail -1)
[ -n "$mx" ] && [ "$mx" -le 4 ] && ok "cellforward --max-clones 3: at most 3 clones with a driver (and the root) in every tumor" || bad "cellforward --max-clones (max clones $mx)"
for k in 1 1000; do rm -rf FC$k; "$FWD" --birth 1 --death 0 --size 1024 --cycle-shape $k -s 10 -n 60 --seed 5 -o FC$k > /dev/null 2>&1; eval "m_$k=\$(awk -F'\t' 'NR>1{s+=\$3;n++}END{printf \"%.2f\",s/n}' FC$k/summary.tsv)"; done
awk -v a="$m_1" -v b="$m_1000" 'BEGIN{exit !(a>6.5 && a<8.5 && b>9.8 && b<10.7)}' && ok "cellforward --cycle-shape: exponential cells reach 1024 cells in ~7.6 generations, cells that divide almost in synchrony in ~10 ($m_1, $m_1000)" || bad "cellforward --cycle-shape ($m_1, $m_1000)"
rm -rf FD; "$FWD" --birth 0.5 --death 0.3 --size 3000 -s 6 -n 20 --driver-rate 0.01 --driver-dist fixed --driver-effect 0.5 --driver-acts-on death --seed 6 -o FD > /dev/null 2>&1; dd=$(awk -F'\t' 'NR>1 && $4==1{print $6" "$9}' FD/drivers.tsv | sort -u | head -1)
[ "$dd" = "1 0.2" ] && ok "cellforward --driver-acts-on death: a driver of effect 0.5 lowers the death rate from 0.3 to 0.2 (net growth x 1.5), birth unchanged" || bad "cellforward driver on death ($dd)"
rm -rf FJ1 FJ3; "$FWD" --birth 0.5 --death 0.1 --size 2000 -s 8 -n 9 --driver-effect 0.3 --seed 9 -j 1 -o FJ1 > /dev/null 2>&1; "$FWD" --birth 0.5 --death 0.1 --size 2000 -s 8 -n 9 --driver-effect 0.3 --seed 9 -j 3 -o FJ3 > /dev/null 2>&1
same=1; for f in trees.nwk clones.tsv drivers.tsv summary.tsv; do cmp -s FJ1/$f FJ3/$f || same=0; done
[ $same -eq 1 ] && [ -z "$(ls -a FJ3 | grep '^\.part')" ] && ok "cellforward -j: the results do not depend on the number of jobs, and no temporary files are left" || bad "cellforward -j"
python3 "$ROOT/tests/reconstructed_check.py" "$FWD" 300 > fwd_recon.out 2>&1 && ok "cellforward --tree-method reconstructed (tree of the sample drawn directly) agrees with the simulation of all the cells: time and divisions of pairs, 5 tumors (time and size, 2e3-3e4 cells)" || { cat fwd_recon.out; bad "cellforward reconstructed sampler"; }
rm -rf FR; "$FWD" --birth 0.5 --death 0.1 --size 100000000000 -s 10 -n 3 --seed 2 -o FR > FR.log 2>&1; [ "$(wc -l < FR/trees.nwk | tr -d ' ')" = "3" ] && ok "cellforward: a neutral tumor of 1e11 cells is sampled directly (no memory per cell)" || bad "cellforward huge tumor"
"$FWD" --birth 0.5 --death 0.1 --size 100000000000 -s 10 --tree-method forward > fe8.out 2>&1; [ $? -ne 0 ] && ok "cellforward: --tree-method forward refuses a tumor too large to keep in memory" || bad "cellforward forward with huge size"
"$FWD" --birth 0.5 --death 0 --size 100 --driver-rate 0.01 --driver-effect 0.5 --driver-acts-on death -s 5 -o FE > FE.log 2>&1; [ $? -ne 0 ] && grep -q "needs a death rate" FE.log && ok "cellforward --driver-acts-on death without death is refused" || bad "cellforward death driver without death"
rm -rf FA; "$FWD" --birth 0.5 --death 0.1 --size 3000 -s 6 -n 40 --driver-effect 0.3 --seed 4 -o FA > FA.log 2>&1; av=$(awk -F'\t' 'NR>1{s+=$5-1;n++}END{printf "%.1f",s/n}' FA/summary.tsv)
awk -v a="$av" 'BEGIN{exit !(a>=2.5 && a<=5)}' && grep -q "driver-rate not given" FA.log && ok "cellforward: without --driver-rate (with --driver-effect) a tumor has about 3-4 clones with a driver ($av)" || bad "cellforward automatic driver rate ($av)"
for f in 1 10; do rm -rf FM$f; "$FWD" --birth 0.5 --death 0.1 --time 25 --clone-time 12 --clone-fitness 1 --clone-mutation-factor $f -s 30 -n 300 --seed 3 -o FM$f > /dev/null 2>&1; done
len1=$(awk -F'[:,();]' '{t=0; for(i=1;i<=NF;i++) if ($i ~ /^[0-9.]+(e[-+][0-9]+)?$/ && $i ~ /[.e]/) t+=$i; s+=t; n++} END{print s/n}' FM1/trees.nwk)
len10=$(awk -F'[:,();]' '{t=0; for(i=1;i<=NF;i++) if ($i ~ /^[0-9.]+(e[-+][0-9]+)?$/ && $i ~ /[.e]/) t+=$i; s+=t; n++} END{print s/n}' FM10/trees.nwk)
awk -v a="$len1" -v b="$len10" 'BEGIN{exit !(b>1.1*a)}' && awk -F'\t' 'NR>1 && $6!=1 && $6!=10{bad=1} NR>1 && $6==10{h=1} END{exit (bad || !h)}' FM10/clones.tsv && ok "cellforward --clone-mutation-factor: a hypermutator clone lengthens the trees (mean total length $len1 -> $len10) and the factor of every cell is in clones.tsv" || bad "cellforward hypermutator clone ($len1, $len10)"
rm -rf FMR; "$FWD" --birth 0.5 --death 0.1 --size 3000 --mutator-rate 0.005 --mutator-factor 5 -s 30 -n 20 --seed 4 -o FMR > /dev/null 2>&1
awk -F'\t' 'NR>1{ if ($6!=1 && $6!=5 && $6!=25 && $6!=125) bad=1; if ($6>1) h=1} END{exit (bad || !h)}' FMR/clones.tsv && ok "cellforward --mutator-rate: random hypermutators are inherited (the factor of a cell is 5 to a power)" || bad "cellforward random hypermutators"
"$FWD" --birth 0.5 --size 100 --mutator-rate 0.1 > fe7.out 2>&1; grep -q "mutator-factor is needed" fe7.out && ok "cellforward: --mutator-rate needs --mutator-factor" || bad "cellforward mutator without factor"
rm -rf FW FW2 FC FT FP FD FG FM1 FM10 FMR
cd "$ROOT"
fi

echo "2g62. Indels in the reads (--indel-rate): CIGAR with insertions and deletions, truth in indels_truth.vcf"
cd "$TMP"
cat > indelchk.py <<'PYEOF'
import gzip,struct,sys,re,collections
d=sys.argv[1]  # reads dir
ref=''.join(l.strip() for l in open(d+'/ref.fa') if not l.startswith('>'))
truth=collections.defaultdict(dict)  # (cell,hap) -> {(pos,type,seq):1}
cells=None; events=[]
for l in open(d+'/indels_truth.vcf'):
    if l.startswith('##'): continue
    f=l.rstrip('\n').split('\t')
    if l.startswith('#'): cells=f[9:]; continue
    pos=int(f[1]); R,A=f[3],f[4]
    key=(pos,'D',R[1:]) if len(R)>len(A) else (pos,'I',A[1:])
    events.append(key)
    for c,g in zip(cells,f[9:]):
        for h,x in enumerate(g.split('|')):
            if x=='1': truth[(c,h+1)][key]=1
nreads=nmis=nbad=0; seen=collections.defaultdict(set); cover=0; miss=0
cig_re=re.compile(r'')
for c in cells:
    data=gzip.open(f'{d}/{c}.bam','rb').read()
    lt=struct.unpack('<i',data[4:8])[0]; o=8+lt; nref=struct.unpack('<i',data[o:o+4])[0]; o+=4
    for _ in range(nref):
        l2=struct.unpack('<i',data[o:o+4])[0]; o+=4+l2+4
    while o<len(data):
        ln=struct.unpack('<i',data[o:o+4])[0]; o+=4
        r=data[o-4:o+ln]
        rid,pos,lrn,mq,bn,nc,fl,lseq,nrid,npos,tl=struct.unpack('<iiBBHHHiiii',data[o:o+32])
        p=o+32+lrn
        cg=struct.unpack('<%dI'%nc,data[p:p+4*nc]); p+=4*nc
        sq=data[p:p+(lseq+1)//2]; p+=(lseq+1)//2
        seq=''.join('=ACMGRSVTWYHKDBN'[(sq[i//2]>>(4 if i%2==0 else 0))&15] for i in range(lseq))
        p+=lseq
        aux=data[p:o+ln]
        m=re.search(rb'HPC(.)',aux); hap=m.group(1)[0]
        o+=ln
        nreads+=1
        rp=pos; qp=0; shown=[]
        for x in cg:
            L,op=x>>4,x&15
            if op==0:
                for i in range(L):
                    if seq[qp+i]!=ref[rp+i]: nmis+=1
                rp+=L; qp+=L
            elif op==1:
                shown.append((rp,'I',seq[qp:qp+L])); qp+=L   # rp = pos(0-based) of next ref base ; anchor 1-based = rp
            elif op==2:
                shown.append((rp,'D',ref[rp:rp+L])); rp+=L
            elif op==4: qp+=L
        for (a,t,s) in shown:
            if (a,t,s) not in truth[(c,hap)]: nbad+=1; print('unexpected',c,hap,a,t,s)
            seen[(c,hap)].add((a,t,s))
        # coverage check: events (1-based anchor pos a) fully spanned by read -> must be shown
        end=rp
        for (a,t,s) in truth[(c,hap)]:
            if t=='D' and pos+1<=a-0 and a+len(s)+1<=end: # anchor base idx a-1 >= pos, and base after deletion within read
                if (a,t,s) not in shown:
                    miss+=1
            if t=='I' and pos+1<=a and a+1<=end:
                if (a,t,s) not in shown:
                    miss+=1
print('reads',nreads,'M mismatches',nmis,'unexpected indels',nbad,'spanned-but-missing',miss,'events',len(events))

assert nmis < 200 and nbad == 0 and miss == 0 and len(events) > 20, (nmis, nbad, miss, len(events))
PYEOF
rm -rf ID1 ID2; "$BIN" -n1 -s8 -l30000 -u1e-9 -e1000 -b1 -C20 -E0 -D0 --bam --fastq --indel-rate 300 -#3 -oID1 >/dev/null 2>&1
python3 indelchk.py ID1/reads_1 >/dev/null 2>&1 && ok "every indel of a read is in the truth, every spanned indel is in the reads, and the aligned bases match the reference" || bad "indels in reads"
grep -q "DEL" ID1/reads_1/indels_truth.vcf && grep -q "INS" ID1/reads_1/indels_truth.vcf && ok "indels_truth.vcf has insertions and deletions" || bad "indels_truth.vcf"
CELLCOAL_READ_TILE=1024 "$BIN" -n1 -s8 -l30000 -u1e-9 -e1000 -b1 -C20 -E0 -D0 --bam --fastq --indel-rate 300 -#3 -oID2 >/dev/null 2>&1
[ "$(zcat ID1/reads_1/cell0003.bam | md5sum)" = "$(zcat ID2/reads_1/cell0003.bam | md5sum)" ] && ok "the BAM with indels does not depend on the size of the tile" || bad "indel BAM depends on tile"
"$BIN" -n1 -s5 -l500 -u1e-6 -C10 -E0.01 --bam --indel-rate 0.1 --indel-max 140 > id1.out 2>&1
grep -q "bad --indel-rate" id1.out && ok "a bad --indel-max is an error" || bad "bad indel option accepted"
"$BIN" -n1 -s5 -l500 -u1e-6 -C10 --indel-rate 0.1 > id2.out 2>&1
grep -q "needs --bam" id2.out && ok "--indel-rate without reads mode is an error" || bad "--indel-rate accepted in the genotypes mode"
"$BIN" -n1 -s6 -l30000 -u1e-8 -e1000 -b1 -C20 -E0 -D0 --bam --indel-rate 30 --seg-gain-rate 20 --seg-gain-len 3000 -#2 -oID3 >/dev/null 2>&1
python3 indelchk.py ID3/reads_1 >/dev/null 2>&1 && ok "indels with segment gains (copies of the haplotype) are in the reads as in the truth" || bad "indels with segment gains"
cat > homochk.py <<'PYEOF'
import gzip,struct,sys,re,glob
d=sys.argv[1]
ref=''.join(l.strip() for l in open(d+'/ref.fa') if not l.startswith('>'))
def run_len(p):  # homopolymer run containing 0-based p in the reference
    a=p
    while a>0 and ref[a-1]==ref[p]: a-=1
    b=p
    while b+1<len(ref) and ref[b+1]==ref[p]: b+=1
    return a,b
nr=nI=nD=nmis=bad=0
for fn in sorted(glob.glob(d+'/cell*.bam')):
    data=gzip.open(fn,'rb').read()
    lt=struct.unpack('<i',data[4:8])[0]; o=8+lt; nref=struct.unpack('<i',data[o:o+4])[0]; o+=4
    for _ in range(nref):
        l2=struct.unpack('<i',data[o:o+4])[0]; o+=4+l2+4
    while o<len(data):
        ln=struct.unpack('<i',data[o:o+4])[0]; o+=4
        rid,pos,lrn,mq,bn,nc,fl,lseq,nrid,npos,tl=struct.unpack('<iiBBHHHiiii',data[o:o+32])
        p=o+32+lrn
        cg=struct.unpack('<%dI'%nc,data[p:p+4*nc]); p+=4*nc
        sq=data[p:p+(lseq+1)//2]
        seq=''.join('=ACMGRSVTWYHKDBN'[(sq[i//2]>>(4 if i%2==0 else 0))&15] for i in range(lseq))
        o+=ln; nr+=1
        rp=pos; qp=0
        for x in cg:
            L,op=x>>4,x&15
            if op==0:
                for i in range(L):
                    if seq[qp+i]!=ref[rp+i]: nmis+=1
                rp+=L; qp+=L
            elif op==1:
                nI+=1; a,b=run_len(rp-1)
                if b-a+1<3 or seq[qp:qp+L]!=ref[rp-1]*L: bad+=1
                qp+=L
            elif op==2:
                nD+=1; a,b=run_len(rp)
                if b-a+1<3 or rp!=a+1 and rp!=a: pass
                if b-a+1<3: bad+=1
                rp+=L
            elif op==4: qp+=L
assert nI > 20 and nD > 20 and bad == 0 and nmis == 0, (nI, nD, bad, nmis)

PYEOF
"$BIN" -n1 -s6 -l30000 -u1e-9 -e1000 -b1 -C20 -E0 -D0 --bam --fastq --homopolymer-error 0.01 -#2 -oID4 >/dev/null 2>&1
python3 homochk.py ID4/reads_1 >/dev/null 2>&1 && ok "--homopolymer-error: insertions and deletions of one base only in homopolymers of 3 or more bases, the rest of the read matches the reference" || bad "homopolymer indel errors"
CELLCOAL_READ_TILE=1024 "$BIN" -n1 -s6 -l30000 -u1e-9 -e1000 -b1 -C20 -E0 -D0 --bam --fastq --homopolymer-error 0.01 -#2 -oID5 >/dev/null 2>&1
[ "$(zcat ID4/reads_1/cell0002.bam | md5sum)" = "$(zcat ID5/reads_1/cell0002.bam | md5sum)" ] && ok "the BAM with homopolymer errors does not depend on the size of the tile" || bad "homopolymer BAM depends on tile"
printf "chr1 1 15000 0.05\nchr1 15001 30000 4\n" > idrr.bed
"$BIN" -n1 -s6 -u1e-9 -e1000 -b1 -C10 -E0 -D0 -l15000 --blocks 2 --contig 1 --bam --indel-rate 300 --region-rates idrr.bed -#2 -oID6 >/dev/null 2>&1
awk '!/^#/{if($2<=15000)a++;else b++}END{exit !(b > 10 * (a + 1))}' ID6/reads_1/indels_truth.vcf && ok "--region-rates: the indels follow the regional multipliers" || bad "indels do not follow --region-rates"
rm -rf ID1 ID2 ID3 ID4 ID5 ID6
cd "$ROOT"

echo "2g64. Indels: long tail (--indel-tail), repeats (--indel-repeat) and germline (--germline-indel-rate)"
cd "$TMP"
cat > replen.py <<'PY'
import sys
d = sys.argv[1]; mode = sys.argv[2]
ref = ''.join(l.strip() for l in open(d + '/ref.fa') if not l.startswith('>'))
def inrep(p):
    for k in range(1, 5):
        m = lambda j: 0 <= j - k and j < len(ref) and ref[j] == ref[j - k]
        j0 = p + 1 if m(p + 1) else (p + 1 + k if m(p + 1 + k) else None)
        if j0 is None: continue
        jl = j0
        while m(jl - 1): jl -= 1
        jr = j0
        while m(jr + 1): jr += 1
        if (jr - jl + 1) + k >= (4 if k == 1 else 3 * k): return True
    return False
n = r = longest = germ = allc = 0
for l in open(d + '/indels_truth.vcf'):
    if l[0] == '#': continue
    f = l.split('\t'); n += 1; r += inrep(int(f[1]) - 1)
    info = dict(kv.split('=') if '=' in kv else (kv, 1) for kv in f[7].split(';'))
    longest = max(longest, int(info['LEN'])); germ += 'GERMLINE' in info; allc += int(info['NC']) == len(f) - 9
if mode == 'tail': sys.exit(0 if n > 30 and 20 <= longest <= 80 else 1)
if mode == 'repeat': sys.exit(0 if n > 60 and r / n > 0.35 else 1)
if mode == 'base': sys.exit(0 if n > 60 and r / n < 0.2 else 1)
if mode == 'germline': sys.exit(0 if n > 10 and germ == n and allc == n else 1)
PY
if command -v python3 >/dev/null; then
  "$BIN" -n1 -s6 -#4 -l40000 -u1e-8 -e1000 -C15 -E0 --bam --indel-rate 20 --indel-max 80 --indel-tail 0.5 -oIT >/dev/null 2>&1
  python3 replen.py IT/reads_1 tail && python3 indelchk.py IT/reads_1 >/dev/null 2>&1 && ok "--indel-tail: long indels (up to --indel-max 80) in the truth and the reads" || bad "indel tail"
  "$BIN" -n1 -s6 -#4 -l60000 -u1e-8 -e1000 -C8 -E0 --bam --indel-rate 12 --indel-repeat 25 -oIR >/dev/null 2>&1
  "$BIN" -n1 -s6 -#4 -l60000 -u1e-8 -e1000 -C8 -E0 --bam --indel-rate 12 -oIB >/dev/null 2>&1
  python3 replen.py IR/reads_1 repeat && python3 replen.py IB/reads_1 base && python3 indelchk.py IR/reads_1 >/dev/null 2>&1 && ok "--indel-repeat 25: most indels are slippages in repeats (and a few without it), reads agree with the truth" || bad "indel repeat"
  "$BIN" -n1 -s6 -#4 -l40000 -u1e-8 -e1000 -C15 -E0 --bulk-depth 15 --bam --germline-indel-rate 1e-3 -oIG >/dev/null 2>&1
  python3 replen.py IG/reads_1 germline && python3 indelchk.py IG/reads_1 >/dev/null 2>&1 && ok "--germline-indel-rate: germline indels (flag GERMLINE) in all the cells, reads agree with the truth" || bad "germline indels"
else echo "  skip (python3 not installed)"; fi
"$BIN" -n1 -s3 -l1000 -u1e-7 -e1000 -C5 --bam --indel-rate 1 --indel-max 101 -oIE >/dev/null 2>&1 && bad "--indel-max 101 accepted" || ok "--indel-max above 100 is rejected"
"$BIN" -n1 -s3 -l1000 -u1e-7 -e1000 -C5 --bam --indel-rate 1 --indel-tail 2 -oIE >/dev/null 2>&1 && bad "--indel-tail 2 accepted" || ok "--indel-tail outside 0 to 1 is rejected"
rm -rf IT IR IB IG IE
cd "$ROOT"


echo "2g54. Documentation: the commands of the manual and README work, every long option is in the manual; a mistyped option gets a suggestion"
cd "$ROOT"
if command -v python3 >/dev/null; then
  python3 tests/manual_check.py "$BIN" > "$TMP/manual_check.txt" 2>&1 && ok "manual and README commands run; every long option of --help is in the manual" || { bad "manual_check.py"; cat "$TMP/manual_check.txt"; }
else echo "  skip (python3 not installed)"; fi
cd "$TMP"
"$BIN" -n1 -s4 -l100 --seg-lenn 5 2>&1 | grep -q "did you mean --seg-len?" && ok "unknown long option: suggests the closest one" || bad "no suggestion for a mistyped option"
"$BIN" --check -n3 -s10 -l500 -u1e-6 --sweep "C=2,5" >Chk.txt 2>&1 && grep -q "not run" Chk.txt && ok "--check with --sweep checks the options and does not run the grid" || bad "--check --sweep"
cd "$TMP" && mkdir -p dd && cd dd
"$BIN" -n1 -s5 -l100 -u1e-7 -e1000 -d1e-7 --blocks 2 >D1.txt 2>&1; r1=$?
"$BIN" -n1 -s5 -l200 -u1e-7 -e1000 -d1e-7 --bam >D2.txt 2>&1; r2=$?
"$BIN" -n1 -s5 -l100 -u1e-7 -e1000 -d1e-7 >D3.txt 2>&1; r3=$?
[ $r1 -ne 0 ] && [ $r2 -ne 0 ] && [ $r3 -eq 0 ] && grep -q "indel-rate" D1.txt && ok "-d is refused with blocks and with reads (it points to --indel-rate), and works for isolated sites" || bad "-d with blocks or reads"
cd "$ROOT"

echo "2g69. --wga mda|pta: hidden shortcut with theoretical starting values; explicit options are kept"
cd "$TMP" && mkdir -p wga && cd wga
"$BIN" -n1 -s4 -l20000 -e1000 -u1e-7 -C10 -#1 --wga mda -ow1 >W1.txt 2>&1 && grep -q "gain-sd 0.8" W1.txt && grep -q -- "-D0.2" W1.txt && ok "--wga mda sets the starting options" || bad "--wga mda"
"$BIN" -n1 -s4 -l20000 -e1000 -u1e-7 -C10 -#1 --wga pta -D0.4 --gain-sd 0.1 -ow2 >W2.txt 2>&1 && ! grep -q -- "-D0.05" W2.txt && ! grep -q "gain-sd 0.2" W2.txt && grep -q "gain-len 20000" W2.txt && ok "options given by the user are not replaced" || bad "--wga with explicit options"
printf 'replicates: 1\ncells: 4\nsites: 20000\npopulation_size: 1000\nmutation_rate: 1e-7\ncoverage: 10\nwga: mda\nseed: 1\n' > w.yaml
"$BIN" --config w.yaml -ow3 >W3.txt 2>&1 && grep -q "^--wga mda" W3.txt && ok "wga: key in a configuration file" || bad "wga key"
"$BIN" -n1 -s4 -l1000 -e1000 -u1e-7 --wga foo >W4.txt 2>&1; [ $? -ne 0 ] && grep -q "mda, pta or none" W4.txt && ok "--wga with an unknown method is refused" || bad "--wga foo"
"$BIN" --help 2>&1 | grep -q -- "--wga" && bad "--wga must stay hidden in --help" || ok "--wga is not shown in --help"
cd "$ROOT"

[ $fail -eq 0 ] && echo "All tests passed" || echo "Some tests FAILED"
exit $fail
