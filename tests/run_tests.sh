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
import re
t=open('$TMP/ad/o.txt',errors='ignore').read()
v=[int(x) for x in re.findall(r'countDropped = (\d+)',t)]
print(round(sum(v)/(len(v)*42.0),4))")
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
for l in open('$TMP/j0/Results/vcf',errors='ignore'):
    if l.startswith('1\t'):
        k+=1
        f=l.rstrip('\n').split('\t')
        for c,x in enumerate(f[9:-1]): n[c]+=int(x.split(':')[1])
m=[x/k for x in n[:100]]
print(round(st.pvariance(m)/st.mean(m)**2,4))")
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
if "$BIN" --cov-ado-corr 0.5 -s5 2>&1 | grep -q 'needs coverage variation among cells (-J) and ADO variation'; then ok "--cov-ado-corr needs -J and -P"; else bad "cov-ado validation"; fi

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
if "$BIN" --bcf -z -s5 2>&1 | grep -q 'cannot be combined' && "$BIN" --bcf -n3 -s5 2>&1 | grep -q 'needs -v'; then ok "--bcf option checks"; else bad "--bcf option checks"; fi

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

[ $fail -eq 0 ] && echo "All tests passed" || echo "Some tests FAILED"
exit $fail
