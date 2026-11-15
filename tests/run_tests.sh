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

[ $fail -eq 0 ] && echo "All tests passed" || echo "Some tests FAILED"
exit $fail
