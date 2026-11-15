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
close() { python3 -c "import sys; sys.exit(0 if abs($1-$2)<=$3 else 1)"; }
M="-n1 -s30 -l3000 -u3e-6 -e10000 -g0 -b1 -1 -Y -#5"
mkdir -p "$TMP/ado" "$TMP/gen" "$TMP/fj"
(cd "$TMP/ado" && "$BIN" $M -D0.3 >/dev/null 2>&1); r=($(python3 "$TMP/rates.py" "$TMP/ado/Results"))
close ${r[0]} 0.30 0.03 && ok "fixed ADO (-D0.3): ${r[0]} changed het genotypes (expected 0.30)" || bad "fixed ADO: ${r[0]} (expected 0.30)"
close ${r[2]} 0.027 0.012 && ok "ADO both alleles: ${r[2]} (expected 0.027)" || bad "ADO both alleles: ${r[2]} (expected 0.027)"
(cd "$TMP/gen" && "$BIN" $M -G0.1 0.001 >/dev/null 2>&1); r=($(python3 "$TMP/rates.py" "$TMP/gen/Results"))
close ${r[1]} 0.10 0.01 && ok "genotyping error (-G0.1): ${r[1]} (expected 0.10)" || bad "genotyping error: ${r[1]} (expected 0.10)"
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
(cd "$TMP/dbl" && "$BIN" -n1 -s100 -l3000 -u1e-6 -e10000 -g0 -b1 -C30 -V1000 -E0 -B0.2 0.001 -v -#5 >/dev/null 2>&1)
r=$(python3 "$TMP/reads.py" dbl "$TMP/dbl/Results/vcf_dir/vcf.0001")
close $r 0.20 0.08 && ok "doublets (-B0.2): fraction of doublet cells $r" || bad "doublets: $r (expected 0.20)"
(cd "$TMP/sig" && "$BIN" -n1 -s20 -l500000 -u1e-8 -e10000 -g0 -b1 -S1 1 1.0 -C20 -V1000 -E0 -v -#5 >/dev/null 2>&1)
r=($(python3 "$TMP/reads.py" sig "$TMP/sig/Results/vcf_dir/vcf.0001"))
close ${r[0]} 0.675 0.05 && ok "signature 1: C>T fraction ${r[0]} (expected 0.675)" || bad "signature 1: C>T ${r[0]} (expected 0.675)"
(cd "$TMP/pq" && "$BIN" $M -P0.2 0.01 >/dev/null 2>&1); r=($(python3 "$TMP/rates.py" "$TMP/pq/Results"))
close ${r[0]} 0.20 0.06 && ok "ADO per cell (-P0.2 0.01): mean ${r[0]} (expected 0.20)" || bad "ADO per cell: ${r[0]} (expected 0.20)"

[ $fail -eq 0 ] && echo "All tests passed" || echo "Some tests FAILED"
exit $fail
