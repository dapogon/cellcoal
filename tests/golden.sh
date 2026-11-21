#!/bin/bash
# golden.sh BIN OUTFILE : runs a battery of configurations and writes md5 of every output file
BIN=$(readlink -f "$1"); OUT=$2; W=$(mktemp -d /home/claude/work/g.XXXXXX); cd $W
norm() { f=$1; case $f in *.bcf) h=$(gzip -dc < "$f" 2>/dev/null | sed -E 's/filedate=[A-Za-z]{3} [A-Za-z]{3} [ 0-9]{2} [0-9:]{8} [0-9]{4}/filedate=X/' | md5sum);; *.bam|*.gz|*.bai|*.csi) h=$(md5sum < "$f");; *) h=$(grep -av '^##filedate' "$f" | md5sum);; esac; echo "${h%% *}  $(basename $f)"; }
run() { name=$1; shift; mkdir $name; (cd $name; "$BIN" "$@" -o$name -#11 >/dev/null 2>&1; echo "$name exit=$?"; find $name -type f ! -name '*.log' ! -name 'settings*' ! -name log ! -name run.json | sort | while read f; do norm "$f"; done ); }
{
run g1  -n3 -s6 -l400 -u1e-6 -e1000 -b1 -1 -2 -3 -4 -5 -6 -7
run g2  -n2 -s8 -l600 -u1e-6 -e1000 -b1 -D0.1 -E0.001 -G0.02 0.0001 -A0.01 0.0002 0 -1 -3
run g3  -n3 -s8 -l500 -u1e-6 -e1000 -b1 -C10 -E0.001 -A0.01 0.0002 0 -D0.1 -1 -3 -4 -5 -8 -9
run g4  -n2 -s10 -l800 -u1e-7 -e1000 -b1 -C8 -E0.002 -A0.01 0.0002 0 -D0.2 --ado-block 200 -g1 -c0.001 -4 -5 -6
run g5  -n1 -s8 -l500 -u1e-6 -e1000 -b1 -C10 -E0.001 -J0.3 -P0.1 0.001 --min-alt-reads 2 -8 -9 -0 --bcf
run g6  -n2 -s8 -l400 -u1e-6 -e1000 -b1 -C10 -E0.001 -A0.01 0.0002 0 --amp-cell-var 0.00005 --gain-sd 0.4 --gain-len 100 -8
run g7  -n1 -s6 -l30000 -u2e-8 -e1000 -b1 -C10 -E0.001 -A0.01 0.0002 0 -D0.1 --ado-block 3000 --gain-sd 0.4 --gain-len 5000 --bam --fastq --bulk-depth 10
run g8  -n1 -s6 -l20000 -u2e-8 -e1000 -b1 -C10 -E0.002 -A0.01 0.0002 0 --qual-sd 0.7 --qual-slope 2 --strand-bias 0.6 --bam --fastq
run g9  -n1 -s6 -u1e-8 -e1000 -b1 -C10 -E0.001 -l20000 --blocks 3 --contig 1 -8
run g10 -n1 -s6 -u1e-8 -e1000 -b1 -C10 -E0.001 -l10000 --blocks 3 --contig 1 --bam
run g11 -n2 -s6 -l300 -u1e-6 -e500 -g1 -H -2 -3 -k2
run g12 -n2 -s6 -l300 -u1e-6 -e500 -b1 -3 -4 -k1
run g13 -n2 -s6 -l300 -u1e-6 -e500 -b1 -m0 -3 -1
run g14 -n2 -s6 -l300 -u1e-6 -e500 -b1 -C5 -E0.001 --fast -8 -9
} > $OUT 2>&1
rm -rf $W
