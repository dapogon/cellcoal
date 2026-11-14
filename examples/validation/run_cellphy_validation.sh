#!/bin/sh
# Simulate single-cell data with CellCoal, infer the tree with CellPhy and compare it with the true tree.
#
# Usage:   ./run_cellphy_validation.sh [clean|noisy] [number of replicates]
# Needs:   a compiled CellCoal (make) and a CellPhy folder (https://github.com/amkozlov/cellphy).
#   CELLCOAL  path to the CellCoal binary   (default: the one in ../../bin)
#   CELLPHY   path to the CellPhy folder    (default: ~/cellphy)
#   THREADS   threads for CellPhy           (default: 2)
# Output: table with one line per replicate (cells, RF distance, normalized RF, correct splits).
set -e

SCENARIO=${1:-clean}
NREPS=${2:-5}
HERE=$(cd "$(dirname "$0")" && pwd)
CELLCOAL=${CELLCOAL:-$(ls "$HERE"/../../bin/cellcoal-* 2>/dev/null | sort | tail -1)}
CELLPHY=${CELLPHY:-$HOME/cellphy}
THREADS=${THREADS:-2}
WORK=${WORKDIR:-validation_$SCENARIO}

case "$(uname -s)" in
  Darwin) RAXML="$CELLPHY/bin/raxml-ng-cellphy-osx" ;;
  *)      RAXML="$CELLPHY/bin/raxml-ng-cellphy-linux" ;;
esac
[ -x "$CELLCOAL" ] || { echo "CellCoal binary not found; set CELLCOAL=/path/to/cellcoal (or run make)" >&2; exit 1; }
[ -x "$RAXML" ]    || { echo "CellPhy not found at $CELLPHY; set CELLPHY=/path/to/cellphy" >&2; exit 1; }

case "$SCENARIO" in
  clean) OPTS="-C30 -V5 -E0.001" ;;                                  # high coverage, almost no errors
  noisy) OPTS="-C15 -V5 -E0.01 -D0.1 -A0.05 0.01 0" ;;               # low coverage, ADO and amplification errors
  *) echo "scenario must be clean or noisy" >&2; exit 1 ;;
esac

rm -rf "$WORK"; mkdir -p "$WORK/sim" "$WORK/cellphy"
echo "Simulating $NREPS replicates ($SCENARIO): 20 cells, 3000 sites"
(cd "$WORK/sim" && "$CELLCOAL" -n"$NREPS" -s20 -l3000 -u1e-6 -e10000 -g1e-4 -b1 -m2 $OPTS -6 -v -#11 > /dev/null 2>&1)

printf 'rep\tcells\tRF\tnormRF\tcorrect_splits\n'
for vcf in "$WORK"/sim/Results/vcf_dir/vcf.*; do
  id=${vcf##*.}
  "$RAXML" --search --msa "$vcf" --model GT10+FO --msa-format VCF --prefix "$WORK/cellphy/r$id" \
           --seed 1 --threads "$THREADS" --tree pars{5},rand{5} --redo > "$WORK/cellphy/log$id.txt" 2>&1
  printf '%s\t%s\n' "$id" "$(python3 "$HERE/rf_distance.py" "$WORK/sim/Results/trees_dir/trees.$id" "$WORK/cellphy/r$id.raxml.bestTree")"
done | tee "$WORK/results.tsv"
echo
awk -F'\t' 'NR>0 && $1!="rep" {s+=$4; n++} END {if (n) printf "mean normalized RF over %d replicates: %.3f\n", n, s/n}' "$WORK/results.tsv"
