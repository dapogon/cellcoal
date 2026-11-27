#!/usr/bin/env python3
"""cellforward: the tree of the sample drawn directly (--tree-method reconstructed, a coalescent point process) against the simulation of all
   the cells (--tree-method forward), for neutral tumors. Compared: the mean time since the coalescence of a pair of sampled cells, the mean
   number of divisions of a pair (with --division-mutation-rate 1) and the depth of the tree.
   usage: reconstructed_check.py CELLFORWARD [replicates=400]"""
import math, os, subprocess, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from newick_tools import newick_pair_times

CF = sys.argv[1]
REPS = int(sys.argv[2]) if len(sys.argv) > 2 else 400
S = 20
tmp = tempfile.mkdtemp()

def run(method, stop, b, d, mode, seed):
    out = os.path.join(tmp, "%s%s%s" % (method, mode, seed))
    mut = ["--mutation-rate", "1", "--division-mutation-rate", "0"] if mode == "time" else ["--mutation-rate", "0", "--division-mutation-rate", "1"]
    subprocess.run([CF, "--birth", str(b), "--death", str(d)] + stop + ["-s", str(S), "-n", str(REPS), "--seed", str(seed), "--tree-method", method, "-o", out] + mut, capture_output=True, check=True)
    pair, root = [], []
    for nwk in open(os.path.join(out, "trees.nwk")).read().split("\n"):
        if nwk.strip():
            v = list(newick_pair_times(nwk.strip(), S, mode == "div").values())
            pair.append(sum(v) / len(v))
    return pair

def stat(x):
    m = sum(x) / len(x)
    return m, math.sqrt(sum((v - m) ** 2 for v in x) / (len(x) - 1) / len(x)), math.sqrt(sum((v - m) ** 2 for v in x) / (len(x) - 1))

ok = True
print("%d replicates, %d sampled cells; mean over pairs, forward / reconstructed (z = difference in standard errors; sd ratio)" % (REPS, S))
print("%-34s %-30s %-30s" % ("tumor", "time since the coalescence", "divisions of the pair"))
for b, d, stop in ((0.5, 0.1, ["--size", "2000"]), (0.5, 0.4, ["--size", "2000"]), (0.5, 0.1, ["--size", "30000"]), (0.5, 0.1, ["--time", "15"]), (0.5, 0.4, ["--time", "40"])):
    cols = []
    for mode in ("time", "div"):
        f, r = stat(run("forward", stop, b, d, mode, 5)), stat(run("reconstructed", stop, b, d, mode, 6))
        z = (f[0] - r[0]) / math.sqrt(f[1] ** 2 + r[1] ** 2)
        sdr = r[2] / f[2]
        good = abs(z) < 3.5 and abs(f[0] / r[0] - 1) < 0.05 and 0.8 < sdr < 1.2
        ok &= good
        cols.append("%.2f / %.2f z=%+.1f sd %.2f%s" % (f[0], r[0], z, r[2] / f[2], "" if good else " FAIL"))
    print("%-34s %-30s %-30s" % ("b=%g d=%g %s" % (b, d, " ".join(stop)), cols[0], cols[1]))
sys.exit(0 if ok else 1)
