#!/usr/bin/env python3
"""Check of cellforward (the tumor simulated forward in time) against CellCoal (the coalescent) and against what it has to do by construction.

1. Neutral tumors, birth-death, the same size: the mean time to the MRCA of a sample of cells in cellforward and in CellCoal (-K birth -L death,
   the model of Ohtsuki and Innan, where -e is the number of cells today) has to agree (the coalescent leaves out the first generations of the tumor,
   which are stochastic, so a few percent are expected).
2. Drivers without effect (--driver-dist fixed --driver-effect 0) in a pure birth tumor: every division has two daughters, and every daughter gets a
   driver with the probability P, so the number of drivers is 2 P (N - 1).

usage: forward_check.py CELLFORWARD CELLCOAL [replicates]
"""
import os
import subprocess
import sys
import tempfile

FWD, CC = sys.argv[1], sys.argv[2]
REPS = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
N, S = 1000, 10
bad = 0


def column(path, name):
    rows = [l.rstrip("\n").split("\t") for l in open(path)]
    k = rows[0].index(name)
    return [float(r[k]) for r in rows[1:]]


with tempfile.TemporaryDirectory() as tmp:
    for birth, death in ((0.1, 0.0), (0.5, 0.1)):
        out1, out2 = os.path.join(tmp, "f"), os.path.join(tmp, "c")
        subprocess.run([FWD, "--birth", str(birth), "--death", str(death), "--size", str(N), "-s", str(S), "-n", str(REPS),
                        "--seed", "5", "-o", out1], check=True, stderr=subprocess.DEVNULL)
        subprocess.run([CC, "-n%d" % REPS, "-s%d" % S, "-l10", "-u1e-9", "-e%d" % N, "-K%g" % birth, "-L%g" % death, "-o" + out2],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        f = column(os.path.join(out1, "summary.tsv"), "tmrca")
        mf = sum(f) / len(f)
        mc = 0.0
        for line in open(os.path.join(out2, "log")):
            if "Mean time to the MRCA" in line:
                mc = float(line.split()[-1]) * N
        ratio = mf / mc
        ok = abs(ratio - 1.0) < 0.08
        bad += not ok
        print("birth %g death %g: mean time to the MRCA of %d cells (generations): cellforward %.1f, CellCoal %.1f, ratio %.3f %s"
              % (birth, death, S, mf, mc, ratio, "ok" if ok else "FAIL"))
    p = 0.01
    out = os.path.join(tmp, "d")
    subprocess.run([FWD, "--birth", "0.1", "--size", str(N), "-s", str(S), "-n", "300", "--seed", "7", "--driver-rate", str(p),
                    "--driver-dist", "fixed", "--driver-effect", "0", "-o", out], check=True, stderr=subprocess.DEVNULL)
    c = column(os.path.join(out, "summary.tsv"), "clones")
    mean = sum(c) / len(c) - 1.0            # the first clone has no driver
    exp = 2 * p * (N - 1)
    ok = abs(mean / exp - 1.0) < 0.05
    bad += not ok
    print("drivers in a pure birth tumor of %d cells: expected 2 P (N - 1) = %.1f, cellforward %.1f %s" % (N, exp, mean, "ok" if ok else "FAIL"))
print("OK" if not bad else "FAILED")
sys.exit(1 if bad else 0)
