#!/usr/bin/env python3
"""Speed regression check: runs fixed workloads with two binaries (the reference, e.g. the last release, and the new one) on the
same machine, one thread, and fails if the new one is slower than the reference by more than the allowed ratio.

Usage: python3 tests/bench.py REFERENCE_BIN NEW_BIN [--max-ratio 1.25] [--repeats 3]

Comparing two builds on the same machine at the same time avoids the noise of absolute timings in shared CI runners.
The best (minimum) of the repeats is used for each workload.
"""
import argparse, os, subprocess, sys, tempfile, time

WORK = {
    "tree+mutations (no reads)": "-n20 -s100 -l50000 -u1e-6 -e10000 -b1 -#1 -y0",
    "reads, VCF (dense)":        "-n3 -s60 -l20000 -u1e-6 -e10000 -b1 -C10 -V5 -E0.001 -P0.2 0.02 -#1 -y0",
    "reads + doublets + -J":     "-n3 -s60 -l20000 -u1e-6 -e10000 -b1 -C10 -V5 -E0.001 -J4 -B0.1 0.01 -#1 -y0",
}

def run(binary, args):
    with tempfile.TemporaryDirectory() as d:
        env = dict(os.environ, OMP_NUM_THREADS="1")
        t = time.perf_counter()
        r = subprocess.run([binary] + args.split(), cwd=d, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        dt = time.perf_counter() - t
    if r.returncode != 0:
        sys.exit("FAIL: %s exited with %d for: %s" % (binary, r.returncode, args))
    return dt

def best(binary, args, repeats):
    return min(run(binary, args) for _ in range(repeats))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("reference"); ap.add_argument("new")
    ap.add_argument("--max-ratio", type=float, default=1.25)
    ap.add_argument("--repeats", type=int, default=3)
    a = ap.parse_args()
    ref, new = os.path.abspath(a.reference), os.path.abspath(a.new)
    bad = 0
    print("%-28s %10s %10s %8s" % ("workload", "reference", "new", "ratio"))
    for name, args in WORK.items():
        # interleave to share any drift of the machine
        tr = best(ref, args, a.repeats); tn = best(new, args, a.repeats)
        ratio = tn / tr
        flag = "" if ratio <= a.max_ratio else "  <-- SLOWER"
        bad += ratio > a.max_ratio
        print("%-28s %9.2fs %9.2fs %8.2f%s" % (name, tr, tn, ratio, flag))
    if bad:
        sys.exit("FAIL: %d workload(s) slower than %.2fx the reference" % (bad, a.max_ratio))
    print("Speed OK (no workload more than %.2fx slower than the reference)" % a.max_ratio)

main()
