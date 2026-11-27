#!/usr/bin/env python3
"""cellforward against exact results of the birth-death process (no other program needed).
   usage: forward_theory.py CELLFORWARD
   1. Yule tree (birth only, all cells sampled): the mean number of divisions from the root to a cell is 2 (H_N - 1).
   2. A clone of fitness F that starts at the time t0 grows, on average over all the replicates (also those where it dies out),
      as exp(F g (T - t0)), with g the net growth rate of the first clone."""
import math, os, re, subprocess, sys, tempfile

F = sys.argv[1]
tmp = tempfile.mkdtemp()

def run(args, out):
    r = subprocess.run([F] + args + ["-o", out], capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stderr); sys.exit(1)

def depths(s):
    sys.setrecursionlimit(100000)
    def parse(i):
        if s[i] == '(':
            res = []
            i += 1
            while True:
                r, i = parse(i)
                res += r
                if s[i] == ',':
                    i += 1
                    continue
                break
            i += 1
            m = re.match(r':([0-9.eE+-]+)', s[i:])
            bl = float(m.group(1)) if m else 0.0
            i += len(m.group(0)) if m else 0
            return [d + bl for d in res], i
        i += len(re.match(r'[^(),:;]+', s[i:]).group(0))
        m = re.match(r':([0-9.eE+-]+)', s[i:])
        return [float(m.group(1)) if m else 0.0], i + (len(m.group(0)) if m else 0)
    return parse(0)[0]

ok = True
N = 200
out = os.path.join(tmp, "yule")
run(["--birth", "1", "--death", "0", "--size", str(N), "-s", str(N), "-n", "80", "--seed", "11", "--mutation-rate", "0", "--division-mutation-rate", "1"], out)
means = []
for l in open(os.path.join(out, "trees.nwk")):
    if l.strip():
        d = depths(l.strip().rstrip(';'))
        means.append(sum(d) / len(d))
theory = 2 * (sum(1.0 / k for k in range(1, N + 1)) - 1)
obs = sum(means) / len(means)
good = abs(obs / theory - 1) < 0.05
print("Yule depth: %.2f divisions to a cell, theory %.2f %s" % (obs, theory, "ok" if good else "FAIL"))
ok &= good

b, d, t0, T = 0.5, 0.1, 10.0, 16.0
for fit in (1.5, 2.0):
    out = os.path.join(tmp, "clone%g" % fit)
    run(["--birth", str(b), "--death", str(d), "--time", str(T), "--clone-time", str(t0), "--clone-fitness", str(fit), "-s", "20", "-n", "800", "--seed", "33"], out)
    tot = 0.0
    for l in open(os.path.join(out, "drivers.tsv")).read().split("\n")[1:]:
        f = l.split("\t")
        if len(f) > 6 and f[1] == "1":
            tot += float(f[6])
    obs = tot / 800
    theory = math.exp(fit * (b - d) * (T - t0))
    good = abs(obs / theory - 1) < 0.15
    print("clone of fitness %g: mean size %.1f, theory %.1f %s" % (fit, obs, theory, "ok" if good else "FAIL"))
    ok &= good
sys.exit(0 if ok else 1)
