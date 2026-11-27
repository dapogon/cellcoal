#!/usr/bin/env python3
"""cellforward against tumopp (Iwasaki and Innan 2017; https://github.com/heavywatal/tumopp), a lattice simulator.
   With the resource competition 'const' (-L const, the default) the birth and death rates of a cell do not depend on its neighbors,
   so the genealogy of the cells is that of a well-mixed birth-death process (the lattice only places the cells): it has to agree with
   cellforward. Compared: the time to reach N cells, the mean time since the coalescence of a pair of 20 sampled cells, and the mean
   number of divisions from the first cell to a sampled cell. -k of tumopp and --cycle-shape of cellforward are the same shape.
   usage: tumopp_compare.py TUMOPP CELLFORWARD [N=2000] [replicates=80]"""
import gzip, math, os, random, subprocess, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from newick_tools import newick_pair_times


TUMOPP, CF = sys.argv[1], sys.argv[2]
N = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
REPS = int(sys.argv[4]) if len(sys.argv) > 4 else 80
S = 20
tmp = tempfile.mkdtemp()

def tumopp(b, d, k, seed):
    out = os.path.join(tmp, "t%d" % seed)
    subprocess.run([TUMOPP, "-D3", "-Cmoore", "-N%d" % N, "-b%g" % b, "-d%g" % d, "-k%g" % k, "--seed", str(seed), "-o", out], capture_output=True, check=True)
    anc, born, died, alive = {}, {}, {}, []
    for l in gzip.open(os.path.join(out, "population.tsv.gz"), "rt").read().split("\n")[1:]:
        f = l.split("\t")
        if len(f) < 8:
            continue
        i = int(f[3])
        anc[i], born[i], died[i] = int(f[4]), float(f[5]), float(f[6])
        if died[i] == 0.0:
            alive.append(i)
    end = max(born.values())
    rng = random.Random(seed)
    sample = rng.sample(alive, S)
    chains = {}
    for c in sample:
        ch, v = [], c
        while v in anc:
            ch.append(v)
            v = anc[v]
        chains[c] = ch
    pairs, depth = [], []
    for a in range(S):
        for c in range(a):
            sa = chains[sample[a]]
            sc = chains[sample[c]]
            seta = set(sa)
            mrca = next(v for v in sc if v in seta)
            pairs.append(end - died[mrca])
            depth.append(sa.index(mrca) + sc.index(mrca))      # divisions of the two lineages from the MRCA
    return end, sum(pairs) / len(pairs), sum(depth) / len(depth)

def cellforward(b, d, k, seed):
    res = []
    for mode in (["--mutation-rate", "1", "--division-mutation-rate", "0"], ["--mutation-rate", "0", "--division-mutation-rate", "1"]):
        out = os.path.join(tmp, "f%d" % seed)
        subprocess.run([CF, "--birth", str(b), "--death", str(d), "--size", str(N), "-s", str(S), "-n", str(REPS), "--seed", str(seed), "--cycle-shape", str(k), "-o", out] + mode, capture_output=True, check=True)
        res.append((open(os.path.join(out, "summary.tsv")).read().split("\n")[1:], open(os.path.join(out, "trees.nwk")).read().split("\n")))
    times = [float(l.split("\t")[2]) for l in res[0][0] if l]
    pt, dp = [], []
    for nwk in res[0][1]:
        if nwk.strip():
            v = list(newick_pair_times(nwk.strip(), S).values())
            pt.append(sum(v) / len(v))
    for nwk in res[1][1]:
        if nwk.strip():
            v = list(newick_pair_times(nwk.strip(), S, True).values())
            dp.append(sum(v) / len(v))      # divisions of the two lineages from their MRCA (the branch lengths are divisions here)
    return times, pt, dp

def stat(x):
    m = sum(x) / len(x)
    return m, math.sqrt(sum((v - m) ** 2 for v in x) / (len(x) - 1) / len(x))

ok = True
print("N=%d cells, %d replicates, %d sampled cells" % (N, REPS, S))
print("%-24s %-18s %-22s %-22s" % ("b, d, k", "time to N", "pair time (since)", "divisions of a pair"))
for b, d, k in ((0.5, 0.1, 1), (0.5, 0.1, 4), (0.5, 0.4, 1)):
    T, P, D = [], [], []
    for seed in range(1, REPS + 1):
        t, p, dd = tumopp(b, d, k, seed)
        T.append(t); P.append(p); D.append(dd)
    ft, fp, fd = cellforward(b, d, k, 7)
    line = []
    for name, x, y in (("time", T, ft), ("pair", P, fp), ("div", D, fd)):
        (m1, e1), (m2, e2) = stat(x), stat(y)
        z = (m1 - m2) / math.sqrt(e1 * e1 + e2 * e2)
        good = abs(z) < 3 and abs(m1 / m2 - 1) < 0.08
        ok &= good
        line.append("%.2f / %.2f (z=%+.1f)%s" % (m1, m2, z, "" if good else " FAIL"))
    print("%-24s %-18s %-22s %-22s" % ("%g, %g, %g" % (b, d, k), line[0], line[1], line[2]))
print("(tumopp / cellforward; z = difference in standard errors)")
sys.exit(0 if ok else 1)
