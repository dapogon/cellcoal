#!/usr/bin/env python3
"""Check of the clones model (--clones) against a forward simulation of the whole tumor.

CellCoal draws the genealogy of a sample of cells from a tumor made of clones of different growth rate with a coalescent in which the size of every
clone grows, N_k(t) = N_k(0) exp(-g_k t) looking back. Here the tumor itself is simulated FORWARD in time, generation by generation, as a haploid
population in which the cells of a generation choose their parent at random among the cells of the previous generation of their clone: the root clone
starts from one cell and grows at the rate g, a clone k is born at its time from founder cells of its parent clone and grows at its own rate. At the end
S cells are sampled at random from the whole tumor, with no knowledge of the clones, and their genealogy is read from the parent of every cell.
CellCoal and the forward simulation are compared in the frequency of every clone in the sample and in the mean time to the coalescence of two cells,
for every pair of clones, over many replicates. Two forward models are used:

  wf    non-overlapping generations (a Wright-Fisher population whose size in every clone is N_k(t)): the coalescence rate of two cells is 1/N_k(t),
        which is -K0.5 -L0.496 here (2 b = 1); it checks the genealogy of the clones, the births and the numerical solution of the times;
  yule  a tumor in which the cells divide and none dies (a pure birth process, birth=g, death=0), simulated by cellforward, the forward program of
        CellCoal: every cell divides with probability g_k per generation, the sizes of the clones are random (as in a real tumor; only the tumors of the
        size given to CellCoal are kept) and the coalescence rate of two cells is 2 g_k / N_k(t).

usage: clones_check.py BIN [replicates] [wf|yule] [CELLFORWARD]  (the yule mode finds cellforward in the folder of BIN)
"""
import glob
import itertools
import math
import os
import re
import subprocess
import sys
import tempfile

import numpy as np

BIN = sys.argv[1]
REPS = int(sys.argv[2]) if len(sys.argv) > 2 else 800
MODE = sys.argv[3] if len(sys.argv) > 3 else "wf"
N, G, S = 500, 4e-3, 8          # cells of the clone 1 today (the whole tumor has these plus the other clones), its growth rate (per generation), cells sampled
TIMES = [2.0, 1.0]              # birth of the clones 2 and 3 (units of N generations here; CellCoal gets generations)
PARENT = [0, 1]                 # parent clones (0-based)
FIT = [1.0, 1.5]                # growth rate relative to the clone 1
FOUNDERS = [1, 1]
K = 3


def cellcoal(tmp, reps):
    model = "birth=0.5,death=%g" % (0.5 - G) if MODE == "wf" else "birth=%g,death=0" % G
    others = sum(FOUNDERS[k - 1] * math.exp(G * FIT[k - 1] * TIMES[k - 1] * N) for k in (1, 2))
    cmd = [BIN, "-n%d" % reps, "-#5", "-l10", "-u1e-9", "-s%d" % S, "-7"] + [
           "--clones", "n=3,cells=%d,%s,sizes=%s,time=%s,parent=%s,fitness=%s,founders=%s" % (
               round(N + others), model, "fixed" if MODE == "wf" else "random", ":".join("%g" % (t * N) for t in TIMES), ":".join(str(p + 1) for p in PARENT),
               ":".join(map(str, FIT)), ":".join(map(str, FOUNDERS))), "-o" + tmp]
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    blocks = open(os.path.join(tmp, "times")).read().split('-------------------- Nodes -------------------')[1:]
    clone = {}
    for line in open(os.path.join(tmp, "clones.tsv")).read().splitlines()[1:]:
        r, c, k = line.split("\t")
        clone.setdefault(int(r), []).append(int(k) - 1)
    out = {}
    for r, blk in enumerate(blocks, 1):
        node = {}
        for line in blk.splitlines():
            m = re.match(r'\s*(tip|internal|ingroupMRCA)\s+\d+\s+(\d+)\s+\(\s*(-?\d+)\s+(-?\d+)\s+(-?\d+)\)\s+\|\s+([-\d.]+)', line)
            if m:
                node[int(m.group(2))] = (int(m.group(5)), float(m.group(6)))

        def path(i):
            p = [i]
            while node[p[-1]][0] in node and node[p[-1]][0] != -1:
                p.append(node[p[-1]][0])
            return p
        for i, j in itertools.combinations(range(1, S + 1), 2):
            pj = set(path(j))
            t = next(node[x][1] for x in path(i) if x in pj)
            key = tuple(sorted((clone[r][i - 1], clone[r][j - 1])))
            out.setdefault(key, []).append(t)
    freq = np.bincount([k for r in clone for k in clone[r]], minlength=K) / float(sum(len(v) for v in clone.values()))
    return out, freq


def forward(rng, gr, cells):
    """one tumor, forward in time; returns the clones of the S sampled cells and the times (generations) to the coalescence of every pair"""
    tb = [0] + [int(round(t * N)) for t in TIMES]               # birth, generations before the present
    start = int(max(math.log(cells[0]) / gr[0], max(tb))) + 2   # the clone 1 has one cell at the start
    size = np.zeros((K, start + 1), dtype=int)                    # size[k][g]: cells of the clone k in the generation g (g = start: the present)
    for k in range(K):
        for g in range(start + 1):
            back = start - g
            if k > 0 and back > tb[k]:
                continue
            size[k, g] = max(FOUNDERS[k - 1] if k > 0 else 1, int(round(cells[k] * math.exp(-gr[k] * back))))
    par = [[None] * (start + 1) for _ in range(K)]               # parents: indexes in the previous generation of the clone (or of its parent clone)
    first = [0] * K
    for k in range(1, K):
        first[k] = start - tb[k]
    for g in range(1, start + 1):
        for k in range(K):
            if size[k, g] == 0:
                continue
            if k > 0 and g == first[k]:
                p = PARENT[k - 1]
                par[k][g] = rng.choice(size[p, g - 1], size=size[k, g], replace=False)
            else:
                par[k][g] = rng.integers(0, size[k, g - 1], size=size[k, g])
    # sample S cells of the whole tumor
    total = size[:, start].sum()
    pick = rng.choice(total, size=S, replace=False)
    edges = np.cumsum(size[:, start])
    lin = []
    for x in pick:
        k = int(np.searchsorted(edges, x, side='right'))
        lin.append((k, int(x - (edges[k - 1] if k else 0))))
    clones = [l[0] for l in lin]
    anc = {i: lin[i] for i in range(S)}
    members = {lin[i]: {i} for i in range(S)}
    times = {}
    cur = dict(members)                 # lineage (clone, index) -> cells below it
    g = start
    while len(cur) > 1 and g > 0:
        nxt = {}
        for (k, idx), cs in cur.items():
            kk, ii = k, idx
            if k > 0 and g == first[k]:
                kk = PARENT[k - 1]
                ii = int(par[k][g][idx])
            else:
                ii = int(par[k][g][idx])
            key = (kk, ii)
            if key in nxt:
                for a in nxt[key]:
                    for b in cs:
                        times[(min(a, b), max(a, b))] = start - (g - 1)
                nxt[key] = nxt[key] | cs
            else:
                nxt[key] = set(cs)
        cur = nxt
        g -= 1
    return clones, times


def newick_pair_times(nwk, n):
    """time to the coalescence of every pair of cells (cell0001 ... celln) of a tree in Newick with branch lengths in generations"""
    parent, length, name = [-1], [0.0], [None]
    cur, i = 0, 0
    while i < len(nwk) and nwk[i] != ";":
        c = nwk[i]
        if c == "(":
            parent.append(cur); length.append(0.0); name.append(None)
            cur = len(parent) - 1
        elif c == ",":
            parent.append(parent[cur]); length.append(0.0); name.append(None)
            cur = len(parent) - 1
        elif c == ")":
            cur = parent[cur]
        elif c == ":":
            j = i + 1
            while nwk[j] not in ",);":
                j += 1
            length[cur] = float(nwk[i + 1:j])
            i = j - 1
        else:
            j = i
            while nwk[j] not in ",):":
                j += 1
            name[cur] = nwk[i:j]
            i = j - 1
        i += 1
    depth = [0.0] * len(parent)             # distance to the root
    for v in range(1, len(parent)):         # parents are numbered before their children
        depth[v] = depth[parent[v]] + length[v]
    leaf = {int(name[v][4:]) - 1: v for v in range(len(parent)) if name[v]}
    anc = {}
    for c, v in leaf.items():
        chain = []
        while v >= 0:
            chain.append(v)
            v = parent[v]
        anc[c] = chain
    times = {}
    end = max(depth[v] for v in leaf.values())
    for x in range(n):
        for y in range(x + 1, n):
            sx = set(anc[x])
            m = next(v for v in anc[y] if v in sx)      # the nearest common ancestor
            times[(x, y)] = end - depth[m]
    return times


def forward_replicates(tmp, gr, cells):
    """tumors simulated forward by cellforward (pure birth, clones started at the times of CellCoal); only the tumors of the size given to CellCoal are kept"""
    fwd = sys.argv[4] if len(sys.argv) > 4 else sorted(glob.glob(os.path.join(os.path.dirname(os.path.abspath(BIN)), "cellforward-*")))[0]
    tend = math.log(N) / G                  # the clone 1 has N cells on average at the end
    out = os.path.join(tmp, "fwd")
    cmd = [fwd, "--birth", "%g" % G, "--time", "%.6f" % tend, "--clone-time", ":".join("%.6f" % (tend - t * N) for t in TIMES),
           "--clone-fitness", ":".join(map(str, FIT)), "--clone-parent", ":".join(map(str, PARENT)), "-s", str(S), "-n", str(40 * REPS),
           "--seed", "11", "-o", out]
    subprocess.run(cmd, check=True, stderr=subprocess.DEVNULL)
    tot = sum(cells)
    keep = {}
    for line in open(os.path.join(out, "summary.tsv")).read().splitlines()[1:]:
        f = line.split("\t")
        if abs(float(f[3]) - tot) <= 0.1 * tot:
            keep[int(f[0])] = None
    clones = {}
    for line in open(os.path.join(out, "clones.tsv")).read().splitlines()[1:]:
        r, c, k = line.split("\t")[:3]
        if int(r) in keep:
            clones.setdefault(int(r), {})[int(c[4:]) - 1] = int(k)
    trees = open(os.path.join(out, "trees.nwk")).read().splitlines()
    res = []
    for r in sorted(clones)[:REPS]:
        lab = [clones[r][i] for i in range(S)]
        res.append((lab, newick_pair_times(trees[r - 1], S)))
    return res


def main():
    gr = [G, G * FIT[0], G * FIT[1]]
    cells = [float(N)]
    for k in (1, 2):
        cells.append(FOUNDERS[k - 1] * math.exp(gr[k] * TIMES[k - 1] * N))
    tot = sum(cells)
    with tempfile.TemporaryDirectory() as tmp:
        a, freqA = cellcoal(tmp, 8 * REPS)
    rng = np.random.default_rng(11)
    b = {}
    cnt = np.zeros(K)
    if MODE == "wf":
        results = []
        tries = 0
        while len(results) < REPS and tries < 60 * REPS:
            tries += 1
            res = forward(rng, gr, cells)
            if res is not None:
                results.append(res)
    else:
        with tempfile.TemporaryDirectory() as tmp:
            results = forward_replicates(tmp, gr, cells)
    for clones, times in results:
        for k in clones:
            cnt[k] += 1
        for (i, j), t in times.items():
            b.setdefault(tuple(sorted((clones[i], clones[j]))), []).append(t)
    print("frequency of the clones: expected %s, CellCoal %s, forward %s" % (
        [round(c / tot, 3) for c in cells], [round(x, 3) for x in freqA], [round(x, 3) for x in cnt / cnt.sum()]))
    ok = True
    for key in sorted(a):
        if key not in b:
            continue
        ma, mb = np.mean(a[key]), np.mean(b[key])
        se = math.sqrt(np.var(a[key]) / len(a[key]) + np.var(b[key]) / len(b[key]))
        print("pair of clones %s: CellCoal %8.1f (%d pairs)  forward %8.1f (%d pairs)  ratio %.3f  z %.1f" % (
            (key[0] + 1, key[1] + 1), ma, len(a[key]), mb, len(b[key]), ma / mb, (ma - mb) / se))
        if len(a[key]) >= 300 and len(b[key]) >= 300 and abs(ma / mb - 1.0) > (0.12 if MODE == "wf" else 0.15):
            ok = False
    for k in range(K):
        if abs(freqA[k] - cells[k] / tot) > 0.03 or abs(cnt[k] / cnt.sum() - cells[k] / tot) > 0.03:
            ok = False
    if MODE != "wf":
        print("INFO: pure birth tumor with random sizes: the model leaves out the stochastic first generations of every clone, so the times are a few percent deeper in CellCoal")
    print("OK" if ok else "FAILED")
    sys.exit(0 if ok else 1)


main()
