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
  yule  a tumor in which the cells divide and none dies (a pure birth process, -K g -L0): every cell divides with probability g_k per
        generation, the sizes of the clones are random (as in a real tumor) and the coalescence rate of two cells is 2 g_k / N_k(t).

usage: clones_check.py BIN [replicates] [wf|yule]
"""
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
           "--clones", "n=3,cells=%d,%s,time=%s,parent=%s,fitness=%s,founders=%s" % (
               round(N + others), model, ":".join("%g" % (t * N) for t in TIMES), ":".join(str(p + 1) for p in PARENT),
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


def forward_yule(rng, gr, cells):
    """a tumor in which every cell divides with probability g_k per generation and none dies (pure birth); the clone 1 starts from one cell"""
    tb = [0] + [int(round(t * N)) for t in TIMES]
    start = int(math.log(cells[0]) / gr[0]) + 1
    first = [0] + [start - tb[k] for k in (1, 2)]
    par = [[None] * (start + 1) for _ in range(K)]
    cur = [1, 0, 0]
    for g in range(1, start + 1):
        new = cur[:]
        for k in range(K):
            if k > 0 and g < first[k]:
                continue
            if k > 0 and g == first[k]:
                par[k][g] = np.array([rng.integers(cur[PARENT[k - 1]])])
                new[k] = 1
                continue
            nk = cur[k]
            b = rng.binomial(nk, gr[k])
            mothers = rng.choice(nk, b, replace=False) if b > 0 else np.zeros(0, dtype=int)
            par[k][g] = np.concatenate([np.arange(nk), mothers]).astype(int)
            new[k] = nk + b
        cur = new
    total = sum(cur)
    pick = rng.choice(total, size=S, replace=False)
    edges = np.cumsum(cur)
    lin = []
    for x in pick:
        k = int(np.searchsorted(edges, x, side='right'))
        lin.append((k, int(x - (edges[k - 1] if k else 0))))
    clones = [l[0] for l in lin]
    curm = {lin[i]: {i} for i in range(S)}
    times = {}
    g = start
    while len(curm) > 1 and g > 0:
        nxt = {}
        for (k, idx), cs in curm.items():
            kk = PARENT[k - 1] if (k > 0 and g == first[k]) else k
            key = (kk, int(par[k][g][idx]))
            if key in nxt:
                for a in nxt[key]:
                    for b_ in cs:
                        times[(min(a, b_), max(a, b_))] = start - (g - 1)
                nxt[key] = nxt[key] | cs
            else:
                nxt[key] = set(cs)
        curm = nxt
        g -= 1
    return clones, times


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
    for _ in range(REPS):
        clones, times = (forward if MODE == "wf" else forward_yule)(rng, gr, cells)
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
        if len(a[key]) >= 300 and len(b[key]) >= 300 and abs(ma / mb - 1.0) > (0.12 if MODE == "wf" else 0.20):
            ok = False
    for k in range(K):
        if abs(freqA[k] - cells[k] / tot) > 0.03 or abs(cnt[k] / cnt.sum() - cells[k] / tot) > (0.03 if MODE == "wf" else 0.06):
            ok = False
    if MODE != "wf":
        print("INFO: the sizes of the clones are random in this tumor (a pure birth process from one cell), but fixed in CellCoal: not a pass/fail test")
        sys.exit(0)
    print("OK" if ok else "FAILED")
    sys.exit(0 if ok else 1)


main()
