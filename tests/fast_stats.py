#!/usr/bin/env python3
"""Statistical test of --fast: the distributions of the read data must not change.

Usage: fast_stats.py BIN [workdir]
Runs the same parameters and seeds with and without --fast (the tree, mutations and true genotypes are the same;
only the random numbers of the reads differ), reads the VCF (DP, read counts RC, genotype likelihoods G10) and compares
  * the distribution of DP and of the read count of each nucleotide (two-sample Kolmogorov-Smirnov test),
  * the mean of DP, of the number of error reads and of each log10 genotype likelihood (z test),
  * the fraction of cells whose maximum-likelihood genotype is the true one.
Thresholds are for p ~ 1e-4 so that a correct implementation essentially never fails, while a change in a
distribution of 1-2% (given the sample sizes) is detected. Needs numpy only."""
import math, os, subprocess, sys, tempfile
import numpy as np

BIN = os.path.abspath(sys.argv[1])
WORK = sys.argv[2] if len(sys.argv) > 2 else tempfile.mkdtemp()
CONFIGS = {
    "default":  "-n1 -s60 -l1500 -u1e-6 -C25 -E0.005 -A0.02 0.01 0 -D0.1 -1 -#{seed}",
    "doublets/dispersion": "-n1 -s60 -l1500 -u1e-6 -C12 -V8 -O3 -J5 -B0.1 0.005 -E0.01 -A0.05 0.02 0 -P0.1 0.01 -Q0.1 0.01 -#{seed}",
    "no amplification error": "-n1 -s60 -l1500 -u1e-6 -C6 -E0.01 -D0.1 -#{seed}",
    "high coverage": "-n1 -s30 -l800 -u1e-6 -C80 -E0.01 -A0.05 0.01 0 -#{seed}",
}
SEEDS = {"exact": [11, 12, 13, 14], "fast": [11, 12, 13, 14]}

KEYED = {}

def run(args, seed, mode):
    d = os.path.join(WORK, "%s_%d" % (mode, seed))
    os.makedirs(d, exist_ok=True)
    cmd = [BIN] + (["--fast"] if mode == "fast" else []) + args.format(seed=seed).split()
    env = dict(os.environ, OMP_NUM_THREADS="2")
    r = subprocess.run(cmd, cwd=d, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, env=env)
    if r.returncode != 0:
        print(r.stderr.decode()[-500:]); sys.exit("run failed: %s" % " ".join(cmd))
    dp, rc, g10, gt_ok = [], [], [], []
    with open(os.path.join(d, "Results", "vcf")) as f:
        for line in f:
            if line[0] == "#" or "\t" not in line:
                continue
            p = line.rstrip("\n").split("\t")
            fmt = p[8].split(":")
            iGT, iDP, iRC, iG = fmt.index("GT"), fmt.index("DP"), fmt.index("RC"), fmt.index("G10")
            for ci, cell in enumerate(p[9:]):
                x = cell.split(":")
                if len(x) < len(fmt):
                    continue
                KEYED[(mode, args, seed, p[2], ci)] = (x[iRC], x[iG])
                dp.append(int(x[iDP]))
                rc.append([int(v) for v in x[iRC].split(",")] if x[iRC] != "." else [0, 0, 0, 0])
                g = [float(v) if v != "." else -999.0 for v in x[iG].split(",")]
                g = (g + [-999.0] * 10)[:10]
                g10.append(g)
    return np.array(dp), np.array(rc), np.array(g10)

def ks_p(a, b):
    a = np.sort(a); b = np.sort(b)
    allv = np.concatenate([a, b])
    d = np.max(np.abs(np.searchsorted(a, allv, side="right") / len(a) - np.searchsorted(b, allv, side="right") / len(b)))
    en = math.sqrt(len(a) * len(b) / (len(a) + len(b)))
    lam = (en + 0.12 + 0.11 / en) * d
    p = 2 * sum((-1) ** (k - 1) * math.exp(-2 * k * k * lam * lam) for k in range(1, 100))
    return d, max(min(p, 1.0), 0.0)

def z_test(a, b):
    se = math.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))
    return 0.0 if se == 0 else (a.mean() - b.mean()) / se

fail = 0
for name, args in CONFIGS.items():
    data = {m: [np.concatenate(v) for v in zip(*[run(args, s, m) for s in SEEDS[m]])] for m in SEEDS}
    (dpE, rcE, gE), (dpF, rcF, gF) = data["exact"], data["fast"]
    print("%s: %d / %d cell-sites" % (name, len(dpE), len(dpF)))
    checks = []
    # same sites and cells in both modes; the two samples are positively correlated, so the unpaired tests are conservative
    checks.append(("DP", ks_p(dpE, dpF)[1], z_test(dpE.astype(float), dpF.astype(float))))
    errE = (rcE.sum(1) - rcE.max(1)).astype(float); errF = (rcF.sum(1) - rcF.max(1)).astype(float)
    checks.append(("non-major reads", ks_p(errE, errF)[1], z_test(errE, errF)))
    zeroE = (dpE == 0).astype(float); zeroF = (dpF == 0).astype(float)
    checks.append(("fraction DP=0", 1.0, z_test(zeroE, zeroF)))
    gmE = gE[dpE > 0].max(1); gmF = gF[dpF > 0].max(1)
    sE = np.sort(gE[dpE > 0], axis=1)[:, ::-1]; sF = np.sort(gF[dpF > 0], axis=1)[:, ::-1]
    gapE = sE[:, 0] - sE[:, 1]; gapF = sF[:, 0] - sF[:, 1]
    checks.append(("GL gap best-second", ks_p(gapE, gapF)[1], z_test(gapE, gapF)))
    for k in range(4):
        checks.append(("reads of the %d-th most frequent nucleotide" % (k + 1),
                       ks_p(np.sort(rcE, 1)[:, 3 - k], np.sort(rcF, 1)[:, 3 - k])[1],
                       z_test(np.sort(rcE, 1)[:, 3 - k].astype(float), np.sort(rcF, 1)[:, 3 - k].astype(float))))
    for what, p, z in checks:
        bad = p < 1e-4 or abs(z) > 4.0
        fail |= bad
        print("  %-45s KS p=%.3g  z=%+.2f  %s" % (what, p, z, "FAIL" if bad else "ok"))
# direct check of the vectorized likelihoods: cells and sites that happen to have the same read counts in both modes
# must have the same G10 values (up to the rounding to 2 decimals of the VCF). Only in the scenario without amplification
# error: with -A the error of every cell and site is random and differs between the modes.
same = worst = 0
for key, (rc, g) in KEYED.items():
    if key[0] != "fast" or rc == "." or "-A" in key[1]:
        continue
    other = KEYED.get(("exact",) + key[1:])
    if other is None or other[0] != rc or g == "." or other[1] == ".":
        continue
    a = [float(v) if v != "." else -999.0 for v in g.split(",")]
    b = [float(v) if v != "." else -999.0 for v in other[1].split(",")]
    same += 1
    worst = max(worst, max(abs(u - v) for u, v in zip(a, b)))
bad = same < 1000 or worst > 0.0101
fail |= bad
print("same read counts in both modes: %d cell-sites, largest difference in G10 = %.3f  %s" % (same, worst, "FAIL" if bad else "ok"))
sys.exit(1 if fail else 0)
