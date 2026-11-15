#!/usr/bin/env python3
"""Biological validation of CellCoal: the simulated data are compared with what theory predicts.

Usage: theory.py BINARY [--quick]
These tests do not check that the code is stable (the regression tests do that) but that the model is right: they
fail if a change alters the distribution of the genealogies, the mutations or the reads. All the runs use fixed seeds,
so the result is reproducible; every statistic is a z-score against its standard error and must be below Z_MAX.

  1. Coalescent genealogy (constant size N): the waiting time with k lineages is exponential with mean 2N/(k(k-1))
     (a pair of cells coalesces at rate 1/N per generation), so the total tree length has mean 2N sum(1/i).
  2. Exponential growth (rate g): the time to the coalescence of two cells has the survival function
     exp(-(exp(g t)-1)/(g N)); its mean is integrated numerically.
  3. Mutations under the diploid infinite-sites model: E[S] = Theta sum(1/i), Theta = 4 N mu L, and the site frequency
     spectrum of the cells is E[xi_i] = Theta / i (Fu 1995), measured on the true SNV genotypes.
  4. Reads: mean and variance of the depth (each chromosome of a cell gets a negative binomial depth with mean C/2 and
     dispersion V), fraction of error reads (= sequencing error E) and fraction of heterozygous genotypes that
     lose one or both alleles by allelic dropout (-D).
Needs only the standard library."""
import math, os, random, re, statistics as st, subprocess, sys, tempfile

binary = os.path.abspath(sys.argv[1])
quick = "--quick" in sys.argv
Z_MAX = 4.0
fails = []

def run(d, *a):
    os.makedirs(d, exist_ok=True)
    p = subprocess.run([binary] + list(a), cwd=d, capture_output=True, text=True, timeout=1500)
    if p.returncode != 0:
        print(p.stdout[-800:], p.stderr[-800:]); sys.exit("simulation failed: " + " ".join(a))
    return p.stdout + p.stderr

def check(name, z, detail):
    ok = abs(z) < Z_MAX
    print("  %s  %-62s z=%+5.2f  %s" % ("ok  " if ok else "FAIL", name, z, detail))
    if not ok: fails.append(name)

def zmean(xs, theory):
    m = st.mean(xs); se = st.stdev(xs) / math.sqrt(len(xs))
    return (m - theory) / se, m

def read_times(path):
    """list of replicates: (sorted coalescence times of the ingroup, total ingroup branch length in generations)"""
    reps = []; cur = None
    for l in open(path):
        if l.startswith("--") and "Nodes" in l:
            cur = []; reps.append(cur); continue
        f = l.replace("|", " ").replace("(", " ").replace(")", " ").split()
        if cur is not None and len(f) == 9 and f[0] in ("tip", "internal", "ingroupMRCA"):
            cur.append((f[0], float(f[6]), float(f[7])))
    out = []
    for r in reps:
        out.append((sorted(t for c, t, _ in r if c != "tip"), sum(l for _, _, l in r)))
    return out

tmp = tempfile.mkdtemp(prefix="cellcoal_theory_")
R = 600 if quick else 3000

# ---------- 1. coalescent ----------
print("1. Coalescent genealogy, constant size")
N, n, u, L = 10000, 6, 1e-6, 1000
d = os.path.join(tmp, "coal")
log = run(d, "-n%d" % R, "-s%d" % n, "-l%d" % L, "-u%g" % u, "-e%d" % N, "-b1", "-k0", "-q0", "-c0", "-m0", "-7", "-#5")
reps = read_times(os.path.join(d, "Results", "times"))
assert len(reps) == R, len(reps)
wait = {k: [] for k in range(2, n + 1)}
for ts, _ in reps:
    prev = 0.0
    for i, t in enumerate(ts):
        wait[n - i].append(t - prev); prev = t
for k in range(2, n + 1):
    th = 2.0 * N / (k * (k - 1))
    z, m = zmean(wait[k], th)
    check("waiting time with k=%d lineages" % k, z, "mean %.0f, theory %.0f" % (m, th))
    cv = st.stdev(wait[k]) / m
    zcv = (cv - 1.0) / (1.0 / math.sqrt(2 * R) * 1.6)   # standard error of the CV of an exponential sample (about 1.6/sqrt(2R))
    check("  ... exponential (CV = 1) for k=%d" % k, zcv, "CV %.3f" % cv)
tot = [t for _, t in reps]
th = 2.0 * N * sum(1.0 / i for i in range(1, n))
z, m = zmean(tot, th)
check("total tree length", z, "mean %.0f, theory %.0f generations" % (m, th))

# ---------- 2. exponential growth ----------
print("2. Exponential growth (-g)")
Ng, g = 10000, 1.0e-4
alpha = g * Ng
# E[T2] = integral of exp(-(exp(g t)-1)/alpha) dt
def surv(t): return math.exp(-(math.exp(g * t) - 1.0) / alpha)
h, th2, t = 1.0, 0.0, 0.0
while surv(t) > 1e-14:
    th2 += 0.5 * h * (surv(t) + surv(t + h)); t += h
d = os.path.join(tmp, "growth")
run(d, "-n%d" % (2 * R), "-s2", "-l100", "-u1e-9", "-e%d" % Ng, "-g%g" % g, "-b1", "-k0", "-q0", "-c0", "-m0", "-7", "-#6")
t2 = [ts[0] for ts, _ in read_times(os.path.join(d, "Results", "times"))]
z, m = zmean(t2, th2)
check("time to coalescence of two cells, N=%d g=%g" % (Ng, g), z, "mean %.0f, theory %.0f" % (m, th2))

# ---------- 3. mutations ----------
print("3. Mutations, diploid infinite-sites model")
theta = 4.0 * N * u * L
m_log = float(re.search(r"Mean number of mutation events\s*=\s*([\d.]+)", log).group(1))
v_log = float(re.search(r"Variance number of mutation events\s*=\s*([\d.]+)", log).group(1))
a_n = sum(1.0 / i for i in range(1, n))
check("mean number of mutations (Theta sum 1/i)", (m_log - theta * a_n) / math.sqrt(v_log / R), "%.2f vs %.2f (Theta=%g)" % (m_log, theta * a_n, theta))

def sfs(d, n, theta):
    """per-replicate counts of SNVs carried by i of the n cells (true genotypes; the outgroup is the reference)"""
    counts = []; lines = open(os.path.join(d, "Results", "snv_gen")).read().split("\n")
    i = 0
    while i < len(lines):
        if lines[i].startswith("[#"):
            ncell, nsnv = map(int, lines[i + 1].split())
            rows = lines[i + 3:i + 3 + ncell] if nsnv > 0 else []
            c = [0] * (n + 1)
            if nsnv > 0:
                geno = [r.split()[1:] for r in rows]; ref = geno[-1]   # outgroup (last row) = ancestral genotypes
                for s in range(nsnv):
                    c[sum(1 for x in range(n) if geno[x][s] != ref[s])] += 1
            counts.append(c); i += 3 + ncell
        else:
            i += 1
    return counts
d = os.path.join(tmp, "sfs")
run(d, "-n%d" % R, "-s%d" % n, "-l%d" % L, "-u%g" % u, "-e%d" % N, "-b1", "-k0", "-q0", "-c0", "-m0", "-1", "-#7")
cnt = sfs(d, n, theta)
for i in range(1, n):
    z, m = zmean([c[i] for c in cnt], theta / i)
    check("site frequency spectrum: SNVs in %d of %d cells (Theta/i)" % (i, n), z, "mean %.2f, theory %.2f" % (m, theta / i))

# ---------- 4. reads ----------
print("4. Reads, sequencing error and allelic dropout")
def bulk(path):
    rows = [l.rstrip("\n").split("\t") for l in open(path)]
    return [r for r in rows if len(r) == 11 and r[0] != "site"]
# The bulk file lists only the sites with alternate reads or true mutations. With a high error rate (0.2) practically
# every site has alternate reads, so all the sites are listed and the depth is not biased by the selection.
cells, C, V = 20, 10.0, 5.0
nsites = 40000 if not quick else 12000
d = os.path.join(tmp, "reads")
run(d, "-n1", "-s%d" % cells, "-l%d" % nsites, "-u1e-8", "-e10000", "-b1", "-C%g" % C, "-V%g" % V, "-E0.2", "-0", "-Z", "-#4")
rows = bulk(os.path.join(d, "Results", "bulk"))
assert len(rows) > 0.999 * nsites, len(rows)
dep = [int(r[6]) for r in rows]
z, m = zmean(dep, cells * C)
check("pooled depth per site, mean (tumour cells x C)", z, "%.2f vs %.2f" % (m, cells * C))
# each of the two chromosomes of a cell gets a negative binomial depth with mean C/2 and dispersion V
var_th = cells * 2.0 * (C / 2 + (C / 2) ** 2 / V)
var = st.variance(dep)
check("pooled depth per site, variance (cells x (C + C^2/2V))", (var - var_th) / (var_th * math.sqrt(2.0 / (len(dep) - 1))), "%.1f vs %.1f" % (var, var_th))
alt = sum(int(r[7]) for r in rows); tdep = sum(dep)
p = alt / tdep
check("fraction of non-reference reads at E=0.2", (p - 0.2) / math.sqrt(0.2 * 0.8 / tdep), "%.5f vs 0.2" % p)

# low error rate: the sites listed have at least one alternate read; maximum likelihood of the zero-truncated binomial
E = 0.01
d = os.path.join(tmp, "err")
run(d, "-n1", "-s%d" % cells, "-l%d" % (60000 if not quick else 20000), "-u1e-8", "-e10000", "-b1", "-C%g" % C, "-V%g" % V, "-E%g" % E, "-0", "-Z", "-#5")
rows = [r for r in bulk(os.path.join(d, "Results", "bulk")) if r[9] == "0"]    # not mutated: every alternate read is an error
data = [(int(r[6]), int(r[7])) for r in rows]
def ll(e): return sum(k * math.log(e) + (D_ - k) * math.log(1 - e) - math.log(1 - (1 - e) ** D_) for D_, k in data)
lo, hi = 1e-4, 0.2
for _ in range(80):
    a1, a2 = lo + (hi - lo) / 3, hi - (hi - lo) / 3
    if ll(a1) < ll(a2): lo = a1
    else: hi = a2
e_hat = 0.5 * (lo + hi); hstep = 1e-5
info = -(ll(e_hat + hstep) - 2 * ll(e_hat) + ll(e_hat - hstep)) / hstep ** 2
check("sequencing error E, from the non-mutated sites (truncated ML)", (e_hat - E) * math.sqrt(info), "%.5f vs %.5f" % (e_hat, E))

# allelic dropout: a genotype loses at least one allele with probability D (no reads: observed genotypes, lost alleles are '?')
D = 0.3
d = os.path.join(tmp, "ado")
run(d, "-n%d" % (15 if quick else 40), "-s30", "-l2000", "-u2e-7", "-e10000", "-b1", "-D%g" % D, "-3", "-Y", "-#8")
def true_haps(path):
    out = []; cur = {}
    for l in open(path):
        f = l.split()
        if l.startswith("[#"):
            if cur: out.append(cur)
            cur = {}
        elif len(f) == 2 and f[0].startswith("cell") and f[0][-1] in "mp":
            cur.setdefault(f[0][:-1], {})[f[0][-1]] = f[1]
    if cur: out.append(cur)
    return out
def full_gen(path):
    out = []; cur = None
    for l in open(path):
        f = l.split()
        if l.startswith("[#"):
            cur = {}; out.append(cur)
        elif len(f) > 2 and f[0].startswith("cell"):
            cur[f[0]] = f[1:]
    return out
th = true_haps(os.path.join(d, "Results", "true_hap")); fg = full_gen(os.path.join(d, "Results", "full_gen"))
het = lost = both = 0
for T, G in zip(th, fg):
    for cname, hp in T.items():
        for j, (a, b) in enumerate(zip(hp["m"], hp["p"])):
            if a != b:
                het += 1
                q = G[cname][j].count("?")          # an allele lost by dropout is written "?"
                lost += q > 0; both += q == 2
p = lost / het
check("heterozygous genotypes that lose an allele (= D)", (p - D) / math.sqrt(D * (1 - D) / het), "%.4f vs %.2f (%d genotypes)" % (p, D, het))
pb = (1 - math.sqrt(1 - D)) ** 2
p = both / het
check("... and lose both alleles (= (1 - sqrt(1-D))^2)", (p - pb) / math.sqrt(pb * (1 - pb) / het), "%.4f vs %.4f" % (p, pb))

print()
if fails:
    print("THEORY: %d check(s) failed: %s" % (len(fails), "; ".join(fails))); sys.exit(1)
print("THEORY ok: all the checks agree with the theory (|z| < %.0f)" % Z_MAX)
