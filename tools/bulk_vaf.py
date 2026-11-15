#!/usr/bin/env python3
"""Bulk-like VAF distribution from a CellCoal SNV genotype file (option -1).

The variant allele frequency (VAF) of each SNV is the fraction of alternate alleles among all the alleles of the
sampled tumour cells (the healthy cell is the reference and is excluded). A bulk sequencing experiment is mimicked by
sampling reads at a given depth (binomial), optionally mixing in healthy cells (purity).

Usage: bulk_vaf.py snv_genotypes_file [--depth 100] [--purity 1.0] [--seed 1] [--out vaf.png] [--tsv vaf.tsv]
"""
import argparse, sys
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("file")
ap.add_argument("--depth", type=float, default=100, help="mean bulk depth (Poisson); 0 = no read sampling")
ap.add_argument("--purity", type=float, default=1.0, help="fraction of tumour cells in the bulk sample")
ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--out", default="vaf.png")
ap.add_argument("--tsv", default=None)
ap.add_argument("--bins", type=int, default=50)
a = ap.parse_args()

lines = [l.rstrip("\n") for l in open(a.file) if l.strip()]
i = 0
while len(lines[i].split()) != 2 or lines[i].startswith("["):   # skip the number of data sets and the [#n] label
    i += 1
ncells1, nsnv = map(int, lines[i].split()); i += 1
pos = lines[i].split(); i += 1
rows = lines[i:i + ncells1]
names = [r.split()[0] for r in rows]
geno = np.array([[g for g in r.split()[1:]] for r in rows])          # (cells+1) x SNVs, each "XY"
ref = geno[-1]                                                       # last row: healthy (outgroup) cell
tumour = geno[:-1]
alt = np.zeros(nsnv)
for j in range(nsnv):
    r = ref[j][0]
    col = tumour[:, j]
    alt[j] = sum((g[0] != r) + (g[1] != r) for g in col)
n_alleles = 2 * tumour.shape[0]
vaf_true = alt / n_alleles
vaf_bulk = a.purity * vaf_true                                       # healthy cells add reference alleles only
rng = np.random.default_rng(a.seed)
if a.depth > 0:
    dp = rng.poisson(a.depth, nsnv)
    vaf_obs = np.where(dp > 0, rng.binomial(dp, vaf_bulk) / np.maximum(dp, 1), np.nan)
else:
    vaf_obs = vaf_bulk

if a.tsv:
    with open(a.tsv, "w") as f:
        f.write("pos\tcells_alt_alleles\ttrue_vaf\tbulk_vaf\n")
        for j in range(nsnv):
            f.write(f"{pos[j]}\t{int(alt[j])}\t{vaf_true[j]:.6f}\t{vaf_obs[j]:.6f}\n")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
bins = np.linspace(0, 0.5 if a.purity <= 1 else 1, a.bins + 1)
ax[0].hist(vaf_true, bins=bins, color="#4c72b0")
ax[0].set_title(f"True VAF ({tumour.shape[0]} cells, {nsnv} SNVs)")
ax[1].hist(vaf_obs[~np.isnan(vaf_obs)], bins=bins, color="#dd8452")
ax[1].set_title(f"Bulk VAF (depth {a.depth:g}x, purity {a.purity:g})" if a.depth > 0 else "Bulk VAF (no read sampling)")
for x in ax:
    x.set_xlabel("VAF")
ax[0].set_ylabel("Number of SNVs")
plt.tight_layout()
plt.savefig(a.out, dpi=150)
print(f"{nsnv} SNVs; mean true VAF {vaf_true.mean():.4f}; plot written to {a.out}")
