#!/usr/bin/env python3
"""Summary plots of CellCoal output.

  cellcoal_plots.py depth   Results/vcf_dir/vcf.0001 [--out depth.png]     read depth per cell, per site and per genotype
  cellcoal_plots.py missing Results/vcf_dir/vcf.0001 [--out missing.png]    genotypes without reads: per cell and per site
  cellcoal_plots.py heatmap Results/snv_genotypes_dir/snv_gen.0001 [--out heatmap.png] [--max-sites 2000]
                                                                          cells x SNVs, coloured by the alleles that differ from the healthy cell
  cellcoal_plots.py all     Results [--replicate 1] [--outdir plots]       every plot that the files of a Results folder allow,
                                                                          plus the bulk VAF, and a self-contained report.html

The VCF can be plain or gzip (-z). 'depth' and 'missing' need reads (option -C); 'heatmap' needs the SNV genotype file (-1).
Every command prints a short summary in the terminal and writes a PNG (and, with --tsv, the numbers behind the plot).
"""
import argparse, gzip, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_vcf_depth(path, with_healthy=False):
    """returns (cell names, site positions, matrix sites x cells of DP); the healthy cell (last column) is dropped by default"""
    op = gzip.open if str(path).endswith(".gz") else open
    names, pos, rows = None, [], []
    with op(path, "rt", errors="ignore") as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                names = line.rstrip("\n").split("\t")[9:]
            elif line[:1].isdigit() and "\t" in line:
                f = line.rstrip("\n").split("\t")
                if len(f) < 10:
                    continue
                k = f[8].split(":").index("DP")
                pos.append(int(f[1]))
                rows.append([int(x.split(":")[k]) if x not in (".", "") else 0 for x in f[9:]])
    if names is None or not rows:
        sys.exit("no sites with reads were found in %s (are reads simulated, option -C, and is it the VCF?)" % path)
    d = np.array(rows, dtype=float)
    if not with_healthy:
        names, d = names[:-1], d[:, :-1]
    return names, np.array(pos), d


def poisson_pmf(k, mu):
    from math import lgamma, log, exp
    return np.array([exp(-mu + x * log(mu) - lgamma(x + 1)) for x in k]) if mu > 0 else np.zeros(len(k))


def cmd_depth(a):
    names, pos, d = read_vcf_depth(a.file)
    per_cell, per_site = d.mean(axis=0), d.sum(axis=1)
    cv2 = lambda x: x.var() / x.mean() ** 2
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    ax[0].hist(per_cell, bins=a.bins, color="#4c78a8")
    ax[0].set(title="Mean depth per cell", xlabel="reads per site", ylabel="cells")
    ax[1].hist(per_site, bins=a.bins, color="#f58518")
    ax[1].set(title="Pooled depth per site", xlabel="reads of all cells", ylabel="sites")
    v = d.ravel().astype(int)
    kmax = int(np.percentile(v, 99.5)) + 1
    ax[2].hist(v, bins=np.arange(kmax + 2) - 0.5, density=True, color="#54a24b", label="simulated")
    ks = np.arange(kmax + 1)
    ax[2].plot(ks, poisson_pmf(ks, v.mean()), "k.-", lw=1, label="Poisson, same mean")
    ax[2].set(title="Depth of a genotype", xlabel="reads", ylabel="fraction"); ax[2].legend(frameon=False)
    for x in ax:
        x.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(a.out, dpi=150)
    print("%d sites, %d cells (healthy cell excluded)" % d.shape)
    print("mean depth %.2f; zero-read genotypes %.1f%%" % (v.mean(), 100 * (v == 0).mean()))
    print("CV^2 between cells (mean depth)  %.3f   (about 1/alpha with -J)" % cv2(per_cell))
    print("CV^2 between sites (pooled depth) %.3f   (about 1/alpha with -O)" % cv2(per_site))
    print("variance/mean of a genotype's depth %.2f   (1 = Poisson; larger with -V, -O, -J)" % (v.var() / v.mean()))
    if a.tsv:
        with open(a.tsv, "w") as o:
            o.write("cell\tmean_depth\tzero_fraction\n")
            for n, m, z in zip(names, per_cell, (d == 0).mean(axis=0)):
                o.write("%s\t%.4f\t%.4f\n" % (n, m, z))
    print("wrote", a.out)


def cmd_missing(a):
    names, pos, d = read_vcf_depth(a.file)
    zero = d == 0
    per_cell, per_site = zero.mean(axis=0), zero.mean(axis=1)
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    ax[0].plot(np.sort(per_cell)[::-1], color="#4c78a8")
    ax[0].set(title="Genotypes without reads, per cell", xlabel="cells (sorted)", ylabel="fraction of sites", ylim=(0, max(0.05, per_cell.max() * 1.1)))
    ax[1].hist(per_site, bins=a.bins, color="#f58518")
    ax[1].set(title="Cells without reads, per site", xlabel="fraction of cells", ylabel="sites")
    ax[2].plot(np.sort(per_site)[::-1], color="#e45756")
    ax[2].set(title="Sites sorted by missing cells", xlabel="sites (sorted)", ylabel="fraction of cells", ylim=(0, 1))
    for x in ax:
        x.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(a.out, dpi=150)
    print("%d sites, %d cells" % d.shape)
    print("genotypes without reads: %.2f%% (cells: mean %.2f%%, max %.2f%%)" % (100 * zero.mean(), 100 * per_cell.mean(), 100 * per_cell.max()))
    print("sites with more than half of the cells without reads: %d" % int((per_site > 0.5).sum()))
    if a.tsv:
        with open(a.tsv, "w") as o:
            o.write("cell\tmissing_fraction\n")
            for n, z in zip(names, per_cell):
                o.write("%s\t%.4f\n" % (n, z))
    print("wrote", a.out)


def read_snv_gen(path):
    """first replicate of an SNV genotype file: (names, sites, list of genotype token lists)"""
    names, rows, sites = [], [], None
    with open(path, errors="ignore") as fh:
        lines = [l.rstrip("\n") for l in fh]
    i = 0
    while i < len(lines) and not lines[i].strip().split()[:1] or (i < len(lines) and lines[i].startswith("[#")):
        i += 1
    n, m = [int(x) for x in lines[i].split()[:2]]
    sites = [int(x) for x in lines[i + 1].split()]
    for l in lines[i + 2:i + 2 + n]:
        f = l.split()
        names.append(f[0]); rows.append(f[1:])
    if not rows or len(rows[0]) != m:
        sys.exit("unexpected format in %s: is it the SNV genotype file (option -1)?" % path)
    return names, sites, rows


def cmd_heatmap(a):
    names, sites, rows = read_snv_gen(a.file)
    ref = rows[-1]                          # the healthy (outgroup) cell
    ncell, nsnv = len(rows) - 1, len(sites)
    code = np.zeros((ncell, nsnv), dtype=np.int8)   # 0 = as healthy, 1 = one allele differs, 2 = both, 3 = missing
    for i in range(ncell):
        for j in range(nsnv):
            g, r = rows[i][j], ref[j]
            if len(g) != 2 or any(c not in "ACGT" for c in g):
                code[i, j] = 3
            else:
                code[i, j] = (g[0] != r[0]) + (g[1] != r[1])
    order = np.arange(nsnv)
    if nsnv > a.max_sites:
        order = np.sort(np.random.default_rng(1).choice(nsnv, a.max_sites, replace=False))
    # cells ordered by similarity: sort by the pattern of the most frequent variants
    freq = (code[:, order] % 3 > 0).mean(axis=0)
    top = order[np.argsort(-freq)[:min(30, len(order))]]
    cell_order = np.lexsort(((code[:, top] % 3 > 0).T[::-1]))
    order = order[np.argsort(-freq, kind="stable")]       # SNVs from the most to the least shared
    cmap = matplotlib.colors.ListedColormap(["#f0f0f0", "#f6b26b", "#b22222", "#7f7f7f"])
    fig, ax = plt.subplots(figsize=(10, max(3, min(12, ncell * 0.04 + 2))))
    ax.imshow(code[np.ix_(cell_order, order)], aspect="auto", cmap=cmap, vmin=0, vmax=3, interpolation="nearest")
    ax.set(xlabel="SNVs, from the most to the least shared (%d of %d)" % (len(order), nsnv), ylabel="cells (%d)" % ncell, title="Genotypes relative to the healthy cell")
    handles = [matplotlib.patches.Patch(color=c, label=l) for c, l in zip(["#f0f0f0", "#f6b26b", "#b22222", "#7f7f7f"], ["as healthy", "one allele differs", "both alleles differ", "missing"])]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=4, frameon=False)
    fig.tight_layout(); fig.savefig(a.out, dpi=150)
    print("%d cells, %d SNVs: %.1f%% as healthy, %.1f%% heterozygous-like, %.1f%% homozygous-like, %.1f%% missing" %
          (ncell, nsnv, *[100 * (code == k).mean() for k in range(4)]))
    print("wrote", a.out)


def find_file(results, names):
    """first existing file among the candidate names (relative to the Results folder)"""
    import os
    for n in names:
        p = os.path.join(results, n)
        if os.path.isfile(p):
            return p
    return None


def cmd_all(a):
    import base64, contextlib, io, os, subprocess, html
    r, rep = a.file, "%04d" % a.replicate
    if not os.path.isdir(r):
        sys.exit("%s is not a folder: give the Results folder of a run" % r)
    os.makedirs(a.outdir, exist_ok=True)
    vcf = find_file(r, ["vcf_dir/vcf." + rep, "vcf_dir/vcf." + rep + ".gz", "vcf", "vcf.gz"])
    gen = find_file(r, ["snv_genotypes_dir/snv_gen." + rep, "snv_gen"])
    bulk = find_file(r, ["bulk_dir/bulk." + rep, "bulk"])
    sections = []   # (title, png path or None, text)

    def run(title, fn, **kw):
        out = os.path.join(a.outdir, kw.pop("png"))
        ns = argparse.Namespace(file=kw.pop("file"), out=out, tsv=None, bins=a.bins, max_sites=a.max_sites)
        buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf):
                fn(ns)
            sections.append((title, out, buf.getvalue()))
        except SystemExit as e:
            sections.append((title, None, str(e)))

    if vcf:
        run("Read depth", cmd_depth, file=vcf, png="depth.png")
        run("Genotypes without reads", cmd_missing, file=vcf, png="missing.png")
    else:
        sections.append(("Read depth, genotypes without reads", None, "no VCF with reads was found (needs -C and a VCF, not -Z)"))
    if gen:
        run("Genotype heatmap", cmd_heatmap, file=gen, png="heatmap.png")
    else:
        sections.append(("Genotype heatmap", None, "no SNV genotype file was found (needs option -1)"))
    vaf_src = [("--bulk", bulk)] if bulk else [("--vcf", vcf)] if vcf else [(None, gen)] if gen else []
    if vaf_src:
        png = os.path.join(a.outdir, "vaf.png")
        flag, path = vaf_src[0]
        cmd = [sys.executable, os.path.join(os.path.dirname(os.path.abspath(__file__)), "bulk_vaf.py")] + ([flag, path] if flag else [path]) + ["--out", png]
        res = subprocess.run(cmd, capture_output=True, text=True)
        sections.append(("Bulk-like VAF distribution", png if os.path.isfile(png) else None, (res.stdout + res.stderr).strip()))
    meta = ""
    rj = find_file(r, ["run.json"])
    if rj:
        import json
        j = json.load(open(rj))
        meta = "%s %s, %s cells x %s sites, %s replicates, seed %s, %s s" % (j.get("program", "CellCoal"), j.get("version", ""), j.get("cells"), j.get("sites"), j.get("replicates"), j.get("seed"), j.get("wall_seconds"))
    parts = ["<!doctype html><meta charset=utf-8><title>CellCoal report</title>",
             "<style>body{font:15px/1.5 system-ui,sans-serif;max-width:1100px;margin:2em auto;padding:0 1em;color:#222}img{max-width:100%}pre{background:#f5f5f5;padding:.6em;overflow:auto}h2{margin-top:1.6em}</style>",
             "<h1>CellCoal report</h1><p>%s</p><p>Folder: <code>%s</code>, replicate %d</p>" % (html.escape(meta), html.escape(os.path.abspath(r)), a.replicate)]
    for title, png, text in sections:
        parts.append("<h2>%s</h2>" % html.escape(title))
        if png and os.path.isfile(png):
            parts.append('<img alt="%s" src="data:image/png;base64,%s">' % (html.escape(title), base64.b64encode(open(png, "rb").read()).decode()))
        if text:
            parts.append("<pre>%s</pre>" % html.escape(text))
    out = os.path.join(a.outdir, "report.html")
    open(out, "w").write("\n".join(parts))
    for title, png, text in sections:
        print("%-38s %s" % (title, "ok" if png else "skipped: " + text.splitlines()[0][:70] if text else "skipped"))
    print("wrote", out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name, fn, out in (("depth", cmd_depth, "depth.png"), ("missing", cmd_missing, "missing.png"), ("heatmap", cmd_heatmap, "heatmap.png")):
        p = sub.add_parser(name)
        p.add_argument("file")
        p.add_argument("--out", default=out)
        p.add_argument("--tsv", default=None, help="write the per-cell numbers to this file")
        p.add_argument("--bins", type=int, default=40)
        p.add_argument("--max-sites", type=int, default=2000, help="heatmap: sample at most this many SNVs")
        p.set_defaults(fn=fn)
    p = sub.add_parser("all")
    p.add_argument("file", help="the Results folder of a run")
    p.add_argument("--outdir", default="plots")
    p.add_argument("--replicate", type=int, default=1, help="replicate whose files are used when the run has one folder per replicate (-v)")
    p.add_argument("--bins", type=int, default=40)
    p.add_argument("--max-sites", type=int, default=2000)
    p.set_defaults(fn=cmd_all)
    a = ap.parse_args()
    a.fn(a)
