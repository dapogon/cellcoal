#!/usr/bin/env python3
"""Makes the figures of the manual (docs/manual/media/fig_*.png) with matplotlib. The panels with data come from real runs of CellCoal.

Usage: python3 tools/make_figures.py [--exe bin/cellcoal-2.7.4]   (from the root of the repository)
"""
import argparse, os, re, subprocess, tempfile, collections
import matplotlib, matplotlib.ticker
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, FancyArrowPatch, Rectangle

ap = argparse.ArgumentParser()
ap.add_argument("--exe", default="bin/cellcoal-2.7.4")
args = ap.parse_args()
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
exe = os.path.join(root, args.exe)
media = os.path.join(root, "docs", "manual", "media")
TMP = tempfile.mkdtemp()
COL = ["#d95f02", "#1b9e77", "#7570b3", "#e7298a"]
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})


def run(*a, out="X"):
    d = os.path.join(TMP, out)
    subprocess.run([exe] + list(a) + ["-o" + d], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True, cwd=TMP)
    return d


# ---------- Newick ----------
def parse(s):
    pos = [0]

    def node():
        if s[pos[0]] == "(":
            pos[0] += 1
            ch = [node()]
            while s[pos[0]] == ",":
                pos[0] += 1
                ch.append(node())
            pos[0] += 1
            name = ""
        else:
            m = re.match(r"[^:,()]+", s[pos[0]:])
            name = m.group(0)
            pos[0] += len(name)
            ch = []
        L = 0.0
        if pos[0] < len(s) and s[pos[0]] == ":":
            m = re.match(r":([-0-9.eE+]+)", s[pos[0]:])
            L = float(m.group(1))
            pos[0] += len(m.group(0))
        return {"name": name, "ch": ch, "len": L}
    return node()


def layout(n, depth=0.0):
    n["x"] = depth + n["len"]
    for c in n["ch"]:
        layout(c, n["x"])
    return n


def tips(n):
    return [n] if not n["ch"] else [t for c in n["ch"] for t in tips(c)]


def draw_tree(ax, tree, deme, lw=1.0):
    ts = tips(tree)
    for i, t in enumerate(ts):
        t["y"] = i
    def assign(n):
        if n["ch"]:
            for c in n["ch"]:
                assign(c)
            n["y"] = sum(c["y"] for c in n["ch"]) / len(n["ch"])
    assign(tree)

    def col(n):
        ds = {deme.get(t["name"], 0) for t in tips(n)}
        return COL[(ds.pop() - 1) % 4] if len(ds) == 1 else "#888888"
    def rec(n, px):
        for c in n["ch"]:
            ax.plot([n["x"], n["x"]], [n["ch"][0]["y"], n["ch"][-1]["y"]], color="#888888", lw=lw, solid_capstyle="butt")
            ax.plot([n["x"], c["x"]], [c["y"], c["y"]], color=col(c), lw=lw + 0.4, solid_capstyle="butt")
            rec(c, n["x"])
    rec(tree, 0)
    H = max(t["x"] for t in ts)
    return H, len(ts)


def demes_of(d):
    deme = {}
    for l in open(os.path.join(d, "demes.tsv")).read().split("\n")[1:]:
        if l:
            a, b = l.split("\t")
            deme[a] = int(b)
    return deme


def tree_panel(ax, d, N, u, times, title, labels):
    t = parse(open(os.path.join(d, "trees")).read().split("\n")[0].strip().rstrip(";"))
    t["ch"] = [c for c in t["ch"] if not c["name"].startswith("outg")]	# the outgroup is not drawn
    if len(t["ch"]) == 1:
        t = t["ch"][0]
    t["len"] = 0.0
    t = layout(t)
    H, n = draw_tree(ax, t, demes_of(d))
    for k, T in enumerate(times):
        x = H - T * N * u
        ax.axvline(x, color=COL[k + 1], ls="--", lw=0.9)
        ax.text(x, n - 0.2 + (0.0 if k == 0 else 1.1), labels[k], color=COL[k + 1], ha="center", va="bottom", fontsize=7.5, bbox=dict(facecolor="white", edgecolor="none", pad=0.6))
    ax.set_title(title, fontsize=9, loc="left", pad=22)
    ax.set_yticks([]); ax.set_xlabel("substitutions per site (x 1000) from the tumor root")
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, p: "%g" % (v * 1000)))
    ax.spines["left"].set_visible(False)
    ax.set_ylim(-0.7, n + 1.6)
    return H


# ================= Figure: structured coalescent =================
def fig_demes():
    N, u = 10000, 1e-7
    base = ["-l200", "-u1e-7", "-e10000", "-6", "-7"]
    spec = "n=6:5:6,size=1:0.5:2,time=1.5:0.5,source=1:2,founders=1:%d"

    def monophyletic(d, k):		# are the cells of the deme k a clade of the genealogy?
        deme = demes_of(d)
        tr = parse(open(os.path.join(d, "trees")).read().split("\n")[0].strip().rstrip(";"))
        want = {c for c, v in deme.items() if v == k}

        def rec(n):
            names = {x["name"] for x in tips(n)}
            if names == want:
                return True
            return any(rec(c) for c in n["ch"])
        return rec(tr)
    dm = dp = None
    for sd in range(1, 60):
        if dm is None:
            d = run(*base, "-#%d" % sd, "--demes", spec % 1, out="dA%d" % sd)
            if monophyletic(d, 3) and not monophyletic(d, 2):
                dm = d
        if dp is None:
            d = run(*base, "-#%d" % sd, "--demes", spec % 3, out="dB%d" % sd)
            if not monophyletic(d, 3):
                dp = d
        if dm and dp:
            break
    fig = plt.figure(figsize=(11, 4.6))
    ax0 = fig.add_axes([0.01, 0.06, 0.27, 0.84]); ax0.axis("off"); ax0.set_xlim(0, 10); ax0.set_ylim(0, 10)
    ax0.set_title("a  The model", fontsize=9, loc="left", fontweight="bold")
    ax0.add_patch(Ellipse((3.6, 7.6), 6.4, 3.6, fc=COL[0], alpha=0.25, ec=COL[0])); ax0.text(3.6, 7.6, "deme 1\nprimary tumor", ha="center", va="center")
    ax0.add_patch(Ellipse((2.3, 3.4), 3.8, 2.4, fc=COL[1], alpha=0.25, ec=COL[1])); ax0.text(2.3, 3.4, "deme 2\nmetastasis", ha="center", va="center")
    ax0.add_patch(Ellipse((7.8, 2.2), 3.4, 2.2, fc=COL[2], alpha=0.25, ec=COL[2])); ax0.text(7.8, 2.2, "deme 3\nmetastasis", ha="center", va="center", fontsize=8.5)
    ax0.add_patch(FancyArrowPatch((2.9, 5.8), (2.4, 4.6), arrowstyle="-|>", mutation_scale=14, color="k"))
    ax0.text(4.4, 5.2, "colonization at T2\nby F2 founders", fontsize=8, va="center")
    ax0.add_patch(FancyArrowPatch((3.5, 3.0), (6.1, 2.4), arrowstyle="-|>", mutation_scale=14, color="k"))
    ax0.text(3.2, 0.5, "colonization from deme 2\nat T3 (later than T2)", fontsize=8, va="center")
    ax0.text(5.2, 4.0, "no exchange of cells\nafter the colonization", fontsize=8, ha="left", style="italic", color="#555555")
    axb = fig.add_axes([0.35, 0.16, 0.30, 0.66]); axc = fig.add_axes([0.70, 0.16, 0.28, 0.66])
    lab = ["deme 2 founded", "deme 3 founded"]
    tree_panel(axb, dm, N, u, [1.5, 0.5], "b  One founder: deme 3 is a clade (nested in deme 2)", lab)
    tree_panel(axc, dp, N, u, [1.5, 0.5], "c  Three founders of deme 3: polyclonal, not a clade", lab)
    for a in (axb, axc):
        a.title.set_fontweight("normal")
    fig.savefig(os.path.join(media, "fig_demes.png"), dpi=170)
    plt.close(fig)


# ================= Figure: blocks and regional rates =================
def fig_blocks():
    fig = plt.figure(figsize=(11, 5.2))
    ax = fig.add_axes([0.02, 0.50, 0.96, 0.44]); ax.axis("off"); ax.set_xlim(0, 100); ax.set_ylim(0, 10)
    ax.set_title("a  A chromosome in blocks: one genealogy, blocks of -l sites that are simulated one after the other", fontsize=9, loc="left", fontweight="bold")
    # one tree on top
    ax.plot([50, 50], [8.6, 9.4], color="#555555"); ax.plot([30, 70], [8.6, 8.6], color="#555555")
    for x in (30, 70):
        ax.plot([x, x], [7.9, 8.6], color="#555555")
    for x0, x1 in ((22, 38), (62, 78)):
        ax.plot([x0, x1], [7.9, 7.9], color="#555555")
        for x in (x0, x1):
            ax.plot([x, x], [7.3, 7.9], color="#555555")
    ax.text(51.5, 9.3, "tree (made in the first block)", fontsize=8, color="#555555")
    xs = [2, 26, 50, 74]
    for k, x in enumerate(xs):
        ax.add_patch(Rectangle((x, 3.6), 23.5, 1.6, fc="#dddddd", ec="k"))
        ax.text(x + 11.7, 5.55, "block %d  (-l sites)" % (k + 1), ha="center", va="center", fontsize=8.5)
        ax.annotate("", xy=(x + 11.7, 6.0), xytext=(50, 7.0), arrowprops=dict(arrowstyle="-", color="#bbbbbb", lw=0.7))
        for m in range(5):
            ax.plot([x + 2 + m * 4.6 + (k * 1.7) % 3] * 2, [3.6, 5.2], color="#c0392b", lw=1.2)
    ax.text(2, 2.7, "SNVs (red) of every block: independent sites given the tree", fontsize=8, color="#c0392b")
    ax.annotate("", xy=(97.5, 1.2), xytext=(2, 1.2), arrowprops=dict(arrowstyle="->"))
    ax.text(50, 0.3, "positions of the chromosome (simulated with --blocks, or real with --fasta and --region)", ha="center", fontsize=8.5)
    # real data: region rates
    exp_mult = [(0, 10000, 4.0), (10000, 25000, 0.25)]
    bed = os.path.join(TMP, "hot.bed")
    open(bed, "w").write("chr1\t0\t10000\t4\nchr1\t10000\t25000\t0.25\n")
    counts = collections.defaultdict(list)
    for sd in range(1, 9):
        d = run("-s10", "-#%d" % sd, "-u4e-7", "-e10000", "-l5000", "--blocks", "8", "--region-rates", bed, "-1", out="rr%d" % sd)
        k = 0
        for l in open(os.path.join(d, "snv_gen")):
            p = l.split()
            if len(p) == 2 and p[0] == "11":
                counts[k].append(int(p[1])); k += 1
    ax1 = fig.add_axes([0.07, 0.09, 0.84, 0.30])
    mean = [sum(counts[k]) / len(counts[k]) for k in sorted(counts)]
    ax1.bar([2500 + 5000 * k for k in range(len(mean))], mean, width=4200, color="#c0392b", alpha=0.85)
    ax1.set_ylabel("SNVs per block (mean of 8 runs)")
    ax1.set_xlabel("position of the block")
    ax1.set_xlim(0, 40000)
    ax2 = ax1.twinx(); ax2.spines["right"].set_visible(True)
    xs2 = [0, 10000, 10000, 25000, 25000, 40000]; ys2 = [4, 4, 0.25, 0.25, 1, 1]
    ax2.plot(xs2, ys2, color="k", lw=1.6); ax2.set_ylabel("multiplier of --region-rates"); ax2.set_ylim(0, 4.6)
    ax1.set_title("b  --region-rates (here 4 and 0.25): the blocks follow the multipliers of their regions; the total number of SNVs does not change", fontsize=9, loc="left", fontweight="bold")
    fig.savefig(os.path.join(media, "fig_blocks.png"), dpi=170)
    plt.close(fig)


# ================= Figure: indels =================
def seqrow(ax, x, y, s, cols=None, dx=0.62, fs=9):
    for i, ch in enumerate(s):
        c = cols[i] if cols else "k"
        ax.text(x + i * dx, y, ch, family="monospace", fontsize=fs, ha="center", va="center", color=c)


def fig_indels():
    fig = plt.figure(figsize=(11, 7.6))
    ax = fig.add_axes([0.01, 0.50, 0.98, 0.46]); ax.axis("off"); ax.set_xlim(0, 100); ax.set_ylim(0, 10)
    ax.set_title("a  An indel is inherited by the cells below its branch and seen in the reads as I or D in the CIGAR", fontsize=9, loc="left", fontweight="bold")
    ref = "ACGTCATTGCAGGTTACCGATCGAATGCCATGGA"
    d0, dl, ins_after, ins = 11, 3, 24, "TT"
    rcols, hcols, rcl, hcl, kind = [], [], [], [], []		# columns of the alignment of the reference and the haplotype
    i = 0
    while i < len(ref):
        if i == d0:
            for k in range(dl):
                rcols.append(ref[i + k]); hcols.append("-"); rcl.append("k"); hcl.append("#c0392b"); kind.append("D")
            i += dl
            continue
        rcols.append(ref[i]); hcols.append(ref[i]); rcl.append("k"); hcl.append("k"); kind.append("M")
        if i == ins_after:
            for ch in ins:
                rcols.append("-"); hcols.append(ch); rcl.append("#2471a3"); hcl.append("#2471a3"); kind.append("I")
        i += 1
    dx = 0.62; X0 = 14
    seqrow(ax, X0, 8.6, "".join(rcols), rcl); ax.text(1, 8.6, "reference", fontsize=8.5, va="center")
    seqrow(ax, X0, 6.6, "".join(hcols), hcl); ax.text(1, 6.6, "haplotype of\nthe cell", fontsize=8.5, va="center")
    ax.text(X0 + (d0 + 1) * dx, 7.5, "deletion (3 bases)", color="#c0392b", fontsize=8, ha="center")
    ax.text(X0 + (ins_after + 2.5) * dx, 7.5, "insertion (2 bases)", color="#2471a3", fontsize=8, ha="center")

    def cigar(c0, c1):
        ops = []
        for k in kind[c0:c1]:
            if ops and ops[-1][1] == k:
                ops[-1][0] += 1
            else:
                ops.append([1, k])
        return "".join("%d%s" % (n, o) for n, o in ops)

    def read(y, c0, L, name):
        c1, q = c0, 0
        while q < L and c1 < len(kind):		# the read has L bases: M and I columns count, D columns do not
            q += kind[c1] != "D"
            c1 += 1
        ax.add_patch(Rectangle((X0 + c0 * dx - dx / 2, y - 0.3), (c1 - c0) * dx, 0.6, fc="#555555", alpha=0.35, ec="#555555"))
        ax.text(X0 + c1 * dx + 0.4, y, "%s  %d bases, CIGAR %s" % (name, q, cigar(c0, c1)), fontsize=8, va="center", family="monospace")
    read(4.7, 0, 16, "read 1:")
    read(3.6, 8, 16, "read 2:")
    read(2.5, 18, 16, "read 3:")
    ax.text(1, 1.2, "A read has L bases of the haplotype: with a deletion it covers more reference sites, with an insertion fewer. The CIGAR counts M (bases that match the reference sites), D and I.\n"
                    "A base inserted at the end of a read is soft-clipped (S). The truth (position, bases, genotype of every cell) is in reads_<n>/indels_truth.vcf.", fontsize=8.5, va="center")
    # b: left alignment
    axb = fig.add_axes([0.03, 0.06, 0.30, 0.36]); axb.axis("off"); axb.set_xlim(0, 10); axb.set_ylim(0, 10)
    axb.set_title("b  Left-aligned, as a caller reports them", fontsize=9, loc="left", fontweight="bold")
    s = "GACTTTTGAC"
    seqrow(axb, 1.5, 8.4, s, dx=0.8); axb.text(0, 9.3, "reference", fontsize=8)
    for k, (lab, delpos) in enumerate((("deleting any T gives the same sequence", None),)):
        pass
    axb.text(0, 6.7, "one T less in the run TTTT:", fontsize=8)
    for k, p in enumerate((3, 4, 5, 6)):
        t = list(s); t[p] = "-"
        seqrow(axb, 1.5, 5.6 - k * 0.95, "".join(t), ["#c0392b" if i == p else "k" for i in range(len(s))], dx=0.8)
    axb.annotate("CellCoal uses this one\n(leftmost position)", xy=(4.3, 5.75), xytext=(5.6, 6.9), fontsize=8, color="#1e8449", arrowprops=dict(arrowstyle="->", color="#1e8449"))
    axb.text(0, 0.2, "Equivalent placements are reported at the left end.", fontsize=8)
    # c: homopolymer error
    axc = fig.add_axes([0.36, 0.06, 0.30, 0.36]); axc.axis("off"); axc.set_xlim(0, 10); axc.set_ylim(0, 10)
    axc.set_title("c  --homopolymer-error", fontsize=9, loc="left", fontweight="bold")
    s0 = "CAGAAAAAGTC"
    seqrow(axc, 1.0, 8.2, s0, dx=0.75); axc.text(0, 9.2, "template (run of 5 A)", fontsize=8)
    seqrow(axc, 1.0, 6.0, "CAGAAAA-GTC", ["k"] * 7 + ["#c0392b"] + ["k"] * 3, dx=0.75); axc.text(0, 6.9, "read with a deletion: a D in the CIGAR", fontsize=8)
    seqrow(axc, 1.0, 3.8, "CAGAAAAAAGTC", ["k"] * 8 + ["#2471a3"] + ["k"] * 3, dx=0.75); axc.text(0, 4.7, "read with an insertion: an I in the CIGAR", fontsize=8)
    axc.text(0, 1.9, "A run of 3 or more equal bases changes by\none base with probability P x (length - 2).\nIndependent in each read; not in the truth.", fontsize=8, va="center")
    # d: length distribution (real)
    axd = fig.add_axes([0.73, 0.10, 0.25, 0.30])
    for k, (ln, mx, c) in enumerate(((2, 10, "#1b9e77"), (4, 10, "#d95f02"))):
        lens = []
        for sd in range(1, 6):
            d = run("-s6", "-#%d" % sd, "-l30000", "-u1e-8", "-e1000", "-b1", "-C10", "-E0", "-D0", "--bam", "--indel-rate", "300", "--indel-len", str(ln), "--indel-max", str(mx), out="il%d_%d" % (ln, sd))
            for dp, _, fs in os.walk(d):
                for f in fs:
                    if f == "indels_truth.vcf":
                        for l in open(os.path.join(dp, f)):
                            if not l.startswith("#"):
                                m = re.search(r"LEN=(\d+)", l)
                                lens.append(int(m.group(1)))
        cnt = collections.Counter(lens); tot = sum(cnt.values())
        axd.bar([i + (k - 0.5) * 0.4 for i in range(1, 11)], [cnt.get(i, 0) / tot for i in range(1, 11)], width=0.4, color=c, label="--indel-len %d" % ln)
    axd.set_xticks(range(1, 11)); axd.set_xlabel("length of the indel (bases)"); axd.set_ylabel("fraction"); axd.legend(frameon=False, fontsize=8)
    axd.set_title("d  Lengths: geometric, up to --indel-max (10)", fontsize=9, loc="left", fontweight="bold")
    fig.savefig(os.path.join(media, "fig_indels.png"), dpi=170)
    plt.close(fig)


# ================= Figure: main flow =================
def fig_flow():
    from matplotlib.patches import FancyBboxPatch
    fig = plt.figure(figsize=(11, 8.4))
    ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off"); ax.set_xlim(0, 100); ax.set_ylim(0, 100)
    G, B, O, P = "#e5e5e5", "#d6e6f5", "#fbe3c8", "#e1d9ef"

    def box(x, y, w, h, text, fc=G, fs=8.6, bold=False):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.25,rounding_size=1.0", fc=fc, ec="#555555", lw=0.9))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, fontweight="bold" if bold else "normal", linespacing=1.25)

    def arrow(x0, y0, x1, y1):
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=12, color="#666666", lw=1.2))

    # inputs
    box(30, 93, 40, 5.2, "command line  |  parameter file  |  configuration file (--config)\n(options in groups: --tree, --demes, --reads ...)", B, 8.6)
    # genealogy
    box(2, 76, 30, 11, "Genealogy\ncoalescent (constant, growth, demographic periods)\nor a user tree (Newick)\ntumor regions: --demes (colonizations)", G)
    box(2, 61, 30, 10, "outgroup (healthy cell)\nvariation of the rate\namong lineages and along the tree", G)
    arrow(50, 92.5, 17, 87.5); arrow(17, 76, 17, 71.2)
    # genome
    box(38, 76, 30, 11, "Genome\nsimulated (-l, --blocks) or real reference (--fasta)\ngermline variation (--germline-vcf, --germline-pop)", G)
    arrow(50, 92.5, 53, 87.5)
    box(38, 61, 30, 10, "ancestral genome\n(maternal and paternal haplotypes)", G)
    arrow(53, 76, 53, 71.2)
    # somatic
    box(2, 40, 66, 15.5, "Somatic variation on the branches of the tree (in random order)\nSNVs: infinite and finite sites models, signatures   |   single-nucleotide deletions   |   cnLOH\nlosses and gains of segments (copy number)   |   small indels (reads mode: --indel-rate)", O, 8.6)
    arrow(17, 61, 17, 55.8); arrow(53, 61, 53, 55.8)
    box(73, 44, 25, 7.5, "true genotypes, haplotypes,\ncopy number, genealogy and\ntruth files (truth.vcf ...)", P, 8.4)
    arrow(68.3, 47.7, 72.7, 47.7)
    # technical layer
    ax.text(32.5, 37.6, "three modes of observation", fontsize=9, fontweight="bold", ha="center")
    box(2, 17.5, 28, 16.5, "1  Genotypes\nallelic dropout (ADO)\ngenotype errors\n(no read counts)", G)
    box(34, 17.5, 30, 16.5, "2  Read counts (-C)\ncoverage, allelic imbalance, ADO,\namplification gain along the genome,\namplification and sequencing errors", G)
    box(68, 17.5, 30, 16.5, "3  Reads (--bam, --fastq)\npaired-end reads of the amplified cells:\nPCR duplicates, qualities, strand bias,\ndamage, doublets, indels, homopolymer errors", G)
    for x in (16, 49, 83):
        arrow(min(max(x, 10), 60), 40, x, 34.3)
    # outputs
    box(2, 2, 28, 10.5, "observed genotypes\n(VCF, genotype files)", P, 8.6)
    box(34, 2, 30, 10.5, "genotype likelihoods (several models)\nML genotypes, VCF with likelihoods,\nbulk read counts, coverage", P, 8.4)
    box(68, 2, 30, 10.5, "BAM + BAI and FASTQ of every cell\ntruth.vcf, indels_truth.vcf, doublets", P, 8.4)
    arrow(16, 17.5, 16, 12.8); arrow(49, 17.5, 49, 12.8); arrow(83, 17.5, 83, 12.8)
    fig.savefig(os.path.join(media, "fig_flow.png"), dpi=170)
    plt.close(fig)


def bam_depth(path, L, minq=0):
    import gzip, struct
    data = gzip.open(path, "rb").read()
    o = 4; lt = struct.unpack_from("<i", data, o)[0]; o += 4 + lt
    nref = struct.unpack_from("<i", data, o)[0]; o += 4
    for _ in range(nref):
        ln = struct.unpack_from("<i", data, o)[0]; o += 4 + ln + 4
    dep = [0] * (L + 2)
    while o < len(data):
        bs = struct.unpack_from("<i", data, o)[0]
        ref, pos, lrn, mapq, bn, ncig, flag, lseq = struct.unpack_from("<iiBBHHHi", data, o + 4)
        co = o + 4 + 32 + lrn; p = pos
        for k in range(ncig):
            v = struct.unpack_from("<I", data, co + 4 * k)[0]; n = v >> 4; op = v & 15
            if op in (0, 7, 8):
                if mapq >= minq:
                    for q in range(p, min(p + n, L)): dep[q] += 1
                p += n
            elif op in (2, 3): p += n
        o += 4 + bs
    return dep


def fig_copynumber():
    L, C = 60000, 20
    d = run("-n1", "-s12", "-#3", "-l%d" % L, "-u1e-8", "-e1000", "-C%d" % C, "-E0.001", "--seg-loss-rate", "3e-4", "--seg-gain-rate", "3e-4", "--seg-len", "5000", "--bam", out="cn")
    seg = collections.defaultdict(list)
    for l in open(os.path.join(d, "cnv_segments")):
        p = l.rstrip("\n").split("\t")
        if len(p) == 6 and p[1] != "start": seg[p[0]].append((int(p[1]), int(p[2]), int(p[5])))
    dep = {c: bam_depth(os.path.join(d, "reads_1", c + ".bam"), L) for c in seg}
    fig = plt.figure(figsize=(11, 6.2))
    gs = fig.add_gridspec(3, 2, width_ratios=[2.2, 1], hspace=0.55, wspace=0.28, left=0.07, right=0.98, top=0.92, bottom=0.1)
    # three cells with the most altered genome
    cells = sorted(seg, key=lambda c: -sum(1 for s in seg[c] if s[2] != 2))[:3]
    B = 400
    for i, c in enumerate(cells):
        ax = fig.add_subplot(gs[i, 0])
        xs = [b * B + B / 2 for b in range(L // B)]
        ys = [sum(dep[c][b * B:(b + 1) * B]) / B for b in range(L // B)]
        ax.plot(xs, ys, color="#555555", lw=0.9, label="mean depth of the reads (400 bp)")
        tx, ty = [], []
        for a, b, cn in seg[c]:
            tx += [a, b]; ty += [cn * C / 2.0] * 2
        ax.plot(tx, ty, color="#d95f02", lw=1.6, label="true copies x C/2")
        ax.set_ylim(0, 2.4 * C); ax.set_xlim(0, L); ax.set_ylabel("depth")
        ax.text(0.005, 0.86, c, transform=ax.transAxes, fontsize=8, fontweight="bold")
        if i == 0:
            ax.set_title("a  Depth along the chromosome follows the true number of copies of every segment", fontsize=9, loc="left", fontweight="bold")
            ax.legend(fontsize=7.5, loc="upper right", frameon=False, ncol=2)
        if i == 2: ax.set_xlabel("position (bp)")
    ax = fig.add_subplot(gs[:, 1])
    pts = collections.defaultdict(list)
    for c, ss in seg.items():
        for a, b, cn in ss:
            if b - a < 1500: continue
            pts[cn].append(sum(dep[c][a + 200:b - 200]) / max(1, b - a - 400))
    for cn in sorted(pts):
        ax.scatter([cn] * len(pts[cn]), pts[cn], s=14, color="#1b9e77", alpha=0.7)
    ax.plot([0, max(pts) + 0.3], [0, (max(pts) + 0.3) * C / 2.0], color="#d95f02", lw=1.4)
    ax.set_xlabel("true number of copies (cnv_segments)"); ax.set_ylabel("mean depth of the segment")
    ax.set_title("b  Every segment", fontsize=9, loc="left", fontweight="bold")
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    fig.savefig(os.path.join(media, "fig_copynumber.png"), dpi=170)
    plt.close(fig)


def fig_mappability():
    L, C = 60000, 30
    bed = os.path.join(TMP, "map.bed")
    regs = [(8000, 14000, 0.25), (24000, 30000, 0.6), (42000, 48000, 0.9)]
    open(bed, "w").write("".join("1\t%d\t%d\t%s\n" % r for r in regs))
    d = run("-n1", "-s4", "-#5", "-l%d" % L, "-u1e-8", "-e1000", "-C%d" % C, "-E0", "--bam", "--mappability", bed, out="mp")
    p = os.path.join(d, "reads_1", "cell0001.bam")
    alld, good = bam_depth(p, L), bam_depth(p, L, 20)
    B = 300
    xs = [b * B + B / 2 for b in range(L // B)]
    sm = lambda v: [sum(v[b * B:(b + 1) * B]) / B for b in range(L // B)]
    fig = plt.figure(figsize=(11, 5.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 2.6], hspace=0.25, left=0.07, right=0.98, top=0.92, bottom=0.11)
    ax = fig.add_subplot(gs[0])
    sx, sy = [0], [1]
    for a, b, m in regs:
        sx += [a, a, b, b]; sy += [1, m, m, 1]
    sx.append(L); sy.append(1)
    ax.fill_between(sx, sy, step=None, color="#7570b3", alpha=0.35); ax.plot(sx, sy, color="#7570b3", lw=1.4)
    ax.set_ylim(0, 1.15); ax.set_xlim(0, L); ax.set_ylabel("mappability"); ax.set_xticklabels([])
    for a, b, m in regs:
        mq = 60 if m >= 1 else (0 if m <= 0.5 else int(60 * (m - 0.5) / 0.5 + 0.5))
        ax.text((a + b) / 2, 0.55 * m, "score %s\nMAPQ %d" % (m, mq), ha="center", va="center", fontsize=8)
    ax.set_title("a  Mappability of the genome (--mappability): the score of each region gives the MAPQ of the reads that cover it", fontsize=9, loc="left", fontweight="bold")
    ax = fig.add_subplot(gs[1])
    ax.plot(xs, sm(alld), color="#999999", lw=1.2, label="all the reads")
    ax.plot(xs, sm(good), color="#d95f02", lw=1.6, label="reads with MAPQ of 20 or more (a usual filter)")
    ax.set_xlim(0, L); ax.set_ylim(0, 1.5 * C); ax.set_xlabel("position (bp)"); ax.set_ylabel("depth")
    ax.legend(fontsize=8, loc="upper right", frameon=False, ncol=2)
    ax.set_title("b  The depth does not change, but a caller that filters by MAPQ loses the reads of the repeated regions", fontsize=9, loc="left", fontweight="bold")
    fig.savefig(os.path.join(media, "fig_mappability.png"), dpi=170)
    plt.close(fig)


def fig_targets():
    L, C = 60000, 30
    bed = os.path.join(TMP, "targets.bed")
    tg = [(8000, 10000), (20000, 21500), (31000, 34000), (44000, 45000), (52000, 55000)]
    open(bed, "w").write("".join("chrT\t%d\t%d\n" % t for t in tg))
    import random
    r = random.Random(5); sq = "".join(r.choice("ACGT") for _ in range(L))
    fa = os.path.join(TMP, "t.fa")
    open(fa, "w").write(">chrT\n" + "\n".join(sq[i:i + 70] for i in range(0, L, 70)) + "\n")
    args = ["-n1", "-s4", "-#5", "-l%d" % L, "--fasta", fa, "--region", "chrT:1-%d" % L, "-u1e-8", "-e1000", "-C%d" % C, "-E0", "--bam"]
    dw = run(*args, out="wgs")
    de = run(*args, "--targets", bed, out="wes")
    B = 250
    xs = [b * B + B / 2 for b in range(L // B)]
    sm = lambda v: [sum(v[b * B:(b + 1) * B]) / B for b in range(L // B)]
    w = bam_depth(os.path.join(dw, "reads_1", "cell0001.bam"), L)
    e = bam_depth(os.path.join(de, "reads_1", "cell0001.bam"), L)
    fig, axs = plt.subplots(2, 1, figsize=(11, 5.4), sharex=True)
    fig.subplots_adjust(left=0.07, right=0.98, top=0.92, bottom=0.11, hspace=0.3)
    for ax, v, col, t in ((axs[0], w, "#999999", "a  Whole genome (scWGS): the depth is even along the chromosome"),
                          (axs[1], e, "#1b9e77", "b  Targeted sequencing (scWES, --targets): the depth is at the targets and at their edges, and 5% elsewhere (--off-target)")):
        ax.plot(xs, sm(v), color=col, lw=1.4)
        ax.set_xlim(0, L); ax.set_ylim(0, 1.5 * C); ax.set_ylabel("depth")
        ax.set_title(t, fontsize=9, loc="left", fontweight="bold")
    for a, b in tg:
        axs[1].axvspan(a, b, color="#1b9e77", alpha=0.15, lw=0)
    axs[1].set_xlabel("position (bp)")
    fig.savefig(os.path.join(media, "fig_targets.png"), dpi=170)
    plt.close(fig)


fig_flow()
fig_demes()
fig_blocks()
fig_indels()
fig_copynumber()
fig_mappability()
fig_targets()
print("figures in", media)
