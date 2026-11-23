#!/usr/bin/env python3
"""Makes the figures of the manual (docs/manual/media/fig_*.png) with matplotlib. The panels with data come from real runs of CellCoal.

Usage: python3 tools/make_figures.py [--exe bin/cellcoal-2.7.0]   (from the root of the repository)
"""
import argparse, os, re, subprocess, tempfile, collections
import matplotlib, matplotlib.ticker
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, FancyArrowPatch, Rectangle

ap = argparse.ArgumentParser()
ap.add_argument("--exe", default="bin/cellcoal-2.7.0")
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
        ax.text(x, n - 0.2 + (0.0 if k == 0 else 1.1), labels[k], color=COL[k + 1], ha="center", va="bottom", fontsize=7.5)
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
    dm = run(*base, "-#4", "--demes", "n=7:6:6,size=1:0.5:0.3,time=1.5:0.6,source=1:2,founders=1:1", out="dA")
    dp = run(*base, "-#4", "--demes", "n=7:6:6,size=1:0.5:0.3,time=1.5:0.6,source=1:2,founders=1:4", out="dB")
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
    tree_panel(axb, dm, N, u, [1.5, 0.6], "b  One founder: every metastasis is a clade", lab)
    tree_panel(axc, dp, N, u, [1.5, 0.6], "c  Four founders of deme 3: polyclonal seeding", lab)
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


fig_flow()
fig_demes()
fig_blocks()
fig_indels()
print("figures in", media)
