#!/usr/bin/env python3
"""chrom.py: simulate a long chromosome with CellCoal by blocks (prototype).

One coalescent tree is generated once and every block of sites is simulated on that same tree (no recombination, so the
sites are independent given the tree).  For each block only the variable sites go to the VCF; the rest of the block is
summarized by the number of sites for each genotype of the healthy root (inv_AA ... inv_TT in the summary table).

Usage (random sequence):
  python3 tools/chrom.py --exe bin/cellcoal-2.3.0 --blocks 50 --block-size 200000 --out chr1 [--jobs 4] [--seed 7]
                         [--keep] [--gzip] -- <CellCoal options: -s -u -e -g -C -D -E --min-alt-reads ...>

Usage (real reference sequence):
  python3 tools/chrom.py --exe bin/cellcoal-2.3.0 --fasta hg38.fa --region chr1:10000001-20000000 --out chr1 -- <options>

  --fasta FILE      reference in standard FASTA (several lines per sequence, upper or lower case; an index FILE.fai is
                    used when it exists).  Sites that are not A, C, G or T (N...) are left out of the simulation: they
                    are not in the VCF and not in the invariant counts (they are counted in chrom_info.txt).
  --region R        chrom:start-end, 1-based and inclusive (default: the whole first sequence)
  --germline-vcf F  germline SNPs given by the user (VCF, optionally .gz; SNVs only; GT 0|1, 1|0, 1|1, 0/1, 1/1).
                    Unphased heterozygous genotypes are phased at random.  It cannot be used together with -c (germline
                    SNPs simulated by CellCoal on the reference).

  --ado-block L (after "--") makes the allelic dropout (-D) come in contiguous segments of mean length L sites; the segments are
  continuous from block to block (the origin of every block is passed with --ado-origin, the same --ado-seed to all).

Do not give -n, -l, -#, -o, -T, -U or --summary after "--": the script sets them.  It also sets --no-outgroup and,
unless you give them, --min-alt-reads 2 --min-alt-vaf 0.3 (a site is written only if a cell has a non-reference call with at least
2 alternative reads and 30% of its reads).  A tree can be given with --tree FILE
(it is then not generated).  Outputs, in the --out folder:
  chrom.vcf[.gz]        variable sites of all the blocks (CHROM and POS of the chromosome)
  chrom_summary.tsv     one row per block, plus a TOTAL row
  chrom_tree.nwk        the tree used for every block
  chrom_coverage.tsv    BED-like: mean depth of every cell over the sites that are NOT in the VCF, per window of
                        --coverage-window sites (default 10000; 0 = off); n_sites = non-variable sites of the window
  chrom_coverage_bed/   one BED file per cell (chrom, start, end, mean depth, n_sites) with the same windows
  chrom_info.txt        region, masked sites, germline sites used
"""
import argparse, concurrent.futures, glob, gzip, os, random, re, shutil, subprocess, sys

COUNT_COLS = ["mutations", "deletions", "cnloh", "snv_true", "snv_observed", "singleton_sites"]
WEIGHTED = {"snv_precision": "snv_observed", "snv_recall": "snv_true", "singleton_recall": "singleton_sites"}
COMP = str.maketrans("ACGT", "TGCA")
COVERAGE_WINDOW = 0


def run(cmd, env=None):
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=env)
    if r.returncode != 0:
        sys.exit("command failed: %s\n%s" % (" ".join(cmd), r.stdout.decode(errors="replace")[-2000:]))


def find_one(folder, patterns):
    for p in patterns:
        hits = sorted(glob.glob(os.path.join(folder, p)))
        if hits:
            return hits[0]
    sys.exit("expected output not found in %s (%s)" % (folder, ", ".join(patterns)))


# ---------------------------------------------------------------------------------------- reference sequence
def parse_region(region):
    m = re.match(r"^([^:]+)(?::([\d,]+)-([\d,]+))?$", region or "")
    if not m:
        sys.exit("--region has to look like chr1:10000001-20000000")
    name = m.group(1)
    if m.group(2) is None:
        return name, None, None
    s, e = int(m.group(2).replace(",", "")), int(m.group(3).replace(",", ""))
    if s < 1 or e < s:
        sys.exit("--region: start has to be >= 1 and end >= start")
    return name, s, e


def first_header(fasta):
    with open(fasta) as fh:
        return fh.readline()[1:].strip().split()[0]


def fetch_region(fasta, name, start, end):
    """Sequence (upper case) of name:start-end (1-based, inclusive) from a FASTA file; start/end None = whole sequence."""
    fai = fasta + ".fai"
    if os.path.exists(fai):
        for line in open(fai):
            f = line.split("\t")
            if f[0] == name:
                length, off, lb, lw = int(f[1]), int(f[2]), int(f[3]), int(f[4])
                s = 1 if start is None else start
                e = length if end is None else min(end, length)
                if s > length:
                    sys.exit("--region starts after the end of %s (%d)" % (name, length))
                a, b = s - 1, e
                with open(fasta, "rb") as fh:
                    fh.seek(off + (a // lb) * lw + a % lb)
                    nbytes = (b - a) + ((b - a) // lb + 2) * (lw - lb) + 1
                    raw = fh.read(nbytes).decode(errors="replace")
                seq = "".join(raw.split())[: b - a]
                return seq.upper(), s
        sys.exit("sequence %s is not in %s" % (name, fai))
    # no index: read the file once
    chunks, pos, inside = [], 0, False
    s = 1 if start is None else start
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                if inside:
                    break
                inside = line[1:].split()[0] == name
                continue
            if not inside:
                continue
            line = line.strip()
            a, b = pos, pos + len(line)
            lo, hi = max(a, s - 1), b if end is None else min(b, end)
            if hi > lo:
                chunks.append(line[lo - a: hi - a])
            pos = b
            if end is not None and pos >= end:
                break
    if not chunks:
        sys.exit("sequence %s (or the region) was not found in %s" % (name, fasta))
    return "".join(chunks).upper(), s


def read_germline(path, name, start, end):
    """{pos: (ref, alt, gt)} of the SNVs of the VCF inside the region."""
    opener = gzip.open if path.endswith(".gz") else open
    alt_name = name[3:] if name.startswith("chr") else "chr" + name
    out = {}
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10 or f[0] not in (name, alt_name):
                continue
            pos = int(f[1])
            if pos < start or pos > end or len(f[3]) != 1 or len(f[4]) != 1 or f[4] not in "ACGT" or f[3] not in "ACGT":
                continue
            gt = f[9].split(":")[0]
            if re.match(r"^[01][|/][01]$", gt) and gt not in ("0|0", "0/0"):
                out[pos] = (f[3], f[4], gt)
    return out


def make_block_genome(seq, offset, germ, rng):
    """Two sequences (maternal, paternal) without the sites that are not ACGT; also the genomic position of every kept site."""
    mat, pat, pos = [], [], []
    nref = 0
    for i, ch in enumerate(seq):
        if ch not in "ACGT":
            continue
        p = offset + i
        m = q = ch
        g = germ.get(p) if germ else None
        if g:
            if g[0] != ch:
                nref += 1
            else:
                gt = g[2]
                if gt[1] == "|":
                    a, b = gt[0] == "1", gt[2] == "1"
                else:
                    a, b = (True, True) if gt == "1/1" else ((True, False) if rng.random() < 0.5 else (False, True))
                if a:
                    m = g[1]
                if b:
                    q = g[1]
        mat.append(m)
        pat.append(q)
        pos.append(p)
    return "".join(mat), "".join(pat), pos, nref


# ---------------------------------------------------------------------------------------- blocks
def one_block(exe, opts, tree, b, size, seed, workdir, genome=None, origin=0, ado_seed=0):
    d = os.path.join(workdir, "block_%04d" % b)
    shutil.rmtree(d, ignore_errors=True)
    env = dict(os.environ, OMP_NUM_THREADS="1")
    cmd = [exe, "-n1", "-l%d" % size, "-T%s" % tree, "--summary", "-#%d" % (seed + b + 1), "-o%s" % d] + opts
    if genome:
        cmd.append("-U%s" % genome)
    if COVERAGE_WINDOW > 0:
        cmd += ["--coverage", str(COVERAGE_WINDOW)]
    if "--ado-block" in opts:  # dropout segments that continue from one block to the next
        cmd += ["--ado-origin", str(origin), "--ado-seed", str(ado_seed)]
    run(cmd, env)
    return b, d


def read_summary(d):
    with open(find_one(d, ["summary.tsv"])) as f:
        rows = [l.rstrip("\n").split("\t") for l in f if l.strip()]
    return rows[0], rows[1]


def vcf_lines(d):
    path = find_one(d, ["vcf", "vcf_dir/vcf.0001", "vcf.0001"])
    with open(path, "rb") as f:
        for raw in f:
            yield raw.decode(errors="replace")


def main():
    ap = argparse.ArgumentParser(usage=__doc__, add_help=True)
    ap.add_argument("--exe", required=True)
    ap.add_argument("--blocks", type=int, default=0)
    ap.add_argument("--block-size", type=int, default=200000)
    ap.add_argument("--out", default="chrom_out")
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--tree", default=None)
    ap.add_argument("--fasta", default=None)
    ap.add_argument("--region", default=None)
    ap.add_argument("--germline-vcf", default=None)
    ap.add_argument("--coverage-window", type=int, default=10000,
                    help="window (sites) of chrom_coverage.tsv, the mean depth per cell of the sites that are not in the VCF (0 = off)")
    ap.add_argument("--keep", action="store_true", help="keep the folders of the blocks")
    ap.add_argument("--gzip", action="store_true", help="write chrom.vcf.gz")
    ap.add_argument("opts", nargs=argparse.REMAINDER)
    a = ap.parse_args()
    opts = [o for o in a.opts if o != "--"]
    for o in opts:
        if re.match(r"^(-n|-l|-#|-o|-T|-U)", o) or o in ("--summary", "--coverage"):
            sys.exit("do not give %s: the script sets -n, -l, -#, -o, -T, -U, --summary and --coverage" % o)
    if a.germline_vcf and any(re.match(r"^-c(\d|\.|$)", o) for o in opts):
        sys.exit("germline SNPs are either simulated (-c) or given by the user (--germline-vcf), not both")
    if a.germline_vcf and not a.fasta:
        sys.exit("--germline-vcf needs --fasta")
    if not a.fasta and a.blocks < 1:
        sys.exit("give --blocks (or --fasta with a region)")
    global COVERAGE_WINDOW
    COVERAGE_WINDOW = max(0, a.coverage_window)
    exe = os.path.abspath(a.exe)
    out = os.path.abspath(a.out)
    os.makedirs(out, exist_ok=True)
    work = os.path.join(out, "blocks")
    os.makedirs(work, exist_ok=True)

    # reference sequence and blocks
    chrom_name, region_seq, region_start, germ = "1", None, 1, None
    info = []
    if a.fasta:
        if a.fasta.endswith(".gz"):  # decompress once (a plain file is needed to read a region)
            plain = os.path.join(work, "reference.fa")
            with gzip.open(a.fasta, "rb") as fi, open(plain, "wb") as fo:
                shutil.copyfileobj(fi, fo)
            a.fasta = plain
        name, s, e = parse_region(a.region) if a.region else (None, None, None)
        hdr = first_header(a.fasta)
        hm = re.match(r"^(.+):(\d+)-(\d+)$", hdr)   # a header such as >chr1:10000001-20000000 (region extracted from a genome)
        if name is None:
            name = hdr
        if hm and not os.path.exists(a.fasta + ".fai") and (name == hdr or name == hm.group(1)):
            h0 = int(hm.group(2))
            seq, _ = fetch_region(a.fasta, hdr, None, None)
            if s is None:
                region_seq, region_start = seq, h0
            else:
                if s < h0 or (e or s) > h0 + len(seq) - 1:
                    sys.exit("--region is outside the sequence %s" % hdr)
                region_seq, region_start = seq[s - h0: e - h0 + 1], s
            name = hm.group(1)
        else:
            region_seq, region_start = fetch_region(a.fasta, name, s, e)
        chrom_name = name
        if len(region_seq) == 0:
            sys.exit("the region is empty")
        nblocks = (len(region_seq) + a.block_size - 1) // a.block_size
        if a.blocks:
            nblocks = min(nblocks, a.blocks)
        total_len = len(region_seq)
        if a.germline_vcf:
            germ = read_germline(a.germline_vcf, name, region_start, region_start + len(region_seq) - 1)
        info.append("reference: %s %s:%d-%d (%d bp)" % (a.fasta, name, region_start, region_start + len(region_seq) - 1, len(region_seq)))
    else:
        nblocks = a.blocks
        total_len = a.blocks * a.block_size

    # 1. one tree for the whole chromosome (tumor only, ingroup tree, as -T reads it)
    tree = os.path.join(out, "chrom_tree.nwk")
    if a.tree:
        shutil.copy(a.tree, tree)
    else:
        tdir = os.path.join(work, "tree")
        shutil.rmtree(tdir, ignore_errors=True)
        topts, skip = [], False
        for o in opts:  # drop the options about reads and outputs, with their values
            if re.match(r"^(-C|-D|-E|-A|-P|-Q|-I|-V|-O|-B|-1|-2|-3|-4|-5|-7|-8|-9|-0|-Y|-z|-Z|-M|-x|-v|--)", o):
                skip = True
            elif o.startswith("-"):
                skip = False
                topts.append(o)
            elif not skip:
                topts.append(o)
        run([exe, "-n1", "-l100", "-6", "--no-outgroup", "-#%d" % a.seed, "-o%s" % tdir] + topts)
        shutil.copy(find_one(tdir, ["trees", "trees_dir/trees.0001"]), tree)
    if "--no-outgroup" not in opts:
        opts.append("--no-outgroup")
    if "--min-alt-reads" not in opts and not any(o.startswith("--min-alt-reads=") for o in opts):
        opts += ["--min-alt-reads", "2"]   # default of the block mode: one alternative read is not believed (see manual 8.23)
    if "--min-alt-vaf" not in opts and not any(o.startswith("--min-alt-vaf=") for o in opts):
        opts += ["--min-alt-vaf", "0.3"]   # and at least 30% of the reads of the cell (separates most single-cell errors)

    # 2. the blocks, in parallel
    rng = random.Random(a.seed)
    poss, genomes, masked, mismatches = {}, {}, 0, 0
    jobs_list = []
    for b in range(nblocks):
        if region_seq is None:
            jobs_list.append((b, a.block_size, None, b * a.block_size))
            continue
        off = b * a.block_size
        seq = region_seq[off: off + a.block_size]
        mat, pat, pos, nref = make_block_genome(seq, region_start + off, germ, rng)
        masked += len(seq) - len(pos)
        mismatches += nref
        if not pos:
            continue
        g = os.path.join(work, "genome_%04d.fa" % b)
        with open(g, "w") as fg:
            fg.write(">maternal\n%s\n>paternal\n%s\n" % (mat, pat))
        poss[b] = pos
        jobs_list.append((b, len(pos), g, off))
    done = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, a.jobs)) as ex:
        futs = [ex.submit(one_block, exe, opts, tree, b, n, a.seed, work, g, o, a.seed + 1) for b, n, g, o in jobs_list]
        for k, fu in enumerate(concurrent.futures.as_completed(futs), 1):
            b, d = fu.result()
            done[b] = d
            print("\rblock %d/%d" % (k, len(jobs_list)), end="", file=sys.stderr, flush=True)
    print(file=sys.stderr)

    # 3. merge: VCF with the positions of the chromosome, summary table with the sums
    opener = (lambda p: gzip.open(p, "wt")) if a.gzip else (lambda p: open(p, "w"))
    vcf_path = os.path.join(out, "chrom.vcf" + (".gz" if a.gzip else ""))
    header, rows, nsites, first = None, [], 0, True
    fcov = open(os.path.join(out, "chrom_coverage.tsv"), "w") if COVERAGE_WINDOW > 0 else None
    bed_dir = os.path.join(out, "chrom_coverage_bed")
    beds = {}
    with opener(vcf_path) as fv:
        for b in sorted(done):
            d = done[b]
            for line in vcf_lines(d):
                if line.startswith("#"):
                    if first:
                        line = re.sub(r"(##contig=<ID=)[^,>]*(,length=)\d+", r"\g<1>%s\g<2>%d" % (chrom_name, region_start + total_len - 1), line)
                        fv.write(line)
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or not f[1].isdigit():
                    continue  # replicate markers
                if region_seq is None:
                    pos = int(f[1]) + b * a.block_size
                else:
                    pos = poss[b][int(f[1]) - 1]
                f[0] = chrom_name
                f[1] = str(pos)
                f[2] = "snv%09d" % pos
                fv.write("\t".join(f) + "\n")
                nsites += 1
            first = False
            if fcov is not None:
                with open(find_one(d, ["coverage.tsv"])) as fc:
                    lines = [l.rstrip("\n").split("\t") for l in fc if l.strip()]
                if b == min(done):
                    fcov.write("chrom\tstart\tend\tn_sites\t" + "\t".join(lines[0][4:]) + "\n")
                    os.makedirs(bed_dir, exist_ok=True)
                    for c in lines[0][4:]:
                        beds[c] = open(os.path.join(bed_dir, c + ".bed"), "w")
                for f in lines[1:]:
                    i0, i1 = int(f[1]), int(f[2])
                    if region_seq is None:
                        g0, g1 = i0 + b * a.block_size, i1 + b * a.block_size
                    else:
                        g0, g1 = poss[b][i0 - 1], poss[b][i1 - 1]
                    fcov.write("%s\t%d\t%d\t%s\t%s\n" % (chrom_name, g0 - 1, g1, f[3], "\t".join(f[4:])))
                    for c, v in zip(lines[0][4:], f[4:]):
                        if v != "NA":
                            beds[c].write("%s\t%d\t%d\t%s\t%s\n" % (chrom_name, g0 - 1, g1, v, f[3]))
            h, r = read_summary(d)
            header = header or h
            rows.append((b, r))
    if fcov is not None:
        fcov.close()
        for fb in beds.values():
            fb.close()
    col = {c: i for i, c in enumerate(header)}
    inv_cols = [c for c in header if c.startswith("inv_")]
    sum_cols = [c for c in COUNT_COLS + inv_cols if c in col]
    total = ["NA"] * len(header)
    total[col["replicate"]] = "TOTAL"
    for c in sum_cols:
        total[col[c]] = str(sum(int(float(r[col[c]])) for _, r in rows))
    for c, wcol in WEIGHTED.items():
        if c in col and wcol in col:
            num = den = 0.0
            for _, r in rows:
                if r[col[c]] != "NA":
                    w = float(r[col[wcol]])
                    num += w * float(r[col[c]])
                    den += w
            if den > 0:
                total[col[c]] = "%.6f" % (num / den)
    for c in ("called_fraction", "genotype_error", "het_error", "hom_error", "het_to_hom", "singlet_error", "ado_simulated"):
        vals = [float(r[col[c]]) for _, r in rows if c in col and r[col[c]] != "NA"]
        if vals:
            total[col[c]] = "%.6f" % (sum(vals) / len(vals))
    with open(os.path.join(out, "chrom_summary.tsv"), "w") as fs:
        fs.write("block\t" + "\t".join(header) + "\n")
        for b, r in rows:
            fs.write("%d\t%s\n" % (b + 1, "\t".join(r)))
        fs.write("TOTAL\t" + "\t".join(total) + "\n")
    inv = sum(int(total[col[c]]) for c in inv_cols) if inv_cols else 0
    if region_seq is not None:
        info.append("sites that are not A, C, G or T (left out): %d" % masked)
        if a.germline_vcf:
            info.append("germline SNVs in the region: %d (REF different from the reference sequence: %d, ignored)" % (len(germ), mismatches))
        elif any(re.match(r"^-c(\d|\.|$)", o) for o in opts):
            info.append("germline SNPs simulated by CellCoal (-c)")
    info.append("variable sites in the VCF: %d; invariant sites (not in the VCF): %d" % (nsites, inv))
    with open(os.path.join(out, "chrom_info.txt"), "w") as fi:
        fi.write("\n".join(info) + "\n")
    if not a.keep:
        shutil.rmtree(work, ignore_errors=True)
    print("chromosome %s: %d bp in %d blocks" % (chrom_name, total_len, len(jobs_list)))
    print("\n".join(info))
    print("written to", out)


if __name__ == "__main__":
    main()
