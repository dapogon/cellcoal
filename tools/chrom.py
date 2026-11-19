#!/usr/bin/env python3
"""chrom.py: simulate a long chromosome with CellCoal by blocks (prototype).

One coalescent tree is generated once and every block of sites is simulated on that same tree (no recombination, so the
sites are independent given the tree).  For each block only the variable sites go to the VCF; the rest of the block is
summarized by the number of sites for each genotype of the healthy root (inv_AA ... inv_TT in the summary table).

Usage:
  python3 tools/chrom.py --exe bin/cellcoal-2.2.7 --blocks 50 --block-size 200000 --out chr1 [--jobs 4] [--seed 7]
                         [--keep] [--gzip] -- <CellCoal options: -s -u -e -g -C -D -E --min-alt-reads ...>

Do not give -n, -l, -#, -o, -T or --summary after "--": the script sets them.  A tree can be given with
--tree FILE (it is then not generated).  Outputs, in the --out folder:
  chrom.vcf[.gz]        variable sites of all the blocks (POS = position in the chromosome)
  chrom_summary.tsv     one row per block, plus a TOTAL row
  chrom_tree.nwk        the tree used for every block
"""
import argparse, concurrent.futures, glob, gzip, os, re, shutil, subprocess, sys

COUNT_COLS = ["mutations", "deletions", "cnloh", "snv_true", "snv_observed", "singleton_sites"]
WEIGHTED = {"snv_precision": "snv_observed", "snv_recall": "snv_true", "singleton_recall": "singleton_sites"}


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


def one_block(exe, opts, tree, b, size, seed, workdir):
    d = os.path.join(workdir, "block_%04d" % b)
    shutil.rmtree(d, ignore_errors=True)
    env = dict(os.environ, OMP_NUM_THREADS="1")
    cmd = [exe, "-n1", "-l%d" % size, "-T%s" % tree, "--summary", "-#%d" % (seed + b + 1), "-o%s" % d] + opts
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
    ap.add_argument("--blocks", type=int, required=True)
    ap.add_argument("--block-size", type=int, default=200000)
    ap.add_argument("--out", default="chrom_out")
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--tree", default=None)
    ap.add_argument("--keep", action="store_true", help="keep the folders of the blocks")
    ap.add_argument("--gzip", action="store_true", help="write chrom.vcf.gz")
    ap.add_argument("opts", nargs=argparse.REMAINDER)
    a = ap.parse_args()
    opts = [o for o in a.opts if o != "--"]
    for o in opts:
        if re.match(r"^(-n|-l|-#|-o|-T)", o) or o == "--summary":
            sys.exit("do not give %s: the script sets -n, -l, -#, -o, -T and --summary" % o)
    exe = os.path.abspath(a.exe)
    out = os.path.abspath(a.out)
    os.makedirs(out, exist_ok=True)
    work = os.path.join(out, "blocks")
    os.makedirs(work, exist_ok=True)

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

    # 2. the blocks, in parallel
    done = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, a.jobs)) as ex:
        futs = [ex.submit(one_block, exe, opts, tree, b, a.block_size, a.seed, work) for b in range(a.blocks)]
        for k, fu in enumerate(concurrent.futures.as_completed(futs), 1):
            b, d = fu.result()
            done[b] = d
            print("\rblock %d/%d" % (k, a.blocks), end="", file=sys.stderr, flush=True)
    print(file=sys.stderr)

    # 3. merge: VCF with the positions of the chromosome, summary table with the sums
    opener = (lambda p: gzip.open(p, "wt")) if a.gzip else (lambda p: open(p, "w"))
    vcf_path = os.path.join(out, "chrom.vcf" + (".gz" if a.gzip else ""))
    total_len = a.blocks * a.block_size
    header, rows = None, []
    nsites = 0
    with opener(vcf_path) as fv:
        for b in range(a.blocks):
            d = done[b]
            offset = b * a.block_size
            for line in vcf_lines(d):
                if line.startswith("#"):
                    if b == 0:
                        line = re.sub(r"(##contig=<ID=[^,>]*,length=)\d+", r"\g<1>%d" % total_len, line)
                        fv.write(line)
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or not f[1].isdigit():
                    continue  # replicate markers
                pos = int(f[1]) + offset
                f[1] = str(pos)
                f[2] = "snv%09d" % pos
                fv.write("\t".join(f) + "\n")
                nsites += 1
            h, r = read_summary(d)
            header = header or h
            rows.append(r)
    col = {c: i for i, c in enumerate(header)}
    inv_cols = [c for c in header if c.startswith("inv_")]
    sum_cols = [c for c in COUNT_COLS + inv_cols if c in col]
    total = ["NA"] * len(header)
    total[col["replicate"]] = "TOTAL"
    for c in sum_cols:
        total[col[c]] = str(sum(int(float(r[col[c]])) for r in rows))
    for c, wcol in WEIGHTED.items():
        if c in col and wcol in col:
            num = den = 0.0
            for r in rows:
                if r[col[c]] != "NA":
                    w = float(r[col[wcol]])
                    num += w * float(r[col[c]])
                    den += w
            if den > 0:
                total[col[c]] = "%.6f" % (num / den)
    for c in ("called_fraction", "genotype_error", "het_error", "hom_error", "het_to_hom", "singlet_error", "ado_simulated"):
        vals = [float(r[col[c]]) for r in rows if c in col and r[col[c]] != "NA"]
        if vals:
            total[col[c]] = "%.6f" % (sum(vals) / len(vals))
    with open(os.path.join(out, "chrom_summary.tsv"), "w") as fs:
        fs.write("block\t" + "\t".join(header) + "\n")
        for b, r in enumerate(rows):
            fs.write("%d\t%s\n" % (b + 1, "\t".join(r)))
        fs.write("TOTAL\t" + "\t".join(total) + "\n")
    if not a.keep:
        shutil.rmtree(work, ignore_errors=True)
    inv = sum(int(total[col[c]]) for c in inv_cols) if inv_cols else 0
    print("chromosome: %d sites in %d blocks; VCF: %d variable sites; invariant (not in VCF): %d" % (total_len, a.blocks, nsites, inv))
    print("written to", out)


if __name__ == "__main__":
    main()
