#!/usr/bin/env python3
"""chrom.py: simulate a long chromosome with CellCoal by blocks (a thin wrapper of the native --blocks mode).

One coalescent tree is generated once and every block of sites is simulated on that same tree (no recombination, so the sites are
independent given the tree).  For each block only the variable sites go to the VCF; the rest of the block is summarized by the number
of sites for each genotype of the healthy root (inv_AA ... inv_TT in the summary table).  The simulation itself is done by the binary
(--blocks, --fasta, --region, --germline-vcf; see the manual, section 8.29); this script sets the options, adds the defaults of the
block mode and writes the files of earlier versions.

Usage (random sequence):
  python3 tools/chrom.py --exe bin/cellcoal-2.5.2 --blocks 50 --block-size 200000 --out chr1 [--seed 7] [--gzip]
                         -- <CellCoal options: -s -u -e -g -C -D -E --min-alt-reads ...>

Usage (real reference sequence):
  python3 tools/chrom.py --exe bin/cellcoal-2.5.2 --fasta hg38.fa.gz --region chr1:10000001-20000000 --out chr1 -- <options>

  --fasta FILE      reference in standard FASTA (plain or .gz, several lines per sequence, upper or lower case).  Sites that are not
                    A, C, G or T (N...) are left out of the simulation: they are not in the VCF and not in the invariant counts.
  --region R        chrom:start-end, 1-based and inclusive (default: the whole first sequence)
  --germline-vcf F  germline SNPs given by the user (VCF, optionally .gz; SNVs only; GT 0|1, 1|0, 1|1, 0/1, 1/1).  Unphased
                    heterozygous genotypes are phased at random.  Not together with -c.  At those sites the VCF has the reference
                    base as REF and the flag GERMLINE.
  --tree FILE       tree to use (otherwise it is generated)
  --jobs N          threads (OpenMP) of the binary

Do not give -n, -l, -#, -o, -T, -U or --summary after "--": the script sets them.  It also sets --no-outgroup and, unless you give
them, --min-alt-reads 2 --min-alt-vaf 0.3 (a site is written only if a cell has a non-reference call with at least 2 alternative reads
and 30% of its reads).  Outputs, in the --out folder:
  chrom.vcf[.gz]        variable sites of all the blocks (CHROM and POS of the chromosome)
  chrom_summary.tsv     one row per block, plus a TOTAL row
  chrom_tree.nwk        the tree used for every block
  chrom_coverage.tsv    BED-like: mean depth of every cell over the sites that are NOT in the VCF, per window of --coverage-window
                        sites (default 10000; 0 = off); n_sites = non-variable sites of the window
  chrom_coverage_bed/   one BED file per cell (chrom, start, end, mean depth, n_sites) with the same windows
  chrom_info.txt        region, masked sites, germline sites used
"""
import argparse, gzip, os, re, shutil, subprocess, sys


def run(cmd, env=None):
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env)
    if r.returncode != 0:
        sys.exit("command failed (%d): %s\n%s" % (r.returncode, " ".join(cmd), r.stdout[-2000:]))
    return r.stdout


def find_one(folder, names):
    for n in names:
        p = os.path.join(folder, n)
        if os.path.exists(p):
            return p
    sys.exit("%s: none of %s was written" % (folder, names))


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
    ap.add_argument("--coverage-window", type=int, default=10000)
    ap.add_argument("--keep", action="store_true", help="keep the folder of the binary (<out>/work)")
    ap.add_argument("--gzip", action="store_true", help="write chrom.vcf.gz")
    a, opts = ap.parse_known_args()
    opts = [o for o in opts if o != "--"]
    if a.fasta is None and a.blocks < 1:
        sys.exit("give --blocks N (random sequence) or --fasta")
    if a.germline_vcf and any(re.match(r"^-c(\d|\.|$)", o) for o in opts):
        sys.exit("germline SNPs are either simulated (-c) or given by the user (--germline-vcf), not both")
    if a.germline_vcf and not a.fasta:
        sys.exit("--germline-vcf needs --fasta")
    exe = os.path.abspath(a.exe)
    out = os.path.abspath(a.out)
    os.makedirs(out, exist_ok=True)
    work = os.path.join(out, "work")
    shutil.rmtree(work, ignore_errors=True)

    cmd = [exe, "-l%d" % a.block_size, "--summary", "-#%d" % a.seed, "-o%s" % work, "--no-outgroup", "-6"]
    if a.fasta:
        cmd += ["--fasta", a.fasta]
        if a.region:
            cmd += ["--region", a.region]
        if a.germline_vcf:
            cmd += ["--germline-vcf", a.germline_vcf]
        if a.blocks:
            cmd += ["--blocks", str(a.blocks)]
    else:
        cmd += ["--blocks", str(a.blocks), "--contig", "1"]
    if a.tree:
        cmd += ["-T%s" % os.path.abspath(a.tree)]
    if "--min-alt-reads" not in opts:
        cmd += ["--min-alt-reads", "2"]   # one alternative read is not believed (manual 8.23)
    if "--min-alt-vaf" not in opts:
        cmd += ["--min-alt-vaf", "0.3"]   # and at least 30% of the reads of the cell
    if a.coverage_window > 0:
        cmd += ["--coverage", str(a.coverage_window)]
    if a.gzip:
        cmd += ["-z"]
    cmd += opts
    run(cmd, dict(os.environ, OMP_NUM_THREADS=str(max(1, a.jobs))))

    # files of the binary -> files of chrom.py
    vcf = find_one(work, ["vcf.gz", "vcf"])
    shutil.copy(vcf, os.path.join(out, "chrom.vcf.gz" if vcf.endswith(".gz") else "chrom.vcf"))
    with open(find_one(work, ["summary.tsv"])) as f, open(os.path.join(out, "chrom_summary.tsv"), "w") as g:
        lines = f.read().splitlines()
        g.write(re.sub(r"^replicate", "block", lines[0]) + "\n")
        g.write("\n".join(lines[1:]) + "\n")
    tree = os.path.join(work, "trees")
    if not os.path.exists(tree):
        tree = find_one(work, ["trees_dir/trees.0001", "trees_dir/trees"])
    shutil.copy(tree, os.path.join(out, "chrom_tree.nwk"))
    shutil.copy(find_one(work, ["blocks_info.txt"]), os.path.join(out, "chrom_info.txt"))
    chrom = "1"
    for l in open(os.path.join(out, "chrom_info.txt")):
        m = re.match(r"blocks: \d+ of \d+ sites; chromosome (\S+),", l)
        if m:
            chrom = m.group(1)
    if a.coverage_window > 0:
        bed_dir = os.path.join(out, "chrom_coverage_bed")
        shutil.rmtree(bed_dir, ignore_errors=True)
        os.makedirs(bed_dir)
        beds = {}
        with open(find_one(work, ["coverage.tsv"])) as f, open(os.path.join(out, "chrom_coverage.tsv"), "w") as g:
            head = f.readline().rstrip("\n").split("\t")
            cells = head[4:]
            g.write("chrom\tstart\tend\tn_sites\t" + "\t".join(cells) + "\n")
            for c in cells:
                beds[c] = open(os.path.join(bed_dir, c + ".bed"), "w")
            for l in f:
                r = l.rstrip("\n").split("\t")
                start, end = int(r[1]) - 1, int(r[2])
                g.write("%s\t%d\t%d\t%s\t%s\n" % (chrom, start, end, r[3], "\t".join(r[4:])))
                for c, v in zip(cells, r[4:]):
                    if v != "NA":
                        beds[c].write("%s\t%d\t%d\t%s\t%s\n" % (chrom, start, end, v, r[3]))
        for b in beds.values():
            b.close()
    if not a.keep:
        shutil.rmtree(work, ignore_errors=True)
    print("\n".join(l.rstrip() for l in open(os.path.join(out, "chrom_info.txt"))))
    print("written to", out)


if __name__ == "__main__":
    main()
