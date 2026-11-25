#!/usr/bin/env python3
"""Small summaries of single-cell sequencing data, to calibrate CellCoal (the --wga presets and the targeted sequencing).

Only the Python standard library is needed (no samtools, no pysam); the output is small, so only the summaries have to leave the cluster.

  BAM (sorted, aligned):   python3 tools/seq_summary.py bam  cell.bam  --out summaries/ [--window 100000] [--stride 1]
  FASTQ (R1 [R2]):         python3 tools/seq_summary.py fastq cell_R1.fastq.gz [cell_R2.fastq.gz] --out summaries/ [--max-reads 5000000]

For a BAM it writes <name>.summary.json and <name>.windows.tsv.gz (reads, bases, duplicates and reads with MAPQ >= 20 per window of
the autosomes). The summary has the depth, the fraction of the genome without reads, the uniformity of the depth (coefficient of variation,
Gini, MAPD), the duplicates, the MAPQ, the insert size, the length of the reads, the soft clipping, the indels and the mismatches (tag NM).
For FASTQ files it writes <name>.fastq.json: the length of the reads, the base quality, the GC, the Ns and the fraction of duplicates
(the same first 40 bases in R1 and the same first 40 in R2, among the first reads).
--stride K reads every K-th record of the BAM (the counts are multiplied by K): faster, for a first look.
"""
import argparse, gzip, json, math, os, re, struct, sys, collections

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("kind", choices=["bam", "fastq"])
ap.add_argument("files", nargs="+")
ap.add_argument("--out", default=".")
ap.add_argument("--window", type=int, default=100000, help="BAM: window of the depth, in bases")
ap.add_argument("--stride", type=int, default=1, help="BAM: use one record of every STRIDE")
ap.add_argument("--max-reads", type=int, default=5000000, help="FASTQ: reads to read")
ap.add_argument("--chroms", default="1-22", help="BAM: chromosomes of the depth, e.g. 1-22 or 1-22,X")
args = ap.parse_args()
os.makedirs(args.out, exist_ok=True)


def chrom_wanted(name):
    n = re.sub(r"^chr", "", name, flags=re.I)
    for part in args.chroms.split(","):
        if "-" in part:
            a, b = part.split("-")
            if n.isdigit() and int(a) <= int(n) <= int(b):
                return True
        elif n == part:
            return True
    return False


def gini(v):
    v = sorted(v); n = len(v); s = sum(v)
    if n == 0 or s == 0:
        return float("nan")
    return (2.0 * sum((i + 1) * x for i, x in enumerate(v)) / (n * s)) - (n + 1.0) / n


def quant(v, q):
    v = sorted(v)
    return v[min(len(v) - 1, int(q * len(v)))] if v else float("nan")


def stats(v):
    n = len(v)
    if n == 0:
        return {}
    m = sum(v) / n
    sd = math.sqrt(sum((x - m) ** 2 for x in v) / n)
    return {"mean": m, "sd": sd, "median": quant(v, 0.5)}


def read_exact(f, n):
    b = f.read(n)
    while len(b) < n:
        c = f.read(n - len(b))
        if not c:
            break
        b += c
    return b


def summarize_bam(path):
    name = re.sub(r"\.bam$", "", os.path.basename(path))
    f = gzip.open(path, "rb")
    if f.read(4) != b"BAM\1":
        sys.exit("%s is not a BAM file" % path)
    lt = struct.unpack("<i", f.read(4))[0]; f.read(lt)
    nref = struct.unpack("<i", f.read(4))[0]
    refs = []
    for _ in range(nref):
        ln = struct.unpack("<i", f.read(4))[0]
        nm = f.read(ln)[:-1].decode(); L = struct.unpack("<i", f.read(4))[0]
        refs.append((nm, L))
    W = args.window
    wins = {}		# ref index -> [reads, bases, dups, q20]
    for i, (nm, L) in enumerate(refs):
        if chrom_wanted(nm):
            nw = (L + W - 1) // W
            wins[i] = [[0, 0, 0, 0] for _ in range(nw)]
    tot = collections.Counter()
    ins = collections.Counter()
    mq = collections.Counter()
    rlen = collections.Counter()
    unp = struct.Struct("<iiBBHHHiiii")
    k = 0
    while True:
        h = f.read(4)
        if len(h) < 4:
            break
        bs = struct.unpack("<i", h)[0]
        rec = read_exact(f, bs)
        k += 1
        if args.stride > 1 and k % args.stride:
            continue
        ref, pos, lrn, mapq, bn, ncig, flag, lseq, nref2, npos, tlen = unp.unpack_from(rec, 0)
        S = args.stride
        tot["records"] += S
        if flag & 4:
            tot["unmapped"] += S
            continue
        if flag & 0x900:
            continue
        tot["primary_mapped"] += S
        if flag & 1:
            tot["paired"] += S
        if flag & 0x400:
            tot["duplicates"] += S
        mq[min(mapq, 60) // 10 * 10] += S
        rlen[lseq] += S
        co = 32 + lrn
        aligned = soft = ins_n = del_n = ins_b = del_b = 0
        for c in range(ncig):
            v = struct.unpack_from("<I", rec, co + 4 * c)[0]; n = v >> 4; op = v & 15
            if op in (0, 7, 8):
                aligned += n
            elif op == 4:
                soft += n
            elif op == 1:
                ins_n += 1; ins_b += n
            elif op == 2:
                del_n += 1; del_b += n
        tot["aligned_bases"] += aligned * S; tot["soft_clipped_bases"] += soft * S
        tot["insertions"] += ins_n * S; tot["deletions"] += del_n * S
        tot["inserted_bases"] += ins_b * S; tot["deleted_bases"] += del_b * S
        # tags: NM
        t = co + 4 * ncig + (lseq + 1) // 2 + lseq
        nm_val = None
        while t + 3 <= len(rec):
            tag = rec[t:t + 2]; ty = chr(rec[t + 2]); t += 3
            if ty in "cC":
                val = rec[t] if ty == "C" else struct.unpack_from("<b", rec, t)[0]; sz = 1
            elif ty in "sS":
                val = struct.unpack_from("<H" if ty == "S" else "<h", rec, t)[0]; sz = 2
            elif ty in "iI":
                val = struct.unpack_from("<I" if ty == "I" else "<i", rec, t)[0]; sz = 4
            elif ty == "f":
                val = None; sz = 4
            elif ty == "A":
                val = None; sz = 1
            elif ty == "Z" or ty == "H":
                e = rec.index(b"\0", t); val = None; sz = e - t + 1
            elif ty == "B":
                sub = chr(rec[t]); cnt = struct.unpack_from("<i", rec, t + 1)[0]
                sz = 5 + cnt * {"c": 1, "C": 1, "s": 2, "S": 2, "i": 4, "I": 4, "f": 4}[sub]; val = None
            else:
                break
            if tag == b"NM":
                nm_val = val; break
            t += sz
        if nm_val is not None:
            tot["nm_sum"] += nm_val * S; tot["nm_reads"] += S
        if (flag & 2) and (flag & 0x40) and 0 < tlen < 3000:
            ins[tlen] += S
        if ref in wins and pos >= 0:
            w = wins[ref][min(pos // W, len(wins[ref]) - 1)]
            w[0] += S; w[1] += lseq * S
            if flag & 0x400:
                w[2] += S
            if mapq >= 20:
                w[3] += S
    f.close()
    # the windows
    rows = []
    for i in sorted(wins):
        nm, L = refs[i]
        for j, w in enumerate(wins[i]):
            span = min(W, L - j * W)
            rows.append((nm, j * W, span, w))
    full = [(r[3][1] / r[2]) for r in rows if r[2] == W]	# depth of the full windows
    med = quant(full, 0.5) if full else 0.0
    cvv = stats(full)
    lr = [math.log2(max(x, 1e-3) / med) for x in full] if med > 0 else []
    mapd = quant([abs(lr[i + 1] - lr[i]) for i in range(len(lr) - 1)], 0.5) if len(lr) > 1 else float("nan")
    out = {
        "file": os.path.basename(path), "window": W, "stride": args.stride, "chroms": args.chroms,
        "records": tot["records"], "primary_mapped": tot["primary_mapped"], "unmapped_fraction": tot["unmapped"] / max(1, tot["records"]),
        "paired_fraction": tot["paired"] / max(1, tot["primary_mapped"]), "duplicate_fraction": tot["duplicates"] / max(1, tot["primary_mapped"]),
        "mapq_histogram": {str(k): v for k, v in sorted(mq.items())},
        "read_length": stats([x for x, c in rlen.items() for _ in range(1)]) if False else {"mode": (rlen.most_common(1)[0][0] if rlen else 0), "mean": sum(x * c for x, c in rlen.items()) / max(1, sum(rlen.values()))},
        "soft_clipped_fraction": tot["soft_clipped_bases"] / max(1, tot["soft_clipped_bases"] + tot["aligned_bases"]),
        "insertions_per_Mb": 1e6 * tot["insertions"] / max(1, tot["aligned_bases"]), "deletions_per_Mb": 1e6 * tot["deletions"] / max(1, tot["aligned_bases"]),
        "mean_insertion_length": tot["inserted_bases"] / max(1, tot["insertions"]), "mean_deletion_length": tot["deleted_bases"] / max(1, tot["deletions"]),
        "mismatches_per_aligned_base(NM)": tot["nm_sum"] / max(1, tot["nm_reads"]) / max(1.0, tot["aligned_bases"] / max(1, tot["primary_mapped"])),
        "insert_size": {},
        "depth_windows": {"n": len(full), "mean": cvv.get("mean"), "median": med, "cv": (cvv["sd"] / cvv["mean"]) if full and cvv["mean"] else None,
                            "gini": gini(full), "mapd": mapd,
                            "fraction_below_10pct_of_median": (sum(1 for x in full if x < 0.1 * med) / len(full)) if full and med > 0 else None,
                            "fraction_zero": (sum(1 for x in full if x == 0) / len(full)) if full else None,
                            "p5_over_median": (quant(full, 0.05) / med) if med > 0 else None, "p95_over_median": (quant(full, 0.95) / med) if med > 0 else None},
    }
    if ins:
        n = sum(ins.values()); m = sum(k * c for k, c in ins.items()) / n
        sd = math.sqrt(sum((k - m) ** 2 * c for k, c in ins.items()) / n)
        acc = 0; med_i = 0
        for k in sorted(ins):
            acc += ins[k]
            if acc >= n / 2:
                med_i = k; break
        out["insert_size"] = {"n": n, "mean": m, "sd": sd, "median": med_i}
    json.dump(out, open(os.path.join(args.out, name + ".summary.json"), "w"), indent=1)
    with gzip.open(os.path.join(args.out, name + ".windows.tsv.gz"), "wt") as g:
        g.write("chrom\tstart\tspan\treads\tbases\tduplicates\treads_q20\n")
        for nm, s, span, w in rows:
            g.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\n" % (nm, s, span, w[0], w[1], w[2], w[3]))
    print(name, "records", tot["records"], "median window depth %.2f" % med, "dup %.3f" % out["duplicate_fraction"])


def fq_records(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        while True:
            h = f.readline()
            if not h:
                return
            s = f.readline().rstrip("\n"); f.readline(); q = f.readline().rstrip("\n")
            yield s, q


def summarize_fastq(paths):
    name = re.sub(r"([_.]R?1)?\.(fastq|fq)(\.gz)?$", "", os.path.basename(paths[0]))
    n = 0; lens = collections.Counter(); gc = ns = bases = 0; qsum = 0; q30 = 0
    seen = {}; dup = 0
    its = [fq_records(p) for p in paths]
    while n < args.max_reads:
        try:
            recs = [next(it) for it in its]
        except StopIteration:
            break
        n += 1
        s, q = recs[0]
        lens[len(s)] += 1; bases += len(s)
        gc += s.count("G") + s.count("C"); ns += s.count("N")
        qs = [ord(c) - 33 for c in q]
        qsum += sum(qs); q30 += sum(1 for x in qs if x >= 30)
        key = tuple(r[0][:40] for r in recs)
        if key in seen:
            dup += 1
        else:
            seen[key] = 1
    out = {"files": [os.path.basename(p) for p in paths], "reads_examined": n, "read_length_mode": lens.most_common(1)[0][0] if lens else 0,
           "gc_fraction": gc / max(1, bases - ns), "n_fraction": ns / max(1, bases), "mean_base_quality": qsum / max(1, bases), "q30_fraction": q30 / max(1, bases),
           "duplicate_fraction_first40": dup / max(1, n), "note": "duplicates counted among the examined reads only: it underestimates the fraction of the whole file"}
    json.dump(out, open(os.path.join(args.out, name + ".fastq.json"), "w"), indent=1)
    print(name, "reads", n, "dup(first40) %.3f" % out["duplicate_fraction_first40"])


if args.kind == "bam":
    for p in args.files:
        summarize_bam(p)
else:
    summarize_fastq(args.files)
