#!/usr/bin/env python3
"""Benchmark single-cell variant callers on the reads simulated by CellCoal (--bam).

  # 1. simulate (the truth goes to reads_1/truth.vcf)
  cellcoal -s50 -l200000 -u2.4e-8 -e10000 -b1 -C10 -A0.01 0.0002 0 -E0.001 -D0.25 --ado-block 3000 \
           --gain-sd 0.4 --gain-len 8000 --amp-signature data/signatures/PTA_artifact_SCAN2.txt --bam --bulk-depth 30 -oout
  # 2a. evaluate any caller: a VCF with one column per cell (cell0001 ...) or one non-empty call per site (--vcf-sites)
  benchmark.py eval  out/reads_1/truth.vcf  caller.vcf
  # 2b. built-in pileup baseline on the BAM files, with and without PCR duplicates, to see what duplicates do
  benchmark.py pileup out/reads_1  --min-alt 2 --min-vaf 0.3

Truth: somatic sites of truth.vcf (germline sites are ignored). A call is a true positive if the cell carries the allele,
a false positive otherwise. Needs numpy (pileup mode only)."""
import argparse, os, sys, collections

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def read_truth(path):
    """-> (cells, {pos: (ref, alts, {cell_index: set(alt bases carried)})}) for somatic sites"""
    cells, sites = [], {}
    for line in open(path):
        if line.startswith('##'):
            continue
        f = line.rstrip('\n').split('\t')
        if line.startswith('#'):
            cells = f[9:]
            continue
        if 'GERMLINE' in f[7]:
            continue
        pos, ref, alts = int(f[1]), f[3], f[4].split(',')
        carry = {}
        for i, gt in enumerate(f[9:]):
            s = set()
            for a in gt.replace('|', '/').split('/'):
                if a.isdigit() and int(a) > 0:
                    s.add(alts[int(a) - 1])
            if s:
                carry[i] = s
        if carry:
            sites[pos] = (ref, alts, carry)
    return cells, sites


def metrics(tp, fp, fn):
    p = tp / (tp + fp) if tp + fp else float('nan')
    r = tp / (tp + fn) if tp + fn else float('nan')
    return p, r


def evaluate(truth_cells, sites, calls, label):
    """calls: set of (pos, cell_index, alt). Cell-level and site-level precision/recall."""
    truth = {(pos, ci, a) for pos, (ref, alts, carry) in sites.items() for ci, s in carry.items() for a in s}
    tp = len(calls & truth)
    fp = len(calls - truth)
    fn = len(truth - calls)
    p, r = metrics(tp, fp, fn)
    cs = {pos for pos, ci, a in calls}
    ts = set(sites)
    sp, sr = metrics(len(cs & ts), len(cs - ts), len(ts - cs))
    # singletons: sites carried by exactly one cell
    single = {(pos, ci, a) for pos, (ref, alts, carry) in sites.items() if len(carry) == 1 for ci, s in carry.items() for a in s}
    sr1 = len(calls & single) / len(single) if single else float('nan')
    print('%-28s cell-level: TP %6d FP %6d FN %6d  precision %.3f recall %.3f | sites: precision %.3f recall %.3f | singleton recall %.3f' % (
        label, tp, fp, fn, p, r, sp, sr, sr1))
    return p, r


def cmd_eval(a):
    cells, sites = read_truth(a.truth)
    calls = set()
    names = None
    for line in open(a.vcf):
        if line.startswith('##'):
            continue
        f = line.rstrip('\n').split('\t')
        if line.startswith('#'):
            names = f[9:]
            continue
        pos = int(f[1])
        alts = f[4].split(',')
        fmt = f[8].split(':')
        gi = fmt.index('GT') if 'GT' in fmt else None
        for k, col in enumerate(f[9:]):
            if gi is None:
                continue
            gt = col.split(':')[gi].replace('|', '/').split('/')
            ci = cells.index(names[k]) if names[k] in cells else k
            for g in gt:
                if g.isdigit() and int(g) > 0 and int(g) <= len(alts):
                    calls.add((pos, ci, alts[int(g) - 1]))
    evaluate(cells, sites, calls, os.path.basename(a.vcf))


def cmd_pileup(a):
    import numpy as np
    import bamcheck as b
    cells, sites = read_truth(os.path.join(a.dir, 'truth.vcf'))
    name, ref = b.read_fasta(os.path.join(a.dir, 'ref.fa'))
    L = len(ref)
    code = np.full(256, 4, np.uint8)
    for i, c in enumerate('ACGT'):
        code[ord(c)] = i
    refc = code[np.frombuffer(ref.encode(), np.uint8)]
    calls = {True: set(), False: set()}      # True = remove duplicates
    for ci, cell in enumerate(cells):
        path = os.path.join(a.dir, cell + '.bam')
        _, _, recs = b.parse_bam(path)
        for dedup in (False, True):
            cnt = np.zeros((5, L), np.int32)
            for r in recs:
                if dedup and r['flag'] & 1024:
                    continue
                cnt[code[np.frombuffer(r['bases'].encode(), np.uint8)], np.arange(r['pos'], r['pos'] + len(r['bases']))] += 1
            depth = cnt[:4].sum(axis=0)
            altc = cnt[:4].copy()
            altc[refc[np.arange(L)] % 4, np.arange(L)] = 0
            altbase = altc.argmax(axis=0)
            altn = altc.max(axis=0)
            sel = np.nonzero((altn >= a.min_alt) & (altn >= a.min_vaf * np.maximum(depth, 1)))[0]
            for pos0 in sel:
                calls[dedup].add((int(pos0) + 1, ci, 'ACGT'[altbase[pos0]]))
        sys.stderr.write('\r%s' % cell)
    sys.stderr.write('\n')
    evaluate(cells, sites, calls[False], 'pileup, all reads')
    evaluate(cells, sites, calls[True], 'pileup, PCR duplicates out')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    e = sub.add_parser('eval')
    e.add_argument('truth')
    e.add_argument('vcf')
    p = sub.add_parser('pileup')
    p.add_argument('dir')
    p.add_argument('--min-alt', type=int, default=2)
    p.add_argument('--min-vaf', type=float, default=0.3)
    a = ap.parse_args()
    (cmd_eval if a.cmd == 'eval' else cmd_pileup)(a)


if __name__ == '__main__':
    main()
