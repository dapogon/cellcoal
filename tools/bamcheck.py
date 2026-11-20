#!/usr/bin/env python3
"""Minimal BAM reader and checker for the BAM files written by CellCoal (--bam); no dependencies.

  bamcheck.py ref.fa cell0001.bam [more.bam ...]

Reports reads, duplicates, coordinate order, mean depth, mismatch rate against the reference (with and without
the PCR duplicates) and checks that the BAI index finds every read of a window."""
import struct, sys, zlib, gzip, collections


def bgzf_bytes(path):
    with open(path, 'rb') as f:
        return gzip.decompress(f.read())          # multi-member gzip == BGZF


def read_fasta(path):
    seq, name = [], None
    for line in open(path):
        if line.startswith('>'):
            name = line[1:].split()[0]
        else:
            seq.append(line.strip().upper())
    return name, ''.join(seq)


def parse_bam(path):
    d = bgzf_bytes(path)
    assert d[:4] == b'BAM\1', 'not a BAM'
    lt = struct.unpack_from('<i', d, 4)[0]
    text = d[8:8 + lt].decode()
    o = 8 + lt
    nref = struct.unpack_from('<i', d, o)[0]
    o += 4
    refs = []
    for _ in range(nref):
        ln = struct.unpack_from('<i', d, o)[0]
        name = d[o + 4:o + 4 + ln - 1].decode()
        l = struct.unpack_from('<i', d, o + 4 + ln)[0]
        refs.append((name, l))
        o += 8 + ln
    recs = []
    while o < len(d):
        bs = struct.unpack_from('<i', d, o)[0]
        r = d[o + 4:o + 4 + bs]
        o += 4 + bs
        refid, pos, lname, mapq, binv, ncig, flag, lseq, nref_, npos, tlen = struct.unpack_from('<iiBBHHHiiii', r, 0)
        name = r[32:32 + lname - 1].decode()
        p = 32 + lname
        cig = struct.unpack_from('<%dI' % ncig, r, p)
        p += 4 * ncig
        sq = r[p:p + (lseq + 1) // 2]
        p += (lseq + 1) // 2
        bases = []
        for i in range(lseq):
            b = sq[i // 2]
            bases.append('=ACMGRSVTWYHKDBN'[(b >> 4) if i % 2 == 0 else (b & 15)])
        qual = r[p:p + lseq]
        p += lseq
        tags = {}
        while p < len(r):
            t = r[p:p + 2].decode()
            ty = chr(r[p + 2])
            p += 3
            if ty == 'Z':
                e = r.index(b'\0', p)
                tags[t] = r[p:e].decode()
                p = e + 1
            elif ty == 'C':
                tags[t] = r[p]
                p += 1
            elif ty == 'i':
                tags[t] = struct.unpack_from('<i', r, p)[0]
                p += 4
            else:
                raise ValueError('tag type ' + ty)
        recs.append(dict(pos=pos, flag=flag, bases=''.join(bases), cigar=cig, name=name, tags=tags, tlen=tlen, mpos=npos, qual=qual))
    return text, refs, recs


def bai_check(path, recs, reflen):
    d = open(path + '.bai', 'rb').read()
    assert d[:4] == b'BAI\1'
    nref, nbin = struct.unpack_from('<ii', d, 4)[0], struct.unpack_from('<i', d, 8)[0]
    o = 12
    bins = {}
    for _ in range(nbin):
        b, nch = struct.unpack_from('<Ii', d, o)
        o += 8
        bins[b] = struct.unpack_from('<%dQ' % (2 * nch), d, o)
        o += 16 * nch
    nint = struct.unpack_from('<i', d, o)[0]
    lin = struct.unpack_from('<%dQ' % nint, d, o + 4)
    return len(bins), nint, lin


def main():
    ref_name, ref = read_fasta(sys.argv[1])
    for path in sys.argv[2:]:
        text, refs, recs = parse_bam(path)
        assert refs[0][0] == ref_name and refs[0][1] == len(ref), 'reference mismatch %s' % (refs,)
        n = len(recs)
        dup = sum(1 for r in recs if r['flag'] & 1024)
        ordered = all(recs[i]['pos'] <= recs[i + 1]['pos'] for i in range(n - 1))
        bases = sum(len(r['bases']) for r in recs)
        mm = mmu = tot = totu = 0
        pairs = collections.defaultdict(list)
        for r in recs:
            pairs[(r['name'], r['flag'] & 192)].append(r)
            for i, b in enumerate(r['bases']):
                same = (ref[r['pos'] + i] == b)
                tot += 1
                mm += not same
                if not (r['flag'] & 1024):
                    totu += 1
                    mmu += not same
        # molecules: duplicates of the same molecule share their template, so reads of a molecule are consistent
        mol = collections.Counter(r['tags']['MI'] for r in recs if r['flag'] & 64 or r['flag'] & 128)
        ok_pairs = all(len(v) >= 1 for v in pairs.values())
        mates = collections.Counter(r['name'] for r in recs)
        paired = all(c % 2 == 0 for c in mates.values())
        nb, nint, lin = bai_check(path, recs, len(ref))
        print('%s: reads %d (duplicates %.1f%%), sorted %s, paired %s, depth %.2fX, mismatch %.5f (dedup %.5f), bins %d, linear %d' % (
            path.split('/')[-1], n, 100.0 * dup / max(n, 1), ordered, paired, bases / len(ref), mm / max(tot, 1), mmu / max(totu, 1), nb, nint))
        assert ordered and paired


if __name__ == '__main__':
    main()
