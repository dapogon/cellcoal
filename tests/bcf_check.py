#!/usr/bin/env python3
"""Checks a BCF file written by CellCoal (--bcf) against the VCF text of the same run.

Usage: bcf_check.py FILE.bcf FILE.vcf
Decodes the BGZF/BCF 2.2 file following the specification (independently of the C converter) and compares every
record, INFO and FORMAT value with the VCF text: strings and integers exactly, floats as float32 of the printed text.
Also checks the BGZF structure (BC extra field, end-of-file marker) and the CSI index FILE.bcf.csi: every record must lie in a chunk
of its bin (virtual offsets computed independently from the BGZF blocks), with consistent loffsets and statistics.
Needs only the standard library."""
import gzip, struct, sys

def bgzf_blocks_ok(path):
    data = open(path, "rb").read()
    pos = n = 0
    while pos < len(data):
        if data[pos:pos+4] != b"\x1f\x8b\x08\x04" or data[pos+12:pos+14] != b"BC":
            sys.exit("bad BGZF block at %d" % pos)
        bsize = struct.unpack_from("<H", data, pos + 16)[0] + 1
        pos += bsize; n += 1
    if not data.endswith(bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")):
        sys.exit("missing BGZF end-of-file marker")
    return n

class R:
    def __init__(s, b): s.b = b; s.p = 0
    def take(s, n): v = s.b[s.p:s.p+n]; s.p += n; return v
    def u32(s): return struct.unpack("<I", s.take(4))[0]
    def i32(s): return struct.unpack("<i", s.take(4))[0]

def typed_len(r):
    d = r.take(1)[0]
    t, n = d & 0xF, d >> 4
    if n == 15:
        n = typed_int(r)
    return t, n

def typed_int(r):
    t, n = typed_len(r); assert n == 1
    return vec(r, t, 1)[0]

MISSING_I = {1: -128, 2: -32768, 3: -2**31}
def vec(r, t, n):
    if t == 0: return []
    if t == 1: v = list(struct.unpack("<%db" % n, r.take(n)))
    elif t == 2: v = list(struct.unpack("<%dh" % n, r.take(2*n)))
    elif t == 3: v = list(struct.unpack("<%di" % n, r.take(4*n)))
    elif t == 5: v = [struct.unpack("<f", r.take(4))[0] if True else 0 for _ in range(n)]
    elif t == 7: return r.take(n).decode().rstrip("\0")
    else: raise ValueError(t)
    return v

def is_end(t, x, raw=None):
    return (t == 1 and x == -127) or (t == 2 and x == -32767) or (t == 3 and x == -2**31 + 1)

def f32(text):
    return struct.unpack("<f", struct.pack("<f", float(text)))[0]

bcf, vcf = sys.argv[1], sys.argv[2]
nblocks = bgzf_blocks_ok(bcf)
r = R(gzip.open(bcf).read())
assert r.take(5) == b"BCF\x02\x02", "magic"
ltext = r.u32(); text = r.take(ltext).decode()
assert text.endswith("\0")
hdr = [l for l in text.rstrip("\0").split("\n") if l]
samples = hdr[-1].split("\t")[9:]
dic = {}
for l in hdr:
    if l.startswith(("##INFO=", "##FORMAT=", "##FILTER=")):
        key = l.split("ID=")[1].split(",")[0].rstrip(">")
        idx = int(l.split("IDX=")[1].rstrip(">"))
        kind = l.split("Type=")[1].split(",")[0] if "Type=" in l else None
        dic[idx] = (key, kind)

vrecs = []
for l in open(vcf):
    if l[0] == "#" or "\t" not in l: continue
    vrecs.append(l.rstrip("\n").split("\t"))

def block_table(path):
    data = open(path, "rb").read()
    pos = 0; tab = []; u = 0
    while pos < len(data):
        bsize = struct.unpack_from("<H", data, pos + 16)[0] + 1
        isize = struct.unpack_from("<I", data, pos + bsize - 4)[0]
        tab.append((pos, u, isize)); u += isize; pos += bsize
    return tab

blocks = block_table(bcf)
def voff(u):
    for caddr, ustart, size in blocks:          # the block holding the byte u (empty blocks hold nothing)
        if size and ustart <= u < ustart + size:
            return (caddr << 16) | (u - ustart)
    return (blocks[-1][0] << 16)                 # end of the data: start of the end-of-file block

size_of = {caddr: size for caddr, _, size in blocks}
next_of = {blocks[i][0]: blocks[i+1][0] for i in range(len(blocks) - 1)}
def norm(v):
    """(block, offset) with offset == block size is the same place as (next block, 0), as htslib reports it"""
    caddr, off = v >> 16, v & 0xFFFF
    if caddr in size_of and off == size_of[caddr] and caddr in next_of:
        return next_of[caddr] << 16
    return v

recs = []
nrec = 0
while r.p < len(r.b):
    ustart = r.p
    ls, li = r.u32(), r.u32()
    sh = R(r.take(ls)); ind = R(r.take(li))
    chrom, pos, rlen = sh.i32(), sh.i32(), sh.i32()
    recs.append((ustart, r.p, chrom, pos, rlen))
    qual = struct.unpack("<I", sh.take(4))[0]
    ninfo = struct.unpack("<H", sh.take(2))[0]; nallele = struct.unpack("<H", sh.take(2))[0]
    w = sh.u32(); nsamp, nfmt = w & 0xFFFFFF, w >> 24
    t, n = typed_len(sh); idv = vec(sh, t, n)
    alleles = []
    for _ in range(nallele):
        t, n = typed_len(sh); alleles.append(vec(sh, t, n))
    t, n = typed_len(sh); filt = vec(sh, t, n) if t else []
    v = vrecs[nrec]; nrec += 1
    assert (chrom, pos + 1, rlen) == (0, int(v[1]), 1), ("pos", v[:3])
    assert (idv or ".") == v[2] and alleles[0] == v[3], ("id/ref", v[:4])
    assert (",".join(alleles[1:]) or ".") == v[4], ("alt", v[:5])
    assert qual == 0x7F800001 and v[5] == ".", "qual"
    assert (";".join(dic[i][0] for i in filt) or ".") == v[6], ("filter", v[6])
    info = {}
    for _ in range(ninfo):
        key = dic[typed_int(sh)][0]; t, n = typed_len(sh)
        info[key] = True if t == 0 else vec(sh, t, n)
    for item in v[7].split(";"):
        k, _, val = item.partition("=")
        if not val:
            assert info[k] is True, ("flag", k)
        elif k == "AF":
            assert [f32(x) for x in val.split(",")] == info[k], ("AF", val, info[k])
        elif k == "NS":
            assert info[k] == [int(val)], ("NS", val)
        else:
            assert info[k] == val, (k, val, info[k])
    assert nsamp == len(samples) and nfmt == len(v[8].split(":")), "counts"
    fkeys = v[8].split(":")
    for fk_i, fk in enumerate(fkeys):
        key = dic[typed_int(ind)][0]; assert key == fk, ("fmt key", key, fk)
        t, n = typed_len(ind)
        for s in range(nsamp):
            txt = v[9 + s].split(":")[fk_i]
            if key == "GT":
                assert t == 1, "GT must be an integer vector"
                g = vec(ind, t, n)
                text_gt = "".join(("." if (x >> 1) == 0 else str((x >> 1) - 1)) + ("" if i == n - 1 else ("|" if g[i + 1] & 1 else "/")) for i, x in enumerate(g))
                assert text_gt == txt, ("GT", s, g, txt)
                continue
            if t == 7:
                got = ind.take(n).decode().rstrip("\0"); assert got == txt, (key, s, got, txt)
                continue
            vals = vec(ind, t, n)
            want = txt.split(",")
            got = []
            for x in vals:
                if t == 5:
                    raw = struct.unpack("<I", struct.pack("<f", x))[0]
                    if raw == 0x7F800002: break
                    got.append(None if raw == 0x7F800001 else x)
                else:
                    if is_end(t, x): break
                    got.append(None if x == MISSING_I[t] else x)
            for g, wv in zip(got, want):
                if wv == ".":
                    assert g is None, (key, s, got, want)
                elif t == 5:
                    assert g == f32(wv), (key, s, g, wv)
                else:
                    assert g == int(wv), (key, s, g, wv)
            assert len(got) == len(want) or (want == ["."] and got == [None]), (key, s, got, want)
assert nrec == len(vrecs), (nrec, len(vrecs))

# ---- CSI index ----
csi = bcf + ".csi"
bgzf_blocks_ok(csi)
c = R(gzip.open(csi).read())
assert c.take(4) == b"CSI\x01", "CSI magic"
min_shift, depth, l_aux, nref = c.i32(), c.i32(), c.i32(), None
c.take(l_aux); nref = c.i32()
ncontig = len([l for l in hdr if l.startswith("##contig=")])
assert nref == ncontig, ("n_ref", nref, ncontig)
assert (1 << (min_shift + 3 * depth)) > max(p + 1 for _, _, _, p, _ in recs), "depth does not cover the contig"
meta_bin = ((1 << (3 * (depth + 1))) - 1) // 7 + 1
bins = []
for _ in range(nref):
    d = {}
    for _ in range(c.i32()):
        b = c.u32(); loff = struct.unpack("<Q", c.take(8))[0]; chunks = []
        for _ in range(c.i32()):
            chunks.append(struct.unpack("<QQ", c.take(16)))
        d[b] = (loff, chunks)
    bins.append(d)
tail = c.take(8) if c.p < len(c.b) else b""     # optional count of records without coordinate
assert c.p == len(c.b), "trailing bytes in the index"

def reg2bin(beg, end):
    end -= 1; s_ = min_shift; t = ((1 << (3 * depth)) - 1) // 7
    for l in range(depth, 0, -1):
        if beg >> s_ == end >> s_: return t + (beg >> s_)
        s_ += 3; t -= 1 << (3 * l)
    return 0

firstvo = {}
for ustart, uend, chrom, pos, rlen in recs:
    b = reg2bin(pos, pos + max(rlen, 1))
    loff, chunks = bins[chrom][b]                 # KeyError = the bin is missing
    vs, ve = norm(voff(ustart)), norm(voff(uend))
    assert any(norm(cb) <= vs and ve <= norm(ce) for cb, ce in chunks), ("record not in a chunk", pos)
    firstvo.setdefault((chrom, b), vs)
for chrom in range(nref):
    for b, (loff, chunks) in bins[chrom].items():
        if b == meta_bin:
            continue
        assert loff <= firstvo[(chrom, b)], ("loffset", b)
        assert chunks == sorted(chunks) and all(chunks[i][1] <= chunks[i+1][0] for i in range(len(chunks) - 1)), ("chunks", b)
    n_recs = [x for x in recs if x[2] == chrom]
    if n_recs:
        (_, (m0, m1)), = [(None, bins[chrom][meta_bin][1][1])]
        assert (m0, m1) == (len(n_recs), 0), ("n_mapped/unmapped", m0, m1)
        assert tuple(norm(x) for x in bins[chrom][meta_bin][1][0]) == (norm(voff(n_recs[0][0])), norm(voff(n_recs[-1][1]))), "offset range"
print("CSI ok: %d bins" % sum(len(d) for d in bins))
print("BCF ok: %d records, %d samples, %d BGZF blocks" % (nrec, len(samples), nblocks))
