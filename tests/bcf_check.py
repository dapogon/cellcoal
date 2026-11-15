#!/usr/bin/env python3
"""Checks a BCF file written by CellCoal (--bcf) against the VCF text of the same run.

Usage: bcf_check.py FILE.bcf FILE.vcf
Decodes the BGZF/BCF 2.2 file following the specification (independently of the C converter) and compares every
record, INFO and FORMAT value with the VCF text: strings and integers exactly, floats as float32 of the printed text.
Also checks the BGZF structure (BC extra field, end-of-file marker). Needs only the standard library."""
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

nrec = 0
while r.p < len(r.b):
    ls, li = r.u32(), r.u32()
    sh = R(r.take(ls)); ind = R(r.take(li))
    chrom, pos, rlen = sh.i32(), sh.i32(), sh.i32()
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
print("BCF ok: %d records, %d samples, %d BGZF blocks" % (nrec, len(samples), nblocks))
