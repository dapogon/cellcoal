#!/usr/bin/env python3
"""Checks the documentation against the program:
 1. every command `cellcoal ...` of the manual and the README is valid (`--check`), and is run for real when it is small;
 2. every long option of `--help` is described in the manual.
Usage: python3 tests/manual_check.py bin/cellcoal-x.y.z   (from the root of the repository)"""
import os, re, shlex, subprocess, sys, tempfile

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
binary = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else sys.exit(__doc__)
docs = {"manual": open(os.path.join(root, "docs/manual/cellcoal-manual.md"), encoding="utf-8").read(), "README": open(os.path.join(root, "README.md"), encoding="utf-8").read()}
SKIP = ("hg38", "chr21.fa", "chr1.fa", ".fa.gz")          # examples that need a real genome that is not in the repository
bad = 0

def commands(text):
    out, lines, i = [], text.split("\n"), 0
    while i < len(lines):
        l = lines[i].strip()
        if l.startswith("cellcoal -") or l.startswith("cellcoal --"):
            cmd = l
            while cmd.endswith("\\") and i + 1 < len(lines):
                i += 1
                cmd = cmd[:-1] + " " + lines[i].strip()
            out.append(cmd)
        i += 1
    return out

def size(args):
    def val(flag, default):
        for a in args:
            if re.fullmatch(r"-%s\d+" % flag, a): return int(a[2:])
        return default
    return val("n", 50) * val("s", 20) * val("l", 1000)

with tempfile.TemporaryDirectory() as d:
    for name in ("examples", "data", "tools"):
        os.symlink(os.path.join(root, name), os.path.join(d, name))
    seen = set()
    for src, text in docs.items():
        for cmd in commands(text):
            cmd = cmd.replace("`", "")
            if cmd in seen or any(k in cmd for k in SKIP): continue
            seen.add(cmd)
            args = shlex.split(cmd)[1:]
            r = subprocess.run([binary, "--check"] + args, cwd=d, capture_output=True, text=True, timeout=60)
            if r.returncode != 0:
                bad += 1; print("INVALID in %s: %s\n   %s" % (src, cmd, (r.stdout + r.stderr).strip().split("\n")[-1][:200])); continue
            if size(args) <= 3_000_000 and "--sweep" not in cmd:
                r = subprocess.run([binary] + args, cwd=d, capture_output=True, text=True, timeout=300)
                if r.returncode != 0:
                    bad += 1; print("FAILS in %s: %s\n   %s" % (src, cmd, (r.stdout + r.stderr).strip().split("\n")[-1][:200]))
    print("commands checked: %d" % len(seen))

help_text = subprocess.run([binary, "--help"], capture_output=True, text=True).stdout
manual = docs["manual"]
opts = sorted(set(re.findall(r"(?<![\w-])--[a-z][a-z0-9-]+", help_text)))
IGNORE = {"--help", "--version"}
missing = [o for o in opts if o not in manual and o not in IGNORE]
for o in missing:
    bad += 1; print("NOT IN THE MANUAL: " + o)
print("long options in --help: %d, missing from the manual: %d" % (len(opts), len(missing)))
sys.exit(1 if bad else 0)
