#!/usr/bin/env python3
"""Fuzzes the command line and the --config parser of CellCoal.

Usage: fuzz.py BINARY [--runs N] [--seed S] [--time SECONDS]
Every run uses a small data set plus random (often invalid) options and values, or a randomly damaged copy of
examples/config.example.yaml / its JSON form. Invalid input must be rejected with a message and a non-zero exit status;
what must never happen is a crash (signal), a sanitizer report (ASan/UBSan "runtime error") or a hang (timeout).
Failing inputs are printed so they can be reproduced. Needs only the standard library."""
import json, os, random, re, subprocess, sys, tempfile, time

args = sys.argv[1:]
binary = args[0]
def opt(name, default):
    return type(default)(args[args.index(name) + 1]) if name in args else default
runs, seed, budget = opt("--runs", 300), opt("--seed", 1), opt("--time", 0.0)
rng = random.Random(seed)
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
cfg_text = open(os.path.join(root, "examples", "config.example.yaml")).read()

# option letters taken from the program itself (the usage line)
usage = subprocess.run([binary, "-?"], capture_output=True, text=True, timeout=20).stdout + \
        subprocess.run([binary, "-?"], capture_output=True, text=True, timeout=20).stderr
letters = sorted(set(re.findall(r"-([A-Za-z0-9#])[#\s]", usage.split("Usage:")[-1].split("\n")[0]))) or list("nslegCVEDA")
long_opts = ["--fast", "--bcf", "--cov-ado-corr", "--config"]

NASTY = ["", "0", "-1", "1", "0.5", "1e-6", "1e308", "-1e308", "nan", "inf", "99999999999", "-99999999999", "abc", "1e", ".", "--", "%s%n", "0x10", "1,2", "ÿ"]
BASE = ["-n1", "-s8", "-l200", "-u1e-6", "-e1000", "-b1", "-C5", "-V1", "-E0.001"]

def value():
    r = rng.random()
    if r < 0.5: return rng.choice(NASTY)
    if r < 0.8: return "%g" % (10 ** rng.uniform(-8, 3) * rng.choice([1, -1, 1, 1]))
    return str(rng.randint(-5, 400))

def random_cmd():
    cmd = list(BASE)
    for _ in range(rng.randint(1, 6)):
        r = rng.random()
        if r < 0.12:
            o = rng.choice(long_opts)
            cmd += [o] + ([value()] if o in ("--cov-ado-corr",) else [])
        else:
            cmd.append("-" + rng.choice(letters) + value())
            if rng.random() < 0.3: cmd.append(value())
    return cmd

def damage(text):
    t = list(text)
    for _ in range(rng.randint(1, 6)):
        k = rng.random(); i = rng.randrange(max(1, len(t)))
        if k < 0.3 and t: del t[i]
        elif k < 0.55: t.insert(i, rng.choice(list(":-{}[],\"'#\n\t ") + [chr(rng.randint(1, 255)), "9999999999", "nan"]))
        elif k < 0.7 and t: t[i] = rng.choice("0123456789.-eE")
        elif k < 0.85: t = t[:i]
        else:
            j = rng.randrange(max(1, len(t))); t[i:i] = t[j:j + rng.randint(1, 40)]
    return "".join(t)

def yaml_to_json(text):
    d = {}
    for l in text.splitlines():
        l = l.split("#")[0].strip()
        if ":" in l:
            k, v = [x.strip() for x in l.split(":", 1)]
            if v in ("true", "false"): d[k] = (v == "true")
            else:
                try: d[k] = float(v) if re.search(r"[.eE]", v) else int(v)
                except ValueError: d[k] = v
    return json.dumps(d)

bad = []
def run(cmd, label, cwd, extra=None):
    try:
        p = subprocess.run(cmd, cwd=cwd, capture_output=True, timeout=30, text=True, errors="replace")
    except subprocess.TimeoutExpired:
        bad.append((label, "timeout (hang?)", cmd, extra)); return
    out = p.stdout + p.stderr
    if p.returncode < 0:
        bad.append((label, "killed by signal %d" % -p.returncode, cmd, extra))
    else:
        # leak reports are ignored: the copy of the arguments made once per process is intentional
        m = re.search(r"^.*(ERROR: AddressSanitizer|runtime error|Segmentation fault|stack smashing).*$", out, re.M)
        if m:
            bad.append((label, "sanitizer report: " + m.group(0)[:200], cmd, extra))

t0 = time.time(); n = 0
with tempfile.TemporaryDirectory() as tmp:
    while n < runs and (budget <= 0 or time.time() - t0 < budget):
        d = os.path.join(tmp, "r%d" % n); os.mkdir(d)
        if n % 3 == 2:    # configuration file
            as_json = rng.random() < 0.4
            body = yaml_to_json(cfg_text) if as_json else cfg_text
            if rng.random() < 0.9: body = damage(body)
            path = os.path.join(d, "c.json" if as_json else "c.yaml")
            open(path, "w", errors="replace").write(body)
            run([binary, "--config", path, "-n1", "-s6", "-l100"], "config", d, body)
        else:
            run([binary] + random_cmd(), "args", d)
        n += 1
if bad:
    print("FUZZ: %d of %d runs misbehaved (seed %d)" % (len(bad), n, seed))
    for label, why, cmd, extra in bad[:10]:
        print(" -", label, why, "\n    ", " ".join(cmd)); 
        if extra: print("    --- config ---\n" + extra[:600])
    sys.exit(1)
print("FUZZ ok: %d runs (seed %d), no crashes, hangs or sanitizer reports" % (n, seed))
