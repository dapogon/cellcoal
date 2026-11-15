#!/usr/bin/env python3
"""Summarizes gcov output: line coverage per source file of src/.

Usage: gcov_summary.py OBJDIR   (the directory with the .gcno/.gcda files of a --coverage build)
Runs gcov on every object, reads the .gcov files and prints a Markdown table (also suitable for $GITHUB_STEP_SUMMARY)
with the executed/executable lines of each file and the total. Standard library only."""
import glob, os, re, subprocess, sys, tempfile

obj = os.path.abspath(sys.argv[1])
tot_exec = tot_lines = 0
rows = []
with tempfile.TemporaryDirectory() as tmp:
    # gcov finds the sources by the paths recorded at compile time (src/x.c, relative to the repository root)
    os.symlink(os.path.abspath("src"), os.path.join(tmp, "src"))
    for gcda in sorted(glob.glob(os.path.join(obj, "*.gcda"))):
        subprocess.run(["gcov", "-o", obj, gcda], cwd=tmp, capture_output=True)
    for g in sorted(glob.glob(os.path.join(tmp, "*.c.gcov"))):
        name = None; ex = ln = 0
        for l in open(g, errors="replace"):
            parts = l.split(":", 2)
            if len(parts) < 3: continue
            cnt, num = parts[0].strip(), parts[1].strip()
            if num == "0":
                m = re.match(r"Source:(.*)", parts[2].strip())
                if m: name = os.path.basename(m.group(1))
                continue
            if cnt == "-": continue
            ln += 1
            if cnt not in ("#####", "====="): ex += 1
        if name and name.endswith(".c") and ln:
            rows.append((name, ex, ln)); tot_exec += ex; tot_lines += ln
print("| file | lines run | lines | coverage |\n|---|---:|---:|---:|")
for name, ex, ln in rows:
    print("| %s | %d | %d | %.1f%% |" % (name, ex, ln, 100.0 * ex / ln))
if tot_lines:
    print("| **total** | %d | %d | **%.1f%%** |" % (tot_exec, tot_lines, 100.0 * tot_exec / tot_lines))
