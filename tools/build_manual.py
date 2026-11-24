#!/usr/bin/env python3
"""Builds the HTML (self-contained) and PDF versions of the manual from its single source, docs/manual/cellcoal-manual.md.

Usage: python3 tools/build_manual.py          (from the root of the repository)

Needs pandoc. The PDF is made with headless Chromium through playwright (python) when available; otherwise only the HTML is built.
"""
import os, re, subprocess, sys, tempfile

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
src = os.path.join(root, "docs", "manual", "cellcoal-manual.md")
html = os.path.join(root, "docs", "cellcoal.manual.v2.0.html")
pdf = os.path.join(root, "docs", "cellcoal.manual.v2.0.pdf")

css = """
body { font-family: -apple-system, 'Segoe UI', Helvetica, Arial, sans-serif; line-height: 1.5; max-width: 62rem; margin: 2rem auto; padding: 0 1rem; color: #222; }
h1, h2, h3, h4 { line-height: 1.25; } h1 { border-bottom: 1px solid #ddd; padding-bottom: .3rem; margin-top: 2.2rem; }
code { background: #f3f3f3; padding: 0 .25em; border-radius: 3px; font-size: 92%; } pre { background: #f6f6f6; padding: .8rem; overflow-x: auto; }
pre code { background: none; padding: 0; } table { border-collapse: collapse; } td, th { border: 1px solid #ccc; padding: .3rem .6rem; }
img { max-width: 100%; height: auto; }
div.caption { font-size: 85%; line-height: 1.4; color: #444; background: #f7f7f9; border-left: 3px solid #999; padding: .4rem .8rem; margin: .3rem 0 1.4rem; }
div.caption p { margin: 0; } nav#TOC { background: #fafafa; border: 1px solid #eee; padding: .5rem 1.5rem; margin-bottom: 2rem; }
@media print { body { max-width: none; margin: 0; } h1 { page-break-before: auto; } img { page-break-inside: avoid; } }
"""
with tempfile.NamedTemporaryFile("w", suffix=".css", delete=False) as f:
    f.write(css)
    cssfile = f.name
# the title, author and date are in the YAML block at the top of the Markdown
subprocess.check_call(["pandoc", src, "-f", "markdown+tex_math_dollars+link_attributes-implicit_figures", "-t", "html5", "--standalone", "--embed-resources",
                       "--resource-path", os.path.dirname(src), "--toc", "--toc-depth=2", "--mathml", "-c", cssfile,
                       "--metadata", "pagetitle=CellCoal manual", "-o", html])
os.unlink(cssfile)
print("wrote", html)
try:
    from playwright.sync_api import sync_playwright
except ImportError:
    print("playwright is not installed: PDF not built")
    sys.exit(0)
with sync_playwright() as p:
    b = p.chromium.launch()
    pg = b.new_page(viewport={"width": 1100, "height": 900})
    pg.goto("file://" + html)
    pg.wait_for_timeout(1500)
    pg.pdf(path=pdf, format="A4", print_background=True, display_header_footer=True, header_template="<span></span>", footer_template='<div style="font-size:8px;color:#666;width:100%;text-align:center;font-family:Helvetica,Arial,sans-serif">CellCoal manual &#8212; <span class="pageNumber"></span> / <span class="totalPages"></span></div>', margin={"top": "15mm", "bottom": "18mm", "left": "12mm", "right": "12mm"})
    b.close()
print("wrote", pdf)
