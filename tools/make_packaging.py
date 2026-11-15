#!/usr/bin/env python3
"""Writes Formula/cellcoal.rb (Homebrew) and conda/meta.yaml from the templates in packaging/.

usage: make_packaging.py VERSION SOURCE_TARBALL
The SHA-256 of the GitHub source archive of the tag vVERSION is computed from SOURCE_TARBALL."""
import hashlib, pathlib, sys

version, tarball = sys.argv[1], sys.argv[2]
sha = hashlib.sha256(open(tarball, 'rb').read()).hexdigest()
root = pathlib.Path(__file__).resolve().parent.parent
for tpl, out in (('packaging/cellcoal.rb.in', 'Formula/cellcoal.rb'), ('packaging/meta.yaml.in', 'conda/meta.yaml')):
    text = (root / tpl).read_text().replace('@VERSION@', version).replace('@SHA256@', sha)
    (root / out).parent.mkdir(exist_ok=True)
    (root / out).write_text(text)
    print('wrote', out)
