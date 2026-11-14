#!/usr/bin/env python3
"""Robinson-Foulds distance between two Newick trees (no external dependencies).

Usage: rf_distance.py TRUE_TREE INFERRED_TREE [--drop NAME]

The leaf NAME (default: outgcell, the CellCoal outgroup) is removed from both trees
before comparing the unrooted topologies. Prints: cells, RF, normalized RF, correct splits.
"""
import re
import sys


def parse(newick):
    s = newick.strip().rstrip(';')
    pos = 0

    def skip_length():
        nonlocal pos
        if pos < len(s) and s[pos] == ':':
            pos += re.match(r':[^,()]*', s[pos:]).end()

    def node():
        nonlocal pos
        if s[pos] == '(':
            pos += 1
            children = [node()]
            while s[pos] == ',':
                pos += 1
                children.append(node())
            assert s[pos] == ')', 'unbalanced parentheses'
            pos += 1
            pos += re.match(r'[^:,()]*', s[pos:]).end()      # optional internal label
            skip_length()
            return children
        name = re.match(r'[^:,()]+', s[pos:])
        pos += name.end()
        skip_length()
        return name.group(0)

    return node()


def leaves(t):
    return {t} if isinstance(t, str) else set().union(*[leaves(c) for c in t])


def prune(t, drop):
    if isinstance(t, str):
        return None if t == drop else t
    kids = [k for k in (prune(c, drop) for c in t) if k is not None]
    if not kids:
        return None
    return kids[0] if len(kids) == 1 else kids


def splits(t, all_leaves):
    """Non-trivial bipartitions, each written as the side that does not contain the smallest leaf name."""
    out, ref = set(), min(all_leaves)

    def rec(t):
        if isinstance(t, str):
            return {t}
        side = set().union(*[rec(c) for c in t])
        if 1 < len(side) < len(all_leaves) - 1:
            out.add(frozenset(side) if ref not in side else frozenset(all_leaves - side))
        return side

    rec(t)
    return out


def rf_distance(true_newick, inferred_newick, drop='outgcell'):
    a = prune(parse(true_newick), drop)
    b = prune(parse(inferred_newick), drop)
    la, lb = leaves(a), leaves(b)
    if la != lb:
        raise ValueError('the trees have different leaves')
    sa, sb = splits(a, la), splits(b, lb)
    rf = len(sa ^ sb)
    return len(la), rf, rf / max(1, len(sa) + len(sb)), len(sa & sb), len(sa)


if __name__ == '__main__':
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    drop = sys.argv[sys.argv.index('--drop') + 1] if '--drop' in sys.argv else 'outgcell'
    n, rf, nrf, ok, tot = rf_distance(open(sys.argv[1]).read(), open(sys.argv[2]).read(), drop)
    print(f'{n}\t{rf}\t{nrf:.3f}\t{ok}/{tot}')
