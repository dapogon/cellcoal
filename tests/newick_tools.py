"""Pair statistics of a tree in Newick (cellNNNN leaves), shared by the tests of cellforward."""


def newick_pair_times(nwk, n, distance=False):
    """time to the coalescence of every pair of cells (cell0001 ... celln) of a tree in Newick with branch lengths in generations"""
    parent, length, name = [-1], [0.0], [None]
    cur, i = 0, 0
    while i < len(nwk) and nwk[i] != ";":
        c = nwk[i]
        if c == "(":
            parent.append(cur); length.append(0.0); name.append(None)
            cur = len(parent) - 1
        elif c == ",":
            parent.append(parent[cur]); length.append(0.0); name.append(None)
            cur = len(parent) - 1
        elif c == ")":
            cur = parent[cur]
        elif c == ":":
            j = i + 1
            while nwk[j] not in ",);":
                j += 1
            length[cur] = float(nwk[i + 1:j])
            i = j - 1
        else:
            j = i
            while nwk[j] not in ",):":
                j += 1
            name[cur] = nwk[i:j]
            i = j - 1
        i += 1
    depth = [0.0] * len(parent)             # distance to the root
    for v in range(1, len(parent)):         # parents are numbered before their children
        depth[v] = depth[parent[v]] + length[v]
    leaf = {int(name[v][4:]) - 1: v for v in range(len(parent)) if name[v]}
    anc = {}
    for c, v in leaf.items():
        chain = []
        while v >= 0:
            chain.append(v)
            v = parent[v]
        anc[c] = chain
    times = {}
    end = max(depth[v] for v in leaf.values())
    for x in range(n):
        for y in range(x + 1, n):
            sx = set(anc[x])
            m = next(v for v in anc[y] if v in sx)      # the nearest common ancestor
            times[(x, y)] = (depth[leaf[x]] + depth[leaf[y]] - 2 * depth[m]) if distance else end - depth[m]
    return times
