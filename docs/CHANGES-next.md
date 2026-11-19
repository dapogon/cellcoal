# Unreleased (after 2.2.7)

## Fixes
- `-T` (user tree): the number of cells is now the number of leaves of the tree, so `-s` no longer has to match it (before, a tree with more leaves than `-s` crashed).
- Branch lengths of the trees written with `-6` use `%.10g` instead of `%8.6f`, so very short branches are no longer rounded to 0 when the tree is reused with `-T`. (Only the text of the tree files changes; simulated data are identical.)
- User tree reader: a branch length shorter than the previous one no longer inherits leftover characters.
