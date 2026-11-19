# Unreleased (after 2.2.7)

## Fixes
- `-T` (user tree): the number of cells is now the number of leaves of the tree, so `-s` no longer has to match it (before, a tree with more leaves than `-s` crashed).
- Branch lengths of the trees written with `-6` use `%.10g` instead of `%8.6f`, so very short branches are no longer rounded to 0 when the tree is reused with `-T`. (Only the text of the tree files changes; simulated data are identical.)
- User tree reader: a branch length shorter than the previous one no longer inherits leftover characters.
- `-T`: when `-s` differs from the number of leaves of the tree, a note says that `-s` is ignored (also shown by `--check`).

## New
- `--no-outgroup` (also `no_outgroup: true` in a configuration file): simulate the tumor only. The healthy outgroup cell is not in the tree file, the VCF/BCF, catg, pileup, genotype files or the statistics, and it does not decide whether a site is variable. The reference is still the healthy root genome (all 0), and the trunk (`-k`) is kept, so truncal mutations are still variants of every cell. The tree written with `-6` has the ingroup only, and `-T` reads it back as it is.
