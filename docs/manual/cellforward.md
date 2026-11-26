# cellforward: a tumor simulated forward in time

`cellforward` is the companion of CellCoal. CellCoal simulates the genealogy of a sample of cells backward in time (the coalescent); `cellforward` simulates the whole tumor forward in time, one cell division at a time, and writes the tree of a random sample of cells. CellCoal then simulates the mutations, copy number changes, amplification, reads and variant calls on that tree (`-T` or `--tree file=`).

Use the coalescent when the tumor is very large or you need thousands of replicates in seconds, and `cellforward` when the model has things that the coalescent cannot do: the fitness of a clone that comes from the driver mutations that it really got, random sizes and extinctions without approximations, or a check of the coalescent. The cost is proportional to the number of cells of the tumor: about 0.14 s for a tumor of 10<sup>6</sup> cells (the limit is 5·10<sup>7</sup> cells, because every division is kept in memory).

## The model

- The tumor starts from **one cell** at time 0. Every cell **divides** at the birth rate `--birth` and **dies** at the rate `--death`, both per generation. A tumor that dies out is simulated again, so the tumors that are kept are those that survive.
- The tumor stops when it has **`--size` cells** (every tumor takes a different time) or after **`--time` generations** (every tumor has a different size). Exactly one of the two is given.
- **Drivers.** In every division, each daughter cell gets a new driver mutation with the probability `--driver-rate`. A driver multiplies the birth rate of the cell, and of all its descendants, by `1 + S` (`--driver-effect`): `S` is the mean of an exponential variable (`--driver-dist exponential`, the default) or the same for all the drivers (`--driver-dist fixed`). The cell with the driver starts a **clone**, and the descendants of the clone that get more drivers start subclones. Nothing about a clone is an input: its fitness, size and age are a result.
- **Clones started by the user** (`--clone-time`, `--clone-fitness`, `--clone-parent`), as `--clones` of CellCoal but forward: at the time given (generations after the start of the tumor) a random cell of the parent clone gets a growth advantage, and its descendants are the clone. The fitness is the **growth rate relative to the first clone** as in CellCoal: the birth rate of the clone is `death + F (birth - death)` (F = 1.5 grows 50% faster). Only the birth and the time are chosen; the size and the frequency of the clone are a result, and they change among replicates. The clones are numbered 1, 2, ... in the order of the lists (the first clone, without drivers, is 0), and can also get random drivers. A tumor in which the parent clone has no cell when the clone has to start, or that reaches `--size` before the time of the clone, is simulated again (so the tumors kept are those in which the clone exists).
- **Sample.** `-s` cells are chosen at random from the tumor, without replacement, and their tree is written.
- Times are in **generations**. The branch lengths of the tree are the generations multiplied by `--mutation-rate` (default 1), because CellCoal reads the branch lengths of a user tree as the expected number of mutations per site: give here the mutation rate per site and generation that you would give to CellCoal with `-u`.

## Use

```
cellforward --birth 0.5 --death 0.1 --size 100000 --driver-rate 1e-5 --driver-effect 0.1 -s 30 -n 100 --mutation-rate 1e-8 -o fwd
cellcoal -n100 -s30 -l10000 -Tfwd/trees.nwk
```

| Option | Meaning | Default |
|---|---|---|
| `--birth B` | birth rate per generation | needed |
| `--death D` | death rate per generation | 0 |
| `--size N` / `--time T` | stop at N cells / at T generations (one of them) | needed |
| `--driver-rate P` | probability that a daughter cell gets a driver | 0 (no drivers) |
| `--driver-effect S` | the birth rate is multiplied by `1 + S` | needed with drivers |
| `--driver-dist exponential\|fixed` | S is a mean or a fixed value | exponential |
| `--clone-time T1:T2:...` | start of each clone, generations after the start (increasing; up to 64 clones) | none |
| `--clone-fitness F1:F2:...` | growth rate of each clone relative to the first clone | needed with `--clone-time` |
| `--clone-parent P1:P2:...` | parent clone of each clone (0 is the first clone, 1 the first of the list) | 0 |
| `-s`, `--sample K` | cells sampled | 20 |
| `-n`, `--replicates R` | number of tumors, each with its tree | 1 |
| `--mutation-rate U` | branch length = generations x U | 1 |
| `--seed X` | seed (written in `log.txt`) | from the clock |
| `-o`, `--output DIR` | folder of the results | `cellforward_out` |

## Results

| File | Content |
|---|---|
| `trees.nwk` | one tree per line (replicate *i* is line *i*), cells `cell0001`, `cell0002`, ...; for `cellcoal -T` |
| `clones.tsv` | replicate, cell, clone, number of drivers, birth rate relative to a cell without drivers |
| `drivers.tsv` | for every clone with cells: replicate, clone, parent clone (-1 for the first), drivers, time of its birth, relative birth rate, cells at the end, cells in the sample |
| `summary.tsv` | per replicate: attempts (tumors that died out or were too small before this one), time, cells, clones, maximum number of drivers, time to the MRCA of the sample (generations before the end) |
| `log.txt` | the command and the seed |

## What is checked

`python3 tests/forward_check.py bin/cellforward-X bin/cellcoal-X` compares neutral tumors of 1000 cells (3000 replicates, samples of 10 cells) with the coalescent of CellCoal under the model of Ohtsuki and Innan (`-K` birth, `-L` death, `-e` the number of cells): the mean time to the MRCA is 64.1 generations in `cellforward` and 66.4 in CellCoal for birth 0.1 and death 0, and 15.3 and 16.0 for birth 0.5 and death 0.1; the coalescent is 3 to 5% deeper because it leaves out the first, stochastic, generations of the tumor. The number of drivers in a pure birth tumor is `2 P (N - 1)` (19.9 against 20.0 in the test). The test suite also checks the outputs, the seed, `--time`, the errors and that CellCoal reads the trees.

## Not in this version

The tumor is a mixed population (no space), the drivers act on the birth rate only, a clone of the user starts from one cell, the mutation rate is the same in all the cells, and the growth is not limited by the size (no competition between clones). These are the next steps.
