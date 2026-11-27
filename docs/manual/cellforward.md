# cellforward: a tumor simulated forward in time

`cellforward` is the companion of CellCoal. CellCoal simulates the genealogy of a sample of cells backward in time (the coalescent); `cellforward` simulates the whole tumor forward in time, one cell division at a time, and writes the tree of a random sample of cells. CellCoal then simulates the mutations, copy number changes, amplification, reads and variant calls on that tree (`-T` or `--tree file=`).

Use the coalescent when the tumor is very large or you need thousands of replicates in seconds, and `cellforward` when the model has things that the coalescent cannot do: the fitness of a clone that comes from the driver mutations that it really got, random sizes and extinctions without approximations, or a check of the coalescent. The cost is proportional to the number of cells of the tumor: about 0.14 s for a tumor of 10<sup>6</sup> cells (the limit is 5·10<sup>7</sup> cells, because every division is kept in memory).

CellCoal and `cellforward` are two separate programs: nothing of one is used in the other. A table that compares them is in the manual of CellCoal (section 5.1.5). The options of the growth model of CellCoal (`-e`, `-g`, `-K`, `-L`, `--clones`) do not exist here, and the options of `cellforward` do not exist in CellCoal; the only things that cross from one to the other are the trees (`trees.nwk`, read by `cellcoal -T`) and the tables of clones.

## The model

- The tumor starts from **one cell** at time 0. Every cell **divides** at the birth rate `--birth` and **dies** at the rate `--death`, both per generation. A tumor that dies out is simulated again, so the tumors that are kept are those that survive.
- The tumor stops when it has **`--size` cells** (every tumor takes a different time) or after **`--time` generations** (every tumor has a different size). Exactly one of the two is given.
- **Drivers.** In every division, each daughter cell gets a new driver mutation with the probability `--driver-rate`. A driver multiplies the **net growth rate** (birth minus death) of the cell, and of all its descendants, by `1 + S` (`--driver-effect`), the same definition as the fitness of the clones of the user and of `--clones` of CellCoal (so the birth rate rises and the death rate does not change; with `--driver-acts-on birth` it is the birth rate that is multiplied by `1 + S`, and then the gain in growth is larger than `S`): `S` is the mean of an exponential variable (`--driver-dist exponential`, the default) or the same for all the drivers (`--driver-dist fixed`). The cell with the driver starts a **clone**, and the descendants of the clone that get more drivers start subclones. Nothing about a clone is an input: its fitness, size and age are a result.
- **Clones started by the user** (`--clone-time`, `--clone-fitness`, `--clone-parent`), as `--clones` of CellCoal but forward: at the time given (generations after the start of the tumor) a random cell of the parent clone gets a growth advantage, and its descendants are the clone. The fitness is the **growth rate relative to the first clone** as in CellCoal: the birth rate of the clone is `death + F (birth - death)` (F = 1.5 grows 50% faster). Only the birth and the time are chosen; the size and the frequency of the clone are a result, and they change among replicates. The clones are numbered 1, 2, ... in the order of the lists (the first clone, without drivers, is 0), and can also get random drivers. A tumor in which the parent clone has no cell when the clone has to start, or that reaches `--size` before the time of the clone, is simulated again (so the tumors kept are those in which the clone exists).
- **Generation, time and division.** In a Wright-Fisher population a generation is one round of reproduction, so the time and the number of divisions of a lineage are the same thing and one mutation rate is enough. In a birth-death process the cells divide at random, with overlapping generations: a cell with birth rate 0.5 per generation divides on average every 2 generations, so the time and the divisions of a lineage are different quantities, and a clone that divides faster makes more divisions in the same time. That is why there are two rates, `U` per generation (mutations that come with time) and `V` per division (errors of the replication of the DNA), and they are added. If you want one rate only, the generation as the division of the Wright-Fisher model is `--birth 1 --death 0` with `--mutation-rate U` (the time unit is then the mean cell cycle), or `--mutation-rate 0 --division-mutation-rate V` for mutations that depend only on the divisions.
- **Hypermutators.** The mutation rate can be different in some cells, independently of the growth: **`--clone-mutation-factor`** gives the mutation rate of every clone of the user relative to normal cells (10 is a hypermutator clone; use `--clone-fitness 1` for a clone that grows as the first one), and **`--mutator-rate`** with **`--mutator-factor`** make a daughter cell a hypermutator at random (the rate of the cell, per generation and per division, is multiplied by the factor, and its descendants inherit it; the factors multiply when it happens again). The mutations per generation of a lineage are `U` times the factor of the cell, integrated along the branch, and the mutations per division are `V` times the factor of the cell that divides. The factor of every sampled cell is in `clones.tsv`.
- **Sample.** `-s` cells are chosen at random from the tumor, without replacement, and their tree is written.
- Times are in **generations**. CellCoal reads the branch lengths of a user tree as the expected number of mutations per site, so the branch lengths of `trees.nwk` are mutations per site: **`--mutation-rate U`** (per site and generation, as `-u` in CellCoal) times the generations of the branch, plus **`--division-mutation-rate V`** (per site and cell division) times the **cell divisions** of the lineage in the branch. With the default (`U` = 1, `V` = 0) the lengths are generations. The mutations linked to the replication of the DNA are the second term: a clone that divides more accumulates more mutations, and the tree in mutations is no longer proportional to the time. Use `--mutation-rate 0` to have the lengths in divisions only.

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
| `--driver-rate P` | probability that a daughter cell gets a driver | none without `--driver-effect`; with it, the rate that gives about 3 or 4 clones with a driver per tumor (printed; it depends on `--size` or `--time` and on the growth) |
| `--driver-effect S` | the birth rate is multiplied by `1 + S` | needed with drivers |
| `--max-clones N` | at most N clones with a random driver per tumor (the clones of `--clone-time` do not count). **It biases the clones towards the early ones**: the N that remain are the first that appear, and later drivers are never assigned, so the natural probabilities change; to keep them, lower `--driver-rate` instead | 0 (no limit) |
| `--driver-dist exponential\|fixed` | S is a mean or a fixed value | exponential |
| `--driver-acts-on growth\|birth\|death` | `growth`: the net growth rate is multiplied by `1 + S` (through the birth rate); `birth`: the birth rate is multiplied by `1 + S`; `death`: the **death rate falls** so that the net growth rate is multiplied by `1 + S` (the birth rate does not change; the death rate cannot go below 0, so the gain is at most birth / (birth - death)); needs `--death` greater than 0. The death rate of every clone is the last column of `drivers.tsv` | growth |
| `--clone-time T1:T2:...` | start of each clone, generations after the start (increasing; up to 64 clones) | none |
| `--clone-fitness F1:F2:...` | growth rate of each clone relative to the first clone | needed with `--clone-time` |
| `--clone-parent P1:P2:...` | parent clone of each clone (0 is the first clone, 1 the first of the list) | 0 |
| `--clone-mutation-factor M1:M2:...` | mutation rate of each clone relative to normal cells | 1 |
| `--mutator-rate P` / `--mutator-factor M` | probability that a daughter becomes a hypermutator / its factor | 0 |
| `-s`, `--sample K` | cells sampled | 20 |
| `-n`, `--replicates R` | number of tumors, each with its tree | 1 |
| `--mutation-rate U` | mutations per site and generation: branch length = generations x U + divisions x V | 1 |
| `--cycle-shape K` | waiting time between divisions: gamma with shape `K` and mean `1 / birth`. `K = 1` is exponential (the cells have no memory of their last division, the default); a large `K` makes the cells divide almost in synchrony, as in a cell cycle, and `K` changes the variance of the number of divisions of a lineage (and so of the mutations linked to the divisions). Deaths stay exponential. A clone of `--clone-time` starts its cycle again | 1 |
| `--division-mutation-rate V` | mutations per site and cell division | 0 |
| `--tree-method auto\|forward\|reconstructed` | `forward`: every cell of the tumor is simulated and the sample is taken at the end. `reconstructed`: only the genealogy of the sampled cells is drawn, directly (a coalescent point process of the birth-death tree, Lambert and Stadler 2013), so the cost does not depend on the size of the tumor (1e9 cells in a few ms per replicate in time mode, tumors of 1e11 cells are possible). Only for neutral tumors that grow (no drivers, clones or hypermutators, `--cycle-shape 1`, birth > death). `auto` uses `reconstructed` whenever it can. With `--time` the number of cells is not simulated (`NA` in `summary.tsv` and `drivers.tsv`). With `--size`, the time to reach the size is drawn exactly up to 2e4 cells and then with a deterministic-growth approximation; with `--size` the divisions that leave no sampled cell are drawn exactly given the tree (the cells that are not sampled add up to N, so they are not independent; extinct side branches are Poisson). Checked against `forward` (mean and standard deviation over replicates of the time and of the divisions of pairs agree within 1-3%, `tests/reconstructed_check.py`) | auto |
| `-j, --jobs J` | replicates run at the same time in `J` processes (0: all the cores). The results do **not** depend on `J`: every replicate has its own random numbers, derived from `--seed` and the number of the replicate | 1 |
| `--seed X` | seed (written in `log.txt`) | from the clock |
| `-o`, `--output DIR` | folder of the results | `cellforward_out` |

## Results

| File | Content |
|---|---|
| `trees.nwk` | one tree per line (replicate *i* is line *i*), cells `cell0001`, `cell0002`, ...; for `cellcoal -T` |
| `clones.tsv` | replicate, cell, clone, number of drivers, birth rate relative to a cell without drivers, mutation factor of the cell |
| `drivers.tsv` | for every clone with cells: replicate, clone, parent clone (-1 for the first), drivers, time of its birth, relative birth rate, cells at the end, cells in the sample |
| `summary.tsv` | per replicate: attempts (tumors that died out or were too small before this one), time, cells, clones, maximum number of drivers, time to the MRCA of the sample (generations before the end) |
| `log.txt` | the command and the seed |

The format of `trees.nwk` and of the tables is the exchange format of the section 5.1.5 of the manual of CellCoal (Newick, one tree per line, branch lengths in expected mutations per site, cells `cell0001`, ...); it is a contract that does not change without a note in the changes of the version.

## What is checked

`python3 tests/forward_check.py bin/cellforward-X bin/cellcoal-X` compares neutral tumors of 1000 cells (3000 replicates, samples of 10 cells) with the coalescent of CellCoal under the model of Ohtsuki and Innan (`-K` birth, `-L` death, `-e` the number of cells): the mean time to the MRCA is 64.1 generations in `cellforward` and 66.4 in CellCoal for birth 0.1 and death 0, and 15.3 and 16.0 for birth 0.5 and death 0.1; the coalescent is 3 to 5% deeper because it leaves out the first, stochastic, generations of the tumor. The number of drivers in a pure birth tumor is `2 P (N - 1)` (19.9 against 20.0 in the test). The test suite also checks the outputs, the seed, `--time`, the errors and that CellCoal reads the trees.

## Not in this version

The tumor is a mixed population (no space), the drivers act on the birth rate only, a clone of the user starts from one cell, hypermutators are inherited but a clone has no other mutation process (no signatures), and the growth is not limited by the size (no competition between clones). These are the next steps.
