## CellCoal 2
CellCoal: coalescent simulation of single-cell NGS genotypes 

&copy;2020 David Posada, University of Vigo, Spain <http://darwin.uvigo.es>

## About CellCoal
CellCoal simulates the somatic evolution of single-cells. CellCoal generates a coalescent genealogy for a sample of somatic cells –no recombination– obtained from a growing population, together with a another cell as outgroup, introduces mutations along this genealogy, and produces single-cell diploid genotypes (single-nuclotide variants or SNVs). CellCoal implements multiple mutations models (0/1, DNA, infinite and finite site models, deletion, copy-neutral LOH, 30 cancer signatures) and is able to generate read counts and genotype likelihoods considering allelic dropout, sequencing and amplification error, plus doublet cells.

## Getting started

1. **Download**: get the program from the GitHub repository at <https://github.com/dapogon/cellcoal>. Under the Code tab you will see a section call *release*. Click on it download the source code in .zip (or tar.gz) format. Then unzip the *cellcoal-x.y.z.zip* file (for example, cellcoal-2.0.0.zip). You should see now a folder called *cellcoal-x.y.z*. Move into the folder to compile and run the program.

2. **Compile**: type `make`. The program should compile without much problem in Linux/MacOSX. A *Makefile* is provided for the gcc compiler. The executable file will be located in the *bin* folder.

3. **Run**: `bin/cellcoal-x.y.z`. In this case, the program arguments will be read from a file called *parameters*. Alternatively, arguments can be entered directly in the command line, or from a specific parameter file with a different name.

4. **Inspect results**: Once the run is finished, different output files will be written by default inside a folder called *results*.

## What's new in 2.0.0

Bug fixes (genotyping error, allelic dropout variation, memory safety), 3–4× faster with about 7× less memory, and parallel read-count generation (`OMP_NUM_THREADS`). Output for a given seed differs from versions ≤ 1.3.4. See [docs/CHANGES-2.0.0.md](docs/CHANGES-2.0.0.md).

## Citation

- If you use CellCoal, please cite: Posada D. (2020) CellCoal: coalescent simulation of single-cell sequencing samples. *Molecular Biology and Evolution* 37: 1535-1542. <a href="https://academic.oup.com/mbe/article/37/5/1535/5728646" target="_blank">https://academic.oup.com/mbe/article/37/5/1535/5728646</a>.

## Manual
A detailed manual is available in [html](https://dapogon.github.io/cellcoal/) or in [pdf](docs/cellcoal.manual.v2.0.pdf).

