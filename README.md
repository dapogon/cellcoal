# CellCoal
Coalescent simulation of somatic cell samples  
(c) 2018, 2019 David Posada, University of Vigo, Spain <http://darwin.uvigo.es>

## About CellCoal
CellCoal simulates the somatic evolution of single-cells. CellCoal generates a coalescent tree and can simulate a sample of diploid genomes from somatic cells –no recombination– from a growing population, together with a another cell as outgroup. For example a growing tumor population together with a healthy cell as outgroup. CellCoal implements multiple mutations models (0/1, DNA, infinite and finite site models, deletion, copy-neutral LOH, 30 cancer signatures) and is able to generate read counts and genotype likelihoods considering allelic dropout, sequencing and amplification error, plus doublet cells.

## Getting started

1. Download: get the program from the GitHub repository at <https://github.com/dapogon/cellcoal>. Under the "Code" tab you will see a section call "release". Click on it download the source code in .zip (or tar.gz) format. Then unzip the *cellcoal-x.y.z.zip* file (for example, cellcoal-1.0.0). You should see now a folder called cellcoal-x.y.z. Move into the folder to compile and run the program.

2. Compile: type make. The program should compile without much problem in Linux/MacOSX. The makefile is provided for the gcc compiler.

3. Run: type ./cellcoal. Arguments will be read from a parameter file called `parameters`. Alternatively, arguments can be entered directly in the command line.


## Citation

- If you use CellCoal, please cite: Posada D. (2019) *CellCoal: coalescent simulation of single-cell genomes*. <https://github.com/dapogon/cellcoal>

## Manual
A detailed manual is available at <https://dapogon.github.io/cellcoal/>. A PDF version of the manual can also be obtained <a href="https://dapogon.github.io/cellcoal/cellcoal.manual.v1.pdf" target="_blank">here</a>.
