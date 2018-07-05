# CoalTumor
Coalescent simulation of somatic cells  
(c) 2018 David Posada, University of Vigo, Spain <http://darwin.uvigo.es>

## About CoalTumor
CoalTumor simulates the somatic evolution of single-cells. CoalTumor generates a coalescent tree and can simulate a sample of diploid genomes from tumoral cells –no recombination– from a growing population, together with a healthy cell as outgroup. CoalTumor implements multiple mutations models (0/1, DNA, infinite and finite site models, cancer signatures) and is able to generate read counts and genotype likelihoods considering allelic dropout, sequencing and amplification error.

## Getting started

1. Download: get the program from the GitHub repository at <https://github.com/dapogon>. An easy way of doing this is to press the 'Download ZIP' button located under the "Clone or Download" drop down green button. Then unzip the *coaltumor-master.zip* file. You should see now a folder called coaltumor-master. Move into the folder to compile and run the program.

2. Compile: type make. The program should compile without much problem in Linux/MacOSX. The makefile is provided for the gcc compiler.

3. Run: type ./coaltumor. Arguments will be read from a parameter file calle. Alternatively, arguments can be entered directly in the command line.


## Citation

- If you use CoalTumor, please cite: Posada D. (2018) *CoalTumor: coalescent simulation of single-cell genomes*. <https://github.com/dapogon/coaltumor>

