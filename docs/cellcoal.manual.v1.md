---
title: "CellCoal.manual.v1"
author: David Posada
date: June, 2019
output:
  html_document:
    fig_caption: yes
    toc: yes
    toc_float:
      collapsed: no
      smooth_scroll: yes
  word_document:
    toc: yes
    toc_depth: '3'
  pdf_document:
    toc: yes
    toc_depth: '3'
---

<style>
p.caption {
  font-size: 0.9em;
   color: gray;
  margin-top: 2%;  
  margin-bottom: 5%;  
  text-align: justify;
}

.small {
   font-size: 0.9em;
}
</style>


#1. About CellCoal
CellCoal simulates the somatic evolution of single-cells including genotypes resulting from single-cell sequencing. CellCoal generates a coalescent tree and can simulate a sample of diploid genomes from somatic cells –no recombination– from a growing population, together with another genome as outgroup. A typical example would be a sample of single-cell genotypes from a growing tumor population together with a healthy cell used as outgroup. CellCoal implements multiple mutation models for single-nucleotide variants (0/1, DNA, infinite and finite site models, point deletions, copy-neutral LOH, cancer mutational signatures) and is able to generate a .VCF file with NGS read counts and genotype likelihoods considering technical artifacts like allelic imbalance, allelic dropout, sequencing error, amplification error, plus doublet cells.


#2. Citation
If you use CellCoal, please cite it as:

Posada D. 2019. CellCoal: coalescent simulation of single-cell genomes.
<https://github.com/dapogon/cellcoal>.


#3. Getting started
1. **Download**: get the program from the GitHub repository at <https://github.com/dapogon/cellcoal>. Under the Code tab you will see a section call *releas*e. Click on it download the source code in .zip (or tar.gz) format. Then unzip the *cellcoal-x.y.z.zip* file (for example, cellcoal-1.0.0.zip). You should see now a folder called *cellcoal-x.y.z*.Move into the folder to compile and run the program.

2. **Compile**: type `make`. The program should compile without much problem in Linux/MacOSX. A Makefile is provided for the gcc compiler. The executable file will be located in the bin folder.

3. **Run**: type `./cellcoal`. In this case, the program arguments will be read from a file called *parameters*. Alternatively, arguments can be entered directly in the command line.

4. **Inspect results**: Once the run is finished, different output files will be written inside a folder called *results*.


#4. Basic usage
CellCoal does not have a Graphical User Interface (GUI). It works on the Linux/Mac command line in a non-interactive fashion. It can run in two main ways parsing its arguments from a file called *parameters*, or parsing its arguments directly from the command line.

If no arguments are passed in the command line, CellCoal will look for a file called *parameters* in the current working directory:

`./cellcoal`

The *parameters* file is a text file that contains the different arguments for the program. The arguments correspond to different letters whose case (uppercase or lowercase) makes a difference. Comments within brackets are allowed. To avoid unexpected problems, users are encouraged to maintain the same order of the arguments as in the example files. For clarity, the example parameters files are divided in different blocks or sections including related arguments.

Otherwise, CellCoal will run under the arguments specified in the command line, using default values for the absent arguments:

`./ cellcoal -n100 -s20 -l1000 -e10000 -g1.0e-04 -k1 -i1 -b1 -j250 -p0.0 -f0.3 0.2 0.2 0.3 -r0.00 0.03 0.12 0.04 0.11 0.00 0.02 0.68 0.68 0.02 0.00 0.11 0.04 0.12 0.03 0.00 -1 -2 -3 -4 -6 -9 -v -x -#200011`

A brief usage guide plus the current default values for the different parameters can be obtained with the argument -h:

`./cellcoal -h`


# 5. CellCoal implementation
The main steps of CellCoal are:

1. Coalescent simulation of the genealogy of a sample of cells (the user can also specify a given tree).

2. Addition of an outgroup and optionally, modification of branches.

3. Simulation of a reference and an ancestral genome, with or without germline variation.

4. Addition of mutations, LOH and deletions along the tree in order to obtain the true genotypes.

5. Addition of technical errors without generating read counts (ADO + genotype errors)
OR Addition of technical errors while generating read counts (ADO + amplification error + sequencing error.

6. Print genotypes/haplotypes to files.

7. Print summary statistics.


![Fig 1. Main flow of CellCoal.](figures/fig1.png){width=100%}


## 5.1 Coalescent genealogy
CellCoal simulates a cell genealogy under the standard neutral coalescent, going backward in time (Fig. 2 left). Note that the coalescent assumes the cell sample comes from a much larger population, which in this case can be constant or growing. In the coalescent we use the *effective population size*, a well-known quantity in population genetics that refers to the number of individuals in an ideal population. CellCoal can consider continuous exponential population growth (Fig. 2 middle) and multiple demographic periods (Fig. 2 right). For the latter, the exponential growth rate during the period (positive or negative) will be deduced from the specified N at the beginning and at the end of the period. The growth rate derived for the last period will continue into the indefinite past. This implementation is borrowed from Hudson (2002).

>If the growth rate is negative, the coalescent time may become
>infinite (i.e., coalescence does not happen). In this case, CellCoal
>will stop and issue an error message.

![Fig. 2. A random coalescent genealogy for a sample of 4 cells under constant (left), and growing population size (middle) and demographic periods (right). The effective population size is indicated as blue background. In the right, from Period 0 to Period 1 the population size experiments a severe bottleneck.](figures/fig2.png){width=100%}



### 5.1.1 Adding an outgroup
After the sample genealogy is simulated, two additional branches are added: a **root branch** and an **outgroup branch** (Fig. 3), in order to include an outgroup cell in the sample. In phylogenetics an outgroup is a distant organism that serves as a reference and is often used for rooting purposes. The length of these two branches is controlled by the user (Fig. 3).

![Fig. 3. Outgroup addition. The length of the root branch will be controlled by the user with the parameter -k, which specifies the length of this branch as a proportion of the sample MRCA mean depth. The length of the outgroup branch depends on the parameter -q, which is the length of this branch as a proportion of the sample/outgroup most recent common ancestor (MRCA) mean depth. In this example, the corresponding values were -p0.5 -k1.0.](figures/fig3.png){width=60%}

In a tumor scenario, the outgroup cell will be naturally a healthy cell (Fig. 4)

![Fig. 4. Addition of a healthy cell in a tumor scenario. In this case we expect that the healthy cell has evolved much slower than the tumoral cells, so a small (or zero) outgroup might make more sense. In this example, the corresponding values where -p1 -k0.05.](figures/fig4.png){width=70%}


### 5.1.2 Rate variation among lineages
The standard coalescent results in an ultrametric tree, where all the tips are at the same distance from the most recent common ancestor (MRCA) of the sample. These branches can be modified by introducing rate variation among linages/branches, here using multipliers sampled from a Gamma distribution (see Fig. 5). This approach was introduced by Yang (for an easy review see Yang 1996) to add substitution rate variation among different nucleotide sites in a DNA sequence. Biologically this is interesting if we want to simulate for example a situation under which some lineages evolve at different rates (i.e., there is not a molecular clock) as if they were under different selective regimes.

![Fig. 5. Rate variation among branches. Rate variation among lineages can be a product for example of selection.](figures/fig5.png){width=60%}


## 5.2 Genotype evolution models
CellCoal simulates the evolution of cell genotypes along the genealogy by adding somatic mutations (single-nucleotide variants -SNVs-, copy-neutral loss of heterozygosity events -cnLOH-, and point deletions), starting from an ancestral genome at the *sample/outgroup MRCA*. Note that CellCoal does not simulate copy number alterations or structural variants. CellCoal considers two possible **alphabets (0/1 or DNA; option -b)** and several mutation models (infinite and finite sites models, including trinucleotide mutational signatures), which can be mixed along the genome. 

>The default SNV model is a diploid infinite site model (ISM).

### 5.2.1 Generation of reference and ancestral genomes
The first step to generate the cell genotypes is the generation of a **reference genome**. Under the binary alphabet, it will consist of all 0s. For DNA, the reference genome will be simulated at random according to the specified nucleotide frequencies, unless genetic signatures are specified, in which case the reference genome will be simulated according to the trinucleotide frequencies in the human genome (version GRCh37). Alternatively, for DNA data, the user can specify its own reference genome in a FASTA file.

The **ancestral genome** at the sample/outgroup MRCA of the genealogy will be simulated by introducing **germline variation (option -c)** (single nucleotide polymorphisms or SNPs) into the reference genome at a rate specified by the user. Logically, if this rate is 0, the reference and the ancestral genomes will be identical. Once the ancestral genome is in place, mutations are introduced according to the specified mutation model.

### 5.2.2 Infinite-site SNV models
CellCoal implements several types of **single-nucleotide mutation models (option -m)**. Infinite-sites models (ISM) allow a given site to change only once. However, in a diploid scenario we can implement this model in at least two ways Under the **diploid ISM**, a given site can change only once in either the maternal and paternal genomes. That is, the change will be for example either 0/0 -> 1/0 or 0/0 -> 0/1 but not both (where 0 is the ancestral allele and 1 the derived one). In consequence, under the diploid ISM, homozygous mutants (1/1) are not possible. A different implementation of the ISM, called here **haploid ISM**, considers instead maternal and paternal sites as different locations, so it allows for two somatic mutations to take place at the same site in both maternal and paternal genomes. Thus, under this model 1/1 genotypes are possible. In both ISMs, the number of somatic mutations is distributed as a Poisson according to the total tree length (which is a function of the mutation rate and time), unless a **fixed number of mutations (option -j)** is requested. For the haploid ISM the expected number of mutations can be easily calculated:

$$
E\left( \text{number of mutations} \right) = \Theta \sum_{i=1}^{i=n-1} \frac{1}{i}
$$


where $\Theta$ = 4*N*μ2*L*, being *N* the **effective population size (option –e)**, μ the **mutation rate per site per generation (option -u)**, and *L* the **number of sites (option -l).** Alternatively, it is also possible to specify a **fixed number of mutations** **(option -j)** to occur along the genealogy.

The specific sites where the ISM mutations occur are selected uniformly along the genome, unless mutational signatures are requested (see below). The particular branches of the cell genealogy where mutations are located are selected according to their lengths. That is, longer branches will receive more mutations on average.

### 5.2.3 Infinite-site mutational signatures
For DNA data, mutational signatures (sensu Alexandrov *et al.* 2013) can be simulated under a **diploid infinite trinucleotide model (ISTM) (option -S)**, that considers the nucleotide context of the mutated site (bases imimagetely 5’ and 3’). Currently, 30 distinct signatures can be incorporated (<https://cancer.sanger.ac.uk/cosmic/signatures>) and mixed in different proportions. As in ISMs, the number of somatic mutations is distributed as a Poisson according to the total tree length

Importantly, in the current implementation we assume that the context is constant along the tree. That is, we select a site based on the context in the ancestral genome, before locating it in a given tree branch. Thus, this assumption will be broken in the case of contiguous mutations.

>In order to get the expected counts under each signature it is important to simulate long genomes,</span> otherwise rare trinucleotide will not be available and, depending on the signature, could bias the observed counts.

### 5.2.4 Infinite-site copy-neutral LOH
**Copy-neutral loss of heterozygosity (cnLOH) (option -H**) events (e.g., A/G -\> A/A, or A/G -\> G/G) can be added assuming a haploid ISM, which means that cnLOH cannot happen in the same site twice unless they occur in a different maternal/paternal genome. As for SNVs, the number of cnLOH events will be distributed as a Poisson according to the total cnLOH tree length (which in this case is a function of the cnLOH rate and time). The specific sites where the deletions occur are selected uniformly along the genome. The particular branches of the cell genealogy where deletions are located are selected according to their
(cnLOH) lengths.

### 5.2.5 Infinite-site single-nucleotide deletions
**Single-nucleotide deletions (option -d)** (e.g., N/N -\> –/N, or N/N -\> N/–) can be added also assuming a haploid ISM. As for SNVs, the number of point deletions will be distributed as a Poisson according to the total deletion tree length (which in this case is a function of the deletion rate and time). The specific sites where the deletions occur are selected uniformly along the genome. The particular branches of the cell genealogy where deletions are located are selected according to their (deletion) lengths.

### 5.2.6 Finite-site SNV models
CellCoal also implement several **finite-sites mutation models (option -m)** (FSM) for SNVs. In this case, multiple mutations can happen along a branch. Typically, these models assume a continuous Markov process and are very typical in phylogenetics, where they are often referred as substitution models (<https://en.wikipedia.org/wiki/Models_of_DNA_evolution>).

For 0/1 data, the model implemented is known as Cavender-Farris-Neyman or Mk2 model (see Lewis 2001), and is equivalent to a Jukes-Cantor (1969) model for two alleles. For DNA data, reversible (assume that the rate X-\>Y is the same as Y-\<X) and non-reversible models, can be built by specifying the **nucleotide frequencies (option -f)** (at equilibrium; there are not expected to change along the genealogy a **transition/transversion rate ratio (option -t)** or a full 4x4 **instantaneous mutation rate matrix (option -r)**, plus **rate variation among sites** **(option -a)** (using a Gamma distribution; see Yang 1996). One can define in this way popular models like Jukes-Cantor (1969), Hasegawa-Kishino-Yano (1985) or General-Time-Reversible (Tavaré 1986). It is also possible to specify non-reversible models (i.e., asymmetrical rates among nucleotides).

### 5.2.7 Mixed ISM-FSM
Both ISM and FSM for SNVs can be incorporated by setting up a **proportion of sites expected to evolve under each type of model (option –p).** Each site is then assigned to each mutation category, ISM/FSM, with probability p and 1-p, respectively. Also, **relative mutation rates (option –w)** can also be specified for ISM and FSM sites

## 5.3 Single-cell genomics
CellCoal can simulate the technical artifacts resulting from single-cell genomics, like allelic dropout (ADO), allelic imbalance, sequencing error, amplification error, plus doublet cells. Here there are two main strategies:

1.  **No sequencing reads are generated**; instead ADO
    and genotyping errors are directly imposed on the evolved genotypes
    at the tip of the genealogy.

2.  **Read counts are simulated** at a variable depth
    along the genome, including allelic imbalance, ADO, amplification
    and sequencing errors. Then, genotype calling is performed under
    different maximum likelihood models for single-cells.

### 5.3.1 Allelic dropout (ADO)
During the amplification of single-cells, one or both alleles at a single site can be lost, so a truly heterozygous genotype might appear as a homozygote (e.g., 0/1 => 0/_ or 0/1 => _/1). CellCoal can remove alleles from single chromosomes at a given **ADO rate per genotype (option -D)**. This ADO rate can vary both **across sites (option -P)** and **across cells (option -Q**) according to Gamma distributions.

### 5.3.2 Genotyping error
**Genotype errors (option -G)** can be directly introduced in the observed genotypes at a given rate per genotype. For DNA data,nucleotides are substituted by a different one according to a user-defined error matrix (option -X). These errors encapsulate into a single class different types of error than can take place along the single-cell genomics pipeline, like amplification, sequencing and calling errors.

### 5.3.3 Sequencing coverage
For each individual genotype the program will generate read counts at a certain sequencing coverage or depth. For a given site the number of reads will follow a Poisson distribution around the **mean coverage (option -C)**, unless some **coverage overdispersion (option -V)** is specified, in which case the coverage per site follows a *Negative Binomial* distribution. Reads will randomly assigned to the maternal or paternal genome, according to a given **allelic imbalance (option -I)** that by default is 0.5 (i.e., the same probability for maternal or paternal read). If the genotype is haploid due to ADO, the expected **haploid coverage reduction (option -R)**, which is 0.5 by default, can also be specified.

### 5.3.4 Amplification error
In order to obtain enough DNA material from single-cells for sequencing library construction, the genomes of single cells are often amplified. During this amplification step the polymerase can introduce nucleotide errors that can appear in a high proportion of the resulting templates for that site if occurring during the first amplification cycles. If many templates contain an error, many reads will be also be wrong for the specific sites. The **amplification error (option -A)** for a given maternal/paternal site follows a Beta distribution with mean and variance specified by the user. The probability of the different types of errors (e.g., A=\>C, A=\>G, A=\>T ,etc.) can be indicated using an **4x4 error matrix (option – X)**. CellCoal considers two different amplification error models:

**Four-template model**: In this case all four bases can be present, due to multiple amplification errors, in the set of amplified templates (Fig. 6).

![Fig 6. Four-template model for amplification (γ) and sequencing (ε) error.](figures/fig6.png){width=40%}


**Two-template model**: In this case only the original base and a single alternative can be present, due to a single amplification error, in the set of amplified templates (Fig. 7).

![Fig 7. Two-template model for amplification (γ) and sequencing (ε) error.](figures/fig7.png){width=40%}

### 5.3.5 Sequencing error
Once the genomes have been amplified, during sequencing, additional errors can add to the reads. **Sequencing error (option -E)** is assumed to be constant across sites and cells. The probability of the different types of errors (e.g., A=\>C, A=\>G, A=\>T, etc.) uses the same 4x4 error matrix **(option – X)** as the amplification error.

### 5.3.6 Doublets 
Sometimes, due to errors during the cell isolation process, a single sequencing library can contain DNA templates from two cells, and therefore reads which are assume to come from a single cell come in fact from two different cells. This is called a “doublet”. CellCoal can consider the effect of such doublets by producing reads from two cells and joining these reads as if they originated from a single cell, given a user-specified **doublet rate per cell (option -B)**.

### 5.3.7 Genotype likelihoods
CellCoal can simulate log10 genotype likelihoods (GLs) for single-cell NGS data, given the read counts simulated, and the sequencing and amplification errors specified by the user, under three different models: GATK standard, 4-template and 2-template models.

#### 5.3.7.1 GATK standard model
The basic model used to calculate the GLs is similar to that implemented in GATK (DePristo et al. 2011; see also Korneliussen et al. 2013). It assumes the same sequencing error rate (ε) for all bases and sites, and no amplification error:

$$
\Pr \left( D|G=\{A_{1},A_{2}\} \right) = \prod_{i=1}^{R}{\Pr\left( r_{i}|G = \{ A_{1},A_{2}\} \right)} = \prod_{i = 1}^{R}\left( \frac{1}{2}\ \text{p}\left( r_{i} | A_{1} \right) + \frac{1}{2}\ \text{p}\left( r_{i}| A_{2} \right) \right)
$$

where
$$
p\left(r|A\right) = \left\{\begin{matrix} \varepsilon/3, r\neq A \\ 1 - \varepsilon, r = A \\ \end{matrix} \right.  
$$

where

<small>
*G* = genotype

*A* = allele

*R* = number of reads

*r<sub>i</sub>* = observed base in read *i*

\(\varepsilon\) = probability of sequencing error


#### 5.3.7.2 Four-template amplification error model
This model extends the previous one in order to also consider the amplification error, which (together with the sequencing error) can be different for distinct bases **(option – X**) and that is sampled for each site from a *Beta* distribution with mean and variance specified by the user **(option -A)**. It allows for multiple amplification errors at a single site, and therefore that all four templates are possible (the correct one and the other three; see Fig. 6 above).

$$
\Pr\left( D \middle| G = \{ A_{1},A_{2}\} \right) = \prod_{i = 1}^{R}{\Pr\left( r_{i} \middle| G = \{ A_{1},A_{2}\} \right)} = \prod_{i = 1}^{R}\left( \frac{1}{2}\ \text{p}\left( r_{i} \middle| A_{1} \right) + \frac{1}{2}\ \text{p}\left( r_{i} \middle| A_{2} \right) \right)
$$

where

$$
p\left(r|A\right) = \sum_{j = 1}^{4}{p\left( T_{j} \middle| A \right)\ p\left( r|T_{j} \right)} = \sum_{j = 1}^{4}\gamma_{A \rightarrow T_{j}}\varepsilon_{T_{j} \rightarrow r}
$$

and
$$
\gamma_{i \rightarrow j} = \left\{ \begin{matrix} 1 - \gamma, i = j \\
\gamma e_{i \rightarrow j}, i \neq j \\
\end{matrix} \right.\ \ \  ,\ \ \  \varepsilon_{i \rightarrow j} = \left\{ \begin{matrix}
1 - \varepsilon, i = j \\
\varepsilon e_{i \rightarrow j}, i \neq j \\
\end{matrix} \right.\
$$

where

<small>
*G* = genotype

*A* = allele

*T* = amplified allele template

*R* = number of reads

*r<sub>i</sub>* = observed base in read *i*

\(\varepsilon\) = probability of sequencing error for a given site

\(\gamma\) = probability of amplification error for a given site

\(\gamma_{i \rightarrow j}\) = probability of amplification of base *i*
into template base *j*

\(\varepsilon_{i \rightarrow j}\) = probability of sequencing error from
template base *i* to read base *j*

\(e_{i \rightarrow j}\) = relative probability of
amplification/sequencing error from base *i* to base *j*
</small>
<br><br>

Note that if the amplification and sequencing errors are constant, the probability of observing a read *r* coming from an allele *A* simplifies to:

$$
p\left( r \middle| A \right) = \left\{ \begin{matrix}
\left( 1 - \gamma \right)\ \varepsilon/3\  + \ \gamma\ \frac{1 - \frac{\varepsilon}{3}}{3},r \neq A \\
\left( 1 - \gamma \right)\left( 1 - \varepsilon \right)\  + \ \gamma\ \varepsilon/3,r = A \\
\end{matrix} \right.\
$$

<br><br>

#### 5.3.7.3 Two-template amplification error model
This model is very similar to the previous one, but it assumes that only a single amplification error can occur at a single site, and that therefore only two templates (the “correct” one and a wrong one) are possible (see Fig. 7 above).

$$
\Pr\left( D \middle| G = \left\{ A_{1},A_{2} \right\},T = \left\{ T_{1},T_{2} \right\} \right) = \prod_{i = 1}^{R}{\Pr\left( r_{i} \middle| G = \left\{ A_{1},A_{2} \right\},T = \left\{ T_{1},T_{2} \right\} \right)} = \prod_{i = 1}^{R}\left( \frac{1}{2}\ \operatorname{p}\left( r_{i} \middle| A_{1},T_{1} \right) + \frac{1}{2}\ \operatorname{p}\left( r_{i} \middle| A_{2},T_{2} \right) \right)
$$

where 
$$
p\left( r \middle| A,T \right) = \sum_{\begin{matrix} T = 1 ,\ T \neq A \\ \end{matrix}}^{4}e_{A \rightarrow T}\lbrack\left( 1 - \gamma \right)\varepsilon_{A \rightarrow r} + \gamma\ \varepsilon_{T \rightarrow r}\rbrack
$$

and
$$
\varepsilon_{i \rightarrow j} = \left\{ \begin{matrix} 1 - \varepsilon,i = j \\ \varepsilon e_{i \rightarrow j},i \neq j \\ \end{matrix} \right.\
$$

where

<small>
*G* = genotype

*A* = allele

*T* = amplified allele template

*R* = number of reads

*r<sub>i</sub>* = observed base in read *i*

\(\varepsilon\) = probability of sequencing error for a given site

\(\gamma\) = probability of amplification error for a given site

\(\gamma_{i \rightarrow j}\) = probability of amplification of base *i*
into template base *j*

\(\varepsilon_{i \rightarrow j}\) = probability of sequencing error from
template base *i* to read base *j*

\(e_{i \rightarrow j}\) = relative probability of
amplification/sequencing error from base *i* to base *j*
</small>
<br><br>


# 6\. Arguments

## 6.1 Coalescent arguments
In this section we describe the arguments that control the simulation of the cell genealogy under the neutral coalescent. Here we can indicate the number of replicates, the number of cells in the sample, the number of sites or loci, and characteristics of the population from which the cells have been sampled, like effective population size, exponential growth or demographic periods.

```
[COALESCENT]
[number of replicates] n100
[number of cells] s20
[number of sites] l1000
[effective population size] e10000
[exponential growth rate] g1e-4 [ln2=0.6931 max]
[demographics] [h2 1000 100 40000 200 30000 20000]
```

**n#: number of replicates; integer \[1-inf)**

The number of samples to be generated. Each sample is an independent realization of the whole process (coalescent, mutations, NGS; see below).

**s#: number of cells; integer \[2-inf)**

The number of cells to be generated for each sample.

**l#: number of sites; integer \[1-inf)**

The total length, in base pairs or nucleotides, of the genomic region simulated. This length includes variable and non-variable sites.

**e#: effective population size; integer \[1-inf)**

The effective size of the population from which the sample was theoretically drawn.

**h# # # # demographic periods; several integers \[1-inf)** 

The number of demographic periods, from present to past, and the effective population size during those periods. The first number specifies the number of periods. Then for each period there should be three consecutive numbers indicating the effective population size at the beginning and at the end of the period, and the duration of the period in generations. This option is incompatible with the *exponential
growth rate option (-g)* below.

>These parameters are looking back in time, so it is not a good idea to imply a negative growth rate for the last period, as the coalescent time could become infinite in the past.

**g\#: exponential growth rate; double (-inf, +inf)**

The rate of exponential growth per individual per generation. This option is incompatible with the 
**demographic periods (-d) option**.

>When simulating cells, that divide in a binary fashion, the maximum exponential growth rate for the population should be ln2 (0.6931). Also, this parameter looks back in time, so it is not a good idea to specify a negative growth rate, as the coalescent time could become infinite in the past.


## 6.2 Genealogy modifiers
```
[GENEALOGY MODIFIERS]
[root branch length ratio] k0.5
[outgroup branch length ratio] q1.0
[rate variation among branches] i1.2
```
**k#: root branch length ratio; double \[0-inf)**

Length of the root branch as a proportion of the mean depth of the *sample MRCA*. Note that this branch is attached before introducing rate variation among branches. See Figs. 1-2.

**q#: outgroup branch length ratio; double \[0-inf)**

Length of the outgroup branch as a proportion of the mean depth of the *sample/outgroup MRCA*. Note that this branch is attached before introducing rate variation among branches. See Figs. 1-2.

**i#: rate variation among branches; double (0-inf)**

Amount of rate variation among lineages, specified with the alpha shape parameter of a continuous gamma distribution with mean 1. The smaller the shape the strongest the heterogeneity introduced. See Fig. 4.


## 6.3 Mutation models
```
[MUTATION MODEL]
[alphabet binary:0 DNA:1] b1
[germline SNP rate] [c1e-5]
[mutation rate] u1e-6
[deletion rate] d1e-6
[CN\_LOH rate] [H1e-5]
[fixed number of mutations - ISM] [j1000]
[trinucleotide genetic signature - ISM] S2 11 0.3 22 0.7
[alternative mutation model ISMhap:0 Mk2:1 finiteDNA:2] [m0]
[proportion of alternative model sites] [p0]
[alternative/default model relative mutation rate] [w1]
[base frequencies] [f0.3 0.2 0.2 0.3]
[transition/transversion ratio] [t1.7]
[rate variation among sites] [a1.0]
[mutation matrix ACGT x ACGT] [r0.00 0.03 0.12 0.04
0.11 0.00 0.02 0.68
0.68 0.02 0.00 0.11
0.04 0.12 0.03 0.00]
```

**b#: alphabet: integer [0/1]**

Alphabet for the simulated characters, which can be binary (0,1; option b0) or DNA nucleotides (A,C,G,T; option b1)


**c#: germline mutation rate: double [0-inf)**

Germline mutation rate per site per generation, which is the same as the probability of introducing a mutation in the ancestral genome (i.e., a SNP in comparison with the reference genome).


**u#: somatic mutation rate: double [0-inf)**

Nucleotide mutation rate per site per generation.


**d#: somatic deletion rate: double [0-inf)**

Deletion rate per site per generation under a haploid ISM.

**H#: copy-neutral loss-of-heterozygosity rate: double [0-inf)**

Copy-neutral LOH rate per site per generation under a haploid ISM.

**j#: fixed number of mutations: integer [0-inf)**

Number of ISM mutations to be added along the genealogy.

**S# # #: trinucleotide genetic signature/s: integer [1-inf), integer [1-30], double [0-1]**

Genetic signatures to be implemented and their proportions. The first
number should be the number of distinct signatures, followed by the
signature/s IDs and the expected proportion of mutations to be drawn
from each particular signature. If these proportions do not add to one
they will be rescaled. For example:

S1 5 1.0 implies that mutations will come from signature 5 only.

S2 11 0.3 22 0.7 implies that mutations will come from two signatures:
from signature 11 with probability 0.3, and from signature 22 with
probability 0.7.

Note that to get the expected counts under each signature it is
important to simulate long genomes, otherwise rare trinucleotide will
not be available and, depending on the signature, could bias the
observed counts.

**m#: alternative mutation model: integer [0-2]**

Alternative mutation model (the default model is ISM diploid). The
options are as follows:
```
0 = ISM haploid
1 = MK2
2 = finite DNA
```

**p#: proportion of alternative model sites: double [0-1]**

Proportion of sites that will be simulated under the alternative model
specified with the -m option.


**w#: alternative/default model relative mutation rate: double [0-inf)**

Ratio of the mutation rate for sites evolving under the alternative
mutation model, relative to the mutation rate for the sites evolving
under the null mutation (option -u).

**f# # # #: base frequencies: double [0-inf)**

Nucleotide ACGT frequencies at equilibrium for the DNA finite-site
models. If they do not add to one they will be rescaled.

**t#: transition/transversion ratio: double (0-inf)**

Ratio of the number of *transitions* to the number of *transversions* for DNA finite-site models. This ratio, often referred as ti/tv, is 0.5 when *transitions*  and *transversions* are equally probable, as there are twice as many possible *transversions* as *transitions* .

**a\#: rate variation among sites: double (0-inf)**

Amount of mutation rate variation among sites, specified with the *alpha shape parameter* of a continuous gamma distribution with mean 1. The smaller the shape the strongest the heterogeneity introduced.

**r#4×4: mutation matrix: double [0-inf)**

This 4×4 matrix specifies the relative mutation rates (instantaneous substitution rates in the case of finite-site models; often called in this case as Q-matrix) among the different nucleotides, in order ACGT for rows and columns. It does not need to be symmetrical, but all diagonals should be zero. If this option is active, the transition/transversion ratio (-t) is ignored. For example, for the Jukes-Cantor (1969) model this matrix would be:

$$\mathbf{Q} = \left[\begin{array}
{rrrr}
1 & 1 & 1 & 1 \\
1 & 1 & 1 & 1 \\
1 & 1 & 1 & 1 \\
1 & 1 & 1 & 1 \\
\end{array}\right]
$$

For the General-Time-Reversible (GTR) (Tavaré 1986) model:
$$\mathbf{Q} = \left[\begin{array}
{rrrr}
0 & 2.1 & 0.3 & 5.6 \\
2.1 & 0 & 1.5 & 7.8 \\
0.3 & 1.5 & 0 & 1.4 \\
5.6 & 7.8 & 1.4 & 0 \\
\end{array}\right]
$$

And for a General-Time-Non-Reversible (GTnR) model:
$$\mathbf{Q} = \left[\begin{array}
{rrrr}
0 & 2.1 & 0.3 & 5.6 \\
1.9 & 0 & 1.5 & 7.8 \\
2.1 & 3.7 & 0 & 1.4 \\
0.8 & 2.4 & 0.9 & 0 \\
\end{array}\right]
$$


## 6.4 NGS parameters

```
[NGS PARAMETERS]
[genotyping error; no reads] [G0.1]
[sequencing coverage; reads simulated] C10
[coverage overdispersion] [V5]
[mat/pat allelic imbalance] [I0.5]
[ADO:allelic dropout] [D0.10]
[ADO variation among sites] [P1.0]
[ADO variation among cells] [Q1.0]
[haploid coverage] [R0.5]
[amplification error] A0.2 0.0001 0
[sequencing error] [E0.00]
[doublet rate per cell] [B0.0]
[error matrix ACGT x ACGT] X0 1 1 1
1 0 1 1
1 1 0 1
1 1 1 0
```

**G#: genotyping error: double [0-1]**

Probability of introducing an error (false positive or false negative;
as a result of amplification/sequencing/calling altogether) per genotype
and cell. In the case of DNA data, the specific error is selected
according to the error matrix (option -X)

If this option is active, no reads are simulated. In fact, it is not
allowed to specify -G and -C (see below) at the same time

**C#: sequencing coverage: integer [1-inf)**

Expected number of reads per site and cell. It is the mean of a Poisson
distribution if there is no overdispersion (option -V). In the case of
overdispersed coverage (option -V > 0), this is the mean of a Negative
Binomial distribution. If there is no allelic imbalance (option -I =
0.5), then the reads are randomly assigned to the maternal/paternal
chromosomes.

**V#: coverage overdispersion: double (0-inf)**

Alpha shape of the gamma distribution in a Negative Binomial
distribution (i.e., *dispersion* parameter) whose mean *mu* is the
sequencing coverage (option -C). The smaller this value is, the more
dispersed is the coverage.

Note: A Negative Binomial distribution can arise as a mixture of Poisson
distributions with mean distributed as a Gamma distribution with scale
parameter (1 - *prob*)/*prob* and shape parameter *dispersion*, where
*prob* = *dispersion*/(*dispersion*+*mu*). The variance is *mu* +
*mu*<sup>2</sup>/*dispersion*.

**I#: maternal / paternal allelic imbalance: double [0-1] **

Probability of a read coming from the maternal chromosome relative to
the probability of coming from the paternal chromosome.

**D#: allelic dropout (ADO) rate: double [0-1]**

Probability of not amplifying one allele, per genotype and cell.

**P#: ADO variation among sites: double (0-inf)**

Amount of ADO rate variation among sites, specified with the alpha shape
parameter of a continuous gamma distribution with mean 1. The smaller
the shape the strongest the heterogeneity introduced.

**Q#: ADO variation among cells: double (0-inf)**

Amount of ADO rate variation among cells, specified with the alpha shape
parameter of a continuous gamma distribution with mean 1. The smaller
the shape the strongest the heterogeneity introduced.

**R#: haploid coverage: double [0-1] **

Expected proportion of reads (*expected* *coverage*, option -C) to be
obtained from a haploid template (e.g., after a deletion or ADO) in
comparison with the expected coverage for a diploid region.

**A# # #: amplification error: double [0-1] double (0-inf) int [0/1]**

Expected probability (*mean*) of introducing an amplification error per
read per site, plus the variance of this error across sites, plus the
amplification model (0: four-template-model; 1: two-template model). For
DNA data, the specific amplification error is selected according to the
error matrix (option -X).

>WARNING: The variance (the second number) has to be bigger than 0 and smaller than *mean* × (1-*mean*).

**E#: sequencing error: double [0-1]**

Probability of introducing a (constant) sequencing error per read site.
For DNA data, the specific sequencing error is selected according to the
error matrix (option -X)

**B#: doublet rate per cell: double [0-1]**

Probability of a doublet per cell.

**X#4×4: error matrix: double [0-inf)**

This 4×4 matrix specifies the relative error rates among the different
nucleotides due to sequencing or amplification errors (the matrix is the
same for both types of errors), in order ACGT for rows and columns. It
does not need to be symmetrical, but all diagonals should be zero.

## 6.5 Output options

All the output files produced by the program are printed inside folder
called by default *results*, located at the same level as the
executable. The program will always output xx data files:

Other outfiles can contain full information on genotypes (-1, -3) or
haplotypes (-2, -4, -9), or about the trees (–6), times (–7). For a
detailed description of these options see below.

```
[OUTPUT]
[print SNV genotypes] 1
[print SNV haplotypes] 2
[print full genotypes] 3
[print full haplotypes] 4
[print ancestors] 5
[print trees] 6
[print times] [7]
[print CATG] [8]
[print true haplotypes] 9
[print replicates in individual folders] v
[print consensus IUPAC haplotypes] [x]
[results folder name] oresults
[user tree file] [Tusertree.0001]
[user genome file] [Uusergenome.fas]
```

**1: print SNV genotypes**

Prints the observed (if no reads are produced) or maximum likelihood
(ML) (if reads are produced, option -C) maternal and paternal genotypes
at variable sites (single nucleotide variants – SNVs), to a file called
*snv_gen* inside the *results* folder, in *Phylip* format
(<http://evolution.genetics.washington.edu/phylip/doc/sequence.html>).
If the option -v is active, each replicate will be printed to an
individual file called *snv_gen.####*, where #### is the
replicate number, in a subfolder *snv_genotypes_dir* inside the
*results* folder.

**2: print SNV haplotypes**

Prints the observed/ML maternal and paternal haplotypes at SNV sites to
a file called *snv_hap* inside the *results* folder. If the option -v
is active, each replicate will be printed to an individual file called
*snv_hap.####*, in a subfolder *snv_haplotypes_dir* inside the
*results* folder. If the option -x is active, a single consensus
sequence with IUPAC codes will be printed (see option -x below).

**3: print full genotypes**

Prints the observed/ML genotypes at all sites to a file called
*full_gen* inside the *results* folder. If the option -v is active,
each replicate will be printed to an individual file called
*full_gen.####*, in a subfolder *full_genotypes_dir* inside the
*results* folder.

**4: print full haplotypes**

Prints the observed/ML maternal and paternal haplotypes at all sites to
a file called *full_hap* inside the *results* folder. If the option -v
is active, each replicate will be printed to an individual file called
*full_hap.####*, in a subfolder *full_haplotypes_dir* inside the
*results* folder. If the option -x is active, a single consensus
sequence with IUPAC codes will be printed (see option -x below).

**5: print MRCA sequences**

Prints genotypes/haplotypes also for the sample MRCA and sample/outgroup
MRCA (see Figure 3).

**6: print trees**

Prints the cell genealogy for each replicate in Newick format.

**7: print times**

Prints for each replicate a list of the coalescent events and branch
lengths in coalescent time.

**8: print CATG**

Prints read counts in CATG format (in addition to VCF). Basically, in the CATG format the first line has the number of samples and the number of sites. Following this is a transposed data matrix with the consensus base calls where each row is a site and each column a different sample. To the right of each site is a tab-separated list of counts of the four bases at that site in the order C, A, T, and G. More information at <http://nbviewer.jupyter.org/gist/dereneaton/d2fd4e70d29f5ee5d195/testing_cat.ipynb#View-the-.cat-results-files>.For example:
```
12 44500
1A0 1B0 1C0 1D0 …
YCCCCCCCCCCC 10,0,10,0 20,0,0,0 20,0,0,0 20,0,0,0 …
GGGGGGGGGGGG 0,0,0,20 0,0,0,20 0,0,0,20 0,0,0,20 …
AAAAAAAAAAAA 0,20,0,0 0,20,0,0 0,20,0,0 0,20,0,0 …
CCCCCCCCCCCC 20,0,0,0 20,0,0,0 20,0,0,0 20,0,0,0 …
AAAAAAAAAAAA 0,20,0,0 0,20,0,0 0, 20,0,0 0,20,0,0 …
```

**9: print true haplotypes**

Prints the true maternal and paternal haplotypes at all sites to a file
called *true_hap* inside the *results* folder. If the option -v is
active, each replicate will be printed to an individual file called
*true_hap.####*, in a subfolder *true_haplotypes_dir* inside the
*results* folder. If the option -x is active, a single consensus
sequence with IUPAC codes will be printed (see option -x below).

**v: print replicates in individual files**

Prints the outputfiles for each replicate

**x: print consensus IUPAC haplotypes**

Prints maternal and paternal sequence as a single consensus sequence
with IUPAC codes as
follows:

| IUPAC codes  |    |   |    |   |   |
| ---------- | ---------- | ---------- | ---------- | ----------- | ---------- |
| A/A => A  | A/C => M  | A/G => R  | A/T => W  | A/_ => a  | A/- => a  |
| C/A => M  | C/C => C  | C/G => S  | C/T => Y  | C/_ => c  | C/- => c  |
| G/A => R  | G/C => S  | G/G => G  | G/T => K  | G/_ => g  | G/- => g  |
| T/A => W  | T/C => Y  | T/G => K  | T/T => T  | T/_ => t  | T/- => t  |
| _/A => a | _/C => c | _/G => g | _/T => t | _/_ => - | _/- => - |
| -/A => a | -/C => c | -/G => g | -/T => t | -/_ => - | -/- => - |

**o$: results folder name**

Name of the folder where the output data files will be stored, at the
same level as the program executable. By default, this folder is called
*results*.

**T$: user tree file**

Name of the file with a rooted tree with branch lengths in Newick format
to be used as the cell genealogy for the simulation of genotypes,
instead of a random coalescent tree.

`((D:0.03,B:0.04):0.02,(C:0.06,A:0.07):0.05);`

**U$: user genome file**

Name of the file with a diploid reference genome in FASTA format, with a
maximum size of 10<sup>6</sup> sites.
```
>maternal genome
AAAAAAACCTGAAACCTTCAGAAGCTGTTTCAGTAGGTAAAGAAAAAGG…
>paternal genome
AAAAAAACCTGAAACCTTCAGAAGCTGTTTCAGTAGGTAAAGAAAAAGG…
```

## 6.6 Other options
```
[OTHER]
[tumor labels] [W]
[noisy] y1
[seed] #1542634309 [if no seed is specified, the current time will be used as seed]
```

**-W: tumor labels**
Use tumor nomenclature for cell labels

| **Meaning**          | **Standard label** | **Tumor label** |
| -------------------- | ------------------ | --------------- |
| sample cell 1        | cell0001           | tumcell0001     |
| outgroup cell        | outgcell           | healthycell     |
| sample MRCA          | ingrroot           | tumoralroot     |
| sample-outgroup MRCA | outgroot           | healthyroot     |

**-y#: noisy**

Amount of information to be printed to the screen:
```
0: does not print anything
1: + simulation summary
2: + replicate information
3: + calculation status and event information.
```

**-##: seed**

Seed for the random number generator. If no seed is specified, the
current time will be used as seed.



#7. Default settings
By default, CellCoal will simulate 10 replicates of samples with 8 cells
(including the outgroup), and 1000 sites, from an exponentially growing
population, with a mutation rate of 1e-07 under a DNA ISM diploid model.
Only the SNV genotypes (“*snv_gen*”) and the genealogies (“*trees*”)
will be printed by default, together with a log file describing the run.
NGS genotype errors/read counts are not simulated by default. To run the
program with the equivalent arguments we would type:

```./cellcoal -n10 -s8 -l100 -e1000 -d0 -b0.0e+00 -r0.0e+00 -k0.0e+00 -h0
-v0.0e+00 -u1.0e-07 -p0.0e+00 -w1 -y1 -g1000```



#8. References
Alexandrov LB et al. 2013. Signatures of mutational processes in human
cancer. Nature 22: 415-421.

DePristo MA et al. 2011. A framework for variation discovery and
genotyping using next-generation DNA sequencing data. Nature Genetics
43: 491–498.

Felsenstein J. 1981. Evolutionary trees from DNA sequences: A maximum
likelihood approach. Journal of Molecular Evolution 17:368-376.

Hasegawa M, Kishino H and Yano T. 1985. Dating the human-ape splitting
by a molecular clock of mitochondrial DNA. Journal of Molecular
Evolution 22:160-174.

Hudson RR. 2002. Generating samples under a Wright-Fisher neutral model
of genetic variation. Bioinformatics 18:337-338.

Jukes TH and Cantor CR. 1969. Evolution of protein molecules. Pp. 21-132
in H. M. Munro, ed. Mammalian Protein Metabolism. Academic Press, New
York, NY.

Lewis PO. 2001. A likelihood approach to estimating phylogeny from
discrete morphological character data. Systematic Biology 50: 913-925.

Korneliussen TS et al. 2013. Calculation of Tajima’s D and other
neutrality test statistics from low depth next-generation sequencing
data. BMC Bioinformatics 14:289.

Tavaré S. 1986. Some probabilistic and statistical problems in the
analysis of DNA sequences. Pp. 57-86 in R. M. Miura, ed. Some
mathematical questions in biology - DNA sequence analysis. Amer. Math.
Soc., Providence, RI.

Yang Z. 1996. Among-site rate variation and its impact on phylogenetic
analysis. Trends in Ecology and Evolution 11:367-372.
