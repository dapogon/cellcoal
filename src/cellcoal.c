/* ****************************************************************************
	GENERAL PUBLIC LICENSE
	This program is free software; you can redistribute it and/or
	modify it under the terms of the GNU General Public License
	as published by the Free Software Foundation; either version 2
	of the License, or (at your option) any later version.
 
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.
 
	You should have received a copy of the GNU General Public License
	along with this program; if not, write to the Free Software
	Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 ******************************************************************************/

/*  cellcoal
//  CellCoal
//
//  Copyright (c) 2016-2022 David Posada. All rights reserved.

	Name:		CellCoal
	Purpose:	To simulate the clonal neutral evolution of single cells (for example within a tumor)
	Code:		David Posada (DP)
	Started:	01 June 2016

 Version 0.1
- started with my code from SNPsim 0.9
- thinnned coalescent without recombination
- attached an additional node to the coalescent MRCA node. This healthy node will be the root for the simulation.
- attached a contemporaneous healthy cell as outgroup
- new data structures
- evolve maternal then paternal genomes over coalescent tree
- redefined output files
- branch length healthy tip and root user defined
- reimplemented ISM diploid (standard) and implemented ISM haploid model
- checked mutation distribution (made a fix for snpsim at the simISMsite function)
- reimplemented fixed SNP model (distributed at random in maternal/paternal genomes)
- reimplemented JC for 0/1 data (=MK2 or CFN model) (made a fix for snpsim at the simJCsite function)
- introduced rate variation among sites
- added rate variation among lineages
- added DNA ISM model
- added rate matrix
- produce read counts
- output format CTAG
- added errors to the NGS reads
- calculate genotype likelihoods (checked to work fine with ANGSD)
- added ADO
- generate reference genome (referenceAllele in struct site)
- implemented HKY+G and nested models (JC69, K80, F81)
- printf FASTA haplotypes with IUPAC codes
- print replicates in different files
- specify results folder name
- add germline variation
- automatically clean previous results directory
- different error probabilities for different bases (Eij)
- can add error directly to genotypes (no reads produced in such a case)
- overdispersed coverage heterogeneity following a negative binomial
- implemented two models for amplification errors. The frequency of wrong templates (max 1; or all possible 3) follows a Beta distribution.
- when ADO we assume there can be a proportion, controlled by a parameter, in the number of reads
- added GTR and GTnR (avoided old redundant calculations)
- print consensus haps for binary data
- simulate point deletions
- added a specific function for JC
- changed letters/symbols for command line arguments
- prints true variants to a file
- can read and use a user tree
- genotype error and ADO now parameterized at the genotype level
- simulate reads for and print all sites in the VCF (changed)
- implemented trinucleotide mutational signatures (alone or mixture)
- implemented user genome
- implemented copy-neutral LOH
- implemented doublets at the NGS reads level
- added option for allelic imbalance
- NGS data: alternate alleles are defined by the reads
- NGS data: SNVs are defined according to the ML genotypes
- NGS data: genotypes and haplotypes show ML genotypes
- VCF: GL is now G10; GL now follows VCF standards
- variable ADO rates among sites and/or cells
- transforming/healthytip branch lengths are not a function of the tumor MRCA / total depth
- option for different taxa names for tumor/non-tumor scenarios
- fixed memory leaks and better printing for VCF (thanks to Alexey Kozlov)
- fixed ISM binary models so germline sites can also accumulate somatic mutations
- new implementation of the 2-template model that corrects the genotype likelihoods
- can use custom name for parameter file

Version 1.0.0 (20/06/2019)
- allow for non integer sequencing coverage

Version 1.1.0 (18/10/2019)
- changed the coalescent time scale from 2N to N
- implemented the Othsuki and Innan 2017 TPB coalescent parameterization with birth and death rates
- genotyping error is now sampled from a Beta-binomial distribution
- allelic imbalance is now sampled from a Beta-binomial distribution
- ADO can now be fixed or sampled from a Beta-binomial distribution, and can vary per cell and/or site
- doublet error is now sampled from a Beta-binomial distribution
- added PL, and also non-normalized versions of G10, GL and PL, to VCF
- compute genotype likelihoods taking into account ADO
- improve genotype likelihood calculations
 
 Version 1.1.1 (21/04/2020)
- print tumor/non-tumor labels in trees
 
 Version 1.1.2 (11/11/2020)
- print PILEUP format
- return 0 after execution
 
 Version 1.2.0 (22/01/2021)
- calculate genotype likelihoods for the phased genotypes and then combine them to obtain the likelihood for the unphased genotypes (important when amplification errors for maternal and paternal sites are distinct)
 
 Version 1.2.1 (22/03/2021)
 - improved memory allocation
 - print all likelihoods for doublets also
 - update example parameter files
 
 Version 1.2.2 (24/03/2021)
 - minor fixes
 
 Version 1.3.0 (11/01/2022)
 - implement tree rate switches
 
 Version 1.3.1 (04/04/2022)
 - fixed bug when using user trees that prevented the calculation of the total tree lenght
 - fixing -P and -Q options when printing help (were interchanged only there)
 
 Version 1.3.3 (01/10/2024)
- fix AllelicDropout so alleleADOrateSite is sampled once and stored for reuse across all cell genotypes. (Thanks to Jerry Barber)
 
 Version 1.3.4 (10/02/2025)
- fix nonISMRelMutRate input so it can be > 1
	(Thanks to Haotian Zhang)

 Version 2.0.3 (14/11/2026)
- faster and lighter genotype handling: allele counting traverses the genomes cell by cell, and genomes are stored with one byte per site (output identical to 2.0.2)

 Version 2.0.2 (14/11/2026)
- fix germline variation (-c): the SNPs were added to the root genome only, so the cells (and the healthy cell) never inherited them; they are now passed down to all nodes

 Version 2.0.1 (14/11/2026)
- memory: the per-cell, per-site blocks are only allocated when an option needs them (reads, haplotypes, full genotypes, ancestors); about 66% less peak memory in runs without reads, with identical output
- fix memory that was not freed at exit (user tree names, genetic signature tables); the example parameter files are now clean under AddressSanitizer/LeakSanitizer

 Version 2.0.0 (13/11/2026)
- fix genotyping error (-G): rows of the error matrix are now normalized, so the erroneous base is drawn with the intended probabilities
- fix allelic dropout (-P/-Q): ADO rates now vary correctly between cells and between sites
- fix memory errors and leaks (user genome/tree reading, signatures), input checks (missing option arguments, buffer bounds, realloc/calloc), VCF G10/G10N/GL/GLN/PL/PLN print '.' when the likelihood is exactly zero (log10 = -inf, e.g. with no sequencing error)
- VCF PL/PLN follow the VCF specification: -10*log10 (non-negative integers, 0 = most likely genotype); before they were 10*log10 (negative)
- much faster and leaner: compact genotype-likelihood storage, budgeted genotype-likelihood cache (env CELLCOAL_GL_CACHE_MB), exact log10 memoization, faster VCF formatting, OpenMP parallel read-count generation and VCF printing
- read counts use independent per-site random streams: results for a given seed are reproducible and independent of the number of threads, but differ from versions <= 1.3.4

 
[TO-DOs]
- add new signatures
- move data structure into cell structure (optional)
- dated tips (optional)
- population structure (optional)
- change random number generator (optional)

[MODELS implemented]
0:ISM_01
1:ISMhap_01
2:M2_01
3:ISM_DNA
4:ISMhap_DNA
5:HKY_JC_DNA
6:GTR_DNA
7:GTRnrDNA
*/

#include "cellcoal.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "definitions.h"

/* The per-cell/per-site blocks (CS(i,j), 48 bytes each) are only needed for read counts
   (NGS, doublets, VCF/CATG/PILEUP) and for haplotype/full-genotype/ancestor outputs. Without them,
   SNV genotypes are printed from the genome data, and allocating (cells+1) x sites blocks anyway costs
   e.g. ~2.4 GB for 1000 cells x 50000 sites. */

/* adds a value to the statistics accumulated over replicates */
void AccumAdd (Accum *a, double x)
	{
	a->sum += x;
	a->sumSq += pow(x,2);
	}

/* mean and (sample) variance after n replicates */
void AccumFinish (Accum *a, int n)
	{
	a->mean = a->sum / (double) n;
	a->var = (1.0 / (double) (n-1)) * (a->sumSq - pow(a->sum,2) / (double) n);
	}

/* default values of all the options */
static void SetDefaultOptions (void)
	{
    numDataSets = 50;			/* the number of samples to simulate */
    numCells = 20;				/* number of ingroup cells in each data set */
    ploidy = 2;                 /* we assume diploid genomes */
    numSites = 1000;			/* number of sites to simulate = length of the genome */
    N = 10000;					/* effective population size */
    numPeriods = 0;				/* number of distinct demographic periods */
    doDemographics = NO;		/* whether to implement demographics */
    doExponential = NO;			/* whether to do exponential growth */
    growthRate = 1.0e-3;		/* rate for the exponential population growth */
	birthRate = 1;				/* birth rate in the Ohtuski and Innan 2017 coalescent model*/
	deathRate = 0;				/* death rate in the Ohtuski and Innan 2017 coalescent model*/
   	mutationRate = 1.0e-7;		/* nucleotide mutation rate per site per generation */
	SNPrate = 0.0;				/* germline variation rate for rooth healthy genome */
	rateVarAmongLineages = NO;	/* modify rate variation among branches (to non-clock) */
	doTreeRateSwitches = NO;	/* do rate switches along the tree */
	numTreeRateSwitches = 0;	/* number of rate switches along the tree */
	alphabet = DNA;          	/* alphabet 0/1 or DNA */
    altModel = 0;				/* by default the alternative model will be ISM haploid */
    propAltModelSites = 0;		/* proportion of sites that will mutate according to alternative model */
    nonISMRelMutRate = 1.0;		/* relative rate alternative/default model for sites */
    equalBaseFreq = YES;        /* DNA base frequencies */
	freq[0] = freq[1] = freq[2] = freq[3] = 0.25;
	titv = 0.5;					/* transition/transversion rate ratio */
	thereIsMij = NO;			/* mutation rate matrix*/
	Mij[0][0] = Mij[1][1] = Mij[2][2] = Mij[3][3] = 0;  /* mutation probabilities */
	Mij[0][1] = Mij[0][2] = Mij[0][3] = 1.0;
	Mij[1][0] = Mij[1][2] = Mij[1][3] = 1.0;
	Mij[2][0] = Mij[2][1] = Mij[2][3] = 1.0;
	Mij[3][0] = Mij[3][1] = Mij[3][2] = 1.0;
	thereIsEij = NO;			/* error rate matrix*/
	Eij[0][0] = Eij[1][1] = Eij[2][2] = Eij[3][3] = 0;  /* sequencing error probabilities */
	Eij[0][1] = Eij[0][2] = Eij[0][3] = 1.0;
	Eij[1][0] = Eij[1][2] = Eij[1][3] = 1.0;
	Eij[2][0] = Eij[2][1] = Eij[2][3] = 1.0;
	Eij[3][0] = Eij[3][1] = Eij[3][2] = 1.0;
	doJC = YES;
	doHKY = NO;
	doGTR = NO;
	doGTnR = NO;
	rateVarAmongSites = NO;     	/* rate variation among different sites along the genome */
    alphaSites = infinity;      	/* alpha shape of the gamma distribution for rate variation among sites */
    alphaBranches = infinity;   	/* alpha shape of the gamma distribution for rate variation among lineages */
    alphaCellCoverage = 0;			/* between-cell coverage dispersion (-J): 0 = none */
    alphaSiteCoverage = 0;			/* between-site coverage dispersion (-O): 0 = none */
    alphaCoverage = infinity;   	/* alpha shape of the gamma distribution for coverage */
	doSimulateFixedNumMutations = NO;	/* whether to simulate a fixed number of mutations */
	doUserTree = NO;				/* whether to assume a user tree instead od making the coalescent */
	doUserGenome = NO;				/* whether to use a user genome instead of a simulated one */
    doPrintSNVgenotypes = NO;		/* whether to print SNVs */
    doPrintSNVhaplotypes = NO;  	/* whether to print haplotypes */
    doPrintTrueHaplotypes = NO;  	/* whether to print the true haplotypes s */
    doPrintFullHaplotypes = NO;		/* whether to print sequences */
    doPrintFullGenotypes = NO;		/* whether to print all genotypes (variable + invariable) */
    doPrintTree = NO;				/* whether to print the coalescent tree */
    doPrintTimes = NO;				/* whether to print coalescent times */
	doPrintAncestors = NO;      	/* whether to print data for ancestral cells */
	doNGS = NO;  					/* produce reads by default */
	doTumorNames = NO;				/* use specific name for taxa in the tumor scenario */
	stringPrecision = 12;			/* precision for taxa names */
	doPrintCATG = NO;				/* whether to print read counts for SNVs in CATG format*/
	doPrintPILEUP = NO;				/* whether to print read counts for SNVs in PILEUP format*/
	doPrintBulk = NO;				/* whether to print pooled (bulk) read counts and VAF per site */
	doPrintVCF = YES;				/* whether to print the VCF file when reads are simulated */
	doZipVCF = NO;					/* whether the VCF is written compressed with gzip */
	minVAFvcf = 0.0;				/* minimum pooled VAF of the sites in the VCF */
	doSimulateData = YES;			/* whether to simulate any data (or just look at the expectations; useful for debugging) */
	doPrintSeparateReplicates = NO; /* whether to put every replica in its own file */
	doPrintIUPAChaplotypes = NO;	/* whether to print IUPAC halotypes */
	doGeneticSignatures = NO;		/* whether to use a genetic signature to model trinucleotide mutations */
	numUserSignatures = 0;			/* by default we do not use a genetic signature */
	transformingBranchLengthRatio = 0.5; /* ratio of the transforming branch length compare to the tumor MRCA length */
	healthyTipBranchLengthRatio = 1.0; 	/* ratio the branch leading to the healthy cell to the (tumor MRCA + transfoming branch) length*/
 	coverage = 5;					/* NGS  depth for read counts */
	rateVarCoverage = NO;			/* there is coverage dispersion */
	fixedADOrate = 0;				/* allelic dropout */
	thereisADO = NO;				/* allelic dropout flag */
	doADOcell = NO;					/* ADO does not change per cell by default */
	doADOsite = NO;					/* ADO does not change per site by default */
	meanADOcell = 0;     			/* mean of beta distribution for ADO expected for a cell */
	varADOcell = 1.0e-06;     		/* var of beta distribution for ADO expected for a cell */
	meanADOsite = 0;				/* mean of beta distribution for ADO expected for a site */
	varADOsite = 1.0e-06; ;			/* var of beta distribution for ADO expected for a site */
	sequencingError = 0;			/* NGS error rate */
	meanGenotypingError = 0;		/* mean of beta distribution for errors added directly in the genotypes */
	varGenotypingError = 1.0e-06;   /* variance of beta distribution for genotype errors */
	meanAmplificationError = 0; 	/* mean of beta distribution for WGA errors */
	varAmplificationError = 1.0e-06; /* variance of beta distribution for WGA errors */
	simulateOnlyTwoTemplates = NO;	/* whether simulating maximum of two templates after single-cell amplification, or there can be all four */
	haploidCoverageReduction = 0.5; /* proportion of reads produced when a single allele is present */
	doAllelicImbalance = NO;		/* by default there is no allelic imbalance */
	meanAllelicImbalance = 0.5;		/* beta mean proportion of maternal/ paternal reads */
	varAllelicImbalance = 1.0e-06;	/* beta var proportion of maternal/ paternal reads */
	meanDoubletRate = 0;			/* mean of beta distribution for doublet errors */
	varDoubletRate = 1.0e-06;		/* variance of beta distribution for doublet errors */
	numNodes = 1000;            	/* initial number of nodes allocated to build the coalescent trees */
	doSpecificParameterFile = NO;	/* by default do not use a custom parameter file */
	userSeed = 0;					/* seed entered by the user */
	noisy = 1;						/* level of information to be printed in the screen (see below) */
	}

/* checks that need several options at once (after reading the command line) */
static void CheckOptionCombinations (void)
	{
	if (doBCF == YES)
		{
		if (doZipVCF == YES || doPrintVCF == NO)
			{
			fprintf (stderr, "\n!!! PARAMETER ERROR: --bcf replaces the VCF file: it cannot be combined with -z or -Z\n\n");
			exit (-1);
			}
		if (numDataSets > 1 && doPrintSeparateReplicates == NO)
			{
			fprintf (stderr, "\n!!! PARAMETER ERROR: --bcf with several replicates needs -v (one BCF file per replicate)\n\n");
			exit (-1);
			}
		}
	if (covAdoCorr != 0.0 && (alphaCellCoverage <= 0.0 || doADOcell == NO))
		{
		fprintf (stderr, "\n!!! PARAMETER ERROR: --cov-ado-corr needs coverage variation among cells (-J) and ADO variation among cells (-P)\n\n");
		exit (-1);
		}
	if (covDoubletCorr != 0.0 && (alphaCellCoverage <= 0.0 || meanDoubletRate <= 0))
		{
		fprintf (stderr, "\n!!! PARAMETER ERROR: --cov-doublet-corr needs coverage variation among cells (-J) and doublets (-B)\n\n");
		exit (-1);
		}
	if (adoDoubletCorr != 0.0 && (doADOcell == NO || meanDoubletRate <= 0))
		{
		fprintf (stderr, "\n!!! PARAMETER ERROR: --ado-doublet-corr needs ADO variation among cells (-P) and doublets (-B)\n\n");
		exit (-1);
		}
	}

/* means and variances over replicates */
static void ComputeSummaryStatistics (void)
	{
	AccumFinish (&runStats.mu, numDataSets);
	AccumFinish (&runStats.snv, numDataSets);
	AccumFinish (&runStats.del, numDataSets);
	AccumFinish (&runStats.tmrca, numDataSets);
	AccumFinish (&runStats.mlErr, numDataSets);
	AccumFinish (&runStats.called, numDataSets);
	runStats.ca.mean = runStats.ca.sum / (double) numDataSets;	/* the variances of these two are not computed */
	runStats.cnloh.mean = runStats.cnloh.sum / (double) numDataSets;
	}

/* reads simulated, summarized and printed (VCF, PILEUP, bulk, CATG) in blocks of sites */
static void GenerateStreamedReads (long int *seed, int dataSetNum)
	{
	/* reads are simulated, summarized and printed (VCF, PILEUP) in blocks of sites */
	if (doPrintVCF == YES && doPrintSeparateReplicates == NO)
		fprintf (fpVCF, "[#%d]\n", dataSetNum+1);
	if (doPrintCATG == YES)
		{
		catgBody = tmpfile ();
		if (!catgBody)
			{
			fprintf (stderr, "Could not open a temporary file for the CATG output\n");
			exit (-1);
			}
		}
	StreamReads (seed, (doPrintVCF == YES) ? fpVCF : NULL, (doPrintPILEUP == YES) ? fpPILEUP : NULL, (doPrintBulk == YES) ? fpBulk : NULL, catgBody);
	if (doPrintCATG == YES)
		{
		char	cbuf[65536];
		size_t	n;

		PrintCATGHeader (fpCATG);
		rewind (catgBody);
		while ((n = fread (cbuf, 1, sizeof(cbuf), catgBody)) > 0)
			fwrite (cbuf, 1, n, fpCATG);
		fprintf (fpCATG, "\n");
		fclose (catgBody);
		catgBody = NULL;
		if (doPrintSeparateReplicates == NO)
			fprintf (fpCATG, "[#%d]\n", dataSetNum+1);
		}
	if (doPrintPILEUP == YES && doPrintSeparateReplicates == NO)
		fprintf (fpPILEUP, "[#%d]\n", dataSetNum+1);
	if (doPrintBulk == YES && doPrintSeparateReplicates == NO)
		fprintf (fpBulk, "[#%d]\n", dataSetNum+1);
	}

/* per-replicate messages and output of genotypes, haplotypes and reads */
static void PrintReplicateOutputs (int dataSetNum)
	{
	if (noisy > 1)
	{
	fprintf (stderr, "\nData set %d",dataSetNum+1);
	fprintf (stderr, "\nNumber of coalescence events   =   %d", numCA);
	fprintf (stderr, "\nNumber of mutational events    =   %d", numMU);
	fprintf (stderr, "\nNumber of SNVs                 =   %d", numSNVs);
	fprintf (stderr, "\nNumber of deletion events      =   %d", numDEL);
	fprintf (stderr, "\nNumber of cnLOH events         =   %d", numCNLOH);
	
	if (doTreeRateSwitches == YES)
		fprintf(stderr, "\nTree rate switches             =   %d %s", numTreeRateSwitches, treeRateInfo);

	fprintf (stderr, "\n\n");
	}

	if (doPrintSNVgenotypes == YES && numSNVs > 0) /* we only print replicates with variation */
		{
		if (doPrintSeparateReplicates == NO)
			fprintf (fpSNVgenotypes, "[#%d]\n", dataSetNum+1);
		PrintSNVGenotypes(fpSNVgenotypes);
		}
	if (doPrintSNVhaplotypes == YES && numSNVs > 0) /* we only print replicates with variation */
		{
		if (doPrintSeparateReplicates == NO)
			fprintf (fpSNVhaplotypes, "[#%d]\n", dataSetNum+1);
		PrintSNVHaplotypes(fpSNVhaplotypes, NO);
		}
	if (doPrintFullGenotypes == YES)
		{
		if (doPrintSeparateReplicates == NO)
			fprintf (fpFullGenotypes, "[#%d]\n", dataSetNum+1);
		PrintFullGenotypes(fpFullGenotypes);
		}
	if (doPrintFullHaplotypes == YES)
		{
		if (doPrintSeparateReplicates == NO)
			fprintf (fpFullHaplotypes, "[#%d]\n", dataSetNum+1);
		PrintFullHaplotypes(fpFullHaplotypes);
		}
	CloseFullStore ();

	if (alphabet == DNA && streamMode == NO && doNGS == YES)
		{
		if (doPrintVCF == YES)
			{
			if (doPrintSeparateReplicates == NO)
				fprintf (fpVCF, "[#%d]\n", dataSetNum+1);
			PrintVCF (fpVCF);
			}
		if (doPrintCATG == YES)
			{
			PrintCATG(fpCATG);
			if (doPrintSeparateReplicates == NO)
				fprintf (fpCATG, "[#%d]\n", dataSetNum+1);
			}
		if (doPrintPILEUP == YES)
			{
			PrintPILEUP(fpPILEUP);
			if (doPrintSeparateReplicates == NO)
				fprintf (fpPILEUP, "[#%d]\n", dataSetNum+1);
			}
		if (doPrintBulk == YES)
			{
			PrintBulk (fpBulk);
			if (doPrintSeparateReplicates == NO)
				fprintf (fpBulk, "[#%d]\n", dataSetNum+1);
			}
		}
	}

/* one replicate: genealogy, mutations, deletions/cnLOH, ADO, errors, reads and all the per-replicate outputs */
static void RunReplicate (int dataSetNum, long int *seed)
	{
	int		i, j;

        numCA = numMU = numDEL = numCNLOH = numProposedMU = TMRCA = numSNVmaternal = 0;

		if (doPrintSeparateReplicates == YES)
			PrepareSeparateFiles(dataSetNum);

		if (noisy == 1)
            {
            fprintf (stdout, "\rReplicate #%3d/%d ", dataSetNum+1, numDataSets);
			fflush (stdout);
			
			/* if stdout is redirected elsewhere, print also to stderr */
			if (!isatty(fileno(stdout)))
				{
				fprintf (stderr, "\rreplicate #%3d/%d", dataSetNum+1, numDataSets);
				fflush (stderr);
				}
			}
			
        /* Generate coalescence tree */
		if (doUserTree == NO)
			{
			if (noisy > 2)
				fprintf (stderr, "\n>> Starting coalescent tree ...");
			MakeCoalescenceTree (numCells, N, seed);
			if (noisy > 2)
				fprintf (stderr, "\n>> Finishing coalescent tree ... DONE");

		   /* Make ingroup cell tree non-clock if needed */
			if (rateVarAmongLineages == YES)
				MakeTreeNonClock (coalTreeMRCA, seed);
			else if (doTreeRateSwitches == YES)
				{
				treeRateInfo = (char *) calloc(MAX_LINE, sizeof(char));
				TreeRateSwitch (coalTreeMRCA, seed);
				}
			
			runStats.ca.sum += numCA;
			AccumAdd (&runStats.tmrca, TMRCA);
			}

		totalTreeLength = SumBranches (healthyRoot);

        if (doPrintTree == YES)
            PrintTree (healthyRoot, fpTrees);

        if (doPrintTimes == YES)
            PrintTimes (0, fpTimes);
 
 
	DecideSparseMode ();

	/* allocate genotype data (see GetG and SetG) */
        data = (Geno ***) calloc (ploidy, sizeof(Geno **));
        if (!data)
            {
            fprintf (stderr, "Could not allocate the data structure\n");
            exit (-1);
            }
        for (i=0; i<ploidy; i++)
            {
            data[i] = (Geno **) calloc (2*numCells+1, sizeof(Geno *));
            if (!data[i])
                {
                fprintf (stderr, "Could not allocate the data[] structure\n");
                exit (-1);
                }
            for (j=0; j<2*numCells+1; j++)
                {
                if (sparseMode == YES && j != OUTGROUP_ROOT)
                    continue;	/* in sparse mode only the genome of the root is stored densely */
                data[i][j] = (Geno *) calloc (numSites, sizeof(Geno));
                if (!data[i][j])
                    {
                    fprintf (stderr, "Could not allocate the data[][] structure\n");
                    exit (-1);
                    }
                }
            }
        if (sparseMode == YES)
            {
            numNodes = 2*numCells+1;
            sparseCols = (Geno **) calloc (numSites, sizeof(Geno *));
            if (!sparseCols)
                {
                fprintf (stderr, "Could not allocate the sparseCols structure\n");
                exit (-1);
                }
            }

		/* allocate cell structure */
		AllocateCellStructure();

        /* allocate memory for site information (equal for maternal and paternal) */
        allSites = (SiteStr*) calloc (numSites, sizeof(SiteStr));
        if (!allSites)
            {
            fprintf (stderr, "Could not allocate the allSites structure\n");
            exit (-1);
            }
		for (i=0; i<numSites; i++)
			{
			allSites[i].alternateAlleles = (int *) calloc (4, sizeof(int));
			if (!allSites[i].alternateAlleles)
                {
                fprintf (stderr, "Could not allocate the allSites[i].alternateAlleles[] structure\n");
                exit (-1);
                }
			}

        /* the arrays below keep the index for different types of sites */
        SNVsites = (int*) calloc (numSites, sizeof(int));
        if (!SNVsites)
            {
            fprintf (stderr, "Could not allocate the SNVsites structure\n");
            exit (-1);
            }

        /* the arrays below keep the index for different types of sites */
        variantSites = (int*) calloc (numSites, sizeof(int));
        if (!variantSites)
            {
            fprintf (stderr, "Could not allocate the variantSites structure\n");
            exit (-1);
            }

        DefaultModelSites = (int*) calloc (numSites, sizeof(int));
        if (!DefaultModelSites)
            {
            fprintf (stderr, "Could not allocate the DefaultModelSites structure\n");
            exit (-1);
            }
			
        AltModelSites = (int*) calloc (numSites, sizeof(int));
        if (!AltModelSites)
            {
            fprintf (stderr, "Could not allocate the AltModelSites structure\n");
            exit (-1);
            }
 
        if (doSimulateData == YES)
            {
            numISMmutations = 0;
			numISMdeletions = 0;
 			numISMCNLOH = 0;
 
			if (sparseMode == YES)
				{
				ResetObservationLayer ();
				ResetSparseColumns ();
				}

             /* assign ancestral states to all genomes */
            InitializeGenomes (healthyRoot, seed);
			
			/* add germline variation if needed */
			if (SNPrate > 0)
				{
				AddGermlineVariation (seed);
				/* germline SNPs were added to the root genome only: pass them down to all the nodes */
				InheritRootGenome (healthyRoot);
				}
	
            /* evolve maternal genome */
            if (noisy > 2)
                fprintf (stderr, "\n>> Evolving maternal genome ... ");
             EvolveSitesOnTree (healthyRoot, MATERNAL, seed);
            if (noisy > 2)
                fprintf (stderr, "DONE");

            /* evolve paternal genome */
             if (noisy > 2)
                fprintf (stderr, "\n>> Evolving paternal genome ... ");
            EvolveSitesOnTree (healthyRoot, PATERNAL, seed);
            if (noisy > 2)
                fprintf (stderr, "DONE\n\n");
	
            AccumAdd (&runStats.mu, numMU);

			if (CNLOHrate > 0)
				{
				/* evolve maternal CN_LOH */
				if (noisy > 2)
					fprintf (stderr, "\n>> Evolving maternal CN_LOH ... ");
				 EvolveCNLOHonTree(healthyRoot, MATERNAL, seed);
				if (noisy > 2)
					fprintf (stderr, "DONE");

				/* evolve paternal CN_LOH  */
				 if (noisy > 2)
					fprintf (stderr, "\n>> Evolving paternal CN_LOH ... ");
				 EvolveCNLOHonTree(healthyRoot, PATERNAL, seed);
				if (noisy > 2)
					fprintf (stderr, "DONE\n\n");

				AccumAdd (&runStats.cnloh, numCNLOH);
				}

			if (deletionRate > 0)
				{
				/* evolve maternal deletions */
				if (noisy > 2)
					fprintf (stderr, "\n>> Evolving maternal deletions ... ");
				 EvolveDeletionsOnTree(healthyRoot, MATERNAL, seed);
				if (noisy > 2)
					fprintf (stderr, "DONE");

				/* evolve paternal deletions  */
				 if (noisy > 2)
					fprintf (stderr, "\n>> Evolving paternal deletions ... ");
				 EvolveDeletionsOnTree(healthyRoot, PATERNAL, seed);
				if (noisy > 2)
					fprintf (stderr, "DONE\n\n");

				AccumAdd (&runStats.del, numDEL);
				}

			/* count how many SNVs we observe before ADO and genotype errors, but after CNLOH and deletion */
			numSNVs = CountTrueVariants();
			numSNVsTrue = numSNVs;

			/* copy true genotypes before doing ADO but after CNLOH and deletion */
			if (NEED_CELL_SITES && streamMode == NO)
			for (i=0; i<numCells+1; i++)
				for (j=0; j<numSites; j++)
					{
					CS(i,j).trueMaternalAllele = GetG(MATERNAL, i, j);
					CS(i,j).truePaternalAllele = GetG(PATERNAL, i, j);
					}

			/* print reference haplotypes without errors */
			if (doPrintTrueHaplotypes == YES)
				{
  				if (doPrintSeparateReplicates == NO)
					fprintf (fpTrueHaplotypes, "[#%d]\n", dataSetNum+1);
                PrintTrueFullHaplotypes(fpTrueHaplotypes);
				}

			/* do alellic dropout */
			adoReplicate = dataSetNum;
			if (fixedADOrate > 0 || doADOcell == YES || doADOsite == YES)
				{
				if (sparseMode == YES)
					SparseAllelicDropout(seed);
				else
					AllelicDropout(seed);
				}
	
			/* introduce errors directly in the genotypes */
			if (meanGenotypingError > 0)
				{
				if (sparseMode == YES)
					SparseGenotypeError(seed);
				else
					GenotypeError(seed);
				}

			if (sparseMode == YES)
				obsActive = YES;	/* from now on the cell genotypes are observed with ADO and errors */
	
			/* generate NGS read counts */
			if (streamMode == YES)
				{
				GenerateStreamedReads (seed, dataSetNum);
				}
			else if (doNGS == YES)
				GenerateReadCounts(seed);

			/* count how many alleles are at each site and how many SNVs we observe now, after ADO and genotype errors/read generation */
			if (streamMode == YES)
				;	/* already counted by StreamReads */
			else if (doNGS == YES)
				numSNVs = CountAllelesInMLGenotypes();
			else
				numSNVs = CountAllelesInObservedGenotypes();

            AccumAdd (&runStats.snv, numSNVs);

			if (fpSummary != NULL)
				{
				double	metric[NUM_TRUTH_METRICS];
				int		k;

				TruthMetrics (metric);
				fprintf (fpSummary, "%d\t%.6f\t%d\t%d\t%d\t%d\t%d\t%.6f", dataSetNum+1, TMRCA, numMU, numDEL, numCNLOH, numSNVsTrue, numSNVs,
					numSNVsTrue > 0 ? numSNVs / (double) numSNVsTrue : 0.0);
				for (k=0; k<NUM_TRUTH_METRICS; k++)	/* NA when there are no reads or the metric is undefined */
					{
					if (metric[k] < 0.0)
						fprintf (fpSummary, "\tNA");
					else
						fprintf (fpSummary, "\t%.6f", metric[k]);
					}
				fprintf (fpSummary, "\n");
				fflush (fpSummary);
				}

			if (numSNVs == 0)
                zeroSNVs++;
		
			PrintReplicateOutputs (dataSetNum);
			} //simdata
 
		#ifdef CHECK_MUT_DISTRIBUTION
        /* Check expectations regarding the position of the mutations in the tree */
        /* Find out how many mutations are per site; we expect theta/i number of sites with i mutations
           This approximation might not be good for the haploid ISM model
           We will count how many sites have 1, 2, 3 ... numCells-1 mutant alleles */
            SiteMut = (int *) calloc(numSites, sizeof(int));
            for (j=0; j<numSites; j++)
                for (i=0; i<numCells; i++)
                    {
                    if (GetG(0, i, j) == 1)
                        SiteMut[j]++;
                    if (GetG(1, i, j) == 1)
                        SiteMut[j]++;
                     }
            for (j=0; j<numSites; j++)
                for (i=0; i<numCells; i++)
                    if (SiteMut[j] == i)
                        MutCount[i]++;
			
            /* Check expectations regarding the distribution of the mutations along the genome */
            /* the average sum of the positions of the segregating sites has to be L/2 */
            sumPos = 0;
            if (numSNVs > 0)
                {
                dataSetsWithSNVs++;
                for (j=0; j<numSites; j++)
                    if (allSites[j].numMutations > 0)
                        sumPos += (j+1);
                meansumPos += (sumPos / (double) numSNVs);
                }
        #endif

		if (noisy > 2)
			for (i=0; i<numSites; i++)
				PrintSiteInfo (stderr, i);
	
		
		#ifdef COUNT_ML_GENOTYPE_ERRORS
			CountMLGenotypingErrors();
		#endif
		
		// free all the necessary stuff for this replicate
		if (doUserTree == NO)
			{
            for (i=0; i<numNodes; i++)
	            free(treeNodes[i].name);
	        free (treeNodes);
	        }
		
		if (doPrintSeparateReplicates == YES)
			{
			if (doPrintTree == YES)
				fclose (fpTrees);
			if (doPrintTimes == YES)
				fclose (fpTimes);
			}
		
        if (doSimulateData == YES)
            {
			/* close file pointers for separate replicate file */
			if (doPrintSeparateReplicates == YES)
				{
				if (doPrintSNVgenotypes == YES)
					fclose (fpSNVgenotypes);
				if (doPrintFullGenotypes == YES)
					fclose (fpFullGenotypes);
				if (doPrintSNVhaplotypes == YES)
					fclose (fpSNVhaplotypes);
				if (doPrintFullHaplotypes == YES)
					fclose (fpFullHaplotypes);
				if (doPrintTrueHaplotypes == YES)
					fclose (fpTrueHaplotypes);
				if (doNGS == YES && doPrintVCF == YES)
					CloseVCF (fpVCF);
				if (doPrintCATG == YES)
					fclose (fpCATG);
				if (doPrintPILEUP == YES)
					fclose (fpPILEUP);
				if (doPrintBulk == YES)
					fclose (fpBulk);
				}
 
			/* free data structures */
			for (i=0; i<ploidy; i++)
                {
                for (j=0; j<2*numCells+1; j++)
                    free (data[i][j]);	/* NULL in sparse mode, except for the root */
                free (data[i]);
                }
			free (data);
			if (sparseMode == YES)
				{
				ResetObservationLayer ();
				ResetSparseColumns ();
				free (sparseCols);
				sparseCols = NULL;
				}

			if (doTreeRateSwitches == YES)
				free (treeRateInfo);
		 
			for (i=0; i<numSites; i++)
              free(allSites[i].alternateAlleles);
            free (allSites);
            free (SNVsites);
            free (mlStore);
            free (cellCovMult);
            cellCovMult = NULL;
            mlStore = NULL;
            mlStoreCap = 0;
			free (variantSites);
            free (DefaultModelSites);
            free (AltModelSites);

			if (doGeneticSignatures == YES)
				FreeTriNucleotides ();
			}
		
		/* free the per-cell structures (they are allocated at every replicate, with or without NGS) */
		for (i=0; i<numCells+1; i++)
			{
			free (cell[i].doublet);
			free (cell[i].site);
			}
		free (cell);
		if (siteGL)
			{
			for (j=0; j<numSites; j++)
				free (siteGL[j]);
			free (siteGL);
			siteGL = NULL;
			}
		
        /* if(((dataSetNum*10)%numDataSets) == 0)
         PrintMemory(stderr); */
	}

int main (int argc, char **argv)
    {
    FILE 		*fp;
    int			i, j;
    int         dataSetNum;
    long int	seed;
    double		a, a2;
    clock_t		start;
    double		secs, wallSecs = 0.0;
		
    /* Default settings */
	seed = time(NULL);				/* seed for random numbers */
	SetDefaultOptions ();

	/*
	 noisy = 0: does not print anything
	 noisy = 1: + simulation summary,
	 noisy = 2: + replicate information,
	 noisy = 3: + calculation status and event information
	*/

	strcpy(parameterFile, "parameters");
	argv = ExpandConfigArgs (&argc, argv);	/* --config file.json|yaml */
	ReadParametersFromCommandLine (argc, argv);
	CheckOptionCombinations ();
	if (argc < 2 || doSpecificParameterFile == YES)
        {
		readingParameterFile = YES;
        if ((fp = freopen(parameterFile, "r", stdin)) != NULL)
            ReadParametersFromFile();
        else
            {
			fprintf (stderr, "\nERROR: Could not open parameter file: %s\n", parameterFile);
            PrintUsage(stderr);
            exit(-1);
            }
        }
	else
		readingParameterFile = NO;
   
   start = clock();
   clock_gettime (CLOCK_MONOTONIC, &wallStart);
	
    if (noisy > 0)
        {
        PrintHeader(stderr);
        PrintDate(stderr);
 		}
		
    if (doSimulateData == NO)
        {
        doPrintSNVgenotypes = NO;
        doPrintSNVhaplotypes = NO;
		doPrintTrueHaplotypes = NO;
		doPrintFullHaplotypes = NO;
        doPrintFullGenotypes = NO;
		doPrintAncestors = NO;
		doNGS = NO;
		doPrintCATG = NO;
		doPrintPILEUP = NO;
		doPrintBulk = NO;
		doPrintVCF = YES;
		doZipVCF = NO;
		minVAFvcf = 0.0;
   		}
	
    /* initialize a few variables */
	numNodes = 2 * numCells + 1;
	maxTreeRateSwitches = numCells;
	a = 0;
    a2 = 0;
    for (i=1; i<numCells; i++) /* ignoring outgroup (healthy cell) */
        {
        a += 1.0/(double)i;
        a2 += 1.0/(double)pow(i,2);
        }
    theta = 2.0 * N * mutationRate * 2.0 * numSites * ((1.0-propAltModelSites) + nonISMRelMutRate*propAltModelSites);
    expNumMU = theta * a;
    expVarNumMU =  theta * a + pow(theta,2) * a2;
    expTMRCA = 2 * (1 - 1.0/numCells) /* 2.0 * N */;
    expVarTMRCA = 0;
	for (i=2; i<=numCells; i++)
		expVarTMRCA += /* 4.0 * pow(N,2) */ 4.0 / (pow(i,2) * pow(i-1,2));

	memset (&runStats, 0, sizeof(runStats));
	zeroSNVs = 0;
		
    altModelMutationRate = mutationRate*nonISMRelMutRate;
	
	if (doSimulateData == YES && alphabet == DNA)
		{
		/* initialize cumfreq */
		cumfreq[0] = freq[0];
		for (i=1; i<4; i++)
			cumfreq[i] = cumfreq[i-1] + freq[i];
	
		/* initialize cumMij */
		for (i=0; i<4; i++)
			{
			cumMij[i][0] = Mij[i][0];
			for (j=1; j<4; j++)
				cumMij[i][j] = cumMij[i][j-1] + Mij[i][j];
			}

		/* initialize cumEij */
		for (i=0; i<4; i++)
			{
			cumEij[i][0] = Eij[i][0];
			for (j=1; j<4; j++)
				cumEij[i][j] = cumEij[i][j-1] + Eij[i][j];
			}
		}

	if (doOthsukiInnanCoal == YES)
		{
		growthRate = birthRate - deathRate;
		if (growthRate < 0)
			{
			fprintf (stderr, "\n!!! PARAMETER ERROR: the birth rate (%f) should be bigger than the death rate %f)", birthRate, deathRate);
			PrintUsage(stderr);
			}
		}
	
    /* Set seed and spin wheels of pseudorandom number generator */
    /* seed = (unsigned int) clock();*/
    if (userSeed > 0)
        seed = userSeed;
	
    originalSeed = seed;
    for (i=0; i<100; i++)
        RandomUniform(&seed);

    if (noisy > 0)
		{
		PrintCommandLine (stderr, argc, argv);
		}

     /* Check the distribution of the mutations */
    #ifdef CHECK_MUT_DISTRIBUTION
        MutCount = (int *) calloc(2*(numCells+1), sizeof(int));
        dataSetsWithSNVs = 0;
        meansumPos = 0;
     #endif
	
	/* set names for taxa */
	if (doTumorNames == YES)
		{
		strcpy(inCellName, "tumcell");
		strcpy(outCellName, "healthycell");
		strcpy(inRootCellName, "tumoralroot");
		strcpy(outRootCellName, "healthyroot");
		stringPrecision = 12;
		}
	else
		{
		strcpy(inCellName, "cell");
		strcpy(outCellName, "outgcell");
		strcpy(inRootCellName, "ingrroot");
		strcpy(outRootCellName, "outgroot");
		stringPrecision = 9;
		}
	
	if (doPrintIUPAChaplotypes == YES)
		stringPrecision += 1;
	
    /* set file dirs and names */
	if (strlen(resultsDir) == 0)
		strcpy(resultsDir, "Results");
	// clean previous results directory
	RemoveDir (resultsDir);
	// create new results directory
	MAKE_DIR (resultsDir);

   	strcpy(treeDir, "trees_dir");
   	strcpy(timesDir, "times_dir");
	strcpy(SNVgenotypesDir, "snv_genotypes_dir");
    strcpy(SNVhaplotypesDir, "snv_haplotypes_dir");
    strcpy(trueHaplotypesDir, "true_haplotypes_dir");
    strcpy(fullHaplotypesDir, "full_haplotypes_dir");
	strcpy(MLhaplotypesDir, "ML_haplotypes_dir");
    strcpy(fullGenotypesDir, "full_genotypes_dir");
    strcpy(CATGdir, "catg_dir");
	strcpy(PILEUPdir, "pileup_dir");
	strcpy(bulkDir, "bulk_dir");
	strcpy(VCFdir, "vcf_dir");
	
	strcpy(SNVgenotypesFile, "snv_gen");
    strcpy(SNVhaplotypesFile, "snv_hap");
    strcpy(trueHaplotypesFile, "true_hap");
    strcpy(fullHaplotypesFile, "full_hap");
	strcpy(MLhaplotypesFile, "ML_hap");
    strcpy(fullGenotypesFile, "full_gen");
    strcpy(treeFile, "trees");
    strcpy(timesFile, "times");
	if (strlen(userTreeFile) == 0)
		strcpy(userTreeFile, "usertree");
	if (strlen(userGenomeFile) == 0)
		strcpy(userGenomeFile, "usergenome");
    strcpy(CATGfile, "catg");
	strcpy(PILEUPfile, "pileup");
	strcpy(bulkFile, "bulk");
	strcpy(VCFfile, "vcf");
    strcpy(logFile, "log");
	strcpy(settingsFile, "log");
 
	if (doPrintSeparateReplicates == NO)
		PrepareGlobalFiles(argc, argv);

	/* Start printing info to log */
	snprintf (File, sizeof(File),"%s/%s", resultsDir, settingsFile);
	if ((fpLog = fopen(File, "w")) == NULL)
		{
		fprintf (stderr, "Can't open \"%s\"\n", File);
		exit(-1);
		}

	PrintHeader(fpLog);
	PrintDate(fpLog);
	snprintf (File, sizeof(File),"%s/%s", resultsDir, settingsFile);
	PrintCommandLine (fpLog, argc, argv);
 
	#ifdef COUNT_ML_GENOTYPE_ERRORS
		snprintf (File, sizeof(File),"%s/%s", resultsDir, "sims");
		 if ((fpSims= fopen(File, "w")) == NULL)
			 {
			 fprintf (stderr, "Can't open \"%s\"\n", File);
			 exit(-1);
			 }
	fprintf(fpSims,"CountCalls, ProbCalls, CountBadCalls, ProbBadCalls, numSNVs\n");
	#endif
	
	/* --summary: one line of metrics per replicate in <results>/summary.tsv */
	if (doSummary == YES)
		{
		snprintf (File, sizeof(File), "%s/%s", resultsDir, "summary.tsv");
		if ((fpSummary = fopen(File, "w")) == NULL)
			{
			fprintf (stderr, "Can't open \"%s\"\n", File);
			exit(-1);
			}
		fprintf (fpSummary, "replicate\ttmrca\tmutations\tdeletions\tcnloh\tsnv_true\tsnv_observed\tsnv_ratio\tcalled_fraction\tgenotype_error\tsnv_precision\tsnv_recall\thet_error\thom_error\thet_to_hom\tdoublet_error\tsinglet_error\tado_simulated\n");
		}

	/* Read user treefile */
	if (doUserTree == YES)
		{
		if (noisy > 2)
			fprintf (stderr, "\n>> Reading user tree ...");
		ReadUserTree(fpUserTree);
		}

	/* Read user genomefile  */
	if (doUserGenome == YES)
		{
		if (noisy > 2)
			fprintf (stderr, "\n>> Reading user genome ...");
		ReadUserGenome(fpUserGenome);
		}

	OUTGROUP_ROOT = 2 * numCells;
	INGROUP_ROOT = (2 * numCells) - 1;

	if (doGeneticSignatures == YES && doSimulateData == YES)
		PrepareGeneticSignatures();

	if (noisy > 1)
		fprintf(stderr, "\n");
	
	for (dataSetNum=0; dataSetNum<numDataSets; dataSetNum++)
		RunReplicate (dataSetNum, &seed);

	ComputeSummaryStatistics ();

	/* close file pointers */
	if (doPrintSeparateReplicates == NO)
		{
		if (doPrintSNVgenotypes == YES)
			fclose(fpSNVgenotypes);
		if (doPrintFullGenotypes == YES)
			fclose(fpFullGenotypes);
		if (doPrintSNVhaplotypes == YES)
			fclose(fpSNVhaplotypes);
		if (doPrintFullHaplotypes == YES)
			fclose(fpFullHaplotypes);
		if (doPrintTrueHaplotypes == YES)
			fclose(fpTrueHaplotypes);
		if (doNGS == YES && doPrintVCF == YES)
			CloseVCF (fpVCF);
		if (doPrintCATG == YES)
			fclose(fpCATG);
		if (doPrintPILEUP == YES)
			fclose(fpPILEUP);
		if (doPrintBulk == YES)
			fclose(fpBulk);
		if (doPrintTree == YES)
			fclose(fpTrees);
		if (doPrintTimes == YES)
			fclose(fpTimes);
		}
	
    if (noisy > 0)
        {
        /* if (noisy < 2)
         PrintHeader (stdout);*/
			
        fprintf (stdout, "\n\n** Simulations finished at ");
        PrintDate (stdout);
        PrintRunInformation (stdout);
			
		/* print run settings to file */
		PrintRunInformation (fpLog);
		}

	#ifdef CHECK_MUT_DISTRIBUTION
	fprintf (stdout,"\n\nChecking distribution of mutations");
	fprintf (stdout,"\n Avg # sites with  0 mutant alelles = %6.2f  (Exp = %6.2f)",(double)MutCount[0]/numDataSets, numSites-expNumMU);
	for (i=1; i<numCells; i++)
		fprintf (stdout,"\n Avg # sites with %2d mutant alelles = %6.2f  (Exp = %6.2f)",i,(double) MutCount[i]/numDataSets, theta/i);
	fprintf (stdout, "\n MeansumPos = %4.2f  (expected (L/2) = %4.2f)\n\n", meansumPos / (double) dataSetsWithSNVs, numSites/2.0);
    #endif
	
	/* elapsed times, and a record of the run (it needs the command line) */
	{
	struct timespec wallEnd;
	clock_gettime (CLOCK_MONOTONIC, &wallEnd);
	wallSecs = (double)(wallEnd.tv_sec - wallStart.tv_sec) + 1e-9 * (double)(wallEnd.tv_nsec - wallStart.tv_nsec);
	}
	secs = (double)(clock() - start) / CLOCKS_PER_SEC;
	WriteRunRecord (wallSecs, secs);

	/* free pending stuff */
	free (CommandLine);
	if (doUserTree == YES)
		FreeUserTree ();
	free (cellNames);
	free (triMutationsCounter);
	free (signatureWeight);
	free (signatureID);
	if (doUserGenome == YES)
		{
		free (paternalUserGenome);
		free (maternalUserGenome);
		}
	if (doTreeRateSwitches == YES)
		free (treeRateSwitchMultiplier);


	if (fpSummary != NULL)
		fclose (fpSummary);
#ifdef COUNT_ML_GENOTYPE_ERRORS
	fprintf (stdout, "\n\n Mean probability for untrue genotype calls     =   %8.6f", runStats.mlErr.mean);
	fprintf (stdout, "\n Variance probability for untrue genotype calls =   %8.6f", runStats.mlErr.var);
	fprintf (stdout, "\n Mean proportion of SNV genotypes called        =   %8.6f", runStats.called.mean);
	fprintf (stdout, "\n Variance proportion of SNV genotypes called    =   %8.6f\n", runStats.called.var);
	fclose(fpSims);
#endif
	

    if (noisy > 0)
        {
        printf("\n\n_________________________________________________________________");
        printf("\nTime processing: %G seconds (%G seconds of CPU time)\n\n", wallSecs, secs);
        }
	return (0);
    }




/********************* WriteRunRecord  ************************/
/* Writes Results/run.json: everything needed to repeat the run and to know how it was run (version, command line, seed,
   threads, environment variables, genome and read storage modes, times). Failing to write it never stops the run. */
static void JsonString (FILE *fp, const char *s)
	{
	fputc ('"', fp);
	for (; s && *s; s++)
		{
		if (*s == '"' || *s == '\\')
			{ fputc ('\\', fp); fputc (*s, fp); }
		else if ((unsigned char) *s < 0x20)
			fprintf (fp, "\\u%04x", (unsigned char) *s);
		else
			fputc (*s, fp);
		}
	fputc ('"', fp);
	}

void WriteRunRecord (double wallSecs, double cpuSecs)
	{
	char		File[2*MAX_NAME+8], when[64];
	FILE		*fp;
	time_t		now = time (NULL);
	struct tm	tmv;
	int			threads = 1, k;
	const char	*envNames[] = {"CELLCOAL_SPARSE", "CELLCOAL_DENSE_MB", "CELLCOAL_STREAM", "CELLCOAL_STREAM_MB", "CELLCOAL_GL_CACHE_MB", "OMP_NUM_THREADS"};
	int			first = YES;

	snprintf (File, sizeof(File), "%s/%s", resultsDir, "run.json");
	if ((fp = fopen (File, "w")) == NULL)
		return;
#ifdef _OPENMP
	threads = omp_get_max_threads ();
#endif
	gmtime_r (&now, &tmv);
	strftime (when, sizeof(when), "%Y-%m-%dT%H:%M:%SZ", &tmv);
	fprintf (fp, "{\n  \"program\": \"%s\",\n  \"version\": ", PROGRAM_NAME);
	JsonString (fp, VERSION_NUMBER);
	fprintf (fp, ",\n  \"finished\": \"%s\",\n  \"command_line\": ", when);
	JsonString (fp, CommandLine ? CommandLine : "");
	{
	/* the options of the run (a valid --config file), with the seed that was used */
	char	*cmd = (char*) calloc (strlen (CommandLine ? CommandLine : "") + 64, 1), *tok, *save = NULL, *copy = strdup (CommandLine ? CommandLine : "");

	for (tok = strtok_r (copy, " ", &save); tok; tok = strtok_r (NULL, " ", &save))
		{
		if (tok[0] == '-' && tok[1] == '#')
			continue;
		strcat (cmd, tok);
		strcat (cmd, " ");
		}
	snprintf (cmd + strlen (cmd), 48, "-#%ld", originalSeed);
	fprintf (fp, ",\n  \"options\": ");
	WriteOptionsJson (fp, cmd, "  ");
	free (cmd);
	free (copy);
	}
	fprintf (fp, ",\n  \"seed\": %ld,\n  \"replicates\": %d,\n  \"cells\": %d,\n  \"sites\": %d,\n", originalSeed, numDataSets, numCells, numSites);
	fprintf (fp, "  \"threads\": %d,\n  \"results_dir\": ", threads);
	JsonString (fp, resultsDir);
	fprintf (fp, ",\n  \"reads_simulated\": %s,\n  \"genome_storage\": \"%s\",\n  \"reads_mode\": \"%s\",\n", (doNGS == YES) ? "true" : "false",
			 (sparseMode == YES) ? "sparse" : "dense", (doNGS == NO) ? "none" : (streamMode == YES) ? ((readsOnly == YES) ? "streamed, reads only" : "streamed") : "stored");
	fprintf (fp, "  \"environment\": {");
	for (k=0; k<(int)(sizeof(envNames)/sizeof(envNames[0])); k++)
		if (getenv (envNames[k]) != NULL)
			{
			fprintf (fp, "%s\n    \"%s\": ", (first == YES) ? "" : ",", envNames[k]);
			JsonString (fp, getenv (envNames[k]));
			first = NO;
			}
	fprintf (fp, "%s},\n  \"compiler\": ", (first == YES) ? "" : "\n  ");
#ifdef __VERSION__
	JsonString (fp, __VERSION__);
#else
	JsonString (fp, "unknown");
#endif
#ifdef _OPENMP
	fprintf (fp, ",\n  \"openmp\": true");
#else
	fprintf (fp, ",\n  \"openmp\": false");
#endif
	fprintf (fp, ",\n  \"wall_seconds\": %.3f,\n  \"cpu_seconds\": %.3f\n}\n", wallSecs, cpuSecs);
	fclose (fp);
	}

