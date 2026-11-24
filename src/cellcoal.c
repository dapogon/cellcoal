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
				Version 2.x: developed by Claude (Anthropic, as a coding agent) in collaboration with DP
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
- produce read counts (the number of reads of each nucleotide per cell and site) or paired-end reads (BAM, FASTQ)
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
			fprintf (stderr, "\n!!! PARAMETER ERROR: --bcf replaces the VCF file: it cannot be combined with --output vcf-gz or --output no-vcf\n\n");
			exit (-1);
			}
		if (numDataSets > 1 && doPrintSeparateReplicates == NO)
			{
			fprintf (stderr, "\n!!! PARAMETER ERROR: --bcf with several replicates needs --output replicate-folders (one BCF file per replicate)\n\n");
			exit (-1);
			}
		}
	if (ampCellVar > 0.0 && (meanAmplificationError <= 0.0 || ampCellVar >= meanAmplificationError * (1.0 - meanAmplificationError)))
		{
		fprintf (stderr, "\n!!! PARAMETER ERROR: --amp-cell-var needs an amplification error (--amplification error) and a variance < mean*(1-mean) of --amplification error\n\n");
		exit (-1);
		}
	if (covAdoCorr != 0.0 && (alphaCellCoverage <= 0.0 || doADOcell == NO))
		{
		fprintf (stderr, "\n!!! PARAMETER ERROR: --cov-ado-corr needs coverage variation among cells (--depth cells) and ADO variation among cells (--ado cells)\n\n");
		exit (-1);
		}
	if (covDoubletCorr != 0.0 && (alphaCellCoverage <= 0.0 || meanDoubletRate <= 0))
		{
		fprintf (stderr, "\n!!! PARAMETER ERROR: --cov-doublet-corr needs coverage variation among cells (--depth cells) and doublets (--depth doublets)\n\n");
		exit (-1);
		}
	if (adoDoubletCorr != 0.0 && (doADOcell == NO || meanDoubletRate <= 0))
		{
		fprintf (stderr, "\n!!! PARAMETER ERROR: --ado-doublet-corr needs ADO variation among cells (--ado cells) and doublets (--depth doublets)\n\n");
		exit (-1);
		}
	}

/* The three ways CellCoal produces data: genotypes (no -C), read counts (-C) and reads (--bam, --fastq) */
int RunMode (void)
	{
	if (doBam || doFastq)
		return 2;
	return (doNGS == YES) ? 1 : 0;
	}

/* options of a mode that the other modes do not use: an error in reads mode, a warning without -C */
static void CheckModeOptions (void)
	{
	int		reads = (RunMode () == 2), counts = (RunMode () == 1), genotypes = (RunMode () == 0);
	struct { const char *name; int given; int countsOnly; } o[] =
		{
		{"-I (allelic imbalance)", doAllelicImbalance == YES, 1},
		{"-V (coverage overdispersion)", rateVarCoverage == YES, 1},
		{"-O (coverage overdispersion between sites)", alphaSiteCoverage > 0.0, 1},
		{"-R (haploid coverage reduction)", haploidCoverageReduction != 0.5, 1},
		{"-M (minimum pooled VAF)", minVAFvcf > 0.0, 1},
		{"--fast", fastMode != 0, 1},
		{"--min-alt-reads", minAltReads > 0, 1},
		{"--min-alt-vaf", minAltVaf > 0.0, 1},
		{"--coverage", coverageWindow > 0, 1},
		{"-8, -9, -0, -Z, -z, --bcf (outputs of read counts)", doPrintCATG == YES || doPrintPILEUP == YES || doPrintBulk == YES || doBCF == YES || doZipVCF == YES, 1},
		{NULL, 0, 0}
		};
	struct { const char *name; int given; } both[] =
		{
		{"-J (coverage overdispersion between cells)", alphaCellCoverage > 0.0},
		{"-A (amplification error)", meanAmplificationError > 0.0},
		{"-E (sequencing error)", sequencingError > 0.0},
		{"-X (error matrix)", thereIsEij == YES},
		{"--gain-sd/--gain-len", gainSd > 0.0 || gainLen > 0.0},
		{"--amp-signature", ampSignatureFile != NULL},
		{NULL, 0}
		};
	struct { const char *name; int given; } readsOnlyArg[] =
		{
		{"--read-len", ExtraOptionGiven ("--read-len")},
		{"--insert-mean", ExtraOptionGiven ("--insert-mean")},
		{"--insert-sd", ExtraOptionGiven ("--insert-sd")},
		{"--dup-mean", ExtraOptionGiven ("--dup-mean")},
		{"--bulk-depth", ExtraOptionGiven ("--bulk-depth")},
		{"--qual-sd", ExtraOptionGiven ("--qual-sd")},
		{"--qual-slope", ExtraOptionGiven ("--qual-slope")},
		{"--strand-bias", ExtraOptionGiven ("--strand-bias")},
		{"--seq-context", ExtraOptionGiven ("--seq-context")},
		{"--damage", ExtraOptionGiven ("--damage")},
		{"--indel-rate", ExtraOptionGiven ("--indel-rate")},
		{"--indel-ins", ExtraOptionGiven ("--indel-ins")},
		{"--indel-len", ExtraOptionGiven ("--indel-len")},
		{"--indel-max", ExtraOptionGiven ("--indel-max")},
		{"--homopolymer-error", ExtraOptionGiven ("--homopolymer-error")},
		{NULL, 0}
		};
	int		k;

	if (reads)
		for (k=0; o[k].name; k++)
			if (o[k].given)
				{
				fprintf (stderr, "\n!!! PARAMETER ERROR: %s is an option of the read counts (-C) and is not used with --bam/--fastq (reads mode): the coverage varies with --depth cells and --gain-sd\n\n", o[k].name);
				exit (-1);
				}
	if (!reads)
		for (k=0; readsOnlyArg[k].name; k++)
			if (readsOnlyArg[k].given)
				{
				fprintf (stderr, "\n!!! PARAMETER ERROR: %s is an option of the reads mode and needs --bam and/or --fastq\n\n", readsOnlyArg[k].name);
				exit (-1);
				}
	if (genotypes)
		{
		for (k=0; o[k].name; k++)
			if (o[k].given)
				fprintf (stderr, "\nNote: %s is ignored without -C (genotypes mode: no reads are simulated)\n", o[k].name);
		for (k=0; both[k].name; k++)
			if (both[k].given)
				fprintf (stderr, "\nNote: %s is ignored without -C (genotypes mode: no reads are simulated)\n", both[k].name);
		}
	(void) counts;
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

/* "[#n]" separates the replicates in the output files (not in the chromosome in blocks, which is one data set) */
static void RepMarker (FILE *fp, int dataSetNum)
	{
	if (blocksOn == 0)
		fprintf (fp, "[#%d]\n", dataSetNum+1);
	}

/* reads simulated, summarized and printed (VCF, PILEUP, bulk, CATG) in blocks of sites */
static void GenerateStreamedReads (long int *seed, int dataSetNum)
	{
	/* reads are simulated, summarized and printed (VCF, PILEUP) in blocks of sites */
	if (doPrintVCF == YES && doPrintSeparateReplicates == NO)
		RepMarker (fpVCF, dataSetNum);
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
			RepMarker (fpCATG, dataSetNum);
		}
	if (doPrintPILEUP == YES && doPrintSeparateReplicates == NO)
		RepMarker (fpPILEUP, dataSetNum);
	if (doPrintBulk == YES && doPrintSeparateReplicates == NO)
		RepMarker (fpBulk, dataSetNum);
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
	if (gainRate > 0)
		fprintf (stderr, "\nNumber of segment gains        =   %d", numSegGain);
	if (lossRate > 0)
		fprintf (stderr, "\nNumber of segment losses       =   %d (variable sites lost in them: %d)", numSegLoss, snvsLostBySegments);
	
	if (doTreeRateSwitches == YES)
		fprintf(stderr, "\nTree rate switches             =   %d %s", numTreeRateSwitches, treeRateInfo);

	fprintf (stderr, "\n\n");
	}

	if (doPrintSNVgenotypes == YES && numSNVs > 0) /* we only print replicates with variation */
		{
		if (doPrintSeparateReplicates == NO)
			RepMarker (fpSNVgenotypes, dataSetNum);
		PrintSNVGenotypes(fpSNVgenotypes);
		}
	if (doPrintSNVhaplotypes == YES && numSNVs > 0) /* we only print replicates with variation */
		{
		if (doPrintSeparateReplicates == NO)
			RepMarker (fpSNVhaplotypes, dataSetNum);
		PrintSNVHaplotypes(fpSNVhaplotypes, NO);
		}
	if (doPrintFullGenotypes == YES)
		{
		if (doPrintSeparateReplicates == NO)
			RepMarker (fpFullGenotypes, dataSetNum);
		PrintFullGenotypes(fpFullGenotypes);
		}
	if (doPrintFullHaplotypes == YES)
		{
		if (doPrintSeparateReplicates == NO)
			RepMarker (fpFullHaplotypes, dataSetNum);
		PrintFullHaplotypes(fpFullHaplotypes);
		}
	CloseFullStore ();

	if (alphabet == DNA && streamMode == NO && doNGS == YES)
		{
		if (doPrintVCF == YES)
			{
			if (doPrintSeparateReplicates == NO)
				RepMarker (fpVCF, dataSetNum);
			PrintVCF (fpVCF);
			}
		if (doPrintCATG == YES)
			{
			PrintCATG(fpCATG);
			if (doPrintSeparateReplicates == NO)
				RepMarker (fpCATG, dataSetNum);
			}
		if (doPrintPILEUP == YES)
			{
			PrintPILEUP(fpPILEUP);
			if (doPrintSeparateReplicates == NO)
				RepMarker (fpPILEUP, dataSetNum);
			}
		if (doPrintBulk == YES)
			{
			PrintBulk (fpBulk);
			if (doPrintSeparateReplicates == NO)
				RepMarker (fpBulk, dataSetNum);
			}
		}
	}

/* chromosome in blocks: the rows of summary.tsv are kept to write the TOTAL row */
static struct { double cnt[5]; double m[NUM_TRUTH_METRICS]; long long inv[10]; } *blockRows = NULL;
static int numBlockRows = 0;

static void StoreBlockRow (int mu, int del, int cnloh, int snvTrue, int snvObs, const double *metric, const long long *inv)
	{
	int		k;

	blockRows = realloc (blockRows, (size_t) (numBlockRows + 1) * sizeof(*blockRows));
	blockRows[numBlockRows].cnt[0] = mu; blockRows[numBlockRows].cnt[1] = del; blockRows[numBlockRows].cnt[2] = cnloh;
	blockRows[numBlockRows].cnt[3] = snvTrue; blockRows[numBlockRows].cnt[4] = snvObs;
	for (k=0; k<NUM_TRUTH_METRICS; k++)
		blockRows[numBlockRows].m[k] = metric[k];
	for (k=0; k<10; k++)
		blockRows[numBlockRows].inv[k] = inv[k];
	numBlockRows++;
	}

/* TOTAL row of the chromosome: counts are added, rates are averaged (precision, recall and singleton recall weighted by their denominator) */
static void WriteBlocksTotal (void)
	{
	int		k, b;
	double	sum[5] = {0,0,0,0,0}, single = 0.0;
	long long	inv[10] = {0,0,0,0,0,0,0,0,0,0};
	static const int meanIdx[] = {0, 1, 4, 5, 6, 8, 9};

	for (b=0; b<numBlockRows; b++)
		{
		for (k=0; k<5; k++)
			sum[k] += blockRows[b].cnt[k];
		for (k=0; k<10; k++)
			inv[k] += blockRows[b].inv[k];
		if (blockRows[b].m[10] > 0.0)
			single += blockRows[b].m[10];
		}
	fprintf (fpSummary, "TOTAL\tNA\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\tNA", sum[0], sum[1], sum[2], sum[3], sum[4]);
	for (k=0; k<NUM_TRUTH_METRICS; k++)
		{
		double	num = 0.0, den = 0.0;
		int		w = (k == 2) ? 4 : (k == 3) ? 3 : -1, isMean = 0, j;

		for (j=0; j<7; j++)
			if (meanIdx[j] == k)
				isMean = 1;
		if (k == 7)
			{ fprintf (fpSummary, "\tNA"); continue; }
		if (k == 10)
			{ fprintf (fpSummary, "\t%.0f", single); continue; }
		for (b=0; b<numBlockRows; b++)
			{
			double v = blockRows[b].m[k], wt = 1.0;

			if (v < 0.0)
				continue;
			if (k == 11)
				wt = blockRows[b].m[10];
			else if (w >= 0)
				wt = blockRows[b].cnt[w];
			if (!isMean && wt <= 0.0)
				continue;
			num += wt * v;
			den += wt;
			}
		if (den > 0.0)
			fprintf (fpSummary, "\t%.6f", num / den);
		else
			fprintf (fpSummary, "\tNA");
		}
	for (k=0; k<10; k++)
		fprintf (fpSummary, "\t%lld", inv[k]);
	if (lossRate > 0)
		fprintf (fpSummary, "\tNA\tNA");
	if (gainRate > 0)
		fprintf (fpSummary, "\tNA");	/* the segment losses are not summed over the blocks */
	fprintf (fpSummary, "\n");
	free (blockRows);
	blockRows = NULL;
	numBlockRows = 0;
	}

/* one replicate: genealogy, mutations, deletions/cnLOH, ADO, errors, reads and all the per-replicate outputs */
/* ---- the phases of a replicate ---- */

static void ReportProgress (int dataSetNum)
	{
	fprintf (stdout, "\rReplicate #%3d/%d ", dataSetNum+1, numDataSets);
	fflush (stdout);
	if (!isatty (fileno (stdout)))		/* if stdout is redirected elsewhere, print also to stderr */
		{
		fprintf (stderr, "\rreplicate #%3d/%d", dataSetNum+1, numDataSets);
		fflush (stderr);
		}
	}

/* coalescent tree (a chromosome in blocks: only for the first block, the others use the same tree) */
static void MakeTreeOfReplicate (int dataSetNum, long int *seed)
	{
	int		sameTree = (blocksOn && dataSetNum > 0);

	if (doUserTree == NO && !sameTree)
		{
		if (noisy > 2)
			fprintf (stderr, "\n>> Starting coalescent tree ...");
		MakeCoalescenceTree (numCells, N, seed);
		if (noisy > 2)
			fprintf (stderr, "\n>> Finishing coalescent tree ... DONE");

		/* make the ingroup cell tree non-clock if needed */
		if (rateVarAmongLineages == YES)
			MakeTreeNonClock (coalTreeMRCA, seed);
		else if (doTreeRateSwitches == YES)
			{
			treeRateInfo = (char *) calloc (MAX_LINE, sizeof(char));
			TreeRateSwitch (coalTreeMRCA, seed);
			}
		runStats.ca.sum += numCA;
		AccumAdd (&runStats.tmrca, TMRCA);
		}
	totalTreeLength = SumBranches (healthyRoot);
	if (doPrintTree == YES && !sameTree)
		PrintTree (healthyRoot, fpTrees);
	if (doPrintTimes == YES && !sameTree)
		PrintTimes (0, fpTimes);
	}

/* memory of the genotypes, of the cells and of the sites of a replicate */
static void AllocateReplicate (void)
	{
	int		i, j;

	DecideSparseMode ();
	data = (Geno ***) Checked (calloc (ploidy, sizeof(Geno **)), "Could not allocate the data structure\n");
	for (i=0; i<ploidy; i++)
		{
		data[i] = (Geno **) Checked (calloc (2*numCells+1, sizeof(Geno *)), "Could not allocate the data[] structure\n");
		for (j=0; j<2*numCells+1; j++)
			{
			if (sparseMode == YES && j != OUTGROUP_ROOT)
				continue;	/* in sparse mode only the genome of the root is stored densely */
			data[i][j] = (Geno *) Checked (calloc (numSites, sizeof(Geno)), "Could not allocate the data[][] structure\n");
			}
		}
	if (sparseMode == YES)
		{
		numNodes = 2*numCells+1;
		sparseCols = (Geno **) Checked (calloc (numSites, sizeof(Geno *)), "Could not allocate the sparseCols structure\n");
		}
	AllocateCellStructure ();
	allSites = (SiteStr*) Checked (calloc (numSites, sizeof(SiteStr)), "Could not allocate the allSites structure\n");
	for (i=0; i<numSites; i++)
		allSites[i].alternateAlleles = (int *) Checked (calloc (4, sizeof(int)), "Could not allocate the allSites[i].alternateAlleles[] structure\n");
	/* the arrays below keep the index for different types of sites */
	SNVsites = (int*) Checked (calloc (numSites, sizeof(int)), "Could not allocate the SNVsites structure\n");
	variantSites = (int*) Checked (calloc (numSites, sizeof(int)), "Could not allocate the variantSites structure\n");
	DefaultModelSites = (int*) Checked (calloc (numSites, sizeof(int)), "Could not allocate the DefaultModelSites structure\n");
	AltModelSites = (int*) Checked (calloc (numSites, sizeof(int)), "Could not allocate the AltModelSites structure\n");
	}

/* one process (mutations, CN-LOH, deletions) on the maternal and then on the paternal genome */
static void EvolveBothGenomes (const char *what, void (*Evolve) (TreeNode *, int, long int *), long int *seed)
	{
	int		genome;

	for (genome=MATERNAL; genome<=PATERNAL; genome++)
		{
		if (noisy > 2)
			fprintf (stderr, "\n>> Evolving %s %s ... ", genome == MATERNAL ? "maternal" : "paternal", what);
		Evolve (healthyRoot, genome, seed);
		if (noisy > 2)
			fprintf (stderr, genome == MATERNAL ? "DONE" : "DONE\n\n");
		}
	}

/* the true genomes of all the cells: mutations, germline SNPs, CN-LOH and deletions */
static void SimulateGenomes (int dataSetNum, long int *seed)
	{
	int		i, j;

	numISMmutations = 0;
	numISMdeletions = 0;
	numISMCNLOH = 0;
	if (sparseMode == YES)
		{
		ResetObservationLayer ();
		ResetSparseColumns ();
		}
	InitializeGenomes (healthyRoot, seed);		/* ancestral states of all the genomes */
	if (SNPrate > 0)
		{
		AddGermlineVariation (seed);
		InheritRootGenome (healthyRoot);		/* the germline SNPs were added to the root genome only */
		}
	snvsLostBySegments = 0;
	ResetSegments ();
	if (lossRate > 0 || gainRate > 0)
		{
		/* copy-number events: SNVs, CNLOH, deletions and segment losses and gains, branch by branch in the order of the tree */
		SimulateCopyNumberGenomes (healthyRoot, seed);
		AccumAdd (&runStats.mu, numMU);
		if (CNLOHrate > 0)
			AccumAdd (&runStats.cnloh, numCNLOH);
		if (deletionRate > 0)
			AccumAdd (&runStats.del, numDEL);
		WriteCnvEvents (dataSetNum);
		WriteCnvSegments (dataSetNum);
		}
	else
		{
		EvolveBothGenomes ("genome", EvolveSitesOnTree, seed);
		AccumAdd (&runStats.mu, numMU);
		if (CNLOHrate > 0)
			{
			EvolveBothGenomes ("CN_LOH", EvolveCNLOHonTree, seed);
			AccumAdd (&runStats.cnloh, numCNLOH);
			}
		if (deletionRate > 0)
			{
			EvolveBothGenomes ("deletions", EvolveDeletionsOnTree, seed);
			AccumAdd (&runStats.del, numDEL);
			}
		}

	/* count how many SNVs we observe before ADO and genotype errors, but after CNLOH and deletion */
	numSNVs = CountTrueVariants ();
	numSNVsTrue = numSNVs;
	if (lossRate > 0 || gainRate > 0)
		snvsLostBySegments = (numPlannedSNV > numSNVs) ? numPlannedSNV - numSNVs : 0;		/* the mutations that were applied and are no longer variable: they were lost with a segment */

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
			RepMarker (fpTrueHaplotypes, dataSetNum);
		PrintTrueFullHaplotypes (fpTrueHaplotypes);
		}
	if (doBam || doFastq)
		WriteTruthVcf (dataSetNum);
	}

/* what the sequencing sees: allelic dropout, amplification gain, genotyping errors, and the reads or read counts */
static void ObserveGenomes (int dataSetNum, long int *seed)
	{
	adoReplicate = dataSetNum;
	if (fixedADOrate > 0 || doADOcell == YES || doADOsite == YES)
		{
		if (sparseMode == YES)
			SparseAllelicDropout (seed);
		else
			AllelicDropout (seed);
		}
	if (gainSd > 0.0)			/* amplification gain along the genome (--gain-sd, --gain-len) */
		BuildGain (seed);
	if (meanGenotypingError > 0)
		{
		if (sparseMode == YES)
			SparseGenotypeError (seed);
		else
			GenotypeError (seed);
		}
	if (sparseMode == YES)
		obsActive = YES;		/* from now on the cell genotypes are observed with ADO and errors */

	if (coverageWindow > 0)		/* the windows of --coverage start empty in every replicate */
		{
		covNumWindows = (numSites + coverageWindow - 1) / coverageWindow;
		free (covSum);
		free (covSites);
		covSum = (long long*) Checked (calloc ((size_t) covNumWindows * numCells, sizeof(long long)), "Could not allocate the coverage windows\n");
		covSites = (long long*) Checked (calloc ((size_t) covNumWindows, sizeof(long long)), "Could not allocate the coverage windows\n");
		}

	/* generate the reads or the read counts, and count the SNVs that we observe now */
	if (doBam || doFastq)
		{
		if (meanDoubletRate > 0.0 && !(blocksOn && blockIndex > 0))
			SetupDoublets (seed, 1);		/* the cells and donors of the doublets (their reads are made in GenerateReadLevel) */
		SetupCellCoverage (seed);
		GenerateReadLevel (seed, dataSetNum);
		numSNVs = numSNVsTrue;	/* the reads are written to files; no genotype calling here */
		}
	else if (streamMode == YES)
		GenerateStreamedReads (seed, dataSetNum);	/* the SNVs are counted by StreamReads */
	else if (doNGS == YES)
		{
		int		j;

		GenerateReadCounts (seed);
		numSNVs = CountAllelesInMLGenotypes ();
		if (coverageWindow > 0)
			for (j=0; j<numSites; j++)
				CoverageAddSite (j);
		}
	else
		numSNVs = CountAllelesInObservedGenotypes ();
	}

/* a row of summary.tsv (and of the table of the blocks) */
static void WriteSummaryRow (int dataSetNum)
	{
	double		metric[NUM_TRUTH_METRICS];
	long long	inv[10];
	int			k;

	TruthMetrics (metric);
	InvariantGenotypeCounts (inv);
	if (blocksOn)
		StoreBlockRow (numMU, numDEL, numCNLOH, numSNVsTrue, numSNVs, metric, inv);
	fprintf (fpSummary, "%d\t%.6f\t%d\t%d\t%d\t%d\t%d\t%.6f", dataSetNum+1, TMRCA, numMU, numDEL, numCNLOH, numSNVsTrue, numSNVs,
		numSNVsTrue > 0 ? numSNVs / (double) numSNVsTrue : 0.0);
	for (k=0; k<NUM_TRUTH_METRICS; k++)	/* NA when there are no reads or the metric is undefined */
		{
		if (metric[k] < 0.0)
			fprintf (fpSummary, "\tNA");
		else
			fprintf (fpSummary, "\t%.6f", metric[k]);
		}
	for (k=0; k<10; k++)
		fprintf (fpSummary, "\t%lld", inv[k]);
	if (lossRate > 0)
		fprintf (fpSummary, "\t%d\t%d", numSegLoss, snvsLostBySegments);
	if (gainRate > 0)
		fprintf (fpSummary, "\t%d", numSegGain);
	fprintf (fpSummary, "\n");
	fflush (fpSummary);
	}

/* closes the files of a replicate that has its own folder, and frees what the replicate allocated */
static void FreeReplicate (int dataSetNum)
	{
	int		i, j;

	if (doUserTree == NO && !(blocksOn && dataSetNum < numDataSets - 1))	/* the tree of a chromosome in blocks is kept until the last block */
		{
		for (i=0; i<numNodes; i++)
			free (treeNodes[i].name);
		free (treeNodes);
		}
	if (doPrintSeparateReplicates == YES)
		{
		if (doPrintTree == YES)
			fclose (fpTrees);
		if (doPrintTimes == YES)
			fclose (fpTimes);
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
		free (allSites[i].alternateAlleles);
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
	/* the per-cell structures (allocated at every replicate, with or without NGS) */
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
	}

static void RunReplicate (int dataSetNum, long int *seed)
	{
	int		i;
	double	keepTMRCA = TMRCA;

	numSegLoss = 0;
	numCA = numMU = numDEL = numCNLOH = numProposedMU = TMRCA = numSNVmaternal = 0;
	if (blocksOn && dataSetNum > 0)
		TMRCA = keepTMRCA;		/* chromosome in blocks: the tree of the first block */
	if (doPrintSeparateReplicates == YES)
		PrepareSeparateFiles (dataSetNum);
	if (noisy == 1)
		ReportProgress (dataSetNum);

	MakeTreeOfReplicate (dataSetNum, seed);
	AllocateReplicate ();
	SimulateGenomes (dataSetNum, seed);
	ObserveGenomes (dataSetNum, seed);
	AccumAdd (&runStats.snv, numSNVs);
	if (fpSummary != NULL)
		WriteSummaryRow (dataSetNum);
	if (coverageWindow > 0)
		CoverageWrite (dataSetNum+1);
	if (numSNVs == 0)
		zeroSNVs++;
	PrintReplicateOutputs (dataSetNum);

	if (noisy > 2)
		for (i=0; i<numSites; i++)
			PrintSiteInfo (stderr, i);
	FreeReplicate (dataSetNum);
	}

static void SetCoalescentExpectations (void)
	{
	int		i;
	double	a = 0, a2 = 0;

	numNodes = 2 * numCells + 1;
	maxTreeRateSwitches = numCells;
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
	}

static void InitCumulativeTables (void)
	{
	int		i, j;

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
			for (j=0; j<4; j++)
				ampSubRow[i][j] = Eij[i][j] / cumEij[i][3];
			}
		}
	}

static void SetNamesAndPaths (void)
	{
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
 
	}

static void CheckReadOptions (void)
	{
	if (gainSd > 0.0 || gainLen > 0.0)
		{
		if (!(gainSd > 0.0 && gainLen > 0.0))
			{
			fprintf (stderr, "\nERROR: --gain-sd and --gain-len go together\n");
			exit (-1);
			}
		if (doNGS == NO)
			{
			fprintf (stderr, "\nERROR: --gain-sd needs read counts (-C)\n");
			exit (-1);
			}
		}

	if (ampSignatureFile != NULL)
		{
		if (doNGS == NO || meanAmplificationError <= 0.0)
			{
			fprintf (stderr, "\nERROR: --amp-signature needs read counts (-C) and an amplification error (--amplification error)\n");
			exit (-1);
			}
		if (alphabet == BINARY)
			{
			fprintf (stderr, "\nERROR: --amp-signature needs the DNA alphabet\n");
			exit (-1);
			}
		SetupAmpSignature (ampSignatureFile);
		}

	if (adoBlockLen > 0.0 && (doADOsite == YES || (fixedADOrate <= 0 && doADOcell == NO)))
		{
		fprintf (stderr, "\nERROR: --ado-block needs allelic dropout (--ado mean or --ado cells) and cannot be combined with ADO among sites (--ado sites)\n");
		exit (-1);
		}

	/* read-level simulator: --bam, --fastq */
	if (doBam || doFastq)
		{
		if (doNGS == NO || alphabet == BINARY || ploidy != 2)
			{
			fprintf (stderr, "\nERROR: --bam/--fastq need the mean coverage (-C) and diploid DNA\n");
			exit (-1);
			}
		if ((fixedADOrate > 0 || doADOcell == YES || doADOsite == YES) && adoBlockLen <= 0.0)
			{
			fprintf (stderr, "\nERROR: in reads mode the allelic dropout has to come in segments (--ado-block L, e.g. 3000): with dropout independent at every site almost no fragment of %d bp survives\n", 2 * readLen);
			exit (-1);
			}
		if (readLen < 20 || readLen > 400 || insertMean < readLen || insertSd < 0.0 || dupMean < 0.0)
			{
			fprintf (stderr, "\nERROR: bad --read-len/--insert-mean/--insert-sd/--dup-mean\n");
			exit (-1);
			}
		if (qualSd < 0.0 || qualSd > 3.0 || qualSlope < -6.0 || qualSlope > 6.0 || strandBias < 0.0 || strandBias > 1.0)
			{
			fprintf (stderr, "\nERROR: bad --qual-sd (0 to 3), --qual-slope (-6 to 6) or --strand-bias (0 to 1)\n");
			exit (-1);
			}
		if (strlen (seqContextFile) > 0 && sequencingError <= 0.0)
			{
			fprintf (stderr, "\nERROR: --seq-context needs a sequencing error (-E) greater than 0\n");
			exit (-1);
			}
		if (indelRate < 0.0 || indelIns < 0.0 || indelIns > 1.0 || indelLen < 1.0 || indelMax < 1 || indelMax > 32)
			{
			fprintf (stderr, "\nERROR: bad --indel-rate (0 or more), --indel-ins (0 to 1), --indel-len (1 or more) or --indel-max (1 to 32)\n");
			exit (-1);
			}
		if (homoError < 0.0 || homoError > 0.1)
			{
			fprintf (stderr, "\nERROR: bad --homopolymer-error (0 to 0.1)\n");
			exit (-1);
			}
		/* the VCF of read counts and genotype likelihoods does not exist in this mode */
		doPrintVCF = NO;
		}
	if (compressLevel < 1 || compressLevel > 9)
		{
		fprintf (stderr, "\nERROR: --compress-level has to be between 1 and 9\n");
		exit (-1);
		}
	}

static void CloseGlobalFiles (void)
	{
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
	
	}

int main (int argc, char **argv)
    {
    FILE 		*fp;
    int         i, dataSetNum;
    long int	seed;
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
	argc = JoinSplitDecimalArgs (argc, argv);	/* -E0.001 split by PowerShell */
	argv = ExpandConfigArgs (&argc, argv);	/* --config file.json|yaml */
	ReadParametersFromCommandLine (argc, argv);
	readingParameterFile = (argc < 2 || doSpecificParameterFile == YES) ? YES : NO;
	if (readingParameterFile == YES)
		{
		if ((fp = freopen (parameterFile, "r", stdin)) == NULL)
			{
			fprintf (stderr, "\nERROR: Could not open parameter file: %s\n", parameterFile);
			PrintUsage (stderr);
			exit (-1);
			}
		ReadParametersFromFile ();
		}
	CheckOptionCombinations ();
	CheckModeOptions ();
   
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
	
	/* with a user tree the number of cells is the number of leaves of the tree (-s is ignored) */
	if (doUserTree == YES)
		{
		int nLeaves = CountUserTreeLeaves (strlen(userTreeFile) == 0 ? "usertree" : userTreeFile);
		if (nLeaves > 0)
			{
			if (nLeaves != numCells)
				{
				ignoredNumCells = numCells;
				if (noisy > 0)
					fprintf (stderr, "\nNOTE: -s%d ignored: the user tree has %d leaves (cells)\n", numCells, nLeaves);
				}
			numCells = nLeaves;
			}
		}

	DemesSetup ();		/* structured coalescent: the number of cells is that of the demes */

	/* tumor only (--no-outgroup): the outgroup stays in the simulation with a branch of length 0 (= the healthy root genome, 0000) but is not in the outputs */
	if (numOutgroup == 0)
		healthyTipBranchLengthRatio = 0.0;

	SetCoalescentExpectations ();

	memset (&runStats, 0, sizeof(runStats));
	zeroSNVs = 0;
		
    altModelMutationRate = mutationRate*nonISMRelMutRate;
    if (checkOnly == YES)
        {
        PrintCheck (stdout);
        exit (0);
        }
	
	InitCumulativeTables ();

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
        {
        RandomUniform(&seed);
        }

	/* single-nucleotide deletions (-d) are a model of isolated sites: with a chromosome (blocks, real reference) or reads, deletions are small indels */
	if (deletionRate > 0.0 && (blocksN > 0 || strlen (fastaFile) > 0 || doBam == YES || doFastq == YES))
		{
		fprintf (stderr, "\nERROR: -d (single-nucleotide deletions of the sites) cannot be used with a chromosome in blocks, a reference (--fasta) or reads.\n");
		fprintf (stderr, "Use --indel-rate (with --indel-ins, --indel-len, --indel-max) to simulate small insertions and deletions along the chromosome.\n");
		exit (-1);
		}

	/* regional mutation rates need a chromosome with coordinates: a simulated one in blocks or a real one (--blocks, --fasta) */
	if (strlen (regionRatesFile) > 0)
		{
		if (blocksN == 0 && strlen (fastaFile) == 0)
			{
			fprintf (stderr, "\nERROR: --region-rates needs a chromosome with coordinates: use --blocks (simulated chromosome) or --fasta (real chromosome).\n");
			fprintf (stderr, "Without them the sites of the simulation have no position, so there are no regions.\n");
			exit (-1);
			}
		if (doGeneticSignatures == YES)
			{
			fprintf (stderr, "\nERROR: --region-rates cannot be combined with mutational signatures\n");
			exit (-1);
			}
		regionRatesOn = YES;
		rateVarAmongSites = YES;
		}

	/* chromosome in blocks (--blocks, --fasta): the reference and the number of blocks */
	if (blocksN > 0 || strlen (fastaFile) > 0)
		{
		if (doUserGenome == YES || (SNPrate > 0 && (strlen (germlineVcf) > 0 || strlen (germlinePop) > 0)))
			{
			fprintf (stderr, "\nERROR: --blocks/--fasta cannot be combined with --mutation genome, and --germline-vcf with --mutation germline-rate\n");
			exit (-1);
			}
		BlocksSetup ();
		if (strlen (fastaFile) > 0)
			{
			doUserGenome = YES;
			snprintf (userGenomeFile, sizeof(userGenomeFile), "%s", fastaFile);
			}
		if (adoSeedOpt == 0)
			adoSeedOpt = originalSeed + 1;
		}

    if (noisy > 0)
		{
		PrintCommandLine (stderr);
		}

     /* Check the distribution of the mutations */
	
	SetNamesAndPaths ();

	if (doPrintSeparateReplicates == NO)
		PrepareGlobalFiles();

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
	PrintCommandLine (fpLog);
 
	
	/* --summary: one line of metrics per replicate in <results>/summary.tsv */
	if (doSummary == YES)
		{
		snprintf (File, sizeof(File), "%s/%s", resultsDir, "summary.tsv");
		if ((fpSummary = fopen(File, "w")) == NULL)
			{
			fprintf (stderr, "Can't open \"%s\"\n", File);
			exit(-1);
			}
		fprintf (fpSummary, "replicate\ttmrca\tmutations\tdeletions\tcnloh\tsnv_true\tsnv_observed\tsnv_ratio\tcalled_fraction\tgenotype_error\tsnv_precision\tsnv_recall\thet_error\thom_error\thet_to_hom\tdoublet_error\tsinglet_error\tado_simulated\tsingleton_sites\tsingleton_recall\tinv_AA\tinv_AC\tinv_AG\tinv_AT\tinv_CC\tinv_CG\tinv_CT\tinv_GG\tinv_GT\tinv_TT%s%s\n", lossRate > 0 ? "\tsegment_losses\tsnv_lost_segments" : "", gainRate > 0 ? "\tsegment_gains" : "");
		}

	/* Read user treefile */
	if (doUserTree == YES)
		{
		if (noisy > 2)
			fprintf (stderr, "\n>> Reading user tree ...");
		ReadUserTree();
		}

	/* Read user genomefile  */
	if (doUserGenome == YES && strlen (fastaFile) == 0)
		{
		if (noisy > 2)
			fprintf (stderr, "\n>> Reading user genome ...");
		ReadUserGenome();
		}

	CheckReadOptions ();

	/* --coverage: depth per cell of the non-variable sites, in <results>/coverage.tsv */
	if (coverageWindow > 0)
		{
		int		ic;

		if (doNGS == NO)
			{
			fprintf (stderr, "\nERROR: --coverage needs read counts (-C)\n");
			exit (-1);
			}
		snprintf (File, sizeof(File), "%s/%s", resultsDir, "coverage.tsv");
		if ((fpCoverage = fopen(File, "w")) == NULL)
			{
			fprintf (stderr, "Can't open \"%s\"\n", File);
			exit(-1);
			}
		fprintf (fpCoverage, "replicate\tstart\tend\tn_sites");
		for (ic=0; ic<numCells; ic++)
			if (doUserTree == NO)
				fprintf (fpCoverage, "\t%s%04d", inCellName, ic+1);
			else
				fprintf (fpCoverage, "\t%s", cellNames[ic]);
		fprintf (fpCoverage, "\n");
		}

	OUTGROUP_ROOT = 2 * numCells;
	INGROUP_ROOT = (2 * numCells) - 1;

	if (doGeneticSignatures == YES && doSimulateData == YES)
		PrepareGeneticSignatures();

	if (noisy > 1)
		fprintf(stderr, "\n");
	
	for (dataSetNum=0; dataSetNum<numDataSets; dataSetNum++)
		{
		if (blocksOn && BlockPrepare (dataSetNum, &seed) == 0)
			continue;	/* a block with no A, C, G, T sites */
		RunReplicate (dataSetNum, &seed);
		}
	if (blocksOn)
		{
		if (fpSummary != NULL)
			WriteBlocksTotal ();
		snprintf (File, sizeof(File), "%s/%s", resultsDir, "blocks_info.txt");
		if ((fp = fopen (File, "w")) != NULL)
			{
			BlocksInfo (fp);
			fclose (fp);
			}
		BlocksFree ();
		FinishReadLevel ();
		}

	ComputeSummaryStatistics ();

	CloseGlobalFiles ();

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

	
	/* elapsed times, and a record of the run (it needs the command line) */
	{
	struct timespec wallEnd;
	clock_gettime (CLOCK_MONOTONIC, &wallEnd);
	wallSecs = (double)(wallEnd.tv_sec - wallStart.tv_sec) + 1e-9 * (double)(wallEnd.tv_nsec - wallStart.tv_nsec);
	}
	secs = (double)(clock() - start) / CLOCKS_PER_SEC;
	WriteDemes ();
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
	if (fpCoverage != NULL)
		fclose (fpCoverage);
	

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
	char		fileArg[2*MAX_NAME+8], when[64];
	FILE		*fp;
	time_t		now = time (NULL);
	struct tm	tmv;
	int			threads = 1, k;
	const char	*envNames[] = {"CELLCOAL_SPARSE", "CELLCOAL_DENSE_MB", "CELLCOAL_STREAM", "CELLCOAL_STREAM_MB", "CELLCOAL_GL_CACHE_MB", "OMP_NUM_THREADS"};
	int			first = YES;

	snprintf (fileArg, sizeof(fileArg), "%s/%s", resultsDir, "run.json");
	if ((fp = fopen (fileArg, "w")) == NULL)
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

