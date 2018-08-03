/******************************************************************************
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
//  Created by David Posada on 01/06/16.
//  Copyright (c) 2016-2018 David Posada. All rights reserved.

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
- simulate reads for and print all sites in the VCF
- implemented trinucleotide mutational signatures (alone or mixture)
- implemented user genome
- implemented copy-neutral LOH
- implemented doublets at the NGS reads level
- added option for allelic imbalance
- print ML haplotypes
 
[TODOs]
//TODO: at some point move data  structure into cell structure
- dated tips (maybe)
- population structure (maybe)
- change random number generator ?


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
#include "definitions.h"

int main (int argc, char **argv)
    {
    FILE 		*fp;
    int			i, j, k;
    int         dataSetNum;
    long int	seed;
    double		a, a2;
    float 		start, secs;
		
    /* Default settings */
    numDataSets = 10;			/* the number of samples to simulate */
    numCells = 8;				/* number of cells in each data set */
    ploidy = 2;                 /* we assume diploid genomes */
    numSites = 200;				/* number of sites (markers, loci) to simulate = length of the chromosomes */
    N = 1000;					/* effective population size */
    numPeriods = 0;				/* number of distinct demographic periods */
    doDemographics = NO;		/* whether to implement demographics */
    doExponential = NO;			/* whether to do exponential growth */
    growthRate = 0;				/* rate for the exponential population growth */
    mutationRate = 1.0e-7;		/* nucleotide mutation rate per site per generation */
	rateVarAmongLineages = NO;	/* modify rate variation among branches  (to non-clock) */
    alphabet = BINARY;          /* alphabet 0/1 or DNA") */
    altModel = 0;				/* by default the alternative model will be ISM haploid */
    propAltModelSites = 0;		/* proportion of sites that will mutate according to alternative model */
    nonISMRelMutRate = 1.0;		/* relative rate alternative/default model for sites */
    equalBaseFreq = YES;        /* DNA base frequencies */
	freq[0] = freq[1] = freq[2] = freq[3] = 0.25;
	titv = 0.5;					/* transition/transversion rate ratio */
	thereIsMij = NO;			/* mutation rate matrix*/
	Mij[0][0] = Mij[1][1] = Mij[2][2] = Mij[3][3] = 0;  /* mutation probabilities */
	Mij[0][1] = Mij[0][2] = Mij[0][3] = 1.0/3;
	Mij[1][0] = Mij[1][2] = Mij[1][3] = 1.0/3;
	Mij[2][0] = Mij[2][1] = Mij[2][3] = 1.0/3;
	Mij[3][0] = Mij[3][1] = Mij[3][2] = 1.0/3;
	thereIsEij = NO;			/* error rate matrix*/
	Eij[0][0] = Eij[1][1] = Eij[2][2] = Eij[3][3] = 0;  /* sequencing error probabilities */
	Eij[0][1] = Eij[0][2] = Eij[0][3] = 1.0/3;
	Eij[1][0] = Eij[1][2] = Eij[1][3] = 1.0/3;
	Eij[2][0] = Eij[2][1] = Eij[2][3] = 1.0/3;
	Eij[3][0] = Eij[3][1] = Eij[3][2] = 1.0/3;
	doJC = YES;
	doHKY = NO;
	doGTR = NO;
	doGTnR = NO;
	rateVarAmongSites = NO;     	/* rate variation among different sites along the genome */
    alphaSites = infinity;      	/* alpha shape of the gamma distribution for rate variation among sites */
    alphaBranches = infinity;   	/* alpha shape of the gamma distribution for rate variation among lineages */
    alphaCoverage = infinity;   	/* alpha shape of the gamma distribution for coverage */
	doSimulateFixedNumMutations = NO;	/* whether to simulate a fixed number of mutations */
	doUserTree = NO;				/* whether to assume a user tree instead od making the coalescent */
	doUserGenome = NO;				/* whether to use a user genome instead of a simulated one */
    doPrintSNVgenotypes = NO;		/* whether to print SNVs */
    doPrintSNVhaplotypes = NO;  	/* whether to print haplotypes */
    doPrintTrueHaplotypes = NO;  	/* whether to print haplotypes without errors */
	doPrintMLhaplotypes = NO;  		/* whether to print ML haplotypes */
    doPrintFullHaplotypes = NO;		/* whether to print sequences */
    doPrintFullGenotypes = NO;		/* whether to print all genotypes (variable + invariable) */
    doPrintTree = NO;				/* whether to print the coalescent tree */
    doPrintTimes = NO;				/* whether to print coalescent times */
	doPrintAncestors = NO;      	/* whether to print data for ancestral cells */
	doSimulateReadCounts = NO;  	/* do not produce reads by default */
	doPrintCATG = NO;				/* whether to print read counts for SNVs in CATG format*/
	doSimulateData = YES;			/* whether to simulate any data (or just look at the expectations; useful for debugging) */
	doPrintSeparateReplicates = NO; /* whether to put every replica in its own file */
	doPrintIUPAChaplotypes = NO;	/* whether to print IUPAC halotypes */
	doGeneticSignatures = NO;		/* whether to use a genetic signature to model trinucleotide mutations */
	numUserSignatures = 0;			/* by default we do not use a genetic signature */
    healthyTipBranchLength = 0; 	/* length of the branch leading to the healthy cell */
    transformingBranchLength = 0; 	/* length of the transforming branch leading to the healthy ancestral cell */
	coverage = 0;					/* NGS  depth for read counts */
	rateVarCoverage = NO;			/* there is coverage dispersion */
	ADOrate = 0;					/* allelic dropout */
	sequencingError = 0;			/* NGS error rate */
	genotypingError = 0;			/* add errors directly in the genotypes */
	SNPrate = 0.0;					/* germline variation rate for rooth healthy genome */
	meanAmplificationError = 0; 	/* mean of beta distribution for WGA errors */
	varAmplificationError = 0;  	/* variance of beta distribution for WGA errors */
	simulateOnlyTwoTemplates = NO;	/* whether simualate maximum of two templates after single-cell amplification, or there can be all four */
	haploidCoverageReduction = 0.5; /* proportion of reads produced when a single allele is present */
	allelicImbalance = 0.5;			/* proportion of maternal/ paternal reads */
	doubletRate = 0.0;				/* no doublets by default */
	numNodes = 3000;            	/* initial number of nodes allocated to build the coalescent trees */
	seed = time(NULL); 				/* seed for random numbers */
    userSeed = 0;					/* seed entered by the user */
	noisy = 1;						/* level of information to be printed in the screen (see below) */

	/*
	 noisy = 0: does not print anything
	 noisy = 1: + simulation summary,
	 noisy = 2: + replicate information,
	 noisy = 3: + calculation status and event information
	*/
		
	ReadParametersFromCommandLine (argc, argv);
    if (argc < 2)
        {
		readingParameterFile = YES;
        if ((fp = freopen("parameters", "r", stdin)) != NULL)
            ReadParametersFromFile();
        else
            {
            fprintf (stderr, "\nERROR: No parameters specified (use command line or parameter file)");
            PrintUsage();
            exit(-1);
            }
        }
	else
		readingParameterFile = NO;

    start = clock();

    if (noisy > 0)
        {
        PrintHeader(stderr);
        PrintDate(stderr);
 		PrintCommandLine (stderr, argc, argv);
		}
		
    if (doSimulateData == NO)
        {
        doPrintSNVgenotypes = NO;
        doPrintSNVhaplotypes = NO;
		doPrintTrueHaplotypes = NO;
		doPrintFullHaplotypes = NO;
		doPrintMLhaplotypes = NO;
        doPrintFullGenotypes = NO;
		doPrintAncestors = NO;
		doSimulateReadCounts = NO;
		doPrintCATG = NO;
   		}
	
    /* initialize a few variables */
    a = 0;
    a2 = 0;
    for (i=1; i<numCells; i++) /* ignoring outgroup (healthy cell) */
        {
        a += 1.0/(double)i;
        a2 += 1.0/(double)pow(i,2);
        }
    theta = 4.0 * N * mutationRate * 2.0 * numSites * ((1.0-propAltModelSites) + nonISMRelMutRate*propAltModelSites);
    expNumMU = theta * a ;
    expVarNumMU =  theta * a + pow(theta,2) * a2;
    expTMRCA = 2 * (1 - 1.0/numCells) /* 2.0 * N */;
    expVarTMRCA = 0;
    for (i=2; i<=numCells; i++)
        expVarTMRCA += /*4.0 * pow(N,2) */ 4.0 / (pow(i,2) * pow(i-1,2));
    cumNumCA = cumNumMU = cumNumMUSq = cumNumSNVs= cumNumSNVsSq = cumCountMLgenotypeErrors = 0;
    zeroSNVs = cumTMRCA = cumTMRCASq = 0;
		
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

    /* Set seed and spin wheels of pseudorandom number generator */
    /* seed = (unsigned int) clock();*/
    if (userSeed > 0)
        seed = userSeed;
	
	originalSeed = seed;
    for (i=0; i<100; i++)
        RandomUniform(&seed);
		
     /* Check the distribution of the mutations */
    #ifdef CHECK_MUT_DISTRIBUTION
        MutCount = (int *) calloc(2*(numCells+1), sizeof(int));
        dataSetsWithSNVs = 0;
        meansumPos = 0;
     #endif
 
    /* set file dirs and names */
	if (strlen(resultsDir) == 0)
		strcpy(resultsDir, "Results");
	// clean previous results directory
	RemoveDir (resultsDir);
	// create new results directory
	mkdir(resultsDir,S_IRWXU);

   	strcpy(treeDir, "trees_dir");
   	strcpy(timesDir, "times_dir");
	strcpy(SNVgenotypesDir, "snv_genotypes_dir");
    strcpy(SNVhaplotypesDir, "snv_haplotypes_dir");
    strcpy(trueHaplotypesDir, "true_haplotypes_dir");
    strcpy(fullHaplotypesDir, "full_haplotypes_dir");
	strcpy(MLhaplotypesDir, "ML_haplotypes_dir");
    strcpy(fullGenotypesDir, "full_genotypes_dir");
    strcpy(CATGdir, "catg_dir");
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
    strcpy(VCFfile, "vcf");
    strcpy(logFile, "log");
	strcpy(settingsFile, "log");
    #ifdef MYDEBUG
        strcpy(mutationsFile, "mutations");
    #endif

	if (doPrintSeparateReplicates == NO)
		PrepareGlobalFiles(argc, argv);

	/* Start printing info to log */
	sprintf(File,"%s/%s", resultsDir, settingsFile);
	if ((fpLog = fopen(File, "w")) == NULL)
		{
		fprintf (stderr, "Can't open \"%s\"\n", File);
		exit(-1);
		}

	PrintHeader(fpLog);
	PrintDate(fpLog);
	sprintf(File,"%s/%s", resultsDir, settingsFile);
	PrintCommandLine (fpLog, argc, argv);
 
	/* Read user treefile  */
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

	HEALTHY_ROOT = 2 * numCells;
	TUMOR_ROOT = (2 * numCells) - 1;

	if (doGeneticSignatures == YES && doSimulateData == YES)
		PrepareGeneticSignatures();

    for (dataSetNum=0; dataSetNum<numDataSets; dataSetNum++)
        {
        numCA = numMU = numDEL = numProposedMU = TMRCA = numSNVmaternal = 0;

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
			MakeCoalescenceTree (numCells, N, &seed);
			if (noisy > 2)
				fprintf (stderr, "\n>> Finishing coalescent tree ... DONE");
	 
		   /* Make cell tree non-clock if needed */
			if (rateVarAmongLineages == YES)
				MakeTreeNonClock (coalTreeMRCA, &seed);
			
			cumNumCA += numCA;
			cumTMRCA += TMRCA;
			cumTMRCASq += pow(TMRCA,2);
			}
 
        if (doPrintTree == YES)
            PrintTree (healthyRoot, fpTrees);
			
        if (doPrintTimes == YES)
            PrintTimes (0, fpTimes);
 
		totalTreeLength = SumBranches(healthyRoot);
 
	/* allocate genotype data to be stored in data[genome][cell][site] */
        data = (int ***) calloc (ploidy, sizeof(int **));
        if (!data)
            {
            fprintf (stderr, "Could not allocate the data structure\n");
            exit (-1);
            }
        for (i=0; i<ploidy; i++)
            {
            data[i] = (int **) calloc (2*numCells+1, sizeof(int *));
            if (!data[i])
                {
                fprintf (stderr, "Could not allocate the data[] structure\n");
                exit (-1);
                }
            for (j=0; j<2*numCells+1; j++)
                {
                data[i][j] = (int *) calloc (numSites, sizeof(int));
                if (!data[i][j])
                    {
                    fprintf (stderr, "Could not allocate the data[][] structure\n");
                    exit (-1);
                    }
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
 
             /* assign ancestral states to all genomes */
            InitializeGenomes (healthyRoot, &seed);
			
			/* add germline variation if needed */
			if (SNPrate > 0)
				AddGermlineVariation (&seed);
				
            /* evolve maternal genome */
            if (noisy > 2)
                fprintf (stderr, "\n>> Evolving maternal genome ... ");
             EvolveSitesOnTree (healthyRoot, MATERNAL, &seed);
            if (noisy > 2)
                fprintf (stderr, "DONE");

            /* evolve paternal genome */
             if (noisy > 2)
                fprintf (stderr, "\n>> Evolving paternal genome ... ");
            EvolveSitesOnTree (healthyRoot, PATERNAL, &seed);
            if (noisy > 2)
                fprintf (stderr, "DONE\n\n");
	
            cumNumMU += numMU;
            cumNumMUSq += pow(numMU,2);

			if (CNLOHrate > 0)
				{
				/* evolve maternal CN_LOH */
				if (noisy > 2)
					fprintf (stderr, "\n>> Evolving maternal CN_LOH ... ");
				 EvolveCNLOHonTree(healthyRoot, MATERNAL, &seed);
				if (noisy > 2)
					fprintf (stderr, "DONE");

				/* evolve paternal CN_LOH  */
				 if (noisy > 2)
					fprintf (stderr, "\n>> Evolving paternal CN_LOH ... ");
				 EvolveCNLOHonTree(healthyRoot, PATERNAL, &seed);
				if (noisy > 2)
					fprintf (stderr, "DONE\n\n");

				cumNumCNLOH += numCNLOH;
				cumNumCNLOHSq += pow(numCNLOH,2);
				}

			if (deletionRate > 0)
				{
				/* evolve maternal deletions */
				if (noisy > 2)
					fprintf (stderr, "\n>> Evolving maternal deletions ... ");
				 EvolveDeletionsOnTree(healthyRoot, MATERNAL, &seed);
				if (noisy > 2)
					fprintf (stderr, "DONE");

				/* evolve paternal deletions  */
				 if (noisy > 2)
					fprintf (stderr, "\n>> Evolving paternal deletions ... ");
				 EvolveDeletionsOnTree(healthyRoot, PATERNAL, &seed);
				if (noisy > 2)
					fprintf (stderr, "DONE\n\n");

				cumNumDEL += numDEL;
				cumNumDELSq += pow(numDEL,2);
				}


			/* count how many SNVs we observe before ADO and genotype errors, but after CNLOH and deletion */
			numSNVs = CountTrueVariants();
				
			/* print reference haplotypes without errors */
			if (doPrintTrueHaplotypes == YES)
				{
  				if (doPrintSeparateReplicates == NO)
					fprintf (fpTrueHaplotypes, "[#%d]\n", dataSetNum+1);
                PrintFullHaplotypes(fpTrueHaplotypes);
				}

			/* copy true genotypes before doing ADO but after CNLOH and deletion */
			for (i=0; i<numCells+1; i++)
				for (j=0; j<numSites; j++)
					{
					cell[i].site[j].trueMaternalAllele = data[MATERNAL][i][j];
					cell[i].site[j].truePaternalAllele = data[PATERNAL][i][j];
					}

			/* do alellic dropout */
			if (ADOrate > 0)
				AllelicDropout(&seed);
			
			/* introduce errors directly in the genotypes */
			if (genotypingError > 0)
				GenotypeError(&seed);

			/* count how many alleles are at each site and how many SNVs we observe now, after ADO and genotype errors */
			numSNVs = CountAlleles();
			
            cumNumSNVs += numSNVs;
            cumNumSNVsSq += pow(numSNVs,2);

			if (numSNVs == 0)
                zeroSNVs++;

			if (noisy > 1)
				{
				fprintf (stderr, "\nData set %d",dataSetNum+1);
				fprintf (stderr, "\nNumber of coalescence events   =   %d", numCA);
				fprintf (stderr, "\nNumber of mutational events    =   %d", numMU);
				fprintf (stderr, "\nNumber of SNVs                 =   %d\n\n", numSNVs);
				fprintf (stderr, "\nNumber of deletion events      =   %d", numDEL);
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
            if (doPrintFullGenotypes == YES )
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
		if (doSimulateReadCounts == YES && doPrintSeparateReplicates == NO)
					fprintf (fpVCF, "[#%d]\n", dataSetNum+1);
			if (doPrintCATG == YES && doPrintSeparateReplicates == NO)
					fprintf (fpCATG, "[#%d]\n", dataSetNum+1);

			if (alphabet == DNA)
				{
				/* Generate NGS read counts */
				if (doSimulateReadCounts)
					{
					GenerateReadCounts(&seed);
					PrintVCF (fpVCF);
					if (doPrintCATG == YES)
						PrintCATG(fpCATG);
	
					if (doPrintMLhaplotypes == YES)
						{
						if (doPrintSeparateReplicates == NO)
							fprintf (fpMLhaplotypes, "[#%d]\n", dataSetNum+1);
						PrintMLhaplotypes(fpMLhaplotypes);
						}
					}
				}
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
                    if (data[0][i][j] == 1)
                        SiteMut[j]++;
                    if (data[1][i][j] == 1)
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
  
			
		// free all the necessary stuff for this replicate
		if (doUserTree == NO)
			free (treeNodes);
		
		if (doPrintSeparateReplicates == YES)
			{
			if (doPrintTree == YES)
				fclose(fpTrees);
			if (doPrintTimes == YES)
				fclose(fpTimes);
			}
		
        if (doSimulateData == YES)
            {
            for (i=0; i<ploidy; i++)
                {
                for (j=0; j<2*numCells+1; j++)
                    free (data[i][j]);
                free (data[i]);
                }
            free (data);
            free (allSites);
            free (SNVsites);
			free (variantSites);
            free (DefaultModelSites);
            free (AltModelSites);

			if (doGeneticSignatures == YES)
				{
				for (i=0; i<64; i++)
					{
					free(triNucleotideMaternal[i].position);
					free(triNucleotidePaternal[i].position);
					}
				free (triNucleotideMaternal);
				free (triNucleotidePaternal);
				}

			if (doPrintSeparateReplicates == YES)
				{
				if (doPrintSNVgenotypes == YES)
					fclose(fpSNVgenotypes);
				if (doPrintFullGenotypes == YES)
					fclose(fpFullGenotypes);
				if (doPrintSNVhaplotypes == YES)
					fclose(fpSNVhaplotypes);
				if (doPrintTrueHaplotypes == YES)
					fclose(fpTrueHaplotypes);
				if (doPrintMLhaplotypes == YES)
					fclose(fpMLhaplotypes);
				if (doPrintFullHaplotypes == YES)
					fclose(fpFullHaplotypes);
				if (doSimulateReadCounts == YES)
					fclose(fpVCF);
				if (doPrintCATG == YES)
					fclose(fpCATG);
				}
			}
		
		if (doSimulateReadCounts == YES)
			{
			for (i=0; i<numCells; i++)
				{
				for (j=0; j<numSites; j++)
					{
					for (k=0; k<4; k++)
						{
						free(cell[i].site[j].genLike[k]);
						free(cell[i].site[j].scaledGenLike[k]);
						free(cell[i].site[j].genLikeDoublet[k]);
						free(cell[i].site[j].scaledGenLikeDoublet[k]);
						}
					free(cell[i].site[j].readCount);
					free(cell[i].site[j].readCountDoublet);
					free(cell[i].site[j].genLike);
					free(cell[i].site[j].genLikeDoublet);
					free(cell[i].site[j].scaledGenLike);
					free(cell[i].site[j].scaledGenLikeDoublet);
					}
				}
			free (cell);
			}
		
        /* if(((dataSetNum*10)%numDataSets) == 0)
         PrintMemory(stderr); */
    } // replicate
		
    meanNumCA  = cumNumCA /  (double) numDataSets;
    meanNumMU  = cumNumMU /  (double) numDataSets;
    meanNumSNVs  = cumNumSNVs / (double) numDataSets;
    meanNumDEL  = cumNumDEL /  (double) numDataSets;
    meanNumCNLOH = cumNumCNLOH/  (double) numDataSets;
    meanTMRCA  =  cumTMRCA /  (double) numDataSets;
		
    varNumMU = (1.0 / (double) (numDataSets-1)) * (cumNumMUSq - pow(cumNumMU,2) / (double) numDataSets);
    varNumSNVs = (1.0 / (double) (numDataSets-1)) * (cumNumSNVsSq - pow(cumNumSNVs,2) / (double) numDataSets);
    varNumDEL = (1.0 / (double) (numDataSets-1)) * (cumNumDELSq - pow(cumNumDEL,2) / (double) numDataSets);
	varTMRCA = (1.0 / (double) (numDataSets-1)) * (cumTMRCASq - pow(cumTMRCA,2) / (double) numDataSets);

    if (noisy > 0)
        {
        /* if (noisy < 2)
         PrintHeader (stdout);*/
			
        fprintf (stdout, "\n\n** Simulations finished at ");
        PrintDate (stdout);
        PrintRunInformation (stdout);
			
		/* print run settings to file */
		PrintRunInformation (fpLog);
	
		fprintf (stdout, "\n\n\nOutput files are in folder \"%s\":", resultsDir);
        if (doPrintTree == YES)
            {
            fprintf (stdout, "\n Tree printed to file \"%s\"", treeFile);
			if (doPrintSeparateReplicates == NO)
				fclose(fpTrees);
			else
				fprintf (stdout, " in folder \"%s\"", treeDir);
	
			}
        if (doPrintTimes == YES)
            {
            fprintf (stdout, "\n Times printed to file \"%s\"", timesFile);
			if (doPrintSeparateReplicates == NO)
				fclose(fpTimes);
			else
				fprintf (stdout, " in folder \"%s\"", timesDir);
            }
			
        if (doSimulateData == YES)
            {
            if (doPrintSNVgenotypes == YES)
                {
                fprintf (stdout, "\n SNV genotypes printed to file \"%s\"", SNVgenotypesFile);
				if (doPrintSeparateReplicates == NO)
					fclose(fpSNVgenotypes);
			else
				fprintf (stdout, " in folder \"%s\"", SNVgenotypesDir);
               }
            if (doPrintFullGenotypes == YES)
                {
                fprintf (stdout, "\n Full genotypes printed to file \"%s\"", fullGenotypesFile);
 				if (doPrintSeparateReplicates == NO)
					fclose(fpFullGenotypes);
				else
					fprintf (stdout, " in folder \"%s\"", fullGenotypesDir);
                }
            if (doPrintSNVhaplotypes == YES)
                {
				if (doPrintIUPAChaplotypes == YES)
					fprintf (stdout, "\n SNV haplotypes (IUPAC codes) printed to file \"%s\"", SNVhaplotypesFile);
				else
					fprintf (stdout, "\n SNV haplotypes printed to file \"%s\"", SNVhaplotypesFile);
 				if (doPrintSeparateReplicates == NO)
					fclose(fpSNVhaplotypes);
				else
					fprintf (stdout, " in folder \"%s\"", SNVhaplotypesDir);
                }
            if (doPrintTrueHaplotypes == YES)
                {
				if (doPrintIUPAChaplotypes == YES)
					fprintf (stdout, "\n True haplotypes (IUPAC codes) printed to file \"%s\"", trueHaplotypesFile);
				else
					fprintf (stdout, "\n True haplotypes printed to file \"%s\"", trueHaplotypesFile);
 				if (doPrintSeparateReplicates == NO)
					fclose(fpTrueHaplotypes);
				else
					fprintf (stdout, " in folder \"%s\"", trueHaplotypesDir);
                }
			if (doPrintMLhaplotypes == YES)
				{
				if (doPrintIUPAChaplotypes == YES)
					fprintf (stdout, "\n ML haplotypes (IUPAC codes) printed to file \"%s\"", MLhaplotypesFile);
				else
					fprintf (stdout, "\n ML haplotypes printed to file \"%s\"", MLhaplotypesFile);
				if (doPrintSeparateReplicates == NO)
					fclose(fpMLhaplotypes);
				else
					fprintf (stdout, " in folder \"%s\"", MLhaplotypesDir);
				}
          if (doPrintFullHaplotypes == YES)
                {
  				if (doPrintIUPAChaplotypes == YES)
					fprintf (stdout, "\n Full haplotypes (IUPAC codes) printed to file \"%s\"", fullHaplotypesFile);
 				else
					fprintf (stdout, "\n Full haplotypes printed to file \"%s\"", fullHaplotypesFile);
 				if (doPrintSeparateReplicates == NO)
					fclose(fpFullHaplotypes);
				else
					fprintf (stdout, " in folder \"%s\"", fullHaplotypesDir);
                }
			if (doSimulateReadCounts == YES)
                {
                fprintf (stdout, "\n Genotype likelihoods printed to file \"%s\"", VCFfile);
 				if (doPrintSeparateReplicates == NO)
					fclose(fpVCF);
				else
					fprintf (stdout, " in folder \"%s\"", VCFdir);
                }
			if (doPrintCATG == YES)
                {
                fprintf (stdout, "\n Read counts printed to file \"%s\"", CATGfile);
 				if (doPrintSeparateReplicates == NO)
					fclose(fpCATG);
				else
					fprintf (stdout, " in folder \"%s\"", CATGdir);
                }
           }
        else
            fprintf (stdout, "\n\n\nNo output data files");
        }

		if (noisy > 2)
			for (i=0; i<numSites; i++)
				PrintSiteInfo (stderr, i);

	#ifdef CHECK_MUT_DISTRIBUTION
        fprintf (stdout,"\n\nChecking distribution of mutations");
        fprintf (stdout,"\n Avg # sites with  0 mutant alelles = %6.2f  (Exp = %6.2f)",(double)MutCount[0]/numDataSets, numSites-expNumMU);
        for (i=1; i<numCells; i++)
            fprintf (stdout,"\n Avg # sites with %2d mutant alelles = %6.2f  (Exp = %6.2f)",i,(double) MutCount[i]/numDataSets, theta/i);
        fprintf (stdout, "\n MeansumPos = %4.2f  (expected (L/2) = %4.2f)\n\n", meansumPos / (double) dataSetsWithSNVs, numSites/2.0);
    #endif
	
	//free pending stuff
	free (CommandLine);
	free (cellNames);
	free (triMutationsCounter);
	free (signatureWeight);
	free (signatureID);
	if (doUserGenome == YES)
		{
		free (paternalUserGenome);
		free (maternalUserGenome);
		}
	
	
    secs = (double)(clock() - start) / CLOCKS_PER_SEC;
		
    if (noisy > 0)
        {
        printf("\n\n_________________________________________________________________");
        printf("\nTime processing: %G seconds\n\n", secs);
        }
	return (1);
    }




/************************* MakeCoalescenceTree  ************************/
/*	Builds a genealogy  under the coalescent with exponential growth and no recombination */

void MakeCoalescenceTree(int numCells, int N, long int *seed)
	{
    int			i, *activeNodes, numActiveNodes,nextAvailableNode,
				firstNode, secondNode, newNode, eventNum,period;
    double		rateCA, timeCA, currentTime, eventTime;
    TreeNode	*p, *q, *r;
	
    p = q = r = NULL;
	
    /* allocate space for tree */
    treeNodes = (TreeNode *) calloc(numNodes, sizeof(TreeNode));
    if (!treeNodes)
        {
        fprintf (stderr, "Could not allocate nodes (%ld)\n", numNodes * sizeof(TreeNode));
        exit (-1);
        }
    activeNodes = (int *) calloc(numNodes, sizeof(int));
    if (!activeNodes)
        {
        fprintf (stderr, "Could not allocate activeNodes (%ld)\n", numNodes * sizeof(int));
        exit (-1);
        }
	
    /* set everything to null */
    for (i=0; i<numNodes; i++)
        {
        treeNodes[i].left = NULL;
        treeNodes[i].right = NULL;
        treeNodes[i].anc = NULL;
        treeNodes[i].time = 0.0;
        treeNodes[i].length = 0.0;
        treeNodes[i].index = 0;
        treeNodes[i].label = 0;
        treeNodes[i].isHealthyRoot = NO;
        treeNodes[i].isHealthyTip = NO;
		treeNodes[i].name = (char*) calloc (MAX_NAME, sizeof(char));
		if (!treeNodes[i].name)
			{
			fprintf (stderr, "Could not allocate the treeNodes[i].name structure\n");
			exit (-1);
			}
		}
	
    /* set up initial nodes */
    numActiveNodes = 0;
    for (i=0; i<numCells; i++)
        {
        p = treeNodes + i;
        activeNodes[numActiveNodes] = i;
        p->index = i;
        numActiveNodes++;
        }
    nextAvailableNode = numActiveNodes;
	
    /* make coalescence tree */
    eventNum = 0;
    period = 1;
    currentTime = 0.0;
    eventTime = 0.0;
    while (numActiveNodes > 1)
        {
        rateCA = numActiveNodes * (numActiveNodes - 1) / 2.0;
			
        if (doDemographics == YES)
            {
            periodGrowth[period] = -log (Nend[period] / (double) Nbegin[period]) / (double) (cumDuration[period] - cumDuration[period-1]);
            if (isnan(periodGrowth[period]) == YES)
                {
                fprintf (stderr, "\nERROR: period growth (%f) is NaN", periodGrowth[period]);
                fprintf (stderr, "\n       This might suggest that the growth rate is too negative");
                fprintf (stderr, "\n       and the coalescent time is therefore infinite.");
                fprintf (stderr, "\n       Try a smaller value");
                exit (-1);
                }
				
            timeCA =  log (1 + RandomExponential (rateCA, seed) * periodGrowth[period] * 2. * Nbegin[period] *
                           exp (-periodGrowth[period] * (currentTime - cumDuration[period-1]))) / periodGrowth[period]; /* time is kept in generations (scaled by 2N)  */
				
            /*	When growth rate is very negative, coalescent time may be infinite
             this results in log (-x) => timCA = NaN. If this not the last period
             just jump to the next. If this is the last period, we have to exit
             the program
             */
            if (isnan(timeCA) == YES)
                {
                if (period < numPeriods)
                    {
                    currentTime = cumDuration[period];
                    period++;
                    periodGrowth[period] = -log (Nend[period] / (double) Nbegin[period]) / (double) cumDuration[period];
                    continue;
                    }
                else
                    {
                    fprintf (stderr, "\nERROR: Coalescent time (%f) is infinite ", timeCA);
                    fprintf (stderr, "\n       This might suggest that the exponential growth rate is too negative");
                    fprintf (stderr, "\n       and the coalescent time is therefore infinite.");
                    fprintf (stderr, "\n       Try a smaller value");
                    exit (-1);
                    }
                }
            }
        else
            {
            timeCA = RandomExponential (rateCA, seed) * 2.0 * N; /* time is kept in 2N generations */
            if (doExponential == YES)
                {
                timeCA = log (exp(growthRate*currentTime) + growthRate * timeCA) / growthRate - currentTime;
					
                /*	When growth rate is very negative, coalescent time may be infinite
                 this results in log (-x) => timCA = NaN. We have to exit
                 the program
                 */
                if (isnan(timeCA) == YES)
                    {
                    fprintf (stderr, "\nERROR: Coalescent time (%f) is infinite ", timeCA);
                    fprintf (stderr, "\n       This might suggest that the exponential growth rate is too negative");
                    fprintf (stderr, "\n       and the coalescent time is therefore infinite.");
                    fprintf (stderr, "\n       Try a smaller value");
                    exit (-1);
                    }
                }
            }
			
        eventTime = timeCA;
			
        /*	if this period is not the last one and if the event time is outside the current interval,
            update period and start again */
        if (doDemographics == YES && period < numPeriods  &&  (currentTime + eventTime) > cumDuration [period])
            {
            currentTime = cumDuration[period];
            period++;
            periodGrowth[period] = -log (Nend[period] - Nbegin[period]) / (double) cumDuration[period];
            continue;
            }
        numCA++;
			
        /* set time */
        currentTime += eventTime;
			
        if (noisy > 2)
            fprintf (stderr, "\nEvent %3d :  rateCA = %lf, currentTime = %lf\n", ++eventNum, rateCA, currentTime);
			
        /* figure out which two nodes will coalesce */
        firstNode = numActiveNodes * RandomUniform(seed);
        if (firstNode >= numActiveNodes)
            {
            fprintf (stderr, "\n\nERROR: firstNode out of range!\n");
            exit (-1);
            }
        do
            {
            secondNode = numActiveNodes * RandomUniform(seed);
            } while (firstNode == secondNode);
			
        newNode = nextAvailableNode;
        if (noisy > 2)
            fprintf (stderr, "  Coalescence involving %d and %d to create node %d\n",
                     activeNodes[firstNode]+1, activeNodes[secondNode]+1, newNode+1);
			
        /* readjust MRCA, pointers and active sites */
        p = treeNodes + activeNodes[firstNode];
        q = treeNodes + activeNodes[secondNode];
        r = treeNodes + newNode;
        r->index = nextAvailableNode;
        r->left = p;
        r->right = q;
        p->anc = r;
        q->anc = r;
        r->time = currentTime;
        p->length = r->time - p->time;
        q->length = r->time - q->time;
        p->branchLength = p->length * mutationRate;
        q->branchLength = q->length * mutationRate;
			
        /* readjust active nodes */
        activeNodes[firstNode] = newNode;
        activeNodes[secondNode] = activeNodes[numActiveNodes-1];
        numActiveNodes--;
        nextAvailableNode++;
        if (nextAvailableNode >= numNodes)
            {
            fprintf (stderr, "\n\nERROR: Too many nodes!\n");
            exit (-1);
            }
    } /* coalescent cell tree finished */
	
    TMRCA = currentTime / (2.0 * N);  /* so TMRCA is in coalescent generations */
    coalTreeMRCA = r;
	
    /* connect the coalescent cell MRCA node with the healthy ancestral cell*/
    if (noisy > 2)
        fprintf (stderr, "\n>> Adding healthy root ... ");
    healthyRoot = treeNodes + nextAvailableNode;
    healthyRoot->index = nextAvailableNode;
    healthyRoot->left = coalTreeMRCA;
    coalTreeMRCA->anc = healthyRoot;
    coalTreeMRCA->length = transformingBranchLength/mutationRate;
    coalTreeMRCA->branchLength = transformingBranchLength;
    healthyRoot->time = currentTime +  transformingBranchLength/mutationRate;
    healthyRoot->length = 0;
    healthyRoot->isHealthyRoot = YES;
    if (noisy > 2)
        fprintf (stderr, "DONE");

    nextAvailableNode++;
	
    /* connect the healthy ancestral cell with the tip healthy cell*/
    if (noisy > 2)
        fprintf (stderr, "\n>> Adding healthy tip ... ");
    healthyTip = treeNodes + nextAvailableNode;
    healthyTip->left = NULL;
    healthyTip->right = NULL;
    healthyTip->anc = healthyRoot;
    healthyRoot->right = healthyTip;
    healthyTip->time = 0;
    healthyTip->length = healthyTipBranchLength/mutationRate;
    healthyTip->branchLength = healthyTipBranchLength;
    healthyTip->isHealthyTip = YES;
    healthyTip->index = nextAvailableNode;
    if (noisy > 2)
        fprintf (stderr, "DONE");
	
    /* Relabel nodes on tree (label tips with indexes) */
    if (noisy > 2)
        fprintf (stderr, "\n>> Relabeling nodes on tree ... ");
    tipLabel = numCells;
    intLabel = numCells+1;
    RelabelNodes(healthyRoot);
    if (noisy > 2)
        fprintf (stderr, "DONE");
	
    free (activeNodes);
	}

/********************************** ReadUserGenome ***********************************/
/* Reads a FASTA with a diploid reference genome */

void ReadUserGenome (FILE *fp)
	{
	char	*header;
	

	/* allocate space for reading */
	header = (char*) calloc (MAX_LINE, sizeof(char));
	if (!header)
		{
		fprintf (stderr, "Could not allocate the header string\n");
		exit (-1);
		}

	maternalUserGenome = (char*) calloc (MAX_GENOME, sizeof(char));
	if (!maternalUserGenome)
		{
		fprintf (stderr, "Could not allocate the maternalUserGenome structure\n");
		exit (-1);
		}

	paternalUserGenome = (char*) calloc (MAX_GENOME, sizeof(char));
	if (!paternalUserGenome)
		{
		fprintf (stderr, "Could not allocate the paternalUserGenome structure\n");
		exit (-1);
		}
	
	/* open user genome */
	if ((fpUserGenome = fopen(userGenomeFile, "r")) == NULL)
		{
		fprintf (stderr, "\nERROR: Can't open user genome file \"%s\"\n", userGenomeFile);
		PrintUsage();
		}

	/* read the FASTA file */
	fgets (header, MAX_LINE, fpUserGenome);
	fgets (maternalUserGenome, MAX_GENOME, fpUserGenome);
	if (maternalUserGenome == NULL)
		{
		fprintf (stderr, "\nERROR: Can't read the user maternal genome");
		PrintUsage();
		}
	
	fgets (header, MAX_LINE, fpUserGenome);
	fgets (paternalUserGenome, MAX_GENOME, fpUserGenome);
	if (paternalUserGenome == NULL)
		{
		fprintf (stderr, "\nERROR: Can't read the user paternal genome");
		PrintUsage();
		}

	/* check lengths */
	if (strlen(maternalUserGenome) != strlen(paternalUserGenome))
		{
		fprintf (stderr, "\nERROR: maternal and paternal user genomes have different lengths");
		PrintUsage();
		}

	numSites = (int) strlen(maternalUserGenome) - 1;
	
	/* realloc memory */
	maternalUserGenome = (char *)realloc(maternalUserGenome, numSites * sizeof(char));
	paternalUserGenome = (char *)realloc(paternalUserGenome, numSites * sizeof(char));


	fprintf (stderr, "\nnumsites = %d", numSites);
	fprintf (stderr, "\nuser maternal genome = %s", maternalUserGenome);
	fprintf (stderr, "\nuser paternal genome = %s", paternalUserGenome);

	free (header);
	fclose(fpUserGenome);
	}


static void ReadUserTree (FILE *fp)
	{
	int			taxonNumber, numDigits;
	int			nextAvailableNode;
	int			i, j;
	int			digit, digit1, digit2, digit3;
	char		*tempString;
	TreeNode	*p, *q;
	
	p = q = NULL;
	
	/* allocate space for reading the tree string */
	treeString = (char*) calloc (MAX_LINE, sizeof(char));
	if (!treeString)
		{
		fprintf (stderr, "Could not allocate the treeString structure\n");
		exit (-1);
		}

	tempString = (char*) calloc (MAX_NAME, sizeof(char));
        if (!tempString)
            {
            fprintf (stderr, "Could not allocate the tempString structure\n");
            exit (-1);
            }
	
	/* allocate space for tree */
    treeNodes = (TreeNode *) calloc(numNodes, sizeof(TreeNode));
    if (!treeNodes)
        {
        fprintf (stderr, "Could not allocate nodes (%ld)\n", numNodes * sizeof(TreeNode));
        exit (-1);
        }
	
	/* allocate space for cellNames */
    cellNames = (char **) calloc(numNodes, sizeof(char *));
    if (!cellNames)
        {
        fprintf (stderr, "Could not allocate cellNames (%ld)\n", numNodes * sizeof(char *));
        exit (-1);
        }
   for (i=0; i<numNodes; i++)
	   {
	   cellNames[i] = (char *) calloc(MAX_NAME, sizeof(char));
	   if (!cellNames[i])
		   {
		   fprintf (stderr, "Could not allocate cellNames[%d] (%ld)\n", i,  MAX_NAME * sizeof(char));
		   exit (-1);
		   }
	   }
	
	/* initialize tree nodes */
    for (i=0; i<numNodes; i++)
        {
        treeNodes[i].left = NULL;
        treeNodes[i].right = NULL;
        treeNodes[i].anc = NULL;
        treeNodes[i].time = 0.0;
        treeNodes[i].length = 0.0;
        treeNodes[i].index = -1;
        treeNodes[i].label = -1;
        treeNodes[i].isHealthyRoot = NO;
        treeNodes[i].isHealthyTip = NO;
		treeNodes[i].name = (char*) calloc (MAX_NAME, sizeof(char));
		if (!treeNodes[i].name)
			{
			fprintf (stderr, "Could not allocate the treeNodes[i].name structure\n");
			PrintUsage();
			}
		}

	/* open treefile */
	if ((fpUserTree = fopen(userTreeFile, "r")) == NULL)
		{
		fprintf (stderr, "\nERROR: Can't open user treefile \"%s\"\n", userTreeFile);
		PrintUsage();
		}

	/* read the tree string */
	fgets (treeString, MAX_LINE, fpUserTree);
	if (treeString == NULL)
		{
		fprintf (stderr, "\nERROR: Can't read tree string");
		PrintUsage();
		}
	fclose(fpUserTree);
	
	CheckTree(treeString);

	/* build the tree linked list */
    numCells = 0;
    nextAvailableNode = 0;
    i = 0;
    do
        {
        if (treeString[i] == '(')
            {
            if (nextAvailableNode == 0)
                {
				p = treeNodes + nextAvailableNode;
                coalTreeMRCA = p;
				p->index = nextAvailableNode;
				nextAvailableNode++;
				}
            else
                {
				q = treeNodes + nextAvailableNode;
				if (p->left == NULL)
					{
					p->left = q;
					q->anc = p;
					}
				else
					{
					p->right = q;
					q->anc = p;
					}
				p = q;
				p->index = nextAvailableNode;
                nextAvailableNode++;
				}
            }
        else if (treeString[i] == ')')
            {
            q = p->anc;
            p = q;
            }
        else if (treeString[i] == ',')
            {
            q = p->anc;
            p = q;
            }
		else if (treeString[i] == ':')
	        {
	        i++;
            for (j=0; j<10; j++)
                tempString[j] = ' ';
            j = 0;
	        while (treeString[i] != ',' && treeString[i] != ')' && treeString[i] != ':' && treeString[i] != ';' && treeString[i] != '(')
	        	tempString[j++] = treeString[i++];
	        i--;
	        sscanf (tempString, "%lf", &p->branchLength);
			}
        else
        	/* now read taxon names*/
        	{
			if (/* DISABLES CODE */ (0)) /* let's read every name as a word */
            /* if (isdigit(treeString[i]) == YES)*/ /* taxon names are numbers */
            	{
	            numDigits = 0;
	            digit1 = digit2 = digit3 = 0;
	            while (treeString[i] != ',' && treeString[i] != ')' && treeString[i] != ':')
	                    {
	                    if      (treeString[i] == '0') digit = 0;
	                    else if (treeString[i] == '1') digit = 1;
	                    else if (treeString[i] == '2') digit = 2;
	                    else if (treeString[i] == '3') digit = 3;
	                    else if (treeString[i] == '4') digit = 4;
	                    else if (treeString[i] == '5') digit = 5;
	                    else if (treeString[i] == '6') digit = 6;
	                    else if (treeString[i] == '7') digit = 7;
	                    else if (treeString[i] == '8') digit = 8;
	                    else if (treeString[i] == '9') digit = 9;
	                    if      (numDigits == 0) digit1 = digit;
	                    else if (numDigits == 1) digit2 = digit;
	                    else if (numDigits == 2) digit3 = digit;
	                    numDigits++;
	                    i++;
	                    }
	            i--;
	            if (numDigits == 1)
	            	taxonNumber = digit1;
	            else if (numDigits == 2)
					taxonNumber = 10*digit1 + digit2;
	            else if (numDigits == 3)
	            	taxonNumber = 100*digit1 + 10*digit2 + digit3;
				}
			else /* taxon names are letters */
	            {
				taxonName = (char*) calloc (MAX_NAME, sizeof(char));
				if (!taxonName)
					{
					fprintf (stderr, "Could not allocate the taxonName structure\n");
					exit (-1);
					}
	           	taxonNamesAreChars = YES;
	            j = 0;
	            while (treeString[i] != ',' && treeString[i] != ')' && treeString[i] != ':')
	            	taxonName[j++] = treeString[i++];
	            i--;
				}
            q = treeNodes + nextAvailableNode;
			if (p->left == NULL)
				{
				p->left = q;
				q->anc = p;
				}
			else
				{
				p->right = q;
				q->anc = p;
				}
			strcpy(q->name,taxonName);
			strcpy(cellNames[numCells],taxonName);
            p = q;
			p->index = nextAvailableNode;
           	nextAvailableNode++;
            numCells++;
 			free (taxonName);
           }
        i++;
        } while (treeString[i] != ';');
	
	/* connect the coalescent cell MRCA node with the outgroup ancestral cell*/
    if (noisy > 2)
        fprintf (stderr, "\n>> Adding healthy root ... ");
    healthyRoot = treeNodes + nextAvailableNode;
    healthyRoot->left = coalTreeMRCA;
    coalTreeMRCA->anc = healthyRoot;
    //coalTreeMRCA->length = transformingBranchLength/mutationRate;
    coalTreeMRCA->branchLength = transformingBranchLength;
    //healthyRoot->time = currentTime +  transformingBranchLength/mutationRate;
    healthyRoot->length = 0;
    healthyRoot->name = "healthyroot";
    healthyRoot->isHealthyRoot = YES;
    healthyRoot->index = nextAvailableNode;
   if (noisy > 2)
        fprintf (stderr, "DONE");

    nextAvailableNode++;
	
    /* connect the healthy ancestral cell with the tip healthy cell*/
    if (noisy > 2)
        fprintf (stderr, "\n>> Adding healthy tip ... ");
    healthyTip = treeNodes + nextAvailableNode;
    healthyTip->left = NULL;
    healthyTip->right = NULL;
    healthyTip->anc = healthyRoot;
    healthyRoot->right = healthyTip;
    healthyTip->time = 0;
    //healthyTip->length = healthyTipBranchLength/mutationRate;
    healthyTip->branchLength = healthyTipBranchLength;
    healthyTip->isHealthyTip = YES;
	healthyTip->name = "healthy";
	healthyTip->index = nextAvailableNode;
   if (noisy > 2)
        fprintf (stderr, "DONE");

    /* Relabel nodes on tree (label tips with consecutive indexes) */
    if (noisy > 2)
        fprintf (stderr, "\n>> Relabeling nodes on tree ... ");

	TipNodeNum = 0;
	IntNodeNum = numCells+1;
	RelabelUserTree (healthyRoot);

    if (noisy > 2)
        fprintf (stderr, "DONE");

	numNodes = nextAvailableNode+1;

	/* reallocate tree nodes */
	treeNodes = (TreeNode *) realloc (treeNodes, numNodes * sizeof (TreeNode));
	if (treeNodes == NULL)
		{
		fprintf (stderr, "Could not reallocate treeNodes (%ld bytes)\n", numNodes * sizeof (TreeNode));
		PrintUsage();
		}

	fprintf (fpLog, "\n\nUser-defined tree: ");
	PrintTree (healthyRoot, fpLog);

	free (tempString);
	}


/***************** CheckTree *******************/
/*  Checks for an error in the input tree string
	Prints message when there is an error in the
	tree and prints correct part of the tree */
static void		CheckTree (char *treeStr)
	{
	int		i, k;
	int		thereIsError, numLeftPar, numRightPar;
	
	thereIsError = NO;
    numLeftPar = 1;
    numRightPar = 0;
	i = 1;

	do
		{
		if (treeStr[i] == '(')
			{
			if (treeStr[i-1] != '(' && treeStr[i-1] != ',')
				thereIsError = YES;
			if (treeStr[i+1] != '(' && isalnum(treeStr[i+1]) == NO)
				thereIsError = YES;
			numLeftPar++;
			}
		else if (treeStr[i] == ')')
			{
			if (treeStr[i-1] != ')' && isalnum(treeStr[i-1]) == NO)
				thereIsError = YES;
			if (treeStr[i+1] != ')' && treeStr[i+1] != ',' && treeStr[i+1] != ':')
				thereIsError = YES;
			if (treeStr[i+1] == ';' && i == strlen(treeStr)-3)
				thereIsError = NO;
			numRightPar++;
			}
		else if (treeStr[i] == ',')
			{
			if (treeStr[i-1] != ')' && isalnum(treeStr[i-1]) == NO)
				thereIsError = YES;
			if (treeStr[i+1] != '(' && isalnum(treeStr[i+1]) == NO)
				thereIsError = YES;
			}
		else if (treeStr[i] == ':')
			{
			if (treeStr[i-1] != ')' && isalnum(treeStr[i-1]) == NO)
				thereIsError = YES;
			if (isdigit(treeStr[i+1]) == NO)
				thereIsError = YES;
			}

		if (thereIsError == YES)
			{
			fprintf (stderr,"\nERROR: There is something wrong in the tree\n");
			for (k=0; k<=i; k++)
				fprintf (stderr,"%c", treeStr[k]);
			fprintf (stderr," <- HERE");
			PrintUsage();
			}
		
		i++;
	
		} while (treeStr[i] != ';');

	if (numLeftPar != numRightPar)
		{
		fprintf (stderr, "Tree seems unbalanced (%d left and %d right parentheses)", numLeftPar, numRightPar);
		PrintUsage();
		}
	if (strrchr(treeString, ':') == NULL)
		{
		fprintf (stderr, "Tree does not have branch lengths");
		PrintUsage();
		}
	}

/****************** RelabelUserTree **************************/
/*	Reindex nodes, count nodes to tips,
	and initialize components in the user tree */

void RelabelUserTree (TreeNode *p)
	{
    if (p != NULL)
        {
		RelabelUserTree (p->left);
		RelabelUserTree (p->right);
		if (p->left == NULL && p->right == NULL)
			p->index = TipNodeNum++;
		else
			p->index = IntNodeNum++;
		p->label = p->index;
		}
	}


/**************** RelabelNodes **************/
/*	Relabel nodes so all tips are consecutive */

void RelabelNodes (TreeNode *p)
    {
    if (p != NULL)
        {
        RelabelNodes (p->left);
        RelabelNodes (p->right);
        if (p->left == NULL && p->right == NULL && p->isHealthyTip == NO) /* is cell tip */
            p->label = p->index;
        else if (p->left == NULL && p->right == NULL && p->isHealthyTip == YES) /* is outgroup tip */
                p->label = tipLabel++;
        else
            p->label = intLabel++; /* is internal */
        }
    }


/**************** MakeTreeNonClock **************/
/*	Introduce rate variation among lineages using gamma rate  */

void MakeTreeNonClock (TreeNode *p, long int *seed)
    {
    if (p != NULL)
        {
        p->branchLength *=  RandomGamma (alphaBranches, seed) / alphaBranches;
        MakeTreeNonClock (p->left, seed);
        MakeTreeNonClock (p->right, seed);
        }
    }


/********************************** InitializeGenomes ***********************************/
/* Initialize all genomes with the reference states  */
void InitializeGenomes (TreeNode *p, long int *seed)
    {
    int     i, cell, anccell, site;
    double  ran;
		
    if (p != NULL)
        {
        cell = p->label;
        if (alphabet == DNA)
            {
            if (p->isHealthyRoot == YES)
                {
				if (doUserGenome == YES)
					{
					for (site=0; site<numSites; site++) /* initialize genome with the user provided sequences */
						{
						data[MATERNAL][cell][site] = WhichNucChar(maternalUserGenome[site]);
						data[PATERNAL][cell][site] = WhichNucChar(paternalUserGenome[site]);
						allSites[site].referenceAllele = data[MATERNAL][cell][site];   // then allSites[site].referenceAllele hosts the reference genome
						}
					}
  				else if (doGeneticSignatures == NO) /* initialize genome with mononucleotide frequencies */
					{
					for (site=0; site<numSites; site++)
						{
						ran = RandomUniform(seed);
						for (i=0; i<4; i++)
							{
							if (ran <= cumfreq[i])
								{
								data[MATERNAL][cell][site] = data[PATERNAL][cell][site] = i;
								allSites[site].referenceAllele = data[MATERNAL][cell][site];   // then allSites[site].referenceAllele hosts the reference genome
								break;
								}
							}
						}
					}
				else /* initialize genome with trinucleotide frequencies */
					{
					SimulateTriNucFreqGenome (cell, seed);
					}
				}
            else
                {
                anccell = p->anc->label;
                for (site=0; site<numSites; site++)
                    {
                    data[MATERNAL][cell][site] = data[MATERNAL][anccell][site];
                    data[PATERNAL][cell][site] = data[PATERNAL][anccell][site];
                    }
                }
             }
        else
            for (site=0; site<numSites; site++)
                data[MATERNAL][cell][site] = data[PATERNAL][cell][site] = 0;

        InitializeGenomes (p->left, seed);
        InitializeGenomes (p->right, seed);
        }
    }



/********************************** SimulateTriNucFreqGenome ***********************************/
/*
	 Simulate a homozygous diploid genome with the trinucleotide frequencies of the human genome

	 We simulate first a random trinucleotide.
	 Then we go site by site taking into account the last two letters  of the
	 previous trinucleotide.  if the first trinucleotide is, for example, CAT
	 then we simulate trinucleotide starting in AT_  (i.e. ATA, ATC, ATG or ATT),
	 and keep going  taking into account the last two letters of the new trinucleotide,
	 and so on.
*/

void SimulateTriNucFreqGenome (int cell, long int *seed)
	{
	int 		chosenTriNucleotide, nextNucleotide;
	int			k, n1, n2, n3, rest, site;
	double 		*prob4, sum;
	
	/* memory allocations */
	prob4 = (double *) calloc (4, sizeof(double));
	if (!prob4)
		{
		fprintf (stderr, "Could not allocate the prob4 vector\n");
		exit (-1);
		}
		
	/* choose first trinucleotide */
	chosenTriNucleotide = ChooseUniformState(triNucFreq, seed);

	/* find bases of the selected trinucleotide */
	n1 = chosenTriNucleotide/16;
	rest = chosenTriNucleotide%16;
	n2 = rest/4;
	n3 = rest%4;
	//fprintf (stderr, "\n%2d %c%c%c ", chosenTriNucleotide, WhichNuc(n1), WhichNuc(n2), WhichNuc(n3));

	site = 0;
	data[MATERNAL][cell][site] = data[PATERNAL][cell][site] = allSites[site].referenceAllele = n1;
	site++;
	data[MATERNAL][cell][site] = data[PATERNAL][cell][site] = allSites[site].referenceAllele = n2;
	site++;
	data[MATERNAL][cell][site] = data[PATERNAL][cell][site] = allSites[site].referenceAllele = n3;
	
	/* fill the rest of the genome */
	/* choose next nucleotide given the last two bases of the previous trinucleotide  */
	for (site=3; site<numSites; site++)
		{
		/* normalize frequencies given the last two bases */
		sum = 0;
		for (k=0; k<4; k++)
			sum += triNucFreq[trinuc(n2,n3,k)];
		for (k=0; k<4; k++)
			prob4[k] = triNucFreq[trinuc(n2,n3,k)] / sum;
		
		nextNucleotide = ChooseUniformState(prob4, seed);
		data[MATERNAL][cell][site] = data[PATERNAL][cell][site] = allSites[site].referenceAllele = nextNucleotide;

		/* move downstream one position */
		n1 = n2;
		n2 = n3;
		n3 = nextNucleotide;
		}
 
	free(prob4);
	}



/********************************** RecordTriNucObservation ***********************************/
/* add a new occurrence of a trinucleotide along the genome */
void RecordTriNucObservation (TriNucStr *trin, int site)
	{
	if (trin->numAvailablePositions == trin->tempLength)
	  {
	  trin->tempLength *= 2;
	  trin->position = (int *)realloc(trin->position, trin->tempLength * sizeof(int));
	  if (!trin->position)
		  {
		  fprintf (stderr, "Could not reallocate trinucleotide position array when adding\n");
		  exit (-1);
		  }
	  }
	trin->position[trin->numAvailablePositions++] = site;
  	}




/********************************** AddGermlineVariation ***********************************/
/* Introduce germline SNPs in the healthy root genome for DNA or 0/1 data  */
void AddGermlineVariation (long int *seed)
    {
	int		j, site, genome, current_state;
	double	ran;
	
	for (site=0; site<numSites; site++)
		{
		if(RandomUniform(seed) < SNPrate) //there is a SNP at this site
			{
			allSites[site].isSNP = YES;
			if(RandomUniform(seed) < 0.5)  // in maternal/paternal genome
				genome = MATERNAL;
			else
				genome = PATERNAL;
		
			if (alphabet == BINARY)
				data[genome][HEALTHY_ROOT][site] = 1;
			else // DNA mutation according to Mij matrix
				{
				current_state = data[genome][HEALTHY_ROOT][site];
				ran = RandomUniform(seed) * cumMij[current_state][3];
                for (j=0; j<4; j++)
                    {
					if (ran <= cumMij[current_state][j])
                        {
                        data[genome][HEALTHY_ROOT][site] = j;
                        break;
                        }
					}
				}
			}
		}
	}


/********************************** EvolveSitesOnTree ***********************************/
/* Evolves all sites (maternal and paternal genomes) on the given tree
   We assume that a site will be ISM or Mk in both maternal and paternal genome
 */
void EvolveSitesOnTree (TreeNode *treeRoot, int genome, long int *seed)
    {
    int i = 0;
		
    if (rateVarAmongSites == YES)
        for (i=0; i<numSites; i++)
            allSites[i].rateMultiplier = RandomGamma (alphaSites, seed) / alphaSites;

	if (doGeneticSignatures == YES)
		{
		SimulateSignatureISM(treeRoot, genome, seed); /* we assume diploid ISM for genetic signatures */
		}
	else
		{
		if (propAltModelSites == 0)  /* only default model (ISM diploid) sites */
			{
			numDefaultModelSites = numSites;
			numAltModelSites = 0;
			for (i=0; i<numSites; i++)
				DefaultModelSites[i] = i;
			SimulateISM (treeRoot, genome, NO, seed);
			}
	   else if (propAltModelSites == 1.0)  /* only alternative model sites */
			{
			numDefaultModelSites = 0;
			numAltModelSites = numSites;
			for (i=0; i<numSites; i++)
				AltModelSites[i] = i;
			if (altModel == ISMhap)
				SimulateISM (treeRoot, genome, YES, seed);
			else if (altModel == Mk2)
				SimulateMk2 (treeRoot, genome, seed);
			else if (altModel == finiteDNA)
				SimulateFiniteDNA (treeRoot, genome, seed);
			else
				{
				fprintf (stderr, "\n\nERROR: Sorry, the specified model is unknown for me");
				exit(0);
				}
			  }
		else /* both ISM and non-ISM sites */
			{
			numDefaultModelSites = 0;
			numAltModelSites = 0;

			for (i=0; i<numSites; i++)
				{
				if (RandomUniform (seed) < propAltModelSites)
					AltModelSites[numAltModelSites++] = i;
				else
					DefaultModelSites[numDefaultModelSites++] = i;
				}
	 
			/* Evolve ISM sites */
			if (numDefaultModelSites > 0)
				SimulateISM (treeRoot, genome, NO, seed);
	 
			/* Evolve non-SIM sites */
			if (numAltModelSites > 0)
				{
				if (altModel == ISMhap)
					SimulateISM (treeRoot, genome, YES, seed);
				else if (altModel == Mk2)
					SimulateMk2 (treeRoot, genome, seed);
				else if (altModel == finiteDNA)
					SimulateFiniteDNA (treeRoot, genome, seed);
				else
					{
					fprintf (stderr, "\n\nERROR: Sorry, the specified model is unknown for me");
					exit(-1);
					}
				}
			}
		}
    }


/*************************** SimulateISM **********************************/
/*	Simulates mutations under an infinite diploid/haploid sites model (ISM).

    Diploid here means that if site 33 mutates in the female genome,
    site 33 in the male genome will not mutate, and viceversa.

    Haploid here means that site 33 in the maternal genome and site 33 in the paternal
    genome will be considered different sites, so both can mutate. The rest is as
    in the diploid, standard ISM
 
    The number of mutations will be distributed as:
 
    (1) a Poisson with parameter the sum of
	branch lengths over all sites  Times are already scaled in
	2N. We will force it to get just one mutation per site
 
	(2) A used-defined number of mutations (segregating sites in the ISM diploid)
 
	Approach 1, which conditions on theta and takes in account
	the branch lengths, is preferred than Hudson's approach 2, which selects per each
	site a random branch and puts a mutation there, and therefore it is
	conditioned on the number of segregating sites. With the coalescent
	you get large or short trees, which imply different number of mutations
	but if you do Hudson's approach you always put the same number of mutations.
*/

void SimulateISM (TreeNode *p, int genome, int doISMhaploid, long int *seed)
    {
    int		i, trials, numMutations, mutationsSoFar;
    double	totalBranchSum;
    int		*modelSites, numModelSites;

	if (doISMhaploid == NO)
		{
		modelSites = DefaultModelSites;
		numModelSites = numDefaultModelSites;
		}
	else
 		{
		modelSites = AltModelSites;
		numModelSites = numAltModelSites;
		}

/*
    totalBranchSum = 0;
    for (i=0; i<numModelSites; i++)
        {
        allSites[modelSites[i]].branchSum = SumBranches (healthyRoot);
        totalBranchSum += allSites[modelSites[i]].branchSum;
        }

*/

	totalBranchSum = totalTreeLength * numModelSites;

    if (doSimulateFixedNumMutations == YES) /* conditioned on a specific number of mutations */
        {
        if (genome == MATERNAL)
			 numMutations = numSNVmaternal =  RandomBinomial (0.5, numFixedMutations, seed);
        else
            numMutations = numFixedMutations - numSNVmaternal;
        }
    else
        numMutations = RandomPoisson (totalBranchSum, seed); /* the number of mutations will be distributed as a Poisson with parameter totalBranchSum */
		
     /* if the number of ISM mutations is bigger than the number of ISM sites quit the program and warn the user about violation of ISM */
     if (doISMhaploid == NO && numISMmutations + numMutations > (numModelSites))
        {
        fprintf (stderr, "\n\nERROR: The diploid infinite sites model (ISM) has been violated. There will be");
        fprintf (stderr, "\nmore mutations (%d existing + %d proposed]) than available sites under this model (%d).", numISMmutations, numMutations, numDefaultModelSites);
        fprintf (stderr, "\nTry using a smaller mutation rate or effective population size");
        fprintf (stderr, "\n(i.e. smaller theta) or check the proportion of JC sites.");
        fprintf (stderr, "\nIf using an user tree, try with smaller branch lengths (including the transforming and healthy branches).\n");
        exit (-1);
        }
    else if (doISMhaploid == YES && numISMmutations + numMutations > (2*numModelSites))
        {
        fprintf (stderr, "\n\nERROR: The haploid infinite sites model (ISM) has been violated. There will be");
        fprintf (stderr, "\nmore mutations (%d existing + %d proposed]) than available sites under this model (%d).", numISMmutations, numMutations, 2*numModelSites);
        fprintf (stderr, "\nTry using a smaller mutation rate or effective population size");
        fprintf (stderr, "\n(i.e. smaller theta) or check the proportion of JC sites.");
        fprintf (stderr, "\nIf using an user tree, try with smaller branch lengths (including the transforming and healthy branches).\n");
        exit (-1);
        }
		
    numISMmutations += numMutations;
    trials = 0;
    mutationsSoFar = 0;
    while (mutationsSoFar < numMutations)
        {
        i = RandomUniformTo(numModelSites, seed);  /* choose a site at random */
        while (((doISMhaploid == NO)  && (allSites[modelSites[i]].numMutations != 0))  ||
               ((doISMhaploid == YES) && (genome == MATERNAL) && (allSites[modelSites[i]].numMutationsMaternal != 0)) ||
               ((doISMhaploid == YES) && (genome == PATERNAL) && (allSites[modelSites[i]].numMutationsPaternal != 0)))
            {
            i = RandomUniformTo(numModelSites, seed);
            if (trials++ > 100*numModelSites)
                {
                fprintf (stderr, "\n\n ERROR: after %d trials cannot find an unmuted site",100*numModelSites);
                fprintf (stderr, "\nmutations = %d   mutations so far = %d\n\n",numMutations, mutationsSoFar);
                for (i=0; i<numModelSites; i++)
                    fprintf (stderr, "\nsite %d has %d mutations",modelSites[i]+1, allSites[modelSites[i]].numMutations);
                exit(-1);
                }
            }
        if (alphabet == DNA)
            SimulateISMDNAforSite (p, genome, modelSites[i], doISMhaploid, seed);
        else
            SimulateISMforSite (p, genome, modelSites[i], doISMhaploid, seed);
			
        mutationsSoFar++;
			
        #ifdef MYDEBUG
            fprintf (stderr, "\nmutations = %d   mutations so far = %d\n",numMutations, mutationsSoFar);
            if (allSites[modelSites[i]].numMutations > 1)
                {
                fprintf (stderr, "\n\n ERROR: %d mutations in site %d",allSites[modelSites[i]].numMutations, modelSites[i]+1);
                exit(-1);
                }
            for (i=0; i<numDefaultModelSites; i++)
                fprintf (stderr, "%2d[%d] ",modelSites[i]+1, allSites[modelSites[i]].numMutations);
        #endif
			
        }
    }


/********************************** SimulateISMForSite ***********************************/
/*	Simulates a 0/1 mutation under an infinite sites model (ISM) for a given site. The branch
	where this mutation is placed is chosen according to its length.
    0 is the reference (healthy) allele
 */
void SimulateISMforSite (TreeNode *p, int genome, int site, int doISMhaploid, long int *seed)
    {
    static double	cumBranchLength, uniform;
    int             cell, anccell;
		
	if (p != NULL)
        {
        cell = p->label;
			
        if (p->isHealthyRoot == YES)
            {
            cumBranchLength = 0;
            uniform = RandomUniform(seed) * totalTreeLength;
            }
        else
            {
            anccell = p->anc->label;
            cumBranchLength += p->branchLength;
				
             if ((cumBranchLength < uniform) || /* => there will be no change */
                ((doISMhaploid == NO)  && (allSites[site].numMutations > 0))  ||
                ((doISMhaploid == YES) && (genome == MATERNAL) && (allSites[site].numMutationsMaternal > 0)) ||
                ((doISMhaploid == YES) && (genome == PATERNAL) && (allSites[site].numMutationsPaternal > 0)))
                {
                data[genome][cell][site] = data[genome][anccell][site];
                }
            else /* => there will be change */
                {
                if (data[genome][anccell][site] == 0)  /* checking all this might be excessive */
                    data[genome][cell][site] = 1;
                else if (data[genome][anccell][site] == 1)
                    {
                    fprintf (stderr, "\n\nERROR: site %d in genome %d of cell %d cannot mutate twice under the ISM model", site, genome, cell);
                    exit(-1);
                    }
                else
                    {
                    fprintf (stderr, "\n\nERROR: site %d in genome %d of cell %d has an unknow state %d under the ISM model", site, genome, anccell, data[genome][anccell][site]);
                    exit(-1);
                    }
					
                if (genome == MATERNAL)
                    allSites[site].numMutationsMaternal++;
                else if (genome == PATERNAL)
                    allSites[site].numMutationsPaternal++;
                allSites[site].numMutations++;
                numMU++;
                }
            }
            SimulateISMforSite (p->left, genome, site, doISMhaploid, seed);
            SimulateISMforSite (p->right, genome, site, doISMhaploid, seed);
        }
    }


/********************************** SimulateISMDNAforSite ***********************************/
/*	Simulates a ACGT mutation under an infinite sites model (ISM) for a given site. The branch
	where this mutation is placed is chosen according to its length.
    The reference (healthy) allele for each site will be determined by the nucleotide frequencies
    A mutation matrix will define the probability of changing to other nucleotides given
    the healthy alelle chose
 */
void SimulateISMDNAforSite (TreeNode *p, int genome, int site, int doISMhaploid, long int *seed)
	{
    static double	cumBranchLength, uniform, ran;
    int             j, cell, anccell, ancstate;
	
    if (p != NULL)
        {
        cell = p->label;
			
        if (p->isHealthyRoot == YES)
            {
            cumBranchLength = 0;
            uniform = RandomUniform(seed) * totalTreeLength;
            }
        else
            {
            anccell = p->anc->label;
            ancstate = data[genome][anccell][site];
            cumBranchLength += p->branchLength;
				
            if ((cumBranchLength < uniform) || /* => there will be no change */
                ((doISMhaploid == NO)  && (allSites[site].numMutations > 0))  ||
                ((doISMhaploid == YES) && (genome == MATERNAL) && (allSites[site].numMutationsMaternal > 0)) ||
                ((doISMhaploid == YES) && (genome == PATERNAL) && (allSites[site].numMutationsPaternal > 0)))
                {
                data[genome][cell][site] = ancstate;
                }
            else /* => there will be change */
                {
                ran = RandomUniform(seed) * cumMij[ancstate][3];
                for (j=0; j<4; j++)
                    {
                    if (ran <= cumMij[ancstate][j])
                        {
                        data[genome][cell][site] = j;
                        break;
                        }
					}
                if (genome == MATERNAL)
                    allSites[site].numMutationsMaternal++;
                else if (genome == PATERNAL)
                    allSites[site].numMutationsPaternal++;
                allSites[site].numMutations++;
                numMU++;
 				}
            }
        SimulateISMDNAforSite (p->left, genome, site, doISMhaploid, seed);
        SimulateISMDNAforSite (p->right, genome, site, doISMhaploid, seed);
		}
	}



/*************************** SimulateSignatureISM **********************************/
/*	Simulates genetic signatures under an infinite diploid trinucleotide model (ISM).

   //MARK: importantly, in the current implementation we assume that the context does not evolve

    Diploid here means that if site 33 mutates in the female genome,
    site 33 in the male genome will not mutate, and viceversa.

    The number of mutations will be distributed as:
 
    (1) a Poisson with parameter the sum of
	branch lengths over all sites  Times are already scaled in
	2N. We will force it to get just one mutation per site
 
	(2) A used-defined number of mutations (segregating sites in the ISM diploid)
 
*/

void SimulateSignatureISM (TreeNode *p, int genome, long int *seed)
    {
    int		i, numMutations, mutationsSoFar;
    double	totalBranchSum;
	int		newState, chosenSite;
	
	newState = -9;
	
	/* allocate memory to store the trinucleotide frequencies */
	triNucleotideMaternal = (TriNucStr *) calloc (64, sizeof(TriNucStr));
	if (!triNucleotideMaternal)
		{
		fprintf (stderr, "Could not allocate the triNucleotideMaternal vector\n");
		exit (-1);
		}
	for (i=0; i<64; i++)
		{
		triNucleotideMaternal[i].tempLength = 100;
		triNucleotideMaternal[i].numAvailablePositions = 0;
		triNucleotideMaternal[i].position = (int *) calloc(100,sizeof(int));
		}

	triNucleotidePaternal = (TriNucStr *) calloc (64, sizeof(TriNucStr));
	if (!triNucleotidePaternal)
		{
		fprintf (stderr, "Could not allocate the triNucleotidePaternal vector\n");
		exit (-1);
		}
	for (i=0; i<64; i++)
		{
		triNucleotidePaternal[i].tempLength = 100;
		triNucleotidePaternal[i].numAvailablePositions = 0;
		triNucleotidePaternal[i].position = (int *) calloc(100,sizeof(int));
		}

	/* count trinucleotide frequencies in the healthy root genome */
	CountTriNucFrequencies (data[genome][HEALTHY_ROOT], genome);
/*
    totalBranchSum = 0;
    for (site=0; site<numSites; site++)
        {
        allSites[site].branchSum = SumBranches (healthyRoot);
        totalBranchSum += allSites[site].branchSum;
        }
*/
	totalBranchSum = totalTreeLength * numSites;

    if (doSimulateFixedNumMutations == YES) /* conditioned on a specified number of mutations */
		{
        if (genome == MATERNAL)
			 numMutations = numSNVmaternal =  RandomBinomial (0.5, numFixedMutations, seed);
        else
            numMutations = numFixedMutations - numSNVmaternal;
        }
    else
        numMutations = RandomPoisson (totalBranchSum, seed); /* the number of mutations will be distributed as a Poisson with parameter totalBranchSum */
		
     /* if the number of ISM mutations is bigger than the number of ISM sites quit the program and warn the user about violation of ISM */
     if ((numISMmutations + numMutations) > numSites)
        {
        fprintf (stderr, "\n\nERROR: The diploid infinite sites model (ISM) has been violated. There will be");
        fprintf (stderr, "\nmore mutations (%d existing + %d proposed]) than available sites under this model (%d).", numISMmutations, numMutations, numDefaultModelSites);
        fprintf (stderr, "\nTry using a smaller mutation rate or effective population size");
        fprintf (stderr, "\n(i.e. smaller theta) or check the proportion of JC sites.");
        fprintf (stderr, "\nIf using an user tree, try with smaller branch lengths (including the transforming and healthy branches).\n");
        exit (-1);
        }

	numISMmutations += numMutations;

	#ifdef PRINT_TRIMUTCOUNTER
	if (genome == MATERNAL)
		{
		triMutationsCounter = (int *) calloc (96, sizeof(int));
		if (!triMutationsCounter)
			{
			fprintf (stderr, "Could not allocate the triMutationsCounter vector\n");
			exit (-1);
			}
		}
	#endif

    mutationsSoFar = 0;
    while (mutationsSoFar < numMutations)
        {
 		chosenSite = ChooseTrinucleotideSite(seed, &newState, genome);
		
		if (allSites[chosenSite].numMutations == 0) /* force ISM diploid */
			{
			SimulateSignatureISMforSite (p, genome, chosenSite, newState, seed);
			mutationsSoFar++;
			}
		else // chose another site as this site has already changed in the maternal chromosome
			{
			//fprintf (stderr, "\nsite %d has already changed", chosenSite+1);
			}
        }
	
	//fprintf (stderr, "\n\nnumISMmutations = %d  numMutations = %d   mutationsSoFar = %d (genome = %d)",numISMmutations, numMutations, mutationsSoFar, genome);
	
	#ifdef PRINT_TRIMUTCOUNTER
	int i, j, k;
	int doDeconstructSigs = YES;
	FILE *fpCounts, *fpRefGenome;

	sprintf(File,"%s/%s", resultsDir, "trimutcounts");
	if ((fpCounts = fopen(File, "w")) == NULL)
		{
		fprintf (stderr, "Can't open \"%s\"\n", File);
		exit(-1);
		}

	if (doDeconstructSigs == NO) /* counter for SomaticSignatures */
		{
		for (i=0; i<96; i++)
			fprintf (fpCounts, "%d\n", triMutationsCounter[i]);
		}
	else /* counter for deconstructSigs */
		{
		for (j=0; j<6; j++)
			for (i=0; i<4; i++)
				for (k=0; k<4; k++)
					fprintf (fpCounts, "%d\n", triMutationsCounter[trimut(i,j,k)]);
		}
	fclose(fpCounts);

	/* print reference genome as genome context for SomaticSignatures */
	sprintf(File,"%s/%s", resultsDir, "contextGenome");
	if ((fpRefGenome = fopen(File, "w")) == NULL)
		{
		fprintf (stderr, "Can't open \"%s\"\n", File);
		exit(-1);
		}
	fprintf (fpRefGenome, ">1\n");
	for (i=0; i<numSites; i++)
		fprintf (fpRefGenome, "%c", WhichNuc(allSites[i].referenceAllele));
	fclose(fpRefGenome);
	#endif
		
    }


/********************************** CountTriNucFrequencies ***********************************/
/*
	 Counts the trinucleotide frequencies in the given genome
*/

void CountTriNucFrequencies (int *genome_array, int genome)
	{
	int i, n1, n2, n3;
	TriNucStr *triNucleotide;

	if (genome == MATERNAL)
		triNucleotide = triNucleotideMaternal;
	else
		triNucleotide = triNucleotidePaternal;

	for (i=0; i<numSites-2; i++)
		{
		n1 = genome_array[i];
		n2 = genome_array[i+1];
		n3 = genome_array[i+2];
		RecordTriNucObservation(&triNucleotide[trinuc(n1,n2,n3)], i);
		}

	#ifdef PRINT_TRINUC_GENOME
		if (genome == MATERNAL)
			{
			fprintf (stderr, "\n\nSimulated trinuc maternal genome: ");
			int rest;
			for (i=0; i<numSites; i++)
				fprintf (stderr, "%c", WhichNuc(allSites[i].referenceAllele));
			
			fprintf (stderr, "\n\nObserved trinuc freqs: ");
			for (i=0; i<64; i++)
				{
				n1 = i/16;
				rest = i%16;
				n2 = rest/4;
				n3 = rest%4;
				fprintf (stderr, "\n%2d %c%c%c [%d] = ", i+1, WhichNuc(n1), WhichNuc(n2), WhichNuc(n3), triNucleotide[i].numAvailablePositions);
				for (int j=0; j<triNucleotide[i].numAvailablePositions; j++)
					fprintf (stderr, "%d ",triNucleotide[i].position[j]+1);
				 }
			fprintf (stderr, "\n\n");
			}
	#endif
	}


/********************************** ChooseTrinucleotideSite ***********************************/
/*
Find a valid trinucleotide site to introduce a mutation
The trinucleotide signatures are observed mutations counts, so we should choose sites according to that.
We cannot just choose at random

Given signatures are referred by the pyrimidine of the mutated Watson–Crick base pair, counts for reverse complements
trinucleotide changes are indistinguishable. For example:
 
	5-AGG-3 => 5-AAG-3 is indistinguishable from
	3-TCC-5 => 3-TTC-5 so is equivalent to
	5-CCT-3 => 5-CTT-3
*/

int ChooseTrinucleotideSite (long int *seed, int *newState, int genome)
	{
	int	n1, n2, n3, r1, r2, r3, rest, index, trindex, trials;
	int	chosenSignature, chosenTriMutation, chosenPosition;
	int signID;
	TriNucStr *triNucleotide;
	
	if (genome == MATERNAL)
		triNucleotide = triNucleotideMaternal;
	else
		triNucleotide = triNucleotidePaternal;

	/* 1: Choose a mutation signature */
	signID = ChooseUniformState(signatureWeight, seed);
	chosenSignature = signatureID[signID];
		
	/* 2: Select a mutation out of the 96 possible ones for the choosen signature */
	//chosenTriMutation = ChooseTrinucMutation(geneticSignature[chosenSignature], seed);
	chosenTriMutation = ChooseUniformState(signatureProbs[chosenSignature-1], seed);

	/* 3: Define in which of the 64 trinucleotide classes this mutation will take place */
	r1 = chosenTriMutation/24;
	rest = chosenTriMutation%24;
	if (rest/4 < 3)
		r2 = C;
	else
		r2 = T;
	r3 = rest%4;

	if (RandomUniform(seed) > 0.5) //plus strand
		{
		n1 = r1;
		n2 = r2;
		n3 = r3;
		*newState = targetTriChange[rest/4];
		//fprintf (stderr, "\n\nplus strand");
		}
	else // minus strand
		{
		n1 = complementBase[r3];
		n2 = complementBase[r2];
		n3 = complementBase[r1];
		*newState = complementBase[targetTriChange[rest/4]];
		//fprintf (stderr, "\n\ncomplement reverse");
		}
		
	trindex = trinuc(n1,n2,n3);

	#ifdef PRINT_TRIMUTATIONS
	fprintf (stderr, "\n\nchosen triMutation (id %d) = %c%c%c",chosenTriMutation, WhichNuc(r1), WhichNuc(r2), WhichNuc(r3));
	fprintf (stderr, "\n%c%c%c [%d] => ", WhichNuc(n1), WhichNuc(n2), WhichNuc(n3), triNucleotide[trindex].numAvailablePositions);
	fprintf (stderr, "%c%c%c", WhichNuc(n1), WhichNuc(*newState), WhichNuc(n3));
	fprintf (stderr, " (%c>%c)", WhichNuc(n2), WhichNuc(*newState));
	#endif

	trials = 0;
	while (triNucleotide[trindex].numAvailablePositions == 0)
		{
		/* 1: Choose a mutation signature */
		signID = ChooseUniformState(signatureWeight, seed);
		chosenSignature = signatureID[signID];
		
		/* 2: Select a mutation out of the 96 possible ones for the choosen signature */
		chosenTriMutation = ChooseUniformState(signatureProbs[chosenSignature-1], seed);

		/* 3: Define the trinucleotide class where this mutation would take place */
		r1 = chosenTriMutation/24;
		rest = chosenTriMutation%24;
		if (rest/4 < 3)
			r2 = C;
		else
			r2 = T;
		r3 = rest%4;

		if (RandomUniform(seed) > 0.5) //plus strand
			{
			n1 = r1;
			n2 = r2;
			n3 = r3;
			*newState = targetTriChange[rest/4];
			//fprintf (stderr, "\nplus strand");
			}
		else // minus strand
			{
			n1 = complementBase[r3];
			n2 = complementBase[r2];
			n3 = complementBase[r1];
			*newState = complementBase[targetTriChange[rest/4]];
			//fprintf (stderr, "\ncomplement reverse");
			}
			
		trindex = trinuc(n1,n2,n3);
		
		#ifdef PRINT_TRIMUTATIONS
		fprintf (stderr, "\nnew chosen triMutation (id %d) = %c%c%c", chosenTriMutation, WhichNuc(r1), WhichNuc(r2), WhichNuc(r3));
		fprintf (stderr, "\n%c%c%c [%d] => ", WhichNuc(n1), WhichNuc(n2), WhichNuc(n3), triNucleotide[trindex].numAvailablePositions);
		fprintf (stderr, "%c%c%c", WhichNuc(n1), WhichNuc(*newState), WhichNuc(n3));
		fprintf (stderr, " (%c>%c)", WhichNuc(n2), WhichNuc(*newState));
		#endif

		if (trials++ > 100*numSites)
			{
			fprintf (stderr, "\n\nERROR: after %d trials cannot find an non-mutated trinuc site, increase the number of sites simulated\n\n",100*numSites);
			exit(-1);
			}
		}

	#ifdef PRINT_TRIMUTATIONS
		fprintf (stderr, "\nfinally chosen trimutation index = %d", chosenTriMutation);
	#endif
	
	/* 3: Select the genomic position where this mutation will take place */
	/* note we assume for convenience that the context is constant along the phylogeny */
	if (triNucleotide[trindex].numAvailablePositions == 1)
		index = 0;
	else if (triNucleotide[trindex].numAvailablePositions > 1)
		index = RandomUniformTo(triNucleotide[trindex].numAvailablePositions, seed);
	else
		{
		fprintf (stderr, "\n\nWe have a problem, no available trinucleotides\n");
		exit(-1);
		}

	chosenPosition = triNucleotide[trindex].position[index] + 1; 
	#ifdef PRINT_TRIMUTATIONS
		fprintf (stderr, "\nposition selected = %d", chosenPosition+1);
	#endif
	
	#ifdef PRINT_TRIMUTCOUNTER
	if (allSites[chosenPosition].numMutations == 0)
		triMutationsCounter[chosenTriMutation]++;
	#endif
	
	/* 4: remove the chosen position from the available list */
	RemoveTriNucObservation(&triNucleotide[trindex], index);

	return chosenPosition;
	}




/********************************** SimulateSignatureISMforSite ***********************************/
/*	Simulates genetic signatures under an infinite sites model (ISM) for a given site. The branch
	where this mutation is placed is chosen according to its length.
 */

void SimulateSignatureISMforSite (TreeNode *p, int genome, int site, int newState, long int *seed)
	{
    static double	cumBranchLength, uniform;
    int             cell, anccell, ancState;
	
    if (p != NULL)
        {
        cell = p->label;
			
        if (p->isHealthyRoot == YES)
            {
            cumBranchLength = 0;
             uniform = RandomUniform(seed) * totalTreeLength;
			//fprintf (stderr, "\nsimulating at root of position %d", site+1);
            }
        else
            {
            anccell = p->anc->label;
            ancState = data[genome][anccell][site];
            cumBranchLength += p->branchLength;
				
            if (cumBranchLength < uniform || allSites[site].numMutations > 0)/* => there will be no change; forcing ISM diploid */
                {
                data[genome][cell][site] = ancState;
 				//fprintf (stderr, "\ncopying ancstate %c", WhichNuc(ancstate));
               }
            else /* => there will be change */
                {
				data[genome][cell][site] = newState;
                if (genome == MATERNAL)
                    allSites[site].numMutationsMaternal++;
                else if (genome == PATERNAL)
                    allSites[site].numMutationsPaternal++;
                allSites[site].numMutations++;
                numMU++;
				//fprintf (stderr, "\nchange at site %d from %c to %c  (numMU = %d)", site+1, WhichNuc(ancState), WhichNuc(newState), numMU);
 				}
            }
        SimulateSignatureISMforSite (p->left, genome, site, newState, seed);
        SimulateSignatureISMforSite (p->right, genome, site, newState, seed);
		}
	}


/********************************** RemoveTriNucObservation ***********************************/
/* Remove the occurrence of a trinucleotide along the genome */
void RemoveTriNucObservation(TriNucStr *trin, int site)
	{
	int		i;
	
	for (i=site; i<trin->tempLength-1; i++)
		 trin->position[i] = trin->position[i+1];
	
	trin->tempLength = --trin->numAvailablePositions;
	
	if (trin->numAvailablePositions > 0 && trin->numAvailablePositions % 100 == 0)
		{
		//fprintf (stderr, "\nreallocating now (available positions = %d)",trin->numAvailablePositions);
		trin->position = (int *)realloc(trin->position, trin->tempLength * sizeof(int));
		if (!trin->position)
			{
			fprintf (stderr, "Could not reallocate trinucleotide position array when removing\n");
			exit (-1);
			}
		}
 /*	else if (trin->numAvailablePositions == 0)
		{
		fprintf (stderr, "\nthis class has zero available positions now");
		}
*/

/* NOTE:  we could put this nucleotide class signature to zero, as it cannot mutate again, always considering the strand
		  But it does not seem necessary 240518 */
  	}



/************************************* SimulateMk2 **********************************************/
/* Simulates the nucleotide substitution process under the Mk2 model (see Lewis 2001),
    also called Cavender-Farris-Neyman CFN  model or Jukes-Cantor (1969) model for two alleles */

void SimulateMk2 (TreeNode *p, int genome, long int *seed)
    {
    int     i;
		
    for (i=0; i<numAltModelSites; i++)
        SimulateMk2forSite (p, genome, AltModelSites[i], seed);
	}

/************************************* SimulateMk2ForSite ***************************************/
/* Simulates the nucleotide substitution process for a given site under Mk2 model (see Lewis 2001)
 with equal rates. 0 is the reference (healthy) allele */

void SimulateMk2forSite (TreeNode *p, int genome, int site, long int *seed)
    {
    double	probOfChange, uniform, branchLength;
    int     cell, anccell;
		
    if (p != NULL)
		{
		if (p->isHealthyRoot == NO)
			{
			cell = p->label;
			anccell = p->anc->label;
			
			if (doUserTree == YES)
				branchLength = p->branchLength;
			else
				{
				if (rateVarAmongSites == YES)
					branchLength = altModelMutationRate * p->length * allSites[site].rateMultiplier;
				else
					branchLength = altModelMutationRate * p->length;
				}
	
			probOfChange = 0.5 - 0.5 * exp (-2.0 * branchLength);

			uniform = RandomUniform(seed);
			if (uniform >= probOfChange) /* => no change */
				data[genome][cell][site] = data[genome][anccell][site];
			else /* => there will be change */
				{
				if (data[genome][anccell][site] == 0)
					data[genome][cell][site] = 1;
				else
					data[genome][cell][site] = 0;
					
				if (genome == MATERNAL)
					allSites[site].numMutationsMaternal++;
				else if (genome == PATERNAL)
					allSites[site].numMutationsPaternal++;
				allSites[site].numMutations++;
				numMU++;
				}
			}
		SimulateMk2forSite (p->left,  genome, site, seed);
		SimulateMk2forSite (p->right, genome, site, seed);
		}
	}


/************************************* SimulateFiniteDNA **********************************************/
/* Simulates the nucleotide substitution process under a 4-state Markov model including JC, HKY, GTR and GTRnr */
/* Note that beta is set such that mean substitution rate will be 1.
 	E.g., for J-C model, beta=4/3, where 12(1/4)(1/4)(4/3) = 1.      */

void SimulateFiniteDNA (TreeNode *p, int genome, long int *seed)
    {
    int     i, j;
		
	if (doJC == YES)
		{
		beta = 4./3;
		}
	else if (doHKY == YES)
        {
		freqR = freq[A] + freq[G];
		freqY = freq[C] + freq[T];
		freqAG = freq[A] * freq[G];
		freqCT = freq[C] * freq[T];
		kappa = (titv*freqR*freqY)/(freqAG+freqCT);
		beta = 0.5 / (freqR*freqY + kappa*(freqAG+freqCT));
		}
    else if (doGTR == YES || doGTnR == YES)
		{
		for (i=0; i<4; i++)
			for (j=0; j<4; j++)
				Qij[i*4+j] = Mij[i][j] / Mij[2][3] * freq[j];
		mr=0;
		for (i=0; i<4; i++)
			{
			Qij[i*4+i]=0;
			Qij[i*4+i]=-(Qij[i*4]+Qij[i*4+1]+Qij[i*4+2]+Qij[i*4+3]);
			mr-=freq[i]*Qij[i*4+i];
			}
		EigenREV(Root, Cijk);
		}
		
    for (i=0; i<numAltModelSites; i++)
        SimulateFiniteDNAforSite (p, genome, AltModelSites[i], seed);
	}



/************************************* SimulateFiniteDNAforSite **********************************************/
/* Simulates JC, HKY, GTR or GTRnr for a given site */
void SimulateFiniteDNAforSite (TreeNode *p, int genome, int site, long int *seed)
	{
    double	branchLength, Pij[4][4];
    int     cell, anccell, ancstate, newstate;
		
    if (p != NULL)
        {
		if (p->isHealthyRoot == NO)
			{
			cell = p->label;
			anccell = p->anc->label;
			ancstate = data[genome][anccell][site];
	
			if (doUserTree == YES)
				branchLength = p->branchLength;
			else
				{
				if (rateVarAmongSites == YES)
					branchLength = altModelMutationRate * p->length * allSites[site].rateMultiplier;
				else
					branchLength = altModelMutationRate * p->length;
				}
				
			FillSubstitutionMatrix (Pij, branchLength);
			newstate = data[genome][cell][site] = ChooseUniformState (Pij[ancstate], seed);
			
			if (newstate != ancstate)
				{
				if (genome == MATERNAL)
					allSites[site].numMutationsMaternal++;
				else if (genome == PATERNAL)
					allSites[site].numMutationsPaternal++;
				allSites[site].numMutations++;
				numMU++;
				}
			}
		SimulateFiniteDNAforSite (p->left,  genome, site, seed);
        SimulateFiniteDNAforSite (p->right, genome, site, seed);
        }
 	}


/********************* FillSubstitutionMatrix **********************/
/* Sets the apropriate model of nucleotide substitution   */
void FillSubstitutionMatrix (double ch_prob[4][4], double branchLength)
	{
	int i, j;
	
	if (branchLength<1e-6)
		{
		for (i=0; i<4; i++)
			{
			for (j=0; j<4; j++)
				{
				if (i == j)
					ch_prob[i][j] = 1.0;
				else
					ch_prob[i][j] = 0.0;
				}
			}
		}
    else if (doJC == YES)
        JCmodel (ch_prob, branchLength);
    else if (doHKY == YES)
        HKYmodel (ch_prob, branchLength);
    else if (doGTR == YES)
		GTRmodel (ch_prob, branchLength);
    else if (doGTnR == YES)
        GTRmodel (ch_prob, branchLength);
	}



/*********************************** JC **************************************/
/*	JC performs Jukes-Cantor 69 correction */
/* Note that beta was set such that mean substitution rate will be 1.
  for the JC model, beta=4/3, where 12(1/4)(1/4)(4/3) = 1.      */
void JCmodel (double Pij[4][4], double branchLength)
	{
	int i, j;
	
	for (i=0; i<4; i++)
		{
		for (j=0; j<4; j++)
			{
			if (i == j)
				Pij[i][j] = 0.25 + 0.75*exp(beta*-branchLength);
			else
				Pij[i][j] = 0.25 - 0.25*exp(beta*-branchLength);
			}
		}
	}

/*********************************** HKY **************************************/
/*	HKY performs Hasegawa-Kishino-Yano 85 correction */

void HKYmodel (double Pij[4][4], double branchLength)
	{
	int			i, j;
	double		AA, t, PIj;

	t = branchLength;

	for (i=0; i<4; i++)
		{
		for (j=0; j<4; j++)
			{
			if (j == A || j == G)	/* purine */
				PIj = freqR;
			else
				PIj = freqY; /* pyrimidine */
				
			AA = 1 + PIj*(kappa-1.0);
			
			if (i==j)
				Pij[i][j] = freq[j] + freq[j]*(1/PIj - 1)*exp(-beta*t) + ((PIj-freq[j])/PIj)*exp(-beta*t*AA);
			else if ((i==A && j==G) || (i==C && j==T) || (i==G && j==A) || (i==T && j==C)) /* transition */
				Pij[i][j] = freq[j] + freq[j]*(1/PIj - 1)*exp(-beta*t) - (freq[j]/PIj)*exp(-beta*t*AA);
			else /* transversion */
				Pij[i][j] = freq[j]*(1-exp(-beta*t));
			}
		}
	}


/*************** GTR **********************/
void GTRmodel (double Pij[4][4], double branchLength)
	{
	int 	i, j, k;
	double	t, expt[4];

	t = branchLength;

	/* P(t)ij = SUM Cijk * exp{Root*t} */
	for (k=1; k<4; k++)
		expt[k]=exp(t*Root[k]);
	for (i=0; i<4; i++)
		for (j=0; j<4; j++)
			{
			Pij[i][j]=Cijk[i*4*4+j*4+0];
			for (k=1; k<4; k++)
				Pij[i][j]+=Cijk[i*4*4+j*4+k]*expt[k];
			}
	}



/********************************** EvolveCNLOHonTree ***********************************/
 /*	Simulates point copy-number neutral loss of heterozygosity events (CN-LOH) under an infinite haploid sites model (ISM).
    Haploid here means that site 33 in the maternal genome and site 33 in the paternal
    genome will be considered different sites, so both can mutate.
 
   This is a very simple mode, where CN_LOH evebts are simulated after the sequences have been evolved
   (i.e., SNVs are already in place)
 
    The number of point CNLOH will be distributed as a Poisson with parameter the sum of
	branch lengths (node length * deletion rate) over all sites.
	
	We will force it to get at most one CNLOH per site
*/
void EvolveCNLOHonTree (TreeNode *p, int genome, long int *seed)
	{
    int		i, trials, numCNLOH, CNLOHsoFar;
    double	totalCNLOHbranchSum;

	totalCNLOHbranchSum = numSites * totalTreeLength / mutationRate * CNLOHrate;
	numCNLOH = RandomPoisson (totalCNLOHbranchSum, seed); /* the number of CNLOH events will be distributed as a Poisson with parameter totalCNLOHbranchSum */
	
     /* if the number of CNLOH events is bigger than the number of sites quit the program and warn the user about violation of ISM */
     if (numISMCNLOH + numCNLOH > (2*numSites))
        {
        fprintf (stderr, "\n\nERROR: The haploid infinite sites model (ISM) for CNLOH events has been violated. There will be");
        fprintf (stderr, "\nmore CNLOH events (%d existing + %d proposed]) than available sites under this model (%d).", numISMCNLOH, numCNLOH, 2*numSites);
        fprintf (stderr, "\nTry using a smaller deletion rate.");
        fprintf (stderr, "\nIf using a user tree, try introducing smaller branch lengths.\n");
        exit (-1);
        }
	
    numISMCNLOH += numCNLOH;
    trials = 0;
    CNLOHsoFar = 0;
    while (CNLOHsoFar < numCNLOH)
        {
        i = RandomUniformTo(numSites, seed);  /* choose a site at random */
        while ((genome == MATERNAL && allSites[i].numCNLOHmaternal != 0) ||
               (genome == PATERNAL && allSites[i].numCNLOHpaternal != 0))
            {
            i = RandomUniformTo(numSites, seed);
            if (trials++ > 100*numSites)
                {
                fprintf (stderr, "\n\n ERROR: after %d trials cannot find a site without CNLOH",100*numSites);
                fprintf (stderr, "\nCNLOH = %d  CNLOH so far = %d\n\n",numCNLOH, CNLOHsoFar);
                for (i=0; i<numSites; i++)
                    fprintf (stderr, "\nsite %d has %d CNLOHs",i+1, allSites[i].numCNLOH);
                 fprintf (stderr, "\n");
                 exit(-1);
                }
            }
	
		SimulateCNLOHforSite (p, genome, i, seed);
			
        CNLOHsoFar++;
			
        #ifdef MYDEBUG
            fprintf (stderr, "\nCNLOH = %d   nCNLOH so far = %d\n",numCNLOH, CNLOHsoFar);
            if (allSites[i]].numDeletions > 1)
                {
                fprintf (stderr, "\n\n ERROR: %d CNLOH in site %d",allSites[i].numCNLOH, i+1);
                exit(-1);
                }
            for (i=0; i<numSites; i++)
                fprintf (stderr, "%2d[%d] ",i+1, allSites[i].numCNLOH);
        #endif
			
        }
	}


/********************************** SimulateCNLOHforSite ***********************************/
/*	Simulates a CNLOH event infinite sites model (ISM) for a given site. The branch
	where this mutation is placed is chosen according to its length.
  */
void SimulateCNLOHforSite (TreeNode *p, int genome, int site, long int *seed)
    {
    static double	cumCNLOHbranchLength, uniform;
    int             cell, anccell;
		
	if (p != NULL)
        {
        cell = p->label;
			
        if (p->isHealthyRoot == YES)
            {
            cumCNLOHbranchLength = 0;
            uniform = RandomUniform(seed) * totalTreeLength / mutationRate * CNLOHrate;;
            }
        else
            {
            anccell = p->anc->label;
            cumCNLOHbranchLength += p->branchLength / mutationRate * CNLOHrate;
				
             if ((cumCNLOHbranchLength < uniform) || /* => there will be no change at this branch */
                ((genome == MATERNAL) && (allSites[site].numCNLOHmaternal > 0)) ||
                ((genome == PATERNAL) && (allSites[site].numCNLOHpaternal > 0)))
                {
                data[genome][cell][site] = data[genome][anccell][site];
                }
			else /* => there will be a CN_LOH event */
                {
				if (genome == MATERNAL)
					{
					data[genome][cell][site] = data[PATERNAL][cell][site];
					allSites[site].numCNLOHmaternal++;
					}
				else
					{
					data[genome][cell][site] = data[MATERNAL][cell][site];
					allSites[site].numCNLOHpaternal++;
					}

                allSites[site].numCNLOH++;
                numCNLOH++;
                }
            }
            SimulateCNLOHforSite (p->left, genome, site, seed);
            SimulateCNLOHforSite (p->right, genome, site, seed);
        }
    }



/********************************** EvolveDeletionsOnTree ***********************************/
 /*	Simulates point deletions under an infinite haploid sites model (ISM).
    Haploid here means that site 33 in the maternal genome and site 33 in the paternal
    genome will be considered different sites, so both can mutate.
 
   This is a very simple mode, where deletions are simulated after the sequences have been evolved
   (i.e., SNVs are already in place)
 
    The number of point deletions will be distributed as a Poisson with parameter the sum of
	branch lengths (node length * deletion rate) over all sites.

	We will force it to get at most one deletion per site
*/
void EvolveDeletionsOnTree (TreeNode *p, int genome, long int *seed)
	{
    int		i, trials, numDeletions, deletionsSoFar;
    double	totalDeletionBranchSum;

    /*
    totalDeletionBranchSum = 0;
    for (i=0; i<numSites; i++)
        {
        allSites[i].deletionBranchSum = SumBranches (healthyRoot) / mutationRate * deletionRate;;
        totalDeletionBranchSum += allSites[i].deletionBranchSum;
        }
	*/
	
	totalDeletionBranchSum = numSites * totalTreeLength / mutationRate * deletionRate;
	numDeletions = RandomPoisson (totalDeletionBranchSum, seed); /* the number of deletions will be distributed as a Poisson with parameter totalDeletionBranchSum */
	
     /* if the number of deletions is bigger than the number of sites quit the program and warn the user about violation of ISM */
     if (numISMdeletions + numDeletions > (2*numSites))
        {
        fprintf (stderr, "\n\nERROR: The haploid infinite sites model (ISM) for deletions has been violated. There will be");
        fprintf (stderr, "\nmore deletions (%d existing + %d proposed]) than available sites under this model (%d).", numISMdeletions, numDeletions, 2*numSites);
        fprintf (stderr, "\nTry using a smaller deletion rate.");
        fprintf (stderr, "\nIf using a user tree, try introducing smaller branch lengths.\n");
        exit (-1);
        }
	
    numISMdeletions += numDeletions;
    trials = 0;
    deletionsSoFar = 0;
    while (deletionsSoFar < numDeletions)
        {
        i = RandomUniformTo(numSites, seed);  /* choose a site at random */
        while ((genome == MATERNAL && allSites[i].numDeletionsMaternal != 0) ||
               (genome == PATERNAL && allSites[i].numDeletionsPaternal != 0))
            {
            i = RandomUniformTo(numSites, seed);
            if (trials++ > 100*numSites)
                {
                fprintf (stderr, "\n\n ERROR: after %d trials cannot find an undeleted site",100*numSites);
                fprintf (stderr, "\ndeletions = %d   deletions so far = %d\n\n",numDeletions, deletionsSoFar);
                for (i=0; i<numSites; i++)
                    fprintf (stderr, "\nsite %d has %d deletions",i+1, allSites[i].numDeletions);
                 fprintf (stderr, "\n");
                 exit(-1);
                }
            }
	
		SimulateDeletionforSite (p, genome, i, seed);
			
        deletionsSoFar++;
			
        #ifdef MYDEBUG
            fprintf (stderr, "\deletions = %d   deletions so far = %d\n",numDeletions, deletionsSoFar);
            if (allSites[i]].numDeletions > 1)
                {
                fprintf (stderr, "\n\n ERROR: %d deletions in site %d",allSites[i].numDeletions, i+1);
                exit(-1);
                }
            for (i=0; i<numSites; i++)
                fprintf (stderr, "%2d[%d] ",i+1, allSites[i].numDeletions);
        #endif
			
        }
	}

/********************************** SimulateDeletionforSite ***********************************/
/*	Simulates a point deletion infinite sites model (ISM) for a given site. The branch
	where this mutation is placed is chosen according to its length.
  */
void SimulateDeletionforSite (TreeNode *p, int genome, int site, long int *seed)
    {
    static double	cumDeletionBranchLength, uniform;
    int             cell, anccell;
		
	if (p != NULL)
        {
        cell = p->label;
			
        if (p->isHealthyRoot == YES)
            {
            cumDeletionBranchLength = 0;
//            uniform = RandomUniform(seed) * allSites[site].deletionBranchSum;
             uniform = RandomUniform(seed) * totalTreeLength / mutationRate * deletionRate;;
           }
        else
            {
            anccell = p->anc->label;
            cumDeletionBranchLength += p->branchLength / mutationRate * deletionRate;
				
             if ((cumDeletionBranchLength < uniform) || /* => there will be no change at this branch */
                ((genome == MATERNAL) && (allSites[site].numDeletionsMaternal > 0)) ||
                ((genome == PATERNAL) && (allSites[site].numDeletionsPaternal > 0)))
                {
                data[genome][cell][site] = data[genome][anccell][site];
                }
			else /* => there will be a deletion */
                {
                if (data[genome][anccell][site] != DELETION)  /* checking all this might be excessive */
                    data[genome][cell][site] = DELETION;
                else if (data[genome][anccell][site] == DELETION)
                    {
                    fprintf (stderr, "\n\n ERROR: site %d in genome %d of cell %d cannot be deleted twice under the ISM model\n", site, genome, cell);
                    exit(-1);
                    }
                else
                    {
                    fprintf (stderr, "\n\n ERROR: site %d in genome %d of cell %d has an unknow state %d under the ISM model\n", site, genome, anccell, data[genome][anccell][site]);
                    exit(-1);
                    }
					
                if (genome == MATERNAL)
                    allSites[site].numDeletionsMaternal++;
                else if (genome == PATERNAL)
                    allSites[site].numDeletionsPaternal++;
                allSites[site].numDeletions++;
                numDEL++;
                }
            }
            SimulateDeletionforSite (p->left, genome, site, seed);
            SimulateDeletionforSite (p->right, genome, site, seed);
        }
    }




/********************* AllelicDropout  ************************/
/*
 Remove alleles from single chromosomes at a given ADO rate per genotype

 We assume that ADOrate is the product of the allele dropout rate as follows:
 A = genotype ADO; a = allele error
 A = 2a - a^2
 a = 1 - sqrt(1-A)
*/

void AllelicDropout (long int *seed)
	{
	int i,j;
	double alleleADOrate;
	
	alleleADOrate = 1.0 - sqrt (1.0 -  ADOrate);

	for (i=0; i<numCells+1; i++)
		{
		for (j=0; j<numSites; j++)
			{
			if (data[MATERNAL][i][j] != DELETION)
				{
				if (RandomUniform(seed) < alleleADOrate)
					{
					data[MATERNAL][i][j] = ADO;
					allSites[i].hasADO = YES;
					}
				}
			if (data[PATERNAL][i][j] != DELETION)
				{
				if (RandomUniform(seed) < alleleADOrate)
					{
					data[PATERNAL][i][j] = ADO;
					allSites[i].hasADO = YES;
					}
				}
			}
		}
	}



/********************* GenotypeError  ************************/
/*
 Introduce false positive errors directly in the genotypes
 (as opposed to errors just in the reads) (note Eij=0 when i=j)
 
 We assume that the genotype error is the product of the allele error as follows:
 G = genotype error; a = allele error
 G = 2a - a^2
 a = 1 - sqrt(1-E)
*/

void GenotypeError (long int *seed)
	{
	int i,j;
	double alleleError;
	
	alleleError = 1.0 - sqrt (1.0 -  genotypingError);
	
	if (alphabet == DNA)
		{
		for (i=0; i<numCells+1; i++)
			{
			for (j=0; j<numSites; j++)
				{
				if (data[MATERNAL][i][j] != ADO && data[MATERNAL][i][j] != DELETION && RandomUniform(seed) < alleleError)
					{
					data[MATERNAL][i][j] = ChooseUniformState (Eij[data[MATERNAL][i][j]], seed);
					allSites[j].hasGenotypeError = YES;
					}
				if (data[PATERNAL][i][j] != ADO && data[PATERNAL][i][j] != DELETION && RandomUniform(seed) < alleleError)
					{
					data[PATERNAL][i][j] = ChooseUniformState (Eij[data[PATERNAL][i][j]], seed);
					allSites[j].hasGenotypeError = YES;
					}
				}
			}
		}
	else  /* for binary data */
		{
		for (i=0; i<numCells+1; i++)
			{
			for (j=0; j<numSites; j++)
				{
				if (data[MATERNAL][i][j] != ADO && data[MATERNAL][i][j] != DELETION && RandomUniform(seed) < alleleError)
					{
					if (data[MATERNAL][i][j] == 0)
						data[MATERNAL][i][j] = 1;
					else if (data[MATERNAL][i][j] == 1)
						data[MATERNAL][i][j] = 0;
					allSites[j].hasGenotypeError = YES;
					}
				if (data[PATERNAL][i][j] != ADO && data[PATERNAL][i][j] != DELETION && RandomUniform(seed) < alleleError)
					{
					if (data[PATERNAL][i][j] == 0)
						data[PATERNAL][i][j] = 1;
					else if (data[PATERNAL][i][j] == 1)
						data[PATERNAL][i][j] = 0;
					allSites[j].hasGenotypeError = YES;
					}
				}
			}
		}
	}



/************************* CountTrueVariants  ************************/
/* Count number of variants in a given data set, assuming there is no ADO (yet).
   A variant is defined as no-deletion genotype different from the healthy root genotype */

int CountTrueVariants () 
	{
    int		cell, site;
    int		nVariants = 0;
		
    for (site=0; site<numSites; site++)
        {
		for (cell=0; cell<numCells+1; cell++)
			{
			if ((data[MATERNAL][cell][site] != DELETION && data[MATERNAL][cell][site] != data[MATERNAL][HEALTHY_ROOT][site]) 	||
			    (data[PATERNAL][cell][site] != DELETION && data[PATERNAL][cell][site] != data[PATERNAL][HEALTHY_ROOT][site]))
                {
				allSites[site].isVariant = YES;
				variantSites[nVariants++] = site;
  				break;
               }
			}
		}
    return nVariants;
	}

/************************* CountAlleles  ************************/
/* Identify reference and alternate alleles plus SNVs, in ingroup plus outgroup cells
   A SNV is defined as no-deletion genotype different from the reference genotype */

int CountAlleles ()
	{
	int		numAltAlleles, cell, site;
	int		countA, countC, countG, countT, countADO, countDEL, nSNVs, countCellswithData;
	
	nSNVs = 0;
	
	for (site=0; site<numSites; site++)
		{
		countA = countC = countG = countT = countADO = countDEL = countCellswithData = 0;
		for (cell=0; cell<numCells+1; cell++)
			{
			if (data[MATERNAL][cell][site] == A)
				countA++;
			else if (data[MATERNAL][cell][site] == C)
				countC++;
			else if (data[MATERNAL][cell][site] == G)
				countG++;
			else if (data[MATERNAL][cell][site] == T)
				countT++;
			else if (data[MATERNAL][cell][site] == ADO)
				countADO++;
			else if (data[MATERNAL][cell][site] == DELETION)
				countDEL++;

			if (data[PATERNAL][cell][site] == A)
				countA++;
			else if (data[PATERNAL][cell][site] == C)
				countC++;
			else if (data[PATERNAL][cell][site] == G)
				countG++;
			else if (data[PATERNAL][cell][site] == T)
				countT++;
			else if (data[PATERNAL][cell][site] == ADO)
				countADO++;
			else if (data[PATERNAL][cell][site] == DELETION)
				countDEL++;
			
			if ((data[MATERNAL][cell][site] != ADO && data[MATERNAL][cell][site] != DELETION) ||
			    (data[PATERNAL][cell][site] != ADO && data[PATERNAL][cell][site] != DELETION))
			    countCellswithData++;
			}

		allSites[site].countA = countA;
		allSites[site].countC = countC;
		allSites[site].countG = countG;
		allSites[site].countT = countT;
		allSites[site].countACGT = countA + countC + countG + countT;
		allSites[site].countDropped = countADO + countDEL;
		allSites[site].countCellswithData = countCellswithData;

		/* count number of alternate alleles, ignoring ADO or DELETION */
		numAltAlleles = 0;
		if (countA > 0 && allSites[site].referenceAllele != A)
			allSites[site].alternateAlleles[numAltAlleles++] = A;
		if (countC > 0 && allSites[site].referenceAllele != C)
			allSites[site].alternateAlleles[numAltAlleles++] = C;
		if (countG > 0 && allSites[site].referenceAllele != G)
			allSites[site].alternateAlleles[numAltAlleles++] = G;
		if (countT > 0 && allSites[site].referenceAllele != T)
			allSites[site].alternateAlleles[numAltAlleles++] = T;
		allSites[site].numAltAlleles = numAltAlleles;
		
		/* find out whether this site is a SNV */
		if (numAltAlleles > 0)
			{
			allSites[site].isSNV = YES;
			SNVsites[nSNVs++] = site;
			}
		}
	return nSNVs;
	}


/*************************** SumBranchlengths **********************************/
/* Returns the sum of the branch lengths for a given tree */
double SumBranches (TreeNode *p)
    {
    static double sum;
		
    if (p != NULL)
        {
        if (p->anc == NULL)
            sum = 0;
        else
            sum += p->branchLength;
        SumBranches (p->left);
        SumBranches (p->right);
        }
		
    return sum;
    }



/********************* CompareGenotypes  ************************/
/*
 Compares two unphased genotypes (i.e, A/T = T/A)
 This is very strict (calling A/A would incorrect if there is a deletion A/_)
 */
int CompareGenotypes (int a1, int a2, int b1, int b2)
	{
	int temp, equal;
	
	temp = 0;
	equal = YES;
	
	/* order the alleles to facilitate the comparisons */
	if (a1 > a2)
		{
		temp = a1;
		a1 = a2;
		a2 = temp;
		}
	
	if (b1 > b2)
		{
		temp = b1;
		b1 = b2;
		b2 = temp;
		}

	if (a1 != b1 || a2 != b2)
		equal = NO;
	
	return equal;
	}



/********************* AllocateCellStructure  ************************/
/*
 This cell structure keeps information basically about read count functions. It was
 build initially for making doublets.
 */
void AllocateCellStructure()
	{
	int i, j, k;
	
	cell = (CellStr*) calloc(numCells+1, sizeof(CellStr));
	if (!cell)
		{
		fprintf (stderr, "Could not allocate cell (%ld)\n", (numCells+1) * sizeof(CellStr));
		exit (-1);
		}

	for (i=0; i<numCells+1; i++)
		{
		cell[i].site  = (CellSiteStr*) calloc (numSites, sizeof(CellSiteStr));
		if (!cell[i].site)
			{
			fprintf (stderr, "Could not allocate the cell[%d].site structure\n", i);
			exit (-1);
			}
		for (j=0; j<numSites; j++)
			{
			cell[i].site[j].readCount = (int*) calloc (4, sizeof(int));
			if (!cell[i].site[j].readCount)
				{
				fprintf (stderr, "Could not allocate the cell[%d].site[%d].readCount structure\n", i, j);
				exit (-1);
				}
			
			cell[i].site[j].readCountDoublet = (int*) calloc (4, sizeof(int));
			if (!cell[i].site[j].readCountDoublet)
				{
				fprintf (stderr, "Could not allocate the cell[%d].site[%d].readCountDoublet structure\n", i, j);
				exit (-1);
				}
			
			cell[i].site[j].genLike = (double**) calloc (4, sizeof(double*));
			if (!cell[i].site[j].genLike)
				{
				fprintf (stderr, "Could not allocate the cell[%d].site[%d].genLike structure\n", i, j);
				exit (-1);
				}
			for (k=0; k<4; k++)
				{
				cell[i].site[j].genLike[k] = (double*) calloc (4, sizeof(double));
				if (! cell[i].site[j].genLike[k] )
					{
					fprintf (stderr, "Could not allocate the cell[%d].site[%d].genLike[%d] structure\n", i,j,k);
					exit (-1);
					}
				}
			
			cell[i].site[j].scaledGenLike = (double**) calloc (4, sizeof(double*));
			if (!cell[i].site[j].scaledGenLike)
				{
				fprintf (stderr, "Could not allocate the cell[%d].site[%d].scaledGenLike structure\n", i, j);
				exit (-1);
				}
			for (k=0; k<4; k++)
				{
				cell[i].site[j].scaledGenLike[k] = (double*) calloc (4, sizeof(double));
				if (!cell[i].site[j].scaledGenLike[k] )
					{
					fprintf (stderr, "Could not allocate the cell[%d].site[%d].scaledGenLike[%d] structure\n", i,j,k);
					exit (-1);
					}
				}
			
			cell[i].site[j].genLikeDoublet = (double**) calloc (4, sizeof(double*));
			if (!cell[i].site[j].genLikeDoublet)
				{
				fprintf (stderr, "Could not allocate the cell[%d].site[%d].genLikeDoublet structure\n", i, j);
				exit (-1);
				}
			for (k=0; k<4; k++)
				{
				cell[i].site[j].genLikeDoublet[k] = (double*) calloc (4, sizeof(double));
				if (!cell[i].site[j].genLikeDoublet[k] )
					{
					fprintf (stderr, "Could not allocate the cell[%d].site[%d].genLikeDoublet[%d] structure\n", i,j,k);
					exit (-1);
					}
				}
			
			cell[i].site[j].scaledGenLikeDoublet = (double**) calloc (4, sizeof(double*));
			if (!cell[i].site[j].scaledGenLikeDoublet)
				{
				fprintf (stderr, "Could not allocate the cell[%d].site[%d].scaledGenLikeDoublet structure\n", i, j);
				exit (-1);
				}
			for (k=0; k<4; k++)
				{
				cell[i].site[j].scaledGenLikeDoublet[k] = (double*) calloc (4, sizeof(double));
				if (!cell[i].site[j].scaledGenLikeDoublet[k] )
					{
					fprintf (stderr, "Could not allocate the cell[%d].site[%d].scaledGenLikeDoublet[%d] structure\n", i,j,k);
					exit (-1);
					}
				}
			}
		}
	}


/********************* GenerateReadCounts  ************************/
/*	For each individual SNV genotype the program will generate read counts
	given some sequencing depth or coverage. The number of reads follow a Poisson
	or a Negative Binomial distribution around mean coverage.
	Reads are randomly assigned to maternal/paternal chromosome, according to an allelic haploid coverage reduction (by default 0.5).
	Reads can contain some errors according to sequencing and amplification error parameters
*/

void GenerateReadCounts (long int *seed)
	{
	int			i, j, snv, equalGenotypes;
	double		*probs, **ngsEij, **ampEijmat, **ampEijpat;
	CellSiteStr *c;
	
	probs = (double*) calloc (4, sizeof(double));
	if (!probs)
		{
		fprintf (stderr, "Could not allocate the probs structure\n");
		exit (-1);
		}
	
	/* allocate error probabilities for NGS */
	ngsEij = (double**) calloc (4, sizeof(double**));
	if (!ngsEij)
		{
		fprintf (stderr, "Could not allocate the ngsEij structure\n");
		exit (-1);
		}
	for (i=0; i<4; i++)
		{
		ngsEij[i] = (double*) calloc (4, sizeof(double*));
		if (!ngsEij[i])
			{
			fprintf (stderr, "Could not allocate the ngsEij[i] structure\n");
			exit (-1);
			}
		}

	/* allocate error probabilities for amplification */
	ampEijmat = (double**) calloc (4, sizeof(double**));
	if (!ampEijmat)
		{
		fprintf (stderr, "Could not allocate the ampEijmat structure\n");
		exit (-1);
		}
	for (i=0; i<4; i++)
		{
		ampEijmat[i] = (double*) calloc (4, sizeof(double*));
		if (!ampEijmat[i])
			{
			fprintf (stderr, "Could not allocate the ampEijmat[i] structure\n");
			exit (-1);
			}
		}

	ampEijpat = (double**) calloc (4, sizeof(double**));
	if (!ampEijpat)
		{
		fprintf (stderr, "Could not allocate the ampEijpat structure\n");
		exit (-1);
		}
	for (i=0; i<4; i++)
		{
		ampEijpat[i] = (double*) calloc (4, sizeof(double*));
		if (!ampEijpat[i])
			{
			fprintf (stderr, "Could not allocate the ampEijpat[i] structure\n");
			exit (-1);
			}
		}

	/* initialize ngsEij now as this will be contant across sites */
	for (i=0; i<4; i++)
		for (j=0; j<4; j++)
			{
			if (i == j)
				ngsEij[i][j] = 1.0 - sequencingError;
			else
				ngsEij[i][j] = Eij[i][j]/cumEij[i][3] * sequencingError;
			}
	
	for (snv=0; snv<numSNVs; snv++) // produce read counts only for observed SNV sites
		{
		j = SNVsites[snv];
		for (i=0; i<numCells+1; i++)
			{
			SiteReadCounts (&cell[i].site[j], i, j, probs, ngsEij, ampEijmat, ampEijpat, seed);
			GenotypeLikelihoods (&cell[i].site[j], i, j, probs, ngsEij, ampEijmat, ampEijpat, seed);
			}
		}
	
	if (doubletRate > 0)
		MakeDoublets(probs, ngsEij, ampEijmat, ampEijpat, seed);

	/* count ML genotyping errors */
	countMLgenotypeErrors = 0;
	for (i=0; i<numCells+1; i++)
		for (snv=0; snv<numSNVs; snv++)
		{
		j = SNVsites[snv];
		c = &cell[i].site[j];
		if (cell[i].hasDoublet == NO)
			equalGenotypes = CompareGenotypes (c->trueMaternalAllele, c->truePaternalAllele, c->MLmatAllele, c->MLpatAllele);
		else
			equalGenotypes = CompareGenotypes (c->trueMaternalAllele, c->truePaternalAllele, c->MLmatAlleleDoublet, c->MLpatAlleleDoublet);

		if (equalGenotypes == NO)
			{
			countMLgenotypeErrors++;
			/*
			if (cell[i].hasDoublet == YES)
				fprintf(stderr,"\nD c%03d s%03d --  %c%c != R[%c%c] T[%c%c] ***", i+1, j+1, WhichNuc(c->MLmatAlleleDoublet), WhichNuc(c->MLpatAlleleDoublet), WhichNuc(c->maternalAllele), WhichNuc(c->paternalAllele), WhichNuc(c->trueMaternalAllele), WhichNuc(c->truePaternalAllele));
			else
				fprintf(stderr,"\nN c%03d s%03d --  %c%c != R[%c%c] T[%c%c]***", i+1, j+1, WhichNuc(c->MLmatAllele), WhichNuc(c->MLpatAllele), WhichNuc(c->maternalAllele), WhichNuc(c->paternalAllele), WhichNuc(c->trueMaternalAllele), WhichNuc(c->truePaternalAllele));
			}
		else
			{
			if (cell[i].hasDoublet == YES)
				fprintf(stderr,"\nD c%03d s%03d --  %c%c == R[%c%c] T[%c%c]", i+1, j+1, WhichNuc(c->MLmatAlleleDoublet), WhichNuc(c->MLpatAlleleDoublet), WhichNuc(c->maternalAllele), WhichNuc(c->paternalAllele), WhichNuc(c->trueMaternalAllele), WhichNuc(c->truePaternalAllele));
			else
				fprintf(stderr,"\nN c%03d s%03d --  %c%c == R[%c%c] T[%c%c]", i+1, j+1, WhichNuc(c->MLmatAllele), WhichNuc(c->MLpatAllele), WhichNuc(c->maternalAllele), WhichNuc(c->paternalAllele), WhichNuc(c->trueMaternalAllele), WhichNuc(c->truePaternalAllele));
			 */
			}
		}
			
	cumCountMLgenotypeErrors += countMLgenotypeErrors;

	free (probs);
	free (ngsEij);
	free (ampEijmat);
	free (ampEijpat);
	}


/********************* SiteReadCounts  ************************/
/*
 Generate read counts for a given genotype (cell and site)
 according to three different models considering sequencing and amplification error
 */
void SiteReadCounts (CellSiteStr *c, int i, int j, double *probs, double **ngsEij, double **ampEijmat, double **ampEijpat, long int *seed)
	{
	int		k, l, m;
	int		read, template;
	int		numMaternalReads, numPaternalReads;
	int		numReads;
	int		*readCount;
	int		maternalAllele, paternalAllele;
	int		thereIsMaternalAllele, thereIsPaternalAllele;
	double  maternalSiteAmplificationError, paternalSiteAmplificationError;
	int		ampErrorMaternalAllele, ampErrorPaternalAllele, goodTemplate, badTemplate;
	int		readsLeft, badReads;
	double	cumProb;
	double	uniform, goodTemplateProb, badTemplateProb, noSeqError;

	readCount = c->readCount;
	
	thereIsMaternalAllele = YES;
	thereIsPaternalAllele = YES;
	
	maternalAllele = data[MATERNAL][i][j];
	if (maternalAllele == ADO || maternalAllele == DELETION)
		thereIsMaternalAllele = NO;
	
	paternalAllele = data[PATERNAL][i][j];
	if (paternalAllele == ADO || paternalAllele == DELETION)
		thereIsPaternalAllele = NO;
	
	/* initialize some variables */
	noSeqError = 0;
	badTemplateProb = 0;
	goodTemplateProb = 0;
	ampErrorMaternalAllele = -9;
	ampErrorPaternalAllele = -9;
	badTemplate = -9;
	goodTemplate = -9;
	maternalSiteAmplificationError = 0;
	paternalSiteAmplificationError = 0;
	
	/* number of reads will follow a Poisson or a Negative Binomial distribution (overdispersed coverage) around mean coverage, and are randomly assigned to maternal/paternal chromosome */
	numMaternalReads = 0;
	numPaternalReads = 0;
	
	if (rateVarCoverage == YES) /* Negative Binomial */
		{
		if (thereIsMaternalAllele == YES && thereIsPaternalAllele == YES) /* both alleles are present */
			{
			numMaternalReads = RandomNegativeBinomial(allelicImbalance*coverage, alphaCoverage, seed);
			numPaternalReads = RandomNegativeBinomial((1.0-allelicImbalance)*coverage, alphaCoverage, seed);
			}
		else if (thereIsMaternalAllele == YES && thereIsPaternalAllele == NO)  /* only maternal allele */
			{
			numMaternalReads = RandomNegativeBinomial(haploidCoverageReduction*coverage, alphaCoverage, seed);
			numPaternalReads = 0;
			}
		else if (thereIsMaternalAllele == NO && thereIsPaternalAllele == YES)  /* only paternal allele */
			{
			numMaternalReads = 0;
			numPaternalReads = RandomNegativeBinomial(haploidCoverageReduction*coverage, alphaCoverage, seed);
			}
		else /* no allele is present (locus dropout/deletion) */
			{
			numMaternalReads = 0;
			numPaternalReads = 0;
			}
		}
	else /* Poisson distributed */
		{
		if (thereIsMaternalAllele == YES && thereIsPaternalAllele == YES) /* both alleles are present */
			{
			numMaternalReads = RandomPoisson(allelicImbalance*coverage, seed);
			numPaternalReads = RandomPoisson((1.0-allelicImbalance)*coverage, seed);
			}
		else if (thereIsMaternalAllele == YES && thereIsPaternalAllele == NO)  /* only maternal allele */
			{
			numMaternalReads = RandomPoisson(haploidCoverageReduction*coverage, seed);
			numPaternalReads = 0;
			}
		else if (thereIsMaternalAllele == NO && thereIsPaternalAllele == YES)  /* only paternal allele */
			{
			numMaternalReads = 0;
			numPaternalReads = RandomPoisson(haploidCoverageReduction*coverage, seed);
			}
		else /* no allele is present (locus dropout/deletion) */
			{
			numMaternalReads = 0;
			numPaternalReads = 0;
			}
		}
	
	numReads = numMaternalReads + numPaternalReads;
	
	/* read count simulation */
	readCount[A] = readCount[C] = readCount[G] = readCount[T] = 0;

	if (numReads > 0)
		{
		/* choose amplification error for this site. In this case we will have a proportion of alternative templates depending on whether the error occurred early or late in the amplification (we use a Beta distribution to modelize this proportion) */
		if (meanAmplificationError > 0)
			{
			maternalSiteAmplificationError = RandomBeta(meanAmplificationError, varAmplificationError, seed);
			paternalSiteAmplificationError = RandomBeta(meanAmplificationError, varAmplificationError, seed);
			}
		
		/* initialize ampEijmat */
		for (k=0; k<4; k++)
			for (l=0; l<4; l++)
				{
				if (k == l)
					ampEijmat[k][l] = 1.0 - maternalSiteAmplificationError;
				else
					ampEijmat[k][l] = Eij[k][l]/cumEij[k][3] * maternalSiteAmplificationError;
				}
		
		/* initialize ampEijpat */
		for (k=0; k<4; k++)
			for (l=0; l<4; l++)
				{
				if (k == l)
					ampEijpat[k][l] = 1.0 - paternalSiteAmplificationError;
				else
					ampEijpat[k][l] = Eij[k][l]/cumEij[k][3] * paternalSiteAmplificationError;
				}
		
		if (simulateOnlyTwoTemplates == YES)
			{
			/* choose base for wrong maternal template */
			ampErrorMaternalAllele = maternalAllele;
			if (maternalSiteAmplificationError > 0 && thereIsMaternalAllele == YES )
				{
				uniform = RandomUniform(seed) * cumEij[maternalAllele][3];
				for (m=0; m<4; m++)
					if (uniform <= cumEij[maternalAllele][m])
						{
						ampErrorMaternalAllele = m;
						break;
						}
				}
			/* choose base for wrong paternal template */
			ampErrorPaternalAllele = paternalAllele;
			if (paternalSiteAmplificationError > 0 && thereIsPaternalAllele == YES )
				{
				uniform = RandomUniform(seed) * cumEij[paternalAllele][3];
				for (m=0; m<4; m++)
					if (uniform <= cumEij[paternalAllele][m])
						{
						ampErrorPaternalAllele = m;
						break;
						}
				}
			}
		
		/* probabilities for the different types of reads  */
		if (sequencingError > 0 || maternalSiteAmplificationError > 0 || paternalSiteAmplificationError > 0)
			{
			if (thereIsMaternalAllele == YES)
				{
				if (simulateOnlyTwoTemplates == YES)
					{
					/* probabilities for maternal reads  (assuming a chosen template )*/
					noSeqError = 1.0 - sequencingError;
					badTemplateProb = maternalSiteAmplificationError;
					goodTemplateProb = 1.0 - maternalSiteAmplificationError;
					badTemplate = ampErrorMaternalAllele;
					goodTemplate = maternalAllele;
					for (read=0; read<4; read++)
						{
						if (read == maternalAllele)
							probs[read] = goodTemplateProb * noSeqError + badTemplateProb * ngsEij[badTemplate][read];
						else if (read == ampErrorMaternalAllele)
							probs[read] = goodTemplateProb * ngsEij[goodTemplate][read] + badTemplateProb * noSeqError;
						else
							probs[read] = goodTemplateProb * ngsEij[goodTemplate][read] + badTemplateProb * ngsEij[badTemplate][read];
						}
					}
				else
					{
					/* probabilities for maternal reads integrating over all potential templates) */
					for (read=0; read<4; read++)
						{
						probs[read] = 0;
						for (template=0; template<4; template++)
							probs[read] += ampEijmat[maternalAllele][template] * ngsEij[template][read];
						}
					}
				
				cumProb = 0;
				readsLeft = numMaternalReads;
				/* multinomial assignment of reads */
				for (read=0; read<4; read++)
					{
					if (read != maternalAllele)
						{
						badReads = RandomBinomial (probs[read] / (1 - cumProb), readsLeft, seed);
						readsLeft -= badReads;
						readCount[read] += badReads;
						cumProb += probs[read];
						}
					}
				readCount[maternalAllele] += readsLeft; /* good reads */
				}
			
			if (thereIsPaternalAllele == YES)
				{
				if (simulateOnlyTwoTemplates == YES)
					{
					/* probabilities for paternal reads  (assuming a chosen template )*/
					badTemplateProb = paternalSiteAmplificationError;
					goodTemplateProb = 1.0 - paternalSiteAmplificationError;
					badTemplate = ampErrorPaternalAllele;
					goodTemplate = paternalAllele;
					
					for (read=0; read<4; read++)
						{
						if (read == paternalAllele)
							probs[read] = goodTemplateProb * noSeqError + badTemplateProb * ngsEij[badTemplate][read];
						else if (read == ampErrorPaternalAllele)
							probs[read] = goodTemplateProb * ngsEij[goodTemplate][read] + badTemplateProb * noSeqError;
						else
							probs[read] = goodTemplateProb * ngsEij[goodTemplate][read] + badTemplateProb * ngsEij[badTemplate][read];
						}
					}
				else
					{
					/* probabilities for paternal reads  (integrating over all potential templates) */
					for (read=0; read<4; read++)
						{
						probs[read] = 0;
						for (template=0; template<4; template++)
							probs[read] += ampEijpat[paternalAllele][template] * ngsEij[template][read];
						}
					}
				
				cumProb = 0;
				readsLeft = numPaternalReads;
				/* multinomial assignment of reads */
				for (read=0; read<4; read++)
					{
					if (read != paternalAllele)
						{
						badReads = RandomBinomial (probs[read] / (1 - cumProb), readsLeft, seed);
						readsLeft -= badReads;
						readCount[read] += badReads;
						cumProb += probs[read];
						}
					}
				readCount[paternalAllele] += readsLeft; // good reads
				}
			}
		else /* no sequencing or amplification error */
			{
			readCount[maternalAllele] += numMaternalReads;
			readCount[paternalAllele] += numPaternalReads;
			}
		} // numReads > 0
	
	c->numReads = numReads;
	c->readCount = readCount;
	c->maternalAllele = maternalAllele;
	c->paternalAllele = paternalAllele;
	c->thereIsMaternalAllele = thereIsMaternalAllele;
	c->thereIsPaternalAllele = thereIsPaternalAllele;
	c->maternalSiteAmplificationError = maternalSiteAmplificationError;
	c->paternalSiteAmplificationError = paternalSiteAmplificationError;
	}


/********************* GenotypeLikelihoods  ************************/
/*
 Calculate genotype likelihoods for a site and cell given the simulated reads
 according to three different models
 
 **	Genotype likelihoods can be calculated given the read counts and the sequencing error
 (see Korneliussen 2013 BMCBioinf)

 00 01 02 03 11 12 13 22 23 33
 AA AC AG AT CC CG CT GG GT TT
 */

void GenotypeLikelihoods (CellSiteStr *c, int i, int j, double *probs, double **ngsEij, double **ampEijmat, double **ampEijpat, long int *seed)
{
	int		k, l, a1, a2;
	int		read, numReads;
	int		*readCount;
	int		maternalAllele, paternalAllele;
	int		thereIsMaternalAllele, thereIsPaternalAllele;
	double  maternalSiteAmplificationError, paternalSiteAmplificationError;
	int		template, template1, template2;
	int		debug_GL, fix_parameters;
	double	maxLike, pReadGivenA1, pReadGivenA2;
	double	homozygote_for_read, heterozygote_for_read, no_allele_for_read;
	double	**genLike, **scaledGenLike;
	double 	MLmatAllele, MLpatAllele;

	readCount = c->readCount;
	numReads = c->numReads;
	genLike = c->genLike;
	scaledGenLike = c->scaledGenLike;
	maternalAllele = c->maternalAllele;
	paternalAllele = c->paternalAllele;
	thereIsMaternalAllele = c->thereIsMaternalAllele;
	thereIsPaternalAllele = c->thereIsPaternalAllele;
	maternalSiteAmplificationError = c->maternalSiteAmplificationError;
	paternalSiteAmplificationError = c->paternalSiteAmplificationError;

	debug_GL = NO;
	fix_parameters = NO;
	if (fix_parameters == YES)
		{
		maternalAllele = c->maternalAllele = A;
		paternalAllele = c->paternalAllele = C;
		thereIsMaternalAllele = c->thereIsMaternalAllele = YES;
		thereIsPaternalAllele = c->thereIsPaternalAllele = YES;
		readCount[A] = 50;
		readCount[C] = 40;
		readCount[G] = 20;
		readCount[T] = 10;
		numReads = c->numReads = 120;
		sequencingError = 0.1;
		maternalSiteAmplificationError = c->maternalSiteAmplificationError = 0.000;
		paternalSiteAmplificationError = c->paternalSiteAmplificationError = 0.000;

		/* reinitialize ngsEij  */
		for (k=0; k<4; k++)
			for (l=0; l<4; l++)
				{
				if (k == l)
					ngsEij[k][l] = 1.0 - sequencingError;
				else
					ngsEij[k][l] = Eij[k][l]/cumEij[k][3] * sequencingError;
				}
		
		/* reinitialize ampEijmat */
		for (k=0; k<4; k++)
			for (l=0; l<4; l++)
				{
				if (k == l)
					ampEijmat[k][l] = 1.0 - maternalSiteAmplificationError;
				else
					ampEijmat[k][l] = Eij[k][l]/cumEij[k][3] * maternalSiteAmplificationError;
				}
		
		/* reinitialize ampEijpat */
		for (k=0; k<4; k++)
			for (l=0; l<4; l++)
				{
				if (k == l)
					ampEijpat[k][l] = 1.0 - paternalSiteAmplificationError;
				else
					ampEijpat[k][l] = Eij[k][l]/cumEij[k][3] * paternalSiteAmplificationError;
				}
		}

	
	/************************** MODEL 0 (GATK) *************************************/
	/* MODEL 0: genotype log10 likelihoods according to the observed read assuming a single sequencing error independent of the nucleotides involved, and no amplification error */
	/* extended from Korneliussen 2013 BMCBioinf
	 p(D|G=A1,A2) = multiply over reads [ 1/2 p(read i|A1) + 1/2 p(read i|A2)]
	 where p(read|A) e/3 if b!=A  or 1-e if b=A
	 N.B: ANGSD seems to use natural logs!
	 */
	if (maternalSiteAmplificationError == 0 && paternalSiteAmplificationError == 0 )
		{
		homozygote_for_read = log10(1.0 - sequencingError);
		no_allele_for_read = log10(sequencingError/3.0);
		heterozygote_for_read = log10((1.0 - sequencingError)/2.0 + sequencingError/6.0);
		
		for (a1=0; a1<4; a1++) //a1 is maternal
			for (a2=a1; a2<4; a2++) //a2 is paternal
				{
				genLike[a1][a2] = -0.0;
				
				for (k=0; k<4; k++)
					{
					if (readCount[k] > 0)
						{
						if (k == a1 && k == a2)
							genLike[a1][a2] += readCount[k] * homozygote_for_read;
						else if (k != a1 && k != a2)
							genLike[a1][a2] += readCount[k] * no_allele_for_read;
						else
							genLike[a1][a2] += readCount[k] * heterozygote_for_read;
							//fprintf (stderr, "\nread=%c gl[%c][%c] = %lf", WhichNuc(k), WhichNuc(a1), WhichNuc(a2), genLike[a1][a2]);
						}
					}
				}
		
		if (debug_GL == YES)
			{
			fprintf (stderr,"\ncell %d site %d: %c%c", i+1, j+1, WhichNuc(maternalAllele), WhichNuc(paternalAllele));
			fprintf (stderr,"|\treads: A:%d C:%d G:%d T:%d  | total:%d ",  readCount[A], readCount[C], readCount[G], readCount[T], numReads);
			fprintf (stderr,"\nmatAmpError = %f", maternalSiteAmplificationError);
			fprintf (stderr,"\npatAmpError = %f", paternalSiteAmplificationError);
			for (a1=0; a1<4; a1++)
				for (a2=a1; a2<4; a2++)
					fprintf (stderr, "\n(1)log10 GL[%c][%c] = %lf", WhichNuc(a1), WhichNuc(a2), genLike[a1][a2]);
			fprintf (stderr, "\n");
			}
		} // model 0
	

	
	/************************** MODEL 1 (4-templates) *************************************/
	/* MODEL 1 (4-templates): genotype log10 likelihoods according to the observed read, assuming a single sequencing error independent of the nucleotides involved, and amplification error with all 4 templates */
	if (simulateOnlyTwoTemplates == NO && (maternalSiteAmplificationError > 0 || paternalSiteAmplificationError > 0))
		{
		for (a1=0; a1<4; a1++) //a1 is maternal
			for (a2=a1; a2<4; a2++) //a2 is paternal
				{
				genLike[a1][a2] = -0.0;
				for (read=0; read<4; read++)
					{
					pReadGivenA1 = pReadGivenA2 = 0.0;
					if (readCount[read] > 0)
						{
						for (template=0; template<4; template++)
							{
							pReadGivenA1 += ampEijmat[a1][template] * ngsEij[template][read];
							pReadGivenA2 += ampEijpat[a2][template] * ngsEij[template][read];
								//fprintf (stderr, "\n ampEijmat[%c]][%c]=%lf ngsEij[%c]][%c]=%lf ", WhichNuc(a1), WhichNuc(template), ampEijmat[a1][template], WhichNuc(a2), WhichNuc(template), ngsEij[template][read]);
							}
						genLike[a1][a2] += readCount[read] * log10(pReadGivenA1*0.5 + pReadGivenA2*0.5);
						}
						//fprintf (stderr, "\n***read=%c gl[%c][%c] = %lf  (p1=%lf p2=%lf)", WhichNuc(read), WhichNuc(a1), WhichNuc(a2), genLike[a1][a2],pReadGivenA1,pReadGivenA2);
					}
				}
		
		if (debug_GL == YES)
			{
			fprintf (stderr,"\ncell %d site %d: %c%c", i+1, j+1, WhichNuc(maternalAllele), WhichNuc(paternalAllele));
			fprintf (stderr,"|\treads: A:%d C:%d G:%d T:%d  | total:%d",  readCount[A], readCount[C], readCount[G], readCount[T], numReads);
			fprintf (stderr,"\nmatAmpError = %f", maternalSiteAmplificationError);
			fprintf (stderr,"\npatAmpError = %f", paternalSiteAmplificationError);
			for (a1=0; a1<4; a1++)
				for (a2=a1; a2<4; a2++)
					fprintf (stderr, "\n(2)log10 GL[%c][%c] = %lf", WhichNuc(a1), WhichNuc(a2), genLike[a1][a2]);
			fprintf (stderr, "\n");
			}
		} //model 1

	
	/************************** MODEL 2 (2-templates) *************************************/
	/* MODEL 2 (2-templates): genotype log10 likelihoods according to the observed read, assuming a single sequencing error independent of the nucleotides involved, and amplification error with up to 2 templates */
	/* this calculation seems OK but I cannot say I am 100% sure */
	if (simulateOnlyTwoTemplates == YES && (maternalSiteAmplificationError > 0 || paternalSiteAmplificationError > 0))
		{
		for (a1=0; a1<4; a1++) //a1 is maternal
			for (a2=a1; a2<4; a2++) //a2 is paternal
				{
				genLike[a1][a2] = -0.0;
				for (template1=0; template1<4; template1++)
					if (template1 != a1)
						{
						for (template2=0; template2<4; template2++)
							if (template2 != a2)
								{
								for (read=0; read<4; read++)
									{
									if (readCount[read] > 0)
										{
										pReadGivenA1 = (1.0 - maternalSiteAmplificationError) * ngsEij[a1][read] + ampEijmat[a1][template1] * ngsEij[template1][read];
										pReadGivenA2 = (1.0 - paternalSiteAmplificationError) * ngsEij[a2][read] + ampEijpat[a2][template2] * ngsEij[template2][read];
										if (pReadGivenA1 > 0 || pReadGivenA2 >0)
											genLike[a1][a2] += Eij[a1][template1]/cumEij[a1][3] * Eij[a2][template2]/cumEij[a2][3] * readCount[read] * log10(pReadGivenA1*0.5 + pReadGivenA2*0.5);
											//fprintf (stderr, "\n***read=%c gl[%c][%c] = %lf  (p1=%lf p2=%lf)  ampEijmat[%c][%c]=%lf  ampEijmat[%c]][%c]=%lf ngsEij[%c]][%c]=%lf ngsEij[%c]][%c]=%lf", WhichNuc(read), WhichNuc(a1), WhichNuc(a2), genLike[a1][a2],pReadGivenA1,pReadGivenA2,
											//WhichNuc(a1), WhichNuc(template1), ampEijmat[a1][template1], WhichNuc(a2), WhichNuc(template2), ampEijmat[a2][template2], WhichNuc(a1), WhichNuc(template1), ngsEij[template1][read], WhichNuc(a2), WhichNuc(template2), ngsEij[template2][read]);
										}
									}
								}
						}
				}
		
		if (debug_GL == YES)
			{
			fprintf (stderr,"\ncell %d site %d: %c%c", i+1, j+1, WhichNuc(maternalAllele), WhichNuc(paternalAllele));
			fprintf (stderr,"|\treads: A:%d C:%d G:%d T:%d  | total:%d",  readCount[A], readCount[C], readCount[G], readCount[T], numReads);
			fprintf (stderr,"\nmatAmpError = %f", maternalSiteAmplificationError);
			fprintf (stderr,"\npatAmpError = %f", paternalSiteAmplificationError);
			for (a1=0; a1<4; a1++)
				for (a2=a1; a2<4; a2++)
					fprintf (stderr, "\n(3)log10 GL[%c][%c] = %lf", WhichNuc(a1), WhichNuc(a2), genLike[a1][a2]);
			fprintf (stderr, "\n");
			}
		} // model 2

	/* find max log10 likelihood  */
	maxLike = genLike[0][0];
	MLmatAllele = MLpatAllele = 0;
	for (a1=0; a1<4; a1++)
		for (a2=a1; a2<4; a2++)
			if (genLike[a1][a2] > maxLike)
				{
				maxLike = genLike[a1][a2];
				MLmatAllele = a1;
				MLpatAllele = a2;
				}

	/* rescale to likelihood ratios */
	for (a1=0; a1<4; a1++)
		for (a2=a1; a2<4; a2++)
			scaledGenLike[a1][a2] = genLike[a1][a2] - maxLike;
	
	// Qphred = -10 log(10) Perror)
	// Perror = 10 ^(-Qphred/10)
	
	if (debug_GL == YES)
		{
		//PrintSiteInfo (stderr, SNVsites[snv]);
		for (a1=0; a1<4; a1++)
			for (a2=a1; a2<4; a2++)
				fprintf (stderr, "\nscalog10 GL[%c][%c] = %lf", WhichNuc(a1), WhichNuc(a2), scaledGenLike[a1][a2]);
		fprintf (stderr, "\n");
		}

	c->MLmatAllele = MLmatAllele;
	c->MLpatAllele = MLpatAllele;

}

/********************* MakeDoublets  ************************/
/*
 We simulate the possibility of amplifying/sequencing two cells at once
 by adding read counts from another sampled cell.
 Only tumor cells are considered.
*/

void MakeDoublets (double *probs, double **ngsEij, double **ampEijmat, double **ampEijpat, long int *seed)
	{
	int		c1, c2;
	int 	j, k, l, snv;
	int		debug_DL;
	int		MLmatAllele, MLpatAllele;
	double	maxLike;
	CellStr celltype;

	debug_DL = NO;
	if (debug_DL == YES)
		fprintf(stderr,"\n\n\n********Doublet read counts and likelihoods");

	/* choose whether a tumor cell will become a doublet */
	for (c1=0; c1<numCells;c1++)
		{
		if (RandomUniform(seed) < doubletRate)
			{
			cell[c1].hasDoublet = YES;
			
			/* select a cell genotype from the tumor sample */
			c2 = RandomUniformTo (numCells,seed);

			/* allocate temporal memory for this cell genotype */
			celltype.site  = (CellSiteStr*) calloc (numSites, sizeof(CellSiteStr));
			if (!celltype.site)
				{
				fprintf (stderr, "Could not allocate the celltype.site structure\n");
				exit (-1);
				}

			for (j=0; j<numSites; j++)
				{
				celltype.site[j].readCount = (int*) calloc (4, sizeof(int));
				if (!celltype.site[j].readCount)
					{
					fprintf (stderr, "Could not allocate the celltype.site[%d].readCount structure\n", j);
					exit (-1);
					}
			
				celltype.site[j].genLike = (double**) calloc (4, sizeof(double*));
				if (!celltype.site[j].genLike)
					{
					fprintf (stderr, "Could not allocate the celltype.site[%d].genLike structure\n", j);
					exit (-1);
					}
				for (k=0; k<4; k++)
					{
					celltype.site[j].genLike[k] = (double*) calloc (4, sizeof(double));
					if (!celltype.site[j].genLike[k] )
						{
						fprintf (stderr, "Could not allocate the celltype.site[%d].genLike[%d] structure\n", j,k);
						exit (-1);
						}
					}

				celltype.site[j].scaledGenLike = (double**) calloc (4, sizeof(double*));
				if (!celltype.site[j].scaledGenLike)
					{
					fprintf (stderr, "Could not allocate the celltype.site[%d].scaledGenLike structure\n", j);
					exit (-1);
					}
				for (k=0; k<4; k++)
					{
					celltype.site[j].scaledGenLike[k] = (double*) calloc (4, sizeof(double));
					if (!celltype.site[j].scaledGenLike[k] )
						{
						fprintf (stderr, "Could not allocate the celltype.site[%d].scaledGenLike[%d] structure\n", j,k);
						exit (-1);
						}
					}
				}
			
			/* get  read counts and likelihoods for this cell genotype for all sites */
			for (snv=0; snv<numSNVs; snv++) // produce read counts only for observed SNV sites
				{
				if (debug_DL == YES)
						fprintf (stderr,"\n\n-- New celltype reads and likelihoods");
				j = SNVsites[snv];
				SiteReadCounts (&celltype.site[j], c2, j, probs, ngsEij, ampEijmat, ampEijpat, seed);
				GenotypeLikelihoods (&celltype.site[j], c2, j, probs, ngsEij, ampEijmat, ampEijpat, seed);

				if (debug_DL == YES)
					{
					fprintf (stderr,"\ncelltype %d site %d: %c%c", c2+1, j+1, WhichNuc(celltype.site[j].maternalAllele), WhichNuc(celltype.site[j].paternalAllele));
					fprintf (stderr,"|\treads: A:%d C:%d G:%d T:%d  | total:%d ",  celltype.site[j].readCount[A], celltype.site[j].readCount[C], celltype.site[j].readCount[G], celltype.site[j].readCount[T], celltype.site[j].numReads);
					fprintf (stderr,"\nmatAmpError = %f", celltype.site[j].maternalSiteAmplificationError);
					fprintf (stderr,"\npatAmpError = %f", celltype.site[j].paternalSiteAmplificationError);
					for (k=0; k<4; k++)
						for (l=k; l<4; l++)
							fprintf (stderr, "\nlog10 GL[%c][%c] = %lf", WhichNuc(k), WhichNuc(l), celltype.site[j].genLike[k][l]);
					fprintf (stderr, "\n");
					}
				
				/* add doublet information */
				cell[c1].site[j].numReadsDoublet = cell[c1].site[j].numReads + celltype.site[j].numReads;
				for (k=0; k<4; k++)
					cell[c1].site[j].readCountDoublet[k] = cell[c1].site[j].readCount[k] + celltype.site[j].readCount[k];

				cell[c1].site[j].maternalAlleleDoublet = celltype.site[j].maternalAllele;
				cell[c1].site[j].paternalAlleleDoublet = celltype.site[j].paternalAllele;

				//fprintf(stderr,"\nln cell\t\tcelltype\tdoublet");
				/* adding the genotype log likelihoods of two sets of reads (c1 + celltype) to get the total doublet likelihoods */
				for (k=0; k<4; k++)
					for (l=k; l<4; l++)
						{
						if (isinf(cell[c1].site[j].genLike[k][l]) || isinf(celltype.site[j].genLike[k][l]))
							cell[c1].site[j].genLikeDoublet[k][l] = -1.0/0.0;
						else
							cell[c1].site[j].genLikeDoublet[k][l] = cell[c1].site[j].genLike[k][l] + celltype.site[j].genLike[k][l];
						}
				
				if (debug_DL == YES)
					{
					fprintf (stderr,"\nCell %d site %d: %c%c", c1+1, j+1, WhichNuc(cell[c1].site[j].maternalAllele), WhichNuc(cell[c1].site[j].paternalAllele));
					fprintf (stderr,"|\treads: A:%d C:%d G:%d T:%d  | total:%d ",   cell[c1].site[j].readCount[A],  cell[c1].site[j].readCount[C],  cell[c1].site[j].readCount[G],  cell[c1].site[j].readCount[T],  cell[c1].site[j].numReads);
					fprintf (stderr,"\n|reads doublet: A:%d C:%d G:%d T:%d  | total:%d ",   cell[c1].site[j].readCountDoublet[A],  cell[c1].site[j].readCountDoublet[C],  cell[c1].site[j].readCountDoublet[G],  cell[c1].site[j].readCountDoublet[T],  cell[c1].site[j].numReadsDoublet);
					fprintf (stderr,"\nmatAmpError = %f",  cell[c1].site[j].maternalSiteAmplificationError);
					fprintf (stderr,"\npatAmpError = %f",  cell[c1].site[j].paternalSiteAmplificationError);
					for (k=0; k<4; k++)
						for (l=k; l<4; l++)
							fprintf (stderr, "\ndoublet log10 GL[%c][%c] = %lf", WhichNuc(k), WhichNuc(l),  cell[c1].site[j].genLikeDoublet[k][l]);
					fprintf (stderr, "\n");
					}

				
				/* reescale doublet genotype likelihood to the maximum likelihood and assign ML alleles */
				/* find max log10 doublet likelihood  */
				maxLike = cell[c1].site[j].genLikeDoublet[0][0];
				MLmatAllele = MLpatAllele = 0;
				for (k=0; k<4; k++)
					for (l=k; l<4; l++)
						if (cell[c1].site[j].genLikeDoublet[k][l] > maxLike)
							{
							maxLike = cell[c1].site[j].genLikeDoublet[k][l];
							MLmatAllele = k;
							MLpatAllele = l;
							}

				/* rescale to likelihood ratios */
				for (k=0; k<4; k++)
					for (l=k; l<4; l++)
						cell[c1].site[j].scaledGenLikeDoublet[k][l] = cell[c1].site[j].genLikeDoublet[k][l] - maxLike;
						
				cell[c1].site[j].MLmatAlleleDoublet = MLmatAllele;
				cell[c1].site[j].MLpatAlleleDoublet = MLpatAllele;
			
				if (debug_DL == YES)
					{
					for (k=0; k<4; k++)
						for (l=k; l<4; l++)
							fprintf (stderr, "\nscaled doublet log10 GL[%c][%c] = %lf", WhichNuc(k), WhichNuc(l), cell[c1].site[j].scaledGenLikeDoublet[k][l]);
					fprintf (stderr, "\n");
					}
				}
			
			//free celltype
			for (j=0; j<numSites; j++)
				{
				free(celltype.site[j].readCount);
				for (k=0; k<4; k++)
					{
					free(celltype.site[j].genLike[k]);
					free(celltype.site[j].scaledGenLike[k]);
					}
				free(celltype.site[j].genLike);
				free(celltype.site[j].scaledGenLike);
				}
			} // if we do doublet for c1
		} //c1
	}



/********************* PrintCATG  ************************/
/**	CATG output format with read counts
	http://nbviewer.jupyter.org/gist/dereneaton/d2fd4e70d29f5ee5d195/testing_cat.ipynb#View-the-.cat-results-files

	The first line has the number of samples and the number of sites.
	Following this is a transposed data matrix with the consensus base calls where each row is a site and each column a different sample.
	To the right of each site is a tab-separated list of counts of the four bases at that site in the order C,A,T,G.

	12 44500
	1A0	1B0	1C0	1D0	2E0	2F0	2G0	2H0	3I0
	YCCCCCCCCCCC	10,0,10,0	20,0,0,0	20,0,0,0	20,0,0, ...
	GGGGGGGGGGGG	0,0,0,20	0,0,0,20	0,0,0,20	0,0,0,20 ...
	AAAAAAAAAAAA	0,20,0,0	0,20,0,0	0,20,0,0	0,20,0,0 ...
	CCCCCCCCCCCC	20,0,0,0	20,0,0,0	20,0,0,0	20,0,0,0 ...
	AAAAAAAAAAAA	0,20,0,0	0,20,0,0	0, 20,0,0	0,20,0,0 ...
	...
*/
void PrintCATG (FILE *fp)
	{
	int	i, j, snv;

	fprintf (fp,"%d %d\n",numCells+1, numSNVs);
	for (i=0; i<numCells; i++)
		fprintf (fp,"tumcell%04d  ", i+1);
	fprintf (fp,"healthycell  ");
		
	for (snv=0; snv<numSNVs; snv++)
		{
		j = SNVsites[snv];
		fprintf (fp,"\n");
		for (i=0; i<numCells+1; i++)
			fprintf (fp,"%c",WhichIUPAC(data[MATERNAL][i][j], data[PATERNAL][i][j]));
		for (i=0; i<numCells+1; i++)
			fprintf (fp,"\t%d,%d,%d,%d",  cell[i].site[j].readCount[C], cell[i].site[j].readCount[A], cell[i].site[j].readCount[T], cell[i].site[j].readCount[G]);
		}

	fprintf (fp,"\n");
	}


/********************* PrintVCF  ************************/
/**	VCF output format with read counts
*/
void PrintVCF (FILE *fp)
	{
	int		i, j, k, l, snv;
	int		maternalAllele, paternalAllele, referenceAllele;
	int		trueMaternalAllele, truePaternalAllele;
	int		maternalAlleleDoublet, paternalAlleleDoublet;

	/* print file header */
	fprintf (fp,"##fileformat=VCFv4.3");
	fprintf (fp,"\n##filedate=");
	PrintDate (fp);
	fprintf (fp,"##source=CellCoal SNV simulation");
	fprintf (fp,"\n##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral allele\">");
	fprintf (fp,"\n##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">");
	fprintf (fp,"\n##INFO=<ID=AF,Number=R,Type=Float,Description=\"True alternate/s allele frequency (ignores doublets)\">");
	fprintf (fp,"\n##FILTER=<ID=s50,Description=\"Less than half of samples have data\">");
	fprintf (fp,"\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"True genotype in the NGS library (considers ADO/DEL; ignores doublets)\">");
	fprintf (fp,"\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">");
	fprintf (fp,"\n##FORMAT=<ID=RC,Number=4,Type=Integer,Description=\"Read counts for AGCT\">");
	fprintf (fp,"\n##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Scaled genotype likelihoods in log10 (AA AC AG AT CC CG CT GG GT TT)\">");
	fprintf (fp,"\n##FORMAT=<ID=ML,Number=1,Type=String,Description=\"Maximum likelihood genotype\">");
	fprintf (fp,"\n##FORMAT=<ID=NG,Number=1,Type=String,Description=\"True SNV genotype in the NGS library (considers ADO/DEL; ignores doublets)\">");
	fprintf (fp,"\n##FORMAT=<ID=DG,Number=1,Type=Integer,Description=\"True SNV genotype in the NGS library for the second cell if there is a doublet\">");
	fprintf (fp,"\n##FORMAT=<ID=TG,Number=1,Type=String,Description=\"True SNV genotype in the original cell (considers DEL; ignores ADO/doublets\">");
	fprintf (fp,"\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	
	for (i=0; i<numCells; i++)
		fprintf (fp,"\ttumcell%04d", i+1);
	fprintf (fp,"\thealthycell");

	for (snv=0; snv<numSNVs; snv++)
		{
		j = SNVsites[snv];
		
		 /* VFC: CHROMOSOME */
		fprintf (fp,"\n%d", 1);
		
		/* VFC: POSITION */
		fprintf (fp,"\t%d", j+1);
		
		/* VFC: ID */
		if (allSites[j].isSNV == NO)
			fprintf (fp,"\tinv%06d", j+1);
		else
			fprintf (fp,"\tsnv%06d", j+1);

		/* VFC: REFERENCE allele(s) => healthy root alleles */
		referenceAllele = allSites[j].referenceAllele;
		fprintf (fp,"\t%c", WhichNuc(referenceAllele));
			
		/* VFC: ALTERNATE allele(s) */
		fprintf (fp,"\t");
		if (allSites[j].countACGT == allSites[j].countA ||
		allSites[j].countACGT == allSites[j].countC ||
		allSites[j].countACGT == allSites[j].countG ||
		allSites[j].countACGT == allSites[j].countT)
			fprintf (fp,".");
		else
			{
			if (allSites[j].countA > 0 && referenceAllele != A)
				fprintf (fp,"A,");
			if (allSites[j].countC > 0 && referenceAllele != C)
				fprintf (fp,"C,");
			if (allSites[j].countG > 0 && referenceAllele != G)
				fprintf (fp,"G,");
			if (allSites[j].countT > 0 && referenceAllele != T)
				fprintf (fp,"T,");
			fseek(fp, -1, SEEK_CUR); 	/* get rid of the last comma */
			}
		
		/* VFC: QUALITY */  /* Qphred probability that a SNV exists at this site; I will put missing info here */
		fprintf (fp, "\t.");

		/* VFC: FILTER */
		if (allSites[j].countACGT < allSites[j].countDropped)  /* if less than half of the samples have missing data */
			fprintf (fp, "\ts50");
		else
			fprintf (fp, "\tPASS");

		/* VFC: INFO: AA (ancestral allele) */
		fprintf (fp, "\tAA=%c", WhichNuc(referenceAllele));

		/* VFC: INFO: NS (number of samples with data) */
		fprintf (fp, ";NS=%d", allSites[j].countCellswithData);

		/* VFC: INFO: AF (alternate allele frequencies) */
		fprintf (fp, ";AF=");
		if (allSites[j].countACGT == allSites[j].countA ||
		allSites[j].countACGT == allSites[j].countC ||
		allSites[j].countACGT == allSites[j].countG ||
		allSites[j].countACGT == allSites[j].countT)
			fprintf (fp,".");
		else
			{
			if (allSites[j].countA > 0 && referenceAllele != A)
				fprintf (fp,"%4.3f,", (double) allSites[j].countA / allSites[j].countACGT);
			if (allSites[j].countC > 0 && referenceAllele != C)
				fprintf (fp,"%4.3f,", (double) allSites[j].countC / allSites[j].countACGT);
			if (allSites[j].countG > 0 && referenceAllele != G)
				fprintf (fp,"%4.3f,", (double) allSites[j].countG / allSites[j].countACGT);
			if (allSites[j].countT > 0 && referenceAllele != T)
				fprintf (fp,"%4.3f,", (double) allSites[j].countT / allSites[j].countACGT);
			fseek(fp, -1, SEEK_CUR); 	/* get rid of the last comma */
			}
	
		/* VFC: INFO : SOMATIC */
		fprintf (fp, ";SOMATIC");
		
		/* VFC: FORMAT */
		fprintf (fp, "\tGT:DP:RC:GL:ML:NG:DG:TG");
		
		for (i=0; i<numCells+1; i++)
			{
			maternalAllele = data[MATERNAL][i][j];
			paternalAllele = data[PATERNAL][i][j];
			trueMaternalAllele = cell[i].site[j].trueMaternalAllele;
			truePaternalAllele = cell[i].site[j].truePaternalAllele;
			maternalAlleleDoublet = cell[i].site[j].maternalAlleleDoublet;
			paternalAlleleDoublet= cell[i].site[j].paternalAlleleDoublet;

			/* VFC: FORMAT: GT (genotypes) */
			if (maternalAllele == ADO)
				fprintf (fp, "\t_");
			else if (maternalAllele == DELETION)
				fprintf (fp, "\t-");
			else if (maternalAllele == referenceAllele)
				fprintf (fp, "\t0");
			else
				{
				fprintf (fp,"\t");
				for (l=0; l<allSites[j].numAltAlleles; l++)
					{
					if (maternalAllele == allSites[j].alternateAlleles[l])
						{
						fprintf (fp,"%d", l+1);
						break;
						}
					}
				}
			fprintf (fp, "|");
			if (paternalAllele == ADO)
				fprintf (fp, "_");
			else if (paternalAllele == DELETION)
				fprintf (fp, "-");
			else if (paternalAllele == referenceAllele)
				fprintf (fp, "0");
			else
				{
				for (l=0; l<allSites[j].numAltAlleles; l++)
					{
					if (paternalAllele == allSites[j].alternateAlleles[l])
						{
						fprintf (fp,"%d", l+1);
						break;
						}
					}
				}
	
			if (cell[i].hasDoublet == NO)
				{
				/* VFC: FORMAT: DP (read depth) */
				fprintf (fp, ":%d", cell[i].site[j].numReads);

				/* VFC: FORMAT: RC (read counts ACGT) */
				fprintf (fp,":%d,%d,%d,%d",  cell[i].site[j].readCount[A], cell[i].site[j].readCount[C], cell[i].site[j].readCount[G], cell[i].site[j].readCount[T]);

				/* VFC: FORMAT: GL (genotype likelihoods) */
				fprintf (fp, ":");
				for (k=0; k<4; k++)
					for (l=k; l<4; l++)
						fprintf (fp, "%3.1f,", cell[i].site[j].scaledGenLike[k][l]);
				fseek(fp, -1, SEEK_CUR); 	/* rewind to get rid of the last comma */

				/* VFC: FORMAT: ML (maximum likelihood genotype) */
				if (cell[i].site[j].numReads > 0)
					{
					fprintf (fp, ":");
					fprintf (fp, "%c/%c",WhichNuc(cell[i].site[j].MLmatAllele), WhichNuc(cell[i].site[j].MLpatAllele));
					}
				else
					fprintf (fp, ":?/?");

				/* VFC: FORMAT: NG (NGS SNV genotype) */
				fprintf (fp, ":");
				fprintf (fp, "%c|%c",WhichNuc(maternalAllele), WhichNuc(paternalAllele));
				
				/* VFC: FORMAT: DG (true SNV second cell genotype, in doublet) */
				fprintf (fp, ":");
				fprintf (fp, ".");

				/* VFC: FORMAT: TG (true SNV genotype) */
				fprintf (fp, ":");
				fprintf (fp, "%c|%c",WhichNuc(trueMaternalAllele), WhichNuc(truePaternalAllele));
				
				}
			else
				{
				/* VFC: FORMAT: DP (read depth) */
				fprintf (fp, ":%d", cell[i].site[j].numReadsDoublet);
				
				/* VFC: FORMAT: RC (read counts ACGT) */
				fprintf (fp,":%d,%d,%d,%d",  cell[i].site[j].readCountDoublet[A], cell[i].site[j].readCountDoublet[C], cell[i].site[j].readCountDoublet[G], cell[i].site[j].readCountDoublet[T]);
				
				/* VFC: FORMAT: GL (genotype likelihoods) */
				fprintf (fp, ":");
				for (k=0; k<4; k++)
					for (l=k; l<4; l++)
						fprintf (fp, "%3.1f,", cell[i].site[j].scaledGenLikeDoublet[k][l]);
				fseek(fp, -1, SEEK_CUR); 	/* rewind to get rid of the last comma */
				
				/* VFC: FORMAT: ML (maximum likelihood genotype) */
				if (cell[i].site[j].numReadsDoublet > 0)
					{
					fprintf (fp, ":");
					fprintf (fp, "%c/%c",WhichNuc(cell[i].site[j].MLmatAlleleDoublet), WhichNuc(cell[i].site[j].MLpatAlleleDoublet));
					}
				else
					fprintf (fp, ":?/?");
				
				/* VFC: FORMAT: NG (NGS SNV genotype) */
				fprintf (fp, ":");
				fprintf (fp, "%c|%c",WhichNuc(maternalAllele), WhichNuc(paternalAllele));
				
				/* VFC: FORMAT: DG (true NGS SNV second cell genotype, in doublet) */
				fprintf (fp, ":");
				fprintf (fp, "%c|%c",WhichNuc(maternalAlleleDoublet), WhichNuc(paternalAlleleDoublet));

				/* VFC: FORMAT: TG (true SNV genotype) */
				fprintf (fp, ":");
				fprintf (fp, "%c|%c",WhichNuc(trueMaternalAllele), WhichNuc(truePaternalAllele));
				
				}


	
			} // cell
		} // snv
	
	fprintf (fp, "\n");
	}

/********************* PrepareGlobalFiles **********************/
/* Open global files to output results */

void PrepareGlobalFiles(int argc, char **argv)
    {
    char File[MAX_NAME];
		
    /* contains the simulated tree in Newick format */
    if (doPrintTree == YES)
        {
		sprintf(File,"%s/%s", resultsDir, treeFile);
        if ((fpTrees = fopen(File, "w")) == NULL)
            {
            fprintf (stderr, "Can't open \"%s\"\n", File);
            exit(-1);
            }
        }
		
    /* contains a list of times and branch lenghts for all nodes in the simulated trees */
    if (doPrintTimes == YES)
        {
		sprintf(File,"%s/%s", resultsDir, timesFile);
        if ((fpTimes = fopen(File, "w")) == NULL)
            {
            fprintf (stderr, "Can't open \"%s\"\n", File);
            exit(-1);
            }
        }
		
    if (doSimulateData == YES)
        {
        /* contains SNV genotypes for every cell */
        if (doPrintSNVgenotypes == YES)
            {
			sprintf(File,"%s/%s", resultsDir, SNVgenotypesFile);
             if ((fpSNVgenotypes = fopen(File, "w")) == NULL)
                {
				fprintf (stderr, "Can't open \"%s\"\n", File);
				exit(-1);
                }
            }
			
        /* contains haplotypes for variable sites for every cell */
        if (doPrintSNVhaplotypes == YES)
            {
			sprintf(File,"%s/%s", resultsDir, SNVhaplotypesFile);
            if ((fpSNVhaplotypes = fopen(File, "w")) == NULL)
                {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
                }
            }

        /* contains error-free haplotypes for variable sites for every cell */
        if (doPrintTrueHaplotypes == YES)
            {
			sprintf(File,"%s/%s", resultsDir, trueHaplotypesFile);
            if ((fpTrueHaplotypes = fopen(File, "w")) == NULL)
                {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
                }
            }

		/* contains ML haplotypes for variable sites for every cell */
		if (doPrintMLhaplotypes == YES)
			{
			sprintf(File,"%s/%s", resultsDir, MLhaplotypesFile);
			if ((fpMLhaplotypes = fopen(File, "w")) == NULL)
				{
				fprintf (stderr, "Can't open \"%s\"\n", File);
				exit(-1);
				}
			}

        /* contains all genotypes (variable or invariable) for every cell */
        if (doPrintFullGenotypes == YES)
            {
			sprintf(File,"%s/%s", resultsDir, fullGenotypesFile);
            if ((fpFullGenotypes = fopen(File, "w")) == NULL)
                {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
                }
            }

        /* contains haplotypes for all sites for every cell */
        if (doPrintFullHaplotypes == YES)
            {
			sprintf(File,"%s/%s", resultsDir, fullHaplotypesFile);
            if ((fpFullHaplotypes = fopen(File, "w")) == NULL)
                {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
                }
            }

        /* contains reads counts and genotype likelihoods for every SNV and cell */
        if (doSimulateReadCounts == YES)
            {
			sprintf(File,"%s/%s", resultsDir, VCFfile);
            if ((fpVCF = fopen(File, "w")) == NULL)
                {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
                }
            }
  
         /* contains reads counts for every SNV and cell */
        if (doPrintCATG == YES)
            {
			sprintf(File,"%s/%s", resultsDir, CATGfile);
            if ((fpCATG = fopen(File, "w")) == NULL)
                {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
                }
            }

       }

   if (doPrintSNVgenotypes == YES)
        {
        fprintf (fpSNVgenotypes, "%s - ",PROGRAM_NAME);
        PrintDate (fpSNVgenotypes);
        fprintf (fpSNVgenotypes, "SNV genotypes\n");
        PrintCommandLine (fpSNVgenotypes, argc, argv);
        fprintf (fpSNVgenotypes,"\n%d\n", numDataSets);
        }

    if (doPrintSNVhaplotypes == YES)
        {
        fprintf (fpSNVhaplotypes, "%s - ",PROGRAM_NAME);
        PrintDate (fpSNVhaplotypes);
        fprintf (fpSNVhaplotypes, "SNV haplotypes\n");
        PrintCommandLine (fpSNVhaplotypes, argc, argv);
        fprintf (fpSNVhaplotypes,"\n%d\n", numDataSets);
        }

     if (doPrintTrueHaplotypes == YES)
        {
        fprintf (fpTrueHaplotypes, "%s - ",PROGRAM_NAME);
        PrintDate (fpTrueHaplotypes);
        fprintf (fpTrueHaplotypes, "True haplotypes\n");
        PrintCommandLine (fpTrueHaplotypes, argc, argv);
        fprintf (fpTrueHaplotypes,"\n%d\n", numDataSets);
        }

	if (doPrintMLhaplotypes == YES)
		{
		fprintf (fpMLhaplotypes, "%s - ",PROGRAM_NAME);
		PrintDate (fpMLhaplotypes);
		fprintf (fpMLhaplotypes, "ML haplotypes\n");
		PrintCommandLine (fpMLhaplotypes, argc, argv);
		fprintf (fpMLhaplotypes,"\n%d\n", numDataSets);
		}

   if (doPrintFullGenotypes == YES)
        {
        fprintf (fpFullGenotypes, "%s - ",PROGRAM_NAME);
        PrintDate (fpFullGenotypes);
        fprintf (fpFullGenotypes, "Full genotypes\n");
        PrintCommandLine (fpFullGenotypes, argc, argv);
        fprintf (fpFullGenotypes,"\n%d\n", numDataSets);
        }

    if (doPrintFullHaplotypes == YES)
        {
        fprintf (fpFullHaplotypes, "%s - ",PROGRAM_NAME);
        PrintDate (fpFullHaplotypes);
        fprintf (fpFullHaplotypes, "Full haplotypes\n");
        PrintCommandLine (fpFullHaplotypes, argc, argv);
        fprintf (fpFullHaplotypes,"\n%d\n", numDataSets);
        }

	if (doSimulateReadCounts == YES)
        {
        fprintf (fpVCF, "%s - ",PROGRAM_NAME);
        PrintDate (fpVCF);
        fprintf (fpVCF, "Read counts and genotype likelihoods\n");
        PrintCommandLine (fpVCF, argc, argv);
        fprintf (fpVCF,"\n%d\n", numDataSets);
        }

	if (doPrintCATG == YES)
        {
        fprintf (fpCATG, "%s - ",PROGRAM_NAME);
        PrintDate (fpCATG);
        fprintf (fpCATG, "Read counts\n");
        PrintCommandLine (fpCATG, argc, argv);
        fprintf (fpCATG,"\n%d\n", numDataSets);
        }
	}



/********************* PrepareSeparateFiles **********************/
/* Open individual files to output results */

void PrepareSeparateFiles(int replicate)
    {
	/* contains the simulated tree in Newick format */
    if (doPrintTree == YES)
        {
		sprintf(File,"%s/%s", resultsDir, treeDir);
		mkdir(File,S_IRWXU);
		sprintf(File,"%s/%s/%s.%04d", resultsDir, treeDir, treeFile, replicate+1);
		if ((fpTrees = fopen(File, "w")) == NULL)
            {
            fprintf (stderr, "Can't open \"%s\"\n", File);
            exit(-1);
            }
        }
		
	if (doPrintTimes == YES)
        {
		sprintf(File,"%s/%s", resultsDir, timesDir);
		mkdir(File,S_IRWXU);
		sprintf(File,"%s/%s/%s.%04d", resultsDir, timesDir, timesFile, replicate+1);
		if ((fpTimes = fopen(File, "w")) == NULL)
            {
            fprintf (stderr, "Can't open \"%s\"\n", File);
            exit(-1);
            }
        }
		
    if (doSimulateData == YES)
        {
        /* contains SNV genotypes for every cell */
        if (doPrintSNVgenotypes == YES)
            {
			sprintf(File,"%s/%s", resultsDir, SNVgenotypesDir);
			mkdir(File,S_IRWXU);
			sprintf(File,"%s/%s/%s.%04d", resultsDir, SNVgenotypesDir, SNVgenotypesFile, replicate+1);
            if ((fpSNVgenotypes = fopen(File, "w")) == NULL)
                {
				fprintf (stderr, "Can't open \"%s\"\n", File);
				exit(-1);
                }
            }
			
        /* contains haplotypes for variable sites for every cell */
        if (doPrintSNVhaplotypes == YES)
            {
			sprintf(File,"%s/%s", resultsDir, SNVhaplotypesDir);
			mkdir(File,S_IRWXU);
			sprintf(File,"%s/%s/%s.%04d", resultsDir, SNVhaplotypesDir, SNVhaplotypesFile, replicate+1);
            if ((fpSNVhaplotypes = fopen(File, "w")) == NULL)
                {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
                }
            }

        /* contains reference haplotypes (before errors)  for every cell */
        if (doPrintTrueHaplotypes == YES)
            {
			sprintf(File,"%s/%s", resultsDir, trueHaplotypesDir);
			mkdir(File,S_IRWXU);
			sprintf(File,"%s/%s/%s.%04d", resultsDir, trueHaplotypesDir, trueHaplotypesFile, replicate+1);
            if ((fpTrueHaplotypes = fopen(File, "w")) == NULL)
                {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
                }
            }

		/* contains ML haplotypes  for every cell */
		if (doPrintMLhaplotypes == YES)
			{
			sprintf(File,"%s/%s", resultsDir, MLhaplotypesDir);
			mkdir(File,S_IRWXU);
			sprintf(File,"%s/%s/%s.%04d", resultsDir, MLhaplotypesDir, MLhaplotypesFile, replicate+1);
			if ((fpMLhaplotypes = fopen(File, "w")) == NULL)
				{
				fprintf (stderr, "Can't open \"%s\"\n", File);
				exit(-1);
				}
			}

        /* contains all genotypes (variable or invariable) for every cell */
        if (doPrintFullGenotypes == YES)
            {
			sprintf(File,"%s/%s", resultsDir, fullGenotypesDir);
			mkdir(File,S_IRWXU);
			sprintf(File,"%s/%s/%s.%04d", resultsDir, fullGenotypesDir, fullGenotypesFile, replicate+1);
            if ((fpFullGenotypes = fopen(File, "w")) == NULL)
                {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
                }
            }

        /* contains haplotypes for all sites for every cell */
        if (doPrintFullHaplotypes == YES)
            {
			sprintf(File,"%s/%s", resultsDir, fullHaplotypesDir);
			mkdir(File,S_IRWXU);
			sprintf(File,"%s/%s/%s.%04d", resultsDir, fullHaplotypesDir, fullHaplotypesFile, replicate+1);
            if ((fpFullHaplotypes = fopen(File, "w")) == NULL)
                {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
                }
            }

        /* contains reads counts and log10 normalized genotype likelihoods for every SNV and cell */
        if (doSimulateReadCounts == YES)
            {
			sprintf(File,"%s/%s", resultsDir, VCFdir);
			mkdir(File,S_IRWXU);
			sprintf(File,"%s/%s/%s.%04d", resultsDir, VCFdir, VCFfile, replicate+1);
            if ((fpVCF = fopen(File, "w")) == NULL)
                {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
                }
            }

        /* contains reads counts for every SNV and cell */
        if (doPrintCATG == YES)
            {
			sprintf(File,"%s/%s", resultsDir, CATGdir);
			mkdir(File,S_IRWXU);
			sprintf(File,"%s/%s/%s.%04d", resultsDir, CATGdir, CATGfile, replicate+1);
            if ((fpCATG = fopen(File, "w")) == NULL)
                {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
                }
            }
        }
	}



/***************************** ReadUserTree  *********************************/
/*	Reads a tree in Newick format, rooted or unroored, binary, and
with or without branch lengths */


/**************** Print Trees ***************/
/*  Print trees to treefile in Newick format */

void PrintTree (TreeNode *treeRoot, FILE *fp)
    {
    WriteTree (treeRoot, fp);
    fprintf (fp,");\n");
    }


/******************* WriteTree ****************/
/* Writes a given tree starting on a particular node (usually the root, but not necessarily) */

void WriteTree (TreeNode *p, FILE *fp)
    {
    if (p != NULL)
        {
        if(p->isHealthyTip == YES)
            {
            fprintf (fp, "healthycell:%8.6f",p->branchLength);
            }
        else if (p->left == NULL && p->right == NULL)
            {
			if (doUserTree == NO)
				fprintf (fp, "tumcell%04d:%8.6f", Index(p),p->branchLength);
			else
				fprintf (fp, "%s:%8.6f", p->name,p->branchLength);
            }
        else
            {
            fprintf (fp, "(");
            WriteTree (p->left, fp);
            fprintf (fp, ",");
            WriteTree (p->right, fp);
            if (p->anc !=NULL)
                fprintf (fp, "):%8.6f",p->branchLength);
            }
        }
    }


/********************* PrintTimes **********************/
/*	Prints to a given file a detailed description of
	the tree as lost of nodes, times, branch lengths.
 */

void PrintTimes(int listPosition, FILE *fp)
    {
    fprintf (fp, "\n-------------------- Nodes -------------------");
    fprintf (fp, "\n       class  label  index  (left right anc) |         time     time length    branch length");
    fprintf (fp, "\n--------------------------------------------------------------------------------------------\n");
    ListTimes (listPosition, fp);
    }


/********************** ListTimes ************************/
/*	Writes a given tree description from ListTimes
	It does not list superfluous nodes
 */

void ListTimes (int position, FILE *fp)
	{
    TreeNode	*p;
    do
    {
        p = treeNodes + position;
 
        if (p->isHealthyTip == YES)
            fprintf (fp, "%12s   %4d   %4d  (%4d %4d %4d) |   %10.2lf      %10.2lf       %10.4lf\n",
                     "healthyTip", Label(p), Index(p), Index(p->left), Index(p->right), Index(p->anc), p->time, p->length, p->branchLength);
        else if (p->anc != NULL && p->left != NULL && p->right != NULL && p->anc->anc == NULL)
            fprintf (fp, "%12s   %4d   %4d  (%4d %4d %4d) |   %10.2lf      %10.2lf       %10.4lf\n",
                     "ingroupMRCA", Label(p), Index(p), Index(p->left), Index(p->right), Index(p->anc), p->time, p->length, p->branchLength);
        else if (p->anc != NULL && p->left != NULL && p->right != NULL)
            fprintf (fp, "%12s   %4d   %4d  (%4d %4d %4d) |   %10.2lf      %10.2lf       %10.4lf\n",
                     "internal", Label(p), Index(p), Index(p->left), Index(p->right), Index(p->anc), p->time, p->length, p->branchLength);
        else if (p->anc != NULL && p->left == NULL && p->right == NULL)
            fprintf (fp, "%12s   %4d   %4d  (%4d %4d %4d) |   %10.2lf      %10.2lf       %10.4lf\n",
                     "tip", Label(p), Index(p), Index(p->left), Index(p->right), Index(p->anc), p->time, p->length, p->branchLength);
        else if (p->anc == NULL && p->left != NULL && p->right != NULL)
            fprintf (fp, "%12s   %4d   %4d  (%4d %4d %4d) |   %10.2lf      %10.2lf       %10.4lf\n",
                     "outgroupRoot", Label(p), Index(p), Index(p->left), Index(p->right), Index(p->anc), p->time, p->length, p->branchLength);
		
        position++;
		
        if (position>500000) /* just in case the tree is weird */
            {
            fprintf (stderr, "\n\nERROR: Problems when listing times in the tree\n\n");
            exit(-1);
            }
		}
    while (p->isHealthyTip == NO);
	}



/************************ PrintSiteInfo ***********************/
/* Prints info for the site struct */

static void PrintSiteInfo (FILE *fp, int i)
	{
	int j;
	fprintf (fp, "\n\nSite %d", i+1);
	fprintf (fp, "\n isSNV = %d", allSites[i].isSNV);
	fprintf (fp, "\n isSNP = %d", allSites[i].isSNP);
	fprintf (fp, "\n numMutations = %d", allSites[i].numMutations);
	fprintf (fp, "\n  numMutationsMaternal = %d", allSites[i].numMutationsMaternal);
	fprintf (fp, "\n  numMutationsPaternal = %d", allSites[i].numMutationsPaternal);
	fprintf (fp, "\n numDeletions = %d", allSites[i].numDeletions);
	fprintf (fp, "\n  numDeletionsMaternal = %d", allSites[i].numDeletionsMaternal);
	fprintf (fp, "\n  numDeletionsPaternal = %d", allSites[i].numDeletionsPaternal);
	fprintf (fp, "\n numCNLOH = %d", allSites[i].numCNLOH);
	fprintf (fp, "\n  numCNLOHmaternal = %d", allSites[i].numCNLOHmaternal);
	fprintf (fp, "\n  numCNLOHpaternal = %d", allSites[i].numCNLOHpaternal);
	fprintf (fp, "\n hasADO = %d", allSites[i].hasADO);
	fprintf (fp, "\n hasGenotypingError = %d", allSites[i].hasGenotypeError);
	fprintf (fp, "\n referenceAllele = %c", WhichNuc(allSites[i].referenceAllele));
	fprintf (fp, "\n numAltAlleles = %d   (", allSites[i].numAltAlleles);
	for (j=0; j <allSites[i].numAltAlleles; j++)
		fprintf (fp, " %c", WhichNuc(allSites[i].alternateAlleles[j]));
	fprintf (fp, " )");
	fprintf (fp, "\n countA = %d", allSites[i].countA);
	fprintf (fp, "\n countC = %d", allSites[i].countC);
	fprintf (fp, "\n countG = %d", allSites[i].countG);
	fprintf (fp, "\n countT = %d", allSites[i].countT);
	fprintf (fp, "\n countDropped = %d", allSites[i].countDropped);
	fprintf (fp, "\n rateMultiplier = %f", allSites[i].rateMultiplier);
}


/************************ PrintSNVGenotypes ***********************/
/* Prints genotypes at variable sites (SNVs) to a file */

static void PrintSNVGenotypes (FILE *fp)
    {
    int		i, j, k;

    if (doPrintAncestors == YES)
        fprintf (fp, "%d %d\n", numCells+3, numSNVs);
    else
        fprintf (fp, "%d %d\n", numCells+1, numSNVs);
		
    /* site information */
    /* fprintf (fp, "%14s", ""); */
	for (i=0; i<numSNVs; i++)
		fprintf (fp, "%d ", SNVsites[i]+1);
    fseek (fp, -1, SEEK_CUR);
    fprintf (fp, "\n");
	k=0;
	
    if (alphabet == DNA)
        {
        for (i=0; i<numCells+1; i++)
            {
            if (i == numCells)
				fprintf (fp,"healthycell ");
			else
				{
				if (doUserTree == NO)
					fprintf (fp,"tumcell%04d ", i+1);
				else
					fprintf (fp,"%-12s", cellNames[i]);
				}
            for (j=0; j<numSNVs; j++)
				fprintf (fp, " %c%c", WhichNuc(data[MATERNAL][i][SNVsites[j]]),WhichNuc(data[PATERNAL][i][SNVsites[j]]));
            fprintf (fp,"\n");
            }
			
        if (doPrintAncestors == YES)
            {
            fprintf (fp,"hearoot%04d ", i+1);
			for (j=0; j<numSNVs; j++)
				fprintf (fp, " %c%c", WhichNuc(data[MATERNAL][HEALTHY_ROOT][SNVsites[j]]),WhichNuc(data[PATERNAL][HEALTHY_ROOT][SNVsites[j]]));
            fprintf (fp,"\ntumroot%04d ", i+1);
			for (j=0; j<numSNVs; j++)
				fprintf (fp, " %c%c", WhichNuc(data[MATERNAL][TUMOR_ROOT][SNVsites[j]]),WhichNuc(data[PATERNAL][TUMOR_ROOT][SNVsites[j]]));
            fprintf (fp,"\n");
            }
        }
    else
        {
        for (i=0; i<numCells+1; i++)
            {
			if (i == numCells)
				fprintf (fp,"healthycell ");
			else
				{
				if (doUserTree == NO)
					fprintf (fp,"tumcell%04d ", i+1);
				else
					fprintf (fp,"%-12s", cellNames[i]);
				}
			for (j=0; j<numSNVs; j++)
				fprintf (fp, " %c%c", WhichMut(data[MATERNAL][i][SNVsites[j]]),WhichMut(data[PATERNAL][i][SNVsites[j]]));
			fprintf (fp,"\n");
            }
			
        if (doPrintAncestors == YES)
            {
			fprintf (fp,"hearoot%04d ", i+1);
			for (j=0; j<numSNVs; j++)
				fprintf (fp, " %d%d", data[MATERNAL][HEALTHY_ROOT][SNVsites[j]],data[PATERNAL][HEALTHY_ROOT][SNVsites[j]]);
			fprintf (fp,"\ntumroot%04d ", i+1);
			for (j=0; j<numSNVs; j++)
				fprintf (fp, " %d%d", data[MATERNAL][TUMOR_ROOT][SNVsites[j]],data[PATERNAL][TUMOR_ROOT][SNVsites[j]]);
			fprintf (fp,"\n");
            }
        }
     }

/***************************** PrintSNVHaplotypes *******************************/
/* Prints haplotypes at variable sites (SNVs), to a file */

static void PrintSNVHaplotypes (FILE *fp, int PrintTrueVariants)
    {
    int		i, j;
    int		*sites;
	
	if (PrintTrueVariants == YES)
		sites = variantSites;
	else
		sites = SNVsites;

	if (doPrintIUPAChaplotypes == YES)
		{
		if (doPrintAncestors == YES)
			fprintf (fp,"%d %d\n", numCells+3, numSNVs);
		else
			fprintf (fp,"%d %d\n", numCells+1, numSNVs);
		}
	else
		{
		if (doPrintAncestors == YES)
			fprintf (fp,"%d %d\n",2*(numCells+3), numSNVs);
		else
			fprintf (fp,"%d %d\n",2*(numCells+1), numSNVs);
		}
		
    /* site information */
	for (i=0; i<numSNVs; i++)
		fprintf (fp, "%d ", sites[i]+1);
	fseek(fp, -1, SEEK_CUR);
	fprintf (fp, "\n");

    if (alphabet == DNA)
        {
		if (doPrintIUPAChaplotypes == YES)
			{
			for (i=0; i<numCells+1; i++)
				{
				/* print haplotype */
				if (i == numCells)
					fprintf (fp,"healthycell  ");
				else
					{
					if (doUserTree == NO)
						fprintf (fp,"tumcell%04d  ", i+1);
					else
						fprintf (fp,"%-12s ", cellNames[i]);
					}
				for (j=0; j<numSNVs; j++)
					fprintf (fp, "%c", WhichIUPAC(data[MATERNAL][i][sites[j]],data[PATERNAL][i][sites[j]]));
				fprintf (fp,"\n");
				}
			
			if (doPrintAncestors == YES)
				{
				fprintf (fp,"hearoot%04d  ", i+1);
				for (j=0; j<numSNVs; j++)
					fprintf (fp, "%c", WhichIUPAC(data[MATERNAL][HEALTHY_ROOT][sites[j]],data[PATERNAL][HEALTHY_ROOT][sites[j]]));
					
				fprintf (fp,"\ntumroot%04d  ", i+1);
				for (j=0; j<numSNVs; j++)
					fprintf (fp, "%c", WhichIUPAC(data[MATERNAL][TUMOR_ROOT][sites[j]],data[PATERNAL][TUMOR_ROOT][sites[j]]));
				fprintf (fp,"\n");
				}
			}
		else // print maternal and paternal DNA haplotypes
			{
			for (i=0; i<numCells+1; i++)
				{
				/* print maternal haplotype */
				if (i == numCells)
					fprintf (fp,"healthycellm ");
				else
					{
					if (doUserTree == NO)
						fprintf (fp,"tumcell%04dm ", i+1);
					else
						fprintf (fp,"m%-12s", cellNames[i]);
					}
				for (j=0; j<numSNVs; j++)
					fprintf (fp, "%c", WhichNuc(data[MATERNAL][i][sites[j]]));
				fprintf (fp,"\n");
				
				/* print paternal haplotype */
				if (i == numCells)
					fprintf (fp,"healthycellp ");
				else
					{
					if (doUserTree == NO)
						fprintf (fp,"tumcell%04dp ", i+1);
					else
						fprintf (fp,"p%-12s", cellNames[i]);
					}
				for (j=0; j<numSNVs; j++)
					fprintf (fp, "%c", WhichNuc(data[PATERNAL][i][sites[j]]));
				fprintf (fp,"\n");
				}
			
			if (doPrintAncestors == YES)
				{
				fprintf (fp,"hearoot%04dm ", i+1);
				for (j=0; j<numSNVs; j++)
					fprintf (fp, "%c", WhichNuc(data[MATERNAL][HEALTHY_ROOT][sites[j]]));
				fprintf (fp,"\nhearoot%04dp ", i+1);
				for (j=0; j<numSNVs; j++)
					fprintf (fp, "%c", WhichNuc(data[PATERNAL][HEALTHY_ROOT][sites[j]]));
					
				fprintf (fp,"\ntumroot%04dm ", i+1);
				for (j=0; j<numSNVs; j++)
					fprintf (fp, "%c", WhichNuc(data[MATERNAL][TUMOR_ROOT][sites[j]]));
				fprintf (fp,"\ntumroot%04dp ", i+1);
				for (j=0; j<numSNVs; j++)
					fprintf (fp, "%c", WhichNuc(data[PATERNAL][TUMOR_ROOT][sites[j]]));
				fprintf (fp,"\n");
				}
			}
		}
    else  //print binary haplotypes
        {
		if (doPrintIUPAChaplotypes == YES) //print binary consensus SNV haplotypes
			{
			for (i=0; i<numCells+1; i++)
				{
				if (i == numCells)
					fprintf (fp,"healthycell  ");
				else
					{
					if (doUserTree == NO)
						fprintf (fp,"tumcell%04d ", i+1);
					else
						fprintf (fp,"%-12s ", cellNames[i]);
					}
				for (j=0; j<numSNVs; j++)
					fprintf (fp, "%c", WhichConsensusBinary(data[MATERNAL][i][sites[j]],data[PATERNAL][i][sites[j]]));
				fprintf (fp,"\n");
				}

			if (doPrintAncestors == YES)
				{
				fprintf (fp,"hearoot%04d  ", i+1);
				for (j=0; j<numSNVs; j++)
					fprintf (fp, "%c", WhichConsensusBinary(data[MATERNAL][HEALTHY_ROOT][sites[j]],data[PATERNAL][HEALTHY_ROOT][sites[j]]));
					
				fprintf (fp,"\ntumroot%04d  ", i+1);
				for (j=0; j<numSNVs; j++)
					fprintf (fp, "%c", WhichConsensusBinary(data[MATERNAL][TUMOR_ROOT][sites[j]],data[PATERNAL][TUMOR_ROOT][sites[j]]));
				fprintf (fp,"\n");
				}
			}
		else // print maternal and paternal binary SNV haplotypes
			{
			for (i=0; i<numCells+1; i++)
				{
				/* print maternal haplotype */
				if (i == numCells)
					fprintf (fp,"healthycellm ");
				else
					{
					if (doUserTree == NO)
						fprintf (fp,"tumcell%04dm ", i+1);
					else
						fprintf (fp,"m%-12s", cellNames[i]);
					}
				for (j=0; j<numSNVs; j++)
					fprintf (fp, "%c", WhichMut(data[MATERNAL][i][sites[j]]));
				fprintf (fp,"\n");
					
				/* print paternal haplotype */
				if (i == numCells)
					fprintf (fp,"healthycellp ");
				else
					{
					if (doUserTree == NO)
						fprintf (fp,"tumcell%04dp ", i+1);
					else
						fprintf (fp,"p%-12s", cellNames[i]);
					}
				for (j=0; j<numSNVs; j++)
					fprintf (fp, "%c", WhichMut(data[PATERNAL][i][sites[j]]));
					
				fprintf (fp,"\n");
				}

			if (doPrintAncestors == YES)
				{
				fprintf (fp,"hearoot%04dm ", i+1);
				for (j=0; j<numSNVs; j++)
					fprintf (fp, "%c", WhichMut(data[MATERNAL][HEALTHY_ROOT][sites[j]]));
				fprintf (fp,"\nhearoot%04dp ", i+1);
				for (j=0; j<numSNVs; j++)
					fprintf (fp, "%c", WhichMut(data[PATERNAL][HEALTHY_ROOT][sites[j]]));
					
				fprintf (fp,"\ntumroot%04dm ", i+1);
				for (j=0; j<numSNVs; j++)
					fprintf (fp, "%c", WhichMut(data[MATERNAL][TUMOR_ROOT][sites[j]]));
				fprintf (fp,"\ntumroot%04dp ", i+1);
				for (j=0; j<numSNVs; j++)
					fprintf (fp, "%c", WhichMut(data[PATERNAL][TUMOR_ROOT][sites[j]]));
				fprintf (fp,"\n");
				}
			}
        }
    }


/********************* PrintFullGenotypes ********************/
/* Prints genotypes at all sites to a file */

static void PrintFullGenotypes (FILE *fp)
    {
    int		i, j;
		
    if (doPrintAncestors == YES)
        fprintf (fp, "%d %d\n", numCells+3, numSites);
    else
        fprintf (fp, "%d %d\n", numCells+1, numSites);

    if (alphabet == DNA)
        {
        for (i=0; i<numCells+1; i++)
            {
            if (i == numCells)
                fprintf (fp,"healthycell ");
			else
				{
				if (doUserTree == NO)
					fprintf (fp,"tumcell%04d ", i+1);
				else
					fprintf (fp,"%-12s", cellNames[i]);
				}
           for (j=0; j<numSites; j++)
                fprintf (fp, " %c%c", WhichNuc(data[0][i][j]), WhichNuc(data[1][i][j]));
            fprintf (fp,"\n");
            }
			
        if (doPrintAncestors == YES)
            {
            fprintf (fp,"hearoot%04d ", 0);
            for (j=0; j<numSites; j++)
                fprintf (fp, " %c%c",  WhichNuc(data[MATERNAL][HEALTHY_ROOT][j]), WhichNuc(data[PATERNAL][HEALTHY_ROOT][j]));
            fprintf (fp,"\ntumroot%04d ", 0);
            for (j=0; j<numSites; j++)
                fprintf (fp, " %c%c",  WhichNuc(data[MATERNAL][TUMOR_ROOT][j]), WhichNuc(data[PATERNAL][TUMOR_ROOT][j]));
            fprintf (fp,"\n");
            }
        }
    else
        {
        for (i=0; i<numCells+1; i++)
            {
            if (i == numCells)
                fprintf (fp,"healthycell ");
 			else
				{
				if (doUserTree == NO)
					fprintf (fp,"tumcell%04d ", i+1);
				else
					fprintf (fp,"%-12s", cellNames[i]);
				}
			for (j=0; j<numSites; j++)
                fprintf (fp, " %c%c", WhichMut(data[0][i][j]),WhichMut(data[1][i][j]));
            fprintf (fp,"\n");
        }
			
        if (doPrintAncestors == YES)
            {
            fprintf (fp,"hearoot%04d ", 0);
            for (j=0; j<numSites; j++)
                fprintf (fp, " %d%d", data[MATERNAL][HEALTHY_ROOT][j],data[PATERNAL][HEALTHY_ROOT][j]);
            fprintf (fp,"\ntumroot%04d ", 0);
            for (j=0; j<numSites; j++)
                fprintf (fp, " %d%d", data[MATERNAL][TUMOR_ROOT][j],data[PATERNAL][TUMOR_ROOT][j]);
            fprintf (fp,"\n");
            }
        }
     }

/***************************** PrintFullHaplotypes *******************************/
/* Prints cells with all sites (variable + invariable) to a file */

static void PrintFullHaplotypes (FILE *fp)
    {
    int		 i, j;
		
	if (doPrintIUPAChaplotypes == NO)
		{
		if (doPrintAncestors == YES)
			fprintf (fp,"%d %d\n",2*(numCells+3), numSites);
		else
			fprintf (fp,"%d %d\n",2*(numCells+1), numSites);
		}
	else
		{
		if (doPrintAncestors == YES)
			fprintf (fp,"%d %d\n", numCells+3, numSites);
		else
			fprintf (fp,"%d %d\n", numCells+1, numSites);
		}
		
    if (alphabet == DNA)
        {
		if (doPrintIUPAChaplotypes == YES)
			{
			for (i=0; i<numCells+1; i++)
				{
				/* print haplotype */
				if (i == numCells)
					fprintf (fp,"healthycell  ");
				else
					{
					if (doUserTree == NO)
						fprintf (fp,"tumcell%04d  ", i+1);
					else
						fprintf (fp,"%-12s ", cellNames[i]);
					}
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichIUPAC(data[MATERNAL][i][j],data[PATERNAL][i][j]));
				fprintf (fp,"\n");
				}

			if (doPrintAncestors == YES)
				{
				fprintf (fp,"hearoot%04d  ", 0);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichIUPAC(data[MATERNAL][HEALTHY_ROOT][j],data[PATERNAL][HEALTHY_ROOT][j]));
					
				fprintf (fp,"\ntumroot%04d  ", 0);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichIUPAC(data[MATERNAL][TUMOR_ROOT][j],data[PATERNAL][TUMOR_ROOT][j]));
				fprintf (fp,"\n");
				}
			}
		else // print maternal and paternal DNA haplotypes
			{
			for (i=0; i<numCells+1; i++)
				{
				/* print maternal haplotype */
				if (i == numCells)
					fprintf (fp,"healthycellm ");
				else
					{
					if (doUserTree == NO)
						fprintf (fp,"tumcell%04dm ", i+1);
					else
						fprintf (fp,"m%-12s", cellNames[i]);
					}
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichNuc(data[MATERNAL][i][j]));
				fprintf (fp,"\n");
				
				/* print paternal haplotype */
				if (i == numCells)
					fprintf (fp,"healthycellp ");
				else
					{
					if (doUserTree == NO)
						fprintf (fp,"tumcell%04dp ", i+1);
					else
						fprintf (fp,"p%-12s", cellNames[i]);
					}
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichNuc(data[PATERNAL][i][j]));
				fprintf (fp,"\n");
				}

			if (doPrintAncestors == YES)
				{
				fprintf (fp,"hearoot%04dm ", 0);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichNuc(data[MATERNAL][HEALTHY_ROOT][j]));
				fprintf (fp,"\nhearoot%04dp ", 0);
					for (j=0; j<numSites; j++)
						fprintf (fp, "%c", WhichNuc(data[PATERNAL][HEALTHY_ROOT][j]));
					
				fprintf (fp,"\ntumroot%04dm ", 0);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichNuc(data[MATERNAL][TUMOR_ROOT][j]));
				fprintf (fp,"\ntumroot%04dp ", 0);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichNuc(data[PATERNAL][TUMOR_ROOT][j]));
				fprintf (fp,"\n");
				}
			}
		}
	else  //print binary haplotypes
        {
		if (doPrintIUPAChaplotypes == YES) //print binary consensus haplotypes
			{
			for (i=0; i<numCells+1; i++)
				{
				if (i == numCells)
					fprintf (fp,"healthycell  ");
				else
					{
					if (doUserTree == NO)
						fprintf (fp,"tumcell%04d ", i+1);
					else
						fprintf (fp,"%-12s ", cellNames[i]);
					}
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichConsensusBinary(data[MATERNAL][i][j],data[PATERNAL][i][j]));
				fprintf (fp,"\n");
				}
			
			if (doPrintAncestors == YES)
				{
				fprintf (fp,"hearoot%04d  ", 0);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichConsensusBinary(data[MATERNAL][HEALTHY_ROOT][j],data[PATERNAL][HEALTHY_ROOT][j]));
					
				fprintf (fp,"\ntumroot%04d  ", 0);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichConsensusBinary(data[MATERNAL][TUMOR_ROOT][j],data[PATERNAL][TUMOR_ROOT][j]));
				fprintf (fp,"\n");
				}
			}
		else // print maternal and paternal binary haplotypes
			{
			for (i=0; i<numCells+1; i++)
				{
				/* print maternal haplotype */
				if (i == numCells)
					fprintf (fp,"healthycellm ");
				else
					{
					if (doUserTree == NO)
						fprintf (fp,"tumcell%04dm ", i+1);
					else
						fprintf (fp,"m%-12s", cellNames[i]);
					}
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichMut(data[MATERNAL][i][j]));
				fprintf (fp,"\n");
					
				/* print paternal haplotype */
				if (i == numCells)
					fprintf (fp,"healthycellp ");
				else
					{
					if (doUserTree == NO)
						fprintf (fp,"tumcell%04dp ", i+1);
					else
						fprintf (fp,"p%-12s", cellNames[i]);
					}
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichMut(data[PATERNAL][i][j]));
				fprintf (fp,"\n");
				}
			
			if (doPrintAncestors == YES)
				{
				fprintf (fp,"hearoot%04dm ", 0);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichMut(data[MATERNAL][HEALTHY_ROOT][j]));
				fprintf (fp,"\nhearoot%04dp ", 0);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichMut(data[PATERNAL][HEALTHY_ROOT][j]));
					
				fprintf (fp,"\ntumroot%04dm ", 0);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichMut(data[MATERNAL][TUMOR_ROOT][j]));
				fprintf (fp,"\ntumroot%04dp ", 0);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichMut(data[PATERNAL][TUMOR_ROOT][j]));
				fprintf (fp,"\n");
				}
			}
        }
    }


/***************************** PrintMLhaplotypes *******************************/
/* Prints cells ML DNA haplotypes to a file (invariable sites are not ML estimates)  */

static void PrintMLhaplotypes (FILE *fp)
{
	int		 i, j;
	
	if (doPrintIUPAChaplotypes == NO)
		fprintf (fp,"%d %d\n",2*(numCells+1), numSites);
	else
		fprintf (fp,"%d %d\n", numCells+1, numSites);
	
	if (doPrintIUPAChaplotypes == YES)
		{
		for (i=0; i<numCells+1; i++)
			{
			if (i == numCells)
				fprintf (fp,"healthycell  ");
			else
				{
				if (doUserTree == NO)
					fprintf (fp,"tumcell%04d  ", i+1);
				else
					fprintf (fp,"%-12s ", cellNames[i]);
				}
			for (j=0; j<numSites; j++)
				{
				if (allSites[j].isSNV == YES)
					fprintf (fp, "%c", WhichIUPAC(cell[i].site[j].MLmatAllele, cell[i].site[j].MLpatAllele));
				else
					fprintf (fp, "%c", WhichIUPAC(data[MATERNAL][i][j],data[PATERNAL][i][j]));
				}
			fprintf (fp,"\n");
			}
		}
	else // print maternal and paternal DNA haplotypes
		{
		for (i=0; i<numCells+1; i++)
			{
			/* print maternal haplotype */
			if (i == numCells)
				fprintf (fp,"healthycellm ");
			else
				{
				if (doUserTree == NO)
					fprintf (fp,"tumcell%04dm ", i+1);
				else
					fprintf (fp,"m%-12s", cellNames[i]);
				}
			for (j=0; j<numSites; j++)
				{
				if (allSites[j].isSNV == YES)
					fprintf (fp, "%c", WhichNuc(cell[i].site[j].MLmatAllele));
				else
					fprintf (fp, "%c", WhichNuc(data[MATERNAL][i][j]));
				}
			fprintf (fp,"\n");
			
			/* print paternal haplotype */
			if (i == numCells)
				fprintf (fp,"healthycellp ");
			else
				{
				if (doUserTree == NO)
					fprintf (fp,"tumcell%04dp ", i+1);
				else
					fprintf (fp,"p%-12s", cellNames[i]);
				}
			for (j=0; j<numSites; j++)
				{
				if (allSites[j].isSNV == YES)
					fprintf (fp, "%c", WhichNuc(cell[i].site[j].MLpatAllele));
				else
					fprintf (fp, "%c", WhichNuc(data[PATERNAL][i][j]));
				}
			fprintf (fp,"\n");
			}
		}
	}



/********************* WhichIUPAC ************************/
/* Returns the IUPAC representation of the genotype */
/*
 UPAC nucleotide code	Base
 A	Adenine
 C	Cytosine
 G	Guanine
 T (or U)	Thymine (or Uracil)
 R	A or G
 Y	C or T
 S	G or C
 W	A or T
 K	G or T
 M	A or C
 B	C or G or T
 D	A or G or T
 H	A or C or T
 V	A or C or G
 N	any base
 . or -	gap

This is what we do:

 A/A => A
 A/C => M
 A/G => R
 A/T => W
 A/_ => a
 
 C/A => M
 C/C => C
 C/G => S
 C/T => Y
 C/_ => c

 G/A => R
 G/C => S
 G/G => G
 G/T => K
 G/_ => g

 T/A => W
 T/C => Y
 T/G => K
 T/T => T
 T/_ => t

 _/A => a
 _/C => c
 _/G => g
 _/T => t
 _/_ => -

 */

char WhichIUPAC (int allele1, int allele2)
	{
	if (allele1 == 0)
		{
		if (allele2 == 0)		//AA
			return ('A');
		else if (allele2 == 1)	//AC
			return ('M');
		else if (allele2 == 2)	//AG
			return ('R');
		else if (allele2 == 3)	//AT
			return ('W');
		else if (allele2 == ADO)	//A?
			return ('a');
		else if (allele2 == DELETION)	//A–
			return ('a');
		else
			return ('N');
		}
	else if (allele1 == 1)
		{
		if (allele2 == 0)		//CA
			return ('M');
		else if (allele2 == 1)	//CC
			return ('C');
		else if (allele2 == 2)	//CG
			return ('S');
		else if (allele2 == 3)	//CT
			return ('Y');
		else if (allele2 == ADO)	//C?
			return ('c');
		else if (allele2 == DELETION)	//C–
			return ('c');
		else
			return ('N');
		}
	else if (allele1 == 2)
		{
		if (allele2 == 0)		//GA
			return ('R');
		else if (allele2 == 1)	//GC
			return ('S');
		else if (allele2 == 2)	//GG
			return ('G');
		else if (allele2 == 3)	//GT
			return ('K');
		else if (allele2 == ADO)	//G?
			return ('g');
		else if (allele2 == DELETION)	//G–
			return ('g');
		else
			return ('N');
		}
	else if (allele1 == 3)
		{
		if (allele2 == 0)		//TA
			return ('W');
		else if (allele2 == 1)	//TC
			return ('Y');
		else if (allele2 == 2)	//TG
			return ('K');
		else if (allele2 == 3)	//TT
			return ('T');
		else if (allele2 == ADO)	//T?
			return ('t');
		else if (allele2 == DELETION)	//T–
			return ('t');
		else
			return ('N');
		}
	else if (allele1 == ADO)
		{
		if (allele2 == 0)		//?A
			return ('a');
		else if (allele2 == 1)	//?C
			return ('c');
		else if (allele2 == 2)	//?G
			return ('g');
		else if (allele2 == 3)	//?T
			return ('t');
		else if (allele2 == ADO)	//??
			return ('-');
		else if (allele2 == DELETION)	//?-
			return ('-');
		else
			return ('N');
		}
	else if (allele1 == DELETION)
		{
		if (allele2 == 0)		//-A
			return ('a');
		else if (allele2 == 1)	//-C
			return ('c');
		else if (allele2 == 2)	//-G
			return ('g');
		else if (allele2 == 3)	//-T
			return ('t');
		else if (allele2 == ADO)	//-?
			return ('-');
		else if (allele2 == DELETION)	//--
			return ('-');
		else
			return ('N');
		}
	else
		return ('N');
	}





/********************* WhichConsensusBinary ************************/
/* Returns a consensus representation of the binary genotype */
/*
 0/0 => 0
 0/1 => 1
 1/0 => 1
 1/1 => 2

 0/_ => 0
 _/0 => 0

 1/_ => 2
 _/1 => 2

 _/_ => -
 */

char WhichConsensusBinary (int allele1, int allele2)
	{
	if (allele1 == 0)
		{
		if (allele2 == 0)		//00
			return ('0');
		else if (allele2 == 1)	//01
			return ('1');
		else if (allele2 == ADO)	//0?
			return ('0');
		else if (allele2 == DELETION)	//0-
			return ('0');
		else
			return ('N');
		}
	else if (allele1 == 1)
		{
		if (allele2 == 0)		//10
			return ('1');
		else if (allele2 == 1)	//11
			return ('2');
		else if (allele2 == ADO)	//1?
			return ('2');
		else if (allele2 == DELETION)	//0-
			return ('2');
		else
			return ('N');
		}
	else if (allele1 == ADO)
		{
		if (allele2 == 0)		//?0
			return ('0');
		else if (allele2 == 1)	//?1
			return ('2');
		else if (allele2 == ADO)	//??
			return ('-');
		else if (allele2 == DELETION)	//?-
			return ('-');
		else
			return ('N');
		}
	else if (allele1 == DELETION)
		{
		if (allele2 == 0)		//-0
			return ('0');
		else if (allele2 == 1)	//-1
			return ('2');
		else if (allele2 == ADO)	//-?
			return ('-');
		else if (allele2 == DELETION)	//--
			return ('-');
		else
			return ('N');
		}
	else
		return ('N');
	}

/********************* WhichNucWhichNucChar ************************/
/* Returns integer representation for character nucleotudes */

int WhichNucChar (char nucleotide)
    {
    if (nucleotide == 'A')
        return (A);
    else if (nucleotide == 'C')
        return (C);
    else if (nucleotide == 'G')
        return (G);
    else if (nucleotide == 'T')
        return (T);
    else if (nucleotide == '?')
        return (ADO);
    else if (nucleotide == '-')
        return (DELETION);
    else
        {
 		fprintf (stderr, "\nERROR in WhicNucChar: nucleotide = %c\n",  nucleotide);
		exit(-1);
        }
    }




/********************* WhichNuc ************************/
/* Returns character representation for nucleotides */

char WhichNuc (int nucleotide)
    {
    if (nucleotide == A)
        return ('A');
    else if (nucleotide == C)
        return ('C');
    else if (nucleotide == G)
        return ('G');
    else if (nucleotide == T)
        return ('T');
    else if (nucleotide == ADO)
        return ('?');
    else if (nucleotide == DELETION)
        return ('-');
    else
        {
 		fprintf (stderr, "\nERROR in WhicNuc: nucleotide = %d\n",  nucleotide);
		exit(-1);
	    return ('N');
        }
    }


/********************* WhichMut ************************/
/* Returns character representation for binary data */

char WhichMut (int state)
    {
    if (state == 0)
        return ('0');
    else if (state == 1)
        return ('1');
    else if (state == ADO)
        return ('?');
     else if (state == DELETION)
        return ('-');
   else
        return ('N');
    }



/***************************** PrintHeader *******************************/
/* Prints program header */

static void PrintHeader(FILE *fp)
    {
    fprintf (fp,"Cell coalescent simulation - %s", PROGRAM_NAME);
    fprintf (fp,"  %s", VERSION_NUMBER);
    fprintf (fp,"\n (c) 2018 David Posada - dposada@uvigo.es");
    fprintf (fp,"\n__________________________________________________\n\n");
    }



/***************************** PrintDate *******************************/
/* Prints current date and time  */

static void PrintDate (FILE *fp)
    {
    time_t now;
    char *date;
		
    now=time(NULL);
    date= ctime(&now);
    fprintf (fp,"%s",date);
    }


/***************************** PrintRunInformation *******************************/
/* Prints a summary of run settings */

static void	PrintRunInformation (FILE *fp)
    {
    int		i;
		
    fprintf (fp, "\n\nRun settings\n------------");
    if (noisy > 1)
		fprintf (fp, "\n[Assumptions in brackets]\n");
    fprintf (fp, "\n Seed                                         =   %-3ld", originalSeed);
    fprintf (fp, "\n Number replicate data sets                   =   %-3d",  numDataSets);
    fprintf (fp, "\n Number of cells                              =   %-3d",  numCells);
	if (doUserGenome == YES)
		fprintf (fp, "\n Number of sites in the user genome           =   %-3d",  numSites);
	else
		fprintf (fp, "\n Number of genomic sites                      =   %-3d",  numSites);
		
    fprintf (fp, "\n\nDemographics");
    fprintf (fp, "\n Effective population size                    =   %-3d",  N);
    if (doExponential == YES)
        fprintf (fp, "\n Exponential growth rate                      =   %2.1e", growthRate);
    else if (doDemographics == YES)
        {
        fprintf (fp, "\n Period   Nbegin     Nend    Duration    Growth rate");
        for (i=1; i<=numPeriods; i++)
            fprintf (fp, "\n   %d    %7d   %7d     %7d       %+2.1e",
                     i, Nbegin[i], Nend[i], cumDuration[i]-cumDuration[i-1], periodGrowth[i]);
        }
	if (doUserTree == YES)
		fprintf (fp, "\n\nUser tree file                                =   %s", userTreeFile);
	else
		{
		fprintf (fp, "\n\nCoalescent");
	   // fprintf (fp, "\n Mean number of coalescence events            =   %3.2f", meanNumCA);
		fprintf (fp, "\n Mean time to the       MRCA (# generations)  =   %3.2f", meanTMRCA);
		if (noisy > 1)
			fprintf (fp, "\n Exp time to the       MRCA [constant Ne]     =   %3.2f", expTMRCA);
		if (numDataSets > 1)
			{
			fprintf (fp, "\n Variance time       MRCA                     =   %3.2f", varTMRCA);
			if (noisy > 1)
				fprintf (fp, "\n Exp variance time       MRCA [constant Ne]   =   %3.2f", expVarTMRCA);
			}
	 
		if (rateVarAmongLineages == YES)
			fprintf (fp, "\n Rate variation among branches (alpha)        =   %-3.2f", alphaBranches);
		else
			fprintf (fp, "\n No rate variation among branches");
		}
		
    fprintf (fp, "\n\nUser-defined branches");
    fprintf (fp, "\n Transforming branch length                   =   %2.1e", transformingBranchLength);
    fprintf (fp, "\n Healthy tip branch length                    =   %2.1e", healthyTipBranchLength);

	if (doUserGenome == YES)
		fprintf (fp, "\n\nUser genome file                              =   %s", userGenomeFile);

    if (doSimulateData == YES)
        {
        fprintf (fp, "\n\nMutation models");
        fprintf (fp, "\n Character alphabet                           =   ");
        if (alphabet == BINARY)
            fprintf (fp, "0/1");
        else if (alphabet == DNA)
            fprintf (fp, "DNA");
		
		if (doGeneticSignatures == YES)
			{
			fprintf (fp, "\n Number of trinucleotide genetic signatures   =   %d", numUserSignatures);
			for (i=0; i<numUserSignatures; i++)
				fprintf (fp, "\n  signature %02d weight                         =   %4.2f", signatureID[i], signatureWeight[i]);
			}
        else if (propAltModelSites == 0)
            fprintf (fp, "\n All sites evolving under the ISM diploid");
			
        else if (propAltModelSites > 0)
            {
			fprintf (fp, "\n Proportion of ISM diploid sites              =   %-3.2f", 1.0 - propAltModelSites);
			fprintf (fp, "\n Proportion of non-ISM diploid sites          =   %-3.2f", propAltModelSites);
			fprintf (fp, "\n Non-ISM diploid mutation model               =   ");
            if (altModel == 0) /* ISM haploid model */
                fprintf (fp, "ISM haploid");
            else if (altModel == 1) /* Markov (Mk2) model */
                fprintf (fp, "Markov 0/1");
             else if (altModel == 2) /* finite DNA model */
                {
                fprintf (fp, "Markov finite DNA: ");
				if (doJC == YES)
					fprintf (fp, "JC");
				else if (doHKY == YES)
					fprintf (fp, "HKY");
				else if (doGTR == YES)
					fprintf (fp, "GTR");
				else if (doGTnR == YES)
					fprintf (fp, "GTnR");
				}
			
			if (rateVarAmongSites == YES)
				fprintf (fp, "\n Rate variation among sites (alpha)           =   %-3.2f", alphaSites);
			else
				fprintf (fp, "\n Equal rates among sites");
		   }

        if (alphabet == DNA && propAltModelSites > 0 && altModel == 2 && doGeneticSignatures == NO)
            {
            fprintf (fp, "\n Base frequencies (ACGT)                      =   %3.2f %3.2f %3.2f %3.2f", freq[0], freq[1], freq[2], freq[3]);
			if (doHKY == YES)
				fprintf (fp, "\n Transition/transversion ratio                =   %3.2f  (kappa = %3.2f)", titv, kappa);
			if (thereIsMij == YES)
				{
				fprintf (fp, "\n Mutation rate matrix                         =   %3.2f %3.2f %3.2f %3.2f", Mij[0][0], Mij[0][1], Mij[0][2], Mij[0][3]);
				fprintf (fp, "\n                                                  %3.2f %3.2f %3.2f %3.2f", Mij[1][0], Mij[1][1], Mij[1][2], Mij[1][3]);
				fprintf (fp, "\n                                                  %3.2f %3.2f %3.2f %3.2f", Mij[2][0], Mij[2][1], Mij[2][2], Mij[2][3]);
				fprintf (fp, "\n                                                  %3.2f %3.2f %3.2f %3.2f", Mij[3][0], Mij[3][1], Mij[3][2], Mij[3][3]);
				}
			}
        if (SNPrate > 0)
            fprintf (fp, "\n SNPrate                                      =   %2.1e", SNPrate);
			
		fprintf (fp, "\n\nSNV rates");
        if (doSimulateFixedNumMutations == YES)
            fprintf (fp, "\n Fixed number of mutations                    =   %-3d", numFixedMutations);
        else
            {
            fprintf (fp, "\n Mutation rate                                =   %2.1e", mutationRate);
            fprintf (fp, "\n Population mutation parameter (4Nu2L)        =   %3.2f", theta);
            if (noisy > 1)
				fprintf (fp, "\n Exp num mutation events [Ne=cte,ISMh]        =   %3.2f", expNumMU);
            }
        fprintf (fp, "\n Mean number of mutation events               =   %3.2f", meanNumMU);
        if (numDataSets > 1)
            {
            if (doSimulateFixedNumMutations == NO)
				{
                if (noisy > 1)
					fprintf (fp, "\n Exp variance num mut events [Ne=cte,ISM]     =   %3.2f", expVarNumMU);
				fprintf (fp, "\n Variance number of mutation events           =   %3.2f", varNumMU);
				}
             }
        fprintf (fp, "\n Mean number of SNVs                          =   %3.2f", meanNumSNVs);
        if (numDataSets > 1)
            {
			fprintf (fp, "\n Variance number of SNVs                      =   %3.2f", varNumSNVs);
			fprintf (fp, "\n Proportion of reps without SNV sites         =   %3.2f", zeroSNVs/(double)numDataSets);
			}

		if (deletionRate > 0)
			{
			fprintf (fp, "\n\nDeletion rates");
			fprintf (fp, "\n Deletion rate                                =   %2.1e", deletionRate);
			fprintf (fp, "\n Mean number of deletion events               =   %3.2f", meanNumDEL);
			if (numDataSets > 1)
				fprintf (fp, "\n Variance number of deletion events           =   %3.2f", varNumDEL);
			}

		if (CNLOHrate > 0)
			{
			fprintf (fp, "\n\nCopy-neutral LOH rates");
			fprintf (fp, "\n CN-LOH rate                                  =   %2.1e", CNLOHrate);
			fprintf (fp, "\n Mean number of CN-LOH events                 =   %3.2f", meanNumCNLOH);
			if (numDataSets > 1)
				fprintf (fp, "\n Variance number of CN-LOH events            =   %3.2f", varNumCNLOH);
			}

		}
 
		if (genotypingError > 0 || doSimulateReadCounts == YES)
			fprintf (fp, "\n\nNGS");
		else
			fprintf (fp, "\n\nNGS genotype errors/read counts are not being simulated");
		
		if (genotypingError > 0)
			{
			fprintf (fp, "\n Genotype error                               =   %2.1e", genotypingError);
			//fprintf (fp, "\n Exp number of FP SNVs due to genotype errors =   %3.2f", (1 - pow(1-genotypingError, numCells+1)) * numSites);
			}
		else if (doSimulateReadCounts == YES)
			{
			fprintf (fp, "\n Sequencing coverage                          =   %dX", coverage);
			if (rateVarCoverage == YES)
				fprintf (fp, "\n Coverage dispersion (alpha)                  =   %-3.2f", alphaCoverage);
			else
				fprintf (fp, "\n  [coverage follows a Poisson distribution]");
			fprintf (fp, "\n Allelic imbalance                            =   %-3.2f", allelicImbalance);
			fprintf (fp, "\n Haploid coverage reduction                   =   %-3.2f", haploidCoverageReduction);
			fprintf (fp, "\n Alellic dropout                              =   %2.1e", ADOrate);
			fprintf (fp, "\n Amplification error");
			if (simulateOnlyTwoTemplates == YES)
				fprintf (fp, " (2-template model)");
			else
				fprintf (fp, " (4-template model)");
			fprintf (fp, "\n  Mean                                        =   %2.1e", meanAmplificationError);
			fprintf (fp, "\n  Variance                                    =   %2.1e", varAmplificationError);
			fprintf (fp, "\n Sequencing error                             =   %2.1e", sequencingError);
			fprintf (fp, "\n Doublet rate per cell                        =   %-3.2f", doubletRate);

			if (thereIsEij == YES)
				{
				fprintf (fp, "\n  NGS error rate matrix                       =   %3.2f %3.2f %3.2f %3.2f", Eij[0][0], Eij[0][1], Eij[0][2], Eij[0][3]);
				fprintf (fp, "\n                                                  %3.2f %3.2f %3.2f %3.2f", Eij[1][0], Eij[1][1], Eij[1][2], Eij[1][3]);
				fprintf (fp, "\n                                                  %3.2f %3.2f %3.2f %3.2f", Eij[2][0], Eij[2][1], Eij[2][2], Eij[2][3]);
				fprintf (fp, "\n                                                  %3.2f %3.2f %3.2f %3.2f", Eij[3][0], Eij[3][1], Eij[3][2], Eij[3][3]);
				}
			
			fprintf (fp, "\n Percentage of untrue genotype calls          =   %3.2f", cumCountMLgenotypeErrors/((numCells+1.0)*numSNVs*numDataSets));
			fprintf (fp, "\n   (see documentation)");
			}
	}

/***************** PrintCommandLine **********************/
/*	Prints the settings used in the simulations */
static void PrintCommandLine (FILE *fp, int argc,char **argv)
    {
	int i;

	if (readingParameterFile == NO)
		fprintf (fp, "Command line arguments (some might not affect) = %s", CommandLine);
	else
		{
		fprintf (fp, "Parameter file arguments (some might not affect) = ");
		
		/* Coalescent */
		fprintf (fp, " -n%d -s%d -l%d -e%d -g%2.1e -h%d", numDataSets, numCells, numSites, N, growthRate, numPeriods);
		for (i=1; i<=numPeriods; i++)
			fprintf (fp, " %d %d %d", Nbegin[i], Nend[i], cumDuration[i]-cumDuration[i-1]);

		/* Post-coalescent */
		fprintf (fp, " -k%2.1e", transformingBranchLength);
		fprintf (fp, " -q%2.1e", healthyTipBranchLength);
		fprintf (fp, " -i%6.4f", alphaBranches);

		/* Mutation model */
		fprintf (fp, " -b%d", alphabet);
		fprintf (fp, " -u%2.1e", mutationRate);
		fprintf (fp, " -d%2.1e", deletionRate);
		fprintf (fp, " -H%2.1e", CNLOHrate);
		fprintf (fp, " -j%d", numFixedMutations);
		fprintf (fp, " -S%d", numUserSignatures);
		for (i=0; i<numUserSignatures; i++)
			fprintf (fp, " %d %4.2f", signatureID[i], signatureWeight[i]);
		fprintf (fp, " -m%d", altModel);
		fprintf (fp, " -p%6.4f", propAltModelSites);
		fprintf (fp, " -w%6.4f", nonISMRelMutRate);
		fprintf (fp, " -c%2.1e", SNPrate);
		fprintf (fp, " -f%3.2f %3.2f %3.2f %3.2f", freq[0], freq[1], freq[2], freq[3]);
		fprintf (fp, " -t%6.4f", titv);
		fprintf (fp, " -a%6.4f", alphaSites);
		fprintf (fp, " -r%3.2f %3.2f %3.2f %3.2f  %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f", Mij[0][0], Mij[0][1], Mij[0][2], Mij[0][3]
		, Mij[1][0], Mij[1][1], Mij[1][2], Mij[1][3],  Mij[2][0], Mij[2][1], Mij[2][2], Mij[2][3],  Mij[3][0], Mij[3][1], Mij[3][2], Mij[3][3]);

		/* NGS */
		fprintf (fp, " -G%2.1e", genotypingError);
		fprintf (fp, " -C%d", coverage);
		fprintf (fp, " -V%6.4f", alphaCoverage);
		fprintf (fp, " -I%2.1e", allelicImbalance);
		fprintf (fp, " -D%2.1e", ADOrate);
		fprintf (fp, " -R%2.1e", haploidCoverageReduction);
		fprintf (fp, " -A%2.1e %2.1e %d", meanAmplificationError, varAmplificationError, simulateOnlyTwoTemplates);
		fprintf (fp, " -E%2.1e", sequencingError);
		fprintf (fp, " -B%2.1e", doubletRate);
		fprintf (fp, " -X%3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f", Eij[0][0], Eij[0][1], Eij[0][2], Eij[0][3]
		, Eij[1][0], Eij[1][1], Eij[1][2], Eij[1][3],  Eij[2][0], Eij[2][1], Eij[2][2], Eij[2][3],  Eij[3][0], Eij[3][1], Eij[3][2], Eij[3][3]);

		/* Output */
		if (doPrintSNVgenotypes == YES)
			fprintf (fp, " -%d", 1);
		if (doPrintSNVhaplotypes == YES)
			fprintf (fp, " -%d", 2);
		if (doPrintFullGenotypes == YES)
			fprintf (fp, " -%d", 3);
		if (doPrintFullHaplotypes == YES)
			fprintf (fp, " -%d", 4);
		if (doPrintAncestors == YES)
			fprintf (fp, " -%d", 5);
		if (doPrintTree == YES)
			fprintf (fp, " -%d", 6);
		if (doPrintTimes == YES)
			fprintf (fp, " -%d", 7);
		if (doPrintCATG == YES)
			fprintf (fp, " -%d", 8);
		if (doPrintTrueHaplotypes == YES)
			fprintf (fp, " -%d", 9);
		if (doPrintMLhaplotypes == YES)
			fprintf (fp, " -%c", '@');
		if (doPrintSeparateReplicates == YES)
			fprintf (fp, " -%c", 'v');
		if (doPrintIUPAChaplotypes == YES)
			fprintf (fp, " -%c", 'x');
		fprintf (fp, " -o%s", resultsDir);
		if (doUserTree == YES)
			fprintf (fp, " -T%s", userTreeFile);
		if (doUserGenome == YES)
			fprintf (fp, " -U%s", userGenomeFile);

		/* Other */
		if (doSimulateData == NO)
			fprintf (fp, " -%d", 0);
		fprintf (fp, " -y%d -z%d -#%ld", noisy, numNodes, originalSeed);
		}
	}


/***************** PrintDefaults **********************/
/*	Prints the default settings */

static void PrintDefaults (FILE *fp)
    {
    int		i;
		
	/* Coalescent */
    fprintf (fp,"\n-n: number of replicates =  %d", numDataSets);
    fprintf (fp,"\n-s: sample size (# cells) =  %d", numCells);
    fprintf (fp,"\n-l: number of sites =  %d", numSites);
    fprintf (fp,"\n-e: effective population size =  %d", N);
    fprintf (fp,"\n-g: growth rate =  %2.1e", growthRate);
    fprintf (fp,"\n-h: number of demographic periods = %d", numPeriods);
 	for (i=1; i<=numPeriods; i++)
			fprintf (fp, "\n  period %d =  %d %d %d", i, Nbegin[i], Nend[i], cumDuration[i]-cumDuration[i-1]);

	/* Post-coalescent */
	fprintf (fp,"\n-k: transforming branch length =  %2.1e", transformingBranchLength);
    fprintf (fp,"\n-q: healthy tip branch length =  %2.1e", healthyTipBranchLength);
    fprintf (fp,"\n-i: shape of the gamma distribution for rate variation among lineages =  %6.4f", alphaBranches);
	
	/* Mutation model */
	fprintf (fp,"\n-b: alphabet [0:binary 1:DNA] =  %d", alphabet);
	fprintf (fp,"\n-u: mutation rate =  %2.1e", mutationRate);
	fprintf (fp,"\n-d: deletion rate =  %2.1e", deletionRate);
	fprintf (fp,"\n-H: CNLOH rate =  %2.1e", CNLOHrate);
    fprintf (fp,"\n-j: fixed number of mutations =  %d", numFixedMutations);
    fprintf (fp,"\n-S: number of trinucleotide genetic signature =  %d", numUserSignatures);
    fprintf (fp,"\n-m: alternative mutation model =  %d", altModel);
    fprintf (fp,"\n-p: proportion of alternative model sites =  %6.4f", propAltModelSites);
    fprintf (fp,"\n-w: alternative/default model relative mutation rate =  %6.4f", nonISMRelMutRate);
	fprintf (fp,"\n-c: germline SNP rate =  %2.1e", SNPrate);
	fprintf (fp,"\n-f: nucleotide base frequencies = f%3.2f %3.2f %3.2f %3.2f", freq[0], freq[1], freq[2], freq[3]);
	fprintf (fp,"\n-t: transition/transversion ratio =  %6.4f", titv);
    fprintf (fp,"\n-a: shape of the gamma distribution for rate variation among sites =  %6.4f", alphaSites);
	fprintf (fp,"\n-r: mutation matrix ACGT x ACGT = %3.2f %3.2f %3.2f %3.2f  %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f", Mij[0][0], Mij[0][1], Mij[0][2], Mij[0][3], Mij[1][0], Mij[1][1], Mij[1][2], Mij[1][3],  Mij[2][0], Mij[2][1], Mij[2][2], Mij[2][3],  Mij[3][0], Mij[3][1], Mij[3][2], Mij[3][3]);

	/* NGS */
	fprintf (fp,"\n-G: genotyping error =  %2.1e", genotypingError);
	fprintf (fp,"\n-C: sequencing coverage =  %d", coverage);
	fprintf (fp,"\n-V: sequencing coverage overdispersion =  %6.4f", alphaCoverage);
	fprintf (fp,"\n-I: allelic imbalance =  %2.1e", allelicImbalance);
	fprintf (fp,"\n-D: allelic dropout =  %2.1e", ADOrate);
	fprintf (fp,"\n-R: haploid coverage reduction =  %2.1e", haploidCoverageReduction);
	fprintf (fp,"\n-A: amplification error, var, 2-template model =  %2.1e, %2.1e, %d)", meanAmplificationError, varAmplificationError, simulateOnlyTwoTemplates);
	fprintf (fp,"\n-E: sequencing error =  %2.1e", sequencingError);
	fprintf (fp,"\n-B: doublet rate per cell =  %2.1e", doubletRate);
	fprintf (fp,"\n-X: error matrix ACGT x ACGT = %3.2f %3.2f %3.2f %3.2f  %3.2f %3.2f %3.2f %3.2f  %3.2f %3.2f %3.2f %3.2f  %3.2f %3.2f %3.2f %3.2f", Eij[0][0], Eij[0][1], Eij[0][2], Eij[0][3], Eij[1][0], Eij[1][1], Eij[1][2], Eij[1][3],  Eij[2][0], Eij[2][1], Eij[2][2], Eij[2][3],  Eij[3][0], Eij[3][1], Eij[3][2], Eij[3][3]);

    /* Output */
    fprintf (fp,"\n-1: print SNV genotypes to a file =  %d", doPrintSNVgenotypes);
    fprintf (fp,"\n-2: print SNV haplotypes to a file =  %d", doPrintSNVhaplotypes);
    fprintf (fp,"\n-3: print full genotypes to a file =  %d", doPrintFullGenotypes);
    fprintf (fp,"\n-4: print full haplotypes to a file =  %d", doPrintFullHaplotypes);
    fprintf (fp,"\n-5: print ancestral genotypes =  %d", doPrintAncestors);
    fprintf (fp,"\n-6: print trees to a file =  %d", doPrintTree);
    fprintf (fp,"\n-7: print times to a file =  %d", doPrintTimes);
    fprintf (fp,"\n-8: print read counts in CATG format =  %d", doPrintCATG);
    fprintf (fp,"\n-9: print true haplotypes to a file =  %d", doPrintTrueHaplotypes);
	fprintf (fp,"\n-@: print ML haplotypes to a file =  %d", doPrintMLhaplotypes);
	fprintf (fp,"\n-v: print replicates in individual folders =  %d", doPrintSeparateReplicates);
	fprintf (fp,"\n-x: print consensus/IUPAC haplotypes =  %d", doPrintIUPAChaplotypes);
	fprintf (fp,"\n-o: results folder name =  %s", resultsDir);
	fprintf (fp,"\n-T: user tree file name =  %s", userTreeFile);
	fprintf (fp,"\n-U: user genome file name =  %s", userGenomeFile);

	/* Other */
	fprintf (fp,"\n-0: simulate just the genealogies =  %d", doSimulateData);
	fprintf (fp,"\n-y: noisy = %d", noisy);
    fprintf (fp,"\n-z: number of nodes to allocate = %d", numNodes);
    fprintf (fp,"\n-#: seed = %ld", originalSeed);
 	}

/***************************** PrintUsage *******************************/
/* Prints a short description of program usage  */
void PrintUsage(void)
{
    fprintf (stderr, "\n\n--------------------------------------------------------------------------------------------------------\n");
    fprintf (stderr, "%s", PROGRAM_NAME_UPPERCASE);
    fprintf (stderr, "\n%s generates a coalescent tree and simulates a sample of diploid genomes from somatic cells (no recombination), together with a (healthy) cell as outgroup.", PROGRAM_NAME);
    fprintf (stderr, "\n\nUsage: %s [-n# -s# -l# -e# -g# -h# (# # #) -k# -q# -i# -b# -u# -d# -H# -j# -S# (# #) -m# -p# -w# -c# -f# # # # -t# -a# -r# # # # # # # # # # # # # # # #  -G# -C# -V# -A# # # -E# -D# -R# -X# # # # # # # # # # # # # # # # -1 -2 -3 -4 -5 -6 -7 -8 -9 -v -x -oTXT -0 -y# -z# -## -?]", PROGRAM_NAME);
	
    fprintf (stderr,"\n-n: number of replicates (e.g. -n1000)");
    fprintf (stderr,"\n-s: sample size (# cells) (e.g. -s8)");
    fprintf (stderr,"\n-l: number of sites (e.g. -l500)");
    fprintf (stderr,"\n-e: effective population size (e.g. -e1000)");
    fprintf (stderr,"\n-g: growth rate (e.g. -g1e-5)");
    fprintf (stderr,"\n-h: number of demographic periods followed by effective  population size at the beginning and at the end of the period, ");
    fprintf (stderr, "\n    and the duration of the period in generations. (e.g. -d1 100 200 30000) (e.g. -d2 100 100 40000 1000 2000 20000)");
	fprintf (stderr,"\n-k: transforming branch length (e.g. -u1e-2)");
    fprintf (stderr,"\n-q: healthy tip branch length (e.g. -q1e-6)");
    fprintf (stderr,"\n-i: shape of the gamma distribution for rate variation among lineages (e.g. -i0.5)");
	fprintf (stderr,"\n-b: alphabet [0:binary 1:DNA] (e.g. -b0)");
	fprintf (stderr,"\n-u: mutation rate (e.g. -u1e-6)");
	fprintf (stderr,"\n-d: deletion rate (e.g. -d1e-5)");
	fprintf (stderr,"\n-H: CNLOH rate (e.g. -H1e-5)");
    fprintf (stderr,"\n-j: fixed number of mutations (e.g. -j100)");
    fprintf (stderr,"\n-S: trinucleotide genetic signatures (e.g. -S2 3 0.8 13 0.2)");
    fprintf (stderr,"\n-m: alternative mutation model [0:ISMhap 1:Mk2 2:finiteDNA] [default mutation model is ISM diploid] (e.g. -m2)");
    fprintf (stderr,"\n-p: proportion of alternative model sites (e.g. -p0.1)");
    fprintf (stderr,"\n-w: alternative/default model relative mutation rate (e.g. -w1)");
	fprintf (stderr,"\n-c: germline SNP rate (e.g. -c1e-5)");
	fprintf (stderr,"\n-f: nucleotide base frequencies (e.g. -f0.4 0.3 0.2 0.1)");
	fprintf (stderr,"\n-t: transition/transversion ratio (e.g. -t2)");
    fprintf (stderr,"\n-a: shape of the gamma distribution for rate variation among sites (e.g. -a0.2)");
	fprintf (stderr,"\n-r: mutation matrix ACGT x ACGT (e.g. -r0 1 2 3 1 0 4 5 2 4 0 1 3 5 1 0)");
	fprintf (stderr,"\n-G: genotyping error [no read counts] (e.g. -G0.01)");
	fprintf (stderr,"\n-C: sequencing coverage [read counts] (e.g. -C60)");
	fprintf (stderr,"\n-V: sequencing coverage overdispersion (e.g. -V5)");
	fprintf (stderr,"\n-I: allelic imbalance (e.g. -I0.5)");
	fprintf (stderr,"\n-D: allelic dropout (e.g. -D0.1)");
	fprintf (stderr,"\n-R: ADO/deletion haploid coverage reduction (e.g. -R0.5)");
	fprintf (stderr,"\n-A: amplification error (e.g. -A0.1 0.01 0)");
	fprintf (stderr,"\n-E: sequencing error (e.g. -E0.001)");
	fprintf (stderr,"\n-B: doublet rate per cell (e.g. -B0.1)");
	fprintf (stderr,"\n-X: error matrix ACGT x ACGT (e.g. -X0 1 1 1 1 0 1 1 1 1 0 1 1 1 1 0)");

    fprintf (stderr,"\n-1: print SNV genotypes to a file (e.g. -1)");
    fprintf (stderr,"\n-2: print SNV haplotypes to a file (e.g. -2)");
    fprintf (stderr,"\n-3: print full genotypes to a file (e.g. -3)");
    fprintf (stderr,"\n-4: print full haplotypes to a file (e.g. -4)");
    fprintf (stderr,"\n-5: print ancestral genotypes (e.g. -5)");
    fprintf (stderr,"\n-6: print trees to a file (e.g. -6)");
    fprintf (stderr,"\n-7: print times to a file (e.g. -7)");
    fprintf (stderr,"\n-8: print read counts in CATG format (e.g. -8)");
    fprintf (stderr,"\n-9: print true SNV haplotypes to a file (e.g. -9)");
	fprintf (stderr,"\n-@: print ML haplotypes to a file (e.g. -@)");
	fprintf (stderr,"\n-v: print replicates in individual folders (e.g. -v)");
	fprintf (stderr,"\n-x: print consensus/IUPAC haplotypes (e.g. -x)");
	fprintf (stderr,"\n-o: results folder name (e.g. -oresultsFolder)");
	fprintf (stderr,"\n-T: user-tree file (e.g. -Susertree)");
	fprintf (stderr,"\n-0: simulate just the genealogies (e.g. -0)");
	fprintf (stderr,"\n-y: noisy (e.g. -y1)");
    fprintf (stderr,"\n      = 0: does not print anything");
    fprintf (stderr,"\n      = 1:  + simulation summary");
    fprintf (stderr,"\n      = 2:  + replicate information");
    fprintf (stderr,"\n      = 3: + calculation status and event information");
    fprintf (stderr,"\n-z: number of nodes to allocate (e.g. -z5000)");
    fprintf (stderr,"\n-#: seed (e.g. -#37864287)");
    fprintf (stderr,"\n-?: Print help");
	
    fprintf (stderr,"\n\nDefaults: ");
    PrintDefaults (stderr);
	fprintf (stderr, "\n--------------------------------------------------------------------------------------------------------\n");
	fprintf (stderr,"\n...now exiting to system...\n\n");
    exit(-1);
}




/***************** Index ***************/
/* Returns index for a given node */

int Index (TreeNode *p)
    {
    return (p == NULL) ? -1 : p->index+1;
    }


/***************** Label ***************/
/* Returns label for a given node */

int Label (TreeNode *p)
    {
    return (p->anc == NULL && p->left == NULL && p->right == NULL) ? -1 : p->label+1;
    }


/************************************* ChooseTrinucMutation **********************************************/
/* Chooses uniformy a random trinucleotide mutation (out of 96) for a given signature 3D matrix with 96 types */

int ChooseTrinucMutation (double ***prob, long int *seed)
	{
	int		i, j = 0, k = 0;
	double	ran, cumProb;
	
	cumProb = 0;
	ran = RandomUniform(seed);
	for (i=0; i<4; i++)
		for (j=0; j<6; j++)
			for (k=0; k<4; k++)
				{
				cumProb += prob[i][j][k];
				if (ran <= cumProb)
					break;
				}
	return trimut(i,j,k);
	}



/************************************* ChooseUniformState **********************************************/
/* Chooses uniformy a random state according to a vector of state probabilities */

int ChooseUniformState (double *prob, long int *seed)
	{
	int			chosenState;
	double		ran, cumProb;
	
	chosenState = 0;
	cumProb = prob[chosenState];
	ran = RandomUniform(seed);

	while (ran > cumProb)
		cumProb += prob[++chosenState];
	
	return chosenState;
	}


/***************************** RandomUniform **********************************/
/*  It returns a random uniform variate in range 0..1. It is described in
 Park, S. K. and K. W. Miller.  1988.  Random number generators: good
 ones are hard to find.  Communications of the ACM, 31(10):1192-1201.
 */

double RandomUniform (long int *seed)
	{
    long int	lo, hi, test;
		
    hi = (*seed) / 127773;
    lo = (*seed) % 127773;
    test = 16807 * lo - 2836 * hi;
    if (test > 0)
        *seed = test;
    else
        *seed = test + 2147483647;
    return (double)(*seed) / (double)2147483647;
	}



/********************** RandomUniformTo ****************************/
/* it returns random uniform in range 0...max-1          */

static int RandomUniformTo (int max, long int *seed)
	{
    double	rd;
    rd = RandomUniform (seed);
    return (floor(rd*max));
	}



/********************* RandomPoisson ********************/
/* Generates a random number from a Poisson distibution with
 mean lambda.
 */

int RandomPoisson (double lambda, long int *seed)
	{
    int		poissonNumber;
    double	sum;
		
    sum = 0;
    poissonNumber = -1;
		
    while (sum <= 1.0)
		{
        sum += RandomExponential (lambda, seed);
        poissonNumber++;
		}
		
    return poissonNumber;
	}


/********************* RandomExponential ********************/
/* Generates a random number from a Exponential distibution with
 mean lambda.
 */

double RandomExponential (double lambda, long int *seed)
	{
		
    double 	exponentialNumber, U;
		
    do
		U = RandomUniform (seed);
	while (U == 0);
		
    exponentialNumber = -log (U) / lambda;
		
    return exponentialNumber;
	}


/********************* RandomBinomial ********************/
/* Generates a random number from a Binomial distibution with sucess probabilty p
	and n trials using the direct method (sum of Bernoulli variables).
 */

int RandomBinomial (double prob, int numTrials, long int *seed)
	{
	int i, sum;
	sum = 0;
	
	for(i=0; i<numTrials; i++)
		{
		if(RandomUniform(seed) < prob)
			sum++;
		}
	return sum;
	}

/********************* RandomNegativeBinomial ********************/
/*
 *    Random variates from the negative binomial distribution.
 *
 *  NOTES
 *
 *    x = the number of failures before the n-th success
 *
 *  REFERENCE
 *
 *    Devroye, L. (1986).
 *    Non-Uniform Random Variate Generation.
 *    New York:Springer-Verlag.  Pages 488 and 543.
 *
 *  METHOD
 *
 *    Generate lambda as gamma with shape parameter "dispersion" (aka size) and scale
 *    parameter "mean/dispersion".  Return a Poisson deviate with mean lambda.

**** NOTE: Extracted from rnbinom.c R code:
	rpois(rgamma(size, (1 - prob) / prob));
	rpois(rgamma(size, mu / size));

The negative binomial distribution with dispersion = n and prob = p has density

p(x) = Gamma(x+n)/(Gamma(n) x!) p^n (1-p)^x

for x = 0, 1, 2, ..., n > 0 and 0 < p <= 1.

A negative binomial distribution can arise as a mixture of Poisson distributions with mean distributed as a Γ (pgamma) distribution with scale parameter (1 - prob)/prob and shape parameter dispersion. In this model prob = scale/(1+scale), and the mean is dispersion * (1 - prob)/prob. The variance in this parametrization is n (1-p)/p^2.

The alternative parameterization, often used in ecology, and the one used here, is by the mean mu, and the dispersion parameter, where prob = dispersion/(dispersion+mu). The variance is mu + mu^2/dispersion in this parametrization.
*/

int RandomNegativeBinomial (double mean, double dispersion, long int *seed)
	{
	int		poissonRand;
	double	gammaRand;
	
	/* the RandomGamma function here has mean 1, so we need to scale it ourselves */
	gammaRand = mean/dispersion * RandomGamma (dispersion, seed);
	poissonRand = RandomPoisson (gammaRand, seed);
		
	return poissonRand;
	}


/**************************** RandomBeta *************************/
/*	Generates a beta random number given mean and variance

An algorithm for generating beta variates B(α,β) is to generate X/(X + Y), where X is a gamma variate with parameters (α, 1) and Y is an independent gamma variate with parameters (β, 1).[53]
From Numerical Recipes 3rd Edition: The Art of Scientific Computing. 2007. Press et al. ISBN-13: 978-0521880688


*/
double RandomBeta (double mean, double var, long int *seed)
	{
	double shape1, shape2, gamma1, gamma2, randBeta;
	
	/* assuming variance < mean (1-mean) */
	shape1 = mean * ((mean*(1.0-mean)/var) - 1.0);
	shape2 = (1.0-mean) * ((mean*(1.0-mean)/var) - 1.0);
	
	gamma1 = RandomGamma(shape1, seed);
	gamma2 = RandomGamma(shape2, seed);
	randBeta = gamma1 / (gamma1 + gamma2);
	
	return randBeta;
	}



/* Gamma functions are taken from Ziheng Yang's PAML package.
 http://abacus.gene.ucl.ac.uk/*/


/**************************** RandomGamma *************************/
/*	Generates a gamma number using routines in Ziheng's
	Yang tools.h in PAML
 
	Random standard gamma (Mean=Var=s,  with shape par=s, scale par=1)
	r^(s-1)*exp(-r)
 
	J. Dagpunar (1988) Principles of random variate generation,
	Clarendon Press, Oxford
 
	Calling rndgamma1() if s<1 or rndgamma2() if s>1 or exponential if s=1
 */

double	RandomGamma (double shape, long int *seed)
{
    double gammaNumber = 0;
	
    if (shape <= 0)
        fprintf (stderr, "ERROR: problems with gamma variable generation, shape < 0");
    else if (shape < 1)
        gammaNumber = RandomGamma1 (shape, seed);
    else if (shape > 1)
        gammaNumber = RandomGamma2 (shape, seed);
    else
        gammaNumber = -log (RandomUniform(seed));
	
    return (gammaNumber);
	
	
}

/*************** RandomGamma1 ***************/
double RandomGamma1 (double s, long int *seed)
{
    /* Random standard gamma for s<1
     switching method
     */
    double			r, x=0.0, small=1e-37, w;
    static double   a, p, uf, ss=10.0, d;
	
    if (s!=ss)
    {
        a  = 1.0-s;
        p  = a/(a+s*exp(-a));
        uf = p*pow(small/a,s);
        d  = a*log(a);
        ss = s;
    }
    for (;;)
    {
        r = RandomUniform(seed);
        if (r > p)
			{
            x = a-log((1.0-r)/(1.0-p));
            w=a*log(x)-d;  /* this was with comma in line above before 270917*/
			}
        else if (r>uf)
            {
            x = a*pow(r/p,1/s);
			w=x; /* this was with comma in line above before 270917*/
			}
        else
            return (0.0);
        r = RandomUniform(seed);
        if (1.0-r <= w && r > 0.0)
            if (r*(w+1.0) >= 1.0 || -log(r) <= w)
                continue;
        break;
    }
    return (x);
}


/*************** RandomGamma2 ***************/
double RandomGamma2 (double s, long int *seed)
{
    /* Random standard gamma for s>1
     Best's (1978) t distribution method
     */
    double			r ,d, f, g, x;
    static double	b, h, ss=0;
	
    if (s!=ss)
    {
        b  = s-1.0;
        h  = sqrt(3.0*s-0.75);
        ss = s;
    }
    for (;;)
    {
        r = RandomUniform(seed);
        g = r-r*r;
        f = (r-0.5)*h/sqrt(g);
        x = b+f;
        if (x <= 0.0)
            continue;
        r = RandomUniform(seed);
        d = 64*r*r*g*g*g;
        if (d*x < x-2.0*f*f || log(d) < 2*(b*log(x/b)-f))
            break;
    }
    return (x);
}



/************************ ReadUntil **************************/
/* Utility to skip comments between brackets. Borrowed from Andy
 Rambaut and Nick Grassly */

void ReadUntil(FILE *fv, char stopChar, char *what)
{
    char ch;
	
    ch=fgetc(fv);
    while (!feof(fv) && ch!=stopChar)
        ch=fgetc(fv);
	
    if (feof(fv) || ch!=stopChar) {
        fprintf (stderr, "%s missing", what);
        exit(0);
    }
}



/************************ Unlink_callback **************************/
/* Call back to remove the passed path .
From http://stackoverflow.com/questions/5467725/how-to-delete-a-directory-and-its-contents-in-posix-c
*/
int Unlink_callback(const char *fpath, const struct stat *sb, int typeflag, struct FTW *ftwbuf)
	{
    int rv = remove(fpath);

    if (rv)
        perror(fpath);

    return rv;
	}


/************************ RemoveDir **************************/
/* Removes directory and its contents
From http://stackoverflow.com/questions/5467725/how-to-delete-a-directory-and-its-contents-in-posix-c
*/
int RemoveDir(char *path)
	{
    return nftw(path, Unlink_callback, 64, FTW_DEPTH | FTW_PHYS);
	}


/************************ CheckMatrixSymmetry **************************/
/* Checks whether a given matrix is symmetric */

int CheckMatrixSymmetry(double matrix[4][4])
	{
	int i,j;
	
	for(i=0; i<4; i++)
		for(j=0; j<4; j++)
			if(matrix[i][j] != matrix[j][i])
				return NO;
	return YES;
	}


/******************** ReadParametersFromCommandLine **************************/
/*
 USED IN ORDER
	n s l e g h k q i b u d H j S m p w c f t a r G C V A E D I R B M 1 2 3 4 5 6 7 8 9 @ v x o T U 0 y z #
 
 USED
	a b c d e f g h i j k l m n o p q r s t u v w x y z
	A B C D E   G H I       M         R S T U V   X
	0 1 2 3 4 5 6 7 8 9 # @
*/
static void ReadParametersFromCommandLine (int argc,char **argv)
{
    int		i, j;
    char	flag, ch, *cat;
    double  sumPi, sum;
    float	argument;
	double  argdouble;
	
	if (argc > 0)
		{
		cat = (char *) calloc(300, sizeof(char));
		CommandLine = (char *) calloc(MAX_LINE, sizeof(char));

		for (i=0; i<argc; i++)
			{
			sprintf (cat, " %s", argv[i]);
			strcat (CommandLine, cat);
			}
		free(cat);
		}
	
    for (i=1; i<argc; i++)
        {
        argv[i]++;
        flag =* argv[i];
        argv[i]++;
        argument = -9999;
			
        switch (flag)
			{
			case 'n':
                argument = atof(argv[i]);
                numDataSets = (int) argument;
                if (numDataSets <1)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad number of replicates (%d)\n\n", numDataSets);
                    PrintUsage();
                    }
                break;
			case 's':
                argument = atof(argv[i]);
                numCells = (int) argument;
                if (numCells < 1)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad sample size (%d)\n\n", numCells);
                    PrintUsage();
                    }
                break;
            case 'l':
                argument = atof(argv[i]);
                numSites = (int) argument;
                if (numSites<1)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad sequence length (%d)\n\n", numSites);
                    PrintUsage();
                    }
                break;
            case 'e':
                argument = atof(argv[i]);
                N = (int) argument;
                if (N < 1)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad effective population size (%d)\n\n", N);
                    PrintUsage();
                    }
                break;
             case 'g':
                growthRate = atof(argv[i]);
                if (growthRate < 0)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad growth rate (%f)\n\n", growthRate);
                    PrintUsage();
                    }
                if (growthRate != 0)
                    {
                    doExponential = YES;
                    if (doDemographics == YES)
                        {
                        fprintf (stderr, "PARAMETER ERROR: Cannot have both exponential (-g) and  demographics (-h)\n\n");
                        exit (-1);
                        }
                    }
                break;
             case 'h':
				argument = atof(argv[i]);
                numPeriods = (int) argument;
				if (numPeriods > 0)
					doDemographics = YES;
                if (doDemographics == YES && doExponential == YES)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Cannot have both demographics periods (-h) and other demographics (-d)\n\n");
                    exit (-1);
                    }
                if (numPeriods <= 0)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad number of periods (%d)\n\n", numPeriods);
                    PrintUsage();
                    }
                Nbegin = 	(int *) calloc(numPeriods+1, sizeof(int));
                Nend =		(int *) calloc(numPeriods+1, sizeof(int));
                cumDuration =	(int *) calloc(numPeriods+1, sizeof(int));
                periodGrowth =	(double *) calloc(numPeriods+1, sizeof(double));
                if (Nbegin == NULL || Nend == NULL || cumDuration == NULL)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Could not allocate demographic vectors (%lu)\n", numPeriods * sizeof(int));
                    exit (-1);
                    }
                for (j=1; j<=numPeriods; j++)
                    {
                    argument = atof(argv[++i]);
                    Nbegin[j] = (int) argument;
                    argument = atof(argv[++i]);
                    Nend[j] = (int) argument;
                    argument = atof(argv[++i]);
                    cumDuration[j] = (int) argument + cumDuration[j-1];
                    }
                break;
            case 'k':
                transformingBranchLength = atof(argv[i]);
              if (transformingBranchLength < 0)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad transforming branch length (%f)\n\n", transformingBranchLength);
                    PrintUsage();
                    }
                break;
			case 'q':
                healthyTipBranchLength = atof(argv[i]);
                if (healthyTipBranchLength < 0)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad healthy tip branch length (%f)\n\n", healthyTipBranchLength);
                    PrintUsage();
                    }
                break;
             case 'i':
                alphaBranches = atof(argv[i]);
                    if (alphaBranches <= 0)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad gamma alphaBranches shape (%f)\n\n", alphaBranches);
                    PrintUsage();
                    }
                rateVarAmongLineages = YES;
                break;
			case 'b':
                argument = atof(argv[i]);
                alphabet = (int) argument;
                   if (alphabet < 0 || alphabet > 1)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad alphabet (%d)\n\n", alphabet);
                    PrintUsage();
                    }
               break;
			case 'u':
                mutationRate = atof(argv[i]);
                if (mutationRate < 0)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad mutation rate (%f)\n\n", mutationRate);
                    PrintUsage();
                    }
                break;
            case 'd':
                deletionRate = atof(argv[i]);
                if (deletionRate < 0)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad deletion rate (%f)\n\n", deletionRate);
                    PrintUsage();
                    }
                break;
            case 'H':
                CNLOHrate = atof(argv[i]);
                if (CNLOHrate < 0)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad CNLOH rate (%f)\n\n", CNLOHrate);
                    PrintUsage();
                    }
                break;
			case 'j':
                argument = atof(argv[i]);
                numFixedMutations = (int) argument;
                if (numFixedMutations <1)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad number of fixed mutations (%d)\n\n", numFixedMutations);
                    PrintUsage();
                    }
                doSimulateFixedNumMutations = YES;
                if (propAltModelSites > 0)
                    {
                    fprintf (stderr, "PARAMETER ERROR: You cannot specify a fixed number of mutations if there is any non-ISM site. Set the proportion of non-ISM sites to zero\n\n");
					PrintUsage();
                    }
                break;
			case 'S':
                argument = atof(argv[i]);
                numUserSignatures = (int) argument;
                if (numUserSignatures < 1 || numUserSignatures > 30)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad  number of genetic signatures (%d), it has to be a single number between 1 and 30\n\n", numUserSignatures);
                    PrintUsage();
                    }
                doGeneticSignatures = YES;

				if (alphabet == BINARY)
					{
					fprintf (stderr, "PARAMETER ERROR: Trinucleotide signatures cannot be simulated with the the binary alphabet");
					PrintUsage();
					}

				signatureWeight = (double*) calloc (numUserSignatures, sizeof(double));
				if (!signatureWeight)
					{
					fprintf (stderr, "Could not allocate the signatureWeight array\n");
					exit (-1);
					}

				signatureID = (int*) calloc (numUserSignatures, sizeof(int));
				if (!signatureID)
					{
					fprintf (stderr, "Could not allocate the signatureID array\n");
					exit (-1);
					}
				sum = 0;
				for (j=0; j<numUserSignatures; j++)
					{
					argument = atof(argv[++i]);
					signatureID[j] = (int) argument;
					argument = atof(argv[++i]);
					signatureWeight[j] = argument;
					if (signatureWeight[j] <= 0)
						{
						fprintf (stderr, "PARAMETER ERROR: Bad signatureWeight (%f)\n\n", signatureWeight[j]);
						PrintUsage();
						}
					sum += signatureWeight[j];
					}
				if (sum != 1.0)
					for (j=0; j<numUserSignatures; j++)
						signatureWeight[j] /= sum;
				    break;
			case 'm':
                argument = atof(argv[i]);
                altModel = (int) argument;
                if (altModel < 0 || altModel > 2)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad alternative mutation model (%d)\n\n", altModel);
                    PrintUsage();
                    }
				if (alphabet == DNA && propAltModelSites > 0)
					{
					if (altModel == Mk)
						{
						fprintf (stderr, "PARAMETER ERROR: The DNA alphabet and the model (%d) specified are incompatible", (int) argument);
						PrintUsage();
						}
					}
				else if (alphabet == BINARY && propAltModelSites > 0)
					{
					if (altModel == finiteDNA)
						{
						fprintf (stderr, "PARAMETER ERROR: The binary alphabet and the model (%d) specified are incompatible", (int) argument);
						PrintUsage();
						}
					}
               break;
		   case 'p':
                propAltModelSites = atof(argv[i]);
                if (propAltModelSites < 0 || propAltModelSites > 1)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad proportion of alternative model sites (%f). It has to be between 0 and 1\n\n", propAltModelSites);
                    PrintUsage();
                    }
				if (propAltModelSites > 0 && altModel != ISMhap && doSimulateFixedNumMutations == YES)
                    {
                    fprintf (stderr, "PARAMETER ERROR: You cannot specify a proportion of non-ISM  sites bigger than zero if the number of mutations is fixed\n\n");
                    PrintUsage();
                    }
 				if (alphabet == DNA && propAltModelSites > 0)
					{
					if (altModel == Mk)
						{
						fprintf (stderr, "PARAMETER ERROR: The DNA alphabet and the alt model (%d) specified are incompatible", (int) altModel);
						PrintUsage();
						}
					}
				else if (alphabet == BINARY && propAltModelSites > 0)
					{
					if (altModel == finiteDNA)
						{
						fprintf (stderr, "PARAMETER ERROR: The binary alphabet and the alt model (%d) specified are incompatible", (int) altModel);
						PrintUsage();
						}
					}
               break;
            case 'w':
                nonISMRelMutRate = atof(argv[i]);
                if (nonISMRelMutRate < 0 || nonISMRelMutRate > 1)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad relative nonISM/ISM mutation rate (%f)\n\n", nonISMRelMutRate);
                    PrintUsage();
                    }
                break;
			case 'c':
                SNPrate = atof(argv[i]);
                if (SNPrate < 0)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad germline SNP rate (%f)\n\n", SNPrate);
                    PrintUsage();
                    }
                break;
			case 'f':
                freq[0] = atof(argv[i]);
                freq[1] = atof(argv[++i]);
                freq[2] = atof(argv[++i]);
                freq[3] = atof(argv[++i]);
					
                if ((freq[0] == freq[1]) && (freq[0] == freq[2]) && (freq[0] == freq[3])) /* TODO fix with loop */
                    equalBaseFreq = YES;
                else
                    equalBaseFreq = NO;
					
                sumPi = freq[0] + freq[1] + freq[2] + freq[3];
                if (sumPi !=1.0)
                    {
                    freq[0]/=sumPi;
                    freq[1]/=sumPi;
                    freq[2]/=sumPi;
                    freq[3]/=sumPi;
                    }
                break;
			case 't':
	        	titv = atof(argv[i]);
				if (titv < 0)
					{
					fprintf (stderr, "PARAMETER ERROR: Bad ti/tv (%f)\n\n", titv);
					PrintUsage();
					}
				else
					{
					doJC = NO;
					doHKY = YES;
					doGTR = NO;
					doGTnR = NO;
					}
				if (thereIsMij == YES)
					{
					fprintf (stderr, "PARAMETER ERROR: Cannot specify a mutation matrix (GTR model) and a ti/tv (HKY model) at the same time\n\n");
					PrintUsage();
					}
				break;
			case 'a':
                alphaSites = atof(argv[i]);
                if (alphaSites < 0)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad gamma alphaSites shape (%f)\n\n", alphaSites);
                    PrintUsage();
                    }
                rateVarAmongSites = YES;
                break;
			case 'r':
				if (doHKY == YES)
					{
					fprintf (stderr, "PARAMETER ERROR: Cannot specify a mutation matrix (GTR model) and a ti/tv (HKY model) at the same time\n\n");
					PrintUsage();
					}
					/*	AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT */
					Mij[0][0]= atof(argv[i]);
					Mij[0][1] = atof(argv[++i]);
					Mij[0][2] = atof(argv[++i]);
					Mij[0][3] = atof(argv[++i]);
					Mij[1][0] = atof(argv[++i]);
					Mij[1][1] = atof(argv[++i]);
					Mij[1][2] = atof(argv[++i]);
					Mij[1][3] = atof(argv[++i]);
					Mij[2][0] = atof(argv[++i]);
					Mij[2][1] = atof(argv[++i]);
					Mij[2][2] = atof(argv[++i]);
					Mij[2][3] = atof(argv[++i]);
					Mij[3][0] = atof(argv[++i]);
					Mij[3][1] = atof(argv[++i]);
					Mij[3][2] = atof(argv[++i]);
					Mij[3][3] = atof(argv[++i]);
					if (Mij[0][0] != 0  || Mij[1][1] != 0 || Mij[2][2] != 0 || Mij[3][3] != 0)
						{
						fprintf(stderr, "PARAMETER ERROR: Bad general rate matrix: diagonals should be 0 \n\n");
						PrintUsage();
						}
					thereIsMij = YES;
					if (CheckMatrixSymmetry (Mij) == YES)
						{
						doJC = NO;
						doHKY = NO;
						doGTR = YES;
						doGTnR = NO;
						}
					else
						{
						doJC = NO;
						doHKY = NO;
						doGTR = NO;
						doGTnR = YES;
						}
				break;
			case 'G':
				genotypingError = atof(argv[i]);
				if (genotypingError < 0 || genotypingError > 1)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad genotyping error (%f)\n\n", genotypingError);
                    PrintUsage();
                    }
  				if (genotypingError > 0 && doSimulateReadCounts == YES)
					{
					fprintf (stderr, "PARAMETER ERROR: Cannot specify a coverage larger than 0, which implies read count generation, and a genotyping error at the same time\n\n");
					PrintUsage();
					}
               break;
			case 'C':
				argument = atof(argv[i]);
				coverage = (int) argument;
				if (coverage < 1)
					{
					fprintf (stderr, "PARAMETER ERROR: Bad coverage (%d)\n\n", coverage);
					PrintUsage();
					}
				if (coverage > 0)
					doSimulateReadCounts = YES;
				if (genotypingError > 0 && doSimulateReadCounts == YES)
					{
					fprintf (stderr, "PARAMETER ERROR: Cannot specify a coverage larger than 0, which implies read count generation, and a genotyping error at the same time\n\n");
					PrintUsage();
					}
				break;
            case 'V':
               alphaCoverage = atof(argv[i]);
				if (alphaCoverage <= 0)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad coverage dispersion (%f)\n\n", alphaCoverage);
                    PrintUsage();
                    }
                rateVarCoverage = YES;
				break;
   			case 'A':
					meanAmplificationError = atof(argv[i]);
					varAmplificationError = atof(argv[++i]);
					simulateOnlyTwoTemplates = (int) atof(argv[++i]);
				if (meanAmplificationError < 0 || meanAmplificationError > 1)
					{
					fprintf(stderr, "PARAMETER ERROR: Bad mean amplification error (%f)\n\n", meanAmplificationError);
					PrintUsage();
					}
				if (varAmplificationError < 0 || varAmplificationError >= (meanAmplificationError * (1.0 - meanAmplificationError)))
					{
					fprintf(stderr, "PARAMETER ERROR: Bad variance amplification error (%f); it has to be < mean*(1-mean)\n\n", meanAmplificationError);
					PrintUsage();
					}
				if (simulateOnlyTwoTemplates < 0 || simulateOnlyTwoTemplates > 1)
					{
					fprintf(stderr, "PARAMETER ERROR: Bad simulateOnlyTwoTemplates error (%d); it has to be 0 (assume 2 templates) or 1 (assume 4 templates)", simulateOnlyTwoTemplates);
					PrintUsage();
					}
				break;
            case 'E':
                sequencingError = atof(argv[i]);
				if (sequencingError < 0 || sequencingError > 1)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad sequencing error (%f)\n\n", sequencingError);
                    PrintUsage();
                    }
                break;
			case 'D':
               ADOrate = atof(argv[i]);
				if (ADOrate < 0 ||ADOrate > 1)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad allelic dropout (%f)\n\n", ADOrate);
                    PrintUsage();
                    }
                break;
			case 'I':
				allelicImbalance = atof(argv[i]);
				if (allelicImbalance < 0 || allelicImbalance > 1)
					{
					fprintf (stderr, "PARAMETER ERROR: Bad allelicImbalance (%f)\n\n", haploidCoverageReduction);
					PrintUsage();
					}
				break;
			case 'R':
				haploidCoverageReduction = atof(argv[i]);
				if (haploidCoverageReduction < 0 || haploidCoverageReduction > 1)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad haploid coverage reduction (%f)\n\n", haploidCoverageReduction);
                    PrintUsage();
                    }
                break;
			case 'B':
				doubletRate = atof(argv[i]);
				if (doubletRate < 0 ||  doubletRate > 1)
					{
					fprintf (stderr, "PARAMETER ERROR: Bad doublet rate (%f)\n\n", doubletRate);
					PrintUsage();
					}
				break;
			case 'X':
					/*	AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT */
					Eij[0][0]= atof(argv[i]);
					Eij[0][1] = atof(argv[++i]);
					Eij[0][2] = atof(argv[++i]);
					Eij[0][3] = atof(argv[++i]);
					Eij[1][0] = atof(argv[++i]);
					Eij[1][1] = atof(argv[++i]);
					Eij[1][2] = atof(argv[++i]);
					Eij[1][3] = atof(argv[++i]);
					Eij[2][0] = atof(argv[++i]);
					Eij[2][1] = atof(argv[++i]);
					Eij[2][2] = atof(argv[++i]);
					Eij[2][3] = atof(argv[++i]);
					Eij[3][0] = atof(argv[++i]);
					Eij[3][1] = atof(argv[++i]);
					Eij[3][2] = atof(argv[++i]);
					Eij[3][3] = atof(argv[++i]);
					if (Eij[0][0] != 0  || Eij[1][1] != 0 || Eij[2][2] != 0 || Eij[3][3] != 0)
						{
						fprintf(stderr, "PARAMETER ERROR: Bad error matrix: diagonals should be 0 \n\n");
						PrintUsage();
						}
				thereIsEij = YES;
				break;
            case '1':
                doPrintSNVgenotypes = YES;
                break;
            case '2':
                doPrintSNVhaplotypes = YES;
                break;
            case '3':
                doPrintFullGenotypes = YES;
                break;
            case '4':
                doPrintFullHaplotypes = YES;
                break;
            case '5':
                doPrintAncestors = YES;
                break;
            case '6':
                doPrintTree = YES;
                break;
            case '7':
                doPrintTimes = YES;
                break;
            case '8':
				doPrintCATG = YES;
   				if (doSimulateReadCounts == NO)
					{
					fprintf (stderr, "PARAMETER ERROR: Cannot print the CATG format when coverage is <= 0\n\n");
					PrintUsage();
					}
               break;
            case '9':
                doPrintTrueHaplotypes = YES;
                break;
			case '@':
				doPrintMLhaplotypes = YES;
				if (doSimulateReadCounts == NO)
					{
					fprintf (stderr, "PARAMETER ERROR: Cannot print ML haplotypes format when coverage is <= 0\n\n");
					PrintUsage();
					}
				break;
 			case 'v':
                doPrintSeparateReplicates = YES;
                break;
			case 'x':
                doPrintIUPAChaplotypes = YES;
                break;
 			case 'o':
				ch = *argv[i];
				if(!isspace(ch))
					strcpy(resultsDir,argv[i]);
			break;
 			case 'T':
				ch = *argv[i];
				if(!isspace(ch))
					strcpy(userTreeFile,argv[i]);
				doUserTree = YES;
			break;
 			case 'U':
				ch = *argv[i];
				if(!isspace(ch))
					strcpy(userGenomeFile,argv[i]);
				doUserGenome = YES;
			break;
			case '0':
                doSimulateData = NO;
                break;
            case 'y':
                argument = atof(argv[i]);
                noisy = (int) argument;
                if (noisy < 0)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad noisy value (%d)\n\n", noisy);
                    PrintUsage();
                    }
                break;
            case 'z':
                argument = atof(argv[i]);
                numNodes = (int) argument;
                if (numNodes<1)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad number of nodes (%d)\n\n", numNodes);
                    PrintUsage();
                    }
                break;
            case '#':
                argdouble = atof(argv[i]);
                userSeed = (long int) argdouble;
                if (userSeed < 0)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad user seed for random number (%ld)\n\n", userSeed);
                    PrintUsage();
                    }
                break;
			case '?':
                PrintUsage();
                break;
            default :
                fprintf (stderr, "PARAMETER ERROR: Incorrect parameter: %c\n\n", flag);
                PrintUsage();
                break;
        }
    }
}



/***************************** ReadParametersFromFile *******************************/
/* Reads parameter values from the parameter file  */
/*
 USED IN ORDER
	n s l e g h k q i b u d H j S m p w c f t a r G C V A E D I R B M 1 2 3 4 5 6 7 8 9 @ v x o T U 0 y z #
 
 USED
	a b c d e f g h i j k l m n o p q r s t u v w x y z
	A B C D E   G H I       M         R S T U V   X
	0 1 2 3 4 5 6 7 8 9 # @
*/

void ReadParametersFromFile ()
{
    int		j;
    char 	ch;
    double	sumPi, sum;
    float	argument;
	double	argdouble;
	
    if(feof(stdin))
        {
        fprintf (stderr, "FILE ERROR: Unable to read parameters from stdin\n");
        exit(-1);
        }
	
    ch=fgetc(stdin);
    while(isspace(ch))
        ch=fgetc(stdin);
    while(ch=='[')
        {
        ReadUntil(stdin, ']', "closing bracket");
        ch=fgetc(stdin);
        while(isspace(ch))
            ch=fgetc(stdin);
        }

    while(!feof(stdin))
        {
        argument = 0;

        switch (ch)
            {
            case 'n':
                if (fscanf(stdin, "%f", &argument) !=1 || argument < 1)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad number of replicates (%d)\n\n", (int) argument);
                    PrintUsage();
                    }
                numDataSets = (int) argument;
                break;
            case 's':
                if (fscanf(stdin, "%f", &argument) !=1 || argument < 1)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad sample size (%d)\n\n", (int) argument);
                    PrintUsage();
                    }
                numCells = (int) argument;
                break;
            case 'l':
                if (fscanf(stdin, "%f", &argument) !=1 || argument < 1)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad sequence length (%d)\n\n", (int) argument);
                    PrintUsage();
                    }
                numSites = (int) argument;
                break;
            case 'e':
                if (fscanf(stdin, "%f", &argument) !=1 || argument < 1)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad effective population size (%d)\n\n", (int) argument);
                    PrintUsage();
                    }
                N = (int) argument;
                break;
             case 'g':
                if (fscanf(stdin, "%lf", &growthRate) !=1)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad exponential growth rate (%f)\n\n", growthRate);
                    PrintUsage();
                    }
                if (growthRate != 0)
                    {
                    doExponential = YES;
                    if (doDemographics == YES)
                        {
                        fprintf (stderr, "PARAMETER ERROR: Cannot have both exponential (-b) and other demographics(-d)\n\n");
                        exit (-1);
                        }
                    }
                break;
             case 'h':
                if (fscanf(stdin, "%f", &argument) !=1 || argument < 1)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad number of periods (%d)\n\n", (int) argument);
                    PrintUsage();
                    }
				numPeriods = (int) argument;
				if (numPeriods > 0)
					doDemographics = YES;
                if (doDemographics == YES && doExponential == YES)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Cannot have both demographics periods (-d) and other demographics (-b)\n\n");
                    exit (-1);
                    }
				Nbegin =	(int *) calloc(numPeriods+1, sizeof(int));
                Nend = 		(int *) calloc(numPeriods+1, sizeof(int));
                cumDuration =	(int *) calloc(numPeriods+1, sizeof(int));
                periodGrowth =	(double *) calloc(numPeriods+1, sizeof(double));
                if (Nbegin == NULL || Nend == NULL || cumDuration == NULL)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Could not allocate demographic vectors (%lu)\n", numPeriods * sizeof(int));
                    exit (-1);
                    }
                for (j=1; j<=numPeriods; j++)
                    {
                    fscanf(stdin, "%f", &argument);
                    Nbegin[j] = (int) argument;
                    fscanf(stdin, "%f", &argument);
                    Nend[j] = (int) argument;
                    fscanf(stdin, "%f", &argument);
                    cumDuration[j] = (int) argument + cumDuration[j-1];
                    }
                break;
            case 'k':
                if (fscanf(stdin, "%lf", &transformingBranchLength)!=1 || transformingBranchLength < 0)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad transforming branch length (%f)\n\n", transformingBranchLength);
                    PrintUsage();
                    }
                 break;
            case 'q':
                if (fscanf(stdin, "%lf", &healthyTipBranchLength)!=1 || healthyTipBranchLength < 0)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad healthy tip branch length (%f)\n\n", healthyTipBranchLength);
                    PrintUsage();
                    }
                break;
            case 'i':
                    if (fscanf(stdin, "%lf", &alphaBranches)!=1 || alphaBranches <= 0)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad alphaBranches shape (%f)\n\n", alphaBranches);
                    PrintUsage();
                    }
                rateVarAmongLineages = YES;
               break;
            case 'b':
                if (fscanf(stdin, "%f", &argument) !=1 || argument < 0  || argument > 1)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad alphabet (%d)\n\n", (int) argument);
                    PrintUsage();
                    }
                alphabet = (int) argument;
                break;
			case 'u':
                if (fscanf(stdin, "%lf", &mutationRate)!=1 || mutationRate < 0 || mutationRate > 1)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad mutation rate (%f)\n\n", mutationRate);
                    PrintUsage();
                    }
                break;
 			case 'd':
                if (fscanf(stdin, "%lf", &deletionRate)!=1 || deletionRate < 0 || deletionRate > 1)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad deletion rate (%f)\n\n", deletionRate);
                    PrintUsage();
                    }
                break;
 			case 'H':
                if (fscanf(stdin, "%lf", &CNLOHrate)!=1 || CNLOHrate < 0 || CNLOHrate > 1)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad CNLOH rate (%f)\n\n", CNLOHrate);
                    PrintUsage();
                    }
                break;
           case 'j':
                if (fscanf(stdin, "%f", &argument) !=1 || argument < 1)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad number of mutations (%d)\n\n", (int) argument);
                    PrintUsage();
                    }
                numFixedMutations = (int) argument;
                doSimulateFixedNumMutations = YES;
                if (propAltModelSites > 0 && altModel != ISMhap)
                    {
                    fprintf (stderr, "PARAMETER ERROR: You cannot specify a fixed number of mutations if there is any non-ISM  site. Set the proportion of non-ISM diploid sites to zero\n\n");
                    PrintUsage();
                    }
                break;
			case 'S':
				if (fscanf(stdin, "%f", &argument) !=1 || argument < 1 || argument > 30)
					{
					fprintf (stderr, "PARAMETER ERROR: Bad number of genetic signatures (%d), it has to be a single number between 1 and 30\n\n", (int) argument);
					PrintUsage();
					}
				numUserSignatures = (int) argument;
				doGeneticSignatures = YES;

				if (alphabet == BINARY)
					{
					fprintf (stderr, "PARAMETER ERROR: Trinucleotide signatures cannot be simulated with the the binary alphabet");
					PrintUsage();
					}

				signatureWeight = (double*) calloc (numUserSignatures, sizeof(double));
				if (!signatureWeight)
					{
					fprintf (stderr, "Could not allocate the signatureWeight array\n");
					exit (-1);
					}

				signatureID = (int*) calloc (numUserSignatures, sizeof(int));
				if (!signatureID)
					{
					fprintf (stderr, "Could not allocate the signatureID array\n");
					exit (-1);
					}

				sum = 0;
				for (j=0; j<numUserSignatures; j++)
					{
					fscanf(stdin, "%f", &argument);
					signatureID[j] = (int) argument;
					fscanf(stdin, "%lf", &signatureWeight[j]);
					if (signatureWeight[j] <= 0)
						{
						fprintf (stderr, "PARAMETER ERROR: Bad signatureWeight (%f)\n\n", signatureWeight[j]);
						PrintUsage();
						}
					sum += signatureWeight[j];
					}
				if (sum != 1.0)
					for (j=0; j<numUserSignatures; j++)
						signatureWeight[j] /= sum;
				break;
            case 'm':
                  if (fscanf(stdin, "%f", &argument) !=1 || argument < 0 || argument > 2)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad alternative mutation model (%d)\n\n", (int) argument);
                    PrintUsage();
                    }
                altModel = (int) argument;
				if (alphabet == DNA && propAltModelSites > 0)
					{
					if (altModel == Mk)
						{
						fprintf (stderr, "PARAMETER ERROR: The DNA alphabet and the model (%d) specified are incompatible", (int) argument);
						PrintUsage();
						}
					}
				else if (alphabet == BINARY && propAltModelSites > 0)
                    {
					if (altModel == finiteDNA)
						{
						fprintf (stderr, "PARAMETER ERROR: The binary alphabet and the model (%d) specified are incompatible", (int) argument);
						PrintUsage();
						}
                    }
				break;
            case 'p':
                if (fscanf(stdin, "%lf", &propAltModelSites) !=1 || propAltModelSites < 0 || propAltModelSites > 1)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad proportion of alternative model sites (%f)\n\n", propAltModelSites);
                    PrintUsage();
                    }
                if (propAltModelSites > 0 && doSimulateFixedNumMutations == YES && altModel != ISMhap)
                    {
                    fprintf (stderr, "PARAMETER ERROR: You cannot specify a proportion of non-ISM diploid sites bigger than zero if the number of mutations is fixed\n\n");
                    PrintUsage();
                    }
				if (alphabet == DNA && propAltModelSites > 0)
					{
					if (altModel == Mk)
						{
						fprintf (stderr, "PARAMETER ERROR: The DNA alphabet and the alt model (%d) specified are incompatible", (int) altModel);
						PrintUsage();
						}
					}
				else if (alphabet == BINARY && propAltModelSites > 0)
                    {
					if (altModel == finiteDNA)
						{
						fprintf (stderr, "PARAMETER ERROR: The binary alphabet and the alt model (%d) specified are incompatible", (int) altModel);
						PrintUsage();
						}
                    }
                break;
            case 'w':
                if (fscanf(stdin, "%lf", &nonISMRelMutRate)!=1 || nonISMRelMutRate < 0 || nonISMRelMutRate > 1)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad relative non-ISM/ISM mutation rate (%f)\n\n", nonISMRelMutRate);
                    PrintUsage();
                    }
                break;
		   case 'c':
                if (fscanf(stdin, "%lf", &SNPrate)!=1 || SNPrate < 0)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad germline SNP rate (%f)\n\n", SNPrate);
                    PrintUsage();
                    }
                break;
            case 'f':
                if (fscanf(stdin, "%lf %lf %lf %lf", &freq[0], &freq[1], &freq[2], &freq[3])!=4)
                    {
                    fprintf(stderr, "PARAMETER ERROR: Bad Base Frequencies\n\n");
                    PrintUsage();
                    }
                else if (freq[0] == freq[1] == freq[2] == freq[3])
                    equalBaseFreq = YES;
                else
                    equalBaseFreq = NO;
                sumPi = freq[0] + freq[1] + freq[2] + freq[3];
                if (sumPi != 1.0)
                    {
                    freq[0]/=sumPi;
                    freq[1]/=sumPi;
                    freq[2]/=sumPi;
                    freq[3]/=sumPi;
                    }
                break;
			case 't':
				if (fscanf(stdin, "%lf", &titv)!=1)
					{
					fprintf(stderr, "Bad ti/tv\n\n");
					PrintUsage();
					}
				else
					{
					doJC = NO;
					doHKY = YES;
					doGTR = NO;
					doGTnR = NO;
					}
				if (thereIsMij == YES)
					{
					fprintf (stderr, "PARAMETER ERROR: Cannot specify a mutation matrix (GTR model) and a ti/tv (HKY model) at the same time\n\n");
					PrintUsage();
					}
				break;
 			case 'a':
                if (fscanf(stdin, "%lf", &alphaSites)!=1 || alphaSites < 0)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad alphaSites shape (%f)\n\n", alphaSites);
                    PrintUsage();
                    }
                rateVarAmongSites = YES;
                break;
			case 'r':
					if (doHKY == YES)
						{
						fprintf (stderr, "PARAMETER ERROR: Cannot specify a mutation matrix (GTR model) and a ti/tv (HKY model) at the same time\n\n");
						PrintUsage();
						}
					if (fscanf(stdin, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf %lf",
							   &Mij[0][0], &Mij[0][1], &Mij[0][2], &Mij[0][3],
							   &Mij[1][0], &Mij[1][1], &Mij[1][2], &Mij[1][3],
							   &Mij[2][0], &Mij[2][1], &Mij[2][2], &Mij[2][3],
							   &Mij[3][0], &Mij[3][1], &Mij[3][2], &Mij[3][3])!=16)
						{
						fprintf(stderr, "PARAMETER ERROR: Bad general rate matrix (-rx x x x x x x x x x x) (AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT)\n\n");
						PrintUsage();
						}
					
					if (Mij[0][0] != 0  || Mij[1][1] != 0 || Mij[2][2] != 0 || Mij[3][3] != 0)
						{
						fprintf(stderr, "PARAMETER ERROR: Bad general rate matrix: diagonals should be 0 \n\n");
						PrintUsage();
						}
					thereIsMij = YES;
					
					if (CheckMatrixSymmetry (Mij) == YES)
						{
						doJC = NO;
						doHKY = NO;
						doGTR = YES;
						doGTnR = NO;
						}
					else
						{
						doJC = NO;
						doHKY = NO;
						doGTR = NO;
						doGTnR = YES;
						}
					break;
			case 'G':
				if (fscanf(stdin, "%lf", &genotypingError) !=1 || genotypingError < 0 || genotypingError > 1)
					{
					fprintf (stderr, "PARAMETER ERROR: Bad genotyping error (%f)\n\n", genotypingError);
					PrintUsage();
					}
  				if (genotypingError > 0 && doSimulateReadCounts == YES)
					{
					fprintf (stderr, "PARAMETER ERROR: Cannot specify a coverage larger than 0, which implies read count generation, and a genotyping error at the same time\n\n");
					PrintUsage();
					}
				break;
			case 'C':
				if (fscanf(stdin, "%f", &argument) !=1 || argument < 1)
					{
					fprintf (stderr, "PARAMETER ERROR: Bad sequencing coverage (%d)\n\n", (int) coverage);
					PrintUsage();
					}
				coverage = (int) argument;
				if (coverage > 0)
					doSimulateReadCounts = YES;
				if (genotypingError > 0 && doSimulateReadCounts == YES)
					{
					fprintf (stderr, "PARAMETER ERROR: Cannot specify a coverage larger than 0, which implies read count generation, and a genotyping error at the same time\n\n");
					PrintUsage();
					}
				break;
            case 'V':
                    if (fscanf(stdin, "%lf", &alphaCoverage)!=1 || alphaCoverage <= 0)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad coverage dispersion (%f)\n\n", alphaCoverage);
                    PrintUsage();
                    }
                rateVarCoverage = YES;
			break;
  			case 'A':
				if (fscanf(stdin, "%lf %lf %d", &meanAmplificationError, &varAmplificationError, &simulateOnlyTwoTemplates) != 3)
					{
					fprintf(stderr, "PARAMETER ERROR: Bad mean/var/model amplification error (%f ; %f ; model=%d)\n\n", meanAmplificationError, varAmplificationError, simulateOnlyTwoTemplates);
					PrintUsage();
					}
				if (meanAmplificationError < 0 || meanAmplificationError > 1)
					{
					fprintf(stderr, "PARAMETER ERROR: Bad mean amplification error (%f)\n\n", meanAmplificationError);
					PrintUsage();
					}
				if (varAmplificationError < 0 || varAmplificationError >= (meanAmplificationError * (1.0 - meanAmplificationError)))
					{
					fprintf(stderr, "PARAMETER ERROR: Bad variance amplification error (%f); it has to be < mean*(1-mean)\n\n", meanAmplificationError);
					PrintUsage();
					}
				if (simulateOnlyTwoTemplates != 0 && simulateOnlyTwoTemplates != 1)
					{
					fprintf(stderr, "PARAMETER ERROR: Bad simulateOnlyTwoTemplates error (%d); it has to be 0 (assume 4 templates) or 1 (assume 2 templates)", simulateOnlyTwoTemplates);
					PrintUsage();
					}
				break;
			case 'E':
				if (fscanf(stdin, "%lf", &sequencingError) !=1 || sequencingError < 0 || sequencingError > 1)
					{
					fprintf (stderr, "PARAMETER ERROR: Bad sequencing error (%f)\n\n", sequencingError);
					PrintUsage();
					}
				break;
			case 'D':
				if (fscanf(stdin, "%lf", &ADOrate) !=1 ||ADOrate < 0 ||ADOrate > 1)
					{
					fprintf (stderr, "PARAMETER ERROR: Bad allelic dropout rate (%f)\n\n", ADOrate);
					PrintUsage();
					}
				break;
			case 'I':
				if (fscanf(stdin, "%lf", &allelicImbalance) !=1 || allelicImbalance < 0 || allelicImbalance > 1)
					{
					fprintf (stderr, "PARAMETER ERROR: Bad allelic imbalance (%f)\n\n", allelicImbalance);
					PrintUsage();
					}
				break;
			case 'R':
				if (fscanf(stdin, "%lf", &haploidCoverageReduction) !=1 ||haploidCoverageReduction < 0 || haploidCoverageReduction > 1)
					{
					fprintf (stderr, "PARAMETER ERROR: Bad haploid coverage reduction (%f)\n\n", haploidCoverageReduction);
					PrintUsage();
					}
				break;
			case 'B':
				if (fscanf(stdin, "%lf", &doubletRate) !=1 || doubletRate < 0 || doubletRate > 1)
					{
					fprintf (stderr, "PARAMETER ERROR: Bad one allele coverage proportion (%f)\n\n", doubletRate);
					PrintUsage();
					}
				break;
			case 'X':
					if (fscanf(stdin, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf %lf",
							   &Eij[0][0], &Eij[0][1], &Eij[0][2], &Eij[0][3],
							   &Eij[1][0], &Eij[1][1], &Eij[1][2], &Eij[1][3],
							   &Eij[2][0], &Eij[2][1], &Eij[2][2], &Eij[2][3],
							   &Eij[3][0], &Eij[3][1], &Eij[3][2], &Eij[3][3])!=16)
						{
						fprintf(stderr, "PARAMETER ERROR: Bad  error matrix (-rx x x x x x x x x x x) (AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT)\n\n");
						PrintUsage();
						}
					
					if (Eij[0][0] != 0  || Eij[1][1] != 0 || Eij[2][2] != 0 || Eij[3][3] != 0)
						{
						fprintf(stderr, "PARAMETER ERROR: Bad error matrix: diagonals should be 0 \n\n");
						PrintUsage();
						}
					thereIsEij = YES;
				break;
            case '1':
                doPrintSNVgenotypes = YES;
                break;
            case '2':
                doPrintSNVhaplotypes = YES;
               break;
		   case '3':
                doPrintFullGenotypes = YES;
                break;
            case '4':
                doPrintFullHaplotypes = YES;
                break;
            case '5':
                doPrintAncestors = YES;
                break;
            case '6':
                doPrintTree = YES;
				break;
            case '7':
                doPrintTimes = YES;
                break;
            case '8':
				doPrintCATG = YES;
				if (doSimulateReadCounts == NO)
					{
					fprintf (stderr, "PARAMETER ERROR: Cannot print the CATG format when coverage is <= 0\n\n");
					PrintUsage();
					}
				break;
			case '9':
                doPrintTrueHaplotypes = YES;
				break;
			case '@':
				doPrintMLhaplotypes = YES;
				if (doSimulateReadCounts == NO)
					{
					fprintf (stderr, "PARAMETER ERROR: Cannot print the ML haplotypes format when coverage is <= 0\n\n");
					PrintUsage();
					}
				break;
			case 'v':
                doPrintSeparateReplicates = YES;
                break;
			case 'x':
                doPrintIUPAChaplotypes = YES;
                break;
			case 'o':
				ch=fgetc(stdin);
				if(!isspace(ch))
					{
					j=0;
					do
						{
						resultsDir[j]=ch;
						j++;
						ch=fgetc(stdin);
						}
					while(!isspace(ch));
					resultsDir[j]='\0';
					}
			break;
			case 'T':
				ch=fgetc(stdin);
				if(!isspace(ch))
					{
					j=0;
					do
						{
						userTreeFile[j]=ch;
						j++;
						ch=fgetc(stdin);
						}
					while(!isspace(ch));
					userTreeFile[j]='\0';
					}
				doUserTree = YES;
			break;
			case 'U':
				ch=fgetc(stdin);
				if(!isspace(ch))
					{
					j=0;
					do
						{
						userGenomeFile[j]=ch;
						j++;
						ch=fgetc(stdin);
						}
					while(!isspace(ch));
					userGenomeFile[j]='\0';
					}
				doUserGenome = YES;
			break;
			case '0':
                doSimulateData = NO;
                break;
            case 'y':
                if (fscanf(stdin, "%f", &argument) !=1 || argument < 0 || argument > 3)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad noisy value (%d)\n\n", (int) argument);
                    PrintUsage();
                    }
                noisy = (int) argument;
                break;
            case 'z':
                if (fscanf(stdin, "%f", &argument) !=1 || argument < 1)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad number of nodes (%d)\n\n", (int) argument);
                    PrintUsage();
                    }
                numNodes = (int) argument;
                break;
			case '#':
                if (fscanf(stdin, "%lf", &argdouble) !=1 || argument < 0)
                    {
                    fprintf (stderr, "PARAMETER ERROR: Bad seed (%ld)\n\n", (long int) argdouble);
                    PrintUsage();
                    }
                userSeed = (long int) argdouble;
                break;
				case '?':
                PrintUsage();
			break;
            default :
                fprintf (stderr, "PARAMETER ERROR: Incorrect parameter: %c\n\n", ch);
                PrintUsage();
                break;
			}
        ch=fgetc(stdin);
        while(isspace(ch) && !feof(stdin))
            ch=fgetc(stdin);
        while(ch=='[')
			{
            ReadUntil(stdin, ']', "closing bracket");
            ch=fgetc(stdin);
            while(isspace(ch))
                ch=fgetc(stdin);
			}
		}
	}
