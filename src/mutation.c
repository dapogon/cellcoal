/*
	CellCoal - genomes and mutation: ISM, Mk2, finite-sites, signatures, cnLOH, deletions, sparse genomes
	See cellcoal.h for the shared declarations. This file was split from cellcoal.c.
*/
#include "cellcoal.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "definitions.h"

/********************************** InitializeGenomes ***********************************/
/* Initialize all genomes with the reference states  */
void InitializeGenomes (TreeNode *p, long int *seed)
    {
    int     i, cell, anccell, site;
    double  ran;
		
    if (p != NULL)
        {
        cell = p->label;
        if (sparseMode == YES && p->isHealthyRoot == NO)
            ;	/* sparse mode: nodes are identical to the root until they are modified */
        else if (alphabet == DNA)
            {
            if (p->isHealthyRoot == YES)
                {
				if (doUserGenome == YES)
					{
					for (site=0; site<numSites; site++) /* initialize genome with the user provided sequences */
						{
						SetG(MATERNAL, cell, site, WhichNucChar(maternalUserGenome[site]));
						SetG(PATERNAL, cell, site, WhichNucChar(paternalUserGenome[site]));
						allSites[site].referenceAllele = (refUserGenome != NULL) ? WhichNucChar(refUserGenome[site]) : GetG(MATERNAL, cell, site);   // then allSites[site].referenceAllele hosts the reference genome
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
								allSites[site].referenceAllele = SetG(MATERNAL, cell, site, SetG(PATERNAL, cell, site, i));
								// then allSites[site].referenceAllele hosts the reference genome
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
                /* here the genomes are dense (sparse mode never gets here): copy the ancestor's rows */
                memcpy (data[MATERNAL][cell], data[MATERNAL][anccell], numSites * sizeof(Geno));
                memcpy (data[PATERNAL][cell], data[PATERNAL][anccell], numSites * sizeof(Geno));
                }
             }
        else /* binary alphabet */
            for (site=0; site<numSites; site++)
                SetG(MATERNAL, cell, site, SetG(PATERNAL, cell, site, 0));

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
	prob4 = (double *) Checked (calloc(4, sizeof(double)), "Could not allocate the prob4 vector\n");
		
	/* choose first trinucleotide */
	chosenTriNucleotide = ChooseUniformState(triNucFreq, seed);

	/* find bases of the selected trinucleotide */
	n1 = chosenTriNucleotide/16;
	rest = chosenTriNucleotide%16;
	n2 = rest/4;
	n3 = rest%4;
	//fprintf (stderr, "\n%2d %c%c%c ", chosenTriNucleotide, WhichNuc(n1), WhichNuc(n2), WhichNuc(n3));

	site = 0;
	SetG(MATERNAL, cell, site, SetG(PATERNAL, cell, site, allSites[site].referenceAllele = n1));
	site++;
	SetG(MATERNAL, cell, site, SetG(PATERNAL, cell, site, allSites[site].referenceAllele = n2));
	site++;
	SetG(MATERNAL, cell, site, SetG(PATERNAL, cell, site, allSites[site].referenceAllele = n3));
	
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
		SetG(MATERNAL, cell, site, SetG(PATERNAL, cell, site, allSites[site].referenceAllele = nextNucleotide));

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
	  int *newPosition = (int *)realloc(trin->position, trin->tempLength * sizeof(int));
	  if (!newPosition)
		  {
		  free (trin->position);	/* realloc leaves the original block untouched when it fails */
		  fprintf (stderr, "Could not reallocate trinucleotide position array when adding\n");
		  exit (-1);
		  }
	  trin->position = newPosition;
	  }
	trin->position[trin->numAvailablePositions++] = site;
  	}




/********************************** Observation layer (sparse mode) ***********************************/
/* In sparse mode ADO and genotyping errors are not stored in the genomes. ADO is a function of a hash of
   (genome, cell, site), compared with the ADO rate of the genotype, and genotyping errors are a sparse list of
   events per site. Only the sites that vary (or have an error) need to be examined: the other sites are identical
   to the root in all the cells. Note that the random numbers differ from those of the dense mode. */

typedef struct { int cell; signed char m, p; } ErrEv;	/* new maternal/paternal alleles (-1: unchanged) */

static double		sparseAdoFixed = 0.0, *sparseAdoCell = NULL, *sparseAdoSite = NULL;
static uint64_t	sparseAdoBase = 0;
static ErrEv		**errEv = NULL;
static int			*errCap = NULL;
static long		totalErrEv = 0;

static inline double HashUniform (uint64_t base, uint64_t a, uint64_t b)
	{
	uint64_t x = base ^ (a * 0x9E3779B97F4A7C15ULL) ^ ((b + 1) * 0xC2B2AE3D27D4EB4FULL);
	(void) SplitMix64 (&x);
	return (double)((SplitMix64 (&x) >> 11) + 0.5) * (1.0 / 9007199254740992.0);
	}

Geno ObservedAllele (int genome, int node, int site, Geno trueAllele)
	{
	double	rate, r;

	if (trueAllele == DELETION)
		return trueAllele;
	if (sparseADO && adoBlocksOn)
		{
		if (AdoBlockLost (node, genome, site))
			return ADO;
		}
	else if (sparseADO)
		{
		if (sparseAdoCell == NULL)
			rate = sparseAdoFixed;
		else
			{
			r = sparseAdoCell[node] + sparseAdoSite[site] - sparseAdoCell[node] * sparseAdoSite[site];
			rate = 1.0 - sqrt (1.0 - r);
			}
		if (HashUniform (sparseAdoBase, (uint64_t) node * 2 + genome, (uint64_t) site) < rate)
			return ADO;
		}
	if (errCount != NULL && errCount[site] > 0)
		{
		int lo = 0, hi = errCount[site] - 1;
		while (lo <= hi)
			{
			int mid = (lo + hi) / 2;
			if (errEv[site][mid].cell == node)
				{
				int ns = (genome == MATERNAL) ? errEv[site][mid].m : errEv[site][mid].p;
				return (ns >= 0) ? (Geno) ns : trueAllele;
				}
			if (errEv[site][mid].cell < node)
				lo = mid + 1;
			else
				hi = mid - 1;
			}
		}
	return trueAllele;
	}

void ResetObservationLayer (void)
	{
	long	site;

	obsActive = NO;
	sparseADO = NO;
	FreeAdoBlocks ();
	free (sparseAdoCell);
	free (sparseAdoSite);
	sparseAdoCell = sparseAdoSite = NULL;
	if (errEv != NULL)
		{
		for (site=0; site<numSites; site++)
			free (errEv[site]);
		free (errEv);
		free (errCount);
		free (errCap);
		errEv = NULL;
		errCount = errCap = NULL;
		}
	totalErrEv = 0;
	}

void SparseAllelicDropout (long int *seed)
	{
	int i, j;

	if (adoBlockLen > 0.0)
		{
		BuildAdoBlocks (seed);
		sparseADO = YES;
		return;
		}
	sparseAdoBase = ((uint64_t)(RandomUniform (seed) * 4294967296.0) << 32) | (uint64_t)(RandomUniform (seed) * 4294967296.0);
	if (doADOcell == NO && doADOsite == NO)
		sparseAdoFixed = 1.0 - sqrt (1.0 - fixedADOrate);	/* ADO rate of an allele */
	else
		{
		sparseAdoCell = (double*) calloc (numCells+1, sizeof(double));
		sparseAdoSite = (double*) calloc (numSites, sizeof(double));
		if (!sparseAdoCell || !sparseAdoSite)
			{
			fprintf (stderr, "Could not allocate the ADO cell/site rate structures\n");
			exit (-1);
			}
		for (i=0; i<numCells+1; i++)
			sparseAdoCell[i] = (doADOcell == YES) ? RandomBetaMeanVar(meanADOcell, varADOcell, seed) : 0;
		if (doADOcell == YES)
			{
			WriteADOCells (sparseAdoCell, adoReplicate);
			KeepCellADO (sparseAdoCell);
			}
		for (j=0; j<numSites; j++)
			sparseAdoSite[j] = (doADOsite == YES) ? RandomBetaMeanVar(meanADOsite, varADOsite, seed) : 0;
		}
	sparseADO = YES;
	}

/* Genotyping errors: a genotype has an error with probability G (its genotyping error), where G ~ Beta(mean, var)
   is drawn for every genotype, and then each allele has an error with probability a = 1 - sqrt(1-G). The genotypes
   with an error are therefore selected with probability E[G] = mean, with G distributed as the size-biased Beta
   (shape1 + 1, shape2) among them, and, given that the genotype has an error, both alleles are changed with
   probability a^2/G and only the maternal (or only the paternal) one with probability a(1-a)/G each. */
void SparseGenotypeError (long int *seed)
	{
	int		i, j, cell, site, k, ns;
	double	normEij[4][4], rowSum, ss, shape1, shape2, g, a, u, pBoth, pOne, g1, g2, N, idx, logq;
	int		newAllele[2];

	for (i=0; i<4; i++)
		{
		rowSum = Eij[i][0] + Eij[i][1] + Eij[i][2] + Eij[i][3];
		for (j=0; j<4; j++)
			normEij[i][j] = (rowSum > 0.0) ? Eij[i][j] / rowSum : 0.0;
		}
	ss = (varGenotypingError > 0.0) ? (meanGenotypingError * (1.0 - meanGenotypingError) / varGenotypingError) - 1.0 : 0.0;
	shape1 = meanGenotypingError * ss;
	shape2 = (1.0 - meanGenotypingError) * ss;

	errEv = (ErrEv **) calloc (numSites, sizeof(ErrEv *));
	errCount = (int *) calloc (numSites, sizeof(int));
	errCap = (int *) calloc (numSites, sizeof(int));
	if (!errEv || !errCount || !errCap)
		{
		fprintf (stderr, "Could not allocate the genotyping error structures\n");
		exit (-1);
		}

	N = (double)(numCells + 1) * (double) numSites;
	logq = log1p (-meanGenotypingError);
	idx = -1.0;
	while (1)
		{
		idx += floor (log (RandomUniform (seed)) / logq) + 1.0;	/* geometric skip to the next genotype with an error */
		if (idx >= N)
			break;
		cell = (int)(idx / numSites);
		site = (int)(idx - (double) cell * numSites);

		if (ss > 0.0)
			{
			g1 = RandomGamma (shape1 + 1.0, seed);
			g2 = RandomGamma (shape2, seed);
			g = g1 / (g1 + g2);
			}
		else
			g = meanGenotypingError;
		a = 1.0 - sqrt (1.0 - g);
		pBoth = a * a / g;
		pOne = a * (1.0 - a) / g;
		u = RandomUniform (seed);
		newAllele[0] = newAllele[1] = -1;
		for (k=0; k<2; k++)
			{
			int change = (k == 0) ? (u < pBoth + pOne) : (u < pBoth || u >= pBoth + pOne);
			Geno b;
			if (!change)
				continue;
			b = GetG ((k == 0) ? MATERNAL : PATERNAL, cell, site);	/* true allele (the layer is not active yet) */
			if (b == DELETION)
				continue;
			if (alphabet == DNA)
				ns = ChooseUniformState (normEij[b], seed);
			else
				ns = (b == 0) ? 1 : (b == 1 ? 0 : -1);
			newAllele[k] = ns;
			}
		if (newAllele[0] < 0 && newAllele[1] < 0)
			continue;
		if (errCount[site] == errCap[site])
			{
			errCap[site] = errCap[site] ? 2 * errCap[site] : 4;
			errEv[site] = (ErrEv *) realloc (errEv[site], errCap[site] * sizeof(ErrEv));
			if (!errEv[site])
				{
				fprintf (stderr, "Could not allocate the genotyping error events\n");
				exit (-1);
				}
			}
		errEv[site][errCount[site]].cell = cell;
		errEv[site][errCount[site]].m = (signed char) newAllele[0];
		errEv[site][errCount[site]].p = (signed char) newAllele[1];
		errCount[site]++;
		totalErrEv++;
		allSites[site].hasGenotypeError = YES;
		}
	}


/********************************** DecideSparseMode ***********************************/
/* Sparse genome storage is used when the dense genome matrices would need too much memory
   (more than CELLCOAL_DENSE_MB megabytes, 2048 by default), or when it is requested with
   the environment variable CELLCOAL_SPARSE=1 (CELLCOAL_SPARSE=0 forces dense storage). */
void DecideSparseMode (void)
	{
	char	*env;
	double	denseMB, limitMB = 2048.0;
	int		force = -1;

	if (doNGS == YES && doBam == 0 && doFastq == 0 && doPrintVCF == NO && doPrintCATG == NO && doPrintPILEUP == NO && doPrintBulk == NO)
		{
		fprintf (stderr, "\nERROR: -Z (no VCF) needs another read output: -8 (CATG), -9 (PILEUP) or -0 (bulk)\n");
		exit (-1);
		}
	sparseMode = NO;
	streamMode = NO;
	readsOnly = NO;
	if (doNGS == YES && doBam == 0 && doFastq == 0)
		{
		/* reads of data sets whose per-cell, per-site data would not fit in memory are streamed in blocks of sites */
		double	cellMB = (double) (numCells+1) * (double) numSites * (double) sizeof(CellSiteStr) / 1048576.0, streamLimitMB = 2048.0;
		int		streamForce = -1, supported;

		if ((env = getenv ("CELLCOAL_STREAM")) != NULL && *env != '\0')
			streamForce = (atoi (env) != 0);
		if ((env = getenv ("CELLCOAL_STREAM_MB")) != NULL && atof (env) > 0)
			streamLimitMB = atof (env);
		supported = YES;	/* -3 and -4 (ML alleles of all the sites) are kept in a temporary file, -5 and -Y read the genomes */
		readsOnly = (doPrintVCF == NO && doPrintCATG == NO && doPrintPILEUP == NO && doPrintBulk == YES && doPrintSNVgenotypes == NO && doPrintSNVhaplotypes == NO &&
					doPrintFullGenotypes == NO && doPrintFullHaplotypes == NO);
		if (streamForce == 1 || (streamForce == -1 && (cellMB > streamLimitMB || readsOnly == YES) && supported == YES))
			{
			streamMode = YES;
			sparseMode = YES;
			return;
			}
		}
	readsOnly = NO;	/* reads-only needs streamed reads */
	if ((env = getenv ("CELLCOAL_SPARSE")) != NULL && *env != '\0')
		force = (atoi (env) != 0);
	if ((env = getenv ("CELLCOAL_DENSE_MB")) != NULL && atof (env) > 0)
		limitMB = atof (env);

	denseMB = (double) ploidy * (2.0 * numCells + 1.0) * (double) numSites / 1048576.0;

	if (force == 1 || (force == -1 && denseMB > limitMB))
		sparseMode = YES;

	if (sparseMode == YES && noisy > 1)
		fprintf (stderr, "\nUsing sparse genome storage (dense genomes would need %.0f MB)\n", denseMB);
	}


/********************************** ResetSparseColumns ***********************************/
void ResetSparseColumns (void)
	{
	long	site;

	if (sparseCols == NULL)
		return;
	for (site=0; site<numSites; site++)
		{
		free (sparseCols[site]);
		sparseCols[site] = NULL;
		}
	numSparseCols = 0;
	}


/********************************** NewSparseColumn ***********************************/
/* Allocate the column of a site that starts to vary in sparse mode: all the nodes get the genotype of the root */
Geno *NewSparseColumn (int site)
	{
	Geno	*c;
	long	n = (long) numNodes;

	c = (Geno *) malloc (2 * n);
	if (!c)
		{
		fprintf (stderr, "\n\nERROR: not enough memory to store the genotypes of site %d (%ld sites vary so far, %.0f MB).\n", site+1, numSparseCols, (double) numSparseCols * 2.0 * n / 1048576.0);
		fprintf (stderr, "Try a smaller mutation rate, fewer cells or sites.\n");
		exit (-1);
		}
	memset (c, data[MATERNAL][OUTGROUP_ROOT][site], n);
	memset (c + n, data[PATERNAL][OUTGROUP_ROOT][site], n);
	sparseCols[site] = c;
	numSparseCols++;
	return c;
	}


/********************************** InheritRootGenome ***********************************/
/* Copy the genomes of the healthy root (with germline variation) to all its descendants */
void InheritRootGenome (TreeNode *p)
    {
    int cell, anccell, site;

    if (p == NULL || sparseMode == YES)
        return;	/* in sparse mode all nodes inherit the root genome implicitly */
    if (p->isHealthyRoot == NO)
        {
        cell = p->label;
        anccell = p->anc->label;
        for (site=0; site<numSites; site++)
            {
            SetG(MATERNAL, cell, site, GetG(MATERNAL, anccell, site));
            SetG(PATERNAL, cell, site, GetG(PATERNAL, anccell, site));
            }
        }
    InheritRootGenome (p->left);
    InheritRootGenome (p->right);
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
				SetG(genome, OUTGROUP_ROOT, site, 1);
			else // DNA mutation according to Mij matrix
				{
				current_state = GetG(genome, OUTGROUP_ROOT, site);
				ran = RandomUniform(seed) * cumMij[current_state][3];
                for (j=0; j<4; j++)
                    {
					if (ran <= cumMij[current_state][j])
                        {
                        SetG(genome, OUTGROUP_ROOT, site, j);
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
                SetG(genome, cell, site, GetG(genome, anccell, site));
                }
            else /* => there will be change */
                {
                if (GetG(genome, anccell, site) == 0)  /* checking all this might be excessive */
                    SetG(genome, cell, site, 1);
				else if (GetG(genome, anccell, site) == 1)
                    {
					if (allSites[site].isSNP == YES) // this is a germline mutation
						SetG(genome, cell, site, 0);
					else
						{
						fprintf (stderr, "\n\nERROR: site %d in genome %d of cell %d cannot mutate twice under the ISM model", site, genome, cell);
                    	exit(-1);
						}
					}
                else
                    {
                    fprintf (stderr, "\n\nERROR: site %d in genome %d of cell %d has an unknow state %d under the ISM model", site, genome, anccell, GetG(genome, anccell, site));
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
            ancstate = GetG(genome, anccell, site);
            cumBranchLength += p->branchLength;
				
            if ((cumBranchLength < uniform) || /* => there will be no change */
                ((doISMhaploid == NO)  && (allSites[site].numMutations > 0))  ||
                ((doISMhaploid == YES) && (genome == MATERNAL) && (allSites[site].numMutationsMaternal > 0)) ||
                ((doISMhaploid == YES) && (genome == PATERNAL) && (allSites[site].numMutationsPaternal > 0)))
                {
                SetG(genome, cell, site, ancstate);
                }
            else /* => there will be change */
                {
                ran = RandomUniform(seed) * cumMij[ancstate][3];
                for (j=0; j<4; j++)
                    {
                    if (ran <= cumMij[ancstate][j])
                        {
                        SetG(genome, cell, site, j);
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

    Diploid here means that if site 33 mutates in the female genome,
    site 33 in the male genome will not mutate, and viceversa.

    The number of mutations will be distributed as:
 
    (1) a Poisson with parameter the sum of
	branch lengths over all sites  Times are already scaled in
	2N. We will force it to get just one mutation per site
 
	(2) A used-defined number of mutations (segregating sites in the ISM diploid)
 
=== Importantly, in the current implementation we assume that the context does not evolve ==
 
 Also, note that to get the expected counts under each signature it is important to simulate
 long genomes, otherwise trinucleotide exhaustion will bias the counts

*/

/* Frees the trinucleotide position vectors. This function is called once per genome
   (maternal and paternal) in every replicate and allocates both vectors each time, so the
   previous ones must be freed first (otherwise they leak). */
void FreeTriNucleotides (void)
	{
	int i;

	if (triNucleotideMaternal)
		{
		for (i=0; i<64; i++)
			free (triNucleotideMaternal[i].position);
		free (triNucleotideMaternal);
		triNucleotideMaternal = NULL;
		}
	if (triNucleotidePaternal)
		{
		for (i=0; i<64; i++)
			free (triNucleotidePaternal[i].position);
		free (triNucleotidePaternal);
		triNucleotidePaternal = NULL;
		}
	}

void SimulateSignatureISM (TreeNode *p, int genome, long int *seed)
    {
    int		i, numMutations, mutationsSoFar;
    double	totalBranchSum;
	int		newState, chosenSite;
	
	newState = -9;

	FreeTriNucleotides ();
	
	/* allocate memory to store the trinucleotide frequencies */
	triNucleotideMaternal = (TriNucStr *) Checked (calloc(64, sizeof(TriNucStr)), "Could not allocate the triNucleotideMaternal vector\n");
	for (i=0; i<64; i++)
		{
		triNucleotideMaternal[i].tempLength = 100;
		triNucleotideMaternal[i].numAvailablePositions = 0;
		triNucleotideMaternal[i].position = (int *) calloc(100,sizeof(int));
		}

	triNucleotidePaternal = (TriNucStr *) Checked (calloc(64, sizeof(TriNucStr)), "Could not allocate the triNucleotidePaternal vector\n");
	for (i=0; i<64; i++)
		{
		triNucleotidePaternal[i].tempLength = 100;
		triNucleotidePaternal[i].numAvailablePositions = 0;
		triNucleotidePaternal[i].position = (int *) calloc(100,sizeof(int));
		}

	/* count trinucleotide frequencies in the healthy root genome */
	CountTriNucFrequencies (data[genome][OUTGROUP_ROOT], genome);

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
		triMutationsCounter = (int *) Checked (calloc(96, sizeof(int)), "Could not allocate the triMutationsCounter vector\n");
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
	int j, k;
	int doDeconstructSigs = YES;
	FILE *fpCounts, *fpRefGenome;

	snprintf (File, sizeof(File),"%s/%s", resultsDir, "trimutcounts");
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
	snprintf (File, sizeof(File),"%s/%s", resultsDir, "contextGenome");
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

void CountTriNucFrequencies (Geno *genome_array, int genome)
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

/* Draws a signature, one of its 96 trinucleotide mutations and a strand; returns the index of the trinucleotide
   class where the mutation takes place and the new state of its central base */
static int DrawTrinucleotideClass (long int *seed, int *newState)
	{
	int		n1, n2, n3, r1, r2, r3, rest, signID, chosenSignature, chosenTriMutation;

	/* 1: Choose a mutation signature */
	signID = ChooseUniformState(signatureWeight, seed);
	chosenSignature = signatureID[signID];

	/* 2: Select a mutation out of the 96 possible ones for the chosen signature */
	chosenTriMutation = ChooseUniformState(signatureProbs[chosenSignature-1], seed);

	/* 3: Define in which of the 64 trinucleotide classes this mutation will take place */
	r1 = chosenTriMutation/24;
	rest = chosenTriMutation%24;
	r2 = (rest/4 < 3) ? C : T;
	r3 = rest%4;

	if (RandomUniform(seed) > 0.5)	/* plus strand */
		{
		n1 = r1;
		n2 = r2;
		n3 = r3;
		*newState = targetTriChange[rest/4];
		}
	else							/* minus strand */
		{
		n1 = complementBase[r3];
		n2 = complementBase[r2];
		n3 = complementBase[r1];
		*newState = complementBase[targetTriChange[rest/4]];
		}
	return trinuc(n1,n2,n3);
	}

int ChooseTrinucleotideSite (long int *seed, int *newState, int genome)
	{
	int			index, trindex, trials = 0;
	TriNucStr	*triNucleotide = (genome == MATERNAL) ? triNucleotideMaternal : triNucleotidePaternal;

	trindex = DrawTrinucleotideClass (seed, newState);
	while (triNucleotide[trindex].numAvailablePositions == 0)
		{
		trindex = DrawTrinucleotideClass (seed, newState);
		if (trials++ > 100*numSites)
			{
			fprintf (stderr, "\n\nERROR: after %d trials cannot find an non-mutated trinuc site, increase the number of sites simulated\n\n", 100*numSites);
			exit(-1);
			}
		}

	/* 3: Select the genomic position where this mutation will take place */
	/* note we assume for convenience that the context is constant along the phylogeny */
	if (triNucleotide[trindex].numAvailablePositions == 1)
		index = 0;
	else
		index = RandomUniformTo(triNucleotide[trindex].numAvailablePositions, seed);

	/* 4: remove the chosen position from the available list */
	{
	int chosenPosition = triNucleotide[trindex].position[index] + 1;

	RemoveTriNucObservation(&triNucleotide[trindex], index);
	return chosenPosition;
	}
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
            ancState = GetG(genome, anccell, site);
            cumBranchLength += p->branchLength;
				
            if (cumBranchLength < uniform || allSites[site].numMutations > 0)/* => there will be no change; forcing ISM diploid */
                {
                SetG(genome, cell, site, ancState);
 				//fprintf (stderr, "\ncopying ancstate %c", WhichNuc(ancstate));
               }
            else /* => there will be change */
                {
				SetG(genome, cell, site, newState);
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
		int *newPosition = (int *)realloc(trin->position, trin->tempLength * sizeof(int));
		if (!newPosition)
			{
			free (trin->position);
			fprintf (stderr, "Could not reallocate trinucleotide position array when removing\n");
			exit (-1);
			}
		trin->position = newPosition;
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
				SetG(genome, cell, site, GetG(genome, anccell, site));
			else /* => there will be change */
				{
				if (GetG(genome, anccell, site) == 0)
					SetG(genome, cell, site, 1);
				else
					SetG(genome, cell, site, 0);
					
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
			ancstate = GetG(genome, anccell, site);
	
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
			newstate = SetG(genome, cell, site, ChooseUniformState (Pij[ancstate], seed));
			
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
                SetG(genome, cell, site, GetG(genome, anccell, site));
                }
			else /* => there will be a CN_LOH event */
                {
				if (genome == MATERNAL)
					{
					SetG(genome, cell, site, GetG(PATERNAL, cell, site));
					allSites[site].numCNLOHmaternal++;
					}
				else
					{
					SetG(genome, cell, site, GetG(MATERNAL, cell, site));
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
                SetG(genome, cell, site, GetG(genome, anccell, site));
                }
			else /* => there will be a deletion */
                {
                if (GetG(genome, anccell, site) != DELETION)  /* checking all this might be excessive */
                    SetG(genome, cell, site, DELETION);
                else if (GetG(genome, anccell, site) == DELETION)
                    {
                    fprintf (stderr, "\n\n ERROR: site %d in genome %d of cell %d cannot be deleted twice under the ISM model\n", site, genome, cell);
                    exit(-1);
                    }
                else
                    {
                    fprintf (stderr, "\n\n ERROR: site %d in genome %d of cell %d has an unknow state %d under the ISM model\n", site, genome, anccell, GetG(genome, anccell, site));
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






/********************************** Segments: copy-number losses and gains ***********************************/
/*	Copy-number losses (--seg-loss-rate) and gains (--seg-gain-rate) of segments, simulated after the SNVs, cnLOH and point deletions.
	The number of events of a haplotype is a Poisson with mean (total tree length in generations) * rate; each event falls on a branch
	chosen according to its length, starts at a uniform site and covers a geometric number of sites (mean --seg-len, clipped at the end
	of the genome). Every cell below the branch carries the event in that haplotype.
	Loss: the sites become DELETION (genotypes, read counts and VCF treat them like point deletions).
	Gain: one extra copy of the haplotype (at most one: a haplotype has 1 or 2 copies, so the genotype of a gained segment is a
	dose 2:1 or 1:2); the genotype is the same, and the read depth of that haplotype is doubled (see CopyGainAt). */
typedef struct { int type, genome, start, len; TreeNode *branch; } CnvEvent;
typedef struct { int *s, *e, n, cap; } SegList;
static CnvEvent	*cnvEv = NULL;
static int		numCnvEv = 0, capCnvEv = 0;
static SegList	*gainSeg = NULL;
static int		gainSegNodes = 0;
int				copyGainOn = 0;

void ResetSegments (void)
	{
	int		k;

	if (gainSeg != NULL)
		{
		for (k=0; k<2*gainSegNodes; k++)
			{
			free (gainSeg[k].s);
			free (gainSeg[k].e);
			}
		free (gainSeg);
		gainSeg = NULL;
		}
	gainSegNodes = 0;
	copyGainOn = NO;
	numCnvEv = 0;
	numSegLoss = numSegGain = 0;
	}

/* 2 inside a gained segment of the haplotype (two copies: twice the reads), 1 elsewhere */
double CopyGainAt (int node, int genome, long site)
	{
	SegList	*v;
	int		k;

	if (gainSeg == NULL || node >= gainSegNodes)
		return 1.0;
	v = &gainSeg[2 * node + genome];
	for (k=0; k<v->n; k++)
		if (site >= v->s[k] && site < v->e[k])
			return 2.0;
	return 1.0;
	}

static void SegCollectNodes (TreeNode *p, TreeNode **list, int *n)
	{
	if (p == NULL)
		return;
	if (p->isHealthyRoot != YES && p->isHealthyTip != YES)	/* the branch of the healthy cell carries no tumour cell */
		list[(*n)++] = p;
	SegCollectNodes (p->left, list, n);
	SegCollectNodes (p->right, list, n);
	}

static void SegLoseBelow (TreeNode *p, int genome, int from, int to)
	{
	int		j;

	if (p == NULL)
		return;
	for (j=from; j<to; j++)
		SetG (genome, p->label, j, DELETION);
	SegLoseBelow (p->left, genome, from, to);
	SegLoseBelow (p->right, genome, from, to);
	}

static void SegGainBelow (TreeNode *p, int genome, int from, int to)
	{
	SegList	*v;

	if (p == NULL)
		return;
	v = &gainSeg[2 * p->label + genome];
	if (v->n == v->cap)
		{
		v->cap = v->cap ? 2 * v->cap : 4;
		v->s = (int *) realloc (v->s, v->cap * sizeof (int));
		v->e = (int *) realloc (v->e, v->cap * sizeof (int));
		if (!v->s || !v->e)
			{
			fprintf (stderr, "Could not allocate the gained segments\n");
			exit (-1);
			}
		}
	v->s[v->n] = from;
	v->e[v->n++] = to;
	SegGainBelow (p->left, genome, from, to);
	SegGainBelow (p->right, genome, from, to);
	}

static void EvolveSegments (TreeNode *p, int genome, long int *seed, int type, double rate)
	{
	TreeNode	**list;
	int			n = 0, e, numEvents, k, start, len;
	double		total = 0.0, u, cum;

	list = (TreeNode **) Checked (calloc (numNodes + 1, sizeof (TreeNode *)), "Could not allocate the list of branches of the segment events\n");
	SegCollectNodes (p, list, &n);
	for (k=0; k<n; k++)
		total += list[k]->branchLength;
	if (type == 1 && gainSeg == NULL)
		{
		gainSegNodes = numNodes + 1;
		gainSeg = (SegList *) Checked (calloc (2 * gainSegNodes, sizeof (SegList)), "Could not allocate the gained segments\n");
		copyGainOn = YES;
		}
	numEvents = RandomPoisson (total / mutationRate * rate, seed);
	for (e=0; e<numEvents; e++)
		{
		u = RandomUniform (seed) * total;
		cum = 0.0;
		for (k=0; k<n-1; k++)
			{
			cum += list[k]->branchLength;
			if (u < cum)
				break;
			}
		start = RandomUniformTo (numSites, seed);
		len = 1 + (int) (log (1.0 - RandomUniform (seed)) / log (1.0 - 1.0 / segLen));
		if (segLen <= 1.0)
			len = 1;
		if (start + len > numSites)
			len = numSites - start;
		if (type == 0)
			{
			SegLoseBelow (list[k], genome, start, start + len);
			numSegLoss++;
			}
		else
			{
			SegGainBelow (list[k], genome, start, start + len);
			numSegGain++;
			}
		if (numCnvEv == capCnvEv)
			{
			capCnvEv = capCnvEv ? 2 * capCnvEv : 16;
			cnvEv = (CnvEvent *) Checked (realloc (cnvEv, capCnvEv * sizeof (CnvEvent)), "Could not allocate the segment events\n");
			}
		cnvEv[numCnvEv].type = type;
		cnvEv[numCnvEv].genome = genome;
		cnvEv[numCnvEv].start = start;
		cnvEv[numCnvEv].len = len;
		cnvEv[numCnvEv++].branch = list[k];
		}
	free (list);
	}

void EvolveSegmentLossesOnTree (TreeNode *p, int genome, long int *seed)
	{
	EvolveSegments (p, genome, seed, 0, lossRate);
	}

void EvolveSegmentGainsOnTree (TreeNode *p, int genome, long int *seed)
	{
	EvolveSegments (p, genome, seed, 1, gainRate);
	}

static void SegLeafNames (TreeNode *p, char *buf, size_t cap, int *count)
	{
	char	name[4*MAX_NAME];

	if (p == NULL)
		return;
	if (p->left == NULL && p->right == NULL && p->label >= 0 && p->label < numCells)
		{
		if (doUserTree == NO)
			snprintf (name, sizeof (name), "%s%04d", inCellName, p->label + 1);
		else
			snprintf (name, sizeof (name), "%s", cellNames[p->label]);
		if (strlen (buf) + strlen (name) + 2 < cap)
			{
			if (*count > 0)
				strcat (buf, ",");
			strcat (buf, name);
			}
		(*count)++;
		}
	SegLeafNames (p->left, buf, cap, count);
	SegLeafNames (p->right, buf, cap, count);
	}

/* Truth of the segment losses and gains of a replicate: file cnv_events (or cnv_events_dir/cnv_events.NNNN with -v), one row per
   event with the type, the haplotype, the segment (1-based, inclusive) and the cells that carry it. */
void WriteCnvEvents (int replicate)
	{
	char	file[4*MAX_NAME], *buf;
	FILE	*fp;
	int		k, count;
	size_t	cap = 64 * (size_t) (numCells + 1) + 64;

	if (doPrintSeparateReplicates == YES)
		{
		snprintf (file, sizeof(file), "%s/cnv_events_dir", resultsDir);
		MAKE_DIR (file);
		snprintf (file, sizeof(file), "%s/cnv_events_dir/cnv_events.%04d", resultsDir, replicate+1);
		fp = fopen (file, "w");
		}
	else
		{
		snprintf (file, sizeof(file), "%s/cnv_events", resultsDir);
		fp = fopen (file, replicate == 0 ? "w" : "a");
		}
	if (!fp)
		{
		fprintf (stderr, "Can't open \"%s\"\n", file);
		exit (-1);
		}
	if (doPrintSeparateReplicates == NO)
		fprintf (fp, "[#%d]\n", replicate+1);
	fprintf (fp, "type\thaplotype\tstart\tend\tnum_cells\tcells\n");
	buf = (char *) Checked (malloc (cap), "Could not allocate the list of cells of the segment events\n");
	for (k=0; k<numCnvEv; k++)
		{
		buf[0] = '\0';
		count = 0;
		SegLeafNames (cnvEv[k].branch, buf, cap, &count);
		fprintf (fp, "%s\t%s\t%d\t%d\t%d\t%s\n", cnvEv[k].type == 0 ? "loss" : "gain", cnvEv[k].genome == MATERNAL ? "maternal" : "paternal",
			cnvEv[k].start + 1, cnvEv[k].start + cnvEv[k].len, count, buf);
		}
	free (buf);
	fclose (fp);
	}
