/*
	CellCoal - NGS read counts, genotype likelihoods and doublets
	See cellcoal.h for the shared declarations. This file was split from cellcoal.c.
*/
#include "cellcoal.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "definitions.h"

/********************* AllocateCellStructure  ************************/
/*
 This cell structure keeps information basically about read count functions. It was
 build initially for making doublets.
 */
void AllocateCellStructure(void)
	{
	int i;
	double limitMB = 256.0;
	char *env;

	/* Genotype likelihoods are only printed for SNV sites, and they are a pure function of the stored reads, so they
	   are kept only for SNV sites and only within a memory budget (CELLCOAL_GL_CACHE_MB, default 256; 0 = keep none);
	   anything not kept is recomputed when printed, with the same result. */
	env = getenv ("CELLCOAL_GL_CACHE_MB");
	if (env && *env)
		limitMB = atof (env);
	glBudget = (size_t) (limitMB * 1048576.0);
	glUsed = 0;
	siteGL = NULL;
	if (doNGS == YES)
		{
		siteGL = (double**) calloc (numSites, sizeof(double*));
		if (!siteGL)
			{
			fprintf (stderr, "Could not allocate the siteGL structure\n");
			exit (-1);
			}
		}

	cell = (CellStr*) calloc(numCells+1, sizeof(CellStr));
	if (!cell)
		{
		fprintf (stderr, "Could not allocate cell (%ld)\n", (numCells+1) * sizeof(CellStr));
		exit (-1);
		}

	/* all the per-site data (read counts, genotype likelihoods) lives inside CellSiteStr, so one block per cell is enough;
	   these blocks are only allocated when an option needs them (see NEED_CELL_SITES) */
	if (!NEED_CELL_SITES)
		return;
	if (streamMode == YES)
		{
		/* only a window of sites is held (see StreamReads): about 64 MB, and enough sites to keep the threads busy */
		long	b = (long) (64.0 * 1048576.0 / (sizeof(CellSiteStr) * (double)(numCells+1)));
		int		minB = 16;
#ifdef _OPENMP
		minB = 8 * omp_get_max_threads();
#endif
		if (b < minB) b = minB;
		if (b > 4096) b = 4096;
		if (b > numSites) b = numSites;
		streamBlock = (int) b;
		}
	for (i=0; i<numCells+1; i++)
		{
		cell[i].site  = (CellSiteStr*) calloc ((streamMode == YES) ? streamBlock : numSites, sizeof(CellSiteStr));
		if (!cell[i].site)
			{
			fprintf (stderr, "Could not allocate the cell[%d].site structure\n", i);
			exit (-1);
			}
		}
	}

/* Total read counts of the four nucleotides at site j (doublets count with their own reads) and the alternate alleles seen in reads */
static void TotalReadCountsAtSite (int j)
	{
	int		i, numAltAlleles;
	int		countReadsA, countReadsC, countReadsG, countReadsT;
	countReadsA = countReadsC = countReadsG = countReadsT = 0;
	for (i=0; i<numCells+1; i++)
		{
		if (cell[i].hasDoublet == YES)
			{
			countReadsA += CD(i,j).readCountDoublet[A];
			countReadsC += CD(i,j).readCountDoublet[C];
			countReadsG += CD(i,j).readCountDoublet[G];
			countReadsT += CD(i,j).readCountDoublet[T];
			}
		else
			{
			countReadsA += CS(i,j).readCount[A];
			countReadsC += CS(i,j).readCount[C];
			countReadsG += CS(i,j).readCount[G];
			countReadsT += CS(i,j).readCount[T];
			}
		}
	
	numAltAlleles = 0;
	if (countReadsA > 0 && allSites[j].referenceAllele != A)
		allSites[j].alternateAlleles[numAltAlleles++] = A;
	if (countReadsC > 0 && allSites[j].referenceAllele != C)
		allSites[j].alternateAlleles[numAltAlleles++] = C;
	if (countReadsG > 0 && allSites[j].referenceAllele != G)
		allSites[j].alternateAlleles[numAltAlleles++] = G;
	if (countReadsT > 0 && allSites[j].referenceAllele != T)
		allSites[j].alternateAlleles[numAltAlleles++] = T;
	allSites[j].numAltAlleles = numAltAlleles;
			
	allSites[j].readCountA = countReadsA;
	allSites[j].readCountC = countReadsC;
	allSites[j].readCountG = countReadsG;
	allSites[j].readCountT = countReadsT;
	allSites[j].readCountACGT = countReadsA + countReadsC + countReadsG + countReadsT;
	}

static void DoubletAtSite (int c1, int c2, int j, const double *gl1, double *probs, double **ngsEij, double **ampEijmat, double **ampEijpat, long int *seed);
static void TotalReadCountsAtSite (int j);
int MLCountsAtSite (int j);


void NewReadScratch (ReadScratch *rs)
	{
	int		t;

	rs->ampMat = (double**) calloc (4, sizeof(double*));
	rs->ampPat = (double**) calloc (4, sizeof(double*));
	rs->probsT = (double*) calloc (4, sizeof(double));
	rs->sbuf = (double*) malloc ((size_t)(numCells+1) * GL_NUM * sizeof(double));
	if (!rs->ampMat || !rs->ampPat || !rs->sbuf || !rs->probsT)
		{
		fprintf (stderr, "Could not allocate the thread amplification matrices\n");
		exit (-1);
		}
	for (t=0; t<4; t++)
		{
		rs->ampMat[t] = (double*) calloc (4, sizeof(double));
		rs->ampPat[t] = (double*) calloc (4, sizeof(double));
		if (!rs->ampMat[t] || !rs->ampPat[t])
			{
			fprintf (stderr, "Could not allocate the thread amplification matrices\n");
			exit (-1);
			}
		}
	}

void FreeReadScratch (ReadScratch *rs)
	{
	int		t;

	for (t=0; t<4; t++)
		{
		free (rs->ampMat[t]);
		free (rs->ampPat[t]);
		}
	free (rs->ampMat);
	free (rs->ampPat);
	free (rs->probsT);
	free (rs->sbuf);
	}

/* Simulates the reads of all the cells at site jj (and their genotype likelihoods) from the stream of the site. When the
   reads are streamed, it also fills the true genotypes, adds the doublets and summarizes the site. */
void SimulateSiteReads (int jj, uint64_t streamBase, double **ngsEij, ReadScratch *rs)
	{
	double	**ampMat = rs->ampMat, **ampPat = rs->ampPat, *sbuf = rs->sbuf, *probsT = rs->probsT;
	long int	dummySeed = 1;	/* never used: the stream of the site is active */
	int		ii;
	int		hasAlt = NO, al;
	size_t	bytes = (size_t)(numCells+1) * GL_NUM * sizeof(double), before;

	/* true genotypes (before ADO and genotyping errors) of the streamed window */
	if (streamMode == YES)
		for (ii=0; ii<numCells+1; ii++)
			{
			CS(ii,jj).trueMaternalAllele = GetTrueG (MATERNAL, ii, jj);
			CS(ii,jj).truePaternalAllele = GetTrueG (PATERNAL, ii, jj);
			}

	/* read counts of all the cells at this site, from the stream of the site */
	StartStream (streamBase, (uint64_t) jj);
	/* coverage varies between sites (-O): the depth of every cell at this site is multiplied by a gamma variate with mean 1 */
	siteCovMult = (alphaSiteCoverage > 0.0) ? RandomGamma (alphaSiteCoverage, &dummySeed) / alphaSiteCoverage : 1.0;
	for (ii=0; ii<numCells+1; ii++)
		{
		SiteReadCounts (&CS(ii,jj), ii, jj, probsT, ngsEij, ampMat, ampPat, &dummySeed);
		if (CS(ii,jj).numReads <= 0)
			CS(ii,jj).MLmatAllele = CS(ii,jj).MLpatAllele = MISSING;
		}

	if (readsOnly == YES)
		{
		/* the bulk file needs reads only: skip the genotype likelihoods (they use no random numbers) */
		if (doubletDonor != NULL)
			for (ii=0; ii<numCells; ii++)
				if (cell[ii].hasDoublet == YES)
					DoubletAtSite (ii, doubletDonor[ii], jj, sbuf, probsT, ngsEij, ampMat, ampPat, &dummySeed);
		streamActive = 0;
		return;
		}

	for (ii=0; ii<numCells+1; ii++)
		{
		double *g = sbuf + GL_NUM*(size_t)ii;
		if (CS(ii,jj).numReads > 0)
			{
			/* rebuild the amplification error matrices of this site exactly as SiteReadCounts did */
			BuildAmpMatrices (&CS(ii,jj), ampMat, ampPat);
			GenotypeLikelihoods (&CS(ii,jj), ii, jj, probsT, ngsEij, ampMat, ampPat, &dummySeed, g);
			}
		else
			memset (g, 0, GL_NUM*sizeof(double));
		/* same criterion as CountAllelesInMLGenotypes: the site is a SNV when an ML allele differs from the reference */
		al = CS(ii,jj).MLmatAllele;
		if (al >= A && al <= T && al != allSites[jj].referenceAllele)
			hasAlt = YES;
		al = CS(ii,jj).MLpatAllele;
		if (al >= A && al <= T && al != allSites[jj].referenceAllele)
			hasAlt = YES;
		}
	/* doublets: the reads of a donor cell are added to those of the doublet cell, from the same stream of the site
	   (stored and streamed reads therefore give the same data) */
	if (doubletDonor != NULL)
		for (ii=0; ii<numCells; ii++)
			if (cell[ii].hasDoublet == YES)
				DoubletAtSite (ii, doubletDonor[ii], jj, sbuf + GL_NUM*(size_t)ii, probsT, ngsEij, ampMat, ampPat, &dummySeed);
	streamActive = 0;

	if (hasAlt == YES)
		{
		#pragma omp atomic capture
			{ before = glUsed; glUsed += bytes; }
		if (before + bytes <= glBudget)
			{
			double *keep = (double*) malloc (bytes);
			if (keep)
				{
				memcpy (keep, sbuf, bytes);
				siteGL[jj] = keep;
				}
			}
		}

	/* summary of the site (what GenerateReadCounts and CountAllelesInMLGenotypes do for all the sites) */
	if (streamMode == YES)
		{
		TotalReadCountsAtSite (jj);
		MLCountsAtSite (jj);
		}
	}

/* Sequencing error probabilities, constant across sites, and the genotype likelihood model */
void SetupReadModel (double **ngsEij)
	{
	int		i, j;

	/* initialize ngsEij now as this will be contant across sites */
	for (i=0; i<4; i++)
		for (j=0; j<4; j++)
			{
			if (i == j)
				ngsEij[i][j] = 1.0 - sequencingError;
			else
				ngsEij[i][j] = Eij[i][j]/cumEij[i][3] * sequencingError;
			}
	
	
	/* define GL model  */
	/* initialize */
	doGATK = NO;
	do2T = NO;
	do4T = NO;
	doGATK_ADO = NO;
	do2T_ADO = NO;
	do4T_ADO = NO;

	if (thereisADO == NO)
		{
		if (meanAmplificationError == 0)
			{
			doGATK = YES;
			strcpy(GLmodel, "GATK");
			}
		else if (simulateOnlyTwoTemplates == YES)
			{
			do2T = YES;
			strcpy(GLmodel, "2T");
			}
		else
			{
			do4T = YES;
			strcpy(GLmodel, "4T");
			}
		}
	else
		{
		if (meanAmplificationError == 0)
			{
			doGATK_ADO = YES;
			strcpy(GLmodel, "GATK_ADO");
			}
		else if (simulateOnlyTwoTemplates == YES)
			{
			do2T_ADO = YES;
			strcpy(GLmodel, "2T_ADO");
			}
		else
			{
			do4T_ADO = YES;
			strcpy(GLmodel, "4T_ADO");
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
	int			i, j;
	double		*probs, **ngsEij, **ampEijmat, **ampEijpat;

	probs = (double*) calloc (4, sizeof(double));
	if (!probs)
		{
		fprintf (stderr, "Could not allocate the probs structure\n");
		exit (-1);
		}
	
	/* allocate error probabilities for NGS */
	ngsEij = (double**) calloc (4, sizeof(double*));
	if (!ngsEij)
		{
		fprintf (stderr, "Could not allocate the ngsEij structure\n");
		exit (-1);
		}
	for (i=0; i<4; i++)
		{
		ngsEij[i] = (double*) calloc (4, sizeof(double));
		if (!ngsEij[i])
			{
			fprintf (stderr, "Could not allocate the ngsEij[i] structure\n");
			exit (-1);
			}
		}

	/* allocate error probabilities for amplification */
	ampEijmat = (double**) calloc (4, sizeof(double*));
	if (!ampEijmat)
		{
		fprintf (stderr, "Could not allocate the ampEijmat structure\n");
		exit (-1);
		}
	for (i=0; i<4; i++)
		{
		ampEijmat[i] = (double*) calloc (4, sizeof(double));
		if (!ampEijmat[i])
			{
			fprintf (stderr, "Could not allocate the ampEijmat[i] structure\n");
			exit (-1);
			}
		}

	ampEijpat = (double**) calloc (4, sizeof(double*));
	if (!ampEijpat)
		{
		fprintf (stderr, "Could not allocate the ampEijpat structure\n");
		exit (-1);
		}
	for (i=0; i<4; i++)
		{
		ampEijpat[i] = (double*) calloc (4, sizeof(double));
		if (!ampEijpat[i])
			{
			fprintf (stderr, "Could not allocate the ampEijpat[i] structure\n");
			exit (-1);
			}
		}

	SetupReadModel (ngsEij);

	/* One stream per site is started from two numbers taken from the main stream (so each replicate
	   differs), which makes the sites independent of each other. For every site the reads of all cells
	   are simulated (in cell order, with the stream of the site) and then the genotype likelihoods are
	   computed, all in parallel when compiled with OpenMP. The result is the same for any number of threads.
	   Each thread needs its own scratch matrices. */
	{
	uint64_t	streamBase;
	long int	s1, s2;

	SetupDoublets (seed, numSites);
	SetupCellCoverage (seed);
	RandomUniform (seed);
	s1 = *seed;
	RandomUniform (seed);
	s2 = *seed;
	streamBase = ((uint64_t) s1 << 32) ^ (uint64_t) s2;

	#pragma omp parallel
		{
		ReadScratch	rs;
		int		jj;

		NewReadScratch (&rs);
		#pragma omp for schedule(static, 8)
		for (jj=0; jj<numSites; jj++)
			SimulateSiteReads (jj, streamBase, ngsEij, &rs);
		FreeReadScratch (&rs);
		}
	}
	free (doubletDonor);	/* the doublets were added site by site above */
	doubletDonor = NULL;

	/* count alternate alleles in total reads */
	for (j=0; j<numSites; j++)
		TotalReadCountsAtSite (j);

	free (probs);
	for (i=0; i<4; i++)
		{
		free (ngsEij[i]);
		free (ampEijmat[i]);
		free (ampEijpat[i]);
		}
	free (ngsEij);
	free (ampEijmat);
	free (ampEijpat);
	}

/* ******************** SiteReadCounts  ************************/
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
	double 	allelicImbalance;
	
	readCount = c->readCount;
	
	thereIsMaternalAllele = YES;
	thereIsPaternalAllele = YES;
	
	maternalAllele = GetG(MATERNAL, i, j);
	if (maternalAllele == ADO || maternalAllele == DELETION)
		thereIsMaternalAllele = NO;
	
	paternalAllele = GetG(PATERNAL, i, j);
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
	
	/* sample maternal allelic imbalance for this genotype */
	if (doAllelicImbalance == YES)
		allelicImbalance = RandomBetaMeanVar(meanAllelicImbalance, varAllelicImbalance, seed);
	else
		allelicImbalance = 0.5;

	/* number of reads will follow a Poisson or a Negative Binomial distribution (overdispersed coverage) around mean coverage, and are randomly assigned to maternal/paternal chromosome */
	numMaternalReads = 0;
	numPaternalReads = 0;
	
	if (rateVarCoverage == YES) /* Negative Binomial */
		{
		if (thereIsMaternalAllele == YES && thereIsPaternalAllele == YES) /* both alleles are present */
			{
			numMaternalReads = RandomNegativeBinomial(allelicImbalance*coverage*siteCovMult*cellCovMult[i], alphaCoverage, seed);
			numPaternalReads = RandomNegativeBinomial((1.0-allelicImbalance)*coverage*siteCovMult*cellCovMult[i], alphaCoverage, seed);
			}
		else if (thereIsMaternalAllele == YES && thereIsPaternalAllele == NO)  /* only maternal allele */
			{
			numMaternalReads = RandomNegativeBinomial(haploidCoverageReduction*coverage*siteCovMult*cellCovMult[i], alphaCoverage, seed);
			numPaternalReads = 0;
			}
		else if (thereIsMaternalAllele == NO && thereIsPaternalAllele == YES)  /* only paternal allele */
			{
			numMaternalReads = 0;
			numPaternalReads = RandomNegativeBinomial(haploidCoverageReduction*coverage*siteCovMult*cellCovMult[i], alphaCoverage, seed);
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
			numMaternalReads = RandomPoisson(allelicImbalance*coverage*siteCovMult*cellCovMult[i], seed);
			numPaternalReads = RandomPoisson((1.0-allelicImbalance)*coverage*siteCovMult*cellCovMult[i], seed);
			}
		else if (thereIsMaternalAllele == YES && thereIsPaternalAllele == NO)  /* only maternal allele */
			{
			numMaternalReads = RandomPoisson(haploidCoverageReduction*coverage*siteCovMult*cellCovMult[i], seed);
			numPaternalReads = 0;
			}
		else if (thereIsMaternalAllele == NO && thereIsPaternalAllele == YES)  /* only paternal allele */
			{
			numMaternalReads = 0;
			numPaternalReads = RandomPoisson(haploidCoverageReduction*coverage*siteCovMult*cellCovMult[i], seed);
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
			maternalSiteAmplificationError = RandomBetaMeanVar(meanAmplificationError, varAmplificationError, seed);
			paternalSiteAmplificationError = RandomBetaMeanVar(meanAmplificationError, varAmplificationError, seed);
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
	c->maternalAllele = maternalAllele;
	c->paternalAllele = paternalAllele;
	c->thereIsMaternalAllele = thereIsMaternalAllele;
	c->thereIsPaternalAllele = thereIsPaternalAllele;
	c->maternalSiteAmplificationError = maternalSiteAmplificationError;
	c->paternalSiteAmplificationError = paternalSiteAmplificationError;
	}


/* amplification error matrices of a site, exactly as SiteReadCounts builds them */
void BuildAmpMatrices (CellSiteStr *c, double **ampMat, double **ampPat)
	{
	int k, l;

	for (k=0; k<4; k++)
		for (l=0; l<4; l++)
			{
			if (k == l)
				{
				ampMat[k][l] = 1.0 - c->maternalSiteAmplificationError;
				ampPat[k][l] = 1.0 - c->paternalSiteAmplificationError;
				}
			else
				{
				ampMat[k][l] = Eij[k][l]/cumEij[k][3] * c->maternalSiteAmplificationError;
				ampPat[k][l] = Eij[k][l]/cumEij[k][3] * c->paternalSiteAmplificationError;
				}
			}
	}

/* Recomputes the genotype likelihoods of a cell and site from its stored data (reads and amplification errors).
   It needs no random numbers, only reads shared data, and is safe to call from several threads. */
void SiteGL (CellSiteStr *c, int i, int j, double *glOut)
	{
	double	ngsRows[4][4], ampMRows[4][4], ampPRows[4][4];
	double	*ngsP[4], *ampMP[4], *ampPP[4], probs[4];
	int		k, l;

	for (k=0; k<4; k++)
		{
		ngsP[k] = ngsRows[k];  ampMP[k] = ampMRows[k];  ampPP[k] = ampPRows[k];
		for (l=0; l<4; l++)
			{
			if (k == l)
				ngsRows[k][l] = 1.0 - sequencingError;
			else
				ngsRows[k][l] = Eij[k][l]/cumEij[k][3] * sequencingError;
			}
		}
	BuildAmpMatrices (c, ampMP, ampPP);
	GenotypeLikelihoods (c, i, j, probs, ngsP, ampMP, ampPP, NULL, glOut);
	}

/* Returns the genotype likelihoods of cell i at site j: the stored ones when kept for this site, otherwise recomputed
   into buf. Sites without reads have all-zero likelihoods, as when they were never computed. */
const double *GetSiteGL (int i, int j, double *buf)
	{
	if (siteGL && siteGL[j])
		return siteGL[j] + GL_NUM*(size_t)i;
	if (CS(i,j).numReads > 0)
		SiteGL (&CS(i,j), i, j, buf);
	else
		memset (buf, 0, GL_NUM*sizeof(double));
	return buf;
	}

/* pow(10,x) without calling the library when the answer is known exactly:
   10^0 = 1, and 10^x underflows to +0 for x < -400 (the smallest double is ~10^-324) */
static inline double Pow10Fast (double x)
	{
	if (x == 0.0)
		return 1.0;
	if (x < -400.0)
		return 0.0;
	if (fastMode)
		return Pow10Approx (x);
	return pow (10.0, x);
	}

/* ******************** GenotypeLikelihoods  ************************/
/*
 Calculate genotype likelihoods (GLs) for a given site and cell given the simulated reads
 according to three different models
 
 **	Genotype likelihoods can be calculated given the read counts and the sequencing error
 (see Korneliussen 2013 BMCBioinf)

 Calculations are done GLs for all 16 phased genotypes, as different errors can occurs for maternal and paternal alleles, so the phase matters. But then we will report the GLs for the 10 unphased genotypes
 00 01 02 03 11 12 13 22 23 33
 AA AC AG AT CC CG CT GG GT TT
 */

/* log10 of x, reusing the result of an earlier call with exactly the same argument (keys/vals hold the
   previous arguments and results, n their number). log10 is a pure function, so the result is bit-identical;
   with symmetric error matrices (the default) many of the arguments of one cell and site coincide, and log10
   is the most expensive operation of the likelihood computation. */
static inline double Log10Memo (double x, double *keys, double *vals, int *n)
	{
	int m;

	for (m=0; m<*n; m++)
		if (keys[m] == x)
			return vals[m];
	keys[*n] = x;
	vals[*n] = fastMode ? Log10Approx (x) : log10 (x);
	return vals[(*n)++];
	}


void GenotypeLikelihoods (CellSiteStr *c, int i, int j, double *probs, double **ngsEij, double **ampEijmat, double **ampEijpat, long int *seed, double *glOut)
	{
	int		k, l, a1, a2;
	int		read, numReads;
	int		*readCount;
	int		maternalAllele, paternalAllele;
	double  maternalSiteAmplificationError, paternalSiteAmplificationError;
	int		template, template1, template2;
	int		debug_GL, debug_GL_EXIT, fix_parameters;
	double	maxLike, pReadGivenA1, pReadGivenA2;
	double	homozygote_for_read, heterozygote_for_read, no_allele_for_read;
	double	(*logGL)[4], (*normLogGL)[4];
	double 	MLmatAllele, MLpatAllele;
	double	alleleDropoutRate;
	long double	term1, term2, term3, maxterm, num1, num2, num3, number1, number2;
	double 	GL[4][4] = {{0}};
	double	pA1tab[4][4], pA2tab[4][4], lgA1tab[4][4], lgA2tab[4][4];
	long double	T2a[4], T3a[4];
	double	logGLm[4][4] = {{0}}, normLogGLm[4][4] = {{0}};
	double	L1tab[4][4][4], memoKeys[32], memoVals[32];
	int		memoN;
	
	numReads = c->numReads;
	readCount = c->readCount;
	logGL = logGLm;
	normLogGL = normLogGLm;
	maternalAllele = c->maternalAllele;
	paternalAllele = c->paternalAllele;
	maternalSiteAmplificationError = c->maternalSiteAmplificationError;
	paternalSiteAmplificationError = c->paternalSiteAmplificationError;

	//!!!!: debug GL
	debug_GL = NO;
	debug_GL_EXIT = YES;
	fix_parameters = NO;
	if (fix_parameters == YES)
		{
		maternalAllele = c->maternalAllele = T;
		paternalAllele = c->paternalAllele = G;
		c->thereIsMaternalAllele = YES;
		c->thereIsPaternalAllele = YES;
		readCount[A] = 0;
		readCount[C] = 0;
		readCount[G] = 14;
		readCount[T] = 0;

		alleleDropoutRate = 0.0;
		thereisADO = NO;
		doADOcell = NO;
		doADOsite = NO;
		fixedADOrate = alleleDropoutRate;
		//doGATK_ADO = YES;
		
		simulateOnlyTwoTemplates = NO;
		numReads = c->numReads = readCount[A] + readCount[C] + readCount[G] + readCount[T];
		sequencingError = 0.00;
		maternalSiteAmplificationError = c->maternalSiteAmplificationError = 0;
		paternalSiteAmplificationError = c->paternalSiteAmplificationError = 0;

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
	
	/* find ADO for this site */
	if (thereisADO == NO)
		alleleDropoutRate = 0;
	else
		{
		if (doADOcell == NO && doADOsite == NO)
			alleleDropoutRate = fixedADOrate;
		else if (doADOcell == YES && doADOsite == NO)
			alleleDropoutRate = meanADOcell;
		else if (doADOcell == NO && doADOsite == YES)
			alleleDropoutRate = meanADOsite;
		else
			alleleDropoutRate = meanADOcell + meanADOsite - (meanADOcell * meanADOsite);
		}
	
	
	/************************** MODEL 0 (GATK) *************************************/
	/* MODEL 0: genotype log10 likelihoods according to the observed read assuming a single sequencing error independent of the nucleotides involved, and no amplification error */
	/* extended from Korneliussen 2013 BMCBioinf
	 p(D|G=A1,A2) = multiply over reads [ 1/2 p(read i|A1) + 1/2 p(read i|A2)]
	 where p(read|A) e/3 if b!=A  or 1-e if b=A
	 N.B: ANGSD seems to use natural logs!
	 */
	if (doGATK == YES)
		{
		homozygote_for_read = log10(1.0 - sequencingError);
		no_allele_for_read = log10(sequencingError/3.0);
		heterozygote_for_read = log10((1.0 - sequencingError)/2.0 + sequencingError/6.0);
		
		for (a1=0; a1<4; a1++) //a1 is maternal
			for (a2=a1; a2<4; a2++) //a2 is paternal
				{
				logGL[a1][a2] = -0.0;
				for (k=0; k<4; k++)
					{
					if (readCount[k] > 0)
						{
						if (k == a1 && k == a2)
							logGL[a1][a2] += readCount[k] * homozygote_for_read;
						else if (k != a1 && k != a2)
							logGL[a1][a2] += readCount[k] * no_allele_for_read;
						else
							logGL[a1][a2] += readCount[k] * heterozygote_for_read;
							//fprintf (stderr, "\nGATK read=%c gl[%c][%c] = %lf", WhichNuc(k), WhichNuc(a1), WhichNuc(a2), logGL[a1][a2]);
						}
					}
				logGL[a2][a1] = logGL[a1][a2];
				}
		
		if (debug_GL == YES)
			{
			fprintf (stderr,"\ncell %d site %d: %c%c", i+1, j+1, WhichNuc(maternalAllele), WhichNuc(paternalAllele));
			fprintf (stderr,"|\treads: A:%d C:%d G:%d T:%d  | total:%d ",  readCount[A], readCount[C], readCount[G], readCount[T], numReads);
			fprintf (stderr,"\nsequencing error = %f",sequencingError);
			fprintf (stderr,"\nmatAmpError = %f", maternalSiteAmplificationError);
			fprintf (stderr,"\npatAmpError = %f", paternalSiteAmplificationError);
			for (a1=0; a1<4; a1++)
				for (a2=a1; a2<4; a2++)
					fprintf (stderr, "\n(GATK)log10 GL[%c][%c] = %lf", WhichNuc(a1), WhichNuc(a2), logGL[a1][a2]);
			fprintf (stderr, "\n");
			}
		} // model 0
	
	
	/************************** MODEL 1 (4-templates) *************************************/
	/* MODEL 1 (4-templates): genotype log10 likelihoods according to the observed read, assuming a single sequencing error independent of the nucleotides involved, and amplification error with all 4 templates */
	if (do4T == YES)
		{
		for (a1=0; a1<4; a1++) //a1 is maternal
			for (a2=0; a2<4; a2++) //a2 is paternal
				{
				logGL[a1][a2] = -0.0;
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
						logGL[a1][a2] += readCount[read] * log10 (meanAllelicImbalance * pReadGivenA1 + (1.0-meanAllelicImbalance) * pReadGivenA2);
						}
					//fprintf (stderr, "\n***read=%c gl[%c][%c] = %lf  (p1=%lf p2=%lf)", WhichNuc(read), WhichNuc(a1), WhichNuc(a2), logGL[a1][a2],pReadGivenA1,pReadGivenA2);
					}
				}
		
		if (debug_GL == YES)
			{
			fprintf (stderr,"\ncell %d site %d: %c%c", i+1, j+1, WhichNuc(maternalAllele), WhichNuc(paternalAllele));
			fprintf (stderr,"|\treads: A:%d C:%d G:%d T:%d  | total:%d",  readCount[A], readCount[C], readCount[G], readCount[T], numReads);
			fprintf (stderr,"\nsequencing error = %f",sequencingError);
			fprintf (stderr,"\nmatAmpError = %f", maternalSiteAmplificationError);
			fprintf (stderr,"\npatAmpError = %f", paternalSiteAmplificationError);
			for (a1=0; a1<4; a1++)
				for (a2=0; a2<4; a2++)
					fprintf (stderr, "\n(4T)log10 GL[%c][%c] = %lf", WhichNuc(a1), WhichNuc(a2), logGL[a1][a2]);
			fprintf (stderr, "\n");
			}
		} //model 1
	
	
	/************************** MODEL 2 (2-templates) *************************************/
	/* MODEL 2 (2-templates): genotype log10 likelihoods according to the observed read, assuming a single sequencing error independent of the nucleotides involved, and amplification error with up to 2 templates */
	if (do2T == YES)
		{
		for (a1=0; a1<4; a1++) //a1 is maternal
			for (a2=0; a2<4; a2++) //a2 is paternal
				{
				logGL[a1][a2] = 0.0;
				for (read=0; read<4; read++)
					{
					if (readCount[read] > 0)
						{
						pReadGivenA1 = pReadGivenA2 = 0.0;
						for (template1=0; template1<4; template1++)
							if (template1 != a1)
								{
								pReadGivenA1 += Eij[a1][template1]/cumEij[a1][3] * ((1.0 - maternalSiteAmplificationError) * ngsEij[a1][read] + maternalSiteAmplificationError * ngsEij[template1][read]);
								//fprintf (stderr, "\nread=%c t=[%c][%c] p1=%lf Eijmat[%c][%c]=%lf", WhichNuc(read), WhichNuc(a1), WhichNuc(template1), pReadGivenA1, WhichNuc(a1), WhichNuc(template1), Eij[a1][template1]/cumEij[a1][3]);
								}
						for (template2=0; template2<4; template2++)
							if (template2 != a2)
								{
								pReadGivenA2 += Eij[a2][template2]/cumEij[a2][3] * ((1.0 - paternalSiteAmplificationError) * ngsEij[a2][read] + paternalSiteAmplificationError * ngsEij[template2][read]);
								//fprintf (stderr, "\nread=%c t=[%c][%c] p2=%lf Eijpat[%c][%c]=%lf", WhichNuc(read), WhichNuc(a2), WhichNuc(template2), pReadGivenA2, WhichNuc(a2), WhichNuc(template2), Eij[a2][template2]/cumEij[a2][3]);
								}
						logGL[a1][a2] += readCount[read] * log10(meanAllelicImbalance * pReadGivenA1 + (1.0 - meanAllelicImbalance) * pReadGivenA2);
						}
					//fprintf (stderr, "\n***read=%c gl[%c][%c] = %lf  (p1=%lf p2=%lf)", WhichNuc(read), WhichNuc(a1), WhichNuc(a2), logGL[a1][a2],pReadGivenA1,pReadGivenA2);
					}//read
				} //a1,a2
		
		if (debug_GL == YES)
			{
			fprintf (stderr,"\ncell %d site %d: %c%c", i+1, j+1, WhichNuc(maternalAllele), WhichNuc(paternalAllele));
			fprintf (stderr,"|\treads: A:%d C:%d G:%d T:%d  | total:%d",  readCount[A], readCount[C], readCount[G], readCount[T], numReads);
			fprintf (stderr,"\nsequencing error = %f",sequencingError);
			fprintf (stderr,"\nmatAmpError = %f", maternalSiteAmplificationError);
			fprintf (stderr,"\npatAmpError = %f", paternalSiteAmplificationError);
			for (a1=0; a1<4; a1++)
				for (a2=0; a2<4; a2++)
					fprintf (stderr, "\n(2T)log10 GL[%c][%c] = %lf", WhichNuc(a1), WhichNuc(a2), logGL[a1][a2]);
			fprintf (stderr, "\n");
			}
		} // model 2
	
	

/* ************************* MODEL GENERAL *************************************/
/* MODEL: genotype log10 likelihoods according to the observed reads, assuming a single sequencing error independent of the nucleotides involved, amplification errors and ADO. Note that ADO does not change the likelihoods

	M = multiply over reads
	b = read

	p1 = p(b_i|A1)
	p2 = p(b_i|A2)
	where p(b_i|A) is defined by the specific amplification model (GATK, 2T or 4T)
	
	P(D|G=A1,A2) = (1-ADO) * M[1/2 p1 + 1/2 p2] + ADO/2 * M[p1] + ADO/2 * M[p2]
 
	term1 = M[1/2 p1 + 1/2 p2]
	term2 = M[p1]
	term3 = M[p2]
 
	M is the product over reads
*/
	if (doGATK_ADO == YES || do4T_ADO == YES || do2T_ADO == YES)
		{
		/* p(read|allele) and its log10 depend only on (allele, read), not on the genotype pair: compute them once */
		for (read=0; read<4; read++)
			{
			if (readCount[read] <= 0)
				continue;
			memoN = 0;
			for (k=0; k<4; k++)
				{
				pReadGivenA1 = pReadGivenA2 = 0.0;
				if (doGATK_ADO == YES)
					{
					pReadGivenA1 += ngsEij[k][read];
					pReadGivenA2 += ngsEij[k][read];
					}
				else if (do2T_ADO == YES)
					{
					for (template1=0; template1<4; template1++)
						if (template1 != k)
							pReadGivenA1 += Eij[k][template1]/cumEij[k][3] * ((1.0 - maternalSiteAmplificationError) * ngsEij[k][read] + maternalSiteAmplificationError * ngsEij[template1][read]);
					for (template2=0; template2<4; template2++)
						if (template2 != k)
							pReadGivenA2 += Eij[k][template2]/cumEij[k][3] * ((1.0 - paternalSiteAmplificationError) * ngsEij[k][read] + paternalSiteAmplificationError * ngsEij[template2][read]);
					}
				else
					{
					for (template=0; template<4; template++)
						{
						pReadGivenA1 += ampEijmat[k][template] * ngsEij[template][read];
						pReadGivenA2 += ampEijpat[k][template] * ngsEij[template][read];
						}
					}
				pA1tab[k][read] = pReadGivenA1;
				pA2tab[k][read] = pReadGivenA2;
				if (alleleDropoutRate != 0)
					{
					lgA1tab[k][read] = Log10Memo (pReadGivenA1, memoKeys, memoVals, &memoN);
					lgA2tab[k][read] = Log10Memo (pReadGivenA2, memoKeys, memoVals, &memoN);
					}
				}
			}

		/* term2 depends only on a1 and term3 only on a2: sum them once per allele (same additions, same order) */
		for (k=0; k<4; k++)
			{
			T2a[k] = 0.0;
			T3a[k] = 0.0;
			if (alleleDropoutRate != 0)
				for (read=0; read<4; read++)
					if (readCount[read] > 0)
						{
						T2a[k] += readCount[read] * lgA1tab[k][read];
						T3a[k] += readCount[read] * lgA2tab[k][read];
						}
			}

		/* log10 (1/2 p1 + 1/2 p2) for every read and genotype, computing each distinct argument only once */
		if (alleleDropoutRate != 1)
			for (read=0; read<4; read++)
				if (readCount[read] > 0)
					{
					memoN = 0;
					for (a1=0; a1<4; a1++)
						for (a2=0; a2<4; a2++)
							L1tab[read][a1][a2] = Log10Memo (meanAllelicImbalance * pA1tab[a1][read] + (1.0 - meanAllelicImbalance) * pA2tab[a2][read], memoKeys, memoVals, &memoN);
					}

		memoN = 0;	/* memo of the final log10 (GL) of each genotype */
		for (a1=0; a1<4; a1++) /* a1 is maternal; a2 is paternal */
			for (a2=0; a2<4; a2++)
				{
				logGL[a1][a2] = 0.0;
				term1 = term2 = term3 = maxterm = 0.0;
				num1 = num2 = num3 = 0.0;
				
				for (read=0; read<4; read++)
					{
					pReadGivenA1 = pReadGivenA2 = 0.0;
					if (readCount[read] > 0)
						{
						pReadGivenA1 = pA1tab[a1][read];
						pReadGivenA2 = pA2tab[a2][read];
						
						/* M [p(b|A)] */
						if (alleleDropoutRate != 1)
							term1 += readCount[read] * L1tab[read][a1][a2];
						}
					//fprintf (stderr, "\n***read=%c count=%3d gl[%c][%c] = %lf  logGL[%c][%c] = %lf (p1=%e p2=%e), term1=%Lf term2=%Lf term3=%Lf", WhichNuc(read), readCount[read], WhichNuc(a1), WhichNuc(a2), GL[a1][a2],  WhichNuc(a1), WhichNuc(a2), logGL[a1][a2], pReadGivenA1,pReadGivenA2, term1, term2, term3);
					} // reads
				if (alleleDropoutRate != 0)
					{
					term2 = T2a[a1];
					term3 = T3a[a2];
					}
				
				/* compute log sum exp; first find maximum term */
				if (alleleDropoutRate == 0)
					maxterm = term1;
				else if (alleleDropoutRate == 1)
					{
					maxterm = term2;
					if (term3 > maxterm)
						maxterm = term3;
					}
				else
					{
					maxterm = term1;
					if (term2 > maxterm)
						maxterm = term2;
					if (term3 > maxterm)
						maxterm = term3;
					}
				
				if (alleleDropoutRate == 0)
					num1 = (1.0 - alleleDropoutRate) * Pow10Fast(term1-maxterm);
				else if (alleleDropoutRate == 1)
					{
					num2 = alleleDropoutRate * meanAllelicImbalance * Pow10Fast(term2-maxterm);
					num3 = alleleDropoutRate * (1.0 - meanAllelicImbalance) * Pow10Fast(term3-maxterm);;
					}
				else
					{
					num1 = (1.0 - alleleDropoutRate) * Pow10Fast(term1-maxterm);
					num2 = alleleDropoutRate * meanAllelicImbalance * Pow10Fast(term2-maxterm);
					num3 = alleleDropoutRate * (1.0 - meanAllelicImbalance) * Pow10Fast(term3-maxterm);;
					}
	
				GL[a1][a2] = num1 + num2 + num3;

				if (isnan(GL[a1][a2]) == NO)
					logGL[a1][a2] = maxterm + Log10Memo (GL[a1][a2], memoKeys, memoVals, &memoN);
				else
					logGL[a1][a2] = -INFINITY;
				/*
				fprintf (stderr, "\nterm1,2,3  = %Lf %Lf %Lf  | max = %Lf", term1, term2, term3, maxterm);
				fprintf (stderr, "\nnum1,2,3 = %Lf %Lf %Lf", num1, num2, num3);
				fprintf (stderr, "\nGL[%c][%c] = %lf", WhichNuc(a1), WhichNuc(a2), GL[a1][a2]);
				fprintf (stderr, "\nlogGL[%c][%c] = %lf\n", WhichNuc(a1), WhichNuc(a2), logGL[a1][a2]);
				*/
				} // genotypes
		
		if (debug_GL == YES)
			{
			fprintf (stderr,"\ncell %d site %d: %c%c", i+1, j+1, WhichNuc(maternalAllele), WhichNuc(paternalAllele));
			fprintf (stderr,"|\treads: A:%d C:%d G:%d T:%d  | total:%d",  readCount[A], readCount[C], readCount[G], readCount[T], numReads);
			fprintf (stderr,"\nmatAmpError = %e", maternalSiteAmplificationError);
			fprintf (stderr,"\npatAmpError = %e", paternalSiteAmplificationError);
			fprintf (stderr,"\nsequencing error = %f",sequencingError);
			fprintf (stderr,"\nallelic dropout rate = %f",alleleDropoutRate);
			for (a1=0; a1<4; a1++)
				for (a2=0; a2<4; a2++) //a2 is paternal
					fprintf (stderr, "\n(general)log10 GL[%c][%c] = %lf", WhichNuc(a1), WhichNuc(a2), logGL[a1][a2]);
			}
		} //model general


	/* obtain unphased GLs by combining GLs from both phases -only heterozygotes */
	for (a1=0; a1<4; a1++)
		for (a2=a1+1; a2<4; a2++)
			if (logGL[a1][a2] != logGL[a2][a1])
				{
				number1 = 0.5 * Pow10Fast(logGL[a1][a2]);
				number2 = 0.5 * Pow10Fast(logGL[a2][a1]);
				logGL[a1][a2] = logGL[a2][a1] = log10(number1 + number2);
				//fprintf (stderr, "\n1/2GL[%c][%c]=%Lf   1/2GL[%c][%c]=%Lf   logGL[%c][%c]=%f  logGL[%c][%c]=%f ", WhichNuc(a1), WhichNuc(a2), number1, WhichNuc(a2), WhichNuc(a1), number2, WhichNuc(a1), WhichNuc(a2), logGL[a1][a2], WhichNuc(a2), WhichNuc(a1), logGL[a2][a1] );
				}

	if (debug_GL == YES)
		{
		//PrintSiteInfo (stderr, SNVsites[snv]);
		fprintf (stderr, "\n");
		for (a1=0; a1<4; a1++)
			for (a2=0; a2<4; a2++)
				fprintf (stderr, "\nunphasedlog10 GL[%c][%c] = %lf", WhichNuc(a1), WhichNuc(a2), logGL[a1][a2]);
		fprintf (stderr, "\n");
		}
	
	/* find max log10 likelihood  */
	maxLike = logGL[A][A];
	MLmatAllele = MLpatAllele = A;
	for (a1=0; a1<4; a1++)
		for (a2=a1; a2<4; a2++)
			if (logGL[a1][a2] > maxLike)
				{
				maxLike = logGL[a1][a2];
				MLmatAllele = a1;
				MLpatAllele = a2;
				}
		
	/* normalize to log likelihood ratios */
	for (a1=0; a1<4; a1++)
		for (a2=a1; a2<4; a2++)
			{
			normLogGL[a1][a2] = logGL[a1][a2] - maxLike;
			normLogGL[a2][a1] = normLogGL[a1][a2];
			}

	//TODO: print GL=-999?
	/* substitute small numbers for printing */
	/*for (a1=0; a1<4; a1++)
		for (a2=a1; a2<4; a2++)
			if (logGL[a1][a2] < 1e-6 )
				{
				logGL[a2][a1] = logGL[a1][a2] = smallNumber;
				normLogGL[a2][a1] = normLogGL[a1][a2] = smallNumber;
				}
	*/
	
	if (debug_GL == YES)
		{
		//PrintSiteInfo (stderr, SNVsites[snv]);
		fprintf (stderr, "\nmaxlike = %lf\n", maxLike);
		for (a1=0; a1<4; a1++)
			for (a2=0; a2<4; a2++)
				fprintf (stderr, "\nunphased_norm_log10 GL[%c][%c] = %lf", WhichNuc(a1), WhichNuc(a2), normLogGL[a1][a2]);
		fprintf (stderr, "\n");
		if (debug_GL_EXIT == YES)
			exit(1);
		}
/*
	if (debug_GL == YES)
		{
		// Check if there are tied ML genotypes
		// We do expect these ties
		k = 0;
		for (a1=0; a1<4; a1++)
			for (a2=a1; a2<4; a2++)
				{
				if (normLogGL[a1][a2] == 0)
					k++;
				}
		if (k > 1)
			{
			MLmatAllele = MLpatAllele = MISSING;
			fprintf (stderr,"\n==>ERR: tied ML genotypes\n\n");
			exit(-1);
			}
		}
 */

	c->MLmatAllele = MLmatAllele;
	c->MLpatAllele = MLpatAllele;

	/* hand the likelihoods over in compact form (not stored when the caller does not want them) */
	if (glOut)
		{
		for (a1=0; a1<4; a1++)
			for (a2=a1; a2<4; a2++)
				glOut[GL_Index(a1,a2)] = logGL[a1][a2];
		glOut[GL_MAX] = maxLike;
		}


	/*
	if (debug_GL == YES)
		//exit(1);
	*/

}

/********************* MakeDoublets  ************************/
/*
 We simulate the possibility of amplifying/sequencing two cells at once
 by adding read counts from another sampled cell.
 Only tumor cells are considered.
*/

/* Read counts and likelihoods of the doublet made by cell c1 and donor cell c2 at site j; gl1 are the genotype likelihoods of c1 */
static void DoubletAtSite (int c1, int c2, int j, const double *gl1, double *probs, double **ngsEij, double **ampEijmat, double **ampEijpat, long int *seed)
	{
	int		k, l;
	int		debug_DL = NO;
	int		MLmatAllele, MLpatAllele;
	double	maxLike;
	long double number1, number2;
	CellSiteStr donor;
	double donorGL[GL_NUM];

	{
		memset (&donor, 0, sizeof(donor));
		memset (donorGL, 0, sizeof(donorGL));
		if (debug_DL == YES)
			fprintf (stderr,"\n\n-- New celldoublet reads and likelihoods");
		SiteReadCounts (&donor, c2, j, probs, ngsEij, ampEijmat, ampEijpat, seed);
		if (readsOnly == NO)
			GenotypeLikelihoods (&donor, c2, j, probs, ngsEij, ampEijmat, ampEijpat, seed, donorGL);

		if (debug_DL == YES)
			{
			fprintf (stderr,"\ncelltype %d site %d: %c%c", c2+1, j+1, WhichNuc(donor.maternalAllele), WhichNuc(donor.paternalAllele));
			fprintf (stderr,"|\treads: A:%d C:%d G:%d T:%d  | total:%d ",  donor.readCount[A], donor.readCount[C], donor.readCount[G], donor.readCount[T], donor.numReads);
			fprintf (stderr,"\nmatAmpError = %f", donor.maternalSiteAmplificationError);
			fprintf (stderr,"\npatAmpError = %f", donor.paternalSiteAmplificationError);
			for (k=0; k<4; k++)
				for (l=0; l<4; l++)
					fprintf (stderr, "\nlog10 GL[%c][%c] = %lf", WhichNuc(k), WhichNuc(l), GLV(donorGL,k,l));
			fprintf (stderr, "\n");
			}
		
		/* add doublet information */
		CD(c1,j).numReadsDoublet = CS(c1,j).numReads + donor.numReads;
		for (k=0; k<4; k++)
			CD(c1,j).readCountDoublet[k] = CS(c1,j).readCount[k] + donor.readCount[k];

		CD(c1,j).maternalAlleleDoublet = donor.maternalAllele;
		CD(c1,j).paternalAlleleDoublet = donor.paternalAllele;
		if (readsOnly == YES)
			return;

		//fprintf(stderr,"\nln cell\t\tcelltype\tdoublet");
		/* adding the genotype log likelihoods of two sets of reads (c1 + celldoublet) to get the total doublet likelihoods */
		for (k=0; k<4; k++)
			for (l=0; l<4; l++)
				{
				if (isinf(GLV(gl1,k,l)) || isinf(GLV(donorGL,k,l)))
					CD(c1,j).logGLdoublet[k][l] = -INFINITY;
				else
					CD(c1,j).logGLdoublet[k][l] = GLV(gl1,k,l) + GLV(donorGL,k,l);
				}
		
		if (debug_DL == YES)
			{
			fprintf (stderr,"\nCell %d site %d: %c%c", c1+1, j+1, WhichNuc(CS(c1,j).maternalAllele), WhichNuc(CS(c1,j).paternalAllele));
			fprintf (stderr,"|\treads: A:%d C:%d G:%d T:%d  | total:%d ",   CS(c1,j).readCount[A],  CS(c1,j).readCount[C],  CS(c1,j).readCount[G],  CS(c1,j).readCount[T],  CS(c1,j).numReads);
			fprintf (stderr,"\n|reads doublet: A:%d C:%d G:%d T:%d  | total:%d ",   CD(c1,j).readCountDoublet[A],  CD(c1,j).readCountDoublet[C],  CD(c1,j).readCountDoublet[G],  CD(c1,j).readCountDoublet[T],  CD(c1,j).numReadsDoublet);
			fprintf (stderr,"\nmatAmpError = %f",  CS(c1,j).maternalSiteAmplificationError);
			fprintf (stderr,"\npatAmpError = %f",  CS(c1,j).paternalSiteAmplificationError);
			for (k=0; k<4; k++)
				for (l=0; l<4; l++)
					fprintf (stderr, "\ndoublet log10 GL[%c][%c] = %lf", WhichNuc(k), WhichNuc(l),  CD(c1,j).logGLdoublet[k][l]);
			fprintf (stderr, "\n");
			}

		/* obtain unphased GLs by combining GLs from both phases -only heterozygotes */
		   for (k=0; k<4; k++)
			   for (l=k+1; l<4; l++)
				   if (CD(c1,j).logGLdoublet[k][l] != CD(c1,j).logGLdoublet[l][k])
					   {
					   number1 = 0.5 * pow(10.0,CD(c1,j).logGLdoublet[k][l]);
					   number2 = 0.5 * pow(10.0,CD(c1,j).logGLdoublet[l][k]);
					   CD(c1,j).logGLdoublet[k][l] = CD(c1,j).logGLdoublet[l][k] = log10(number1 + number2);
					   //fprintf (stderr, "\n1/2GL[%c][%c]=%Lf   1/2GL[%c][%c]=%Lf   logGL[%c][%c]=%f  logGL[%c][%c]=%f ", WhichNuc(a1), WhichNuc(a2), number1, WhichNuc(a2), WhichNuc(a1), number2, WhichNuc(a1), WhichNuc(a2), logGL[a1][a2], WhichNuc(a2), WhichNuc(a1), logGL[a2][a1] );
					   }

		/* find max log10 doublet likelihood  */
		maxLike = CD(c1,j).logGLdoublet[A][A];
		MLmatAllele = MLpatAllele = A;
		for (k=0; k<4; k++)
			for (l=k; l<4; l++)
				if (CD(c1,j).logGLdoublet[k][l] > maxLike)
					{
					maxLike = CD(c1,j).logGLdoublet[k][l];
					MLmatAllele = k;
					MLpatAllele = l;
					}
		
		/* normalize to likelihood ratios */
		for (k=0; k<4; k++)
			for (l=k; l<4; l++)
			{
			CD(c1,j).normLogGLdoublet[k][l] = CD(c1,j).logGLdoublet[k][l] - maxLike;
			CD(c1,j).normLogGLdoublet[l][k] = CD(c1,j).normLogGLdoublet[k][l];
			}

		//TODO: print GL=-999?
		/* substitute small numbers for printing */
		/*for (k=0; k<4; k++)
			for (l=k; l<4; l++)
				if (CD(c1,j).logGLdoublet[k][l] < 1e-6 )
					{
					CD(c1,j).logGLdoublet[l][k] = CD(c1,j).logGLdoublet[k][l] = smallNumber;
					CD(c1,j).normLogGLdoublet[l][k] = CD(c1,j).normLogGLdoublet[k][l] = smallNumber;
					}
		 */
		CD(c1,j).MLmatAlleleDoublet = MLmatAllele;
		CD(c1,j).MLpatAlleleDoublet = MLpatAllele;
			
		if (debug_DL == YES)
			{
			for (k=0; k<4; k++)
				for (l=0; l<4; l++)
					fprintf (stderr, "\nscaled doublet log10 GL[%c][%c] = %lf", WhichNuc(k), WhichNuc(l), CD(c1,j).normLogGLdoublet[k][l]);
			fprintf (stderr, "\n");
			}
		}
	}


