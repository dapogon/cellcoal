/*
	CellCoal - allelic dropout, genotyping errors and counting of alleles
	See cellcoal.h for the shared declarations. This file was split from cellcoal.c.
*/
#include "cellcoal.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "definitions.h"

/********************* AllelicDropout  ************************/
/*
 Remove alleles from single chromosomes at a given ADO rate per genotype

 We assume that ADO is the product of the allele dropout rate as follows:
 A = genotype ADO; a = allele ADO
 A = 2a - a^2
 a = 1 - sqrt(1-A)
*/

void AllelicDropout (long int *seed)
	{
	int i,j;
	double alleleADOrateMean = 0.0, *ADOrateCell = NULL, *ADOrateSite = NULL;
	double alleleRate, cellSiteRate;

	/* The rate of each genotype is computed on the fly from the rate of its cell and the rate of its site,
	   so that no cells x sites matrix of rates has to be stored. */
	if (doADOcell == NO && doADOsite == NO)
		{
		// fixed ADO rate
		alleleADOrateMean = 1.0 - sqrt (1.0 - fixedADOrate);
		}
	else // variable ADO rates
		{
		/* The rate of each cell and of each site must be sampled once and stored, so that every
		   genotype uses the rate of its own cell and its own site (previously only the last
		   sampled value of each was kept and applied to all cells and sites). The sampling
		   order (cells first, then sites) is unchanged. */
		ADOrateCell = (double*) calloc (numCells+1, sizeof(double));
		ADOrateSite = (double*) calloc (numSites, sizeof(double));
		if (!ADOrateCell || !ADOrateSite)
			{
			fprintf (stderr, "Could not allocate the ADO cell/site rate structures\n");
			exit (-1);
			}
		for (i=0; i<numCells+1; i++)
			ADOrateCell[i] = (doADOcell == YES) ? RandomBetaMeanVar(meanADOcell, varADOcell, seed) : 0;
		if (doADOcell == YES)
			{
			WriteADOCells (ADOrateCell, adoReplicate);
			KeepCellADO (ADOrateCell);
			}
		for (j=0; j<numSites; j++)
			ADOrateSite[j] = (doADOsite == YES) ? RandomBetaMeanVar(meanADOsite, varADOsite, seed) : 0;
		}

	for (i=0; i<numCells+1; i++)
		{
		for (j=0; j<numSites; j++)
			{
			if (ADOrateCell == NULL)
				alleleRate = alleleADOrateMean;
			else
				{
				cellSiteRate = ADOrateCell[i] + ADOrateSite[j] - (ADOrateCell[i] * ADOrateSite[j]);
				alleleRate = 1.0 - sqrt (1.0 - cellSiteRate); // at the allele level
				}
			if (GetG(MATERNAL, i, j) != DELETION)
				{
				if (RandomUniform(seed) < alleleRate)
					{
					SetG(MATERNAL, i, j, ADO);
					allSites[j].hasADO = YES;
					}
				}
			if (GetG(PATERNAL, i, j) != DELETION)
				{
				if (RandomUniform(seed) < alleleRate)
					{
					SetG(PATERNAL, i, j, ADO);
					allSites[j].hasADO = YES;
					}
				}
			}
		}

	free (ADOrateCell);
	free (ADOrateSite);
	}


/*
void AllelicDropoutGamma (long int *seed)
	{
	int i,j;
	double alleleADOrateMean, alleleADOrateMod;
	double **alleleADOrate;
	
	// calculate ADO rate for a given allele
	alleleADOrateMean = 1.0 - sqrt (1.0 -  ADOrate);

	// allocate space for the ADO rates
	alleleADOrate = (double**) calloc (numCells+1, sizeof(double*));
	if (!alleleADOrate)
		{
		fprintf (stderr, "Could not allocate the alleleADOrate structure\n");
		exit (-1);
		}
	for (i=0; i<numCells+1; i++)
		{
		alleleADOrate[i] = (double*) calloc (numSites, sizeof(double));
		if (!alleleADOrate[i])
			{
			fprintf (stderr, "Could not allocate the alleleADOrate[i] structure\n");
			exit (-1);
			}
		}

	// sample rate modifiers
	for (i=0; i<numCells+1; i++)
		for (j=0; j<numSites; j++)
			alleleADOrate[i][j] = alleleADOrateMean;

	if (ADOvarAmongCells == YES)
		for (i=0; i<numCells+1; i++)
			{
			alleleADOrateMod = RandomGamma (alphaADOcells, seed) / alphaADOcells;
			for (j=0; j<numSites; j++)
				alleleADOrate[i][j] *= alleleADOrateMod;
			}
	if (ADOvarAmongSites == YES)
		for (j=0; j<numSites; j++)
			{
			alleleADOrateMod = RandomGamma (alphaADOsites, seed) / alphaADOsites;
			for (i=0; i<numCells+1; i++)
				alleleADOrate[i][j] *= alleleADOrateMod;
			}

//	for (i=0; i<numCells+1; i++)
//		for (j=0; j<numSites; j++)
//			fprintf (stderr, "\ncell %d  site %d  alleleADOrate[i][j] = %6.4f", i+1, j+1, alleleADOrate[i][j]);

	
	for (i=0; i<numCells+1; i++)
		{
		for (j=0; j<numSites; j++)
			{
			if (GetG(MATERNAL, i, j) != DELETION)
				{
				if (RandomUniform(seed) < alleleADOrate[i][j])
					{
					SetG(MATERNAL, i, j, ADO);
					allSites[j].hasADO = YES;
					}
				}
			if (GetG(PATERNAL, i, j) != DELETION)
				{
				if (RandomUniform(seed) < alleleADOrate[i][j])
					{
					SetG(PATERNAL, i, j, ADO);
					allSites[j].hasADO = YES;
					}
				}
			}
		}
	
	// free allele ADO rates
	for (i=0; i<numCells+1; i++)
		free(alleleADOrate[i]);
	free(alleleADOrate);
	}
*/

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
	int 	i,j;
	double	genotypingError, alleleError, normEij[4][4], rowSum;
	
	/* ChooseUniformState expects probabilities that sum to 1, but the rows of Eij are
	   relative weights (e.g. 0 1 1 1 by default), so normalize them first. Without this
	   the first non-zero entry always wins and the wrong base is not drawn at random. */
	for (i=0; i<4; i++)
		{
		rowSum = Eij[i][0] + Eij[i][1] + Eij[i][2] + Eij[i][3];
		for (j=0; j<4; j++)
			normEij[i][j] = (rowSum > 0.0) ? Eij[i][j] / rowSum : 0.0;
		}
	
	if (alphabet == DNA)
		{
		for (i=0; i<numCells+1; i++)
			{
			for (j=0; j<numSites; j++)
				{
				genotypingError = RandomBetaMeanVar(meanGenotypingError, varGenotypingError, seed);
				alleleError = 1.0 - sqrt (1.0 -  genotypingError);

				if (GetG(MATERNAL, i, j) != ADO && GetG(MATERNAL, i, j) != DELETION && RandomUniform(seed) < alleleError)
					{
					SetG(MATERNAL, i, j, ChooseUniformState (normEij[GetG(MATERNAL, i, j)], seed));
					allSites[j].hasGenotypeError = YES;
					}
				if (GetG(PATERNAL, i, j) != ADO && GetG(PATERNAL, i, j) != DELETION && RandomUniform(seed) < alleleError)
					{
					SetG(PATERNAL, i, j, ChooseUniformState (normEij[GetG(PATERNAL, i, j)], seed));
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
				genotypingError = RandomBetaMeanVar(meanGenotypingError, varGenotypingError, seed);
				alleleError = 1.0 - sqrt (1.0 -  genotypingError);
				if (GetG(MATERNAL, i, j) != ADO && GetG(MATERNAL, i, j) != DELETION && RandomUniform(seed) < alleleError)
					{
					if (GetG(MATERNAL, i, j) == 0)
						SetG(MATERNAL, i, j, 1);
					else if (GetG(MATERNAL, i, j) == 1)
						SetG(MATERNAL, i, j, 0);
					allSites[j].hasGenotypeError = YES;
					}
				if (GetG(PATERNAL, i, j) != ADO && GetG(PATERNAL, i, j) != DELETION && RandomUniform(seed) < alleleError)
					{
					if (GetG(PATERNAL, i, j) == 0)
						SetG(PATERNAL, i, j, 1);
					else if (GetG(PATERNAL, i, j) == 1)
						SetG(PATERNAL, i, j, 0);
					allSites[j].hasGenotypeError = YES;
					}
				}
			}
		}
	}



/************************* CountTrueVariants  ************************/
/* Count number of variants in a given data set, assuming there is no ADO (yet).
   A variant is defined as no-deletion genotype different from the healthy root genotype */

int CountTrueVariants (void)
	{
    int		cell, site;
    int		nVariants = 0;
    char	*isVar;
    Geno	*matRoot, *patRoot, *mat, *pat;

    if (sparseMode == YES)
        {
        /* only the sites with a column can differ from the root */
        for (site=0; site<numSites; site++)
            {
            Geno *c = sparseCols[site];
            Geno mr, pr;
            if (c == NULL)
                continue;
            mr = data[MATERNAL][OUTGROUP_ROOT][site];
            pr = data[PATERNAL][OUTGROUP_ROOT][site];
            for (cell=0; cell<numCells+1; cell++)
                {
                if ((c[cell] != DELETION && c[cell] != mr) || (c[numNodes+cell] != DELETION && c[numNodes+cell] != pr))
                    {
                    allSites[site].isVariant = YES;
                    variantSites[nVariants++] = site;
                    break;
                    }
                }
            }
        return nVariants;
        }

    /* cell-major traversal (genomes are stored contiguously per cell) */
    isVar = (char *) calloc (numSites, sizeof(char));
    if (!isVar)
        {
        fprintf (stderr, "Could not allocate the isVar vector\n");
        exit (-1);
        }
    matRoot = data[MATERNAL][OUTGROUP_ROOT];
    patRoot = data[PATERNAL][OUTGROUP_ROOT];
    for (cell=0; cell<numCells+1; cell++)
        {
        mat = data[MATERNAL][cell];
        pat = data[PATERNAL][cell];
        for (site=0; site<numSites; site++)
            {
            if (isVar[site] == 0 &&
                ((mat[site] != DELETION && mat[site] != matRoot[site]) ||
                 (pat[site] != DELETION && pat[site] != patRoot[site])))
                isVar[site] = 1;
            }
        }
    for (site=0; site<numSites; site++)
        {
        if (isVar[site])
            {
            allSites[site].isVariant = YES;
            variantSites[nVariants++] = site;
            }
        }
    free (isVar);
    return nVariants;
	}


/************************* CountAllelesInObservedGenotypes  ************************/
/* Identify reference and alternate alleles plus SNVs, in ingroup plus outgroup cell genotypes
 	after introducing mutations, cnLOH, deletions, ADO and genotype errors
   A SNV is defined as no-deletion genotype different from the reference genotype */

int CountAllelesInObservedGenotypes (void)
	{
	int		numAltAlleles, cell, site, s0, s1, k, ch;
	int		countA, countC, countG, countT, countADO, countDEL, nSNVs;
	int		*cnt[7];
	Geno	*row;
	unsigned char slot[19];
	enum { BLK = 4096 };
	
	nSNVs = 0;
	
	/* counts for a block of sites are accumulated cell by cell (genomes are stored contiguously per cell) */
	/* slot of each genotype value (offset by 9 to allow MISSING = -9); slot 6 collects everything else */
	for (k=0; k<19; k++)
		slot[k] = 6;
	slot[A+9] = 0; slot[C+9] = 1; slot[G+9] = 2; slot[T+9] = 3; slot[ADO+9] = 4; slot[DELETION+9] = 5;
	for (k=0; k<7; k++)
		{
		cnt[k] = (int *) calloc (BLK, sizeof(int));
		if (!cnt[k])
			{
			fprintf (stderr, "Could not allocate the allele counts vector\n");
			exit (-1);
			}
		}

	for (s0=0; s0<numSites; s0+=BLK)
		{
		s1 = (s0 + BLK < numSites) ? s0 + BLK : numSites;
		for (k=0; k<7; k++)
			memset (cnt[k], 0, BLK * sizeof(int));
		if (sparseMode == YES)
			{
			for (site=s0; site<s1; site++)
				{
				Geno *c = sparseCols[site];
				if (c == NULL && !(errCount != NULL && errCount[site] > 0))
					{
					/* all the cells are identical to the root (ADO at these sites is counted by PrintSiteInfo, when needed) */
					cnt[slot[data[MATERNAL][OUTGROUP_ROOT][site]+9]][site-s0] += numCells+1;
					cnt[slot[data[PATERNAL][OUTGROUP_ROOT][site]+9]][site-s0] += numCells+1;
					}
				else
					{
					for (cell=0; cell<numCells+1; cell++)
						{
						cnt[slot[GetG(MATERNAL, cell, site)+9]][site-s0]++;
						cnt[slot[GetG(PATERNAL, cell, site)+9]][site-s0]++;
						}
					}
				}
			}
		else
		for (cell=0; cell<numCells+1; cell++)
			{
			for (ch=0; ch<2; ch++)
				{
				row = data[ch == 0 ? MATERNAL : PATERNAL][cell];
				for (site=s0; site<s1; site++)
					{
					cnt[slot[row[site]+9]][site-s0]++;
					}
				}
			}

		for (site=s0; site<s1; site++)
			{
			countA = cnt[0][site-s0];
			countC = cnt[1][site-s0];
			countG = cnt[2][site-s0];
			countT = cnt[3][site-s0];
			countADO = cnt[4][site-s0];
			countDEL = cnt[5][site-s0];

		allSites[site].countA = countA;
		allSites[site].countC = countC;
		allSites[site].countG = countG;
		allSites[site].countT = countT;
		allSites[site].countACGT = countA + countC + countG + countT;
		allSites[site].countDropped = countADO + countDEL;

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
		}

	for (k=0; k<7; k++)
		free (cnt[k]);
	return nSNVs;
	}


#ifdef COUNT_ML_GENOTYPE_ERRORS

/************************* CountMLGenotypingErrors  ************************/
/*
 Counts how many SNV genotypes are called wrongly
*/
void CountMLGenotypingErrors()
	{
	int			i, j, equalGenotypes, countCalledGenotypes, countMLgenotypeErrors;
	int			trueMaternalAllele, truePaternalAllele, MLmatAllele, MLpatAllele;
	
	CellSiteStr *c;

	countMLgenotypeErrors = 0;
	countCalledGenotypes = 0;
	for (i=0; i<numCells+1; i++)
		{
	//	for (snv=0; snv<numSNVs; snv++) // otherwise, with low coverage, increasing amplification error increases accuracy, as it results in many SNVS where most cells do not contain errors...
			//{
			//j = SNVsites[snv];*/
		for (j=0; j<numSites; j++)
			{
			c = &CS(i,j);
			trueMaternalAllele = c->trueMaternalAllele;
			truePaternalAllele = c->truePaternalAllele;
			if (cell[i].hasDoublet == NO)
				{
				MLmatAllele = c->MLmatAllele;
				MLpatAllele = c->MLpatAllele;
				}
			else
				{
				MLmatAllele = CD(i,j).MLmatAlleleDoublet;
				MLpatAllele = CD(i,j).MLpatAlleleDoublet;
				}
			
			if (MLmatAllele == MISSING) // no reads for this genotype
				equalGenotypes = YES;
			else
				{
				countCalledGenotypes++;
				equalGenotypes = CompareGenotypes (trueMaternalAllele, truePaternalAllele, MLmatAllele, MLpatAllele);
				}
			
			if (equalGenotypes == NO)
				countMLgenotypeErrors++;
			}
		}

	fprintf(fpSims, "%d, %10.8f, %d, %10.8f, %d\n", countCalledGenotypes, countCalledGenotypes / ((numCells+1) * (double)numSites), countMLgenotypeErrors, countMLgenotypeErrors/(double)countCalledGenotypes, numSNVs);
	cumCountCalledGenotypes += countCalledGenotypes / ((numCells+1) * (double)numSites);
	cumCountCalledGenotypesSq += pow (countCalledGenotypes / ((numCells+1) * (double)numSites), 2);
	cumCountMLgenotypeErrors += countMLgenotypeErrors/(double)countCalledGenotypes;
	cumCountMLgenotypeErrorsSq += pow (countMLgenotypeErrors/(double)countCalledGenotypes, 2);
	}


/********************* CompareGenotypes  ************************/
/*
 Compares two unphased genotypes (i.e, A/T = T/A)
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

#endif


/************************* CountAllelesInMLGenotypes  ************************/
/*
 Used when NGS reads are generated and ML genotypes called.
 A SNV is defined as an ML genotype different from the reference genotype
 Note that ML genotypes cannot be ADO or DELETION
 */
/* Counts the ML alleles of the cells at site j (before doublets), and returns YES if the site is an apparent SNV */
int MLCountsAtSite (int j)
	{
	int		numAltAllelesInCalls, i;
	int		countA, countC, countG, countT, countN;

	countA = countC = countG = countT = countN = 0;
	for (i=0; i<numCells+1; i++)
		{
		int m = CS(i,j).MLmatAllele, q = CS(i,j).MLpatAllele;
		if (m == A) countA++; else if (m == C) countC++; else if (m == G) countG++; else if (m == T) countT++; else countN++;
		if (q == A) countA++; else if (q == C) countC++; else if (q == G) countG++; else if (q == T) countT++; else countN++;
		}

	/* keep counts on ML genotypes */
	allSites[j].countA = countA;
	allSites[j].countC = countC;
	allSites[j].countG = countG;
	allSites[j].countT = countT;
	allSites[j].countN = countN;
	allSites[j].countACGT = countA + countC + countG + countT;
	allSites[j].countCellswithData = allSites[j].countACGT / 2.0;

	/* count number of alternate alleles in reads that appear in the ML genotypes */
	numAltAllelesInCalls = 0;
	if (countA > 0 && allSites[j].referenceAllele != A)
		numAltAllelesInCalls++;
	if (countC > 0 && allSites[j].referenceAllele != C)
		numAltAllelesInCalls++;
	if (countG > 0 && allSites[j].referenceAllele != G)
		numAltAllelesInCalls++;
	if (countT > 0 && allSites[j].referenceAllele != T)
		numAltAllelesInCalls++;

	/* find out whether this site is a SNV */
	if (numAltAllelesInCalls > 0)
		{
		allSites[j].isSNV = YES;
		return YES;
		}
	return NO;
	}

int CountAllelesInMLGenotypes(void)
	{
	int		j, nSNVs;

	nSNVs = 0;
	for (j=0; j<numSites; j++)
		if (MLCountsAtSite (j) == YES)
			SNVsites[nSNVs++] = j;
	return nSNVs;
}



