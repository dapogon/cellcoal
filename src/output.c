/*
	CellCoal - result files and the printing of genotypes, haplotypes, trees, times and information
	See cellcoal.h for the shared declarations. This file was split from cellcoal.c.
*/
#include "cellcoal.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "definitions.h"

/* The files of the results: a file is written when its flag is on (and `and`, if there is one, too). The first column of the
   file is the number of replicates (header = 1); the VCF is opened by OpenVCF. */
typedef struct
	{
	const int	*on, *and;
	FILE		**fp;
	const char	*dir, *file;
	int			header, vcf;
	} OutputFile;

static const OutputFile OUTPUT_FILES[] =
	{
	{&doPrintTree,			NULL,		&fpTrees,			treeDir,			treeFile,			0, 0},
	{&doPrintTimes,			NULL,		&fpTimes,			timesDir,			timesFile,			0, 0},
	{&doPrintSNVgenotypes,	NULL,		&fpSNVgenotypes,	SNVgenotypesDir,	SNVgenotypesFile,	1, 0},
	{&doPrintSNVhaplotypes,	NULL,		&fpSNVhaplotypes,	SNVhaplotypesDir,	SNVhaplotypesFile,	1, 0},
	{&doPrintTrueHaplotypes,NULL,		&fpTrueHaplotypes,	trueHaplotypesDir,	trueHaplotypesFile,	1, 0},
	{&doPrintFullGenotypes,	NULL,		&fpFullGenotypes,	fullGenotypesDir,	fullGenotypesFile,	1, 0},
	{&doPrintFullHaplotypes,NULL,		&fpFullHaplotypes,	fullHaplotypesDir,	fullHaplotypesFile,	1, 0},
	{&doPrintVCF,			&doNGS,		&fpVCF,				VCFdir,				VCFfile,			1, 1},
	{&doPrintCATG,			NULL,		&fpCATG,			CATGdir,			CATGfile,			1, 0},
	{&doPrintPILEUP,		NULL,		&fpPILEUP,			PILEUPdir,			PILEUPfile,			1, 0},
	{&doPrintBulk,			NULL,		&fpBulk,			bulkDir,			bulkFile,			1, 0},
	{NULL, NULL, NULL, NULL, NULL, 0, 0}
	};

static int OutputWanted (const OutputFile *o)
	{
	return *o->on == YES && (o->and == NULL || *o->and == YES);
	}

/* opens the files of the results: all the replicates in one file (replicate < 0) or the file of one replicate in its folder */
static void OpenOutputFiles (int replicate)
	{
	const OutputFile	*o;
	char				path[4*MAX_NAME];

	for (o=OUTPUT_FILES; o->on != NULL; o++)
		{
		if (!OutputWanted (o))
			continue;
		if (replicate < 0)
			snprintf (path, sizeof(path), "%s/%s", resultsDir, o->file);
		else
			{
			snprintf (path, sizeof(path), "%s/%s", resultsDir, o->dir);
			MAKE_DIR (path);
			snprintf (path, sizeof(path), "%s/%s/%s.%04d", resultsDir, o->dir, o->file, replicate + 1);
			}
		if ((*o->fp = o->vcf ? OpenVCF (path) : fopen (path, "w")) == NULL)
			{
			fprintf (stderr, "Can't open \"%s\"\n", path);
			exit (-1);
			}
		}
	}

/********************* PrepareGlobalFiles **********************/
/* Open global files to output results */
void PrepareGlobalFiles (void)
	{
	const OutputFile	*o;

	OpenOutputFiles (-1);
	for (o=OUTPUT_FILES; o->on != NULL; o++)
		if (OutputWanted (o) && o->header && !(o->vcf && blocksOn))	/* a chromosome in blocks is one data set: a plain VCF */
			fprintf (*o->fp, "%d\n", numDataSets);
	if (doPrintBulk == YES)
		fprintf (fpBulk, "site\tref\treadsA\treadsC\treadsG\treadsT\tdepth\talt_reads\tbulk_vaf\ttrue_alt_alleles\ttrue_vaf\n");
	}

/********************* PrepareSeparateFiles **********************/
/* Open individual files to output results */
void PrepareSeparateFiles (int replicate)
	{
	OpenOutputFiles (replicate);
	if (doPrintBulk == YES)
		fprintf (fpBulk, "site\tref\treadsA\treadsC\treadsG\treadsT\tdepth\talt_reads\tbulk_vaf\ttrue_alt_alleles\ttrue_vaf\n");
	}


/**************** Print Trees ***************/
/*  Print trees to treefile in Newick format */

void PrintTree (TreeNode *treeRoot, FILE *fp)
    {
    if (numOutgroup == 0 && treeRoot == healthyRoot && coalTreeMRCA != NULL)
        {
        /* tumor only: the ingroup tree, as -T reads it (without the trunk and the outgroup) */
        fprintf (fp, "(");
        WriteTree (coalTreeMRCA->left, fp);
        fprintf (fp, ",");
        WriteTree (coalTreeMRCA->right, fp);
        fprintf (fp,");\n");
        return;
        }
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
           // fprintf (fp, "healthycell:%.10g",p->branchLength);
			fprintf (fp, "%s:%.10g",outCellName, p->branchLength);
			}
        else if (p->left == NULL && p->right == NULL)
            {
			if (doUserTree == NO)
				//fprintf (fp, "tumcell%04d:%.10g", Index(p),p->branchLength);
				fprintf (fp, "%s%04d:%.10g", inCellName, Index(p),p->branchLength);
			else
				fprintf (fp, "%s:%.10g", p->name,p->branchLength);
            }
        else
            {
            fprintf (fp, "(");
            WriteTree (p->left, fp);
            fprintf (fp, ",");
            WriteTree (p->right, fp);
            if (p->anc !=NULL)
				{
				fprintf (fp, "):%.10g",p->branchLength);
				//TODO: anotate tree?
				/*if (p->rateSwitch > 0)
					fprintf (fp, "[&rate_switch=%e]",p->rateSwitch);*/
				}
            }
        }
    }


/********************* PrintTimes **********************/
/*	Prints to a given file a detailed description of
	the tree as lost of nodes, times, branch lengths.
	Note that only the branch lengths can be modified
	after the coalescent (by rate variation among lineages
	or by tree rate switches)
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
                     "outgroupTip", Label(p), Index(p), Index(p->left), Index(p->right), Index(p->anc), p->time, p->length, p->branchLength);
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



/* the label of the row of a cell (haplotype 'm' or 'p'): the name of the cell and the spaces that follow it; LABEL_TIGHT does not
   put a space after the name of a cell of a user tree */
enum { LABEL_TIGHT = 0, LABEL_SPACED = 1 };

static void PrintRowLabel (FILE *fp, int i, int kind)
	{
	if (kind == 'm' || kind == 'p')
		{
		if (i == numCells)
			fprintf (fp, "%s%c  ", outCellName, kind);
		else if (doUserTree == NO)
			fprintf (fp, "%s%04d%c  ", inCellName, i+1, kind);
		else
			fprintf (fp, "%c%-*s ", kind, stringPrecision, cellNames[i]);
		}
	else if (i == numCells)
		fprintf (fp, "%s ", outCellName);
	else if (doUserTree == NO)
		fprintf (fp, "%s%04d ", inCellName, i+1);
	else
		fprintf (fp, kind == LABEL_SPACED ? "%-*s " : "%-*s", stringPrecision, cellNames[i]);
	}

/* Genotype of one cell at one site with extra copies (copy-number gains): the alleles of the maternal copies, '|', and those of the
   paternal copies; the first copy of each haplotype comes first. Without gains it is the usual two letters. */
static void PrintCopyGenotype (FILE *fp, int i, int site, int binary)
	{
	Geno	ex[MAX_COPIES];
	int		n, e;
	Geno	a[2] = {GetG(MATERNAL, i, site), GetG(PATERNAL, i, site)};
	int		h;

	for (h=0; h<2; h++)
		{
		if (h == 1 && CN_EVENTS && doPhase)
			fputc ('|', fp);
		fputc (binary ? WhichMut (a[h]) : WhichNuc (a[h]), fp);
		if (copyGainOn == YES && a[h] != DELETION)
			{
			n = ExtraCopies (i, h, site, ex);
			for (e=0; e<n; e++)
				fputc (binary ? WhichMut (ex[e]) : WhichNuc (ex[e]), fp);
			}
		}
	}

/************************ PrintSiteInfo ***********************/
/* Prints info for the site struct */

void PrintSiteInfo (FILE *fp, int i)
	{
	int j;

	/* In sparse mode the ADO of the sites that do not vary is not counted while the alleles are counted (it would cost
	   cells x sites operations); it is computed here, only for the sites that are printed. */
	if (sparseMode == YES && sparseADO == YES)
		{
		int	cc, dropped = 0;

		for (cc=0; cc<NUM_TAXA; cc++)
			{
			Geno gm = GetG(MATERNAL, cc, i), gp = GetG(PATERNAL, cc, i);

			if (gm == ADO || gm == DELETION)
				dropped++;
			if (gp == ADO || gp == DELETION)
				dropped++;
			if (gm == ADO || gp == ADO)
				allSites[i].hasADO = YES;
			}
		allSites[i].countDropped = dropped;
		}
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
	fprintf (fp, "\n numAltAlleles (in read counts) = %d   (", allSites[i].numAltAlleles);
	for (j=0; j <allSites[i].numAltAlleles; j++)
		fprintf (fp, " %c", WhichNuc(allSites[i].alternateAlleles[j]));
	fprintf (fp, " )");
	fprintf (fp, "\n alleleA = %d", allSites[i].countA);
	fprintf (fp, "\n alleleC = %d", allSites[i].countC);
	fprintf (fp, "\n alleleG = %d", allSites[i].countG);
	fprintf (fp, "\n alleleT = %d", allSites[i].countT);
	fprintf (fp, "\n alleleN = %d", allSites[i].countN);
	fprintf (fp, "\n countDropped = %d", allSites[i].countDropped);
	fprintf (fp, "\n readCountA = %d", allSites[i].readCountA);
	fprintf (fp, "\n readCountC = %d", allSites[i].readCountC);
	fprintf (fp, "\n readCountG = %d", allSites[i].readCountG);
	fprintf (fp, "\n readCountT = %d", allSites[i].readCountT);
	fprintf (fp, "\n rateMultiplier = %f", allSites[i].rateMultiplier);
	}

/************************ PrintSNVGenotypes ***********************/
/* Prints observed/ML genotypes at variable sites (SNVs) to a file */
void PrintSNVGenotypes (FILE *fp)
    {
    int		i, j;

    if (doPrintAncestors == YES && doNGS == NO)
        {
        fprintf (fp, "%d %d\n", NUM_TAXA+2, numSNVs);
        }
    else
        {
        fprintf (fp, "%d %d\n", NUM_TAXA, numSNVs);
        }

    /* site information */
	fprintf (fp, "%d", SNVsites[0]+1);
	for (i=1; i<numSNVs; i++)
		fprintf (fp, " %d", SNVsites[i]+1);
    fprintf (fp, "\n");
	
    if (alphabet == DNA)
        {
        if (doNGS == NO && !CN_EVENTS)
            {
            /* Genotypes are gathered for a block of cells, site by site (the genomes are stored by site in sparse
               mode), in a text buffer that is then written row by row. */
            enum { CELLBLOCK = 128 };
            long	rowLen = 3L * numSNVs;
            char	*buf = (char *) malloc ((size_t) CELLBLOCK * rowLen + 1);
            int		i0, i1;

            if (!buf)
                {
                fprintf (stderr, "Could not allocate the genotype print buffer\n");
                exit (-1);
                }
            for (i0=0; i0<NUM_TAXA; i0+=CELLBLOCK)
                {
                i1 = (i0 + CELLBLOCK < NUM_TAXA) ? i0 + CELLBLOCK : NUM_TAXA;
                for (j=0; j<numSNVs; j++)
                    {
                    int site = SNVsites[j];
                    for (i=i0; i<i1; i++)
                        {
                        char *q = buf + (long long)(i-i0) * rowLen + 3L * j;
                        q[0] = ' ';
                        q[1] = WhichNuc(GetG(MATERNAL, i, site));
                        q[2] = WhichNuc(GetG(PATERNAL, i, site));
                        }
                    }
                for (i=i0; i<i1; i++)
                    {
                    PrintRowLabel (fp, i, LABEL_TIGHT);fwrite (buf + (long long)(i-i0) * rowLen, 1, (size_t) rowLen, fp);
                    fprintf (fp,"\n");
                    }
                }
            free (buf);
            }
        else
        for (i=0; i<NUM_TAXA; i++)
            {
            PrintRowLabel (fp, i, LABEL_TIGHT);
			for (j=0; j<numSNVs; j++)
				{
				int	mm, pp;

				if (doNGS == NO && CN_EVENTS)
					{
					fputc (' ', fp);
					PrintCopyGenotype (fp, i, SNVsites[j], NO);
					}
				else if (doNGS == NO)
					fprintf (fp, " %c%c", WhichNuc(GetG(MATERNAL, i, SNVsites[j])), WhichNuc(GetG(PATERNAL, i, SNVsites[j])));
				else
					{
					MLPairAt (i, j, &mm, &pp);
					fprintf (fp, " %c%c", WhichNuc(mm), WhichNuc(pp));
					}
				}
           fprintf (fp,"\n");
            }
			
        if (doPrintAncestors == YES && doNGS == NO)
            {
            fprintf (fp,"%s ", outRootCellName);
			for (j=0; j<numSNVs; j++)
				fprintf (fp, " %c%c", WhichNuc(GetG(MATERNAL, OUTGROUP_ROOT, SNVsites[j])),WhichNuc(GetG(PATERNAL, OUTGROUP_ROOT, SNVsites[j])));
            fprintf (fp,"\n%s ", inRootCellName);
			for (j=0; j<numSNVs; j++)
				fprintf (fp, " %c%c", WhichNuc(GetG(MATERNAL, INGROUP_ROOT, SNVsites[j])),WhichNuc(GetG(PATERNAL, INGROUP_ROOT, SNVsites[j])));
            fprintf (fp,"\n");
            }
        }
    else // binary data
        {
        for (i=0; i<NUM_TAXA; i++)
            {
			PrintRowLabel (fp, i, LABEL_TIGHT);
			for (j=0; j<numSNVs; j++)
				{
				if (CN_EVENTS)
					{
					fputc (' ', fp);
					PrintCopyGenotype (fp, i, SNVsites[j], YES);
					}
				else
					fprintf (fp, " %c%c", WhichMut(GetG(MATERNAL, i, SNVsites[j])),WhichMut(GetG(PATERNAL, i, SNVsites[j])));
				}
			fprintf (fp,"\n");
            }
			
        if (doPrintAncestors == YES)
            {
			fprintf (fp,"%s ", outRootCellName);
			for (j=0; j<numSNVs; j++)
				fprintf (fp, " %d%d", GetG(MATERNAL, OUTGROUP_ROOT, SNVsites[j]),GetG(PATERNAL, OUTGROUP_ROOT, SNVsites[j]));
			fprintf (fp,"\n%s ", inRootCellName);
			for (j=0; j<numSNVs; j++)
				fprintf (fp, " %d%d", GetG(MATERNAL, INGROUP_ROOT, SNVsites[j]),GetG(PATERNAL, INGROUP_ROOT, SNVsites[j]));
			fprintf (fp,"\n");
            }
        }
     }

/********************* PrintFullGenotypes ********************/
/* Prints observed/ML genotypes at all sites to a file */

void PrintFullGenotypes (FILE *fp)
    {
    int		i, j, ml1 = MISSING, ml2 = MISSING;
		
	if (doPrintAncestors == YES && doNGS == NO)
        fprintf (fp, "%d %d\n", NUM_TAXA+2, numSites);
    else
        fprintf (fp, "%d %d\n", NUM_TAXA, numSites);
	
    if (alphabet == DNA)
        {
        for (i=0; i<NUM_TAXA; i++)
            {
	LoadFullRow (i);
            PrintRowLabel (fp, i, LABEL_TIGHT);
			for (j=0; j<numSites; j++)
				{
				if (doNGS == NO && CN_EVENTS)
					{
					fputc (' ', fp);
					PrintCopyGenotype (fp, i, j, NO);
					}
				else if (doNGS == NO)
					fprintf (fp, " %c%c", WhichNuc(GetG(MATERNAL, i, j)), WhichNuc(GetG(PATERNAL, i, j)));
				else
					{
					FullMLPair (i, j, &ml1, &ml2);
					fprintf (fp, " %c%c", WhichNuc(ml1), WhichNuc(ml2));
					}
				}
			fprintf (fp,"\n");
            }
		
		/* print ancestral sequences, but not for the NGS case */
        if (doPrintAncestors == YES && doNGS == NO)
            {
            fprintf (fp,"%s ", outRootCellName);
            for (j=0; j<numSites; j++)
				fprintf (fp, " %c%c",  WhichNuc(GetG(MATERNAL, OUTGROUP_ROOT, j)), WhichNuc(GetG(PATERNAL, OUTGROUP_ROOT, j)));
            fprintf (fp,"\n%s ", inRootCellName);

            for (j=0; j<numSites; j++)
                fprintf (fp, " %c%c",  WhichNuc(GetG(MATERNAL, INGROUP_ROOT, j)), WhichNuc(GetG(PATERNAL, INGROUP_ROOT, j)));
            fprintf (fp,"\n");
            }
        }
    else // binary alphabet
        {
        for (i=0; i<NUM_TAXA; i++)
            {
            PrintRowLabel (fp, i, LABEL_TIGHT);
			for (j=0; j<numSites; j++)
                {
                if (CN_EVENTS)
                    {
                    fputc (' ', fp);
                    PrintCopyGenotype (fp, i, j, YES);
                    }
                else
                    fprintf (fp, " %c%c", WhichMut(GetG(0, i, j)),WhichMut(GetG(1, i, j)));
                }
            fprintf (fp,"\n");
			}
			
        if (doPrintAncestors == YES)
            {
            fprintf (fp,"%s ", outRootCellName);
            for (j=0; j<numSites; j++)
                fprintf (fp, " %d%d", GetG(MATERNAL, OUTGROUP_ROOT, j),GetG(PATERNAL, OUTGROUP_ROOT, j));
            fprintf (fp,"\n%s ", inRootCellName);
            for (j=0; j<numSites; j++)
                fprintf (fp, " %d%d", GetG(MATERNAL, INGROUP_ROOT, j),GetG(PATERNAL, INGROUP_ROOT, j));
            fprintf (fp,"\n");
            }
		}
     }

/* what a haplotype row shows at each site */
enum {ROW_IUPAC, ROW_MAT, ROW_PAT, ROW_BIN_CONSENSUS, ROW_BIN_MAT, ROW_BIN_PAT};

/* character shown at site j of a row; ml is the pair of maximum-likelihood
   alleles of the cell when haplotypes are inferred from read counts */
static char HaplotypeChar (int row, int idx, const int *sites, int kind, int useML)
	{
	int		ml1 = MISSING, ml2 = MISSING, j = (sites == NULL) ? idx : sites[idx];

	if (useML && sites == NULL)
		FullMLPair (row, j, &ml1, &ml2);
	else if (useML)
		MLPairAt (row, idx, &ml1, &ml2);
	else
		{
		ml1 = GetG(MATERNAL, row, j);
		ml2 = GetG(PATERNAL, row, j);
		}
	switch (kind)
		{
		case ROW_IUPAC:         return WhichIUPAC (ml1, ml2);
		case ROW_MAT:           return WhichNuc (ml1);
		case ROW_PAT:           return WhichNuc (ml2);
		case ROW_BIN_CONSENSUS: return WhichConsensusBinary (ml1, ml2);
		case ROW_BIN_MAT:       return WhichMut (ml1);
		default:                return WhichMut (ml2);
		}
	}

static void PrintHaplotypeRow (FILE *fp, int row, const int *sites, int n, int kind, int useML)
	{
	int		j;

	for (j=0; j<n; j++)
		fputc (HaplotypeChar (row, j, sites, kind, useML), fp);
	fputc ('\n', fp);
	}

/* label of the row of the rank-th copy (rank >= 2) of a haplotype of a cell: cell0001m2, cell0001p3 */
static void PrintCopyLabel (FILE *fp, int i, char kind, int rank)
	{
	if (doUserTree == NO)
		fprintf (fp, "%s%04d%c%d ", inCellName, i+1, kind, rank);
	else
		fprintf (fp, "%c%d%-*s ", kind, rank, stringPrecision, cellNames[i]);
	}

/* row of the rank-th copy (rank >= 2) of a haplotype: the allele of the copy at each site, or a gap where the cell has fewer copies */
static void PrintCopyRow (FILE *fp, int row, const int *sites, int n, int h, int rank)
	{
	Geno	ex[MAX_COPIES];
	int		j, site, k;

	for (j=0; j<n; j++)
		{
		site = (sites == NULL) ? j : sites[j];
		k = ExtraCopies (row, h, site, ex);
		if (GetG (h, row, site) == DELETION || k < rank - 1)
			fputc ('-', fp);
		else
			fputc (alphabet == DNA ? WhichNuc (ex[rank - 2]) : WhichMut (ex[rank - 2]), fp);
		}
	fputc ('\n', fp);
	}

/* Prints the haplotypes of a list of sites (all of them when sites is NULL): the true ones, or those
   inferred from the read counts when fromReads is YES (then the root rows are not printed) */
static void PrintSiteHaplotypes (FILE *fp, const int *sites, int n, int fromReads, int loadRows)
	{
	int		i, consensus = doPrintIUPAChaplotypes, ancestors = (doPrintAncestors == YES && fromReads == NO);
	int		isDNA = (alphabet == DNA), useML = (fromReads && isDNA);
	int		kindM = isDNA ? ROW_MAT : ROW_BIN_MAT, kindP = isDNA ? ROW_PAT : ROW_BIN_PAT;
	int		kindC = isDNA ? ROW_IUPAC : ROW_BIN_CONSENSUS;
	int		roots[2] = {OUTGROUP_ROOT, INGROUP_ROOT};
	const char	*rootNames[2] = {outRootCellName, inRootCellName};

	int		extraRows = 0, *extraM = NULL, *extraP = NULL, r;

	if (CN_EVENTS && fromReads == NO && consensus == NO && copyGainOn == YES)		/* one row for every copy of every haplotype */
		{
		extraM = (int *) Checked (calloc (NUM_TAXA, sizeof (int)), "Could not allocate the copy rows\n");
		extraP = (int *) Checked (calloc (NUM_TAXA, sizeof (int)), "Could not allocate the copy rows\n");
		for (i=0; i<NUM_TAXA; i++)
			if (i < numCells)
				{
				extraM[i] = MaxExtraCopies (i, MATERNAL);
				extraP[i] = MaxExtraCopies (i, PATERNAL);
				extraRows += extraM[i] + extraP[i];
				}
		}
	fprintf (fp, "%d %d\n", (consensus == YES ? 1 : 2) * (NUM_TAXA + (ancestors ? 2 : 0)) + extraRows, n);

	if (sites != NULL)		/* positions of the sites */
		{
		fprintf (fp, "%d", sites[0]+1);
		for (i=1; i<n; i++)
			fprintf (fp, " %d", sites[i]+1);
		fprintf (fp, "\n");
		}

	for (i=0; i<NUM_TAXA; i++)
		{
		if (loadRows && isDNA)
			LoadFullRow (i);
		if (consensus == YES)
			{
			PrintRowLabel (fp, i, LABEL_SPACED);
			PrintHaplotypeRow (fp, i, sites, n, kindC, useML);
			}
		else
			{
			PrintRowLabel (fp, i, 'm');
			PrintHaplotypeRow (fp, i, sites, n, kindM, useML);
			if (extraM != NULL && i < numCells)
				for (r=0; r<extraM[i]; r++)
					{
					PrintCopyLabel (fp, i, 'm', r + 2);
					PrintCopyRow (fp, i, sites, n, MATERNAL, r + 2);
					}
			PrintRowLabel (fp, i, 'p');
			PrintHaplotypeRow (fp, i, sites, n, kindP, useML);
			if (extraP != NULL && i < numCells)
				for (r=0; r<extraP[i]; r++)
					{
					PrintCopyLabel (fp, i, 'p', r + 2);
					PrintCopyRow (fp, i, sites, n, PATERNAL, r + 2);
					}
			}
		}

	if (ancestors)
		for (i=0; i<2; i++)
			{
			if (consensus == YES)
				{
				fprintf (fp, "%s ", rootNames[i]);
				PrintHaplotypeRow (fp, roots[i], sites, n, kindC, NO);
				}
			else
				{
				fprintf (fp, "%sm  ", rootNames[i]);
				PrintHaplotypeRow (fp, roots[i], sites, n, kindM, NO);
				fprintf (fp, "%sp  ", rootNames[i]);
				PrintHaplotypeRow (fp, roots[i], sites, n, kindP, NO);
				}
			}
	free (extraM);
	free (extraP);
	}

/***************************** PrintSNVHaplotypes *******************************/
/* Prints observed/ML (or true) haplotypes at variable sites (SNVs), to a file */
void PrintSNVHaplotypes (FILE *fp, int PrintTrueVariants)
	{
	PrintSiteHaplotypes (fp, (PrintTrueVariants == YES) ? variantSites : SNVsites, numSNVs, doNGS, NO);
	}

/***************************** PrintFullHaplotypes *******************************/
/* Prints the haplotypes for all sites (variable + invariable) inferred from the read counts */
void PrintFullHaplotypes (FILE *fp)
	{
	PrintSiteHaplotypes (fp, NULL, numSites, doNGS, YES);
	}

/***************************** PrintTrueFullHaplotypes *******************************/
/* Prints true haplotypes for all sites (variable + invariable) to a file */
void PrintTrueFullHaplotypes (FILE *fp)
	{
	PrintSiteHaplotypes (fp, NULL, numSites, NO, NO);
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
 N	unknown state
 . or -	gap

This is what we do:

 A/A => A
 A/C => M
 A/G => R
 A/T => W
 A/_ => a
 A/- => a
 
 C/A => M
 C/C => C
 C/G => S
 C/T => Y
 C/_ => c
 C/- => c
 
 G/A => R
 G/C => S
 G/G => G
 G/T => K
 G/_ => g
 G/- => g
 
 T/A => W
 T/C => Y
 T/G => K
 T/T => T
 T/_ => t
 T/- => t

 _/A => a
 _/C => c
 _/G => g
 _/T => t
 _/_ => -
 _/- => -


 -/A => a
 -/C => c
 -/G => g
 -/T => t
 -/_ => -
 -/- => -

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
		else if (allele2 == ADO)	//A_
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
		if (allele2 == 0)		//_A
			return ('a');
		else if (allele2 == 1)	//_C
			return ('c');
		else if (allele2 == 2)	//_G
			return ('g');
		else if (allele2 == 3)	//_T
			return ('t');
		else if (allele2 == ADO)	//__
			return ('-');
		else if (allele2 == DELETION)	//_-
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
		else if (allele2 == ADO)	//-_
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
		else if (allele2 == DELETION)	//1-
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
/* Returns integer representation for character nucleotides */

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
	    return ('N');
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

void PrintHeader(FILE *fp)
    {
	fprintf (fp, "______________________________________________________________________\n\n");
    fprintf (fp,"Cell coalescent simulation - %s", PROGRAM_NAME);
    fprintf (fp,"  %s", VERSION_NUMBER);
    fprintf (fp,"\n(c) 2016-2022 David Posada - dposada@uvigo.es");
	fprintf (fp, "\n______________________________________________________________________\n\n");
    }



/***************************** PrintDate *******************************/
/* Prints current date and time  */

void PrintDate (FILE *fp)
    {
    time_t now;
    char *date;
		
    now=time(NULL);
    date= ctime(&now);
    fprintf (fp,"%s",date);
    }


/***************************** PrintRunInformation *******************************/
/* Prints a summary of run settings */

void	PrintRunInformation (FILE *fp)
    {
    int		i;
		
    fprintf (fp, "\n\nRun settings\n------------");
	fprintf (fp, "\n[Parameter file = %s]\n", parameterFile);
	if (noisy > 1)
		fprintf (fp, "\n[Assumptions in brackets]\n");
	fprintf (fp, "\nGeneral");
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
    if (doOthsukiInnanCoal == YES)
		{
		fprintf (fp, "\n Cell birth rate                              =   %2.1e", birthRate);
		fprintf (fp, "\n Cell death rate                              =   %2.1e", deathRate);
		fprintf (fp, "\n Exponential growth rate                      =   %2.1e", growthRate);
		}
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
		fprintf (fp, "\n\nGenealogy");
		fprintf (fp, "\n Mean time to the MRCA (Ne generations)       =   %6.4f", runStats.tmrca.mean);
		if (noisy > 1)
			fprintf (fp, "\n Exp time to the MRCA [constant Ne]           =   %6.4f", expTMRCA);
		if (numDataSets > 1)
			{
			fprintf (fp, "\n Variance time MRCA                           =   %6.4f", runStats.tmrca.var);
			if (noisy > 1)
				fprintf (fp, "\n Exp variance time MRCA [constant Ne]         =   %6.4f", expVarTMRCA);
			}
	 
		if (rateVarAmongLineages == YES)
			fprintf (fp, "\n Rate variation among branches (alpha)        =   %-3.2f", alphaBranches);
		else if (doTreeRateSwitches == YES)
			fprintf (fp, "\n Number of tree rate switches                 =   %-3d",  numTreeRateSwitches);
		else
			fprintf (fp, "\n No rate variation among branches");
		}
    fprintf (fp, "\n\nUser-defined branches");
    fprintf (fp, "\n Transforming branch length ratio             =   %2.1e", transformingBranchLengthRatio);
    fprintf (fp, "\n Healthy tip branch length ratio              =   %2.1e", healthyTipBranchLengthRatio);

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
            fprintf (fp, "\n Population mutation parameter (4NuL)         =   %3.2f", theta);
            if (noisy > 1)
				fprintf (fp, "\n Exp num mutation events [Ne=cte,ISMh]        =   %3.2f", expNumMU);
            }
        fprintf (fp, "\n Mean number of mutation events               =   %3.2f", runStats.mu.mean);
        if (numDataSets > 1)
            {
            if (doSimulateFixedNumMutations == NO)
				{
                if (noisy > 1)
					fprintf (fp, "\n Exp variance num mut events [Ne=cte,ISM]     =   %3.2f", expVarNumMU);
				fprintf (fp, "\n Variance number of mutation events           =   %3.2f", runStats.mu.var);
				}
             }
        fprintf (fp, "\n Mean number of SNVs                          =   %3.2f", runStats.snv.mean);
        if (numDataSets > 1)
            {
			fprintf (fp, "\n Variance number of SNVs                      =   %3.2f", runStats.snv.var);
			fprintf (fp, "\n Proportion of reps without SNV sites         =   %3.2f", zeroSNVs/(double)numDataSets);
			}

		if (deletionRate > 0)
			{
			fprintf (fp, "\n\nDeletion rates");
			fprintf (fp, "\n Deletion rate                                =   %2.1e", deletionRate);
			fprintf (fp, "\n Mean number of deletion events               =   %3.2f", runStats.del.mean);
			if (numDataSets > 1)
				fprintf (fp, "\n Variance number of deletion events           =   %3.2f", runStats.del.var);
			}

		if (CNLOHrate > 0)
			{
			fprintf (fp, "\n\nCopy-neutral LOH rates");
			fprintf (fp, "\n CN-LOH rate                                  =   %2.1e", CNLOHrate);
			fprintf (fp, "\n Mean number of CN-LOH events                 =   %3.2f", runStats.cnloh.mean);
			if (numDataSets > 1)
				fprintf (fp, "\n Variance number of CN-LOH events             =   %3.2f", runStats.cnloh.var);
			}
 
		fprintf (fp, "\n\nSingle-cell WGA");
		if (fixedADOrate == 0 && doADOcell == NO && doADOsite == NO)
			fprintf (fp, "\n Alellic dropout                              =   %2.1f", fixedADOrate);
		if (fixedADOrate > 0)
			fprintf (fp, "\n Alellic dropout fixed                        =   %2.1e", fixedADOrate);
		if (doADOcell == YES)
			fprintf (fp, "\n Alellic dropout per cell (mean, var)         =   %2.1e, %2.1e", meanADOcell, varADOcell);
		if (doADOsite == YES)
			fprintf (fp, "\n Alellic dropout per site (mean, var)         =   %2.1e, %2.1e", meanADOsite, varADOsite);

		if (doNGS == YES)
			{
			if (doAllelicImbalance == YES)
				fprintf (fp, "\n Maternal allelic imbalance (mean, var)       =   %2.1e, %2.1e", meanAllelicImbalance, varAllelicImbalance);
			else
				fprintf (fp, "\n Maternal allelic imbalance                   =   %2.1f", meanAllelicImbalance);

			if (meanAmplificationError > 0)
				{
				fprintf (fp, "\n Amplification error ");
				if (simulateOnlyTwoTemplates == YES)
					fprintf (fp, "\n  (2-template model)");
				else
					fprintf (fp, "\n  (4-template model)");
				fprintf (fp, " (mean, var)              =   %2.1e, %2.1e", meanAmplificationError, varAmplificationError);
				if (ampCellVar > 0)
					fprintf (fp, "\n  variance among cells of the mean           =   %2.1e", ampCellVar);
				}
			else
				fprintf (fp, "\n Amplification error                          =   %2.1f", meanAmplificationError);

			if (meanDoubletRate > 0)
				fprintf (fp, "\n Doublet rate per cell (mean, var)            =   %2.1e, %2.1e", meanDoubletRate, varDoubletRate);
			else
				fprintf (fp, "\n Doublet rate per cell                        =   %2.1f", meanDoubletRate);
			}
	
		if (meanGenotypingError == 0 && doNGS == NO)
			fprintf (fp, "\n\nNGS genotype errors are not being simulated");

		if (meanGenotypingError > 0)
			{
			fprintf (fp, "\n\nGenotyping");
			fprintf (fp, "\n Genotype error (mean, var)                   =   %2.1e, %2.1e", meanGenotypingError, varGenotypingError);
			//fprintf (fp, "\n Exp number of FP SNVs due to genotype errors =   %3.2f", (1 - pow(1-genotypingError, NUM_TAXA)) * numSites);
			}

		fprintf (fp, "\n\nNGS");
		if (doNGS == YES)
			{
			if (rateVarCoverage == NO)
				fprintf (fp, "\n Sequencing coverage (Poisson)                =   %-3.2fX", coverage);
			else
				{
				fprintf (fp, "\n Sequencing coverage (Negative binomial)");
				fprintf (fp, "\n  mean                                        =   %-3.2fX", coverage);
				fprintf (fp, "\n  dispersion (alpha)                          =   %-3.2f", alphaCoverage);
				}
			if (alphaSiteCoverage > 0.0)
				fprintf (fp, "\n Coverage dispersion between sites (alpha)    =   %-3.2f", alphaSiteCoverage);
			if (alphaCellCoverage > 0.0)
				fprintf (fp, "\n Coverage dispersion between cells (alpha)    =   %-3.2f", alphaCellCoverage);
			fprintf (fp, "\n Haploid coverage reduction                   =   %-3.2f", haploidCoverageReduction);
			fprintf (fp, "\n Sequencing error                             =   %2.1e", sequencingError);
			if (thereIsEij == YES)
				{
				fprintf (fp, "\n NGS error rate matrix                        =   %3.2f %3.2f %3.2f %3.2f", Eij[0][0], Eij[0][1], Eij[0][2], Eij[0][3]);
				fprintf (fp, "\n                                                  %3.2f %3.2f %3.2f %3.2f", Eij[1][0], Eij[1][1], Eij[1][2], Eij[1][3]);
				fprintf (fp, "\n                                                  %3.2f %3.2f %3.2f %3.2f", Eij[2][0], Eij[2][1], Eij[2][2], Eij[2][3]);
				fprintf (fp, "\n                                                  %3.2f %3.2f %3.2f %3.2f", Eij[3][0], Eij[3][1], Eij[3][2], Eij[3][3]);
				}
			}
		else
			fprintf (fp, "\n No read counts were produced");
		}

	   fprintf (fp, "\n\n\nOutput files are in folder \"%s\":", resultsDir);
		   
		   if (doPrintTree == YES)
			   {
			   fprintf (fp, "\n Tree printed to file \"%s\"", treeFile);
			   if (doPrintSeparateReplicates == YES)
				   fprintf (fp, " in folder \"%s\"", treeDir);
	   
			   }
	
		   if (doPrintTimes == YES)
			   {
			   fprintf (fp, "\n Times printed to file \"%s\"", timesFile);
			   if (doPrintSeparateReplicates == YES)
				   fprintf (fp, " in folder \"%s\"", timesDir);
			   }
		   
		   if (doSimulateData == YES)
			   {
			   if (doPrintSNVgenotypes == YES)
				   {
				   fprintf (fp, "\n SNV genotypes printed to file \"%s\"", SNVgenotypesFile);
				   if (doPrintSeparateReplicates == YES)
					   fprintf (fp, " in folder \"%s\"", SNVgenotypesDir);
				  }
	
			   if (doPrintFullGenotypes == YES)
				   {
				   fprintf (fp, "\n Full genotypes printed to file \"%s\"", fullGenotypesFile);
				   if (doPrintSeparateReplicates == YES)
					   fprintf (fp, " in folder \"%s\"", fullGenotypesDir);
				   }
	
			   if (doPrintSNVhaplotypes == YES)
				   {
				   if (doPrintIUPAChaplotypes == YES)
					   fprintf (fp, "\n SNV haplotypes (IUPAC codes) printed to file \"%s\"", SNVhaplotypesFile);
				   else
					   fprintf (fp, "\n SNV haplotypes printed to file \"%s\"", SNVhaplotypesFile);
				   if (doPrintSeparateReplicates == YES)
					   fprintf (fp, " in folder \"%s\"", SNVhaplotypesDir);
				   }

			   if (doPrintFullHaplotypes == YES)
				   {
				   if (doPrintIUPAChaplotypes == YES)
					   fprintf (fp, "\n Full haplotypes (IUPAC codes) printed to file \"%s\"", fullHaplotypesFile);
				   else
					   fprintf (fp, "\n Full haplotypes printed to file \"%s\"", fullHaplotypesFile);
				   if (doPrintSeparateReplicates == YES)
					   fprintf (fp, " in folder \"%s\"", fullHaplotypesDir);
				   }

			   if (doPrintTrueHaplotypes == YES)
				   {
				   if (doPrintIUPAChaplotypes == YES)
					   fprintf (fp, "\n True haplotypes (IUPAC codes) printed to file \"%s\"", trueHaplotypesFile);
				   else
					   fprintf (fp, "\n True haplotypes printed to file \"%s\"", trueHaplotypesFile);
				   if (doPrintSeparateReplicates == YES)
					   fprintf (fp, " in folder \"%s\"", trueHaplotypesDir);
				   }

			   if (doNGS == YES && doPrintVCF == YES)
				   {
				   fprintf (fp, "\n Genotype likelihoods (%s model) printed to file \"%s%s\"", GLmodel, VCFfile, (doZipVCF == YES) ? ".gz" : "");
				   if (doPrintSeparateReplicates == YES)
					   fprintf (fp, " in folder \"%s\"", VCFdir);
				   }

			   if (doPrintCATG == YES)
				   {
				   fprintf (fp, "\n Read counts printed in CATG format to file \"%s\"", CATGfile);
				   if (doPrintSeparateReplicates == YES)
					   fprintf (fp, " in folder \"%s\"", CATGdir);
				   }

			   if (doPrintPILEUP == YES)
				   {
				   fprintf (fp, "\n Read counts printed in PILEUP format to file \"%s\"", PILEUPfile);
				   if (doPrintSeparateReplicates == YES)
					   fprintf (fp, " in folder \"%s\"", PILEUPdir);
				   }

			   if (doPrintBulk == YES)
				   {
				   fprintf (fp, "\n Pooled (bulk) read counts and VAF printed to file \"%s\"", bulkFile);
				   if (doPrintSeparateReplicates == YES)
					   fprintf (fp, " in folder \"%s\"", bulkDir);
				   }
			  }
		   else
			   fprintf (fp, "\n\n\nNo output data files");

	}

/***************** PrintCommandLine **********************/
/*	Prints the options of the run: the command line, or the equivalent command line of the parameter file */
void PrintCommandLine (FILE *fp)
	{
	if (readingParameterFile == NO)
		{
		fprintf (fp, "Command line arguments (some might not affect) = %s", CommandLine);
		return;
		}
	fprintf (fp, "Parameter file [\"%s\"] arguments (some might not affect) = %s", parameterFile, paramFileArgs ? paramFileArgs : "");
	if (paramFileArgs == NULL || strstr (paramFileArgs, " -#") == NULL)
		fprintf (fp, " -#%ld", originalSeed);
	}


/***************************** PrintUsage *******************************/
/* Prints a short description of program usage  */
/* every option (--help, -?): the basic ones, the groups of related parameters, the run options and the older single options */
void PrintUsageAll (FILE *fp)
{
	PrintHeader(fp);
    fprintf (fp, "\n%s generates a coalescent tree and simulates a sample of diploid genomes from somatic cells (no recombination), together with a (healthy) cell as outgroup.", PROGRAM_NAME);
	fprintf (fp, "\n\nUsage: %s [-n# -s# -l# -e# -u# -C# -#SEED] [--GROUP KEY=VALUE,KEY,...] [--config FILE] [--check]", PROGRAM_NAME);
	fprintf (fp,"\n\nThree modes:");
	fprintf (fp,"\n  genotypes    no coverage (-C): genotypes with allelic dropout and genotyping error; no reads.");
	fprintf (fp,"\n  read counts  coverage (-C): read counts per base, cell and site (CellCoal does not simulate individual reads); VCF with genotype likelihoods.");
	fprintf (fp,"\n  reads        coverage with --reads bam,fastq: paired-end reads (BAM, FASTQ) of the amplified cells, with truth.vcf.");
	fprintf (fp,"\nAn option that the mode does not use is an error (reads) or a note (genotypes). --check shows the mode.");

	fprintf (fp,"\n\n*** BASIC OPTIONS ***");
	fprintf (fp,"\n-n: number of replicates (e.g. -n1000)");
	fprintf (fp,"\n-s: sample size (# cells) (e.g. -s8)");
	fprintf (fp,"\n-l: number of sites (e.g. -l500)");
	fprintf (fp,"\n-e: effective population size (e.g. -e1000)");
	fprintf (fp,"\n-u: somatic mutation rate (e.g. -u1e-6)  [also --mutation rate=]");
	fprintf (fp,"\n-C: mean sequencing coverage; without it no reads are simulated (e.g. -C60)  [also --depth mean=]");
	fprintf (fp,"\n-#: seed (e.g. -#37864287)");

	fprintf (fp,"\n\n*** GROUPS OF PARAMETERS ***");
	fprintf (fp,"\nThe parameters that belong together are given in one option, as a list of keys separated by commas; a key has a value (key=value;");
	fprintf (fp,"\nseveral numbers are separated by ':') or is a switch. Example:");
	fprintf (fp,"\n  %s -n10 -s20 -l1000 -u1e-6 -C10 --ado mean=0.1,block=3000 --amplification error=0.05:0.01:0 --reads bam,fastq --output snv-genotypes,vcf-gz", PROGRAM_NAME);
	PrintGroupHelp (fp);

	fprintf (fp,"\n\n*** RUN ***");
	fprintf (fp,"\n--config FILE: read the options from a JSON or YAML file (the keys are the option names, e.g. cells: 20, or a group as a string, e.g. ado: mean=0.1,block=3000)");
	fprintf (fp,"\n--check: check the options and show what would run, without simulating");
	fprintf (fp,"\n--fast: faster read simulation (same distributions, different random numbers; see the manual)");
	fprintf (fp,"\n--sweep \"KEY=V1,V2 ...\": run the grid of values (e.g. --sweep \"C=2,5,10 E=0.001,0.01\") and write sweep.tsv with mean and standard error of the metrics; --sweep-keep keeps the runs; --sweep-jobs N runs N at a time; --sweep-resume skips the finished runs");
	fprintf (fp,"\n-F: user-defined parameter file (older: --config) (e.g. -Fmyparameters)");
	fprintf (fp,"\n--version: print the version and exit");
	fprintf (fp,"\n-?, --help: print this help");

	fprintf (fp,"\n\n*** SINGLE OPTIONS (older commands) ***");
	fprintf (fp,"\nEvery key of a group is also an option of its own; this is the table of equivalences (group key = option):");
	PrintGroupEquivalents (fp);

	PrintConfigKeys (fp);
	fprintf (fp,"\n\nDefaults: ");
	PrintDefaults (fp);
	fprintf (fp, "\n-------------------------------------------------------------------------------------------------------");
	fprintf (fp,"\n");
}

/* the whole list of options (--help, -?) */
void PrintHelp (FILE *fp)
{
	PrintUsageAll (fp);
}

/* a wrong option: the list of options and exit */
void PrintUsage (FILE *fp)
{
	PrintUsageAll (fp);
	fprintf (fp, "\n...now exiting to system...\n\n");
	exit(-1);
}


/***************** PrintDefaults **********************/
/*	Prints the default settings */

void PrintDefaults (FILE *fp)
	{
	int		i;

	fprintf (fp,"\n-F: parameter file =  %s", "parameters");

	/* Coalescent */
	fprintf (fp,"\n-n: number of replicates =  %d", numDataSets);
	fprintf (fp,"\n-s: sample size (# cells) =  %d", numCells);
	fprintf (fp,"\n-l: number of sites =  %d", numSites);
	fprintf (fp,"\n-e: effective population size =  %d", N);
	fprintf (fp,"\n-g: exponential growth rate =  %2.1e", growthRate);
	fprintf (fp,"\n-K: cell birth rate =  %2.1e", birthRate);
	fprintf (fp,"\n-L: cell death rate =  %2.1e", deathRate);
	fprintf (fp,"\n-h: number of demographic periods = %d", numPeriods);
	for (i=1; i<=numPeriods; i++)
		fprintf (fp, "\n  period %d =  %d %d %lld", i, Nbegin[i], Nend[i], (long long) cumDuration[i] - (long long) cumDuration[i-1]);


	/* Post-coalescent */
	fprintf (fp,"\n-k: transforming branch length ratio =  %2.1e", transformingBranchLengthRatio);
	fprintf (fp,"\n-q: healthy tip branch length ratio =  %2.1e", healthyTipBranchLengthRatio);
	fprintf (fp,"\n-i: shape of the gamma distribution for rate variation among lineages =  %6.4f", alphaBranches);
	
	/* Mutation model */
	fprintf (fp,"\n-b: alphabet [0:binary 1:DNA] =  %d", alphabet);
	fprintf (fp,"\n-c: germline SNP rate =  %2.1e", SNPrate);
	fprintf (fp,"\n-u: mutation rate =  %2.1e", mutationRate);
	fprintf (fp,"\n-d: deletion rate =  %2.1e", deletionRate);
	fprintf (fp,"\n-H: CNLOH rate =  %2.1e", CNLOHrate);
	fprintf (fp,"\n-j: fixed number of mutations =  %d", numFixedMutations);
	fprintf (fp,"\n-S: number of trinucleotide genetic signature =  %d", numUserSignatures);
	fprintf (fp,"\n-m: alternative mutation model =  %d", altModel);
	fprintf (fp,"\n-p: proportion of alternative model sites =  %6.4f", propAltModelSites);
	fprintf (fp,"\n-w: alternative/default model relative mutation rate =  %6.4f", nonISMRelMutRate);
	fprintf (fp,"\n-f: nucleotide base frequencies = f%3.2f %3.2f %3.2f %3.2f", freq[0], freq[1], freq[2], freq[3]);
	fprintf (fp,"\n-t: transition/transversion ratio =  %6.4f", titv);
	fprintf (fp,"\n-a: shape of the gamma distribution for rate variation among sites =  %6.4f", alphaSites);
	fprintf (fp,"\n-r: mutation matrix ACGT x ACGT = %3.2f %3.2f %3.2f %3.2f  %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f", Mij[0][0], Mij[0][1], Mij[0][2], Mij[0][3], Mij[1][0], Mij[1][1], Mij[1][2], Mij[1][3],  Mij[2][0], Mij[2][1], Mij[2][2], Mij[2][3],  Mij[3][0], Mij[3][1], Mij[3][2], Mij[3][3]);
	
	/* scWGA */
	fprintf (fp,"\n-D: fixed allelic dropout (ADO) =  %2.1e", fixedADOrate);
	fprintf (fp,"\n-P: ADO per cell mean, var  =  %2.1e %2.1e", meanADOcell, varADOcell);
	fprintf (fp,"\n-Q: ADO per site mean, var =  %2.1e %2.1e", meanADOsite, varADOsite);
	fprintf (fp,"\n-I: allelic imbalance =  %2.1e %2.1e", meanAllelicImbalance, varAllelicImbalance);
	fprintf (fp,"\n-A: amplification error: mean, var, 2-template model =  %2.1e, %2.1e, %d", meanAmplificationError, varAmplificationError, simulateOnlyTwoTemplates);
	fprintf (fp,"\n-B: doublet rate per cell; mean, var =  %2.1e, %2.1e", meanDoubletRate, varDoubletRate);

	/* GENOTYPING */
	fprintf (fp,"\n-G: genotyping error: mean, var =  %2.1e %2.1e", meanGenotypingError, varGenotypingError);

	/* NGS */
	fprintf (fp,"\n-C: sequencing coverage =  %3.2f", coverage);
	fprintf (fp,"\n-V: sequencing coverage overdispersion =  %6.4f", alphaCoverage);
	fprintf (fp,"\n-O: coverage overdispersion between sites =  %6.4f", alphaSiteCoverage);
	fprintf (fp,"\n-J: coverage overdispersion between cells =  %6.4f", alphaCellCoverage);
	fprintf (fp,"\n-R: haploid coverage reduction =  %2.1e", haploidCoverageReduction);
	fprintf (fp,"\n-E: sequencing error =  %2.1e", sequencingError);
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
	fprintf (fp,"\n-9: print read counts in PILEUP format =  %d", doPrintPILEUP);
	fprintf (fp,"\n-0: print pooled (bulk) read counts and VAF per site =  %d", doPrintBulk);
	fprintf (fp,"\n-Z: print the VCF file =  %d", doPrintVCF);
	fprintf (fp,"\n-z: compress the VCF file with gzip =  %d", doZipVCF);
	fprintf (fp,"\n-M: minimum pooled VAF of the sites in the VCF =  %g", minVAFvcf);
	fprintf (fp,"\n-Y: print true haplotypes to a file =  %d", doPrintTrueHaplotypes);
	fprintf (fp,"\n-v: print replicates in individual folders =  %d", doPrintSeparateReplicates);
	fprintf (fp,"\n-x: print consensus/IUPAC haplotypes =  %d", doPrintIUPAChaplotypes);
	if (strlen(userTreeFile) == 0)
		strcpy(userTreeFile, "none");
	fprintf (fp,"\n-T: user tree file name =  %s", userTreeFile);
	if (strlen(userGenomeFile) == 0)
		strcpy(userGenomeFile, "none");
	fprintf (fp,"\n-U: user genome file name =  %s", userGenomeFile);
	if (strlen(resultsDir) == 0)
		strcpy(resultsDir, "Results");
	fprintf (fp,"\n-o: results folder name =  %s", resultsDir);

	/* Other */
	fprintf (fp,"\n-W: use tumor nomenclature for cells =  %d", doTumorNames);
	fprintf (fp,"\n-y: noisy = %d", noisy);
	fprintf (fp,"\n-#: seed = %ld", originalSeed);
	}




/* Writes the allelic dropout rate drawn for every cell (-P) to the file ado_cells (or ado_cells_dir/ado_cells.NNNN with -v).
   This is the rate of the cell itself; the ADO of a genotype also depends on the rate of its site (-Q). */
void WriteADOCells (const double *rate, int replicate)
	{
	char	file[4*MAX_NAME];
	FILE	*fp;
	int		i;

	if (doPrintSeparateReplicates == YES)
		{
		snprintf (file, sizeof(file), "%s/ado_cells_dir", resultsDir);
		MAKE_DIR (file);
		snprintf (file, sizeof(file), "%s/ado_cells_dir/ado_cells.%04d", resultsDir, replicate+1);
		fp = fopen (file, "w");
		}
	else
		{
		snprintf (file, sizeof(file), "%s/ado_cells", resultsDir);
		fp = fopen (file, replicate == 0 ? "w" : "a");
		}
	if (!fp)
		{
		fprintf (stderr, "Can't open \"%s\"\n", file);
		exit (-1);
		}
	if (doPrintSeparateReplicates == NO)
		fprintf (fp, "[#%d]\n", replicate+1);
	fprintf (fp, "cell\tado_rate\n");
	for (i=0; i<NUM_TAXA; i++)
		{
		if (i == numCells)
			fprintf (fp, "%s", outCellName);
		else if (doUserTree == NO)
			fprintf (fp, "%s%04d", inCellName, i+1);
		else
			fprintf (fp, "%s", cellNames[i]);
		fprintf (fp, "\t%.6f\n", rate[i]);
		}
	fclose (fp);
	}


/***************************** PrintCheck *******************************/
/* --check: the options have been read and validated; says what the run would do without running it */

void PrintCheck (FILE *fp)
	{
	struct stat	st;
	int			any = 0;

	fprintf (fp, "CellCoal %s: the options are valid. Nothing was simulated.\n\n", VERSION_NUMBER);
	fprintf (fp, "  mode:                %s\n", RunMode () == 2 ? "reads (paired-end reads of the amplified cells: --bam/--fastq)" : RunMode () == 1 ? "read counts (-C: read counts per base, cell and site; no individual reads)" : "genotypes (no -C: genotypes with ADO and genotyping error, no reads)");
	fprintf (fp, "  replicates:          %d\n", numDataSets);
	fprintf (fp, "  cells:               %d%s\n", numCells, numOutgroup == 0 ? " (tumor only: no outgroup cell)" : " (and one healthy outgroup cell)");
	fprintf (fp, "  sites:               %d\n", numSites);
	if (doSimulateData == NO)
		fprintf (fp, "  data:                only the genealogy is simulated (no mutations)\n");
	else if (doSimulateFixedNumMutations == YES)
		fprintf (fp, "  mutations:           %d fixed in each replicate\n", numFixedMutations);
	else
		fprintf (fp, "  mutations:           about %.0f per replicate expected under the model (the number varies among replicates)\n", expNumMU);
	fprintf (fp, "  expected tree depth: %.3f (in units of N generations, mean TMRCA of the tumour cells)\n", expTMRCA);
	if (doSimulateData == YES)
		{
		if (doNGS == YES)
			fprintf (fp, "  read counts:         yes (coverage %.1f; CellCoal simulates read counts, not individual reads)\n", coverage);
		else
			fprintf (fp, "  read counts:         no\n");
		}
	fprintf (fp, "  seed:                %s\n", userSeed > 0 ? "fixed by the user (-#)" : "from the clock (it is written in the log)");
	fprintf (fp, "  threads:             ");
#ifdef _OPENMP
	fprintf (fp, "up to %d (OMP_NUM_THREADS)\n", omp_get_max_threads ());
#else
	fprintf (fp, "1 (built without OpenMP)\n");
#endif
	fprintf (fp, "  outputs:            ");
	if (doPrintTree == YES) { fprintf (fp, " trees"); any = 1; }
	if (doPrintSNVgenotypes == YES) { fprintf (fp, " snv_gen"); any = 1; }
	if (doPrintSNVhaplotypes == YES) { fprintf (fp, " snv_hap"); any = 1; }
	if (doPrintTrueHaplotypes == YES || doPrintFullHaplotypes == YES || doPrintFullGenotypes == YES) { fprintf (fp, " full genomes"); any = 1; }
	if (doPrintCATG == YES) { fprintf (fp, " catg"); any = 1; }
	if (doPrintPILEUP == YES) { fprintf (fp, " pileup"); any = 1; }
	if (doPrintBulk == YES) { fprintf (fp, " bulk"); any = 1; }
	if (doPrintVCF == YES && doNGS == YES && readsOnly == 0) { fprintf (fp, " vcf%s", doZipVCF == YES ? " (compressed)" : ""); any = 1; }
	if (doSummary == YES) { fprintf (fp, " summary.tsv"); any = 1; }
	fprintf (fp, "%s plus log and run.json\n", any ? "," : "");
	if (doPrintVCF == YES && doNGS == YES && readsOnly == 0)
		fprintf (fp, "  note:                the VCF is the largest file (roughly 200-350 bytes per cell and variable site; it grows with the coverage)\n");
	fprintf (fp, "  results folder:      %s", strlen (resultsDir) == 0 ? "Results" : resultsDir);
	if (stat (strlen (resultsDir) == 0 ? "Results" : resultsDir, &st) == 0)
		fprintf (fp, "  (it exists: a run would delete it and write it again)");
	fprintf (fp, "\n");
	if (doUserTree == YES)
		{
		fprintf (fp, "  user tree file:      %s%s\n", userTreeFile, stat (userTreeFile, &st) == 0 ? "" : "  !!! NOT FOUND");
		if (ignoredNumCells > 0)
			fprintf (fp, "  note:                -s%d is ignored: the user tree has %d leaves (cells)\n", ignoredNumCells, numCells);
		}
	if (doUserGenome == YES)
		fprintf (fp, "  user genome file:    %s%s\n", userGenomeFile, stat (userGenomeFile, &st) == 0 ? "" : "  !!! NOT FOUND");
	}
