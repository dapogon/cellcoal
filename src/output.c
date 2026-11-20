/*
	CellCoal - result files and the printing of genotypes, haplotypes, trees, times and information
	See cellcoal.h for the shared declarations. This file was split from cellcoal.c.
*/
#include "cellcoal.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "definitions.h"

/********************* PrepareGlobalFiles **********************/
/* Open global files to output results */
void PrepareGlobalFiles(int argc, char **argv)
    {
    char File[4*MAX_NAME];
		
    /* contains the simulated tree in Newick format */
    if (doPrintTree == YES)
        {
		snprintf (File, sizeof(File),"%s/%s", resultsDir, treeFile);
        if ((fpTrees = fopen(File, "w")) == NULL)
            {
            fprintf (stderr, "Can't open \"%s\"\n", File);
            exit(-1);
            }
        }
		
    /* contains a list of times and branch lenghts for all nodes in the simulated trees */
    if (doPrintTimes == YES)
        {
		snprintf (File, sizeof(File),"%s/%s", resultsDir, timesFile);
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
			snprintf (File, sizeof(File),"%s/%s", resultsDir, SNVgenotypesFile);
             if ((fpSNVgenotypes = fopen(File, "w")) == NULL)
                {
				fprintf (stderr, "Can't open \"%s\"\n", File);
				exit(-1);
                }
            }
			
        /* contains haplotypes for variable sites for every cell */
        if (doPrintSNVhaplotypes == YES)
            {
			snprintf (File, sizeof(File),"%s/%s", resultsDir, SNVhaplotypesFile);
            if ((fpSNVhaplotypes = fopen(File, "w")) == NULL)
                {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
                }
            }

        /* contains error-free haplotypes for variable sites for every cell */
        if (doPrintTrueHaplotypes == YES)
            {
			snprintf (File, sizeof(File),"%s/%s", resultsDir, trueHaplotypesFile);
            if ((fpTrueHaplotypes = fopen(File, "w")) == NULL)
                {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
                }
            }

        /* contains all genotypes (variable or invariable) for every cell */
        if (doPrintFullGenotypes == YES)
            {
			snprintf (File, sizeof(File),"%s/%s", resultsDir, fullGenotypesFile);
            if ((fpFullGenotypes = fopen(File, "w")) == NULL)
                {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
                }
            }

        /* contains haplotypes for all sites for every cell */
        if (doPrintFullHaplotypes == YES)
            {
			snprintf (File, sizeof(File),"%s/%s", resultsDir, fullHaplotypesFile);
            if ((fpFullHaplotypes = fopen(File, "w")) == NULL)
                {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
                }
            }

        /* contains reads counts and genotype likelihoods for every SNV and cell */
        if (doNGS == YES && doPrintVCF == YES)
            {
			snprintf (File, sizeof(File),"%s/%s", resultsDir, VCFfile);
            if ((fpVCF = OpenVCF (File)) == NULL)
                {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
                }
            }
  
         /* contains reads counts for every SNV and cell */
        if (doPrintCATG == YES)
            {
			snprintf (File, sizeof(File),"%s/%s", resultsDir, CATGfile);
            if ((fpCATG = fopen(File, "w")) == NULL)
                {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
                }
            }

		if (doPrintPILEUP == YES)
			{
			snprintf (File, sizeof(File),"%s/%s", resultsDir, PILEUPfile);
			if ((fpPILEUP = fopen(File, "w")) == NULL)
				{
				fprintf (stderr, "Can't open \"%s\"\n", File);
				exit(-1);
				}
			}

		if (doPrintBulk == YES)
			{
			snprintf (File, sizeof(File),"%s/%s", resultsDir, bulkFile);
			if ((fpBulk = fopen(File, "w")) == NULL)
				{
				fprintf (stderr, "Can't open \"%s\"\n", File);
				exit(-1);
				}
			}

       }

   if (doPrintSNVgenotypes == YES)
        {
       /*fprintf (fpSNVgenotypes, "%s - ",PROGRAM_NAME);
        PrintDate (fpSNVgenotypes);
        fprintf (fpSNVgenotypes, "SNV genotypes\n");
        PrintCommandLine (fpSNVgenotypes, argc, argv);*/
        fprintf (fpSNVgenotypes,"%d\n", numDataSets);
        }

    if (doPrintSNVhaplotypes == YES)
        {
        /*fprintf (fpSNVhaplotypes, "%s - ",PROGRAM_NAME);
        PrintDate (fpSNVhaplotypes);
        fprintf (fpSNVhaplotypes, "SNV haplotypes\n");
        PrintCommandLine (fpSNVhaplotypes, argc, argv);*/
        fprintf (fpSNVhaplotypes,"%d\n", numDataSets);
        }

   if (doPrintFullGenotypes == YES)
        {
        /*fprintf (fpFullGenotypes, "%s - ",PROGRAM_NAME);
        PrintDate (fpFullGenotypes);
        fprintf (fpFullGenotypes, "Full genotypes\n");
        PrintCommandLine (fpFullGenotypes, argc, argv);*/
        fprintf (fpFullGenotypes,"%d\n", numDataSets);
        }

    if (doPrintFullHaplotypes == YES)
        {
        /*fprintf (fpFullHaplotypes, "%s - ",PROGRAM_NAME);
        PrintDate (fpFullHaplotypes);
        fprintf (fpFullHaplotypes, "Full haplotypes\n");
        PrintCommandLine (fpFullHaplotypes, argc, argv);*/
        fprintf (fpFullHaplotypes,"%d\n", numDataSets);
        }

	if (doPrintTrueHaplotypes == YES)
		{
		/*fprintf (fpTrueHaplotypes, "%s - ",PROGRAM_NAME);
		 PrintDate (fpTrueHaplotypes);
		 fprintf (fpTrueHaplotypes, "True haplotypes\n");
		 PrintCommandLine (fpTrueHaplotypes, argc, argv);*/
		fprintf (fpTrueHaplotypes,"%d\n", numDataSets);
		}
	
	if (doNGS == YES && doPrintVCF == YES)
        {
        /*fprintf (fpVCF, "%s - ",PROGRAM_NAME);
        PrintDate (fpVCF);
        fprintf (fpVCF, "Read counts and genotype likelihoods\n");
        PrintCommandLine (fpVCF, argc, argv);*/
        fprintf (fpVCF,"%d\n", numDataSets);
        }

	if (doPrintCATG == YES)
        {
        /*fprintf (fpCATG, "%s - ",PROGRAM_NAME);
        PrintDate (fpCATG);
        fprintf (fpCATG, "Read counts\n");
        PrintCommandLine (fpCATG, argc, argv);*/
        fprintf (fpCATG,"%d\n", numDataSets);
        }

	if (doPrintPILEUP == YES)
		{
		/*fprintf (fpPILEUP, "%s - ",PROGRAM_NAME);
		PrintDate (fpPILEUP);
		fprintf (fpPILEUP, "Read counts\n");
		PrintCommandLine (fpPILEUP, argc, argv);*/
		fprintf (fpPILEUP,"%d\n", numDataSets);
		}

	if (doPrintBulk == YES)
		{
		fprintf (fpBulk,"%d\n", numDataSets);
		fprintf (fpBulk,"site\tref\treadsA\treadsC\treadsG\treadsT\tdepth\talt_reads\tbulk_vaf\ttrue_alt_alleles\ttrue_vaf\n");
		}
	}



/********************* PrepareSeparateFiles **********************/
/* Open individual files to output results */

void PrepareSeparateFiles(int replicate)
    {
	/* contains the simulated tree in Newick format */
    if (doPrintTree == YES)
        {
		snprintf (File, sizeof(File),"%s/%s", resultsDir, treeDir);
		MAKE_DIR (File);
		snprintf (File, sizeof(File),"%s/%s/%s.%04d", resultsDir, treeDir, treeFile, replicate+1);
		if ((fpTrees = fopen(File, "w")) == NULL)
            {
            fprintf (stderr, "Can't open \"%s\"\n", File);
            exit(-1);
            }
        }
		
	if (doPrintTimes == YES)
        {
		snprintf (File, sizeof(File),"%s/%s", resultsDir, timesDir);
		MAKE_DIR (File);
		snprintf (File, sizeof(File),"%s/%s/%s.%04d", resultsDir, timesDir, timesFile, replicate+1);
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
			snprintf (File, sizeof(File),"%s/%s", resultsDir, SNVgenotypesDir);
			MAKE_DIR (File);
			snprintf (File, sizeof(File),"%s/%s/%s.%04d", resultsDir, SNVgenotypesDir, SNVgenotypesFile, replicate+1);
            if ((fpSNVgenotypes = fopen(File, "w")) == NULL)
                {
				fprintf (stderr, "Can't open \"%s\"\n", File);
				exit(-1);
                }
            }
			
        /* contains haplotypes for variable sites for every cell */
        if (doPrintSNVhaplotypes == YES)
            {
			snprintf (File, sizeof(File),"%s/%s", resultsDir, SNVhaplotypesDir);
			MAKE_DIR (File);
			snprintf (File, sizeof(File),"%s/%s/%s.%04d", resultsDir, SNVhaplotypesDir, SNVhaplotypesFile, replicate+1);
            if ((fpSNVhaplotypes = fopen(File, "w")) == NULL)
                {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
                }
            }

        /* contains reference haplotypes (before errors)  for every cell */
        if (doPrintTrueHaplotypes == YES)
            {
			snprintf (File, sizeof(File),"%s/%s", resultsDir, trueHaplotypesDir);
			MAKE_DIR (File);
			snprintf (File, sizeof(File),"%s/%s/%s.%04d", resultsDir, trueHaplotypesDir, trueHaplotypesFile, replicate+1);
            if ((fpTrueHaplotypes = fopen(File, "w")) == NULL)
                {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
                }
            }

        /* contains all genotypes (variable or invariable) for every cell */
        if (doPrintFullGenotypes == YES)
            {
			snprintf (File, sizeof(File),"%s/%s", resultsDir, fullGenotypesDir);
			MAKE_DIR (File);
			snprintf (File, sizeof(File),"%s/%s/%s.%04d", resultsDir, fullGenotypesDir, fullGenotypesFile, replicate+1);
            if ((fpFullGenotypes = fopen(File, "w")) == NULL)
                {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
                }
            }

        /* contains haplotypes for all sites for every cell */
        if (doPrintFullHaplotypes == YES)
            {
			snprintf (File, sizeof(File),"%s/%s", resultsDir, fullHaplotypesDir);
			MAKE_DIR (File);
			snprintf (File, sizeof(File),"%s/%s/%s.%04d", resultsDir, fullHaplotypesDir, fullHaplotypesFile, replicate+1);
            if ((fpFullHaplotypes = fopen(File, "w")) == NULL)
                {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
                }
            }

        /* contains reads counts and log10 normalized genotype likelihoods for every SNV and cell */
        if (doNGS == YES && doPrintVCF == YES)
            {
			snprintf (File, sizeof(File),"%s/%s", resultsDir, VCFdir);
			MAKE_DIR (File);
			snprintf (File, sizeof(File),"%s/%s/%s.%04d", resultsDir, VCFdir, VCFfile, replicate+1);
            if ((fpVCF = OpenVCF (File)) == NULL)
                {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
                }
            }

        /* contains reads counts for every SNV and cell */
        if (doPrintCATG == YES)
            {
			snprintf (File, sizeof(File),"%s/%s", resultsDir, CATGdir);
			MAKE_DIR (File);
			snprintf (File, sizeof(File),"%s/%s/%s.%04d", resultsDir, CATGdir, CATGfile, replicate+1);
            if ((fpCATG = fopen(File, "w")) == NULL)
                {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
                }
            }

		/* contains reads counts for every SNV and cell */
		if (doPrintPILEUP == YES)
			{
			snprintf (File, sizeof(File),"%s/%s", resultsDir, PILEUPdir);
			MAKE_DIR (File);
			snprintf (File, sizeof(File),"%s/%s/%s.%04d", resultsDir, PILEUPdir, PILEUPfile, replicate+1);
			if ((fpPILEUP = fopen(File, "w")) == NULL)
				{
				fprintf (stderr, "Can't open \"%s\"\n", File);
				exit(-1);
				}
			}

		/* pooled (bulk) read counts and VAF for every site */
		if (doPrintBulk == YES)
			{
			snprintf (File, sizeof(File),"%s/%s", resultsDir, bulkDir);
			MAKE_DIR (File);
			snprintf (File, sizeof(File),"%s/%s/%s.%04d", resultsDir, bulkDir, bulkFile, replicate+1);
			if ((fpBulk = fopen(File, "w")) == NULL)
				{
				fprintf (stderr, "Can't open \"%s\"\n", File);
				exit(-1);
				}
			fprintf (fpBulk,"site\tref\treadsA\treadsC\treadsG\treadsT\tdepth\talt_reads\tbulk_vaf\ttrue_alt_alleles\ttrue_vaf\n");
			}
        }
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
        if (doNGS == NO)
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
                    if (i == numCells)
                        fprintf (fp,"%s ", outCellName);
                    else if (doUserTree == NO)
                        fprintf (fp,"%s%04d ", inCellName,i+1);
                    else
                        fprintf (fp,"%-*s", stringPrecision, cellNames[i]);
                    fwrite (buf + (long long)(i-i0) * rowLen, 1, (size_t) rowLen, fp);
                    fprintf (fp,"\n");
                    }
                }
            free (buf);
            }
        else
        for (i=0; i<NUM_TAXA; i++)
            {
            if (i == numCells)
				fprintf (fp,"%s ", outCellName);
			else
				{
				if (doUserTree == NO)
					fprintf (fp,"%s%04d ", inCellName,i+1);
				else
					fprintf (fp,"%-*s", stringPrecision, cellNames[i]);
				}
			for (j=0; j<numSNVs; j++)
				{
				int	mm, pp;

				if (doNGS == NO)
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
			if (i == numCells)
				fprintf (fp,"%s ", outCellName);
			else
				{
				if (doUserTree == NO)
					fprintf (fp,"%s%04d ", inCellName,i+1);
				else
					fprintf (fp,"%-*s", stringPrecision, cellNames[i]);
				}
			for (j=0; j<numSNVs; j++)
				fprintf (fp, " %c%c", WhichMut(GetG(MATERNAL, i, SNVsites[j])),WhichMut(GetG(PATERNAL, i, SNVsites[j])));
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

/***************************** PrintSNVHaplotypes *******************************/
/* Prints obsserved/ML haplotypes at variable sites (SNVs), to a file */

void PrintSNVHaplotypes (FILE *fp, int PrintTrueVariants)
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
			fprintf (fp,"%d %d\n", NUM_TAXA+2, numSNVs);
		else
			fprintf (fp,"%d %d\n", NUM_TAXA, numSNVs);
		}
	else
		{
		if (doPrintAncestors == YES && doNGS == NO)
			fprintf (fp,"%d %d\n",2*(NUM_TAXA+2), numSNVs);
		else
			fprintf (fp,"%d %d\n",2*(NUM_TAXA), numSNVs);
		}
		
    /* site information */
	fprintf (fp, "%d", sites[0]+1);
	for (i=1; i<numSNVs; i++)
		fprintf (fp, " %d", sites[i]+1);
	fprintf (fp, "\n");

    if (alphabet == DNA)
        {
		if (doPrintIUPAChaplotypes == YES)
			{
			for (i=0; i<NUM_TAXA; i++)
				{
				/* print haplotype */
				if (i == numCells)
					fprintf (fp,"%s ", outCellName);
				else
					{
					if (doUserTree == NO)
						fprintf (fp,"%s%04d ", inCellName,i+1);
					else
						fprintf (fp,"%-*s ", stringPrecision, cellNames[i]);
					}
				for (j=0; j<numSNVs; j++)
					{
					if (doNGS == NO)
						fprintf (fp, "%c", WhichIUPAC(GetG(MATERNAL, i, sites[j]),GetG(PATERNAL, i, sites[j])));
					else
						{
						int mm, pp;
						MLPairAt (i, j, &mm, &pp);
						fprintf (fp, "%c", WhichIUPAC(mm, pp));
						}
					}
				fprintf (fp,"\n");
				}
			
			if (doPrintAncestors == YES && doNGS == NO)
				{
				fprintf (fp,"%s ", outRootCellName);
				for (j=0; j<numSNVs; j++)
					fprintf (fp, "%c", WhichIUPAC(GetG(MATERNAL, OUTGROUP_ROOT, sites[j]),GetG(PATERNAL, OUTGROUP_ROOT, sites[j])));
				fprintf (fp,"\n%s ", inRootCellName);
				for (j=0; j<numSNVs; j++)
					fprintf (fp, "%c", WhichIUPAC(GetG(MATERNAL, INGROUP_ROOT, sites[j]),GetG(PATERNAL, INGROUP_ROOT, sites[j])));
				fprintf (fp,"\n");
				}
			}
		else // print maternal and paternal DNA haplotypes
			{
			for (i=0; i<NUM_TAXA; i++)
				{
				/* print maternal haplotype */
				if (i == numCells)
					fprintf (fp,"%sm  ", outCellName);
				else
					{
					if (doUserTree == NO)
						fprintf (fp,"%s%04dm  ", inCellName, i+1);
					else
						fprintf (fp,"m%-*s ", stringPrecision, cellNames[i]);
					}
				for (j=0; j<numSNVs; j++)
					{
					if (doNGS == NO)
						fprintf (fp, "%c", WhichNuc(GetG(MATERNAL, i, SNVsites[j])));
					else
						{
						int mm, pp;
						MLPairAt (i, j, &mm, &pp);
						fprintf (fp, "%c", WhichNuc(mm));
						}
					}
				fprintf (fp,"\n");
				
				/* print paternal haplotype */
				if (i == numCells)
					fprintf (fp,"%sp  ", outCellName);
				else
					{
					if (doUserTree == NO)
						fprintf (fp,"%s%04dp  ", inCellName, i+1);
					else
						fprintf (fp,"p%-*s ", stringPrecision, cellNames[i]);
					}
				for (j=0; j<numSNVs; j++)
					{
					if (doNGS == NO)
						fprintf (fp, "%c", WhichNuc(GetG(PATERNAL, i, SNVsites[j])));
					else
						{
						int mm, pp;
						MLPairAt (i, j, &mm, &pp);
						fprintf (fp, "%c", WhichNuc(pp));
						}
					}
				fprintf (fp,"\n");
				}
			
			if (doPrintAncestors == YES && doNGS == NO)
				{
				fprintf (fp,"%sm  ", outRootCellName);
				for (j=0; j<numSNVs; j++)
					fprintf (fp, "%c", WhichNuc(GetG(MATERNAL, OUTGROUP_ROOT, sites[j])));
				fprintf (fp,"\n%sp  ", outRootCellName);
				for (j=0; j<numSNVs; j++)
					fprintf (fp, "%c", WhichNuc(GetG(PATERNAL, OUTGROUP_ROOT, sites[j])));
				fprintf (fp,"\n%sm  ", inRootCellName);
				for (j=0; j<numSNVs; j++)
					fprintf (fp, "%c", WhichNuc(GetG(MATERNAL, INGROUP_ROOT, sites[j])));
				fprintf (fp,"\n%sp  ", inRootCellName);
				for (j=0; j<numSNVs; j++)
					fprintf (fp, "%c", WhichNuc(GetG(PATERNAL, INGROUP_ROOT, sites[j])));
				fprintf (fp,"\n");
				}
			}
		}
    else  //print binary haplotypes
        {
		if (doPrintIUPAChaplotypes == YES) //print binary consensus SNV haplotypes
			{
			for (i=0; i<NUM_TAXA; i++)
				{
				if (i == numCells)
					fprintf (fp,"%s ", outCellName);
				else
					{
					if (doUserTree == NO)
						fprintf (fp,"%s%04d ", inCellName,i+1);
					else
						fprintf (fp,"%-*s ", stringPrecision, cellNames[i]);
					}
				for (j=0; j<numSNVs; j++)
					fprintf (fp, "%c", WhichConsensusBinary(GetG(MATERNAL, i, sites[j]),GetG(PATERNAL, i, sites[j])));
				fprintf (fp,"\n");
				}

			if (doPrintAncestors == YES)
				{
				fprintf (fp,"%s ", outRootCellName);
				for (j=0; j<numSNVs; j++)
					fprintf (fp, "%c", WhichConsensusBinary(GetG(MATERNAL, OUTGROUP_ROOT, sites[j]),GetG(PATERNAL, OUTGROUP_ROOT, sites[j])));
				fprintf (fp,"\n%s ", inRootCellName);
				for (j=0; j<numSNVs; j++)
					fprintf (fp, "%c", WhichConsensusBinary(GetG(MATERNAL, INGROUP_ROOT, sites[j]),GetG(PATERNAL, INGROUP_ROOT, sites[j])));
				fprintf (fp,"\n");
				}
			}
		else // print maternal and paternal binary SNV haplotypes
			{
			for (i=0; i<NUM_TAXA; i++)
				{
				/* print maternal haplotype */
				if (i == numCells)
					fprintf (fp,"%sm  ", outCellName);
				else
					{
					if (doUserTree == NO)
						fprintf (fp,"%s%04dm  ", inCellName, i+1);
					else
						fprintf (fp,"m%-*s ", stringPrecision, cellNames[i]);
					}
				for (j=0; j<numSNVs; j++)
					fprintf (fp, "%c", WhichMut(GetG(MATERNAL, i, sites[j])));
				fprintf (fp,"\n");
					
				/* print paternal haplotype */
				if (i == numCells)
					fprintf (fp,"%sp  ", outCellName);
				else
					{
					if (doUserTree == NO)
						fprintf (fp,"%s%04dp  ", inCellName, i+1);
					else
						fprintf (fp,"p%-*s ", stringPrecision, cellNames[i]);
					}
				for (j=0; j<numSNVs; j++)
					fprintf (fp, "%c", WhichMut(GetG(PATERNAL, i, sites[j])));
				fprintf (fp,"\n");
				}

			if (doPrintAncestors == YES)
				{
				fprintf (fp,"%sm  ", outRootCellName);
				for (j=0; j<numSNVs; j++)
					fprintf (fp, "%c", WhichMut(GetG(MATERNAL, OUTGROUP_ROOT, sites[j])));
				fprintf (fp,"\n%sp  ", outRootCellName);
				for (j=0; j<numSNVs; j++)
					fprintf (fp, "%c", WhichMut(GetG(PATERNAL, OUTGROUP_ROOT, sites[j])));
				fprintf (fp,"\n%sm  ", inRootCellName);
				for (j=0; j<numSNVs; j++)
					fprintf (fp, "%c", WhichMut(GetG(MATERNAL, INGROUP_ROOT, sites[j])));
				fprintf (fp,"\n%sp  ", inRootCellName);
				for (j=0; j<numSNVs; j++)
					fprintf (fp, "%c", WhichMut(GetG(PATERNAL, INGROUP_ROOT, sites[j])));
				fprintf (fp,"\n");
				}
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
            if (i == numCells)
                fprintf (fp,"%s ", outCellName);
			else
				{
				if (doUserTree == NO)
					fprintf (fp,"%s%04d ", inCellName,i+1);
				else
					fprintf (fp,"%-*s", stringPrecision, cellNames[i]);
				}
			for (j=0; j<numSites; j++)
				{
				if (doNGS == NO)
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
            if (i == numCells)
                fprintf (fp,"%s ", outCellName);
 			else
				{
				if (doUserTree == NO)
					fprintf (fp,"%s%04d ", inCellName,i+1);
				else
					fprintf (fp,"%-*s", stringPrecision, cellNames[i]);
				}
			for (j=0; j<numSites; j++)
                fprintf (fp, " %c%c", WhichMut(GetG(0, i, j)),WhichMut(GetG(1, i, j)));
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

/***************************** PrintFullHaplotypes *******************************/
/* Prints observed/ML haplotypes for all sites (variable + invariable) to a file */

void PrintFullHaplotypes (FILE *fp)
    {
    int		 i, j, ml1 = MISSING, ml2 = MISSING;
		
	if (doPrintIUPAChaplotypes == NO)
		{
		if (doPrintAncestors == YES && doNGS == NO)
			fprintf (fp,"%d %d\n",2*(NUM_TAXA+2), numSites);
		else
			fprintf (fp,"%d %d\n",2*(NUM_TAXA), numSites);
		}
	else
		{
		if (doPrintAncestors == YES && doNGS == NO)
			fprintf (fp,"%d %d\n", NUM_TAXA+2, numSites);
		else
			fprintf (fp,"%d %d\n", NUM_TAXA, numSites);
		}
		
    if (alphabet == DNA)
        {
		if (doPrintIUPAChaplotypes == YES)
			{
			for (i=0; i<NUM_TAXA; i++)
				{
	LoadFullRow (i);
				/* print IUPAC haplotype */
				if (i == numCells)
					fprintf (fp,"%s ", outCellName);
				else
					{
					if (doUserTree == NO)
						fprintf (fp,"%s%04d ", inCellName,i+1);

					else
						fprintf (fp,"%-*s ", stringPrecision, cellNames[i]);
					}
				for (j=0; j<numSites; j++)
					{
					if (doNGS == NO)
						fprintf (fp, "%c", WhichIUPAC(GetG(MATERNAL, i, j),GetG(PATERNAL, i, j)));
					else
						{
						FullMLPair (i, j, &ml1, &ml2);
						fprintf (fp, "%c", WhichIUPAC(ml1, ml2));
						}
					}
				fprintf (fp,"\n");
				}

			if (doPrintAncestors == YES && doNGS == NO)
				{
				fprintf (fp,"%s ", outRootCellName);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichIUPAC(GetG(MATERNAL, OUTGROUP_ROOT, j),GetG(PATERNAL, OUTGROUP_ROOT, j)));
				fprintf (fp,"\n%s ", inRootCellName);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichIUPAC(GetG(MATERNAL, INGROUP_ROOT, j),GetG(PATERNAL, INGROUP_ROOT, j)));
				fprintf (fp,"\n");
				}
			}
		else // print maternal and paternal DNA haplotypes
			{
			for (i=0; i<NUM_TAXA; i++)
				{
	LoadFullRow (i);
				/* print maternal haplotype */
				if (i == numCells)
					fprintf (fp,"%sm  ", outCellName);
				else
					{
					if (doUserTree == NO)
						fprintf (fp,"%s%04dm  ", inCellName, i+1);

					else
						fprintf (fp,"m%-*s ", stringPrecision, cellNames[i]);
					}
	
				for (j=0; j<numSites; j++)
					{
					if (doNGS == NO)
						fprintf (fp, "%c", WhichNuc(GetG(MATERNAL, i, j)));
					else
						{
						FullMLPair (i, j, &ml1, &ml2);
						fprintf (fp, "%c", WhichNuc(ml1));
						}
					}
				fprintf (fp,"\n");
				
				/* print paternal haplotype */
				if (i == numCells)
					fprintf (fp,"%sp  ", outCellName);
				else
					{
					if (doUserTree == NO)
						fprintf (fp,"%s%04dp  ", inCellName, i+1);
					else
						fprintf (fp,"p%-*s ", stringPrecision, cellNames[i]);
					}

				for (j=0; j<numSites; j++)
					{
					if (doNGS == NO)
						fprintf (fp, "%c", WhichNuc(GetG(PATERNAL, i, j)));
					else
						{
						FullMLPair (i, j, &ml1, &ml2);
						fprintf (fp, "%c", WhichNuc(ml2));
						}
					}
				fprintf (fp,"\n");
				}

			if (doPrintAncestors == YES && doNGS == NO)
				{
				fprintf (fp,"%sm  ", outRootCellName);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichNuc(GetG(MATERNAL, OUTGROUP_ROOT, j)));
				fprintf (fp,"\n%sp  ", outRootCellName);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichNuc(GetG(PATERNAL, OUTGROUP_ROOT, j)));
				fprintf (fp,"\n%sm  ", inRootCellName);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichNuc(GetG(MATERNAL, INGROUP_ROOT, j)));
				fprintf (fp,"\n%sp  ", inRootCellName);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichNuc(GetG(PATERNAL, INGROUP_ROOT, j)));
				fprintf (fp,"\n");
				}
			}
		}
	else  //print binary haplotypes
        {
		if (doPrintIUPAChaplotypes == YES) // print binary consensus haplotypes
			{
			for (i=0; i<NUM_TAXA; i++)
				{
				if (i == numCells)
					fprintf (fp,"%s ", outCellName);
				else
					{
					if (doUserTree == NO)
						fprintf (fp,"%s%04d ", inCellName,i+1);
					else
						fprintf (fp,"%-*s ", stringPrecision, cellNames[i]);
					}
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichConsensusBinary(GetG(MATERNAL, i, j),GetG(PATERNAL, i, j)));
				fprintf (fp,"\n");
				}
			
			if (doPrintAncestors == YES)
				{
				fprintf (fp,"%s ", outRootCellName);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichConsensusBinary(GetG(MATERNAL, OUTGROUP_ROOT, j),GetG(PATERNAL, OUTGROUP_ROOT, j)));
				fprintf (fp,"\n%s ", inRootCellName);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichConsensusBinary(GetG(MATERNAL, INGROUP_ROOT, j),GetG(PATERNAL, INGROUP_ROOT, j)));
				fprintf (fp,"\n");
				}
			}
		else // print maternal and paternal binary haplotypes
			{
			for (i=0; i<NUM_TAXA; i++)
				{
				/* print maternal haplotype */
				if (i == numCells)
					fprintf (fp,"%sm  ", outCellName);
				else
					{
					if (doUserTree == NO)
						fprintf (fp,"%s%04dm  ", inCellName, i+1);
					else
						fprintf (fp,"m%-*s ", stringPrecision, cellNames[i]);
					}
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichMut(GetG(MATERNAL, i, j)));
				fprintf (fp,"\n");
					
				/* print paternal haplotype */
				if (i == numCells)
					fprintf (fp,"%sp  ", outCellName);
				else
					{
					if (doUserTree == NO)
						fprintf (fp,"%s%04dp  ", inCellName, i+1);
					else
						fprintf (fp,"p%-*s ", stringPrecision, cellNames[i]);
					}
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichMut(GetG(PATERNAL, i, j)));
				fprintf (fp,"\n");
				}
			
			if (doPrintAncestors == YES)
				{
				fprintf (fp,"%sm  ", outRootCellName);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichMut(GetG(MATERNAL, OUTGROUP_ROOT, j)));
				fprintf (fp,"\n%sp  ", outRootCellName);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichMut(GetG(PATERNAL, OUTGROUP_ROOT, j)));
				fprintf (fp,"\n%sm  ", inRootCellName);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichMut(GetG(MATERNAL, INGROUP_ROOT, j)));
				fprintf (fp,"\n%sp  ", inRootCellName);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichMut(GetG(PATERNAL, INGROUP_ROOT, j)));
				fprintf (fp,"\n");
				}
			}
        }
    }

/***************************** PrintTrueFullHaplotypes *******************************/
/* Prints true haplotypes for all sites (variable + invariable) to a file */

void PrintTrueFullHaplotypes (FILE *fp)
	{
	int		 i, j;
	
	if (doPrintIUPAChaplotypes == NO)
		{
		if (doPrintAncestors == YES)
			fprintf (fp,"%d %d\n",2*(NUM_TAXA+2), numSites);
		else
			fprintf (fp,"%d %d\n",2*(NUM_TAXA), numSites);
		}
	else
		{
		if (doPrintAncestors == YES)
			fprintf (fp,"%d %d\n", NUM_TAXA+2, numSites);
		else
			fprintf (fp,"%d %d\n", NUM_TAXA, numSites);
		}
	
	if (alphabet == DNA)
		{
		if (doPrintIUPAChaplotypes == YES)
			{
			for (i=0; i<NUM_TAXA; i++)
				{
				/* print IUPAC haplotype */
				if (i == numCells)
					fprintf (fp,"%s ", outCellName);
				else
					{
					if (doUserTree == NO)
						fprintf (fp,"%s%04d ", inCellName,i+1);
					else
						fprintf (fp,"%-*s ", stringPrecision, cellNames[i]);
					}
				for (j=0; j<numSites; j++)
					{
					fprintf (fp, "%c", WhichIUPAC(GetG(MATERNAL, i, j),GetG(PATERNAL, i, j)));
					}
				fprintf (fp,"\n");
				}
			
			if (doPrintAncestors == YES)
				{
				fprintf (fp,"%s ", outRootCellName);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichIUPAC(GetG(MATERNAL, OUTGROUP_ROOT, j),GetG(PATERNAL, OUTGROUP_ROOT, j)));
				fprintf (fp,"\n%s ", inRootCellName);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichIUPAC(GetG(MATERNAL, INGROUP_ROOT, j),GetG(PATERNAL, INGROUP_ROOT, j)));
				fprintf (fp,"\n");
				}
			}
		else // print maternal and paternal DNA haplotypes
			{
			for (i=0; i<NUM_TAXA; i++)
				{
				/* print maternal haplotype */
				if (i == numCells)
					fprintf (fp,"%sm  ", outCellName);
				else
					{
					if (doUserTree == NO)
						fprintf (fp,"%s%04dm  ", inCellName, i+1);
					else
						fprintf (fp,"m%-*s ", stringPrecision, cellNames[i]);
					}
				for (j=0; j<numSites; j++)
					{
					fprintf (fp, "%c", WhichNuc(GetG(MATERNAL, i, j)));
					}
				fprintf (fp,"\n");
				
				/* print paternal haplotype */
				if (i == numCells)
					fprintf (fp,"%sp  ", outCellName);
				else
					{
					if (doUserTree == NO)
						fprintf (fp,"%s%04dp  ", inCellName, i+1);
					else
						fprintf (fp,"p%-*s ", stringPrecision, cellNames[i]);
					}
				for (j=0; j<numSites; j++)
					{
					fprintf (fp, "%c", WhichNuc(GetG(PATERNAL, i, j)));
					}
				fprintf (fp,"\n");
				}
			
			if (doPrintAncestors == YES)
				{
				fprintf (fp,"%sm  ", outRootCellName);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichNuc(GetG(MATERNAL, OUTGROUP_ROOT, j)));
				fprintf (fp,"\n%sp  ", outRootCellName);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichNuc(GetG(PATERNAL, OUTGROUP_ROOT, j)));
				fprintf (fp,"\n%sm  ", inRootCellName);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichNuc(GetG(MATERNAL, INGROUP_ROOT, j)));
				fprintf (fp,"\n%sp  ", inRootCellName);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichNuc(GetG(PATERNAL, INGROUP_ROOT, j)));
				fprintf (fp,"\n");
				}
			}
		}
	else  //print binary haplotypes
		{
		if (doPrintIUPAChaplotypes == YES) // print binary consensus haplotypes
			{
			for (i=0; i<NUM_TAXA; i++)
				{
				if (i == numCells)
					fprintf (fp,"%s ", outCellName);
				else
					{
					if (doUserTree == NO)
						fprintf (fp,"%s%04d ", inCellName,i+1);
					else
						fprintf (fp,"%-*s ", stringPrecision, cellNames[i]);
					}
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichConsensusBinary(GetG(MATERNAL, i, j),GetG(PATERNAL, i, j)));
				fprintf (fp,"\n");
				}

			if (doPrintAncestors == YES)
				{
				fprintf (fp,"%s ", outRootCellName);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichConsensusBinary(GetG(MATERNAL, OUTGROUP_ROOT, j),GetG(PATERNAL, OUTGROUP_ROOT, j)));
				fprintf (fp,"\n%s ", inRootCellName);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichConsensusBinary(GetG(MATERNAL, INGROUP_ROOT, j),GetG(PATERNAL, INGROUP_ROOT, j)));
				fprintf (fp,"\n");
				}
			}
		else // print maternal and paternal binary haplotypes
			{
			for (i=0; i<NUM_TAXA; i++)
				{
				/* print maternal haplotype */
				if (i == numCells)
					fprintf (fp,"%sm  ", outCellName);
				else
					{
					if (doUserTree == NO)
						fprintf (fp,"%s%04dm  ", inCellName, i+1);
					else
						fprintf (fp,"m%-*s ", stringPrecision, cellNames[i]);
					}
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichMut(GetG(MATERNAL, i, j)));
				fprintf (fp,"\n");
				
				/* print paternal haplotype */
				if (i == numCells)
					fprintf (fp,"%sp  ", outCellName);
				else
					{
					if (doUserTree == NO)
						fprintf (fp,"%s%04dp  ", inCellName, i+1);
					else
						fprintf (fp,"p%-*s ", stringPrecision, cellNames[i]);
					}
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichMut(GetG(PATERNAL, i, j)));
				fprintf (fp,"\n");
				}
			
			if (doPrintAncestors == YES)
				{
				fprintf (fp,"%sm  ", outRootCellName);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichMut(GetG(MATERNAL, OUTGROUP_ROOT, j)));
				fprintf (fp,"\n%sp  ", outRootCellName);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichMut(GetG(PATERNAL, OUTGROUP_ROOT, j)));
				fprintf (fp,"\n%sm  ", inRootCellName);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichMut(GetG(MATERNAL, INGROUP_ROOT, j)));
				fprintf (fp,"\n%sp  ", inRootCellName);
				for (j=0; j<numSites; j++)
					fprintf (fp, "%c", WhichMut(GetG(PATERNAL, INGROUP_ROOT, j)));
				fprintf (fp,"\n");
				}
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
/*	Prints the settings used in the simulations */
void PrintCommandLine (FILE *fp, int argc,char **argv)
    {
	int i;

	if (readingParameterFile == NO)
		fprintf (fp, "Command line arguments (some might not affect) = %s", CommandLine);
	else
		{
		fprintf (fp, "Parameter file [\"%s\"] arguments (some might not affect) = ", parameterFile);

		/* Coalescent */
		fprintf (fp, " -n%d -s%d -l%d -e%d -g%2.1e -K%2.1e -L%2.1e -h%d", numDataSets, numCells, numSites, N, growthRate, birthRate, deathRate, numPeriods);
		for (i=1; i<=numPeriods; i++)
			fprintf (fp, " %d %d %d", Nbegin[i], Nend[i], cumDuration[i]-cumDuration[i-1]);

		/* Post-coalescent */
		fprintf (fp, " -k%2.1e", transformingBranchLengthRatio);
		fprintf (fp, " -q%2.1e", healthyTipBranchLengthRatio);
		fprintf (fp, " -i%6.4f", alphaBranches);
		if (doTreeRateSwitches == YES)
			{
			fprintf (fp, " -N%d", numTreeRateSwitches);
			for (i=0; i<numTreeRateSwitches; i++)
				fprintf (fp, " %4.2f", treeRateSwitchMultiplier[i]);
			}
		
		/* Mutation model */
		fprintf (fp, " -b%d", alphabet);
		fprintf (fp, " -c%2.1e", SNPrate);
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
		fprintf (fp, " -f%3.2f %3.2f %3.2f %3.2f", freq[0], freq[1], freq[2], freq[3]);
		fprintf (fp, " -t%6.4f", titv);
		fprintf (fp, " -a%6.4f", alphaSites);
		fprintf (fp, " -r%3.2f %3.2f %3.2f %3.2f  %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f", Mij[0][0], Mij[0][1], Mij[0][2], Mij[0][3]
		, Mij[1][0], Mij[1][1], Mij[1][2], Mij[1][3],  Mij[2][0], Mij[2][1], Mij[2][2], Mij[2][3],  Mij[3][0], Mij[3][1], Mij[3][2], Mij[3][3]);

		/* scWGA */
		fprintf (fp, " -D%2.1e", fixedADOrate);
		fprintf (fp, " -P%2.1e %2.1e", meanADOcell, varADOcell);
		fprintf (fp, " -Q%2.1e %2.1e", meanADOsite, varADOsite);
		fprintf (fp, " -I%2.1e %2.1e", meanAllelicImbalance, varAllelicImbalance);
		fprintf (fp, " -A%2.1e %2.1e %d", meanAmplificationError, varAmplificationError, simulateOnlyTwoTemplates);
		fprintf (fp, " -B%2.1e %2.1e", meanDoubletRate, varDoubletRate);

		/* GENOTYPING */
		fprintf (fp, " -G%2.1e %2.1e", meanGenotypingError, varGenotypingError);

		/* NGS */
		fprintf (fp, " -C%3.2f", coverage);
		fprintf (fp, " -V%6.4f", alphaCoverage);
		if (alphaSiteCoverage > 0.0)
			fprintf (fp, " -O%6.4f", alphaSiteCoverage);
		if (alphaCellCoverage > 0.0)
			fprintf (fp, " -J%6.4f", alphaCellCoverage);
		fprintf (fp, " -R%2.1e", haploidCoverageReduction);
		fprintf (fp, " -E%2.1e", sequencingError);
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
		if (doPrintPILEUP == YES)
			fprintf (fp, " -%d", 9);
		if (doPrintBulk == YES)
			fprintf (fp, " -%d", 0);
		if (doPrintVCF == NO)
			fprintf (fp, " -%c", 'Z');
		if (doZipVCF == YES)
			fprintf (fp, " -%c", 'z');
		if (minVAFvcf > 0.0)
			fprintf (fp, " -M%g", minVAFvcf);
		if (doPrintTrueHaplotypes == YES)
			fprintf (fp, " -%c", 'Y');
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
		if (doTumorNames == YES)
			fprintf (fp, " -%c", 'W');
		fprintf (fp, " -y%d -#%ld", noisy, originalSeed);
		}
	}


/***************************** PrintUsage *******************************/
/* Prints a short description of program usage  */
void PrintUsage(FILE *fp)
{
	PrintHeader(fp);
    fprintf (fp, "\n%s generates a coalescent tree and simulates a sample of diploid genomes from somatic cells (no recombination), together with a (healthy) cell as outgroup.", PROGRAM_NAME);
	fprintf (fp, "\n\nUsage: %s [-FTXT -n# -s# -l# -e# -g# -h# (# # #) -k# -q# -i# -b# -c# -u# -d# -H# -j# -S# (# #) -m# -p# -w# -f# # # # -t# -a# -r# # # # # # # # # # # # # # # # -G# -C# -V# -I# -D# -P# -Q# -R# -A# # # -E# -B# -X# # # # # # # # # # # # # # # # -1 -2 -3 -4 -5 -6 -7 -8 -9 -Y -v -x -oTXT -TTXT -UTXT -W# -y# -## -?]", PROGRAM_NAME_LOWERCASE);
	fprintf (fp,"\n");

	fprintf (fp,"\n-F: user-defined parameter file (e.g. -Fmyparameters)");
	fprintf (fp,"\n-n: number of replicates (e.g. -n1000)");
    fprintf (fp,"\n-s: sample size (# cells) (e.g. -s8)");
    fprintf (fp,"\n-l: number of sites (e.g. -l500)");
    fprintf (fp,"\n-e: effective population size (e.g. -e1000)");
	fprintf (fp,"\n-g: exponential growth rate (e.g. -g1e-5)");
	fprintf (fp,"\n-K: cell birth rate (e.g. -K1)");
	fprintf (fp,"\n-L: cell death rate (e.g. -L0.001)");
   	fprintf (fp,"\n-h: number of demographic periods followed by effective population size at the beginning and at the end of the period, ");
    fprintf (fp, "\n    and the duration of the period in generations. (e.g. -d1 100 200 30000) (e.g. -d2 100 100 40000 1000 2000 20000)");
	fprintf (fp,"\n-k: transforming branch length (e.g. -u1e-2)");
    fprintf (fp,"\n-q: healthy tip branch length (e.g. -q1e-6)");
    fprintf (fp,"\n-i: shape of the gamma distribution for rate variation among lineages (e.g. -i0.5)");
	fprintf (fp,"\n-N: number of tree rate switches followed by the rate multipliers  (e.g. -N2 10 0.2)");
	fprintf (fp,"\n-b: alphabet [0:binary 1:DNA] (e.g. -b0)");
	fprintf (fp,"\n-c: germline SNP rate (e.g. -c1e-5)");
	fprintf (fp,"\n-u: somatic mutation rate (e.g. -u1e-6)");
	fprintf (fp,"\n-d: deletion rate (e.g. -d1e-5)");
	fprintf (fp,"\n-H: CNLOH rate (e.g. -H1e-5)");
    fprintf (fp,"\n-j: fixed number of mutations (e.g. -j100)");
    fprintf (fp,"\n-S: trinucleotide genetic signatures (e.g. -S2 3 0.8 13 0.2)");
    fprintf (fp,"\n-m: alternative mutation model [0:ISMhap 1:Mk2 2:finiteDNA] [default mutation model is ISM diploid] (e.g. -m2)");
    fprintf (fp,"\n-p: proportion of alternative model sites (e.g. -p0.1)");
    fprintf (fp,"\n-w: alternative/default model relative mutation rate (e.g. -w1)");
	fprintf (fp,"\n-f: nucleotide base frequencies (e.g. -f0.4 0.3 0.2 0.1)");
	fprintf (fp,"\n-t: transition/transversion ratio (e.g. -t2)");
    fprintf (fp,"\n-a: shape of the gamma distribution for rate variation among sites (e.g. -a0.2)");
	fprintf (fp,"\n-r: mutation matrix ACGT x ACGT (e.g. -r0 1 2 3 1 0 4 5 2 4 0 1 3 5 1 0)");

	fprintf (fp,"\n-D: allelic dropout (e.g. -D0.1)");
	fprintf (fp,"\n-P: allelic dropout variation among cells (e.g. -P1)");
	fprintf (fp,"\n-Q: allelic dropout variation among sites (e.g. -Q1)");
	fprintf (fp,"\n-I: maternal allelic imbalance (e.g. -I0.5)");
	fprintf (fp,"\n-A: amplification error (e.g. -A0.1 0.01 0)");
	fprintf (fp,"\n-B: doublet rate per cell (e.g. -B0.1)");

	fprintf (fp,"\n-G: genotyping error [no read counts]: mean and variance (e.g. -G0.01 0.001)");

	fprintf (fp,"\n-C: sequencing coverage [read counts, not individual reads] (e.g. -C60)");
	fprintf (fp,"\n-V: sequencing coverage overdispersion (e.g. -V5)");
	fprintf (fp,"\n-O: coverage overdispersion between sites, shared by all cells (e.g. -O10)");
	fprintf (fp,"\n-J: coverage overdispersion between cells, shared by all sites (e.g. -J10)");
	fprintf (fp,"\n-R: ADO/deletion haploid coverage reduction (e.g. -R0.5)");
	fprintf (fp,"\n-E: sequencing error (e.g. -E0.001)");
	fprintf (fp,"\n-X: error matrix ACGT x ACGT (e.g. -X0 1 1 1 1 0 1 1 1 1 0 1 1 1 1 0)");

    fprintf (fp,"\n-1: print SNV genotypes to a file (e.g. -1)");
    fprintf (fp,"\n-2: print SNV haplotypes to a file (e.g. -2)");
    fprintf (fp,"\n-3: print full genotypes to a file (e.g. -3)");
    fprintf (fp,"\n-4: print full haplotypes to a file (e.g. -4)");
    fprintf (fp,"\n-5: print ancestral genotypes (e.g. -5)");
    fprintf (fp,"\n-6: print trees to a file (e.g. -6)");
    fprintf (fp,"\n-7: print times to a file (e.g. -7)");
    fprintf (fp,"\n-8: print read counts in CATG format (e.g. -8)");
	fprintf (fp,"\n-9: print read counts in PILEUP format (e.g. -9)");
	fprintf (fp,"\n-0: print pooled (bulk) read counts and VAF per site (e.g. -0)");
	fprintf (fp,"\n-Z: do not print the VCF file (needs -8, -9 or -0) (e.g. -Z)");
	fprintf (fp,"\n-z: write the VCF compressed with gzip, as vcf.gz (e.g. -z)");
	fprintf (fp,"\n-M: minimum pooled VAF in the read counts of the sites listed in the VCF (e.g. -M0.02)");
    fprintf (fp,"\n-Y: print true SNV haplotypes to a file (e.g. -Y)");
	fprintf (fp,"\n-v: print replicates in individual folders (e.g. -v)");
	fprintf (fp,"\n-x: print consensus/IUPAC haplotypes (e.g. -x)");
	fprintf (fp,"\n-o: results folder name (e.g. -oresultsFolder)");
	fprintf (fp,"\n--amp-cell-var V: the mean amplification error (-A) varies among cells, beta with the mean of -A and variance V (e.g. -A0.05 0.002 0 --amp-cell-var 0.0007)");
	fprintf (fp,"\n--cov-ado-corr RHO: correlation (-1..1) between the coverage (-J) and the ADO rate (-P) of the cells (e.g. --cov-ado-corr -0.5)");
	fprintf (fp,"\n--cov-doublet-corr RHO: correlation (-1..1) between the coverage (-J) and being a doublet (-B)");
	fprintf (fp,"\n--ado-doublet-corr RHO: correlation (-1..1) between the ADO rate (-P) and being a doublet (-B)");
	fprintf (fp,"\n--version: print the version and exit");
	fprintf (fp,"\n--check: check the options and show what would run, without simulating");
	fprintf (fp,"\n--summary: write summary.tsv, one row of metrics per replicate");
	fprintf (fp,"\n--sweep \"KEY=V1,V2 ...\": run the grid of values (e.g. --sweep \"C=2,5,10 E=0.001,0.01\") and write sweep.tsv with mean and standard error of the metrics; --sweep-keep keeps the runs");
	fprintf (fp,"\n--bcf: write the genotype likelihoods as a BCF file (vcf.bcf) instead of a VCF file");
	fprintf (fp,"\n--min-alt-reads N: a site is written (VCF, catg, pileup) only if a cell has a non-reference call with at least N reads of that allele (e.g. --min-alt-reads 2)");
	fprintf (fp,"\n--min-alt-vaf F: ... and at least a fraction F of the reads of that cell (e.g. --min-alt-vaf 0.3)");
	fprintf (fp,"\n--ado-block L: allelic dropout in contiguous segments of mean length L sites (a Markov chain along each allele, with the ADO rate of -D/-P); --ado-origin O: position of the first site on the chromosome; --ado-seed S: seed of the segments");
	fprintf (fp,"\n--coverage W: mean read depth of every cell over the sites that are not variable sites, in windows of W sites (coverage.tsv; needs -C)");
	fprintf (fp,"\n--no-outgroup: simulate the tumor only: no healthy outgroup cell in the outputs (the reference is still the healthy root genome)");
	fprintf (fp,"\n--fast: faster read simulation (same distributions, different random numbers; see the manual)");
	fprintf (fp,"\n--config FILE: read the options from a JSON or YAML file");
	fprintf (fp,"\n-T: user-tree file (e.g. -Tusertree)");
	fprintf (fp,"\n-U: user genome (e.g. -Uusergenome)");
	fprintf (fp,"\n-W: use tumor nomenclature for cells (e.g. -W)");
	fprintf (fp,"\n-y: noisy (e.g. -y1)");
    fprintf (fp,"\n      = 0: does not print anything");
    fprintf (fp,"\n      = 1:  + simulation summary");
    fprintf (fp,"\n      = 2:  + replicate information");
    fprintf (fp,"\n      = 3: + calculation status and event information");
    fprintf (fp,"\n-#: seed (e.g. -#37864287)");
    fprintf (fp,"\n-?: Print help");
	
    fprintf (fp,"\n\nDefaults: ");
    PrintDefaults (fp);
	fprintf (fp, "\n--------------------------------------------------------------------------------------------------------\n");
	fprintf (fp,"\n...now exiting to system...\n\n");
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
	fprintf (fp,"\n-Y: print true haplotypes to a file =  %c", doPrintTrueHaplotypes);
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
