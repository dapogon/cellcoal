#include "cellcoal.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "definitions.h"

/* Output of every replicate and of the summary: messages, genotype files, summary.tsv and its TOTAL row */

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

/* means and variances over replicates */
void ComputeSummaryStatistics (void)
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
void RepMarker (FILE *fp, int dataSetNum)
	{
	if (blocksOn == 0)
		fprintf (fp, "[#%d]\n", dataSetNum+1);
	}

/* per-replicate messages and output of genotypes, haplotypes and reads */
void PrintReplicateOutputs (int dataSetNum)
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

void StoreBlockRow (int mu, int del, int cnloh, int snvTrue, int snvObs, const double *metric, const long long *inv)
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
void WriteBlocksTotal (void)
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

void ReportProgress (int dataSetNum)
	{
	fprintf (stdout, "\rReplicate #%3d/%d ", dataSetNum+1, numDataSets);
	fflush (stdout);
	if (!isatty (fileno (stdout)))		/* if stdout is redirected elsewhere, print also to stderr */
		{
		fprintf (stderr, "\rreplicate #%3d/%d", dataSetNum+1, numDataSets);
		fflush (stderr);
		}
	}

/* a row of summary.tsv (and of the table of the blocks) */
void WriteSummaryRow (int dataSetNum)
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
