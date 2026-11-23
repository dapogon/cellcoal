/*
	CellCoal - coalescent tree: simulation, user trees, rate variation
	See cellcoal.h for the shared declarations. This file was split from cellcoal.c.
*/
#include "cellcoal.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "definitions.h"

/* Structured coalescent with colonizations (--deme-n, --deme-size, --deme-source, --deme-time, --deme-founders; group --demes). The sampled cells
   belong to regions (demes) of a tumor: the primary tumor (deme 1) and its metastases. There is no exchange of cells between the demes: a deme is
   founded once, at a time (in N generations before the present), by a number of cells of its source deme, and after that it grows apart. Looking
   back in time, two lineages of the same deme coalesce at the rate 1/size (size relative to N), until the time of the colonization of the deme; then
   its lineages go through the bottleneck of the founders (they coalesce until at most F remain) and continue in the source deme. */
char	*demeN = "", *demeSize = "", *demeSource = "", *demeTime = "", *demeFounders = "";
int		numDemes = 0, *cellDeme = NULL;
static double	*demeSz = NULL, *demeT = NULL;		/* relative sizes and times of colonization (the first deme, the primary, has no time) */
static int		*demeSrc = NULL, *demeFnd = NULL;	/* source deme (0-based) and number of founders */

static int ParseNumbers (const char *s, double *v, int max)
	{
	int		n = 0;
	char	*copy = strdup (s), *tok, *save = NULL;

	for (tok = strtok_r (copy, ":", &save); tok; tok = strtok_r (NULL, ":", &save))
		{
		char *e;

		if (n >= max)
			{
			free (copy);
			return -1;
			}
		v[n] = strtod (tok, &e);
		if (e == tok || *e != '\0')
			{
			free (copy);
			return -1;
			}
		n++;
		}
	free (copy);
	return n;
	}

void DemesSetup (void)
	{
	double	v[256];
	int		k, n, K, total = 0;

	if (strlen (demeN) == 0)
		{
		if (strlen (demeSize) > 0 || strlen (demeSource) > 0 || strlen (demeTime) > 0 || strlen (demeFounders) > 0)
			{
			fprintf (stderr, "\nERROR: --deme-size, --deme-source, --deme-time and --deme-founders need --deme-n (the number of cells of every deme)\n");
			exit (-1);
			}
		return;
		}
	if (doUserTree == YES || doDemographics == YES || doExponential == YES || doOthsukiInnanCoal == YES)
		{
		fprintf (stderr, "\nERROR: the structured coalescent (--demes) cannot be combined with a user tree, growth or demographic periods\n");
		exit (-1);
		}
	K = ParseNumbers (demeN, v, 64);
	if (K < 2)
		{
		fprintf (stderr, "\nERROR: --deme-n needs the number of cells of at least two demes, separated by : (for example 20:10:10)\n");
		exit (-1);
		}
	for (k=0; k<K; k++)
		{
		if (v[k] < 0.0 || v[k] != floor (v[k]) || (v[k] < 1.0 && k > 0))
			{
			fprintf (stderr, "\nERROR: --deme-n: the number of cells of a deme has to be a whole number of at least 1 (0 is allowed only for the first deme, a primary tumor that is not sampled)\n");
			exit (-1);
			}
		total += (int) v[k];
		}
	if (total < 2)
		{
		fprintf (stderr, "\nERROR: --deme-n: at least two cells are needed\n");
		exit (-1);
		}
	cellDeme = (int*) Checked (calloc (total, sizeof(int)), "Could not allocate the demes\n");
	for (k=0, n=0; k<K; k++)
		{
		int c, m = (int) v[k];

		for (c=0; c<m; c++)
			cellDeme[n++] = k;
		}
	demeSz = (double*) Checked (calloc (K, sizeof(double)), "Could not allocate the demes\n");
	demeT = (double*) Checked (calloc (K, sizeof(double)), "Could not allocate the demes\n");
	demeSrc = (int*) Checked (calloc (K, sizeof(int)), "Could not allocate the demes\n");
	demeFnd = (int*) Checked (calloc (K, sizeof(int)), "Could not allocate the demes\n");
	for (k=0; k<K; k++)
		{
		demeSz[k] = 1.0;
		demeFnd[k] = 1;
		}
	if (strlen (demeSize) > 0)
		{
		n = ParseNumbers (demeSize, v, 64);
		if (n != K)
			{
			fprintf (stderr, "\nERROR: --deme-size needs %d sizes (one per deme), separated by :\n", K);
			exit (-1);
			}
		for (k=0; k<K; k++)
			{
			if (v[k] <= 0.0)
				{
				fprintf (stderr, "\nERROR: --deme-size: the sizes have to be greater than 0\n");
				exit (-1);
				}
			demeSz[k] = v[k];
			}
		}
	if (strlen (demeTime) == 0)
		{
		fprintf (stderr, "\nERROR: --deme-time is needed: the time of colonization of the demes 2 to %d, in N generations before the present (%d numbers separated by :)\n", K, K - 1);
		exit (-1);
		}
	n = ParseNumbers (demeTime, v, 64);
	if (n != K - 1)
		{
		fprintf (stderr, "\nERROR: --deme-time needs %d times (demes 2 to %d), separated by :\n", K - 1, K);
		exit (-1);
		}
	for (k=1; k<K; k++)
		{
		if (v[k-1] <= 0.0)
			{
			fprintf (stderr, "\nERROR: --deme-time: the times have to be greater than 0\n");
			exit (-1);
			}
		demeT[k] = v[k-1];
		demeSrc[k] = 0;
		}
	if (strlen (demeSource) > 0)
		{
		n = ParseNumbers (demeSource, v, 64);
		if (n != K - 1)
			{
			fprintf (stderr, "\nERROR: --deme-source needs %d demes (the source of the demes 2 to %d), separated by :\n", K - 1, K);
			exit (-1);
			}
		for (k=1; k<K; k++)
			{
			int sdeme = (int) v[k-1];

			if (v[k-1] != sdeme || sdeme < 1 || sdeme > K || sdeme == k + 1)
				{
				fprintf (stderr, "\nERROR: --deme-source: the source of the deme %d has to be another deme, 1 to %d\n", k + 1, K);
				exit (-1);
				}
			demeSrc[k] = sdeme - 1;
			}
		}
	if (strlen (demeFounders) > 0)
		{
		n = ParseNumbers (demeFounders, v, 64);
		if (n != K - 1)
			{
			fprintf (stderr, "\nERROR: --deme-founders needs %d numbers (the founder cells of the demes 2 to %d), separated by :\n", K - 1, K);
			exit (-1);
			}
		for (k=1; k<K; k++)
			{
			if (v[k-1] < 1.0 || v[k-1] != floor (v[k-1]))
				{
				fprintf (stderr, "\nERROR: --deme-founders: a deme is founded by a whole number of cells, at least 1\n");
				exit (-1);
				}
			demeFnd[k] = (int) v[k-1];
			}
		}
	/* a deme can only be colonized from a deme that exists then: the source is older (and the chain of sources reaches the first deme) */
	for (k=1; k<K; k++)
		if (demeSrc[k] != 0 && !(demeT[demeSrc[k]] > demeT[k]))
			{
			fprintf (stderr, "\nERROR: the deme %d is colonized (time %g) from the deme %d, which is colonized later (time %g): the source has to be older\n", k + 1, demeT[k], demeSrc[k] + 1, demeT[demeSrc[k]]);
			exit (-1);
			}
	numDemes = K;
	numCells = total;			/* -s is ignored */
	}

/* demes.tsv: the deme of every cell */
void WriteDemes (void)
	{
	char	file[2*MAX_NAME+16];
	FILE	*fp;
	int		i;

	if (numDemes == 0)
		return;
	snprintf (file, sizeof(file), "%s/demes.tsv", resultsDir);
	if ((fp = fopen (file, "w")) == NULL)
		return;
	fprintf (fp, "cell\tdeme\n");
	for (i=0; i<numCells; i++)
		fprintf (fp, "%s%04d\t%d\n", inCellName, i + 1, cellDeme[i] + 1);
	fclose (fp);
	}

/************************* MakeCoalescenceTree  ************************/
/*	Builds a genealogy  under the coalescent with exponential growth without recombination */

/* the coalescent time is infinite when the growth rate is too negative */
static void FatalInfiniteCoalescentTime (double t)
	{
	fprintf (stderr, "\nERROR: Coalescent time (%f) is infinite ", t);
	fprintf (stderr, "\n       This might suggest that the exponential growth rate is too negative");
	fprintf (stderr, "\n       and the coalescent time is therefore infinite.");
	fprintf (stderr, "\n       Try a smaller value");
	exit (-1);
	}

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
        fprintf (stderr, "Could not allocate nodes (%zu)\n", numNodes * sizeof(TreeNode));
        exit (-1);
        }
    activeNodes = (int *) calloc(numNodes, sizeof(int));
    if (!activeNodes)
        {
        fprintf (stderr, "Could not allocate activeNodes (%zu)\n", numNodes * sizeof(int));
        exit (-1);
        }
	
    /* the nodes come zeroed (calloc): only the names need space */
    for (i=0; i<numNodes; i++)
        {
		treeNodes[i].name = (char*) Checked (calloc(MAX_NAME, sizeof(char)), "Could not allocate the treeNodes[i].name structure\n");
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
	
    /* structured coalescent with colonizations: coalescences inside the demes that exist, and the colonization events */
    if (numDemes > 0)
        {
        int		K = numDemes, **lin, *cnt, *done, total = numCells, k, j, a, b;
        double	*coal, h, u, tNext, wait;
        int		kNext;

        lin = (int**) Checked (calloc (K, sizeof(int*)), "Could not allocate the demes\n");
        cnt = (int*) Checked (calloc (2 * K, sizeof(int)), "Could not allocate the demes\n");
        done = cnt + K;
        coal = (double*) Checked (calloc (K, sizeof(double)), "Could not allocate the demes\n");
        for (k=0; k<K; k++)
            lin[k] = (int*) Checked (calloc (numCells, sizeof(int)), "Could not allocate the demes\n");
        for (i=0; i<numCells; i++)
            lin[cellDeme[i]][cnt[cellDeme[i]]++] = i;
        currentTime = 0.0;
        while (total > 1)
            {
            h = 0.0;
            tNext = -1.0;
            kNext = -1;
            for (k=0; k<K; k++)		/* the demes that exist now: the ones not yet colonized (looking back) */
                {
                coal[k] = done[k] ? 0.0 : cnt[k] * (cnt[k] - 1) / 2.0 / demeSz[k];
                h += coal[k];
                if (k > 0 && !done[k] && (kNext < 0 || demeT[k] * N < tNext))
                    {
                    kNext = k;
                    tNext = demeT[k] * N;
                    }
                }
            wait = (h > 0.0) ? RandomExponential (h, seed) * N : -1.0;
            if (kNext >= 0 && (wait < 0.0 || currentTime + wait >= tNext))
                {
                /* colonization of the deme kNext: its lineages go through the bottleneck of the founders and continue in the source deme */
                int		F = demeFnd[kNext], m = cnt[kNext], *slot, c, src = demeSrc[kNext];
                int		ncoal = 0;

                currentTime = tNext;
                slot = (int*) malloc ((m > 0 ? m : 1) * sizeof(int));
                for (c=0; c<m; c++)
                    slot[c] = (int) RandomUniformTo (F, seed);	/* every lineage descends from one of the F founders, at random */
                for (c=0; c<m; c++)					/* the lineages of the same founder coalesce at the time of the colonization */
                    {
                    int d, first = -1;

                    for (d=0; d<c; d++)
                        if (slot[d] == slot[c])
                            {
                            first = d;
                            break;
                            }
                    if (first >= 0 || slot[c] < 0)
                        continue;
                    /* c is the first lineage of its founder: coalesce the following ones into it */
                    for (d=c+1; d<m; d++)
                        if (slot[d] == slot[c])
                            {
                            double t = currentTime * (1.0 + 1e-9 * (++ncoal));
                            TreeNode *pp = treeNodes + lin[kNext][c], *qq = treeNodes + lin[kNext][d], *rr = treeNodes + nextAvailableNode;

                            numCA++;
                            rr->index = nextAvailableNode;
                            rr->left = pp;
                            rr->right = qq;
                            pp->anc = rr;
                            qq->anc = rr;
                            rr->time = t;
                            pp->length = rr->time - pp->time;
                            qq->length = rr->time - qq->time;
                            pp->branchLength = pp->length * mutationRate;
                            qq->branchLength = qq->length * mutationRate;
                            lin[kNext][c] = nextAvailableNode;
                            r = rr;
                            nextAvailableNode++;
                            total--;
                            slot[d] = -1;		/* absorbed */
                            if (nextAvailableNode >= numNodes)
                                {
                                fprintf (stderr, "\n\nERROR: Too many nodes!\n");
                                exit (-1);
                                }
                            }
                    }
                for (c=0; c<m; c++)
                    if (slot[c] >= 0)
                        lin[src][cnt[src]++] = lin[kNext][c];
                free (slot);
                cnt[kNext] = 0;
                done[kNext] = 1;
                continue;
                }
            if (h <= 0.0)
                {
                fprintf (stderr, "\n\nERROR: the structured coalescent cannot finish (check --deme-source and --deme-time)\n");
                exit (-1);
                }
            currentTime += wait;
            u = RandomUniform (seed) * h;
            for (k=0; k<K-1; k++)
                {
                if (u < coal[k])
                    break;
                u -= coal[k];
                }
            while (coal[k] <= 0.0 && k > 0)
                k--;
            a = (int) (cnt[k] * RandomUniform (seed));
            if (a >= cnt[k]) a = cnt[k] - 1;
            do
                {
                b = (int) (cnt[k] * RandomUniform (seed));
                } while (a == b || b >= cnt[k]);
            numCA++;
            p = treeNodes + lin[k][a];
            q = treeNodes + lin[k][b];
            r = treeNodes + nextAvailableNode;
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
            lin[k][a] = nextAvailableNode;
            lin[k][b] = lin[k][cnt[k] - 1];
            cnt[k]--;
            total--;
            nextAvailableNode++;
            if (nextAvailableNode >= numNodes)
                {
                fprintf (stderr, "\n\nERROR: Too many nodes!\n");
                exit (-1);
                }
            }
        (void) j;
        for (k=0; k<K; k++)
            free (lin[k]);
        free (lin);
        free (cnt);
        free (coal);
        }

    /* make coalescence tree */
    eventNum = 0;
    period = 1;
    currentTime = 0.0;
    eventTime = 0.0;
    while (numDemes == 0 && numActiveNodes > 1)
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
				
            timeCA =  log (1 + RandomExponential (rateCA, seed) * periodGrowth[period] * 2.0 * Nbegin[period] *
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
                    FatalInfiniteCoalescentTime (timeCA);
                    }
                }
            }
        else
            {
            timeCA = RandomExponential (rateCA, seed) * N; /* time is kept in N generations */
            if (doExponential == YES)
                {
                timeCA = log (exp(growthRate*currentTime) + growthRate * timeCA) / growthRate - currentTime;
					
                /*	When growth rate is very negative, coalescent time may be infinite
                 this results in log (-x) => timCA = NaN. We have to exit
                 the program
                 */
                if (isnan(timeCA) == YES)
                    {
                    FatalInfiniteCoalescentTime (timeCA);
                    }
                }
			else if (doOthsukiInnanCoal == YES)
				{
				timeCA = RandomExponential (rateCA * 2.0 * birthRate, seed) * N;
				timeCA = log (exp(growthRate*currentTime) + growthRate * timeCA) / growthRate - currentTime;

				/*	When growth rate is very negative, coalescent time may be infinite
				 this results in log (-x) => timCA = NaN. We have to exit
				 the program
				 */
				if (isnan(timeCA) == YES)
					{
					FatalInfiniteCoalescentTime (timeCA);
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
	
  	TMRCA = currentTime / N ;  /* TMRCA is now is back in N coalescent generations */
	coalTreeMRCA = r;
	
    /* Connect the coalescent cell MRCA node with the healthy ancestral cell */
    if (noisy > 2)
        fprintf (stderr, "\n>> Adding healthy root ... ");
    healthyRoot = treeNodes + nextAvailableNode;
    healthyRoot->index = nextAvailableNode;
    healthyRoot->left = coalTreeMRCA;
    coalTreeMRCA->anc = healthyRoot;
    healthyRoot->time = currentTime + (currentTime * transformingBranchLengthRatio);
    healthyRoot->length = 0;
    healthyRoot->isHealthyRoot = YES;
	coalTreeMRCA->length = healthyRoot->time - coalTreeMRCA->time;
	coalTreeMRCA->branchLength = coalTreeMRCA->length * mutationRate;
   if (noisy > 2)
        fprintf (stderr, "DONE");

    nextAvailableNode++;
	
    /* Connect the healthy ancestral cell with the tip healthy cell*/
    if (noisy > 2)
        fprintf (stderr, "\n>> Adding healthy tip ... ");
    healthyTip = treeNodes + nextAvailableNode;
	healthyTip->index = nextAvailableNode;
  	healthyTip->left = NULL;
    healthyTip->right = NULL;
    healthyTip->anc = healthyRoot;
    healthyRoot->right = healthyTip;
    healthyTip->time = 0;
    healthyTip->length = (healthyRoot->time - healthyTip->time) * healthyTipBranchLengthRatio;
	healthyTipBranchLength = healthyTip->branchLength = healthyTip->length * mutationRate;
	healthyTip->isHealthyTip = YES;
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

void ReadUserGenome (void)
	{
	char	*header, *newGenome;
	

	/* allocate space for reading */
	header = (char*) Checked (calloc(MAX_LINE, sizeof(char)), "Could not allocate the header string\n");

	maternalUserGenome = (char*) Checked (calloc(MAX_GENOME, sizeof(char)), "Could not allocate the maternalUserGenome structure\n");

	paternalUserGenome = (char*) Checked (calloc(MAX_GENOME, sizeof(char)), "Could not allocate the paternalUserGenome structure\n");
	
	/* open user genome */
	if ((fpUserGenome = fopen(userGenomeFile, "r")) == NULL)
		{
		fprintf (stderr, "\nERROR: Can't open user genome file \"%s\"\n", userGenomeFile);
		PrintUsage(stderr);
		}

	/* read the FASTA file */
	if (fgets (header, MAX_LINE, fpUserGenome) == NULL || fgets (maternalUserGenome, MAX_GENOME, fpUserGenome) == NULL)
		{
		fprintf (stderr, "\nERROR: Can't read the user maternal genome");
		PrintUsage(stderr);
		}
	
	if (fgets (header, MAX_LINE, fpUserGenome) == NULL || fgets (paternalUserGenome, MAX_GENOME, fpUserGenome) == NULL)
		{
		fprintf (stderr, "\nERROR: Can't read the user paternal genome");
		PrintUsage(stderr);
		}

	/* remove the end-of-line characters (the last line may not have one) */
	maternalUserGenome[strcspn (maternalUserGenome, "\r\n")] = '\0';
	paternalUserGenome[strcspn (paternalUserGenome, "\r\n")] = '\0';

	/* check lengths */
	if (strlen(maternalUserGenome) != strlen(paternalUserGenome))
		{
		fprintf (stderr, "\nERROR: maternal and paternal user genomes have different lengths");
		PrintUsage(stderr);
		}

	numSites = (int) strlen(maternalUserGenome);
	
	/* realloc memory (keep room for the string terminator) */
	newGenome = (char *)realloc(maternalUserGenome, (numSites + 1) * sizeof(char));
	if (newGenome == NULL)
		{
		fprintf (stderr, "\nERROR: Could not reallocate the user genomes\n");
		exit (-1);
		}
	maternalUserGenome = newGenome;
	newGenome = (char *)realloc(paternalUserGenome, (numSites + 1) * sizeof(char));
	if (newGenome == NULL)
		{
		fprintf (stderr, "\nERROR: Could not reallocate the user genomes\n");
		exit (-1);
		}
	paternalUserGenome = newGenome;


	fprintf (stderr, "\nnumsites = %d", numSites);
	fprintf (stderr, "\nuser maternal genome = %s", maternalUserGenome);
	fprintf (stderr, "\nuser paternal genome = %s", paternalUserGenome);

	free (header);
	fclose(fpUserGenome);
	}


/* sizes of the arrays allocated by ReadUserTree, needed to free them at the end of the run
   (numNodes and numCells change while the tree is read) */
static int	userTreeAllocNodes = 0, userTreeAllocCells = 0;

void FreeUserTree (void)
	{
	int		i;

	if (treeNodes)
		{
		for (i=0; i<userTreeAllocNodes; i++)
			free (treeNodes[i].name);
		free (treeNodes);
		treeNodes = NULL;
		}
	if (cellNames)
		{
		for (i=0; i<userTreeAllocCells; i++)
			free (cellNames[i]);
		}
	}

/********************************** CountUserTreeLeaves ***********************************/
/* Counts the leaves (cells) of the user tree: number of commas + 1 in its first line.
   It is called before any array is allocated, so that -s does not need to match the tree. */
int CountUserTreeLeaves (const char *fileName)
	{
	FILE	*fp;
	char	*line;
	int		i, n = 0;

	if ((fp = fopen (fileName, "r")) == NULL)
		return 0;	/* ReadUserTree will report the error */
	line = (char*) calloc (MAX_LINE, sizeof(char));
	if (line && fgets (line, MAX_LINE, fp) != NULL)
		{
		for (i=0; line[i] != '\0' && line[i] != ';'; i++)
			if (line[i] == ',')
				n++;
		n++;
		}
	free (line);
	fclose (fp);
	return n;
	}


/********************************** ReadUserTree ***********************************/
/* Reads a rooted phyogenetic tree with branch lengths in Newick format */
/* Hangs q below p as its first free child */
static void AttachChild (TreeNode *p, TreeNode *q)
	{
	if (p->left == NULL)
		p->left = q;
	else
		p->right = q;
	q->anc = p;
	}

/* Descendant counts and mean depth of an internal node whose children are complete */
static void SummarizeDescendants (TreeNode *p)
	{
	p->numLeftDescendants  += p->left->numLeftDescendants  + p->left->numRightDescendants;
	p->numRightDescendants += p->right->numLeftDescendants + p->right->numRightDescendants;
	p->meanDepth = ((p->left->branchLength  + p->left->meanDepth)  * p->numLeftDescendants
				  + (p->right->branchLength + p->right->meanDepth) * p->numRightDescendants)
				  / (p->numLeftDescendants + p->numRightDescendants);
	}

/* Copies the characters of treeString from position *i up to (not including) a
   delimiter into buf, leaving *i on the last character read */
static void ReadToken (char *buf, int *i, const char *delims, const char *what)
	{
	int		j = 0;

	memset (buf, 0, MAX_NAME);
	while (strchr (delims, treeString[*i]) == NULL)
		{
		if (j >= MAX_NAME-1)
			{
			fprintf (stderr, "\nERROR: %s in the user tree is too long\n", what);
			exit (-1);
			}
		buf[j++] = treeString[(*i)++];
		}
	(*i)--;
	}

void ReadUserTree (void)
	{
	int			nextAvailableNode, i;
	char		*tempString;
	TreeNode	*p, *q;
	
	p = q = NULL;
	
	treeString = (char*) Checked (calloc(MAX_LINE, sizeof(char)), "Could not allocate the treeString structure\n");
	tempString = (char*) Checked (calloc(MAX_NAME, sizeof(char)), "Could not allocate the tempString structure\n");
	treeNodes  = (TreeNode *) Checked (calloc(numNodes, sizeof(TreeNode)), "Could not allocate the tree nodes\n");
	cellNames  = (char **) Checked (calloc(numCells, sizeof(char *)), "Could not allocate cellNames\n");
	userTreeAllocNodes = numNodes;
	userTreeAllocCells = numCells;
	for (i=0; i<numCells; i++)
		cellNames[i] = (char *) Checked (calloc(MAX_NAME, sizeof(char)), "Could not allocate cellNames[i]\n");
	
	/* initialize tree nodes */
	for (i=0; i<numNodes; i++)
		{
		treeNodes[i].index = -1;
		treeNodes[i].label = -1;
		treeNodes[i].isHealthyRoot = NO;
		treeNodes[i].isHealthyTip = NO;
		treeNodes[i].name = (char*) Checked (calloc (MAX_NAME, sizeof(char)), "Could not allocate the treeNodes[i].name structure\n");
		}

	/* read the tree string */
	if ((fpUserTree = fopen(userTreeFile, "r")) == NULL)
		{
		fprintf (stderr, "\nERROR: Can't open user treefile \"%s\"\n", userTreeFile);
		PrintUsage(stderr);
		}
	if (fgets (treeString, MAX_LINE, fpUserTree) == NULL)
		{
		fprintf (stderr, "\nERROR: Can't read tree string");
		PrintUsage(stderr);
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
			q = treeNodes + nextAvailableNode;
			if (nextAvailableNode == 0)
				coalTreeMRCA = q;
			else
				AttachChild (p, q);
			p = q;
			p->index = nextAvailableNode++;
			}
		else if (treeString[i] == ')' || treeString[i] == ',')
			{
			if (p == NULL || p->anc == NULL)
				{
				fprintf (stderr, "\nERROR: The user tree is not correctly formed (unbalanced parentheses or misplaced '%c')\n", treeString[i]);
				exit (-1);
				}
			p = p->anc;
			}
		else if (treeString[i] == ':')
			{
			i++;
			ReadToken (tempString, &i, ",):;(", "A branch length");
			sscanf (tempString, "%lf", &p->branchLength);
			if (p->left == NULL)
				p->numLeftDescendants = 1;
			if (p->right == NULL)
				p->numRightDescendants = 1;
			if (p->left != NULL && p->right != NULL)
				SummarizeDescendants (p);
			}
		else
			{
			/* a taxon name: every name is read as a word */
			ReadToken (tempString, &i, ",):", "A taxon name");
			taxonNamesAreChars = YES;
			q = treeNodes + nextAvailableNode;
			AttachChild (p, q);
			strcpy (q->name, tempString);
			strcpy (cellNames[numCells], tempString);
			p = q;
			p->index = nextAvailableNode++;
			numCells++;
			}
		i++;
		} while (treeString[i] != ';');
	
	/* Connect the coalescent cell MRCA node with the outgroup ancestral cell */
	p = coalTreeMRCA;
	if (p->left != NULL && p->right != NULL)
		SummarizeDescendants (p);
	
	if (noisy > 2)
		fprintf (stderr, "\n>> Adding healthy root to user tree... ");
	healthyRoot = treeNodes + nextAvailableNode;
	healthyRoot->index = nextAvailableNode;
	healthyRoot->left = coalTreeMRCA;
	healthyRoot->length = 0;
	healthyRoot->isHealthyRoot = YES;
	strcpy (healthyRoot->name, "healthyroot");
	coalTreeMRCA->anc = healthyRoot;
	coalTreeMRCA->branchLength = coalTreeMRCA->meanDepth * transformingBranchLengthRatio;
	if (noisy > 2)
		fprintf (stderr, "DONE");
	nextAvailableNode++;
	
	/* Connect the healthy ancestral cell with the tip healthy cell*/
	if (noisy > 2)
		fprintf (stderr, "\n>> Adding healthy tip to user tree... ");
	healthyTip = treeNodes + nextAvailableNode;
	healthyTip->anc = healthyRoot;
	healthyRoot->right = healthyTip;
	healthyTip->time = 0;
	healthyTip->branchLength = (coalTreeMRCA->meanDepth + coalTreeMRCA->branchLength) * healthyTipBranchLengthRatio;
	healthyTip->isHealthyTip = YES;
	strcpy (healthyTip->name, "healthy");
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

	/* Note: treeNodes is NOT shrunk here. The tree nodes point to each other (and healthyRoot
	   points into this array), so a realloc that moved the block would leave all of them
	   dangling (use-after-free). The unused tail of the array is harmless. */

	fprintf (fpLog, "\n\nUser-defined tree: ");
	PrintTree (healthyRoot, fpLog);

	free (tempString);
	}


/***************** CheckTree *******************/
/*  Checks for an error in the input tree string
	Prints message when there is an error in the
	tree and prints the correct part of the tree */
void		CheckTree (char *treeStr)
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
			if (treeStr[i+1] == ';' && (size_t) i == strlen(treeStr)-3)
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
			PrintUsage(stderr);
			}
		
		i++;
	
		} while (treeStr[i] != ';');

	if (numLeftPar != numRightPar)
		{
		fprintf (stderr, "Tree seems unbalanced (%d left and %d right parentheses)", numLeftPar, numRightPar);
		PrintUsage(stderr);
		}
	if (strrchr(treeString, ':') == NULL)
		{
		fprintf (stderr, "Tree does not have branch lengths");
		PrintUsage(stderr);
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



/*************************** SumTreeHeights **********************************/
/* Returns the total height (sum of all paths from root to tips)
	and the number of tips for for a given tree/subtree
*/
/*
 double SumTreeHeights (TreeNode *p, double total, int *numTips)
	{
	if (p == NULL)
		return 0;

	if (p->left == NULL)
		{
		(*numTips)++;
		return total;
		}
	else return SumTreeHeights(p->left, total + p->left->branchLength, numTips) +
				SumTreeHeights(p->right, total + p->right->branchLength, numTips);
	}
*/

/*************************** SumBranchlengths **********************************/
/* Returns the sum of the branch lengths for a given tree/subtree */
//because we have a static sum, we only want to call this function once
/*
 double SumBranchesOnce (TreeNode *p)
	 {
	 static double sum;
		 
	 if (p != NULL)
		 {
		 if (p->anc == NULL)
			 sum = 0;
		 else
			 sum += p->branchLength;
		 SumBranchesOnce (p->left);
		 SumBranchesOnce (p->right);
		 }
		 
	 return sum;
	 }
*/


/************************** SumBranchlengths **********************************/
/* Returns the sum of the branch lengths for a given tree/subtree */
double SumBranches (TreeNode *p)
   {
	double totalSum, leftSum, rightSum;
	totalSum = leftSum = rightSum = 0;
	
	if (p == NULL)
	   return 0;
	else
		{
		leftSum  = SumBranches(p->left);
		rightSum = SumBranches(p->right);
		totalSum =  p->branchLength + leftSum + rightSum;
		return totalSum;
		}
   }


/**************** MakeTreeNonClock **************/
/*	Introduce rate variation among lineages using a gamma multiplier  */

void MakeTreeNonClock (TreeNode *p, long int *seed)
    {
    if (p != NULL)
        {
        p->branchLength *=  RandomGamma (alphaBranches, seed) / alphaBranches;
        MakeTreeNonClock (p->left, seed);
        MakeTreeNonClock (p->right, seed);
        }
    }


/**************** TreeRateSwitch **************/
/*	Introduce rate switches on subtrees
Chooses a random branch and from there on multiplies
the branch lengths of descendants
*/

void TreeRateSwitch (TreeNode *p, long int *seed)
	{
	double		uniform, inTreeLength, cumBranchLength;
	int			i, j, numTipsSwitched, equal;
	TreeNode	**switchNodes;
	char		*concat;
	
	concat = (char *) calloc(300, sizeof(char));
	
	switchNodes = (TreeNode **) calloc(numTreeRateSwitches, sizeof(TreeNode*));
	if (!switchNodes)
		{
		fprintf (stderr, "Could not allocate switchNode (%zu)\n", numTreeRateSwitches * sizeof(TreeNode*));
		exit (-1);
		}
	
	inTreeLength = SumBranches(p) - p->branchLength;
	//fprintf(stderr,"\ninTreeBL = %6.4f;  pBL = %6.4f ; htBL = %6.4f;  total = %6.4f\n", inTreeLength, p->branchLength, healthyTipBranchLength, inTreeLength + p->branchLength + healthyTipBranchLength);
	
	// choose all switches
	for (i=0; i<numTreeRateSwitches; i++)
		{
		uniform = RandomUniform(seed) * inTreeLength;
		cumBranchLength = 0;

		for (j=0; j<numNodes; j++)
			{
			p = treeNodes + j;
			cumBranchLength += p->branchLength;
			if (cumBranchLength >= uniform)
				break;
			}
		
		// avoid switches in the same branch
		equal = NO;
		if (i>0)
			{
			for (j=i; j>=0; j--)
				if (p == switchNodes[j])
					equal = YES;
			}

		if (equal == NO)
			switchNodes[i] = p;
		else
			i--;
		}

	// do all switches
	for (i=0; i<numTreeRateSwitches; i++)
		{
		p = switchNodes[i];
		sprintf (concat, "\n  node %d  *%4.2f  ( ", p->label+1, treeRateSwitchMultiplier[i]);
		if (strlen (treeRateInfo) + strlen (concat) < MAX_LINE)
			strcat (treeRateInfo, concat);
		numTipsSwitched = 0;
		BranchRateSwitch (p, treeRateSwitchMultiplier[i], &numTipsSwitched, seed);
		if (strlen (treeRateInfo) + 2 < MAX_LINE)
			strcat (treeRateInfo, ")");
		p->rateSwitch = treeRateSwitchMultiplier[i];
		}
	free (switchNodes);
	free (concat);
	}


/**************** BranchRateSwitch **************/
/* multiplies branch lengths of the given tree/subtree */

void BranchRateSwitch (TreeNode *p, float multiplier, int *numTips, long int *seed)
	{
	if (p != NULL)
		{
		p->branchLength *= multiplier;
		if (p->left == NULL && p->right == NULL)
			{
			(*numTips)++;
			char *concat;
			concat = (char *) calloc(20, sizeof(char));
			sprintf (concat, "%d ", p->label+1);
			if (strlen (treeRateInfo) + strlen (concat) < MAX_LINE)
				strcat (treeRateInfo, concat);
			free (concat);
			}
		BranchRateSwitch (p->left, multiplier, numTips, seed);
		BranchRateSwitch (p->right, multiplier, numTips, seed);
		}
	}


