/*
	CellCoal - coalescent tree: simulation, user trees, rate variation
	See cellcoal.h for the shared declarations. This file was split from cellcoal.c.
*/
#include "cellcoal.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "definitions.h"

/************************* MakeCoalescenceTree  ************************/
/*	Builds a genealogy  under the coalescent with exponential growth without recombination */

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
		treeNodes[i].rateSwitch = 0.0;
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
                    fprintf (stderr, "\nERROR: Coalescent time (%f) is infinite ", timeCA);
                    fprintf (stderr, "\n       This might suggest that the exponential growth rate is too negative");
                    fprintf (stderr, "\n       and the coalescent time is therefore infinite.");
                    fprintf (stderr, "\n       Try a smaller value");
                    exit (-1);
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

void ReadUserGenome (FILE *fp)
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

void ReadUserTree (FILE *fp)
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


