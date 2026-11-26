/*
 * cellforward: simulation of a growing tumor forward in time, one cell division at a time
 * (the companion of CellCoal, which simulates the genealogy of a tumor backward in time).
 *
 * The tumor starts from one cell and grows as a birth-death process; mutations of the "driver" kind appear at random in the daughters of a
 * division and raise the birth rate of the cell and of all its descendants, so the fitness of a clone is the result of its drivers and
 * not an input. When the tumor reaches the size (--size) or the time (--time) asked, a random sample of cells is taken and the tree of
 * the sample is written in Newick format, one tree per replicate, ready for CellCoal (-T / --tree file=), which simulates the mutations,
 * copy number changes and reads on it.
 *
 * Time is in generations. The time between events comes from a thinning of the total rate: a cell is chosen at random, and it divides, dies
 * or does nothing with the probabilities birth / (bmax + death), death / (bmax + death) and the rest, where bmax is the largest birth
 * rate of the clones. The cost is proportional to the number of events.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <sys/stat.h>
#include <sys/types.h>

#define PROGRAM  "cellforward"
#define VERSION  "0.1.0"

/* ------------------------------ random numbers (xoshiro256**, seeded with splitmix64) ------------------------------ */

static uint64_t rngState[4];

static uint64_t SplitMix (uint64_t *x)
	{
	uint64_t z = (*x += 0x9e3779b97f4a7c15ULL);

	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
	z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
	return z ^ (z >> 31);
	}

static void SeedRandom (uint64_t seed)
	{
	int		i;

	for (i=0; i<4; i++)
		rngState[i] = SplitMix (&seed);
	}

static uint64_t Rotl (uint64_t x, int k)
	{
	return (x << k) | (x >> (64 - k));
	}

static uint64_t NextRandom (void)
	{
	uint64_t result = Rotl (rngState[1] * 5, 7) * 9, t = rngState[1] << 17;

	rngState[2] ^= rngState[0];
	rngState[3] ^= rngState[1];
	rngState[1] ^= rngState[2];
	rngState[0] ^= rngState[3];
	rngState[2] ^= t;
	rngState[3] = Rotl (rngState[3], 45);
	return result;
	}

/* uniform in [0,1) */
static double Uniform (void)
	{
	return (NextRandom () >> 11) * (1.0 / 9007199254740992.0);
	}

static double Exponential (double rate)
	{
	return -log (1.0 - Uniform ()) / rate;
	}

/* ------------------------------ parameters ------------------------------ */

static int		replicates = 1, sampleSize = 20;
static double	birthRate = -1.0, deathRate = 0.0, stopSize = 0.0, stopTime = 0.0;
static double	driverRate = 0.0, driverEffect = -1.0, mutationRate = 1.0;
static int		driverExponential = 1;
static uint64_t	seed = 0;
static const char *outDir = "cellforward_out";

static void Fail (const char *fmt, const char *a)
	{
	fprintf (stderr, "\nERROR: ");
	fprintf (stderr, fmt, a);
	fprintf (stderr, "\n");
	exit (-1);
	}

static void *Checked (void *p)
	{
	if (p == NULL)
		Fail ("out of memory (a tumor of this size does not fit: lower --size or --time)%s", "");
	return p;
	}

static void Usage (FILE *fp)
	{
	fprintf (fp,
"%s %s: simulation of a growing tumor forward in time (the tree of a sample of cells, for CellCoal)\n\n"
"usage: %s --birth B (--size N | --time T) [options]\n\n"
"growth (the tumor starts from one cell; the rates are per generation)\n"
"  --birth B           birth (division) rate of a cell, needed; the clones without drivers have this rate\n"
"  --death D           death rate of a cell (default 0)\n"
"  --size N            stop when the tumor has N cells\n"
"  --time T            stop after T generations (the tumor can have any size)\n"
"                      exactly one of --size and --time is needed; a tumor that dies out is simulated again\n\n"
"driver mutations (they raise the birth rate; the fitness of a clone is not an input but the result of its drivers)\n"
"  --driver-rate P     probability that a daughter cell of a division gets a new driver (default 0: no drivers)\n"
"  --driver-effect S   effect of a driver on the birth rate: the rate of the cell is multiplied by 1 + S (needed with --driver-rate)\n"
"  --driver-dist D     exponential (default: S is the mean of an exponential variable) or fixed (every driver has the effect S)\n\n"
"sample and output\n"
"  -s, --sample K      number of cells sampled at random from the tumor (default 20)\n"
"  -n, --replicates R  number of tumors, each with its own tree (default 1)\n"
"  --mutation-rate U   the branch lengths are the generations times U, the expected mutations per site (CellCoal reads the\n"
"                      branch lengths of a user tree as such); default 1 (the lengths are generations)\n"
"  --seed X            seed of the random numbers (default: from the clock; it is written in log.txt)\n"
"  -o, --output DIR    folder of the results (default cellforward_out)\n"
"  -h, --help          this help\n\n"
"results (in the folder): trees.nwk (one tree per line, the cells are cell0001, cell0002, ...), clones.tsv (the clone of every sampled\n"
"cell), drivers.tsv (every clone: its parent, birth time, birth rate and cells), summary.tsv (one line per replicate), log.txt\n\n"
"example, with CellCoal simulating the genotypes on the trees:\n"
"  %s --birth 0.5 --death 0.1 --size 100000 --driver-rate 1e-5 --driver-effect 0.1 -s 30 -n 100 -o fwd\n"
"  cellcoal -n100 -s30 -l10000 -Tfwd/trees.nwk\n",
		PROGRAM, VERSION, PROGRAM, PROGRAM);
	}

static double Number (const char *opt, const char *text)
	{
	char	*end;
	double	v = strtod (text, &end);

	if (end == text || *end != '\0' || !isfinite (v))
		Fail ("%s needs a number", opt);
	return v;
	}

static void ReadOptions (int argc, char **argv)
	{
	int		i;
	char	*opt, *val, name[64];
	int		haveSeed = 0;

	for (i=1; i<argc; i++)
		{
		opt = argv[i];
		val = NULL;
		if (strncmp (opt, "--", 2) == 0 && strchr (opt, '=') != NULL)
			{
			val = strchr (opt, '=') + 1;
			snprintf (name, sizeof (name), "%.*s", (int) (strchr (opt, '=') - opt), opt);
			opt = name;
			}
		if (strcmp (opt, "-h") == 0 || strcmp (opt, "--help") == 0)
			{
			Usage (stdout);
			exit (0);
			}
		if (val == NULL)
			{
			if (i + 1 >= argc)
				Fail ("%s needs a value (see --help)", opt);
			val = argv[++i];
			}
		if (strcmp (opt, "--birth") == 0)
			birthRate = Number (opt, val);
		else if (strcmp (opt, "--death") == 0)
			deathRate = Number (opt, val);
		else if (strcmp (opt, "--size") == 0)
			stopSize = Number (opt, val);
		else if (strcmp (opt, "--time") == 0)
			stopTime = Number (opt, val);
		else if (strcmp (opt, "--driver-rate") == 0)
			driverRate = Number (opt, val);
		else if (strcmp (opt, "--driver-effect") == 0)
			driverEffect = Number (opt, val);
		else if (strcmp (opt, "--driver-dist") == 0)
			{
			if (strcmp (val, "exponential") == 0)
				driverExponential = 1;
			else if (strcmp (val, "fixed") == 0)
				driverExponential = 0;
			else
				Fail ("--driver-dist is exponential or fixed (not %s)", val);
			}
		else if (strcmp (opt, "-s") == 0 || strcmp (opt, "--sample") == 0)
			sampleSize = (int) Number (opt, val);
		else if (strcmp (opt, "-n") == 0 || strcmp (opt, "--replicates") == 0)
			replicates = (int) Number (opt, val);
		else if (strcmp (opt, "--mutation-rate") == 0)
			mutationRate = Number (opt, val);
		else if (strcmp (opt, "--seed") == 0)
			{
			seed = (uint64_t) Number (opt, val);
			haveSeed = 1;
			}
		else if (strcmp (opt, "-o") == 0 || strcmp (opt, "--output") == 0)
			outDir = val;
		else
			Fail ("unknown option %s (see --help)", opt);
		}
	if (!haveSeed)
		seed = (uint64_t) time (NULL) * 2654435761ULL + (uint64_t) clock ();
	if (!(birthRate > 0.0))
		Fail ("--birth is needed (the birth rate per generation, greater than 0)%s", "");
	if (deathRate < 0.0)
		Fail ("--death cannot be negative%s", "");
	if ((stopSize > 0.0) == (stopTime > 0.0))
		Fail ("exactly one of --size (cells) and --time (generations) is needed%s", "");
	if (stopSize > 0.0 && (stopSize < 2.0 || stopSize > 1.0e9))
		Fail ("--size is between 2 and 1e9 cells%s", "");
	if (sampleSize < 2)
		Fail ("-s needs at least 2 cells%s", "");
	if (stopSize > 0.0 && sampleSize > stopSize)
		Fail ("the sample (-s) is larger than the tumor (--size)%s", "");
	if (replicates < 1)
		Fail ("-n needs at least 1 replicate%s", "");
	if (driverRate < 0.0 || driverRate > 1.0)
		Fail ("--driver-rate is a probability, between 0 and 1%s", "");
	if (driverRate > 0.0 && !(driverEffect > -1.0) )
		Fail ("--driver-effect is needed with --driver-rate (the birth rate is multiplied by 1 + S, so S > -1)%s", "");
	if (driverRate > 0.0 && driverExponential && !(driverEffect > 0.0))
		Fail ("with --driver-dist exponential the effect --driver-effect (the mean) has to be greater than 0%s", "");
	if (!(mutationRate > 0.0))
		Fail ("--mutation-rate has to be greater than 0%s", "");
	}

/* ------------------------------ the tumor ------------------------------ */

/* the lineages of all the cells ever born: a binary tree of divisions. A node is born at the time of the division of its parent. */
static int32_t	*nodeParent = NULL;
static double	*nodeTime = NULL;
static long		numNodes = 0, capNodes = 0;

/* the cells that are alive */
static int32_t	*cellNode = NULL, *cellClone = NULL;
static long		numAlive = 0, capAlive = 0;

/* the clones: a clone starts with a driver and has the drivers of its parent plus that one */
static int32_t	*cloneParent = NULL, *cloneDepth = NULL;
static double	*cloneTime = NULL, *cloneRate = NULL;		/* birth rate relative to the clones without drivers */
static long		*cloneAlive = NULL;
static long		numClones = 0, capClones = 0;

static double	bMax;

static void GrowNodes (void)
	{
	if (numNodes + 2 > capNodes)
		{
		capNodes = capNodes ? 2 * capNodes : 1024;
		if (capNodes > 2000000000L)
			Fail ("the tumor needs more than 2e9 divisions: lower --size or --time%s", "");
		nodeParent = (int32_t*) Checked (realloc (nodeParent, capNodes * sizeof (int32_t)));
		nodeTime = (double*) Checked (realloc (nodeTime, capNodes * sizeof (double)));
		}
	}

static void GrowAlive (void)
	{
	if (numAlive + 1 > capAlive)
		{
		capAlive = capAlive ? 2 * capAlive : 1024;
		cellNode = (int32_t*) Checked (realloc (cellNode, capAlive * sizeof (int32_t)));
		cellClone = (int32_t*) Checked (realloc (cellClone, capAlive * sizeof (int32_t)));
		}
	}

static int NewClone (int parent, double t)
	{
	double	s;

	if (numClones + 1 > capClones)
		{
		capClones = capClones ? 2 * capClones : 256;
		cloneParent = (int32_t*) Checked (realloc (cloneParent, capClones * sizeof (int32_t)));
		cloneDepth = (int32_t*) Checked (realloc (cloneDepth, capClones * sizeof (int32_t)));
		cloneTime = (double*) Checked (realloc (cloneTime, capClones * sizeof (double)));
		cloneRate = (double*) Checked (realloc (cloneRate, capClones * sizeof (double)));
		cloneAlive = (long*) Checked (realloc (cloneAlive, capClones * sizeof (long)));
		}
	cloneParent[numClones] = parent;
	cloneTime[numClones] = t;
	cloneAlive[numClones] = 0;
	if (parent < 0)
		{
		cloneDepth[numClones] = 0;
		cloneRate[numClones] = 1.0;
		}
	else
		{
		s = driverExponential ? driverEffect * Exponential (1.0) : driverEffect;
		cloneDepth[numClones] = cloneDepth[parent] + 1;
		cloneRate[numClones] = cloneRate[parent] * (1.0 + s);
		if (birthRate * cloneRate[numClones] > bMax)
			bMax = birthRate * cloneRate[numClones];
		}
	return (int) numClones++;
	}

/* a daughter cell of a division: it can get a driver, which makes a new clone */
static int Daughter (int clone, double t)
	{
	if (driverRate > 0.0 && Uniform () < driverRate)
		return NewClone (clone, t);
	return clone;
	}

/* one tumor, from one cell until the size or the time; returns 1 if the tumor survives to the end with enough cells, 0 if it has to be simulated again */
static int GrowTumor (double *endTime)
	{
	double	t = 0.0, rate, r, bk;
	long	i, n;
	int		k, c1, c2;

	numNodes = 0;
	numAlive = 0;
	numClones = 0;
	bMax = birthRate;
	GrowNodes ();
	GrowAlive ();
	nodeParent[0] = -1;
	nodeTime[0] = 0.0;
	numNodes = 1;
	NewClone (-1, 0.0);
	cellNode[0] = 0;
	cellClone[0] = 0;
	cloneAlive[0] = 1;
	numAlive = 1;
	for (;;)
		{
		rate = bMax + deathRate;
		t += Exponential ((double) numAlive * rate);
		if (stopTime > 0.0 && t >= stopTime)
			{
			t = stopTime;
			break;
			}
		i = (long) (Uniform () * (double) numAlive);
		k = cellClone[i];
		bk = birthRate * cloneRate[k];
		r = Uniform () * rate;
		if (r < bk)					/* division: both daughters are new lineages, and can get a driver */
			{
			GrowNodes ();
			GrowAlive ();
			c1 = (int) numNodes++;
			c2 = (int) numNodes++;
			nodeParent[c1] = nodeParent[c2] = cellNode[i];
			nodeTime[c1] = nodeTime[c2] = t;
			cloneAlive[k]--;
			n = numAlive++;
			cellNode[i] = c1;
			cellNode[n] = c2;
			cellClone[i] = Daughter (k, t);
			cellClone[n] = Daughter (k, t);
			cloneAlive[cellClone[i]]++;
			cloneAlive[cellClone[n]]++;
			if (stopSize > 0.0 && (double) numAlive >= stopSize)
				break;
			}
		else if (r < bk + deathRate)		/* death: the lineage stays in the tree of divisions, as a branch with no sampled cell */
			{
			cloneAlive[k]--;
			numAlive--;
			cellNode[i] = cellNode[numAlive];
			cellClone[i] = cellClone[numAlive];
			if (numAlive == 0)
				return 0;
			}
		}
	*endTime = t;
	return numAlive >= sampleSize;
	}

/* ------------------------------ the tree of the sample ------------------------------ */

static int32_t	*sCount, *sChildren, *sKeptAnc, *sKid0, *sKid1;
static double	*sDiv;
static unsigned char *sLeaf;
static int32_t	*leafOfNode;		/* number of the sampled cell (1 to K) of a leaf node */
static double	treeEnd;

static void WriteNewick (FILE *fp, int node, int parentKept)
	{
	double	h, hp;

	if (sKid0[node] >= 0)
		{
		fputc ('(', fp);
		WriteNewick (fp, sKid0[node], node);
		fputc (',', fp);
		WriteNewick (fp, sKid1[node], node);
		fputc (')', fp);
		h = sDiv[node];
		}
	else
		{
		fprintf (fp, "cell%04d", leafOfNode[node]);
		h = treeEnd;
		}
	if (parentKept >= 0)
		{
		hp = sDiv[parentKept];
		fprintf (fp, ":%.10g", (h - hp) * mutationRate);
		}
	}

/* samples the cells, finds the tree of the sample inside the tree of divisions and writes it; returns the time to the MRCA and the sampled cells of every clone */
static double SampleAndWriteTree (FILE *fpTrees, FILE *fpClones, int rep, int *sampledClones)
	{
	long	i, j, n = numNodes, *pick;
	int		root = -1, v, k;
	double	tmrca;

	sCount = (int32_t*) Checked (calloc (n, sizeof (int32_t)));
	sChildren = (int32_t*) Checked (calloc (n, sizeof (int32_t)));
	sKeptAnc = (int32_t*) Checked (malloc (n * sizeof (int32_t)));
	sKid0 = (int32_t*) Checked (malloc (n * sizeof (int32_t)));
	sKid1 = (int32_t*) Checked (malloc (n * sizeof (int32_t)));
	leafOfNode = (int32_t*) Checked (calloc (n, sizeof (int32_t)));
	sDiv = (double*) Checked (calloc (n, sizeof (double)));
	sLeaf = (unsigned char*) Checked (calloc (n, 1));
	pick = (long*) Checked (malloc (numAlive * sizeof (long)));
	for (i=0; i<n; i++)
		sKid0[i] = sKid1[i] = -1;
	for (i=0; i<numAlive; i++)
		pick[i] = i;
	for (i=0; i<sampleSize; i++)		/* a random sample without replacement (partial shuffle) */
		{
		j = i + (long) (Uniform () * (double) (numAlive - i));
		v = (int) pick[j];
		pick[j] = pick[i];
		pick[i] = v;
		sLeaf[cellNode[v]] = 1;
		sCount[cellNode[v]] = 1;
		leafOfNode[cellNode[v]] = (int32_t) (i + 1);
		if (fpClones != NULL)
			{
			k = cellClone[v];
			fprintf (fpClones, "%d\tcell%04ld\t%d\t%d\t%.6g\n", rep, i + 1, k, cloneDepth[k], cloneRate[k]);
			sampledClones[k]++;
			}
		}
	for (i=n-1; i>0; i--)				/* the parents have a smaller number than their children */
		if (sCount[i] > 0)
			{
			sCount[nodeParent[i]] += sCount[i];
			sChildren[nodeParent[i]]++;
			sDiv[nodeParent[i]] = nodeTime[i];	/* the time of the division of the parent: when its children are born */
			}
	for (i=0; i<n; i++)					/* nearest ancestor that is a leaf or a branching node */
		{
		if (i == 0)
			sKeptAnc[i] = -1;
		else
			sKeptAnc[i] = (sLeaf[nodeParent[i]] || sChildren[nodeParent[i]] >= 2) ? nodeParent[i] : sKeptAnc[nodeParent[i]];
		}
	for (i=0; i<n; i++)
		if (sCount[i] > 0 && (sLeaf[i] || sChildren[i] >= 2))
			{
			if (sKeptAnc[i] < 0)
				root = (int) i;
			else if (sKid0[sKeptAnc[i]] < 0)
				sKid0[sKeptAnc[i]] = (int32_t) i;
			else
				sKid1[sKeptAnc[i]] = (int32_t) i;
			}
	if (root < 0)
		Fail ("internal error: no root in the tree of the sample%s", "");
	WriteNewick (fpTrees, root, -1);
	fputs (";\n", fpTrees);
	tmrca = sDiv[root];
	free (sCount);
	free (sChildren);
	free (sKeptAnc);
	free (sKid0);
	free (sKid1);
	free (leafOfNode);
	free (sDiv);
	free (sLeaf);
	free (pick);
	return tmrca;
	}

/* ------------------------------ main ------------------------------ */

int main (int argc, char **argv)
	{
	FILE	*fpTrees, *fpClones, *fpDrivers, *fpSummary, *fpLog;
	char	path[4096];
	int		rep, attempts, i, usedClones, *sampledClones, drivers;
	double	t, tmrca, mean;
	long	totalAttempts = 0;

	if (argc < 2)
		{
		Usage (stderr);
		return 1;
		}
	ReadOptions (argc, argv);
	SeedRandom (seed);
	mkdir (outDir, 0777);
	snprintf (path, sizeof (path), "%s/trees.nwk", outDir);
	fpTrees = fopen (path, "w");
	snprintf (path, sizeof (path), "%s/clones.tsv", outDir);
	fpClones = fopen (path, "w");
	snprintf (path, sizeof (path), "%s/drivers.tsv", outDir);
	fpDrivers = fopen (path, "w");
	snprintf (path, sizeof (path), "%s/summary.tsv", outDir);
	fpSummary = fopen (path, "w");
	snprintf (path, sizeof (path), "%s/log.txt", outDir);
	fpLog = fopen (path, "w");
	if (!fpTrees || !fpClones || !fpDrivers || !fpSummary || !fpLog)
		Fail ("cannot write in the folder %s", outDir);
	fprintf (fpLog, "%s %s\ncommand:", PROGRAM, VERSION);
	for (i=1; i<argc; i++)
		fprintf (fpLog, " %s", argv[i]);
	fprintf (fpLog, "\nseed: %llu\n", (unsigned long long) seed);
	fclose (fpLog);
	fprintf (fpClones, "replicate\tcell\tclone\tdrivers\tbirth_rate_relative\n");
	fprintf (fpDrivers, "replicate\tclone\tparent_clone\tdrivers\ttime\tbirth_rate_relative\tcells_at_end\tcells_sampled\n");
	fprintf (fpSummary, "replicate\tattempts\ttime\tcells\tclones\tmax_drivers\ttmrca\n");
	for (rep=1; rep<=replicates; rep++)
		{
		attempts = 0;
		while (!GrowTumor (&t))
			if (++attempts >= 100000)
				Fail ("the tumor dies out, or has fewer cells than the sample, in 100000 attempts: raise --birth, or lower --death or -s%s", "");
		totalAttempts += attempts + 1;
		treeEnd = t;
		usedClones = (int) numClones;
		sampledClones = (int*) Checked (calloc (usedClones, sizeof (int)));
		tmrca = SampleAndWriteTree (fpTrees, fpClones, rep, sampledClones);
		drivers = 0;
		for (i=0; i<usedClones; i++)
			{
			if (cloneDepth[i] > drivers)
				drivers = cloneDepth[i];
			if (cloneAlive[i] > 0 || sampledClones[i] > 0)
				fprintf (fpDrivers, "%d\t%d\t%d\t%d\t%.6g\t%.6g\t%ld\t%d\n", rep, i, cloneParent[i], cloneDepth[i], cloneTime[i], cloneRate[i], cloneAlive[i], sampledClones[i]);
			}
		mean = t - tmrca;
		fprintf (fpSummary, "%d\t%d\t%.6g\t%ld\t%ld\t%d\t%.6g\n", rep, attempts + 1, t, numAlive, numClones, drivers, mean);
		free (sampledClones);
		}
	fclose (fpTrees);
	fclose (fpClones);
	fclose (fpDrivers);
	fclose (fpSummary);
	fprintf (stderr, "%s: %d replicates (%ld tumors simulated), results in %s/\n", PROGRAM, replicates, totalAttempts, outDir);
	return 0;
	}
