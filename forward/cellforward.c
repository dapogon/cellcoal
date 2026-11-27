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
#include <sys/wait.h>
#include <unistd.h>

#define MAX_CELLS 5.0e7
#define MAX_PRESET 64
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

static double Normal (void)
	{
	return sqrt (-2.0 * log (1.0 - Uniform ())) * cos (6.283185307179586 * Uniform ());
	}

/* gamma with shape k and scale 1 (Marsaglia and Tsang) */
static double Gamma (double k)
	{
	double	d, c, x, v, u;

	if (k < 1.0)
		return Gamma (k + 1.0) * pow (1.0 - Uniform (), 1.0 / k);
	d = k - 1.0 / 3.0;
	c = 1.0 / sqrt (9.0 * d);
	for (;;)
		{
		do
			{
			x = Normal ();
			v = 1.0 + c * x;
			}
		while (v <= 0.0);
		v = v * v * v;
		u = Uniform ();
		if (u < 1.0 - 0.0331 * x * x * x * x || log (1.0 - u) < 0.5 * x * x + d * (1.0 - v + log (v)))
			return d * v;
		}
	}

/* ------------------------------ parameters ------------------------------ */

static int		replicates = 1, sampleSize = 20;
static double	birthRate = -1.0, deathRate = 0.0, stopSize = 0.0, stopTime = 0.0;
static double	driverRate = 0.0, driverEffect = -1.0, mutationRate = 1.0, divisionMutationRate = 0.0;
static int		driverExponential = 1, driverOnGrowth = 1;
static int		driverRateGiven = 0, driverOnDeath = 0, jobs = 1;
static int		treeMethod = 0;				/* 0 auto, 1 forward, 2 reconstructed */
static double	cycleShape = 1.0;						/* shape k of the gamma waiting time between divisions (1: exponential) */
static double	driverTarget = 3.5;						/* expected number of clones with a driver per tumor when --driver-rate is not given */
static long		maxClones = 0, driverClones = 0;		/* limit of the clones with a random driver per tumor (0 = no limit) */
static double	mutatorRate = 0.0, mutatorFactor = -1.0, presetMutation[MAX_PRESET];	/* hypermutators */
static int		mutatorsOn = 0;
static int		presetN = 0, presetParent[MAX_PRESET];			/* clones that the user starts at a given time (clone 1, 2, ...; the first clone, without drivers, is 0) */
static double	presetTime[MAX_PRESET], presetFitness[MAX_PRESET];
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
"  --driver-rate P     probability that a daughter cell of a division gets a new driver (default: none without --driver-effect; with it, a rate\n"
"                      that gives about 3 or 4 clones with a driver per tumor, from --size or --time)\n"
"  --driver-effect S   effect of a driver on the birth rate: the rate of the cell is multiplied by 1 + S (needed with --driver-rate)\n"
"  --max-clones N      at most N clones with a random driver per tumor (default 0: no limit; the clones of --clone-time do not count).\n"
"                      The N clones are the first that appear; to keep the natural probabilities lower --driver-rate instead\n"
"  --driver-dist D     exponential (default: S is the mean of an exponential variable) or fixed (every driver has the effect S)\n"
"  --driver-acts-on W  growth (default): the net growth rate (birth - death) is multiplied by 1 + S, as the fitness of --clone-fitness;\n"
"                      birth: the birth rate is multiplied by 1 + S and the death rate does not change (the gain in growth is then larger than S);\n"
"                      death: the DEATH rate falls so that the net growth rate is multiplied by 1 + S (the birth rate does not change; the death rate\n"
"                      cannot go below 0, so the gain is at most birth / (birth - death)); needs --death greater than 0\n\n"
"clones started by the user (a cell of the parent clone gets a growth advantage at a time; the size and the frequency of the clone are a result)\n"
"  --clone-time T1:T2:...   generations after the start of the tumor at which each clone starts (increasing; at most %d clones)\n"
"  --clone-fitness F1:F2:.. growth rate of each clone relative to the first clone, as in CellCoal --clones: the birth rate of the clone\n"
"                      is death + F (birth - death); 1.5 grows 50%% faster\n"
"  --clone-parent P1:P2:... clone (0 is the first clone, 1 the first started by the user, ...) from which each clone starts (default 0)\n"
"                      a clone that has no cell when it has to start, or a tumor that reaches --size before the time of a clone, is simulated again\n\n"
"  --clone-mutation-factor M1:M2:...  mutation rate of each clone relative to the normal cells (a hypermutator clone has M > 1); default 1.\n"
"                      The mutation rate does not need to follow the growth: use --clone-fitness 1 for a clone that grows as the first one\n\n"
"hypermutators that appear at random (independent of the drivers of growth; the mutation rate is inherited by the descendants)\n"
"  --mutator-rate P    probability that a daughter cell of a division gets a mutator mutation (default 0)\n"
"  --mutator-factor M  the mutation rate of the cell (per generation and per division) is multiplied by M (needed with --mutator-rate)\n\n"
"sample and output\n"
"  -s, --sample K      number of cells sampled at random from the tumor (default 20)\n"
"  -n, --replicates R  number of tumors, each with its own tree (default 1)\n"
"  -j, --jobs J       replicates run at the same time in J processes (default 1; 0: all the cores). The results do not depend on J: every replicate has\n"
"                      its own random numbers, from --seed and the number of the replicate\n"
"  --tree-method M    auto (default), forward or reconstructed. A neutral tumor (no drivers, clones or hypermutators, --cycle-shape 1) is a birth-death\n"
"                      process: with auto its tree is drawn directly (coalescent point process), in a time that does not depend on the size of the tumor;\n"
"                      forward simulates all the cells. The two give the same distribution (see the manual)\n"
"  --mutation-rate U   mutations per site and generation: the branch length is the generations times U plus the divisions times V\n"
"                      (CellCoal reads the branch lengths of a user tree as expected mutations per site); default 1, which gives generations\n"
"  --cycle-shape K     waiting time between divisions: gamma with shape K and mean 1 / birth rate (default 1: exponential, the cells have no memory of their\n"
"                      last division; a large K makes the cells divide almost in synchrony, as in a cell cycle). Deaths stay exponential. A clone of\n"
"                      --clone-time starts its cycle again\n"
"  --division-mutation-rate V  mutations per site and cell division (those linked to the replication of the DNA; default 0). With V the number of\n"
"                      mutations of a branch depends on the divisions of the lineage and not only on the time, so a clone that divides more\n"
"                      accumulates more mutations. Use --mutation-rate 0 for divisions only\n"
"  --seed X            seed of the random numbers (default: from the clock; it is written in log.txt)\n"
"  -o, --output DIR    folder of the results (default cellforward_out)\n"
"  -h, --help          this help\n\n"
"results (in the folder): trees.nwk (one tree per line, the cells are cell0001, cell0002, ...), clones.tsv (the clone of every sampled\n"
"cell), drivers.tsv (every clone: its parent, birth time, birth rate and cells), summary.tsv (one line per replicate), log.txt\n\n"
"example, with CellCoal simulating the genotypes on the trees:\n"
"  %s --birth 0.5 --death 0.1 --size 100000 --driver-rate 1e-5 --driver-effect 0.1 -s 30 -n 100 -o fwd\n"
"  cellcoal -n100 -s30 -l10000 -Tfwd/trees.nwk\n",
		PROGRAM, VERSION, PROGRAM, MAX_PRESET, PROGRAM);
	}

static double Number (const char *opt, const char *text)
	{
	char	*end;
	double	v = strtod (text, &end);

	if (end == text || *end != '\0' || !isfinite (v))
		Fail ("%s needs a number", opt);
	return v;
	}

/* a list of numbers separated by : */
static int NumberList (const char *opt, const char *text, double *v, int max)
	{
	int		n = 0;
	char	*copy = strdup (text), *tok, *save = NULL;

	for (tok = strtok_r (copy, ":", &save); tok != NULL; tok = strtok_r (NULL, ":", &save))
		{
		if (n >= max)
			Fail ("%s has too many values", opt);
		v[n++] = Number (opt, tok);
		}
	free (copy);
	if (n == 0)
		Fail ("%s needs numbers separated by :", opt);
	return n;
	}

static void ReadOptions (int argc, char **argv)
	{
	int		i;
	char	*opt, *val, name[64];
	int		haveSeed = 0, nTime = 0, nFit = 0, nPar = 0, nMut = 0, k;
	double	parents[MAX_PRESET];

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
			{
			driverRate = Number (opt, val);
			driverRateGiven = 1;
			}
		else if (strcmp (opt, "--cycle-shape") == 0)
			cycleShape = Number (opt, val);
		else if (strcmp (opt, "--max-clones") == 0)
			maxClones = (long) Number (opt, val);
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
		else if (strcmp (opt, "--clone-time") == 0)
			nTime = NumberList (opt, val, presetTime, MAX_PRESET);
		else if (strcmp (opt, "--clone-fitness") == 0)
			nFit = NumberList (opt, val, presetFitness, MAX_PRESET);
		else if (strcmp (opt, "--clone-mutation-factor") == 0)
			nMut = NumberList (opt, val, presetMutation, MAX_PRESET);
		else if (strcmp (opt, "--mutator-rate") == 0)
			mutatorRate = Number (opt, val);
		else if (strcmp (opt, "--mutator-factor") == 0)
			mutatorFactor = Number (opt, val);
		else if (strcmp (opt, "--clone-parent") == 0)
			nPar = NumberList (opt, val, parents, MAX_PRESET);
		else if (strcmp (opt, "-s") == 0 || strcmp (opt, "--sample") == 0)
			sampleSize = (int) Number (opt, val);
		else if (strcmp (opt, "--tree-method") == 0)
			{
			if (strcmp (val, "auto") == 0)
				treeMethod = 0;
			else if (strcmp (val, "forward") == 0)
				treeMethod = 1;
			else if (strcmp (val, "reconstructed") == 0)
				treeMethod = 2;
			else
				Fail ("--tree-method is auto, forward or reconstructed (not %s)", val);
			}
		else if (strcmp (opt, "-j") == 0 || strcmp (opt, "--jobs") == 0)
			jobs = (int) Number (opt, val);
		else if (strcmp (opt, "-n") == 0 || strcmp (opt, "--replicates") == 0)
			replicates = (int) Number (opt, val);
		else if (strcmp (opt, "--mutation-rate") == 0)
			mutationRate = Number (opt, val);
		else if (strcmp (opt, "--division-mutation-rate") == 0)
			divisionMutationRate = Number (opt, val);
		else if (strcmp (opt, "--driver-acts-on") == 0)
			{
			if (strcmp (val, "growth") == 0)
				driverOnGrowth = 1;
			else if (strcmp (val, "birth") == 0)
				driverOnGrowth = 0;
			else if (strcmp (val, "death") == 0)
				{
				driverOnGrowth = 0;
				driverOnDeath = 1;
				}
			else
				Fail ("--driver-acts-on is growth, birth or death (not %s)", val);
			}
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
	if (stopSize > 0.0 && (stopSize < 2.0 || stopSize > 1e15))
		Fail ("--size is between 2 and 1e15 cells%s", "");
	if (sampleSize < 2)
		Fail ("-s needs at least 2 cells%s", "");
	if (stopSize > 0.0 && sampleSize > stopSize)
		Fail ("the sample (-s) is larger than the tumor (--size)%s", "");
	if (replicates < 1)
		Fail ("-n needs at least 1 replicate%s", "");
	if (!driverRateGiven && driverEffect > -1.0)		/* the rate of the drivers is chosen so that a tumor has a few clones (about 3 or 4) */
		{
		double	g = birthRate - deathRate, cells;

		if (!(g > 0.0))
			Fail ("--driver-rate is needed when the tumor does not grow (--birth not greater than --death)%s", "");
		cells = (stopSize > 0.0) ? stopSize : exp (g * stopTime) * birthRate / g;		/* with --time: the size of the tumors that survive */
		driverRate = driverTarget * g / (2.0 * birthRate * cells);		/* the tumor makes about cells * birth / growth divisions, with two daughters each */
		if (driverRate > 1.0)
			driverRate = 1.0;
		fprintf (stderr, "--driver-rate not given: %.3g, for about %.1f clones with a driver per tumor\n", driverRate, driverTarget);
		}
	if (!(cycleShape > 0.0) || cycleShape > 1e6)
		Fail ("--cycle-shape is the shape of a gamma distribution, greater than 0 (and not above 1e6)%s", "");
	if (driverRate > 0.0 && driverOnDeath && !(deathRate > 0.0))
		Fail ("--driver-acts-on death needs a death rate (--death greater than 0)%s", "");
	if (driverEffect > -1.0 && driverOnDeath && !(birthRate > deathRate))
		Fail ("--driver-acts-on death needs a tumor that grows (--birth greater than --death)%s", "");
	if (maxClones < 0)
		Fail ("--max-clones is a number of clones, 0 or more%s", "");
	if (driverRate < 0.0 || driverRate > 1.0)
		Fail ("--driver-rate is a probability, between 0 and 1%s", "");
	if (driverRate > 0.0 && !(driverEffect > -1.0) )
		Fail ("--driver-effect is needed with --driver-rate (the birth rate is multiplied by 1 + S, so S > -1)%s", "");
	if (driverRate > 0.0 && driverExponential && !(driverEffect > 0.0))
		Fail ("with --driver-dist exponential the effect --driver-effect (the mean) has to be greater than 0%s", "");
	if (mutationRate < 0.0 || divisionMutationRate < 0.0 || (mutationRate == 0.0 && divisionMutationRate == 0.0))
		Fail ("--mutation-rate and --division-mutation-rate cannot be negative, and one of them has to be greater than 0%s", "");
	if (driverRate > 0.0 && driverOnGrowth && !(birthRate > deathRate))
		Fail ("--driver-acts-on growth needs a tumor that grows (--birth greater than --death); use --driver-acts-on birth%s", "");
	presetN = nTime;
	if (mutatorRate > 0.0)
		mutatorsOn = 1;
	if (nTime != nFit || (nPar > 0 && nPar != nTime))
		Fail ("--clone-time, --clone-fitness and --clone-parent need the same number of values (--clone-parent is optional)%s", "");
	if (nMut > 0 && nMut != nTime)
		Fail ("--clone-mutation-factor needs one value for every clone of --clone-time%s", "");
	if (mutatorRate < 0.0 || mutatorRate > 1.0)
		Fail ("--mutator-rate is a probability, between 0 and 1%s", "");
	if (mutatorRate > 0.0 && !(mutatorFactor > 0.0))
		Fail ("--mutator-factor is needed with --mutator-rate (greater than 0)%s", "");
	for (k=0; k<presetN; k++)
		{
		if (nMut == 0)
			presetMutation[k] = 1.0;
		else if (!(presetMutation[k] > 0.0))
			Fail ("--clone-mutation-factor: the factors have to be greater than 0%s", "");
		if (presetMutation[k] != 1.0)
			mutatorsOn = 1;
		presetParent[k] = (nPar > 0) ? (int) parents[k] : 0;
		if (!(presetTime[k] > 0.0) || (k > 0 && !(presetTime[k] > presetTime[k-1])))
			Fail ("--clone-time: the times are greater than 0 and increasing%s", "");
		if (stopTime > 0.0 && presetTime[k] >= stopTime)
			Fail ("--clone-time: a clone starts after the end of the tumor (--time)%s", "");
		if (!(presetFitness[k] > 0.0))
			Fail ("--clone-fitness: the fitness (growth rate relative to the first clone) has to be greater than 0%s", "");
		if (presetParent[k] < 0 || presetParent[k] > k || (nPar > 0 && parents[k] != (double) presetParent[k]))
			Fail ("--clone-parent: the parent of a clone is a clone that starts before it (0 is the first clone)%s", "");
		if (deathRate + presetFitness[k] * (birthRate - deathRate) <= 0.0)
			Fail ("--clone-fitness: the birth rate of a clone would not be positive%s", "");
		}
	if (presetN > 0 && !(birthRate > deathRate))
		Fail ("clones of different fitness need a tumor that grows (--birth greater than --death)%s", "");
	}

/* ------------------------------ the tumor ------------------------------ */

/* the lineages of all the cells ever born: a binary tree of divisions. A node is born at the time of the division of its parent. */
static int32_t	*nodeParent = NULL;
static double	*nodeTime = NULL;
static int32_t	*nodeDivs = NULL;			/* divisions of the lineage from the start of the tumor */
static double	*nodeMutStart = NULL, *nodeMutEnd = NULL;	/* with hypermutators: mutations (per site) of the lineage when the node starts and ends */
static long		numNodes = 0, capNodes = 0;
static double	*nodeNext = NULL;			/* gamma waiting times: time of the next event of the cell of a node */
static int32_t	*nodeCell = NULL;			/* gamma waiting times: alive cell of a node (-1: the node has divided or died) */
static unsigned char *nodeDeath = NULL;		/* gamma waiting times: the next event is a death */

/* the cells that are alive */
static int32_t	*cellNode = NULL, *cellClone = NULL;
static long		numAlive = 0, capAlive = 0;
static double	*cellBase = NULL, *cellLast = NULL, *cellMult = NULL;	/* with hypermutators: mutations of the lineage at the last update, time of it, mutation factor */

/* mutations accumulated by the lineage of the cell i up to the time t, at its mutation factor */
static void SyncMutations (long i, double t)
	{
	cellBase[i] += (t - cellLast[i]) * mutationRate * cellMult[i];
	cellLast[i] = t;
	}

/* the clones: a clone starts with a driver and has the drivers of its parent plus that one */
static int32_t	*cloneParent = NULL, *cloneDepth = NULL;
static double	*cloneTime = NULL, *cloneRate = NULL;		/* birth rate relative to the clones without drivers */
static long		*cloneAlive = NULL;
static long		numClones = 0, capClones = 0;

static double	bMax, dMax;
static double	*cloneDeath = NULL;						/* death rate of every clone */

static void GrowNodes (void)
	{
	if (numNodes + 2 > capNodes)
		{
		capNodes = capNodes ? 2 * capNodes : 1024;
		if (capNodes > 2000000000L)
			Fail ("the tumor needs more than 2e9 divisions: lower --size or --time%s", "");
		nodeParent = (int32_t*) Checked (realloc (nodeParent, capNodes * sizeof (int32_t)));
		nodeTime = (double*) Checked (realloc (nodeTime, capNodes * sizeof (double)));
		nodeDivs = (int32_t*) Checked (realloc (nodeDivs, capNodes * sizeof (int32_t)));
		if (cycleShape != 1.0)
			{
			nodeNext = (double*) Checked (realloc (nodeNext, capNodes * sizeof (double)));
			nodeCell = (int32_t*) Checked (realloc (nodeCell, capNodes * sizeof (int32_t)));
			nodeDeath = (unsigned char*) Checked (realloc (nodeDeath, capNodes * sizeof (unsigned char)));
			}
		if (mutatorsOn)
			{
			nodeMutStart = (double*) Checked (realloc (nodeMutStart, capNodes * sizeof (double)));
			nodeMutEnd = (double*) Checked (realloc (nodeMutEnd, capNodes * sizeof (double)));
			}
		}
	}

static void GrowAlive (void)
	{
	if (numAlive + 1 > capAlive)
		{
		capAlive = capAlive ? 2 * capAlive : 1024;
		cellNode = (int32_t*) Checked (realloc (cellNode, capAlive * sizeof (int32_t)));
		cellClone = (int32_t*) Checked (realloc (cellClone, capAlive * sizeof (int32_t)));
		if (mutatorsOn)
			{
			cellBase = (double*) Checked (realloc (cellBase, capAlive * sizeof (double)));
			cellLast = (double*) Checked (realloc (cellLast, capAlive * sizeof (double)));
			cellMult = (double*) Checked (realloc (cellMult, capAlive * sizeof (double)));
			}
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
		cloneDeath = (double*) Checked (realloc (cloneDeath, capClones * sizeof (double)));
		cloneAlive = (long*) Checked (realloc (cloneAlive, capClones * sizeof (long)));
		}
	cloneParent[numClones] = parent;
	cloneTime[numClones] = t;
	cloneAlive[numClones] = 0;
	if (parent < 0)
		{
		cloneDepth[numClones] = 0;
		cloneRate[numClones] = 1.0;
		cloneDeath[numClones] = deathRate;
		}
	else
		{
		s = driverExponential ? driverEffect * Exponential (1.0) : driverEffect;
		cloneDepth[numClones] = cloneDepth[parent] + 1;
		cloneDeath[numClones] = cloneDeath[parent];
		if (driverOnDeath)			/* the death rate falls so that the net growth rate is multiplied by 1 + s */
			{
			cloneRate[numClones] = cloneRate[parent];
			cloneDeath[numClones] = birthRate * cloneRate[parent] - (1.0 + s) * (birthRate * cloneRate[parent] - cloneDeath[parent]);
			if (cloneDeath[numClones] < 0.0)
				cloneDeath[numClones] = 0.0;
			}
		else if (driverOnGrowth)			/* the net growth rate is multiplied by 1 + s */
			cloneRate[numClones] = (deathRate + (1.0 + s) * (birthRate * cloneRate[parent] - deathRate)) / birthRate;
		else						/* the birth rate is multiplied by 1 + s */
			cloneRate[numClones] = cloneRate[parent] * (1.0 + s);
		if (birthRate * cloneRate[numClones] > bMax)
			bMax = birthRate * cloneRate[numClones];
		if (cloneDeath[numClones] > dMax)
			dMax = cloneDeath[numClones];
		}
	return (int) numClones++;
	}

/* a daughter cell of a division: it can get a driver, which makes a new clone */
static int Daughter (int clone, double t)
	{
	if (driverRate > 0.0 && (maxClones == 0 || driverClones < maxClones) && Uniform () < driverRate)
		{
		driverClones++;
		return NewClone (clone, t);
		}
	return clone;
	}

/* a division of the alive cell i at the time t: both daughters are new lineages, and can get a driver */
static void Division (long i, double t)
	{
	long	n;
	int		k = cellClone[i], c1, c2;

	GrowNodes ();
	GrowAlive ();
	c1 = (int) numNodes++;
	c2 = (int) numNodes++;
	nodeParent[c1] = nodeParent[c2] = cellNode[i];
	nodeTime[c1] = nodeTime[c2] = t;
	nodeDivs[c1] = nodeDivs[c2] = nodeDivs[cellNode[i]] + 1;
	if (cycleShape != 1.0)
		nodeCell[cellNode[i]] = -1;
	if (mutatorsOn)			/* mutations per generation until now, and per division: both at the mutation factor of the dividing cell */
		{
		SyncMutations (i, t);
		nodeMutEnd[cellNode[i]] = cellBase[i];
		cellBase[i] += divisionMutationRate * cellMult[i];
		nodeMutStart[c1] = nodeMutStart[c2] = cellBase[i];
		cellBase[numAlive] = cellBase[i];
		cellLast[numAlive] = t;
		cellMult[numAlive] = cellMult[i];
		if (mutatorRate > 0.0 && Uniform () < mutatorRate)
			cellMult[i] *= mutatorFactor;
		if (mutatorRate > 0.0 && Uniform () < mutatorRate)
			cellMult[numAlive] *= mutatorFactor;
		}
	cloneAlive[k]--;
	n = numAlive++;
	cellNode[i] = c1;
	cellNode[n] = c2;
	cellClone[i] = Daughter (k, t);
	cellClone[n] = Daughter (k, t);
	cloneAlive[cellClone[i]]++;
	cloneAlive[cellClone[n]]++;
	}

/* a death of the alive cell i: the lineage stays in the tree of divisions, as a branch with no sampled cell */
static void Death (long i)
	{
	cloneAlive[cellClone[i]]--;
	if (cycleShape != 1.0)
		nodeCell[cellNode[i]] = -1;
	numAlive--;
	cellNode[i] = cellNode[numAlive];
	cellClone[i] = cellClone[numAlive];
	if (cycleShape != 1.0 && i < numAlive)
		nodeCell[cellNode[i]] = (int32_t) i;
	if (mutatorsOn)
		{
		cellBase[i] = cellBase[numAlive];
		cellLast[i] = cellLast[numAlive];
		cellMult[i] = cellMult[numAlive];
		}
	}

/* a clone of the user starts: a cell of its parent clone gets the advantage; returns the cell, or -1 if the parent clone has no cell */
static long StartPresetClone (int nextPreset, double t)
	{
	long	i, tries;
	int		j = nextPreset + 1, p = cloneParent[j];

	if (cloneAlive[p] <= 0)
		return -1;
	for (tries=0, i=0; tries<100000; tries++)
		{
		i = (long) (Uniform () * (double) numAlive);
		if (cellClone[i] == p)
			break;
		}
	if (tries == 100000)			/* the parent is rare: look for one of its cells */
		for (i=(long) (Uniform () * (double) numAlive), tries=0; tries<numAlive && cellClone[i] != p; tries++)
			i = (i + 1) % numAlive;
	if (mutatorsOn)
		{
		SyncMutations (i, t);
		cellMult[i] = presetMutation[j-1];
		}
	cloneAlive[p]--;
	cellClone[i] = j;
	cloneAlive[j]++;
	return i;
	}

/* gamma waiting times: the events of the cells in a heap ordered by time (an event is valid if its node still has the same time) */
typedef struct { double t; int32_t node; } Event;
static Event	*heap = NULL;
static long		heapN = 0, heapCap = 0;

static void HeapPush (double t, int32_t node)
	{
	long	c, p;

	if (heapN + 1 > heapCap)
		{
		heapCap = heapCap ? 2 * heapCap : 1024;
		heap = (Event*) Checked (realloc (heap, heapCap * sizeof (Event)));
		}
	for (c=heapN++; c>0; c=p)
		{
		p = (c - 1) / 2;
		if (heap[p].t <= t)
			break;
		heap[c] = heap[p];
		}
	heap[c].t = t;
	heap[c].node = node;
	}

static void HeapPop (void)
	{
	Event	last;
	long	c, p = 0;

	last = heap[--heapN];
	if (heapN == 0)
		return;
	/* bottom-up: the hole goes down to a leaf following the smaller child (one comparison per level), and the last event climbs from there */
	for (;;)
		{
		c = 2 * p + 1;
		if (c >= heapN)
			break;
		if (c + 1 < heapN && heap[c+1].t < heap[c].t)
			c++;
		heap[p] = heap[c];
		p = c;
		}
	for (; p>0; p=c)
		{
		c = (p - 1) / 2;
		if (heap[c].t <= last.t)
			break;
		heap[p] = heap[c];
		}
	heap[p] = last;
	}

/* the next event of the cell i, which is new or has changed clone at the time t: its division (gamma) or its death (exponential), whichever comes first */
static void Schedule (long i, double t)
	{
	int32_t	node = cellNode[i];
	int		k = cellClone[i];
	double	bk = birthRate * cloneRate[k], td, tx = 1e300;

	td = t + Gamma (cycleShape) / (cycleShape * bk);
	if (cloneDeath[k] > 0.0)
		tx = t + Exponential (cloneDeath[k]);
	nodeCell[node] = (int32_t) i;
	nodeDeath[node] = (tx < td);
	nodeNext[node] = (tx < td) ? tx : td;
	HeapPush (nodeNext[node], node);
	}

/* one tumor, from one cell until the size or the time; returns 1 if the tumor survives to the end with enough cells, 0 if it has to be simulated again */
static int GrowTumor (double *endTime)
	{
	double	t = 0.0, rate, r, bk;
	long	i, n;
	int		k, nextPreset = 0, j;

	numNodes = 0;
	numAlive = 0;
	numClones = 0;
	driverClones = 0;
	heapN = 0;
	bMax = birthRate;
	dMax = deathRate;
	GrowNodes ();
	GrowAlive ();
	nodeParent[0] = -1;
	nodeTime[0] = 0.0;
	nodeDivs[0] = 0;
	if (mutatorsOn)
		{
		nodeMutStart[0] = 0.0;
		cellBase[0] = 0.0;
		cellLast[0] = 0.0;
		cellMult[0] = 1.0;
		}
	numNodes = 1;
	NewClone (-1, 0.0);
	for (j=0; j<presetN; j++)			/* the clones of the user, with the numbers 1, 2, ... */
		{
		NewClone (-1, presetTime[j]);
		cloneParent[j+1] = presetParent[j];
		cloneDepth[j+1] = cloneDepth[presetParent[j]] + 1;
		cloneRate[j+1] = (deathRate + presetFitness[j] * (birthRate - deathRate)) / birthRate;
		if (birthRate * cloneRate[j+1] > bMax)
			bMax = birthRate * cloneRate[j+1];
		}
	cellNode[0] = 0;
	cellClone[0] = 0;
	cloneAlive[0] = 1;
	numAlive = 1;
	if (cycleShape != 1.0)
		{
		Schedule (0, 0.0);
		for (;;)
			{
			while (heapN > 0 && !(nodeCell[heap[0].node] >= 0 && nodeNext[heap[0].node] == heap[0].t))	/* events of cells that have changed */
				HeapPop ();
			if (heapN == 0)
				return 0;
			t = heap[0].t;
			if (nextPreset < presetN && t >= presetTime[nextPreset])
				{
				t = presetTime[nextPreset];
				i = StartPresetClone (nextPreset, t);
				if (i < 0)
					return 0;
				nextPreset++;
				Schedule (i, t);				/* the cell starts its cycle again with the rates of the clone */
				continue;
				}
			if (stopTime > 0.0 && t >= stopTime)
				{
				t = stopTime;
				break;
				}
			i = nodeCell[heap[0].node];
			k = nodeDeath[heap[0].node];
			HeapPop ();
			if (k)
				{
				Death (i);
				if (numAlive == 0)
					return 0;
				}
			else
				{
				Division (i, t);
				n = numAlive - 1;
				if (stopSize > 0.0 && (double) numAlive >= stopSize)
					break;
				if ((double) numAlive > MAX_CELLS)
					Fail ("the tumor has more than 5e7 cells (every division is kept in memory): lower --time or the growth rate%s", "");
				Schedule (i, t);
				Schedule (n, t);
				}
			}
		}
	else
	for (;;)
		{
		rate = bMax + dMax;
		t += Exponential ((double) numAlive * rate);
		if (nextPreset < presetN && t >= presetTime[nextPreset])	/* a clone of the user starts here */
			{
			t = presetTime[nextPreset];
			if (StartPresetClone (nextPreset, t) < 0)
				return 0;
			nextPreset++;
			continue;
			}
		if (stopTime > 0.0 && t >= stopTime)
			{
			t = stopTime;
			break;
			}
		i = (long) (Uniform () * (double) numAlive);
		k = cellClone[i];
		bk = birthRate * cloneRate[k];
		r = Uniform () * rate;
		if (r < bk)					/* division */
			{
			Division (i, t);
			if (stopSize > 0.0 && (double) numAlive >= stopSize)
				break;
			if ((double) numAlive > MAX_CELLS)
				Fail ("the tumor has more than 5e7 cells (every division is kept in memory): lower --time or the growth rate%s", "");
			}
		else if (r < bk + cloneDeath[k])		/* death */
			{
			Death (i);
			if (numAlive == 0)
				return 0;
			}
		}
	*endTime = t;
	if (nextPreset < presetN)			/* the tumor reached its size before a clone of the user started */
		return 0;
	return numAlive >= sampleSize;
	}

/* ------------------------------ the tree of the sample ------------------------------ */

static int32_t	*sCount, *sChildren, *sKeptAnc, *sKid0, *sKid1;
static double	*sDiv;
static unsigned char *sLeaf;
static int32_t	*leafOfNode;		/* number of the sampled cell (1 to K) of a leaf node */
static double	treeEnd;

static int32_t	*sTop;			/* for a kept node: the first node of its branch (the child of the kept ancestor) */

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
		if (mutatorsOn)
			fprintf (fp, ":%.10g", nodeMutEnd[node] - nodeMutStart[sTop[node]]);
		else
			fprintf (fp, ":%.10g", (h - hp) * mutationRate + (double) (nodeDivs[node] - nodeDivs[parentKept]) * divisionMutationRate);
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
	sTop = (int32_t*) Checked (calloc (n, sizeof (int32_t)));
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
		if (mutatorsOn)
			{
			SyncMutations (v, treeEnd);
			nodeMutEnd[cellNode[v]] = cellBase[v];
			}
		if (fpClones != NULL)
			{
			k = cellClone[v];
			fprintf (fpClones, "%d\tcell%04ld\t%d\t%d\t%.6g\t%.6g\n", rep, i + 1, k, cloneDepth[k], cloneRate[k], mutatorsOn ? cellMult[v] : 1.0);
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
	for (i=1; i<n; i++)
		sTop[i] = (sLeaf[nodeParent[i]] || sChildren[nodeParent[i]] >= 2) ? (int32_t) i : sTop[nodeParent[i]];
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
	free (sTop);
	free (pick);
	return tmrca;
	}

/* ------------------------------ main ------------------------------ */

/* ------------------------------ the tree of the sample, without the tumor ------------------------------ */

/* A neutral tumor (birth b, death d, exponential waiting times, no drivers) is a birth-death process, and the tree of the cells that are
   sampled can be drawn directly, with a cost that does not depend on the size of the tumor:
   - the reconstructed tree of the cells alive at the time T is a coalescent point process (Lambert and Stadler 2013): the depths of the
     nodes are independent, with P(H > t) = 1 / F(t), F(t) = 1 + (b / g) (e^{g t} - 1), g = b - d;
   - at a fixed time T (--time), the tree of s cells chosen at random is exact (Lambert 2018): Y is drawn from a density s (1 - a) y^(s-1) /
     (1 - a (1 - y))^(s+1), with a = 1 - 1 / F(T), and then s - 1 independent depths with density f_y / c on (0, T), where F_y = 1 - y + y F;
   - when the tumor stops at N cells (--size): the time T is that of a chain of the number of cells (exact up to 2e4 cells, and then the
     exponential growth), and the s cells are chosen among the N tips of the coalescent point process: the depth between two chosen tips is
     the maximum of the depths of the tips between them;
   - the divisions of the lineages that do not leave a sampled cell (the other daughter has no sampled cells, or dies out) happen along the
     branches at the rate 2 b (1 - p(t)), p(t) = 1 / (b / g + (1 / y - b / g) e^{-g t}) being the probability that a cell at the time t before the
     end has a sampled descendant, when every cell is sampled with probability y (s / N for --size). Every branch has these divisions and the
     one of its first node. */
/* can the tree of the sample be drawn without simulating the tumor? */
static int FastOK (void)
	{
	int		ok = (driverRate == 0.0 && presetN == 0 && !mutatorsOn && cycleShape == 1.0 && birthRate > deathRate && sampleSize >= 2);

	if (treeMethod == 2 && !ok)
		Fail ("--tree-method reconstructed needs a neutral tumor that grows: no drivers, no clones, no hypermutators, --cycle-shape 1 and --birth greater than --death%s", "");
	return treeMethod != 1 && ok;
	}

static double Poisson (double mean)
	{
	double	l, p;
	long	k, total = 0;

	while (mean > 20.0)
		{
		l = exp (-20.0);
		for (k=0, p=1.0; ; k++)
			{
			p *= Uniform ();
			if (p <= l)
				break;
			}
		total += k;
		mean -= 20.0;
		}
	l = exp (-mean);
	for (k=0, p=1.0; ; k++)
		{
		p *= Uniform ();
		if (p <= l)
			break;
		}
	return (double) (total + k);
	}

/* time of the first passage of the number of cells to N (chain of the births and deaths of the cells) */
static double TimeToSize (double N)
	{
	double	t, M = (N < 2e4) ? N : 2e4, g = birthRate - deathRate, tot = birthRate + deathRate;
	long	n;

	for (;;)
		{
		n = 1;
		t = 0.0;
		while (n > 0 && (double) n < M)
			{
			t += Exponential ((double) n * tot);
			n += (Uniform () * tot < birthRate) ? 1 : -1;
			}
		if (n > 0)
			break;
		}
	if (N > M)
		t += log (N / M) / g;
	return t;
	}

static double FofT (double t)			/* F(t) - 1 */
	{
	return birthRate / (birthRate - deathRate) * expm1 ((birthRate - deathRate) * t);
	}

static double TofF1 (double f1)			/* t such that F(t) - 1 = f1 */
	{
	double	g = birthRate - deathRate;

	return log1p (g / birthRate * f1) / g;
	}

static int		*fKid0, *fKid1, *fLeafLabel;
static double	*fDepth, *fHidden;

static void WriteNewickFast (FILE *fp, int node, int parent)
	{
	if (node < sampleSize)
		fprintf (fp, "cell%04d", fLeafLabel[node] + 1);
	else
		{
		fputc ('(', fp);
		WriteNewickFast (fp, fKid0[node - sampleSize], node);
		fputc (',', fp);
		WriteNewickFast (fp, fKid1[node - sampleSize], node);
		fputc (')', fp);
		}
	if (parent >= 0)
		fprintf (fp, ":%.10g", (fDepth[parent] - fDepth[node]) * mutationRate + (1.0 + fHidden[node]) * divisionMutationRate);
	}

/* one replicate without simulating the tumor */
static void FastTree (FILE *fpTrees, FILE *fpClones, FILE *fpDrivers, FILE *fpSummary, int rep, int attempts)
	{
	int		k = sampleSize, i, j, top, *stack, *left, *right, root, m, tmp;
	double	T, N = -1.0, g = birthRate - deathRate, b = birthRate, F1T, oneMinusA, a, y, c, w, v, u, *x, A, B, tau1, tau2, mean, rootDepth;
	long	*pos, count, cell;

	if (stopSize > 0.0)
		{
		N = stopSize;
		T = TimeToSize (N);
		}
	else
		T = stopTime;
	F1T = FofT (T);
	oneMinusA = 1.0 / (1.0 + F1T);					/* 1 - a, with a = 1 - 1 / F(T) */
	a = F1T / (1.0 + F1T);
	x = (double*) Checked (calloc (k, sizeof (double)));
	if (N < 0.0)							/* a fixed time: the sample is exact (Lambert 2018) */
		{
		v = pow (Uniform (), 1.0 / (double) k);
		y = v * oneMinusA / ((1.0 - v) + v * oneMinusA);
		w = oneMinusA + a * y;
		c = a * y / w;
		for (i=0; i<k-1; i++)
			{
			u = Uniform ();
			x[i] = TofF1 (u * c / (y * (1.0 - u * c)));
			}
		}
	else								/* the number of cells at the end is N: k tips chosen among the N tips of the coalescent point process */
		{
		y = (double) k / N;
		pos = (long*) Checked (calloc (k, sizeof (long)));
		for (count=0, cell=(long) N - k; cell<(long) N; cell++)		/* Floyd's sample of k of the N positions, kept in order */
			{
			long	r = (long) (Uniform () * (double) (cell + 1)), lo = 0, hi = count, mid;

			while (lo < hi)
				{
				mid = (lo + hi) / 2;
				if (pos[mid] < r)
					lo = mid + 1;
				else
					hi = mid;
				}
			if (lo < count && pos[lo] == r)			/* taken: the largest position is added */
				lo = count, r = cell;
			memmove (pos + lo + 1, pos + lo, (size_t) (count - lo) * sizeof (long));
			pos[lo] = r;
			count++;
			}
		for (i=0; i<k-1; i++)
			{
			double	gap = (double) (pos[i+1] - pos[i]), u1 = pow (Uniform (), 1.0 / gap);	/* maximum of gap depths: the distribution function to the power gap */

			/* P(H < t | H < T) = (1 - 1/F(t)) / (1 - 1/F(T)) = u1, so F(t) - 1 = u1 a / (1 - u1 a) */
			x[i] = TofF1 (u1 * a / ((1.0 - u1) + u1 * oneMinusA));
			}
		free (pos);
		}
	/* the tree: a Cartesian tree on the depths (leaves 0 ... k-1, internal nodes k ... 2k-2) */
	stack = (int*) Checked (calloc (k, sizeof (int)));
	left = (int*) Checked (calloc (k, sizeof (int)));
	right = (int*) Checked (calloc (k, sizeof (int)));
	fKid0 = (int*) Checked (calloc (k, sizeof (int)));
	fKid1 = (int*) Checked (calloc (k, sizeof (int)));
	fDepth = (double*) Checked (calloc (2 * k, sizeof (double)));
	fHidden = (double*) Checked (calloc (2 * k, sizeof (double)));
	fLeafLabel = (int*) Checked (calloc (k, sizeof (int)));
	for (i=0; i<k-1; i++)
		left[i] = right[i] = -1;
	for (i=0, top=0; i<k-1; i++)
		{
		j = -1;
		while (top > 0 && x[stack[top-1]] < x[i])
			j = stack[--top];
		left[i] = j;
		if (top > 0)
			right[stack[top-1]] = i;
		stack[top++] = i;
		}
	root = stack[0];
	for (i=0; i<k-1; i++)
		{
		fKid0[i] = (left[i] >= 0) ? k + left[i] : i;
		fKid1[i] = (right[i] >= 0) ? k + right[i] : i + 1;
		fDepth[k + i] = x[i];
		}
	rootDepth = x[root];
	/* the divisions that leave no sampled cell, on every branch */
	A = b / g;
	B = 1.0 / y - A;
	for (i=0; i<k-1; i++)
		for (j=0; j<2; j++)
			{
			m = (j == 0) ? fKid0[i] : fKid1[i];
			tau1 = fDepth[m];
			tau2 = x[i];
			mean = 2.0 * b * (tau2 - tau1) - 2.0 * ((g * tau2 + log (A + B * exp (-g * tau2))) - (g * tau1 + log (A + B * exp (-g * tau1))));
			fHidden[m] = (mean > 0.0) ? Poisson (mean) : 0.0;
			}
	for (i=0; i<k; i++)						/* the leaves get the labels in random order */
		fLeafLabel[i] = i;
	for (i=k-1; i>0; i--)
		{
		j = (int) (Uniform () * (double) (i + 1));
		tmp = fLeafLabel[i];
		fLeafLabel[i] = fLeafLabel[j];
		fLeafLabel[j] = tmp;
		}
	WriteNewickFast (fpTrees, k + root, -1);
	fputs (";\n", fpTrees);
	for (i=0; i<k; i++)
		fprintf (fpClones, "%d\tcell%04d\t0\t0\t1\t1\n", rep, i + 1);
	if (N >= 0.0)
		fprintf (fpDrivers, "%d\t0\t-1\t0\t0\t1\t%.0f\t%d\t%.6g\n", rep, N, k, deathRate);
	else
		fprintf (fpDrivers, "%d\t0\t-1\t0\t0\t1\tNA\t%d\t%.6g\n", rep, k, deathRate);
	if (N >= 0.0)
		fprintf (fpSummary, "%d\t%d\t%.6g\t%.0f\t1\t0\t%.6g\n", rep, attempts, T, N, rootDepth);
	else
		fprintf (fpSummary, "%d\t%d\t%.6g\tNA\t1\t0\t%.6g\n", rep, attempts, T, rootDepth);
	free (stack); free (left); free (right); free (fKid0); free (fKid1); free (fDepth); free (fHidden); free (fLeafLabel); free (x);
	}


/* one replicate: simulates the tumor and writes its tree and tables; returns the number of tumors that were simulated (more than one if they die out) */
static long RunReplicate (int rep, FILE *fpTrees, FILE *fpClones, FILE *fpDrivers, FILE *fpSummary)
	{
	int		attempts = 0, i, usedClones, *sampledClones, drivers;
	double	t, tmrca, mean;

	SeedRandom (seed + (uint64_t) rep * 0x9E3779B97F4A7C15ULL);		/* the replicate does not depend on the others, nor on the number of jobs */
	if (FastOK ())
		{
		FastTree (fpTrees, fpClones, fpDrivers, fpSummary, rep, 1);
		return 1;
		}
	if (stopSize > MAX_CELLS)
		Fail ("--size above 5e7 cells (every division is kept in memory) is only possible for a neutral tumor that grows, with --tree-method auto or reconstructed%s", "");
	while (!GrowTumor (&t))
		if (++attempts >= 100000)
			Fail ("the tumor dies out, or has fewer cells than the sample, in 100000 attempts: raise --birth, or lower --death or -s%s", "");
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
			fprintf (fpDrivers, "%d\t%d\t%d\t%d\t%.6g\t%.6g\t%ld\t%d\t%.6g\n", rep, i, cloneParent[i], cloneDepth[i], cloneTime[i], cloneRate[i], cloneAlive[i], sampledClones[i], cloneDeath[i]);
		}
	mean = t - tmrca;
	fprintf (fpSummary, "%d\t%d\t%.6g\t%ld\t%ld\t%d\t%.6g\n", rep, attempts + 1, t, numAlive, numClones, drivers, mean);
	free (sampledClones);
	return attempts + 1;
	}

static const char *kinds[4] = {"trees.nwk", "clones.tsv", "drivers.tsv", "summary.tsv"};

static FILE *Open (const char *name, const char *mode)
	{
	char	path[4096];
	FILE	*fp;

	snprintf (path, sizeof (path), "%s/%s", outDir, name);
	fp = fopen (path, mode);
	if (fp == NULL)
		Fail ("cannot write in the folder %s", outDir);
	return fp;
	}

int main (int argc, char **argv)
	{
	FILE	*fp[4], *fpLog, *part[4], *in;
	char	name[4096], buf[65536];
	int		rep, i, w, status, failed = 0;
	long	totalAttempts = 0, got;
	pid_t	*pid;

	if (argc < 2)
		{
		Usage (stderr);
		return 1;
		}
	ReadOptions (argc, argv);
	if (jobs < 0)
		Fail ("-j is a number of jobs, 0 or more%s", "");
	if (jobs == 0)
		jobs = (int) sysconf (_SC_NPROCESSORS_ONLN);
	if (jobs < 1)
		jobs = 1;
	if (jobs > replicates)
		jobs = replicates;
	mkdir (outDir, 0777);
	for (i=0; i<4; i++)
		fp[i] = Open (kinds[i], "w");
	fpLog = Open ("log.txt", "w");
	fprintf (fpLog, "%s %s\ncommand:", PROGRAM, VERSION);
	for (i=1; i<argc; i++)
		fprintf (fpLog, " %s", argv[i]);
	fprintf (fpLog, "\nseed: %llu\n", (unsigned long long) seed);
	fclose (fpLog);
	fprintf (fp[1], "replicate\tcell\tclone\tdrivers\tbirth_rate_relative\tmutation_factor\n");
	fprintf (fp[2], "replicate\tclone\tparent_clone\tdrivers\ttime\tbirth_rate_relative\tcells_at_end\tcells_sampled\tdeath_rate\n");
	fprintf (fp[3], "replicate\tattempts\ttime\tcells\tclones\tmax_drivers\ttmrca\n");
	if (jobs == 1)
		{
		for (rep=1; rep<=replicates; rep++)
			totalAttempts += RunReplicate (rep, fp[0], fp[1], fp[2], fp[3]);
		}
	else
		{
		/* every process takes the replicates w+1, w+1+jobs, ... and writes each one in its own small files, which are joined in order */
		fflush (NULL);
		pid = (pid_t*) Checked (calloc (jobs, sizeof (pid_t)));
		for (w=0; w<jobs; w++)
			{
			pid[w] = fork ();
			if (pid[w] < 0)
				Fail ("cannot start the processes of -j%s", "");
			if (pid[w] == 0)
				{
				for (i=0; i<4; i++)
					fclose (fp[i]);
				for (rep=w+1; rep<=replicates; rep+=jobs)
					{
					for (i=0; i<4; i++)
						{
						snprintf (name, sizeof (name), ".part%d.%s", rep, kinds[i]);
						part[i] = Open (name, "w");
						}
					RunReplicate (rep, part[0], part[1], part[2], part[3]);
					for (i=0; i<4; i++)
						fclose (part[i]);
					}
				_exit (0);
				}
			}
		for (w=0; w<jobs; w++)
			{
			waitpid (pid[w], &status, 0);
			if (!WIFEXITED (status) || WEXITSTATUS (status) != 0)
				failed = 1;
			}
		if (failed)
			Fail ("a process of -j failed (see the message above)%s", "");
		for (rep=1; rep<=replicates; rep++)
			for (i=0; i<4; i++)
				{
				snprintf (name, sizeof (name), "%s/.part%d.%s", outDir, rep, kinds[i]);
				in = fopen (name, "r");
				if (in == NULL)
					Fail ("a replicate is missing (a process of -j failed)%s", "");
				while ((got = (long) fread (buf, 1, sizeof (buf), in)) > 0)
					fwrite (buf, 1, got, fp[i]);
				fclose (in);
				if (i == 3)
					{
					int		r, a;
					in = fopen (name, "r");
					if (in && fscanf (in, "%d %d", &r, &a) == 2)
						totalAttempts += a;
					if (in)
						fclose (in);
					}
				remove (name);
				}
		free (pid);
		}
	for (i=0; i<4; i++)
		fclose (fp[i]);
	fprintf (stderr, "%s: %d replicates (%ld tumors simulated), results in %s/\n", PROGRAM, replicates, totalAttempts, outDir);
	return 0;
	}
