/* Regional mutation rates (--region-rates FILE).
   The file has one region per line: chromosome, start, end and a multiplier of the mutation rate (BED coordinates: the start is
   0-based and the end is not included; fields are separated by blanks or tabs; lines that start with # are comments). The sites that
   are not in a region have multiplier 1. It needs a chromosome with coordinates (--blocks or --fasta). The multipliers are divided by their mean over the chromosome, so the average rate does not change, only where the mutations are: a region with multiplier 2 has twice the
   mutations of a site outside it. With blocks, a block in a hot region has more mutations than the others. */
#include "cellcoal.h"

char	*regionRatesFile = "";
int		regionRatesOn = NO;
double	regionRmax = 1.0;

typedef struct { long s, e; double m; } Region;		/* first and last position (1-based, included) and multiplier */

static Region	*reg = NULL;
static int		rn = 0, loaded = 0;
static double	regMean = 1.0;					/* mean multiplier over the range of the simulation */

static int CmpRegion (const void *a, const void *b)
	{
	long x = ((const Region*) a)->s, y = ((const Region*) b)->s;

	return (x > y) - (x < y);
	}

static void Load (void)
	{
	FILE	*fp;
	char	line[4096], chrom[256], want[256], extra[64];
	int		isPipe, cap = 0, ln = 0, k;
	int		filter = (strlen (fastaFile) > 0);

	snprintf (want, sizeof(want), "%s", BlocksChrom ());
	fp = OpenText (regionRatesFile, &isPipe);
	while (fgets (line, sizeof(line), fp))
		{
		long	a, b;
		double	m;
		char	*p = line;

		ln++;
		while (*p == ' ' || *p == '\t')
			p++;
		if (*p == '#' || *p == '\n' || *p == '\r' || *p == '\0' || strncmp (p, "track", 5) == 0 || strncmp (p, "browser", 7) == 0)
			continue;
		if (sscanf (p, "%255s %ld %ld %lf %63s", chrom, &a, &b, &m, extra) < 4)
			{
			fprintf (stderr, "\nERROR: %s, line %d: expected \"chromosome start end multiplier\"\n", regionRatesFile, ln);
			exit (-1);
			}
		if (a < 0 || b <= a || m < 0.0)
			{
			fprintf (stderr, "\nERROR: %s, line %d: need 0 <= start < end and multiplier >= 0\n", regionRatesFile, ln);
			exit (-1);
			}
		if (filter)
			{
			const char *c1 = chrom, *c2 = want;

			if (strncmp (c1, "chr", 3) == 0) c1 += 3;
			if (strncmp (c2, "chr", 3) == 0) c2 += 3;
			if (strcmp (c1, c2) != 0)
				continue;
			}
		if (rn == cap)
			{
			cap = cap ? 2 * cap : 256;
			reg = (Region*) Checked (realloc (reg, cap * sizeof(Region)), "Could not allocate the regional rates\n");
			}
		reg[rn].s = a + 1;
		reg[rn].e = b;
		reg[rn].m = m;
		rn++;
		}
	CloseText (fp, isPipe);
	if (rn == 0)
		{
		fprintf (stderr, "\nERROR: %s has no regions%s\n", regionRatesFile, filter ? " of the chromosome of the FASTA" : "");
		exit (-1);
		}
	qsort (reg, rn, sizeof(Region), CmpRegion);
	for (k=1; k<rn; k++)
		if (reg[k].s <= reg[k-1].e)
			{
			fprintf (stderr, "\nERROR: %s: the regions %ld-%ld and %ld-%ld overlap\n", regionRatesFile, reg[k-1].s, reg[k-1].e, reg[k].s, reg[k].e);
			exit (-1);
			}
	}

/* mean multiplier over the positions lo..hi (the sites outside the regions have 1) */
static double MeanOver (long lo, long hi)
	{
	double	sum = 0.0;
	long	len = hi - lo + 1;
	int		k;

	for (k=0; k<rn; k++)
		{
		long a = reg[k].s > lo ? reg[k].s : lo, b = reg[k].e < hi ? reg[k].e : hi;

		if (b >= a)
			sum += (double) (b - a + 1) * (reg[k].m - 1.0);
		}
	return 1.0 + sum / (double) len;
	}

/* to be called at the start of every data set (the number of sites and the block can change) */
void RegionRatesSetup (void)
	{
	double	mx = 1.0;
	int		k;

	if (regionRatesOn == NO)
		return;
	if (!loaded)
		{
		Load ();
		loaded = 1;
		}
	regMean = MeanOver (BlocksRegionStart (), chromLength);
	if (regMean <= 0.0)
		{
		fprintf (stderr, "\nERROR: %s: the mean mutation rate of the chromosome is 0 (all the multipliers are 0)\n", regionRatesFile);
		exit (-1);
		}
	for (k=0; k<rn; k++)
		if (reg[k].m > mx)
			mx = reg[k].m;
	regionRmax = mx / regMean;
	}

/* relative rate of a site (mean 1 over the chromosome) */
double RegionRel (int site)
	{
	long	pos = BlockSitePos (site);
	int		lo = 0, hi = rn - 1;

	while (lo <= hi)
		{
		int mid = (lo + hi) / 2;

		if (pos < reg[mid].s)
			hi = mid - 1;
		else if (pos > reg[mid].e)
			lo = mid + 1;
		else
			return reg[mid].m / regMean;
		}
	return 1.0 / regMean;
	}

/* mean relative rate of a list of sites: the number of mutations of a model is proportional to it */
double RegionMean (const int *sites, int n)
	{
	double	sum = 0.0;
	int		i;

	if (n == 0)
		return 1.0;
	for (i=0; i<n; i++)
		sum += RegionRel (sites[i]);
	return sum / (double) n;
	}

/* a site of the list with probability proportional to its relative rate (rejection sampling) */
int RegionPickSite (const int *sites, int n, long int *seed)
	{
	long	trials = 0;

	for (;;)
		{
		int i = RandomUniformTo (n, seed);

		if (regionRatesOn == NO || RandomUniform (seed) * regionRmax < RegionRel (sites[i]))
			return i;
		if (++trials > 100000000L)
			{
			fprintf (stderr, "\nERROR: cannot find a site with a positive mutation rate (--region-rates)\n");
			exit (-1);
			}
		}
	}

/* rate multipliers of the sites: gamma (-a) times the regional rate */
void SiteRateMultipliers (long int *seed)
	{
	int		i, gam = (alphaSites != infinity);

	if (rateVarAmongSites == NO)
		return;
	RegionRatesSetup ();
	for (i=0; i<numSites; i++)
		{
		double v = gam ? RandomGamma (alphaSites, seed) / alphaSites : 1.0;

		if (regionRatesOn == YES)
			v *= RegionRel (i);
		allSites[i].rateMultiplier = v;
		}
	}

void RegionRatesInfo (FILE *fp)
	{
	if (regionRatesOn == YES)
		fprintf (fp, "\n Regional mutation rates                      =   %s", regionRatesFile);
	}
