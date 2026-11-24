/* Mappability of the genome (--mappability FILE, --map-drop F), in the reads mode.
   The file has one region per line: chromosome, start, end and a mappability score from 0 to 1 (BED coordinates: the start is 0-based and the
   end is not included; fields are separated by blanks or tabs; lines that start with # are comments). A score of 1 is a sequence that is found
   once in the genome and a score of 1/n a sequence with n copies (the scores of GEM, Umap or the tracks of the UCSC genome browser). The sites that
   are not in a region have score 1. The mapping quality (MAPQ) of a read is a function of the mean score of the bases that it covers:
   60 for 1, 0 for 0.5 or less (the read is ambiguous: an aligner cannot decide among two or more places), and a straight line in between. */
#include "cellcoal.h"

char	*mappabilityFile = "";
double	mapDrop = 0.0;

typedef struct { long s, e; double m; } MapReg;		/* first and last position (0-based, included) and score */

static MapReg	*mapReg = NULL;
static int		mapN = 0, mloaded = 0;

static int CmpMap (const void *a, const void *b)
	{
	long x = ((const MapReg*) a)->s, y = ((const MapReg*) b)->s;

	return (x > y) - (x < y);
	}

static void MapLoad (void)
	{
	FILE	*fp;
	char	line[4096], chrom[256], want[256], extra[64];
	int		isPipe, cap = 0, ln = 0, k;
	int		filter = (strlen (fastaFile) > 0);

	snprintf (want, sizeof(want), "%s", BlocksChrom ());
	fp = OpenText (mappabilityFile, &isPipe);
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
			fprintf (stderr, "\nERROR: %s, line %d: expected \"chromosome start end score\"\n", mappabilityFile, ln);
			exit (-1);
			}
		if (a < 0 || b <= a || m < 0.0 || m > 1.0)
			{
			fprintf (stderr, "\nERROR: %s, line %d: need 0 <= start < end and a score from 0 to 1\n", mappabilityFile, ln);
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
		if (mapN == cap)
			{
			cap = cap ? 2 * cap : 256;
			mapReg = (MapReg*) Checked (realloc (mapReg, cap * sizeof(MapReg)), "Could not allocate the mappability\n");
			}
		mapReg[mapN].s = a;
		mapReg[mapN].e = b - 1;
		mapReg[mapN].m = m;
		mapN++;
		}
	CloseText (fp, isPipe);
	if (mapN == 0)
		{
		fprintf (stderr, "\nERROR: %s has no regions%s\n", mappabilityFile, filter ? " of the chromosome of the FASTA" : "");
		exit (-1);
		}
	qsort (mapReg, mapN, sizeof(MapReg), CmpMap);
	for (k=1; k<mapN; k++)
		if (mapReg[k].s <= mapReg[k-1].e)
			{
			fprintf (stderr, "\nERROR: %s: the regions %ld-%ld and %ld-%ld overlap\n", mappabilityFile, mapReg[k-1].s, mapReg[k-1].e, mapReg[k].s, mapReg[k].e);
			exit (-1);
			}
	}

/* to be called once before the reads */
void MappabilitySetup (void)
	{
	if (strlen (mappabilityFile) == 0 || mloaded)
		return;
	MapLoad ();
	mloaded = 1;
	}

/* mean score over the positions a .. b (0-based, included) of the chromosome */
static double MapMean (long a, long b)
	{
	double	sum = 0.0;
	long	len = b - a + 1;
	int		lo = 0, hi = mapN - 1, k;

	while (lo < hi)
		{
		int mid = (lo + hi) / 2;

		if (mapReg[mid].e < a)
			lo = mid + 1;
		else
			hi = mid;
		}
	for (k=lo; k<mapN && mapReg[k].s <= b; k++)
		{
		long x = mapReg[k].s > a ? mapReg[k].s : a, y = mapReg[k].e < b ? mapReg[k].e : b;

		if (y >= x)
			sum += (double) (y - x + 1) * (mapReg[k].m - 1.0);
		}
	return 1.0 + sum / (double) len;
	}

/* mapping quality of a read that covers the positions a .. b (0-based) of the chromosome */
int MapQual (long a, long b)
	{
	double	m;

	if (!mloaded)
		return 60;
	m = MapMean (a, b);
	if (m >= 1.0)
		return 60;
	if (m <= 0.5)
		return 0;
	return (int) (60.0 * (m - 0.5) / 0.5 + 0.5);
	}
