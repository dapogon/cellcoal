/* Targeted sequencing (--targets FILE, --off-target F), in the reads mode: a whole-exome (or panel) experiment instead of a whole-genome one.
   The file has one region per line: chromosome, start and end (BED coordinates: the start is 0-based and the end is not included; extra fields are
   ignored; lines that start with # are comments). The depth -C is the depth of the target bases; a fragment that overlaps no target is kept with
   probability F (the depth off target, relative to the targets: 0.05 by default, 0 for a perfect capture). The fragments that only partly overlap
   a target are kept, so the depth falls at the edges of the targets as it does with a capture. Without --targets the whole genome is sequenced. */
#include "cellcoal.h"

char	*targetsFile = "";
double	offTarget = 0.05;

typedef struct { long s, e; } Target;		/* first and last position (0-based, included) */

static Target	*tg = NULL;
static int		tgN = 0, tgLoaded = 0;
static long		tgReach = 0;

static int CmpTarget (const void *a, const void *b)
	{
	long x = ((const Target*) a)->s, y = ((const Target*) b)->s;

	return (x > y) - (x < y);
	}

static void TargetsLoad (void)
	{
	FILE	*fp;
	char	line[4096], chrom[256], want[256];
	int		isPipe, cap = 0, ln = 0, k, merged = 0;
	int		filter = (strlen (fastaFile) > 0);

	snprintf (want, sizeof(want), "%s", BlocksChrom ());
	fp = OpenText (targetsFile, &isPipe);
	while (fgets (line, sizeof(line), fp))
		{
		long	a, b;
		char	*p = line;

		ln++;
		while (*p == ' ' || *p == '\t')
			p++;
		if (*p == '#' || *p == '\n' || *p == '\r' || *p == '\0' || strncmp (p, "track", 5) == 0 || strncmp (p, "browser", 7) == 0)
			continue;
		if (sscanf (p, "%255s %ld %ld", chrom, &a, &b) < 3)
			{
			fprintf (stderr, "\nERROR: %s, line %d: expected \"chromosome start end\"\n", targetsFile, ln);
			exit (-1);
			}
		if (a < 0 || b <= a)
			{
			fprintf (stderr, "\nERROR: %s, line %d: need 0 <= start < end\n", targetsFile, ln);
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
		if (tgN == cap)
			{
			cap = cap ? 2 * cap : 256;
			tg = (Target*) Checked (realloc (tg, cap * sizeof(Target)), "Could not allocate the targets\n");
			}
		tg[tgN].s = a;
		tg[tgN].e = b - 1;
		tgN++;
		}
	CloseText (fp, isPipe);
	if (tgN == 0)
		{
		fprintf (stderr, "\nERROR: %s has no regions%s\n", targetsFile, filter ? " of the chromosome of the FASTA" : "");
		exit (-1);
		}
	qsort (tg, tgN, sizeof(Target), CmpTarget);
	for (k=1; k<tgN; k++)		/* overlapping targets are joined */
		if (tg[k].s <= tg[merged].e + 1)
			{
			if (tg[k].e > tg[merged].e)
				tg[merged].e = tg[k].e;
			}
		else
			tg[++merged] = tg[k];
	tgN = merged + 1;
	}

/* to be called once before the reads */
void TargetsSetup (void)
	{
	if (strlen (targetsFile) == 0 || tgLoaded)
		return;
	TargetsLoad ();
	tgLoaded = 1;
	tgReach = (long) (insertMean + 6.0 * insertSd) + readLen + 1;
	}

/* first target that ends at or after position a */
static int FirstFrom (long a)
	{
	int lo = 0, hi = tgN;

	while (lo < hi)
		{
		int mid = (lo + hi) / 2;

		if (tg[mid].e < a)
			lo = mid + 1;
		else
			hi = mid;
		}
	return lo;
	}

/* does a target overlap the positions a .. b (0-based, included) of the chromosome? */
int TargetHit (long a, long b)
	{
	int k = FirstFrom (a);

	return k < tgN && tg[k].s <= b;
	}

/* can a fragment that starts in a .. b hit a target? (false: the whole window has the off-target depth) */
int TargetNear (long a, long b)
	{
	return TargetHit (a - tgReach, b + tgReach);
	}
