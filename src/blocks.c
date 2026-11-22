/* Chromosome in blocks (--blocks, --fasta, --region, --germline-vcf).

   One coalescent tree is generated once (first block) and every block of sites is simulated on that same tree: with no recombination the
   sites are independent given the tree. The output is that of one chromosome: one VCF with the positions of the chromosome and one summary
   row per block plus a TOTAL row. The processes along the genome (--ado-block, --gain-sd) continue from one block to the next.

   With --fasta the reference is a real sequence (plain or gzip FASTA, several lines per sequence): the sites that are not A, C, G or T are left
   out (they are not simulated, not in the VCF and not in the invariant counts) and the positions of the VCF are the real ones.
   --germline-vcf gives the germline SNVs (maternal/paternal alleles of the healthy genome); unphased heterozygous genotypes are phased at random. */

#include "cellcoal.h"
#include <ctype.h>
#include <errno.h>

int			blocksN = 0, blockIndex = 0, blocksOn = 0;
char		*fastaFile = "", *regionSpec = "", *germlineVcf = "";
long		blockOffset = 0, blockSize = 0, chromLength = 0, maskedSites = 0, germlineUsed = 0, germlineRefMismatch = 0;
int			*posMap = NULL;
char		*refUserGenome = NULL;	/* reference base of every site of the block (the maternal genome may carry the germline allele) */
char		*germFlag = NULL;		/* the site of the block is a germline SNV */

typedef struct { long pos; char ref, alt; char gt[4]; } GermSnv;

static char		*regSeq = NULL;
static long		regStart = 1, regLen = 0, nBlocksTotal = 0;
static char		chromName[MAX_NAME] = "1";
static GermSnv	*germ = NULL;
static long		numGerm = 0;

FILE *OpenText (const char *path, int *isPipe)
	{
	size_t	n = strlen (path);
	FILE	*fp;

	*isPipe = 0;
	if (n > 3 && strcmp (path + n - 3, ".gz") == 0)
		{
		char	cmd[2*MAX_NAME+32];

		if (strchr (path, '\'') != NULL)
			{
			fprintf (stderr, "\nERROR: cannot read a compressed file whose name has a quote: %s\n", path);
			exit (-1);
			}
		snprintf (cmd, sizeof(cmd), "gzip -dc '%s'", path);
		fp = popen (cmd, "r");
		*isPipe = 1;
		}
	else
		fp = fopen (path, "r");
	if (fp == NULL)
		{
		fprintf (stderr, "\nERROR: cannot open %s\n", path);
		exit (-1);
		}
	return fp;
	}

void CloseText (FILE *fp, int isPipe)
	{
	if (isPipe)
		pclose (fp);
	else
		fclose (fp);
	}

static int SameChrom (const char *a, const char *name)
	{
	if (strcmp (a, name) == 0)
		return 1;
	if (strncmp (a, "chr", 3) == 0 && strcmp (a + 3, name) == 0)
		return 1;
	if (strncmp (name, "chr", 3) == 0 && strcmp (a, name + 3) == 0)
		return 1;
	return 0;
	}

/* name:start-end (1-based, inclusive) or just a name; returns 1 if there are coordinates */
static int ParseRegion (const char *spec, char *name, long *s, long *e)
	{
	const char	*c = strrchr (spec, ':');
	long		a, b;

	if (c != NULL && sscanf (c + 1, "%ld-%ld", &a, &b) == 2)
		{
		size_t n = (size_t) (c - spec);

		if (n >= MAX_NAME)
			n = MAX_NAME - 1;
		memcpy (name, spec, n);
		name[n] = '\0';
		*s = a; *e = b;
		return 1;
		}
	snprintf (name, MAX_NAME, "%s", spec);
	return 0;
	}

static int CmpGerm (const void *a, const void *b)
	{
	long x = ((const GermSnv*) a)->pos, y = ((const GermSnv*) b)->pos;

	return (x > y) - (x < y);
	}

static void LoadGermline (void)
	{
	FILE	*fp;
	int		isPipe;
	char	*line = NULL;
	size_t	cap = 0;
	long	cnt = 0, alloc = 0;

	fp = OpenText (germlineVcf, &isPipe);
	while (getline (&line, &cap, fp) > 0)
		{
		char	*f[10], *p = line, *gt;
		int		k;
		long	pos;

		if (line[0] == '#')
			continue;
		for (k=0; k<10; k++)
			{
			f[k] = p;
			if (k < 9)
				{
				p = strchr (p, '\t');
				if (p == NULL)
					break;
				*p++ = '\0';
				}
			}
		if (k < 9)
			continue;
		if (!SameChrom (f[0], chromName))
			continue;
		pos = atol (f[1]);
		if (pos < regStart || pos >= regStart + regLen)
			continue;
		if (strlen (f[3]) != 1 || strlen (f[4]) != 1 || !strchr ("ACGT", f[3][0]) || !strchr ("ACGT", f[4][0]))
			continue;
		gt = f[9];
		if (!((gt[0] == '0' || gt[0] == '1') && (gt[1] == '|' || gt[1] == '/') && (gt[2] == '0' || gt[2] == '1')))
			continue;
		if (gt[0] == '0' && gt[2] == '0')
			continue;
		if (cnt == alloc)
			{
			alloc = alloc ? alloc * 2 : 1024;
			germ = (GermSnv*) realloc (germ, alloc * sizeof(GermSnv));
			}
		germ[cnt].pos = pos;
		germ[cnt].ref = f[3][0];
		germ[cnt].alt = f[4][0];
		memcpy (germ[cnt].gt, gt, 3);
		germ[cnt].gt[3] = '\0';
		cnt++;
		}
	free (line);
	CloseText (fp, isPipe);
	if (cnt > 1)
		qsort (germ, cnt, sizeof(GermSnv), CmpGerm);
	numGerm = cnt;
	}

/* reads the reference region and decides the number of blocks; called once, after the options are read */
void BlocksSetup (void)
	{
	blocksOn = (blocksN > 0 || strlen (fastaFile) > 0);
	if (!blocksOn)
		return;
	blockSize = numSites;
	if (blockSize < 1)
		{
		fprintf (stderr, "\nERROR: the block size (-l) has to be positive\n");
		exit (-1);
		}
	if (strlen (fastaFile) > 0)
		{
		FILE	*fp;
		int		isPipe, found = 0, inSeq = 0;
		char	*line = NULL, want[MAX_NAME] = "", name[MAX_NAME];
		size_t	cap = 0;
		long	s = 0, e = 0, h0 = 1, pos = 1, n = 0, alloc = 0;
		int		haveRegion = 0;

		if (strlen (regionSpec) > 0)
			haveRegion = ParseRegion (regionSpec, want, &s, &e);
		fp = OpenText (fastaFile, &isPipe);
		while (getline (&line, &cap, fp) > 0)
			{
			size_t	len = strlen (line);

			while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r'))
				line[--len] = '\0';
			if (line[0] == '>')
				{
				char	*sp = line + 1, *colon;
				long	a, b;

				if (inSeq)
					break;	/* the next sequence: the region was in the previous one */
				for (len=0; sp[len] && !isspace ((unsigned char) sp[len]); len++)
					;
				snprintf (name, sizeof(name), "%.*s", (int) (len < MAX_NAME-1 ? len : MAX_NAME-1), sp);
				h0 = 1;
				colon = strrchr (name, ':');
				if (colon != NULL && sscanf (colon + 1, "%ld-%ld", &a, &b) == 2)
					{
					*colon = '\0';
					h0 = a;
					}
				if (want[0] == '\0' || SameChrom (name, want))
					{
					inSeq = 1;
					found = 1;
					pos = h0;
					snprintf (chromName, sizeof(chromName), "%s", name);
					if (!haveRegion)
						{ s = h0; e = (long) 1 << 40; }
					}
				continue;
				}
			if (!inSeq)
				continue;
			for (len=0; line[len]; len++)
				{
				if (pos >= s && pos <= e)
					{
					if (n == alloc)
						{
						alloc = alloc ? alloc * 2 : 1 << 20;
						regSeq = (char*) realloc (regSeq, alloc);
						if (!regSeq)
							{
							fprintf (stderr, "\nERROR: not enough memory for the reference\n");
							exit (-1);
							}
						}
					regSeq[n++] = (char) toupper ((unsigned char) line[len]);
					}
				pos++;
				}
			if (pos > e)
				break;
			}
		free (line);
		CloseText (fp, isPipe);
		if (!found || n == 0)
			{
			fprintf (stderr, "\nERROR: sequence/region %s not found in %s\n", strlen (regionSpec) ? regionSpec : "(first sequence)", fastaFile);
			exit (-1);
			}
		regStart = haveRegion ? s : h0;
		regLen = n;
		nBlocksTotal = (n + blockSize - 1) / blockSize;
		if (blocksN > 0 && blocksN < nBlocksTotal)
			nBlocksTotal = blocksN;
		blocksN = (int) nBlocksTotal;
		chromLength = regStart + regLen - 1;
		if (strlen (germlineVcf) > 0)
			LoadGermline ();
		}
	else
		{
		if (strlen (germlineVcf) > 0 || strlen (regionSpec) > 0)
			{
			fprintf (stderr, "\nERROR: --region and --germline-vcf need --fasta\n");
			exit (-1);
			}
		chromLength = (long) blocksN * blockSize;
		}
	numDataSets = blocksN;
	}

/* prepares block b: its genome (with --fasta) and its coordinates. Returns the number of sites (0 = nothing to simulate). */
int BlockPrepare (int b, long int *seed)
	{
	blockIndex = b;
	blockOffset = (long) b * blockSize;
	adoOrigin = blockOffset;
	if (strlen (fastaFile) == 0)
		{
		numSites = (int) blockSize;
		return numSites;
		}
	{
	long	off = blockOffset, len = regLen - off, i, kept = 0, g = 0;
	char	*mat, *pat;

	if (len > blockSize)
		len = blockSize;
	free (posMap); free (refUserGenome); free (germFlag);
	posMap = (int*) malloc ((size_t) len * sizeof(int));
	mat = (char*) malloc ((size_t) len + 1);
	pat = (char*) malloc ((size_t) len + 1);
	refUserGenome = (char*) malloc ((size_t) len + 1);
	germFlag = (char*) calloc ((size_t) len + 1, 1);
	if (!posMap || !mat || !pat || !refUserGenome || !germFlag)
		{
		fprintf (stderr, "\nERROR: not enough memory for the block %d\n", b + 1);
		exit (-1);
		}
	while (g < numGerm && germ[g].pos < regStart + off)
		g++;
	for (i=0; i<len; i++)
		{
		char	ch = regSeq[off + i];
		long	p = regStart + off + i;
		char	m = ch, q = ch;

		while (g < numGerm && germ[g].pos < p)
			g++;
		if (!strchr ("ACGT", ch) || ch == '\0')
			continue;
		if (g < numGerm && germ[g].pos == p)
			{
			if (germ[g].ref != ch)
				germlineRefMismatch++;
			else
				{
				const char *gt = germ[g].gt;
				int a, c;

				if (gt[1] == '|')
					{ a = (gt[0] == '1'); c = (gt[2] == '1'); }
				else if (gt[0] == '1' && gt[2] == '1')
					{ a = c = 1; }
				else if (RandomUniform (seed) < 0.5)
					{ a = 1; c = 0; }
				else
					{ a = 0; c = 1; }
				if (a) m = germ[g].alt;
				if (c) q = germ[g].alt;
				germFlag[kept] = 1;
				germlineUsed++;
				}
			}
		mat[kept] = m;
		pat[kept] = q;
		refUserGenome[kept] = ch;
		posMap[kept] = (int) p;
		kept++;
		}
	mat[kept] = '\0';
	pat[kept] = '\0';
	maskedSites += len - kept;
	free (maternalUserGenome);
	free (paternalUserGenome);
	maternalUserGenome = mat;
	paternalUserGenome = pat;
	numSites = (int) kept;
	return numSites;
	}
	}

/* position of the site j of the current block in the chromosome */
long BlockSitePos (int j)
	{
	if (posMap != NULL)
		return posMap[j];
	return blockOffset + j + 1;
	}

void BlocksInfo (FILE *fp)
	{
	fprintf (fp, "blocks: %d of %ld sites; chromosome %s, length %ld\n", blocksN, blockSize, strlen (fastaFile) ? chromName : contigName, chromLength);
	if (strlen (fastaFile) > 0)
		{
		fprintf (fp, "reference: %s %s:%ld-%ld (%ld bp)\n", fastaFile, chromName, regStart, regStart + regLen - 1, regLen);
		fprintf (fp, "sites that are not A, C, G or T (left out): %ld\n", maskedSites);
		if (strlen (germlineVcf) > 0)
			fprintf (fp, "germline SNVs in the region: %ld; used: %ld; REF different from the reference (ignored): %ld\n", numGerm, germlineUsed, germlineRefMismatch);
		}
	}

void BlocksFree (void)
	{
	free (regSeq); regSeq = NULL;
	free (germ); germ = NULL;
	free (posMap); posMap = NULL;
	free (refUserGenome); refUserGenome = NULL;
	free (germFlag); germFlag = NULL;
	}

long BlocksRegionStart (void)
	{
	return regStart;
	}

const char *BlocksChrom (void)
	{
	return (strlen (fastaFile) > 0) ? chromName : contigName;
	}
