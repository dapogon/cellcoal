/* Small insertions and deletions in the read-level simulation (--indel-rate ...). The indels are events of the genealogy: they happen on the
   branches of the tree, on one haplotype, at a position of the chromosome of the block, and every cell below the branch inherits them. They live in
   a layer of their own, apart from the site alleles: the read simulator applies the indels of a cell and haplotype to the fragments that overlap them
   (the bases that are deleted disappear, the inserted bases are new) and writes the CIGAR of the reads. The truth of the indels is written to
   indels_truth.vcf. Every indel is left-aligned against the reference, as a variant caller would report it. */
#include "cellcoal.h"

double		indelRate = 0.0, indelIns = 0.5, indelLen = 2.0;
int			indelMax = 10;

static IndelEv	*master = NULL;			/* all the events of the block */
static long		nMaster = 0, capMaster = 0;
static IndelEv	**cellEv[2] = {NULL, NULL};	/* resolved events of every cell and haplotype, sorted by position */
static long		*cellN[2] = {NULL, NULL};
static int		nCellsAlloc = 0;
static long		nextId = 0;

static char RefAt (long j)
	{
	int a;

	if (j < 0 || j >= numSites)
		return 'N';
	if (refUserGenome != NULL)
		return refUserGenome[j];
	a = data[MATERNAL][OUTGROUP_ROOT][j];
	return (a >= 0 && a < 4) ? "ACGT"[a] : 'N';
	}

static long CP (long j)
	{
	return blocksOn ? BlockSitePos ((int) j) - 1 : j;
	}

static int ValidBase (char c)
	{
	return c == 'A' || c == 'C' || c == 'G' || c == 'T';
	}

static int CmpEv (const void *a, const void *b)
	{
	const IndelEv *x = (const IndelEv*) a, *y = (const IndelEv*) b;

	if (x->pos != y->pos)
		return (x->pos < y->pos) ? -1 : 1;
	return (x->id < y->id) ? -1 : (x->id > y->id);
	}

/* one new event on a branch; returns 1 when it was made */
static int NewEvent (IndelEv *e, long int *seed)
	{
	long	p, k;
	int		i, len;
	double	q = (indelLen > 1.0) ? 1.0 - 1.0 / indelLen : 0.0;

	len = 1;
	while (len < indelMax && RandomUniform (seed) < q)
		len++;
	e->type = (RandomUniform (seed) < indelIns) ? 0 : 1;
	e->len = len;
	p = (long) (RandomUniform (seed) * (numSites - 1));
	for (i=0; i<len; i++)
		e->seq[i] = "ACGT"[(int) (RandomUniform (seed) * 4.0) & 3];
	e->seq[len] = '\0';
	if (p < 0 || p >= numSites - 1 || (e->type == 1 && p + len >= numSites))
		return 0;
	/* left-align: the same indel is described by the leftmost position (a deletion in a repeat, an insertion in a homopolymer) */
	if (e->type == 1)
		{
		for (k=0; k<len; k++)
			e->seq[k] = RefAt (p + 1 + k);
		while (p > 0 && RefAt (p) == e->seq[len - 1] && ValidBase (RefAt (p)))
			{
			memmove (e->seq + 1, e->seq, len - 1);
			e->seq[0] = RefAt (p);
			p--;
			}
		}
	else
		while (p > 0 && ValidBase (RefAt (p)) && RefAt (p) == e->seq[len - 1])
			{
			char last = e->seq[len - 1];

			memmove (e->seq + 1, e->seq, len - 1);
			e->seq[0] = last;
			p--;
			}
	e->pos = p;
	if (!ValidBase (RefAt (p)))
		return 0;
	if (e->type == 1)
		for (k=1; k<=len; k++)
			if (!ValidBase (RefAt (p + k)))
				return 0;
	e->id = nextId++;
	return 1;
	}

/* the events of a haplotype that overlap are resolved: an event whose anchor lies in a deleted segment is lost with it */
static long Resolve (IndelEv *l, long n)
	{
	long	i, m = 0, from = 1, to = 0;

	qsort (l, n, sizeof(IndelEv), CmpEv);
	for (i=0; i<n; i++)
		{
		if (l[i].pos >= from && l[i].pos <= to)
			continue;
		if (l[i].type == 1)
			{
			from = l[i].pos + 1;
			to = l[i].pos + l[i].len;
			}
		l[m++] = l[i];
		}
	return m;
	}

static void Walk (TreeNode *p, IndelEv *par[2], long npar[2], long int *seed)
	{
	IndelEv	*cur[2];
	long	n[2], cap[2];
	int		h;
	int		own = (p->anc != NULL && !p->isHealthyTip && !p->isHealthyRoot);

	for (h=0; h<2; h++)
		{
		long m = 0, k;

		if (own)
			m = RandomPoisson (p->branchLength * (double) numSites * indelRate, seed);
		cap[h] = npar[h] + m + 1;
		cur[h] = (IndelEv*) Checked (malloc (cap[h] * sizeof(IndelEv)), "Could not allocate the indels\n");
		if (npar[h] > 0)
			memcpy (cur[h], par[h], npar[h] * sizeof(IndelEv));
		n[h] = npar[h];
		for (k=0; k<m; k++)
			if (NewEvent (&cur[h][n[h]], seed))
				{
				if (nMaster == capMaster)
					{
					capMaster = capMaster ? 2 * capMaster : 1024;
					master = (IndelEv*) Checked (realloc (master, capMaster * sizeof(IndelEv)), "Could not allocate the indels\n");
					}
				cur[h][n[h]].hap = h;
				master[nMaster++] = cur[h][n[h]];
				n[h]++;
				}
		}
	if (p->left == NULL && p->right == NULL)
		{
		if (!p->isHealthyTip && p->label >= 0 && p->label < numCells)
			for (h=0; h<2; h++)
				{
				n[h] = Resolve (cur[h], n[h]);
				cellEv[h][p->label] = cur[h];
				cellN[h][p->label] = n[h];
				cur[h] = NULL;
				}
		}
	else
		{
		if (p->left != NULL)
			Walk (p->left, cur, n, seed);
		if (p->right != NULL)
			Walk (p->right, cur, n, seed);
		}
	free (cur[0]);
	free (cur[1]);
	}

void IndelsFree (void)
	{
	int h, i;

	for (h=0; h<2; h++)
		{
		for (i=0; i<nCellsAlloc && cellEv[h] != NULL; i++)
			free (cellEv[h][i]);
		free (cellEv[h]);
		free (cellN[h]);
		cellEv[h] = NULL;
		cellN[h] = NULL;
		}
	free (master);
	master = NULL;
	nMaster = capMaster = 0;
	nCellsAlloc = 0;
	}

/* the events of the block, on the branches of the tree (before the reads of the block) */
void IndelsGenerate (long int *seed)
	{
	int		h;
	IndelEv	*none[2] = {NULL, NULL};
	long	nn[2] = {0, 0};
	long	s1;

	IndelsFree ();
	if (indelRate <= 0.0 || numSites < 3)
		return;
	nCellsAlloc = numCells;
	for (h=0; h<2; h++)
		{
		cellEv[h] = (IndelEv**) Checked (calloc (numCells, sizeof(IndelEv*)), "Could not allocate the indels\n");
		cellN[h] = (long*) Checked (calloc (numCells, sizeof(long)), "Could not allocate the indels\n");
		}
	RandomUniform (seed);
	s1 = *seed;
	nextId = 0;
	{
	long int sd = s1;

	Walk (healthyRoot, none, nn, &sd);
	}
	qsort (master, nMaster, sizeof(IndelEv), CmpEv);
	}

const IndelEv *IndelsOf (int cellIdx, int hap, long *n)
	{
	if (cellIdx < 0 || cellEv[hap] == NULL || cellIdx >= nCellsAlloc)
		{
		*n = 0;
		return NULL;
		}
	*n = cellN[hap][cellIdx];
	return cellEv[hap][cellIdx];
	}

/* indels_truth.vcf: the indels of the block (the left-aligned variant of the reference) and the cells that carry them */
void IndelsWriteVcf (const char *dir)
	{
	static int started = 0;
	char	path[4*MAX_NAME];
	FILE	*fp;
	long	e, *ptr[2];
	int		i, h;

	if (nMaster == 0 && started)
		return;
	snprintf (path, sizeof(path), "%s/indels_truth.vcf", dir);
	if ((fp = fopen (path, (blocksOn && started) ? "a" : "w")) == NULL)
		return;
	if (!(blocksOn && started))
		{
		fprintf (fp, "##fileformat=VCFv4.2\n##source=cellcoal indels truth\n##contig=<ID=%s,length=%ld>\n", ReadContig (), ReadRefLen ());
		fprintf (fp, "##INFO=<ID=NC,Number=1,Type=Integer,Description=\"Cells that carry the indel\">\n");
		fprintf (fp, "##INFO=<ID=TYPE,Number=1,Type=String,Description=\"INS or DEL\">\n");
		fprintf (fp, "##INFO=<ID=LEN,Number=1,Type=Integer,Description=\"Length of the indel\">\n");
		fprintf (fp, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"True genotype (maternal|paternal): 1 when the haplotype carries the indel\">\n");
		fprintf (fp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
		for (i=0; i<numCells; i++)
			fprintf (fp, "\tcell%04d", i + 1);
		fprintf (fp, "\n");
		}
	started = 1;
	for (h=0; h<2; h++)
		ptr[h] = (long*) Checked (calloc (numCells, sizeof(long)), "Could not allocate the indels\n");
	for (e=0; e<nMaster; e++)
		{
		const IndelEv *ev = &master[e];
		int		nc = 0;
		long	k;
		char	*gts = (char*) Checked (malloc (4 * (size_t) numCells + 4), "Could not allocate the indels\n");
		char	*g = gts;

		fprintf (fp, "%s\t%ld\t.\t", ReadContig (), CP (ev->pos) + 1);
		if (ev->type == 1)
			{
			for (k=0; k<=ev->len; k++)
				fputc (RefAt (ev->pos + k), fp);
			fprintf (fp, "\t%c", RefAt (ev->pos));
			}
		else
			fprintf (fp, "%c\t%c%s", RefAt (ev->pos), RefAt (ev->pos), ev->seq);
		for (i=0; i<numCells; i++)
			{
			int has[2] = {0, 0};

			for (h=0; h<2; h++)
				{
				while (ptr[h][i] < cellN[h][i] && CmpEv (&cellEv[h][i][ptr[h][i]], ev) < 0)
					ptr[h][i]++;
				if (ptr[h][i] < cellN[h][i] && cellEv[h][i][ptr[h][i]].id == ev->id)
					has[h] = 1;
				}
			nc += (has[0] || has[1]);
			g += sprintf (g, "\t%d|%d", has[0], has[1]);
			}
		fprintf (fp, "\t.\tPASS\tNC=%d;TYPE=%s;LEN=%d\tGT%s\n", nc, ev->type == 1 ? "DEL" : "INS", ev->len, gts);
		free (gts);
		}
	for (h=0; h<2; h++)
		free (ptr[h]);
	fclose (fp);
	}
