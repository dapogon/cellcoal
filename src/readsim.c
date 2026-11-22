/* Read-level simulator (--bam, --fastq): paired-end reads from the amplified genome of every cell.

   For every cellIdx and haplotype it draws fragments along the genome with a density given by the cellIdx coverage (-C, -J), the amplification
   gain (--gain-sd) and the allelic dropout (the lost segments are not covered).
   The fragments are molecules: errors of the amplification are carried by the molecule and are copied to its PCR duplicates (--dup-mean);
   sequencing errors (-E) are drawn independently in every read. Output: one BAM (+ BAI) and/or FASTQ.gz pair per cellIdx, the reference (ref.fa)
   and, with --bulk-depth, the BAM of the healthy genome. The BAM flag 0x400 marks the PCR duplicates, the tag MI is the molecule and HP the haplotype. */

#include "cellcoal.h"
#include "bgzf.h"
#include <sys/stat.h>
#include <math.h>
#include <string.h>

typedef struct { long pos; size_t off; int len; } RecIdx;

typedef struct
	{
	unsigned char	*buf;
	size_t			n, cap;
	RecIdx			*idx;
	long			ni, capi;
	} RecSet;

typedef struct { long site; int alt; double f; int strand; } AmpEvent;	/* strand: -1 both, 0 forward, 1 reverse */

static const char BASES[] = "ACGT";

static void Grow (unsigned char **b, size_t *cap, size_t need)
	{
	if (need <= *cap)
		return;
	while (*cap < need)
		*cap = (*cap == 0) ? 1 << 20 : *cap * 2;
	*b = (unsigned char*) realloc (*b, *cap);
	if (!*b)
		{
		fprintf (stderr, "\nERROR: not enough memory for the reads of a cell\n");
		exit (-1);
		}
	}

static inline int Allele (int hap, int cellIdx, long j)
	{
	if (cellIdx < 0)
		return data[hap][OUTGROUP_ROOT][j];
	return GetG (hap, cellIdx, (int) j);
	}

static double Normal (long int *seed)
	{
	double u1 = RandomUniform (seed), u2 = RandomUniform (seed);

	if (u1 < 1e-300)
		u1 = 1e-300;
	return sqrt (-2.0 * log (u1)) * cos (6.283185307179586 * u2);
	}

static int Reg2Bin (long beg, long end)
	{
	end--;
	if (beg >> 14 == end >> 14) return (int) (((1 << 15) - 1) / 7 + (beg >> 14));
	if (beg >> 17 == end >> 17) return (int) (((1 << 12) - 1) / 7 + (beg >> 17));
	if (beg >> 20 == end >> 20) return (int) (((1 << 9) - 1) / 7 + (beg >> 20));
	if (beg >> 23 == end >> 23) return (int) (((1 << 6) - 1) / 7 + (beg >> 23));
	if (beg >> 26 == end >> 26) return (int) (((1 << 3) - 1) / 7 + (beg >> 26));
	return 0;
	}

static void Put16 (unsigned char *p, int v) { p[0] = v & 0xFF; p[1] = (v >> 8) & 0xFF; }

/* BAM record of one read (sequence in the orientation of the reference) */
static void AddRecord (RecSet *rs, const char *name, int flag, long pos, long mpos, long tlen, const char *seq, int len, int q, const unsigned char *qv, const char *rg, int hap, int mol)
	{
	size_t	need, o;
	int		nl = (int) strlen (name) + 1, rgl = (int) strlen (rg) + 1, i;
	unsigned char *p;

	need = 4 + 32 + nl + 4 + (len + 1) / 2 + len + 3 + rgl + 3 + 8 + 16;
	Grow (&rs->buf, &rs->cap, rs->n + need);
	if (rs->ni == rs->capi)
		{
		rs->capi = rs->capi ? rs->capi * 2 : 4096;
		rs->idx = (RecIdx*) realloc (rs->idx, rs->capi * sizeof(RecIdx));
		}
	p = rs->buf + rs->n;
	o = 4;
	Put32le (p + o, 0); o += 4;												/* refID */
	Put32le (p + o, (uint32_t) pos); o += 4;
	p[o++] = (unsigned char) nl; p[o++] = 60;								/* l_read_name, mapq */
	Put16 (p + o, Reg2Bin (pos, pos + len)); o += 2;
	Put16 (p + o, 1); o += 2;												/* n_cigar_op */
	Put16 (p + o, flag); o += 2;
	Put32le (p + o, (uint32_t) len); o += 4;
	Put32le (p + o, 0); o += 4;												/* next refID */
	Put32le (p + o, (uint32_t) mpos); o += 4;
	Put32le (p + o, (uint32_t) (int32_t) tlen); o += 4;
	memcpy (p + o, name, nl); o += nl;
	Put32le (p + o, (uint32_t) len << 4); o += 4;							/* CIGAR: len M */
	for (i=0; i<len; i+=2)
		{
		int a = (int) (strchr (BASES, seq[i]) - BASES), b = (i+1 < len) ? (int) (strchr (BASES, seq[i+1]) - BASES) : -1;
		static const int code[4] = {1, 2, 4, 8};

		p[o++] = (unsigned char) ((code[a] << 4) | (b >= 0 ? code[b] : 0));
		}
	if (qv != NULL)
		memcpy (p + o, qv, len);
	else
		memset (p + o, q, len);
	o += len;
	p[o++] = 'R'; p[o++] = 'G'; p[o++] = 'Z'; memcpy (p + o, rg, rgl); o += rgl;
	p[o++] = 'H'; p[o++] = 'P'; p[o++] = 'C'; p[o++] = (unsigned char) hap;
	p[o++] = 'M'; p[o++] = 'I'; p[o++] = 'i'; Put32le (p + o, (uint32_t) mol); o += 4;
	Put32le (p, (uint32_t) (o - 4));
	rs->idx[rs->ni].pos = pos;
	rs->idx[rs->ni].off = rs->n;
	rs->idx[rs->ni].len = (int) o;
	rs->ni++;
	rs->n += o;
	}

static int CmpIdx (const void *a, const void *b)
	{
	const RecIdx *x = (const RecIdx*) a, *y = (const RecIdx*) b;

	if (x->pos != y->pos)
		return x->pos < y->pos ? -1 : 1;
	return x->off < y->off ? -1 : (x->off > y->off);
	}

/* ---------------------------- BAM writer with its index ---------------------------- */

typedef struct { uint64_t beg, end; } Chunk;

/* a BAM file that is written in pieces (the blocks of a chromosome), with its index */
typedef struct
	{
	Bgzf		bz;
	FILE		*fp;
	Chunk		*bins;
	uint64_t	*lin;
	long		nlin, ni, reflen;
	uint64_t	firstOff, lastEnd;
	char		path[4*MAX_NAME+8];
	int			isOpen;
	} BamOut;

static long ReadRefLen (void)		/* length of the reference sequence of the BAM header */
	{
	return blocksOn ? chromLength : (long) numSites;
	}

static const char *ReadContig (void)
	{
	return blocksOn ? BlocksChrom () : contigName;
	}

static void BamOpen (BamOut *bo, const char *path, const char *rg)
	{
	unsigned char l4[4], hdr[8];
	char		text[1024];
	int			tl;
	const char	*cn = ReadContig ();

	memset (bo, 0, sizeof(*bo));
	snprintf (bo->path, sizeof(bo->path), "%s", path);
	if ((bo->fp = fopen (path, "wb")) == NULL)
		{
		fprintf (stderr, "\nERROR: cannot write %s\n", path);
		exit (-1);
		}
	bo->bz.fp = bo->fp;
	CrcInit ();
	bo->reflen = ReadRefLen ();
	bo->nlin = (bo->reflen >> 14) + 2;
	bo->bins = (Chunk*) calloc (37450 + 1, sizeof(Chunk));
	bo->lin = (uint64_t*) calloc (bo->nlin, sizeof(uint64_t));
	tl = snprintf (text, sizeof(text), "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:%s\tLN:%ld\n@RG\tID:%s\tSM:%s\tPL:ILLUMINA\n@PG\tID:cellcoal\tPN:cellcoal\n", cn, bo->reflen, rg, rg);
	memcpy (hdr, "BAM\1", 4);
	BgzfWrite (&bo->bz, hdr, 4);
	Put32le (l4, (uint32_t) tl); BgzfWrite (&bo->bz, l4, 4);
	BgzfWrite (&bo->bz, text, tl);
	Put32le (l4, 1); BgzfWrite (&bo->bz, l4, 4);
	Put32le (l4, (uint32_t) strlen (cn) + 1); BgzfWrite (&bo->bz, l4, 4);
	BgzfWrite (&bo->bz, cn, strlen (cn) + 1);
	Put32le (l4, (uint32_t) bo->reflen); BgzfWrite (&bo->bz, l4, 4);
	bo->isOpen = 1;
	}

/* appends the records of rs (sorted here; their positions have to follow those already written) */
static void BamAppend (BamOut *bo, RecSet *rs)
	{
	unsigned char *buf = rs->buf;
	long		k;

	qsort (rs->idx, rs->ni, sizeof(RecIdx), CmpIdx);
	for (k=0; k<rs->ni; k++)
		{
		RecIdx	*r = &rs->idx[k];
		uint64_t b0 = BgzfTell (&bo->bz), b1;
		long	end = r->pos + (long) (buf[r->off+20] | (buf[r->off+21]<<8) | (buf[r->off+22]<<16) | ((long)buf[r->off+23]<<24));
		int		bin = Reg2Bin (r->pos, end);
		long	w;

		BgzfWrite (&bo->bz, buf + r->off, r->len);
		b1 = BgzfTell (&bo->bz);
		if (bo->ni == 0)
			bo->firstOff = b0;
		bo->ni++;
		bo->lastEnd = b1;
		if (bo->bins[bin].end == 0)
			bo->bins[bin].beg = b0;
		bo->bins[bin].end = b1;
		for (w = r->pos >> 14; w <= (end - 1) >> 14; w++)
			if (bo->lin[w] == 0)
				bo->lin[w] = b0;
		}
	}

static void BamClose (BamOut *bo)
	{
	FILE		*fp;
	int			nbin = 37450 + 1, i;
	char		bai[4*MAX_NAME+16];
	long		k;

	if (!bo->isOpen)
		return;
	BgzfClose (&bo->bz);
	fclose (bo->fp);
	bo->isOpen = 0;
	snprintf (bai, sizeof(bai), "%s.bai", bo->path);
	if ((fp = fopen (bai, "wb")) != NULL)
		{
		int nb = 0;
		unsigned char w4[4];
		uint64_t w8;
		long last = (bo->ni > 0) ? ((bo->reflen - 1) >> 14) + 1 : 0;
		uint64_t prev = 0;

		for (i=0; i<nbin-1; i++)
			if (bo->bins[i].end)
				nb++;
		fwrite ("BAI\1", 1, 4, fp);
		Put32le (w4, 1); fwrite (w4, 1, 4, fp);
		Put32le (w4, (uint32_t) (nb + 1)); fwrite (w4, 1, 4, fp);
		for (i=0; i<nbin-1; i++)
			if (bo->bins[i].end)
				{
				Put32le (w4, (uint32_t) i); fwrite (w4, 1, 4, fp);
				Put32le (w4, 1); fwrite (w4, 1, 4, fp);
				w8 = bo->bins[i].beg; fwrite (&w8, 1, 8, fp);
				w8 = bo->bins[i].end; fwrite (&w8, 1, 8, fp);
				}
		Put32le (w4, 37450); fwrite (w4, 1, 4, fp);
		Put32le (w4, 2); fwrite (w4, 1, 4, fp);
		w8 = bo->firstOff; fwrite (&w8, 1, 8, fp);
		w8 = bo->lastEnd; fwrite (&w8, 1, 8, fp);
		w8 = (uint64_t) bo->ni; fwrite (&w8, 1, 8, fp);
		w8 = 0; fwrite (&w8, 1, 8, fp);
		Put32le (w4, (uint32_t) last); fwrite (w4, 1, 4, fp);
		for (k=0; k<last; k++)
			{
			if (bo->lin[k] == 0)
				bo->lin[k] = prev;
			prev = bo->lin[k];
			fwrite (&bo->lin[k], 1, 8, fp);
			}
		fclose (fp);
		}
	free (bo->bins);
	free (bo->lin);
	bo->bins = NULL; bo->lin = NULL;
	}

/* ------------------------------------ fragments ------------------------------------ */

static int CmpEvent (const void *a, const void *b)
	{
	const AmpEvent *x = (const AmpEvent*) a, *y = (const AmpEvent*) b;

	return (x->site > y->site) - (x->site < y->site);
	}

/* Errors of the amplification of one haplotype: sites where a fraction f of the molecules carry a wrong base */
static AmpEvent *BuildAmpEvents (int cellIdx, int hap, long int *seed, long *nev)
	{
	double	eps = meanAmplificationError, v = varAmplificationError, conc, a, meanSb, rho, maxMult = 1.0;
	long	n, i, cnt = 0, cap = 0;
	AmpEvent *ev = NULL;
	int		x, y, z;

	*nev = 0;
	if (cellIdx < 0 || eps <= 0.0)
		return NULL;
	if (v <= 0.0 || v >= eps * (1.0 - eps))
		v = 0.5 * eps * (1.0 - eps);
	if (cellAmpMean != NULL)	/* --amp-cell-var: the mean of this cell; the concentration of -A is kept */
		{
		double e0 = eps;

		eps = cellAmpMean[cellIdx];
		if (eps <= 0.0)
			return NULL;
		if (eps > 0.95)
			eps = 0.95;
		v = v * eps * (1.0 - eps) / (e0 * (1.0 - e0));
		}
	conc = eps * (1.0 - eps) / v - 1.0;
	a = eps * conc;
	meanSb = (a + 1.0) / (conc + 1.0);			/* an error at a site is seen in a fraction of molecules distributed as Beta(a+1,b) */
	rho = eps / meanSb;
	if (ampSigOn)
		for (x=0; x<4; x++)
			for (y=0; y<4; y++)
				for (z=0; z<4; z++)
					if (ampSigMult[x][y][z] > maxMult)
						maxMult = ampSigMult[x][y][z];
	n = RandomPoisson (rho * maxMult * (double) numSites, seed);
	for (i=0; i<n; i++)
		{
		long	j = (long) (RandomUniform (seed) * numSites);
		int		al = Allele (hap, cellIdx, j), alt, l;
		double	u, cum;

		if (al < 0 || al > 3)
			continue;
		if (RandomUniform (seed) * maxMult > AmpSiteMult (j, al))
			continue;
		u = RandomUniform (seed);
		cum = 0.0;
		alt = -1;
		for (l=0; l<4; l++)
			if (l != al)
				{
				cum += AmpSub (j, al, l);
				if (u <= cum || alt < 0)
					alt = l;
				if (u <= cum)
					break;
				}
		if (cnt == cap)
			{
			cap = cap ? cap * 2 : 1024;
			ev = (AmpEvent*) realloc (ev, cap * sizeof(AmpEvent));
			}
		ev[cnt].site = j;
		ev[cnt].alt = alt;
		ev[cnt].f = RandomBetaMeanVar (meanSb, meanSb * (1.0 - meanSb) / (conc + 2.0), seed);
		ev[cnt].strand = -1;
		if (strandBias > 0.0 && RandomUniform (seed) < strandBias)
			ev[cnt].strand = (RandomUniform (seed) < 0.5) ? 0 : 1;
		cnt++;
		}
	if (cnt > 1)
		qsort (ev, cnt, sizeof(AmpEvent), CmpEvent);
	*nev = cnt;
	return ev;
	}

static void SeqErrors (char *s, int len, long int *seed)
	{
	int		k, ne;

	if (sequencingError <= 0.0)
		return;
	ne = RandomPoisson (sequencingError * len, seed);
	for (k=0; k<ne; k++)
		{
		int		p = RandomUniformTo (len, seed), a = (int) (strchr (BASES, s[p]) - BASES), l;
		double	u = RandomUniform (seed), cum = 0.0;

		for (l=0; l<4; l++)
			if (l != a)
				{
				cum += ampSubRow[a][l];
				if (u <= cum)
					break;
				}
		if (l > 3)
			l = (a + 1) & 3;
		s[p] = BASES[l];
		}
	}

/* Base qualities (--qual-sd, --qual-slope): the error probability of every base is eps * w(cycle) * m, where w rises along the read (log-slope
   --qual-slope between the first and the last cycle) and m is a lognormal multiplier of the base (log-sd --qual-sd); both have mean 1, so -E stays the
   mean error. The quality written is the one of that probability (calibrated: Q = -10 log10 p). The qualities are in the orientation of the BAM. */
static double	qualMult[1024], qualW[512];
static int		qualReady = 0;

static void QualSetup (int len, long int *seed)
	{
	int		k, i;
	double	sum = 0.0, sd = qualSd;

	if (qualReady)
		return;
	for (k=0; k<1024; k++)
		{
		qualMult[k] = (sd > 0.0) ? exp (sd * Normal (seed)) : 1.0;
		sum += qualMult[k];
		}
	for (k=0; k<1024; k++)
		qualMult[k] *= 1024.0 / sum;
	sum = 0.0;
	for (i=0; i<len; i++)
		{
		qualW[i] = exp (qualSlope * ((len > 1 ? (double) i / (len - 1) : 0.0) - 0.5));
		sum += qualW[i];
		}
	for (i=0; i<len; i++)
		qualW[i] *= len / sum;
	qualReady = 1;
	}

static void SeqErrorsQ (char *s, int len, int reverse, unsigned char *qv, double eps, long int *seed)
	{
	int		i;

	for (i=0; i<len; i++)
		{
		double	e = eps * qualW[reverse ? len - 1 - i : i] * qualMult[RandomUniformTo (1024, seed)];
		int		qq;

		if (e > 0.75) e = 0.75;
		qq = (int) (-10.0 * log10 (e) + 0.5);
		qv[i] = (unsigned char) (qq > 41 ? 41 : (qq < 2 ? 2 : qq));
		if (sequencingError > 0.0 && RandomUniform (seed) < e)
			{
			int		a = (int) (strchr (BASES, s[i]) - BASES), l;
			double	u = RandomUniform (seed), cum = 0.0;

			for (l=0; l<4; l++)
				if (l != a)
					{
					cum += ampSubRow[a][l];
					if (u <= cum)
						break;
					}
			if (l > 3)
				l = (a + 1) & 3;
			s[i] = BASES[l];
			}
		}
	}

static void RevComp (const char *s, char *o, int len)
	{
	int	i;

	for (i=0; i<len; i++)
		{
		char c = s[len - 1 - i];

		o[i] = (c == 'A') ? 'T' : (c == 'C') ? 'G' : (c == 'G') ? 'C' : 'A';
		}
	}

typedef struct { Bgzf r1, r2; FILE *f1, *f2; } FastqOut;

/* the output files of one cell (or of the healthy genome); in a chromosome in blocks they stay open until the last block */
typedef struct { BamOut bo; FastqOut fq; long molecule; } ReadOut;

/* position (0-based) in the BAM of the site j of the current block */
static inline long CP (long j)
	{
	return blocksOn ? BlockSitePos ((int) j) - 1 : j;
	}

static void FastqPut (Bgzf *z, const char *name, int mate, const char *seq, int len, int q, const unsigned char *qv, int reverse)
	{
	char	line[512], qs[1024];
	int		n;

	n = snprintf (line, sizeof(line), "@%s/%d\n", name, mate);
	BgzfWrite (z, line, n);
	BgzfWrite (z, seq, len);
	BgzfWrite (z, "\n+\n", 3);
	if (qv == NULL)
		memset (qs, q + 33, len);
	else
		{
		int i;

		for (i=0; i<len; i++)
			qs[i] = (char) (qv[reverse ? len - 1 - i : i] + 33);
		}
	BgzfWrite (z, qs, len);
	BgzfWrite (z, "\n", 1);
	}

/* the two templates of a fragment (first and last L bases): 1 when a site is lost (deletion or ADO), else 0 */
static int LoadTemplates (int hap, int cellIdx, long s, long ins, int L, char *ta, char *tb, int copyRank)
	{
	int		j, a, b;

	for (j=0; j<L; j++)
		{
		a = Allele (hap, cellIdx, s + j);
		b = Allele (hap, cellIdx, s + ins - L + j);
		if (copyRank > 0)		/* the fragment comes from an extra copy of the haplotype (--seg-gain-rate) where there is one */
			{
			Geno	ex[MAX_COPIES];

			if (ExtraCopies (cellIdx, hap, s + j, ex) >= copyRank)
				a = ex[copyRank-1];
			if (ExtraCopies (cellIdx, hap, s + ins - L + j, ex) >= copyRank)
				b = ex[copyRank-1];
			}
		if (a < 0 || a > 3 || b < 0 || b > 3)
			return 1;
		ta[j] = BASES[a];
		tb[j] = BASES[b];
		}
	return 0;
	}

/* Errors of the amplification carried by a molecule (shared by its PCR duplicates): each event within reach of the
   reads of the fragment [s, s+ins) is inherited with probability f, on the strand it happened on */
static void ApplyAmpErrors (const AmpEvent *ev, long nev, long s, long ins, int L, char *ta, char *tb, long int *seed)
	{
	long	lo = 0, hi = nev, mid, k, e;
	long	errSite[16];
	int		errAlt[16], errStr[16], nerr = 0;

	while (lo < hi)
		{
		mid = (lo + hi) / 2;
		if (ev[mid].site < s)
			lo = mid + 1;
		else
			hi = mid;
		}
	for (k=lo; k<nev && ev[k].site < s + ins && nerr < 16; k++)
		{
		e = ev[k].site;
		if (e < s + L || e >= s + ins - L)
			if (RandomUniform (seed) < ev[k].f)
				{
				errSite[nerr] = e; errAlt[nerr] = ev[k].alt; errStr[nerr] = ev[k].strand; nerr++;
				}
		}
	for (k=0; k<nerr; k++)
		{
		if (errSite[k] < s + L && errStr[k] != 1)
			ta[errSite[k] - s] = BASES[errAlt[k]];
		if (errSite[k] >= s + ins - L && errStr[k] != 0)
			tb[errSite[k] - (s + ins - L)] = BASES[errAlt[k]];
		}
	}

/* number of copies of a molecule after PCR duplication (geometric, mean 1 + dupMean) */
static int DuplicateCopies (long int *seed)
	{
	int		copies = 1;
	double	qd = dupMean / (1.0 + dupMean);

	if (dupMean > 0.0)
		while (RandomUniform (seed) < qd && copies < 50)
			copies++;
	return copies;
	}

/* writes one read pair (copy c of a molecule) to the BAM records and/or the FASTQ files; ra and rb are the first and last reads of the fragment */
static void EmitPair (RecSet *rs, FastqOut *fq, const char *name, int fwdFrag, int c, long s, long ins, const char *ra, const char *rb, int q,
					  const unsigned char *qa, const unsigned char *qb, const char *rg, int hap, long molecule)
	{
	int		L = readLen, fa, fb;
	char	rc[512];

	if (fwdFrag)
		{ fa = 99; fb = 147; }
	else
		{ fa = 163; fb = 83; }
	if (c > 0)
		{ fa |= 1024; fb |= 1024; }
	if (doBam)
		{
		AddRecord (rs, name, fa, CP (s), CP (s) + ins - L, ins, ra, L, q, qa, rg, hap + 1, (int) molecule);
		AddRecord (rs, name, fb, CP (s) + ins - L, CP (s), -ins, rb, L, q, qb, rg, hap + 1, (int) molecule);
		}
	if (doFastq)
		{
		RevComp (rb, rc, L);
		if (fwdFrag)
			{
			FastqPut (&fq->r1, name, 1, ra, L, q, qa, 0);
			FastqPut (&fq->r2, name, 2, rc, L, q, qb, 1);
			}
		else
			{
			FastqPut (&fq->r1, name, 1, rc, L, q, qb, 1);
			FastqPut (&fq->r2, name, 2, ra, L, q, qa, 0);
			}
		}
	}

/* Reads of one cell (cellIdx >= 0) or of the healthy genome (cellIdx = -1). Returns the number of read pairs. */
static long SimulateGenome (int cellIdx0, double depth, int amplified, const char *rg, long int *seed, ReadOut *ro, int donor)
	{
	RecSet	rs;
	FastqOut *fq = &ro->fq;
	AmpEvent *ev[2] = {NULL, NULL};
	long	nev[2] = {0, 0}, pairs = 0, molecule = ro->molecule, pos;
	int		hap, q, L = readLen, i;
	char	*ta, *tb, *ra, *rb, name[96];
	unsigned char qa[512], qb[512];
	int		qualOn = (qualSd > 0.0 || qualSlope != 0.0);
	double	qerr = (sequencingError > 0.0) ? sequencingError : 1e-4;
	const long W = 256;
	int		srcs[2] = {cellIdx0, donor}, nsrc = (donor >= 0) ? 2 : 1, si;

	memset (&rs, 0, sizeof(rs));
	q = (int) (-10.0 * log10 (qerr) + 0.5);
	if (q > 41) q = 41;
	if (q < 2) q = 2;
	if (qualOn)
		QualSetup (L, seed);
	ta = (char*) malloc (4 * (size_t) (L + 2));
	tb = ta + L + 2; ra = tb + L + 2; rb = ra + L + 2;
	/* a doublet (a cell captured together with its donor) has the reads of both genomes in the same files, in one sorted set of records */
	for (si=0; si<nsrc; si++)
	{
	int		cellIdx = srcs[si];
	double	cm = (cellIdx >= 0 && cellCovMult != NULL) ? cellCovMult[cellIdx] : 1.0;

	nev[0] = nev[1] = 0;
	for (hap=0; hap<2; hap++)
		{
		if (amplified)
			ev[hap] = BuildAmpEvents (cellIdx, hap, seed, &nev[hap]);
		}

	for (hap=0; hap<2; hap++)
		{
		for (pos = 0; pos < numSites; pos += W)
			{
			double	g = 1.0, lambda;
			int		nfrag, f;

			if (amplified)
				{
				if (gainOn == YES && cellIdx >= 0)
					g = GainAt (cellIdx, hap, pos + W/2);
				if (copyGainOn == YES && cellIdx >= 0)
					{
					Geno	ex[MAX_COPIES];

					g *= 1.0 + ExtraCopies (cellIdx, hap, pos + W/2, ex);	/* the extra copies (--seg-gain-rate) add fragments */
					}
				}
			lambda = depth * cm * 0.5 * g / (2.0 * L * (1.0 + dupMean)) * (double) W;
			nfrag = RandomPoisson (lambda, seed);
			for (f=0; f<nfrag; f++)
				{
				long	s = pos + (long) (RandomUniform (seed) * W), ins;
				int		copies, c, fwdFrag, copyRank = 0;

				ins = (long) (insertMean + insertSd * Normal (seed) + 0.5);
				if (ins < L)
					ins = L;
				if (s + ins > numSites)
					continue;
				if (posMap != NULL && CP (s + ins - 1) - CP (s) != ins - 1)
					continue;	/* the fragment would span masked sites (N) */
				if (copyGainOn == YES && cellIdx >= 0)
					{
					Geno	ex[MAX_COPIES];

					copyRank = RandomUniformTo (1 + ExtraCopies (cellIdx, hap, s, ex), seed);	/* the fragment comes from one of the copies */
					}
				if (LoadTemplates (hap, cellIdx, s, ins, L, ta, tb, copyRank))
					continue;
				if (nev[hap] > 0)
					ApplyAmpErrors (ev[hap], nev[hap], s, ins, L, ta, tb, seed);
				copies = DuplicateCopies (seed);
				fwdFrag = (RandomUniform (seed) < 0.5);
				molecule++;
				for (c=0; c<copies; c++)
					{
					memcpy (ra, ta, L); memcpy (rb, tb, L);
					if (qualOn)
						{
						SeqErrorsQ (ra, L, 0, qa, qerr, seed);
						SeqErrorsQ (rb, L, 1, qb, qerr, seed);
						}
					else
						{
						SeqErrors (ra, L, seed);
						SeqErrors (rb, L, seed);
						}
					snprintf (name, sizeof(name), "%s:%ld", rg, molecule);
					EmitPair (&rs, fq, name, fwdFrag, c, s, ins, ra, rb, q, qualOn ? qa : NULL, qualOn ? qb : NULL, rg, hap, molecule);
					pairs++;
					}
				}
			}
		}
	for (i=0; i<2; i++)
		{
		free (ev[i]);
		ev[i] = NULL;
		}
	}
	if (doBam)
		BamAppend (&ro->bo, &rs);
	ro->molecule = molecule;
	free (ta);
	free (rs.buf);
	free (rs.idx);
	return pairs;
	}

/********************************** GenerateReadLevel ***********************************/
/* Output slots: one per cell and one for the healthy genome (--bulk-depth). Without --blocks they are opened and closed here;
   with --blocks they stay open from the first block to FinishReadLevel. */
static ReadOut	*routs = NULL;
static int		nrouts = 0;
static FILE		*refFp = NULL;
static long		refPos = 0;			/* bases of the chromosome already written to ref.fa */
static int		*readDonor = NULL;	/* donor cell of every doublet (-1: not a doublet); chosen in the first block of a chromosome */
static char		readDir[3*MAX_NAME];

static void RefPut (int c)
	{
	fputc (c, refFp);
	refPos++;
	if (refPos % 60 == 0)
		fputc ('\n', refFp);
	}

static void ReadOutOpen (const char *name, ReadOut *ro)
	{
	char	path[4*MAX_NAME];

	memset (ro, 0, sizeof(*ro));
	if (doBam)
		{
		snprintf (path, sizeof(path), "%s/%s.bam", readDir, name);
		BamOpen (&ro->bo, path, name);
		}
	if (doFastq)
		{
		snprintf (path, sizeof(path), "%s/%s.R1.fastq.gz", readDir, name);
		ro->fq.f1 = fopen (path, "wb");
		snprintf (path, sizeof(path), "%s/%s.R2.fastq.gz", readDir, name);
		ro->fq.f2 = fopen (path, "wb");
		if (!ro->fq.f1 || !ro->fq.f2)
			{
			fprintf (stderr, "\nERROR: cannot write the FASTQ files %s/%s\n", readDir, name);
			exit (-1);
			}
		ro->fq.r1.fp = ro->fq.f1;
		ro->fq.r2.fp = ro->fq.f2;
		CrcInit ();
		}
	}

static void ReadOutClose (ReadOut *ro)
	{
	if (doBam)
		BamClose (&ro->bo);
	if (doFastq)
		{
		BgzfClose (&ro->fq.r1); BgzfClose (&ro->fq.r2);
		fclose (ro->fq.f1); fclose (ro->fq.f2);
		}
	}

void FinishReadLevel (void)
	{
	char	path[4*MAX_NAME];
	FILE	*fp;
	int		i;
	long	len = ReadRefLen ();

	if (routs == NULL)
		return;
	for (i=0; i<nrouts; i++)
		ReadOutClose (&routs[i]);
	free (routs);
	routs = NULL;
	free (readDonor);
	readDonor = NULL;
	while (refPos < len)
		RefPut ('N');
	if (refPos % 60)
		fputc ('\n', refFp);
	fclose (refFp);
	refFp = NULL;
	/* index of the reference (samtools faidx format) */
	snprintf (path, sizeof(path), "%s/ref.fa.fai", readDir);
	if ((fp = fopen (path, "w")) != NULL)
		{
		fprintf (fp, "%s\t%ld\t%d\t60\t61\n", ReadContig (), len, (int) strlen (ReadContig ()) + 2);
		fclose (fp);
		}
	}

void GenerateReadLevel (long int *seed, int dataSetNum)
	{
	char	path[4*MAX_NAME], rg[64];
	long	j, total = 0;
	int		i;

	if (routs == NULL)
		{
		snprintf (readDir, sizeof(readDir), "%s/reads_%d", resultsDir, blocksOn ? 1 : dataSetNum + 1);
		MAKE_DIR (readDir);
		snprintf (path, sizeof(path), "%s/ref.fa", readDir);
		if ((refFp = fopen (path, "w")) == NULL)
			{
			fprintf (stderr, "\nERROR: cannot write %s\n", path);
			exit (-1);
			}
		fprintf (refFp, ">%s\n", ReadContig ());
		refPos = 0;
		nrouts = numCells + (bulkDepth > 0.0 ? 1 : 0);
		routs = (ReadOut*) calloc (nrouts, sizeof(ReadOut));
		for (i=0; i<nrouts; i++)
			{
			if (i < numCells)
				snprintf (rg, sizeof(rg), "cell%04d", i + 1);
			else
				snprintf (rg, sizeof(rg), "normal");
			ReadOutOpen (rg, &routs[i]);
			}
		}

	/* doublets: the cells and donors were chosen (SetupDoublets, with the correlations of SetupCellCoverage); their reads are those of the cell and of its donor */
	if (meanDoubletRate > 0.0 && readDonor == NULL && doubletDonor != NULL)
		{
		int		nd = 0;

		readDonor = (int*) malloc (numCells * sizeof(int));
		for (i=0; i<numCells; i++)
			{
			readDonor[i] = (cell[i].hasDoublet == YES) ? doubletDonor[i] : -1;
			if (readDonor[i] == i && numCells > 1)		/* a cell cannot be its own partner */
				{
				int d = RandomUniformTo (numCells - 1, seed);

				readDonor[i] = (d >= i) ? d + 1 : d;
				}
			nd += (readDonor[i] >= 0);
			free (cell[i].doublet);
			cell[i].doublet = NULL;
			cell[i].hasDoublet = NO;
			}
		free (doubletDonor);
		doubletDonor = NULL;
		if (nd > 0)
			{
			FILE *fd;

			snprintf (path, sizeof(path), "%s/doublets.tsv", readDir);
			if ((fd = fopen (path, "w")) != NULL)
				{
				fprintf (fd, "cell\tdonor\n");
				for (i=0; i<numCells; i++)
					if (readDonor[i] >= 0)
						{
						if (doUserTree == NO)
							fprintf (fd, "%s%04d\t%s%04d\n", inCellName, i + 1, inCellName, readDonor[i] + 1);
						else
							fprintf (fd, "%s\t%s\n", cellNames[i], cellNames[readDonor[i]]);
						}
				fclose (fd);
				}
			}
		}

	/* the reference is the healthy root genome (the FASTA with --fasta); sites that are not in the simulation are N */
	for (j=0; j<numSites; j++)
		{
		long p0 = CP (j);
		int a = data[MATERNAL][OUTGROUP_ROOT][j];
		char c = (a >= 0 && a < 4) ? BASES[a] : 'N';

		if (refUserGenome != NULL)
			c = refUserGenome[j];
		while (refPos < p0)
			RefPut ('N');
		RefPut (c);
		}

	/* every cell (and the healthy bulk) has its own random stream, taken from the main stream, so the cells are simulated in parallel and
	   the files are the same for any number of threads */
	{
	uint64_t	base;
	long		s1, s2;

	if (qualSd > 0.0 || qualSlope != 0.0)
		QualSetup (readLen, seed);		/* the multipliers of the bases are shared by all the cells: drawn once, before the threads start */
	RandomUniform (seed);
	s1 = *seed;
	RandomUniform (seed);
	s2 = *seed;
	base = ((uint64_t) s1 << 32) ^ (uint64_t) s2;

	#pragma omp parallel for schedule(dynamic,1) reduction(+:total)
	for (i=0; i<nrouts; i++)
		{
		char	rgi[64];

		if (i < numCells)
			snprintf (rgi, sizeof(rgi), "cell%04d", i + 1);
		else
			snprintf (rgi, sizeof(rgi), "normal");
		StartStream (base, (uint64_t) i + 1);
		if (i < numCells)
			total += SimulateGenome (i, coverage, 1, rgi, seed, &routs[i], readDonor != NULL ? readDonor[i] : -1);
		else
			SimulateGenome (-1, bulkDepth, 0, rgi, seed, &routs[i], -1);
		streamActive = 0;
		}
	}
	if (noisy > 0)
		fprintf (stderr, "\rreads: %ld pairs in %d cells -> %s\n", total, numCells, readDir);
	if (blocksOn == 0)
		{
		FinishReadLevel ();
		free (readDonor);
		readDonor = NULL;
		}
	}


/********************************** WriteTruthVcf ***********************************/
/* True genotypes of the cells (before dropout and genotyping errors) at the sites where any cell or the germline differs from the reference
   (ref.fa, the maternal healthy genome). Called before the observation layer is active. GT is maternal|paternal; a deleted haplotype is '.'.
   INFO: NC = cells with a non-reference allele, GERMLINE = the healthy genome is heterozygous or differs from the reference, TYPE = SNV. */
void WriteTruthVcf (int dataSetNum)
	{
	char	dir[3*MAX_NAME], path[4*MAX_NAME];
	FILE	*fp;
	long	j;
	int		i, h;

	static int started = 0;		/* with --blocks the file is written block after block */

	snprintf (dir, sizeof(dir), "%s/reads_%d", resultsDir, blocksOn ? 1 : dataSetNum + 1);
	MAKE_DIR (dir);
	snprintf (path, sizeof(path), "%s/truth.vcf", dir);
	if ((fp = fopen (path, (blocksOn && started) ? "a" : "w")) == NULL)
		{
		fprintf (stderr, "\nERROR: cannot write %s\n", path);
		exit (-1);
		}
	if (!(blocksOn && started))
		{
	fprintf (fp, "##fileformat=VCFv4.2\n##source=cellcoal truth\n##contig=<ID=%s,length=%ld>\n", ReadContig (), ReadRefLen ());
	fprintf (fp, "##INFO=<ID=NC,Number=1,Type=Integer,Description=\"Cells with a non-reference allele\">\n");
	fprintf (fp, "##INFO=<ID=GERMLINE,Number=0,Type=Flag,Description=\"The healthy genome differs from the reference or is heterozygous\">\n");
	fprintf (fp, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"True genotype (maternal|paternal)\">\n");
	fprintf (fp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	for (i=0; i<numCells; i++)
		fprintf (fp, "\tcell%04d", i + 1);
	fprintf (fp, "\n");
		}
	started = 1;
	for (j=0; j<numSites; j++)
		{
		int		mr = data[MATERNAL][OUTGROUP_ROOT][j], ref = mr, pr = data[PATERNAL][OUTGROUP_ROOT][j], var, nc = 0, nalt = 0, altBase[4] = {-1,-1,-1,-1}, g[2];
		char	gt[8];

		if (refUserGenome != NULL)	/* the reference is the FASTA; the healthy genome may carry germline alleles */
			ref = (int) (strchr (BASES, refUserGenome[j]) - BASES);
		var = (pr != ref || mr != ref);
		if (ref < 0 || ref > 3)
			continue;
		if (sparseMode == YES && sparseCols[j] == NULL && !var)
			continue;
		for (i=0; i<numCells; i++)
			{
			int any = 0;

			for (h=0; h<2; h++)
				{
				g[h] = GetG (h, i, (int) j);
				if (g[h] >= 0 && g[h] < 4 && g[h] != ref)
					{
					int k;

					any = 1;
					for (k=0; k<nalt && altBase[k] != g[h]; k++)
						;
					if (k == nalt && nalt < 4)
						altBase[nalt++] = g[h];
					}
				}
			nc += any;
			}
		for (h=0; h<2; h++)
			{
			int k, gh = (h == 0) ? mr : pr;

			if (gh >= 0 && gh < 4 && gh != ref)
				{
				for (k=0; k<nalt && altBase[k] != gh; k++)
					;
				if (k == nalt && nalt < 4)
					altBase[nalt++] = gh;
				}
			}
		if (nalt == 0)
			continue;
		fprintf (fp, "%s\t%ld\t.\t%c\t", ReadContig (), CP (j) + 1, BASES[ref]);
		{
		int k;

		for (k=0; k<nalt; k++)
			fprintf (fp, "%s%c", k ? "," : "", BASES[altBase[k]]);
		}
		fprintf (fp, "\t.\tPASS\tNC=%d%s\tGT", nc, var ? ";GERMLINE" : "");
		for (i=0; i<numCells; i++)
			{
			for (h=0; h<2; h++)
				{
				int a = GetG (h, i, (int) j), k, code = -1;

				if (a == ref)
					code = 0;
				else
					for (k=0; k<nalt; k++)
						if (altBase[k] == a)
							code = k + 1;
				if (code < 0)
					strcpy (gt + 2*h, ".");
				else
					sprintf (gt + 2*h, "%d", code);
				if (h == 0)
					gt[1] = '|';
				}
			gt[3] = '\0';
			fprintf (fp, "\t%s", gt);
			}
		fprintf (fp, "\n");
		}
	fclose (fp);
	}
