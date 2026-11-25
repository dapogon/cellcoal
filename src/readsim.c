/* Read-level simulator (--bam, --fastq): paired-end reads from the amplified genome of every cell.

   For every cellIdx and haplotype it draws fragments along the genome with a density given by the cellIdx coverage (-C, -J), the amplification
   gain (--gain-sd) and the allelic dropout (the lost segments are not covered).
   The fragments are molecules: errors of the amplification are carried by the molecule and are copied to its PCR duplicates (--dup-mean);
   sequencing errors (-E) are drawn independently in every read. Output: one BAM (+ BAI) and/or FASTQ.gz pair per cellIdx, the reference (ref.fa)
   and, with --bulk-depth, the BAM of the healthy genome. The BAM flag 0x400 marks the PCR duplicates, the tag MI is the molecule and HP the haplotype. */

#include <limits.h>
#include "cellcoal.h"
#include "bgzf.h"
#include "bam.h"
#include <sys/stat.h>
#include <math.h>
#include <string.h>


typedef struct { long site; int alt; double f; int strand; } AmpEvent;	/* strand: -1 both, 0 forward, 1 reverse */



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

/* Sequencing errors that depend on the sequence (--seq-context FILE). A line of the file is a context of three bases in the direction of the
   sequencing (the two bases read before and the base itself: GGC is a C that follows GG) and a multiplier of the error probability of that
   base; N stands for any base. The lines are applied in order (a later line overrides the earlier ones on the contexts that both match);
   the contexts that no line matches have multiplier 1. The multipliers are divided by their mean over the trinucleotides of the reference,
   so -E stays the mean error of the reads. */
static double	ctxMult[64];
static int		ctxOn = 0;

static int BaseCode (int c)
	{
	return (c == 'A') ? 0 : (c == 'C') ? 1 : (c == 'G') ? 2 : (c == 'T') ? 3 : -1;
	}

static void LoadSeqContext (void)
	{
	FILE	*fp = fopen (seqContextFile, "r");
	char	line[512], pat[16];
	int		ln = 0, k, i;
	double	m, freq3[64], sum, wsum, mean;
	long	tot = 0;

	if (!fp)
		{
		fprintf (stderr, "\nERROR: cannot open %s\n", seqContextFile);
		exit (-1);
		}
	for (k=0; k<64; k++)
		ctxMult[k] = 1.0;
	while (fgets (line, sizeof(line), fp))
		{
		int ok = 1, a[3], x, y, z;

		ln++;
		if (line[0] == '#' || line[0] == '\n' || line[0] == '\r' || line[0] == '\0')
			continue;
		if (sscanf (line, "%15s %lf", pat, &m) != 2 || strlen (pat) != 3 || m < 0.0)
			{
			fprintf (stderr, "\nERROR: %s, line %d: expected \"CONTEXT MULTIPLIER\" (three bases, A C G T or N, and a number >= 0)\n", seqContextFile, ln);
			exit (-1);
			}
		for (i=0; i<3; i++)
			{
			int c = toupper ((unsigned char) pat[i]);

			a[i] = (c == 'N') ? -1 : BaseCode (c);
			if (c != 'N' && a[i] < 0)
				ok = 0;
			}
		if (!ok)
			{
			fprintf (stderr, "\nERROR: %s, line %d: the context has to be made of A, C, G, T or N\n", seqContextFile, ln);
			exit (-1);
			}
		for (x=0; x<4; x++)
			for (y=0; y<4; y++)
				for (z=0; z<4; z++)
					if ((a[0] < 0 || a[0] == x) && (a[1] < 0 || a[1] == y) && (a[2] < 0 || a[2] == z))
						ctxMult[x*16 + y*4 + z] = m;
		}
	fclose (fp);
	/* trinucleotides of the reference: the mean multiplier has to be 1 */
	for (k=0; k<64; k++)
		freq3[k] = 0.0;
	for (i=0; i+2<numSites && tot < 2000000; i++)
		{
		int x = data[MATERNAL][OUTGROUP_ROOT][i], y = data[MATERNAL][OUTGROUP_ROOT][i+1], z = data[MATERNAL][OUTGROUP_ROOT][i+2];

		if (refUserGenome != NULL)
			{
			x = BaseCode (refUserGenome[i]); y = BaseCode (refUserGenome[i+1]); z = BaseCode (refUserGenome[i+2]);
			}
		if (x < 0 || x > 3 || y < 0 || y > 3 || z < 0 || z > 3)
			continue;
		freq3[x*16 + y*4 + z] += 1.0;
		tot++;
		}
	sum = wsum = 0.0;
	for (k=0; k<64; k++)
		{
		sum += freq3[k];
		wsum += freq3[k] * ctxMult[k];
		}
	mean = (sum > 0.0) ? wsum / sum : 1.0;
	if (mean <= 0.0)
		{
		fprintf (stderr, "\nERROR: %s: the sequencing error of the contexts of the reference is 0 everywhere\n", seqContextFile);
		exit (-1);
		}
	for (k=0; k<64; k++)
		ctxMult[k] /= mean;
	ctxOn = 1;
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
	char	orig[512];

	memcpy (orig, s, len);		/* the context is that of the template, not of the bases that already have an error */
	for (i=0; i<len; i++)
		{
		double	e = eps * qualW[reverse ? len - 1 - i : i] * qualMult[RandomUniformTo (1024, seed)];
		int		qq;

		if (ctxOn)
			{
			int x, y, z;

			if (!reverse)			/* read in the direction of the genome: the bases before are the previous ones */
				{
				x = (i >= 2) ? BaseCode (orig[i-2]) : -1;
				y = (i >= 1) ? BaseCode (orig[i-1]) : -1;
				z = BaseCode (orig[i]);
				}
			else					/* the second read is the reverse complement: it reads the next bases first */
				{
				x = (i + 2 < len) ? 3 - BaseCode (orig[i+2]) : -1;
				y = (i + 1 < len) ? 3 - BaseCode (orig[i+1]) : -1;
				z = 3 - BaseCode (orig[i]);
				if (x > 3) x = -1;
				if (y > 3) y = -1;
				if (z > 3) z = -1;
				}
			if (x >= 0 && y >= 0 && z >= 0)
				e *= ctxMult[x*16 + y*4 + z];
			}

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
static void ApplyAmpErrors (const AmpEvent *ev, long nev, long s, long ins, int L, char *ta, char *tb, long int *seed, const long *ma, const long *mb)
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
	if (ma != NULL)		/* the reads of a fragment with indels: the bases of the reads are found by their site */
		{
		for (k=0; k<nerr; k++)
			{
			int i;

			for (i=0; i<L; i++)
				{
				if (ma[i] == errSite[k] && errStr[k] != 1)
					ta[i] = BASES[errAlt[k]];
				if (mb[i] == errSite[k] && errStr[k] != 0)
					tb[i] = BASES[errAlt[k]];
				}
			}
		return;
		}
	for (k=0; k<nerr; k++)
		{
		if (errSite[k] < s + L && errStr[k] != 1)
			ta[errSite[k] - s] = BASES[errAlt[k]];
		if (errSite[k] >= s + ins - L && errStr[k] != 0)
			tb[errSite[k] - (s + ins - L)] = BASES[errAlt[k]];
		}
	}

/* Damage of the library (--damage "G>T:0.002,C>T:0.001"): a base of one strand of the molecule changes before the PCR, so the PCR duplicates
   inherit it, and it is seen only in the read that sequences that strand (the orientation bias of the oxidation, G>T, and of the deamination,
   C>T). Every item is FROM>TO:RATE, with RATE the probability that a FROM base of a strand is damaged. A damaged G of the top strand is a T of
   the first read; a damaged G of the bottom strand is a C of the genome that is seen as an A in the second read. */
typedef struct { char from, to; double rate; } DamageItem;
static DamageItem	dmg[16];
static int			ndmg = 0;

static char Comp (char c)
	{
	return (c == 'A') ? 'T' : (c == 'C') ? 'G' : (c == 'G') ? 'C' : (c == 'T') ? 'A' : c;
	}

static void LoadDamage (void)
	{
	char	*copy = strdup (damageSpec), *tok, *save = NULL;

	for (tok = strtok_r (copy, ",", &save); tok; tok = strtok_r (NULL, ",", &save))
		{
		char	f, t;
		double	r;
		int		n = 0;

		if (sscanf (tok, " %c>%c:%lf%n", &f, &t, &r, &n) != 3 || tok[n] != '\0')
			{
			fprintf (stderr, "\nERROR: --damage: \"%s\" is not FROM>TO:RATE (for example G>T:0.002)\n", tok);
			exit (-1);
			}
		f = toupper ((unsigned char) f); t = toupper ((unsigned char) t);
		if (!strchr ("ACGT", f) || !strchr ("ACGT", t) || f == t || r < 0.0 || r > 1.0)
			{
			fprintf (stderr, "\nERROR: --damage: \"%s\": the bases are A, C, G or T (different ones) and the rate is between 0 and 1\n", tok);
			exit (-1);
			}
		if (ndmg == 16)
			{
			fprintf (stderr, "\nERROR: --damage: at most 16 items\n");
			exit (-1);
			}
		dmg[ndmg].from = f; dmg[ndmg].to = t; dmg[ndmg].rate = r;
		ndmg++;
		}
	free (copy);
	if (ndmg == 0)
		{
		fprintf (stderr, "\nERROR: --damage: no items (for example G>T:0.002)\n");
		exit (-1);
		}
	}

static void ApplyDamage (int L, char *ta, char *tb, long int *seed)
	{
	int		k, i, n, p;

	for (k=0; k<ndmg; k++)
		{
		n = RandomPoisson (dmg[k].rate * L, seed);		/* every position is a candidate with probability rate, and changes if it is a FROM base */
		for (i=0; i<n; i++)
			{
			p = RandomUniformTo (L, seed);
			if (ta[p] == dmg[k].from)
				ta[p] = dmg[k].to;
			}
		n = RandomPoisson (dmg[k].rate * L, seed);
		for (i=0; i<n; i++)
			{
			p = RandomUniformTo (L, seed);
			if (tb[p] == Comp (dmg[k].from))
				tb[p] = Comp (dmg[k].to);
			}
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

/* The reads of a fragment that carries indels: the haplotype of the fragment (the sites s .. s+ins-1 with its insertions and deletions) is built, the first
   and the last L bases are the reads, and every base remembers its site (-1 for an inserted base), which gives the CIGAR and the positions. */
#define EXT 8		/* extra bases of the templates: a deletion error of the sequencer takes bases that follow the read */
typedef struct
	{
	long		ma[520], mb[520];			/* site of every base of the templates of the two reads (-1: inserted base); the template of the second read has extB bases before it */
	int			la, lb, extB;				/* length of the templates, and bases of the second template before its read */
	long		mra[512], mrb[512];			/* site of every base of the reads */
	uint32_t	ca[64], cb[64];				/* CIGAR of the reads */
	int			nca, ncb;
	long		posA, endA, posB, endB;		/* first and last site aligned by every read */
	} Aln;

static int SiteBase (int hap, int cellIdx, long j, int copyRank)
	{
	int a = Allele (hap, cellIdx, j);

	if (copyRank > 0)
		{
		Geno	ex[MAX_COPIES];

		if (ExtraCopies (cellIdx, hap, j, ex) >= copyRank)
			a = ex[copyRank-1];
		}
	return (a < 0 || a > 3) ? -1 : a;
	}

static int MakeCigar (const long *pm, int L, uint32_t *cig, long *first, long *last)
	{
	int		n = 0, i, fm = 0, lm = L - 1;
	long	prev;
	uint32_t op = 99, len = 0;

#define PUSHOP(O,N) do { if ((uint32_t)(O) == op) len += (N); else { if (len > 0) { if (n >= 63) return -1; cig[n++] = (len << 4) | op; } op = (O); len = (N); } } while (0)
	while (fm < L && pm[fm] < 0)
		fm++;
	while (lm >= 0 && pm[lm] < 0)
		lm--;
	if (fm >= L || lm < 0)
		return -1;
	if (fm > 0)
		PUSHOP (4, fm);
	prev = -1;
	for (i=fm; i<=lm; i++)
		{
		if (pm[i] < 0)
			PUSHOP (1, 1);
		else
			{
			if (prev >= 0 && pm[i] - prev > 1)
				PUSHOP (2, (uint32_t) (pm[i] - prev - 1));
			PUSHOP (0, 1);
			prev = pm[i];
			}
		}
	if (lm < L - 1)
		PUSHOP (4, L - 1 - lm);
	if (len > 0)
		{
		if (n >= 64)
			return -1;
		cig[n++] = (len << 4) | op;
		}
#undef PUSHOP
	*first = pm[fm];
	*last = pm[lm];
	return n;
	}

/* The CIGARs and positions of the reads from the sites of their bases; 0 when they can be aligned */
static int FinishAln (Aln *aln, int L)
	{
	aln->nca = MakeCigar (aln->mra, L, aln->ca, &aln->posA, &aln->endA);
	aln->ncb = MakeCigar (aln->mrb, L, aln->cb, &aln->posB, &aln->endB);
	return (aln->nca < 0 || aln->ncb < 0) ? -1 : 0;
	}

/* Templates of a fragment from the haplotype of the cell with its indels (the sites s .. s+ins-1 with the insertions and deletions that fall in it);
   with force (sequencing indel errors) the template is built even if no indel falls in the fragment. The template of the first read is the first
   L bases (and up to EXT more) and that of the second read the last L bases (with up to EXT before them). Returns 0 when there is nothing to do
   (no indel), 2 when the templates and their maps were built, and 3 when the fragment is lost (a site lost by deletion or dropout in the reads). */
static int IndelTemplates (const IndelEv *ie, long nie, int hap, int cellIdx, long s, long ins, int L, int copyRank, char *ta, char *tb, Aln *aln, int force)
	{
	long	lo = 0, hi = nie, mid, k, j, n = 0, cap;
	int		*hb;
	long	*hp;
	int		i, ok = 1, ext;

	while (lo < hi)
		{
		mid = (lo + hi) / 2;
		if (ie[mid].pos < s)
			lo = mid + 1;
		else
			hi = mid;
		}
	if (!force && (lo >= nie || ie[lo].pos >= s + ins - 1))
		return 0;
	cap = ins + 64 * (indelMax + 1) + 8;
	cap += cap & 1;						/* even, so that the long values that follow the ints are aligned */
	hb = (int*) Checked (malloc (cap * (sizeof(int) + sizeof(long))), "Could not allocate the reads\n");
	hp = (long*) (hb + cap);
	k = lo;
	for (j=s; j<s+ins; j++)
		{
		hb[n] = SiteBase (hap, cellIdx, j, copyRank);
		hp[n] = j;
		n++;
		while (k < nie && ie[k].pos == j)
			{
			const IndelEv *e = &ie[k++];

			if (e->type == 0)
				{
				for (i=0; i<e->len && n < cap - 1; i++)
					{
					hb[n] = (int) (strchr (BASES, e->seq[i]) - BASES);
					hp[n] = -1;
					n++;
					}
				}
			else
				j += (j + e->len < s + ins) ? e->len : (s + ins - 1 - j);
			}
		if (n >= cap - 1)
			break;
		}
	if (n < L)
		ok = 0;
	ext = (homoError > 0.0 && ok) ? (int) ((n - L < EXT) ? n - L : EXT) : 0;
	aln->la = L + ext;
	aln->lb = L + ext;
	aln->extB = ext;
	for (i=0; ok && i<L + ext; i++)
		{
		int ia = i, ib = (int) (n - L - ext + i);

		if (hb[ia] < 0 && i < L)
			ok = 0;
		if (hb[ib] < 0 && i >= ext)
			ok = 0;
		ta[i] = (hb[ia] >= 0) ? BASES[hb[ia]] : 'N';
		tb[i] = (hb[ib] >= 0) ? BASES[hb[ib]] : 'N';
		aln->ma[i] = hp[ia];
		aln->mb[i] = hp[ib];
		}
	free (hb);
	if (!ok)
		return 3;
	for (i=0; i<L; i++)
		{
		aln->mra[i] = aln->ma[i];
		aln->mrb[i] = aln->mb[ext + i];
		}
	return (FinishAln (aln, L) < 0) ? 3 : 2;
	}

/* Sequencing indel errors in homopolymers (--homopolymer-error P): a run of 3 or more equal bases of a read gets one base more or less with
   probability P times (length - 2). The deleted base is the first of the run and the inserted one follows it (left-aligned). The templates are
   slipped and the reads are their first L bases (first read) and last L bases (second read); the sites of the bases give the CIGAR. */
static int Slip (const char *t, const long *m, int n, char *out, long *mo, long int *seed)
	{
	int		i = 0, o = 0, j, c;

	while (i < n)
		{
		int len, kind = 0;

		for (j=i; j<n && t[j] == t[i]; j++)
			;
		len = j - i;
		if (len >= 3 && t[i] != 'N' && m[i] >= 0)
			{
			double p = homoError * (len - 2);

			if (p > 0.5)
				p = 0.5;
			if (RandomUniform (seed) < p)
				kind = (RandomUniform (seed) < 0.5) ? 1 : 2;
			}
		for (c=i; c<j; c++)
			{
			if (kind == 1 && c == i)
				continue;
			out[o] = t[c];
			mo[o++] = m[c];
			if (kind == 2 && c == i)
				{
				out[o] = t[c];
				mo[o++] = -1;
				}
			}
		i = j;
		}
	return o;
	}

/* the reads of a copy of the molecule from the templates (with or without sequencing indel errors); -1 when the reads cannot be aligned */
static int SlipReads (Aln *aln, const char *ta, const char *tb, int L, char *ra, char *rb, long int *seed)
	{
	char	oa[2*(512+EXT)+8], ob[2*(512+EXT)+8];
	long	ma[2*(512+EXT)+8], mb[2*(512+EXT)+8];
	int		i, na = 0, nb = 0, ok = 0;

	memcpy (ra, ta, L);
	memcpy (rb, tb + aln->extB, L);
	for (i=0; i<L; i++)
		{
		aln->mra[i] = aln->ma[i];
		aln->mrb[i] = aln->mb[aln->extB + i];
		}
	if (homoError > 0.0)
		{
		na = Slip (ta, aln->ma, aln->la, oa, ma, seed);
		nb = Slip (tb, aln->mb, aln->lb, ob, mb, seed);
		ok = (na >= L && nb >= L);
		for (i=0; ok && i<L; i++)
			if (oa[i] == 'N' || ob[nb - L + i] == 'N')
				ok = 0;
		}
	if (ok)
		{
		memcpy (ra, oa, L);
		memcpy (rb, ob + nb - L, L);
		for (i=0; i<L; i++)
			{
			aln->mra[i] = ma[i];
			aln->mrb[i] = mb[nb - L + i];
			}
		}
	return FinishAln (aln, L);
	}

/* writes one read pair (copy c of a molecule) to the BAM records and/or the FASTQ files; ra and rb are the first and last reads of the fragment */
static int		mappabilityOn = 0;		/* --mappability: the MAPQ of the reads follows the mappability of the genome */

static void EmitPair (RecSet *rs, FastqOut *fq, const char *name, int fwdFrag, int c, long s, long ins, const char *ra, const char *rb, int q,
					  const unsigned char *qa, const unsigned char *qb, const char *rg, int hap, long molecule, const Aln *aln)
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
		int	mqa = 60, mqb = 60;

		if (mappabilityOn)
			{
			if (aln != NULL)
				{
				mqa = MapQual (CP (aln->posA), CP (aln->endA));
				mqb = MapQual (CP (aln->posB), CP (aln->endB));
				}
			else
				{
				mqa = MapQual (CP (s), CP (s) + L - 1);
				mqb = MapQual (CP (s) + ins - L, CP (s) + ins - 1);
				}
			/* ambiguous pairs (MAPQ 0) are left out of the BAM with probability --map-drop; the decision depends on the molecule, so the
			   PCR duplicates of a read go together */
			if (mapDrop > 0.0 && (mqa == 0 || mqb == 0))
				{
				unsigned long h = (unsigned long) molecule * 2654435761UL + (unsigned long) hap * 40503UL + 12345UL;

				h ^= h >> 15; h *= 2246822519UL; h ^= h >> 13;
				if ((double) (h % 1000003UL) / 1000003.0 < mapDrop)
					goto fastq;
				}
			}
		if (aln != NULL)
			{
			long pa = CP (aln->posA), pb = CP (aln->posB), tl = CP (aln->endB) - pa + 1;

			AddRecord (rs, name, fa, pa, pb, tl, ra, L, q, qa, rg, hap + 1, (int) molecule, aln->ca, aln->nca, (int) (aln->endA - aln->posA + 1), mqa);
			AddRecord (rs, name, fb, pb, pa, -tl, rb, L, q, qb, rg, hap + 1, (int) molecule, aln->cb, aln->ncb, (int) (aln->endB - aln->posB + 1), mqb);
			}
		else
			{
			AddRecord (rs, name, fa, CP (s), CP (s) + ins - L, ins, ra, L, q, qa, rg, hap + 1, (int) molecule, NULL, 0, L, mqa);
			AddRecord (rs, name, fb, CP (s) + ins - L, CP (s), -ins, rb, L, q, qb, rg, hap + 1, (int) molecule, NULL, 0, L, mqb);
			}
		}
fastq:
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
static long		readTile = 262144;		/* bases of a tile (CELLCOAL_READ_TILE) */

static long SimulateGenome (int cellIdx0, double depth, int amplified, const char *rg, long int *seed, ReadOut *ro, int donor)
	{
	RecSet	rs;
	FastqOut *fq = &ro->fq;
	AmpEvent *ev[2][2] = {{NULL, NULL}, {NULL, NULL}};
	long	nev[2][2] = {{0, 0}, {0, 0}}, pairs = 0, molecule = ro->molecule, pos;
	int		hap, q, L = readLen, i;
	char	*ta, *tb, *ra, *rb, name[96];
	unsigned char qa[512], qb[512];
	int		qualOn = (qualSd > 0.0 || qualSlope != 0.0 || ctxOn);
	double	qerr = (sequencingError > 0.0) ? sequencingError : 1e-4;
	const long W = 256;
	const int indelsOn = (indelRate > 0.0 || homoError > 0.0 || germIndelRate > 0.0);
	const int homoOn = (homoError > 0.0);
	int		srcs[2] = {cellIdx0, donor}, nsrc = (donor >= 0) ? 2 : 1, si;
	Aln		aln;
	int		indelAln = 0;

	memset (&rs, 0, sizeof(rs));
	q = (int) (-10.0 * log10 (qerr) + 0.5);
	if (q > 41) q = 41;
	if (q < 2) q = 2;
	if (qualOn)
		QualSetup (L, seed);
	ta = (char*) malloc (4 * (size_t) (L + EXT + 2));
	tb = ta + L + EXT + 2; ra = tb + L + EXT + 2; rb = ra + L + EXT + 2;
	/* A doublet (a cell captured together with its donor) has the reads of both genomes in the same files, in one sorted set of records.
	   The genome is simulated in tiles of CELLCOAL_READ_TILE bases (262144 by default) and the BAM records that are complete after each tile
	   are written, so the memory does not depend on the length of the genome. Every source genome and haplotype has its own random stream, which
	   goes on from tile to tile, so the reads do not depend on the size of the tile (or on the number of threads). */
	{
	int			ns = nsrc * 2, k;
	uint64_t	st[4];
	long		cnt[4] = {0, 0, 0, 0}, tile0, tile1, maxcnt = 0;
	double		cmv[2];

	for (si=0; si<nsrc; si++)
		{
		int cellIdx = srcs[si];

		cmv[si] = (cellIdx >= 0 && cellCovMult != NULL) ? cellCovMult[cellIdx] : 1.0;
		for (hap=0; hap<2; hap++)
			ev[si][hap] = amplified ? BuildAmpEvents (cellIdx, hap, seed, &nev[si][hap]) : NULL;
		}
	for (k=0; k<ns; k++)
		{
		uint64_t x = streamState + (uint64_t) (k + 1) * 0xA24BAED4963EE407ULL;

		st[k] = SplitMix64 (&x);
		}
	for (tile0 = 0; tile0 < numSites; tile0 += readTile)
		{
		tile1 = (tile0 + readTile < numSites) ? tile0 + readTile : numSites;
		for (si=0; si<nsrc; si++)
		for (hap=0; hap<2; hap++)
			{
			int		cellIdx = srcs[si], sk = si * 2 + hap;
			double	cm = cmv[si];

			streamState = st[sk];
			rs.key = sk;
			for (pos = tile0; pos < tile1; pos += W)
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
					long	mol;

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
					indelAln = 0;
					if (indelsOn)
						{
						long		nie;
						const IndelEv *ie = IndelsOf (cellIdx, hap, &nie);
						int			r = (nie > 0 || homoOn) ? IndelTemplates (ie, nie, hap, cellIdx, s, ins, L, copyRank, ta, tb, &aln, homoOn) : 0;

						if (r == 3)
							continue;
						indelAln = (r == 2);
						}
					if (!indelAln && LoadTemplates (hap, cellIdx, s, ins, L, ta, tb, copyRank))
						continue;
					if (nev[si][hap] > 0)
						ApplyAmpErrors (ev[si][hap], nev[si][hap], s, ins, L, ta, tb + (indelAln ? aln.extB : 0), seed, indelAln ? aln.ma : NULL, indelAln ? aln.mb + aln.extB : NULL);
					if (ndmg > 0)
						ApplyDamage (L, ta, tb + (indelAln ? aln.extB : 0), seed);
					copies = DuplicateCopies (seed);
					fwdFrag = (RandomUniform (seed) < 0.5);
					mol = molecule + 1 + 4 * cnt[sk] + sk;		/* unique in the cell, and the same for any size of tile */
					cnt[sk]++;
					for (c=0; c<copies; c++)
						{
						if (indelAln)
							{
							if (SlipReads (&aln, ta, tb, L, ra, rb, seed) < 0)
								continue;
							}
						else
							{
							memcpy (ra, ta, L); memcpy (rb, tb, L);
							}
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
						snprintf (name, sizeof(name), "%s:%ld", rg, mol);
						EmitPair (&rs, fq, name, fwdFrag, c, s, ins, ra, rb, q, qualOn ? qa : NULL, qualOn ? qb : NULL, rg, hap, mol, indelAln ? &aln : NULL);
						pairs++;
						}
					}
				}
			st[sk] = streamState;
			}
		if (doBam)
			BamFlush (&ro->bo, &rs, (tile1 < numSites) ? CP (tile1) : LONG_MAX);
		}
	for (k=0; k<ns; k++)
		if (cnt[k] > maxcnt)
			maxcnt = cnt[k];
	molecule += 4 * maxcnt + 4;
	for (si=0; si<2; si++)
		for (i=0; i<2; i++)
			{
			free (ev[si][i]);
			ev[si][i] = NULL;
			}
	}
	if (doBam)
		BamFlush (&ro->bo, &rs, LONG_MAX);
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

	MappabilitySetup ();
	mappabilityOn = (strlen (mappabilityFile) > 0);
	if (indelRate > 0.0 || germIndelRate > 0.0)
		{
		IndelsGenerate (seed);
		IndelsWriteVcf (readDir);
		}

	/* every cell (and the healthy bulk) has its own random stream, taken from the main stream, so the cells are simulated in parallel and
	   the files are the same for any number of threads */
	{
	uint64_t	base;
	long		s1, s2;

	{
	const char *env = getenv ("CELLCOAL_READ_TILE");

	if (env != NULL && atol (env) > 0)
		readTile = ((atol (env) + 255) / 256) * 256;
	}
	if (strlen (damageSpec) > 0 && ndmg == 0)
		LoadDamage ();
	if (strlen (seqContextFile) > 0 && !ctxOn)
		LoadSeqContext ();		/* before the threads start; the table is read only */
	if (qualSd > 0.0 || qualSlope != 0.0 || ctxOn)
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
	IndelsFree ();
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
		int		mrate = data[MATERNAL][OUTGROUP_ROOT][j], ref = mrate, pr = data[PATERNAL][OUTGROUP_ROOT][j], var, nc = 0, nalt = 0, altBase[4] = {-1,-1,-1,-1}, g[2];
		char	gt[8];

		if (refUserGenome != NULL)	/* the reference is the FASTA; the healthy genome may carry germline alleles */
			ref = (int) (strchr (BASES, refUserGenome[j]) - BASES);
		var = (pr != ref || mrate != ref);
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
			int k, gh = (h == 0) ? mrate : pr;

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
