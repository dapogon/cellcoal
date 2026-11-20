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

typedef struct { long site; int alt; double f; } AmpEvent;

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
static void AddRecord (RecSet *rs, const char *name, int flag, long pos, long mpos, long tlen, const char *seq, int len, int q, const char *rg, int hap, int mol)
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
	memset (p + o, q, len); o += len;
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

static void WriteBam (const char *path, RecSet *rs, const char *rg, int sorted)
	{
	Bgzf		bz;
	FILE		*fp;
	unsigned char l4[4], hdr[1024];
	char		text[1024];
	int			tl, nbin = 37450 + 1, i;
	Chunk		*bins;
	uint64_t	*lin;
	long		nlin = (numSites >> 14) + 2, k;
	char		bai[4*MAX_NAME+8];
	unsigned char *buf;
	uint64_t	firstOff = 0, lastEnd = 0;

	if ((fp = fopen (path, "wb")) == NULL)
		{
		fprintf (stderr, "\nERROR: cannot write %s\n", path);
		exit (-1);
		}
	memset (&bz, 0, sizeof(bz));
	bz.fp = fp;
	CrcInit ();
	tl = snprintf (text, sizeof(text), "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:%s\tLN:%d\n@RG\tID:%s\tSM:%s\tPL:ILLUMINA\n@PG\tID:cellcoal\tPN:cellcoal\n", contigName, numSites, rg, rg);
	memcpy (hdr, "BAM\1", 4);
	BgzfWrite (&bz, hdr, 4);
	Put32le (l4, (uint32_t) tl); BgzfWrite (&bz, l4, 4);
	BgzfWrite (&bz, text, tl);
	Put32le (l4, 1); BgzfWrite (&bz, l4, 4);
	Put32le (l4, (uint32_t) strlen (contigName) + 1); BgzfWrite (&bz, l4, 4);
	BgzfWrite (&bz, contigName, strlen (contigName) + 1);
	Put32le (l4, (uint32_t) numSites); BgzfWrite (&bz, l4, 4);

	if (sorted)
		qsort (rs->idx, rs->ni, sizeof(RecIdx), CmpIdx);
	bins = (Chunk*) calloc (nbin, sizeof(Chunk));
	lin = (uint64_t*) calloc (nlin, sizeof(uint64_t));
	buf = rs->buf;
	for (k=0; k<rs->ni; k++)
		{
		RecIdx	*r = &rs->idx[k];
		uint64_t b0 = BgzfTell (&bz), b1;
		long	end = r->pos + (long) (buf[r->off+20] | (buf[r->off+21]<<8) | (buf[r->off+22]<<16) | ((long)buf[r->off+23]<<24));
		int		bin = Reg2Bin (r->pos, end);
		long	w;

		BgzfWrite (&bz, buf + r->off, r->len);
		b1 = BgzfTell (&bz);
		if (k == 0)
			firstOff = b0;
		lastEnd = b1;
		if (bins[bin].end == 0)
			bins[bin].beg = b0;
		bins[bin].end = b1;
		for (w = r->pos >> 14; w <= (end - 1) >> 14; w++)
			if (lin[w] == 0)
				lin[w] = b0;
		}
	BgzfClose (&bz);
	fclose (fp);

	/* BAI */
	snprintf (bai, sizeof(bai), "%s.bai", path);
	if ((fp = fopen (bai, "wb")) == NULL)
		return;
	{
	int nb = 0;
	unsigned char w4[4];
	uint64_t w8;

	for (i=0; i<nbin-1; i++)
		if (bins[i].end)
			nb++;
	fwrite ("BAI\1", 1, 4, fp);
	Put32le (w4, 1); fwrite (w4, 1, 4, fp);
	Put32le (w4, (uint32_t) (nb + 1)); fwrite (w4, 1, 4, fp);
	for (i=0; i<nbin-1; i++)
		if (bins[i].end)
			{
			Put32le (w4, (uint32_t) i); fwrite (w4, 1, 4, fp);
			Put32le (w4, 1); fwrite (w4, 1, 4, fp);
			w8 = bins[i].beg; fwrite (&w8, 1, 8, fp);
			w8 = bins[i].end; fwrite (&w8, 1, 8, fp);
			}
	Put32le (w4, 37450); fwrite (w4, 1, 4, fp);
	Put32le (w4, 2); fwrite (w4, 1, 4, fp);
	w8 = firstOff; fwrite (&w8, 1, 8, fp);
	w8 = lastEnd; fwrite (&w8, 1, 8, fp);
	w8 = (uint64_t) rs->ni; fwrite (&w8, 1, 8, fp);
	w8 = 0; fwrite (&w8, 1, 8, fp);
	{
	long last = (rs->ni > 0) ? ((numSites - 1) >> 14) + 1 : 0;
	uint64_t prev = 0;

	Put32le (w4, (uint32_t) last); fwrite (w4, 1, 4, fp);
	for (k=0; k<last; k++)
		{
		if (lin[k] == 0)
			lin[k] = prev;
		prev = lin[k];
		fwrite (&lin[k], 1, 8, fp);
		}
	}
	}
	fclose (fp);
	free (bins);
	free (lin);
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
				cum += Eij[a][l] / cumEij[a][3];
				if (u <= cum)
					break;
				}
		if (l > 3)
			l = (a + 1) & 3;
		s[p] = BASES[l];
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

static void FastqPut (Bgzf *z, const char *name, int mate, const char *seq, int len, int q)
	{
	char	line[512], qs[1024];
	int		n;

	n = snprintf (line, sizeof(line), "@%s/%d\n", name, mate);
	BgzfWrite (z, line, n);
	BgzfWrite (z, seq, len);
	BgzfWrite (z, "\n+\n", 3);
	memset (qs, q, len);
	BgzfWrite (z, qs, len);
	BgzfWrite (z, "\n", 1);
	}

/* Reads of one cell (cellIdx >= 0) or of the healthy genome (cellIdx = -1). Returns the number of read pairs. */
static long SimulateGenome (int cellIdx, double depth, int amplified, const char *prefix, const char *rg, long int *seed)
	{
	RecSet	rs;
	FastqOut fq;
	AmpEvent *ev[2] = {NULL, NULL};
	long	nev[2] = {0, 0}, pairs = 0, molecule = 0, pos;
	int		hap, q, L = readLen, i;
	char	*ta, *tb, *ra, *rb, name[96], path[4*MAX_NAME];
	double	qerr = (sequencingError > 0.0) ? sequencingError : 1e-4;
	const long W = 256;
	double	cm = (cellIdx >= 0 && cellCovMult != NULL) ? cellCovMult[cellIdx] : 1.0;

	memset (&rs, 0, sizeof(rs));
	memset (&fq, 0, sizeof(fq));
	q = (int) (-10.0 * log10 (qerr) + 0.5);
	if (q > 41) q = 41;
	if (q < 2) q = 2;
	q += 0;	/* BAM stores the raw Phred value */
	ta = (char*) malloc (4 * (size_t) (L + 2));
	tb = ta + L + 2; ra = tb + L + 2; rb = ra + L + 2;
	if (doFastq)
		{
		snprintf (path, sizeof(path), "%s.R1.fastq.gz", prefix);
		fq.f1 = fopen (path, "wb");
		snprintf (path, sizeof(path), "%s.R2.fastq.gz", prefix);
		fq.f2 = fopen (path, "wb");
		if (!fq.f1 || !fq.f2)
			{
			fprintf (stderr, "\nERROR: cannot write the FASTQ files %s\n", prefix);
			exit (-1);
			}
		memset (&fq.r1, 0, sizeof(Bgzf)); memset (&fq.r2, 0, sizeof(Bgzf));
		fq.r1.fp = fq.f1; fq.r2.fp = fq.f2;
		CrcInit ();
		}
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
				}
			lambda = depth * cm * 0.5 * g / (2.0 * L * (1.0 + dupMean)) * (double) W;
			nfrag = RandomPoisson (lambda, seed);
			for (f=0; f<nfrag; f++)
				{
				long	s = pos + (long) (RandomUniform (seed) * W), ins, j, e, k;
				int		lost = 0, copies, c, fwdFrag, nerr = 0;
				long	errSite[16];
				int		errAlt[16];

				ins = (long) (insertMean + insertSd * Normal (seed) + 0.5);
				if (ins < L)
					ins = L;
				if (s + ins > numSites)
					continue;
				for (j=0; j<L && !lost; j++)
					{
					int a = Allele (hap, cellIdx, s + j), b = Allele (hap, cellIdx, s + ins - L + j);

					if (a < 0 || a > 3 || b < 0 || b > 3)
						lost = 1;
					else
						{
						ta[j] = BASES[a];
						tb[j] = BASES[b];
						}
					}
				if (lost)
					continue;
				/* errors of the amplification carried by this molecule (shared by its PCR duplicates) */
				if (nev[hap] > 0)
					{
					long lo = 0, hi = nev[hap];

					while (lo < hi)
						{
						long mid = (lo + hi) / 2;

						if (ev[hap][mid].site < s)
							lo = mid + 1;
						else
							hi = mid;
						}
					for (k=lo; k<nev[hap] && ev[hap][k].site < s + ins && nerr < 16; k++)
						{
						e = ev[hap][k].site;
						if (e < s + L || e >= s + ins - L)
							if (RandomUniform (seed) < ev[hap][k].f)
								{
								errSite[nerr] = e; errAlt[nerr] = ev[hap][k].alt; nerr++;
								}
						}
					for (k=0; k<nerr; k++)
						{
						if (errSite[k] < s + L)
							ta[errSite[k] - s] = BASES[errAlt[k]];
						if (errSite[k] >= s + ins - L)
							tb[errSite[k] - (s + ins - L)] = BASES[errAlt[k]];
						}
					}
				copies = 1;
				if (dupMean > 0.0)
					{
					double qd = dupMean / (1.0 + dupMean);

					while (RandomUniform (seed) < qd && copies < 50)
						copies++;
					}
				fwdFrag = (RandomUniform (seed) < 0.5);
				molecule++;
				for (c=0; c<copies; c++)
					{
					int fa, fb;

					memcpy (ra, ta, L); memcpy (rb, tb, L);
					SeqErrors (ra, L, seed);
					SeqErrors (rb, L, seed);
					snprintf (name, sizeof(name), "%s:%ld", rg, molecule);
					if (fwdFrag)
						{ fa = 99; fb = 147; }
					else
						{ fa = 163; fb = 83; }
					if (c > 0)
						{ fa |= 1024; fb |= 1024; }
					if (doBam)
						{
						AddRecord (&rs, name, fa, s, s + ins - L, ins, ra, L, q, rg, hap + 1, (int) molecule);
						AddRecord (&rs, name, fb, s + ins - L, s, -ins, rb, L, q, rg, hap + 1, (int) molecule);
						}
					if (doFastq)
						{
						char rc[512];

						if (fwdFrag)
							{
							RevComp (rb, rc, L);
							FastqPut (&fq.r1, name, 1, ra, L, q + 33);
							FastqPut (&fq.r2, name, 2, rc, L, q + 33);
							}
						else
							{
							RevComp (rb, rc, L);
							FastqPut (&fq.r1, name, 1, rc, L, q + 33);
							FastqPut (&fq.r2, name, 2, ra, L, q + 33);
							}
						}
					pairs++;
					}
				}
			}
		}
	if (doBam)
		{
		snprintf (path, sizeof(path), "%s.bam", prefix);
		WriteBam (path, &rs, rg, 1);
		}
	if (doFastq)
		{
		BgzfClose (&fq.r1); BgzfClose (&fq.r2);
		fclose (fq.f1); fclose (fq.f2);
		}
	for (i=0; i<2; i++)
		{
		free (ev[i]);
		}
	free (ta);
	free (rs.buf);
	free (rs.idx);
	return pairs;
	}

/********************************** GenerateReadLevel ***********************************/
void GenerateReadLevel (long int *seed, int dataSetNum)
	{
	char	dir[3*MAX_NAME], path[4*MAX_NAME], rg[64];
	FILE	*fp;
	long	j, total = 0;
	int		i;

	snprintf (dir, sizeof(dir), "%s/reads_%d", resultsDir, dataSetNum + 1);
	mkdir (dir, 0755);

	/* the reference is the healthy root genome */
	snprintf (path, sizeof(path), "%s/ref.fa", dir);
	if ((fp = fopen (path, "w")) == NULL)
		{
		fprintf (stderr, "\nERROR: cannot write %s\n", path);
		exit (-1);
		}
	fprintf (fp, ">%s\n", contigName);
	for (j=0; j<numSites; j++)
		{
		int a = data[MATERNAL][OUTGROUP_ROOT][j];

		fputc ((a >= 0 && a < 4) ? BASES[a] : 'N', fp);
		if ((j + 1) % 60 == 0)
			fputc ('\n', fp);
		}
	if (numSites % 60)
		fputc ('\n', fp);
	fclose (fp);
	/* index of the reference (samtools faidx format) */
	snprintf (path, sizeof(path), "%s/ref.fa.fai", dir);
	if ((fp = fopen (path, "w")) != NULL)
		{
		fprintf (fp, "%s\t%d\t%d\t60\t61\n", contigName, numSites, (int) strlen (contigName) + 2);
		fclose (fp);
		}

	for (i=0; i<numCells; i++)
		{
		snprintf (rg, sizeof(rg), "cell%04d", i + 1);
		snprintf (path, sizeof(path), "%s/%s", dir, rg);
		total += SimulateGenome (i, coverage, 1, path, rg, seed);
		if (noisy > 0)
			fprintf (stderr, "\rreads: cell %d/%d", i + 1, numCells);
		}
	if (bulkDepth > 0.0)
		{
		snprintf (path, sizeof(path), "%s/normal", dir);
		SimulateGenome (-1, bulkDepth, 0, path, "normal", seed);
		}
	if (noisy > 0)
		fprintf (stderr, "\rreads: %ld pairs in %d cells -> %s\n", total, numCells, dir);
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

	snprintf (dir, sizeof(dir), "%s/reads_%d", resultsDir, dataSetNum + 1);
	mkdir (dir, 0755);
	snprintf (path, sizeof(path), "%s/truth.vcf", dir);
	if ((fp = fopen (path, "w")) == NULL)
		{
		fprintf (stderr, "\nERROR: cannot write %s\n", path);
		exit (-1);
		}
	fprintf (fp, "##fileformat=VCFv4.2\n##source=cellcoal truth\n##contig=<ID=%s,length=%d>\n", contigName, numSites);
	fprintf (fp, "##INFO=<ID=NC,Number=1,Type=Integer,Description=\"Cells with a non-reference allele\">\n");
	fprintf (fp, "##INFO=<ID=GERMLINE,Number=0,Type=Flag,Description=\"The healthy genome differs from the reference or is heterozygous\">\n");
	fprintf (fp, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"True genotype (maternal|paternal)\">\n");
	fprintf (fp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	for (i=0; i<numCells; i++)
		fprintf (fp, "\tcell%04d", i + 1);
	fprintf (fp, "\n");
	for (j=0; j<numSites; j++)
		{
		int		ref = data[MATERNAL][OUTGROUP_ROOT][j], pr = data[PATERNAL][OUTGROUP_ROOT][j], var = (pr != ref), nc = 0, nalt = 0, altBase[4] = {-1,-1,-1,-1}, g[2];
		char	gt[8];

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
		if (pr >= 0 && pr < 4 && pr != ref)
			{
			int k;

			for (k=0; k<nalt && altBase[k] != pr; k++)
				;
			if (k == nalt && nalt < 4)
				altBase[nalt++] = pr;
			}
		if (nalt == 0)
			continue;
		fprintf (fp, "%s\t%ld\t.\t%c\t", contigName, j + 1, BASES[ref]);
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
