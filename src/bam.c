/* BAM writer: records, sorted flush in BGZF blocks, and the BAI index */

#include <limits.h>
#include <string.h>
#include <math.h>
#include "cellcoal.h"
#include "bam.h"

void Grow (unsigned char **b, size_t *cap, size_t need)
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
void AddRecord (RecSet *rs, const char *name, int flag, long pos, long mpos, long tlen, const char *seq, int len, int q, const unsigned char *qv, const char *rg, int hap, int mol, const uint32_t *cig, int ncig, int reflen, int mapq)
	{
	size_t	need, o;
	int		nl = (int) strlen (name) + 1, rgl = (int) strlen (rg) + 1, i;
	unsigned char *p;

	need = 4 + 32 + nl + 4 * (ncig > 0 ? ncig : 1) + (len + 1) / 2 + len + 3 + rgl + 3 + 8 + 16;
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
	p[o++] = (unsigned char) nl; p[o++] = (unsigned char) mapq;								/* l_read_name, mapq */
	Put16 (p + o, Reg2Bin (pos, pos + reflen)); o += 2;
	Put16 (p + o, ncig > 0 ? ncig : 1); o += 2;												/* n_cigar_op */
	Put16 (p + o, flag); o += 2;
	Put32le (p + o, (uint32_t) len); o += 4;
	Put32le (p + o, 0); o += 4;												/* next refID */
	Put32le (p + o, (uint32_t) mpos); o += 4;
	Put32le (p + o, (uint32_t) (int32_t) tlen); o += 4;
	memcpy (p + o, name, nl); o += nl;
	if (ncig > 0)														/* CIGAR */
		{
		int ci;

		for (ci=0; ci<ncig; ci++)
			{
			Put32le (p + o, cig[ci]); o += 4;
			}
		}
	else
		{
		Put32le (p + o, (uint32_t) len << 4); o += 4;						/* len M */
		}
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
	rs->idx[rs->ni].key = rs->key;
	rs->ni++;
	rs->n += o;
	}

static int CmpIdx (const void *a, const void *b)
	{
	const RecIdx *x = (const RecIdx*) a, *y = (const RecIdx*) b;

	if (x->pos != y->pos)
		return x->pos < y->pos ? -1 : 1;
	if (x->key != y->key)
		return x->key < y->key ? -1 : 1;
	return x->off < y->off ? -1 : (x->off > y->off);
	}

/* ---------------------------- BAM writer with its index ---------------------------- */

long ReadRefLen (void)		/* length of the reference sequence of the BAM header */
	{
	return blocksOn ? chromLength : (long) numSites;
	}

const char *ReadContig (void)
	{
	return blocksOn ? BlocksChrom () : contigName;
	}

void BamOpen (BamOut *bo, const char *path, const char *rg)
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
/* Writes the records of the set that start before limit (all of them if limit is LONG_MAX), sorted, and keeps the others for the next call. The records
   that start at or after the limit are those of the molecules that have not been simulated yet or the second reads of molecules that have. */
void BamFlush (BamOut *bo, RecSet *rs, long limit)
	{
	unsigned char *buf = rs->buf;
	long		k, nw = 0;

	qsort (rs->idx, rs->ni, sizeof(RecIdx), CmpIdx);
	for (k=0; k<rs->ni; k++)
		{
		RecIdx	*r = &rs->idx[k];
		uint64_t b0, b1;
		long	end, w;
		int		bin;

		if (r->pos >= limit)
			break;
		b0 = BgzfTell (&bo->bz);
		end = r->pos + (long) (buf[r->off+20] | (buf[r->off+21]<<8) | (buf[r->off+22]<<16) | ((long)buf[r->off+23]<<24));
		bin = Reg2Bin (r->pos, end);
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
		nw++;
		}
	if (nw == rs->ni)
		{
		rs->n = 0;
		rs->ni = 0;
		return;
		}
	{							/* the records that stay move to the start of the buffer */
	size_t	keep = 0, o = 0;
	unsigned char *tmp;
	long	m = rs->ni - nw;

	for (k=nw; k<rs->ni; k++)
		keep += rs->idx[k].len;
	tmp = (unsigned char*) malloc (keep);
	for (k=nw; k<rs->ni; k++)
		{
		memcpy (tmp + o, buf + rs->idx[k].off, rs->idx[k].len);
		rs->idx[k].off = o;
		o += rs->idx[k].len;
		}
	memcpy (buf, tmp, keep);
	free (tmp);
	memmove (rs->idx, rs->idx + nw, m * sizeof(RecIdx));
	rs->ni = m;
	rs->n = keep;
	}
	}

void BamClose (BamOut *bo)
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
