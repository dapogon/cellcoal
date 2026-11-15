/*
 *  bcf.c
 *  CellCoal: conversion of the VCF written by the program to BCF (--bcf)
 *
 *  The VCF is written as text (to a temporary file) and converted here, record by record, to BCF 2.2 inside a BGZF
 *  container (the format of bcftools/htslib), so that no external library is needed except zlib, which is used when it is
 *  available (HAVE_ZLIB); otherwise the BGZF blocks are stored without compression. The converter only needs to
 *  understand the VCF that CellCoal writes (one contig, INFO and FORMAT fields declared in its header).
 *  Floats are written as float32 (the VCF text has 2-3 decimals); integers use the smallest type that holds all the
 *  samples of a field.
 */
#include "cellcoal.h"
#ifdef HAVE_ZLIB
#include <zlib.h>
#endif

/* ---------- BGZF writer ---------- */
typedef struct { FILE *fp; unsigned char buf[0xff00]; int n; } Bgzf;

static uint32_t crcTable[256];
static void CrcInit (void)
	{
	uint32_t c, i, k;

	for (i=0; i<256; i++)
		{
		c = i;
		for (k=0; k<8; k++)
			c = (c & 1) ? 0xEDB88320u ^ (c >> 1) : c >> 1;
		crcTable[i] = c;
		}
	}

static uint32_t Crc32 (const unsigned char *d, size_t n)
	{
	uint32_t c = 0xFFFFFFFFu;
	size_t i;

	for (i=0; i<n; i++)
		c = crcTable[(c ^ d[i]) & 0xFF] ^ (c >> 8);
	return c ^ 0xFFFFFFFFu;
	}

static void Put32le (unsigned char *p, uint32_t v)
	{
	p[0] = v & 0xFF; p[1] = (v >> 8) & 0xFF; p[2] = (v >> 16) & 0xFF; p[3] = (v >> 24) & 0xFF;
	}

/* writes one BGZF block (a gzip member with the BC extra field) with the n bytes of d */
static void BgzfBlock (FILE *fp, const unsigned char *d, int n)
	{
	static const unsigned char head[16] = {0x1f,0x8b,8,4, 0,0,0,0, 0,0xff, 6,0, 'B','C', 2,0};
	unsigned char out[0x10000 + 64], tail[8];
	int clen = -1, total;

#ifdef HAVE_ZLIB
	if (n > 0)
		{
		z_stream zs;

		memset (&zs, 0, sizeof(zs));
		if (deflateInit2 (&zs, 6, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) == Z_OK)
			{
			zs.next_in = (Bytef*) d;
			zs.avail_in = n;
			zs.next_out = out + 18;
			zs.avail_out = 0x10000 - 26;
			if (deflate (&zs, Z_FINISH) == Z_STREAM_END)
				clen = (int) zs.total_out;
			deflateEnd (&zs);
			}
		}
#endif
	if (clen < 0)	/* stored deflate block (also used for the empty end-of-file block) */
		{
		out[18] = 1;
		out[19] = n & 0xFF; out[20] = (n >> 8) & 0xFF;
		out[21] = (~n) & 0xFF; out[22] = ((~n) >> 8) & 0xFF;
		memcpy (out + 23, d, n);
		clen = 5 + n;
		}
	memcpy (out, head, 16);
	total = 18 + clen + 8;
	out[16] = (total - 1) & 0xFF;
	out[17] = ((total - 1) >> 8) & 0xFF;
	Put32le (tail, Crc32 (d, n));
	Put32le (tail + 4, (uint32_t) n);
	memcpy (out + 18 + clen, tail, 8);
	fwrite (out, 1, total, fp);
	}

static void BgzfWrite (Bgzf *b, const void *data, size_t len)
	{
	const unsigned char *p = (const unsigned char*) data;

	while (len > 0)
		{
		size_t room = sizeof(b->buf) - b->n, k = len < room ? len : room;

		memcpy (b->buf + b->n, p, k);
		b->n += (int) k; p += k; len -= k;
		if (b->n == (int) sizeof(b->buf))
			{
			BgzfBlock (b->fp, b->buf, b->n);
			b->n = 0;
			}
		}
	}

static void BgzfClose (Bgzf *b)
	{
	if (b->n > 0)
		BgzfBlock (b->fp, b->buf, b->n);
	{
	/* the standard 28-byte BGZF end-of-file marker (an empty block) */
	static const unsigned char eof[28] = {0x1f,0x8b,0x08,0x04,0,0,0,0,0,0xff,0x06,0,0x42,0x43,0x02,0,0x1b,0,0x03,0,0,0,0,0,0,0,0,0};

	fwrite (eof, 1, sizeof(eof), b->fp);
	}
	}

/* ---------- growable byte buffer ---------- */
typedef struct { unsigned char *d; size_t n, cap; } Bytes;

static void BNeed (Bytes *b, size_t k)
	{
	if (b->n + k > b->cap)
		{
		b->cap = (b->n + k) * 2 + 1024;
		b->d = (unsigned char*) realloc (b->d, b->cap);
		if (!b->d)
			{
			fprintf (stderr, "Out of memory writing the BCF\n");
			exit (-1);
			}
		}
	}
static void BPut8 (Bytes *b, unsigned v) { BNeed (b, 1); b->d[b->n++] = (unsigned char) v; }
static void BPut16 (Bytes *b, unsigned v) { BNeed (b, 2); b->d[b->n] = v & 0xFF; b->d[b->n+1] = (v >> 8) & 0xFF; b->n += 2; }
static void BPut32 (Bytes *b, uint32_t v) { BNeed (b, 4); Put32le (b->d + b->n, v); b->n += 4; }
static void BPutFloat (Bytes *b, float f) { uint32_t u; memcpy (&u, &f, 4); BPut32 (b, u); }
static void BPutRaw (Bytes *b, const void *p, size_t k) { BNeed (b, k); memcpy (b->d + b->n, p, k); b->n += k; }

/* BCF types */
enum { T_NULL = 0, T_INT8 = 1, T_INT16 = 2, T_INT32 = 3, T_FLOAT = 5, T_CHAR = 7 };
#define F_MISSING	0x7F800001u
#define F_END		0x7F800002u

/* type descriptor: length (>= 15 goes in a following typed integer) and type */
static void PutTypedInt (Bytes *b, long v);
static void PutDesc (Bytes *b, int type, long len)
	{
	if (len < 15)
		BPut8 (b, (unsigned)(len << 4) | type);
	else
		{
		BPut8 (b, 0xF0 | type);
		PutTypedInt (b, len);
		}
	}

/* one typed integer scalar, in the smallest type */
static void PutTypedInt (Bytes *b, long v)
	{
	if (v >= -120 && v <= 127)
		{ BPut8 (b, 0x10 | T_INT8); BPut8 (b, (unsigned)(v & 0xFF)); }
	else if (v >= -32760 && v <= 32767)
		{ BPut8 (b, 0x10 | T_INT16); BPut16 (b, (unsigned)(v & 0xFFFF)); }
	else
		{ BPut8 (b, 0x10 | T_INT32); BPut32 (b, (uint32_t) v); }
	}

static void PutTypedString (Bytes *b, const char *s, size_t len)
	{
	PutDesc (b, T_CHAR, (long) len);
	BPutRaw (b, s, len);
	}

/* ---------- dictionaries (from the VCF header) ---------- */
#define MAXDICT 128
enum { K_INT = 1, K_FLOAT = 2, K_STRING = 3, K_FLAG = 4 };
typedef struct { char id[40]; int kind; } DictEntry;
static DictEntry	dict[MAXDICT];
static int			ndict = 0;
static char			contigName[8][64];
static int			ncontig = 0;

static int DictFind (const char *id)
	{
	int i;

	for (i=0; i<ndict; i++)
		if (strcmp (dict[i].id, id) == 0)
			return i;
	return -1;
	}

/* reads ID=..., Type=... from a header line like ##INFO=<ID=AF,Number=R,Type=Float,Description="..."> */
static void HeaderEntry (const char *line, int hasType)
	{
	const char *p = strstr (line, "ID="), *t;
	int k = 0;

	if (!p || ndict >= MAXDICT)
		return;
	p += 3;
	if (ndict == 0)
		{
		strcpy (dict[0].id, "PASS");
		dict[0].kind = 0;
		ndict = 1;
		}
	while (*p && *p != ',' && *p != '>' && k < 39)
		dict[ndict].id[k++] = *p++;
	dict[ndict].id[k] = '\0';
	dict[ndict].kind = 0;
	if (hasType && (t = strstr (line, "Type=")) != NULL)
		{
		t += 5;
		if (strncmp (t, "Integer", 7) == 0) dict[ndict].kind = K_INT;
		else if (strncmp (t, "Float", 5) == 0) dict[ndict].kind = K_FLOAT;
		else if (strncmp (t, "Flag", 4) == 0) dict[ndict].kind = K_FLAG;
		else dict[ndict].kind = K_STRING;
		}
	if (DictFind (dict[ndict].id) < 0)
		ndict++;
	}

/* writes the header line with its IDX appended (before the closing '>') into text */
static void HeaderLineWithIdx (Bytes *text, const char *line, int idx)
	{
	const char *gt = strrchr (line, '>');
	char add[32];

	if (!gt)
		{
		BPutRaw (text, line, strlen (line));
		return;
		}
	snprintf (add, sizeof(add), ",IDX=%d", idx);
	BPutRaw (text, line, (size_t)(gt - line));
	BPutRaw (text, add, strlen (add));
	BPutRaw (text, gt, strlen (gt));
	}

/* ---------- one VCF record to BCF ---------- */
static int Split (char *s, char sep, char **out, int max)
	{
	int n = 0;

	out[n++] = s;
	while (*s && n < max)
		{
		if (*s == sep)
			{
			*s = '\0';
			out[n++] = s + 1;
			}
		s++;
		}
	return n;
	}

static void FailLine (long lineNo, const char *msg)
	{
	fprintf (stderr, "\n!!! BCF conversion: %s (line %ld of the VCF)\n\n", msg, lineNo);
	exit (-1);
	}

/* integers of one FORMAT field of all samples: writes the type descriptor and the values */
static void PutIntVectors (Bytes *b, long *vals, unsigned char *miss, int ns, int len)
	{
	long lo = 0, hi = 0;
	int i, type;

	for (i=0; i<ns*len; i++)
		if (!miss[i])
			{
			if (vals[i] < lo) lo = vals[i];
			if (vals[i] > hi) hi = vals[i];
			}
	type = (lo >= -120 && hi <= 127) ? T_INT8 : (lo >= -32760 && hi <= 32767) ? T_INT16 : T_INT32;
	PutDesc (b, type, len);
	for (i=0; i<ns*len; i++)
		{
		/* miss: 0 value, 1 missing, 2 vector end */
		if (type == T_INT8)
			BPut8 (b, miss[i] == 1 ? 0x80 : miss[i] == 2 ? 0x81 : (unsigned)(vals[i] & 0xFF));
		else if (type == T_INT16)
			BPut16 (b, miss[i] == 1 ? 0x8000 : miss[i] == 2 ? 0x8001 : (unsigned)(vals[i] & 0xFFFF));
		else
			BPut32 (b, miss[i] == 1 ? 0x80000000u : miss[i] == 2 ? 0x80000001u : (uint32_t) vals[i]);
		}
	}

static void ConvertRecord (Bytes *rec, char *line, int nsamples, long lineNo)
	{
	enum { MAXKEYS = 32 };
	static char		**tok = NULL;
	char			*alts[8], *filt[8], *info[16], *keys[MAXKEYS], *vals[64];
	Bytes			shared = {0,0,0}, indiv = {0,0,0};
	int				ntok, nalt, nallele, ninfo, nfilt, nkeys, k, s, v, idx, nv;
	long			pos;

	if (!tok)
		tok = (char**) malloc ((size_t)(nsamples + 16) * sizeof(char*));
	ntok = Split (line, '\t', tok, nsamples + 16);
	if (ntok < 9 + nsamples)
		FailLine (lineNo, "too few fields");

	/* CHROM, POS, rlen, QUAL, counts */
	idx = -1;
	for (k=0; k<ncontig; k++)
		if (strcmp (contigName[k], tok[0]) == 0)
			idx = k;
	if (idx < 0)
		FailLine (lineNo, "unknown contig");
	BPut32 (&shared, (uint32_t) idx);
	pos = atol (tok[1]);
	BPut32 (&shared, (uint32_t)(pos - 1));
	BPut32 (&shared, (uint32_t) strlen (tok[3]));
	if (strcmp (tok[5], ".") == 0)
		BPut32 (&shared, F_MISSING);
	else
		BPutFloat (&shared, strtof (tok[5], NULL));

	/* alleles */
	if (strcmp (tok[4], ".") == 0)
		nalt = 0;
	else
		nalt = Split (tok[4], ',', alts, 8);
	nallele = 1 + nalt;
	nfilt = (strcmp (tok[6], ".") == 0) ? 0 : Split (tok[6], ';', filt, 8);
	ninfo = (strcmp (tok[7], ".") == 0) ? 0 : Split (tok[7], ';', info, 16);

	{
	Bytes	body = {0,0,0};

	/* ID */
	if (strcmp (tok[2], ".") == 0)
		PutTypedString (&body, "", 0);
	else
		PutTypedString (&body, tok[2], strlen (tok[2]));
	PutTypedString (&body, tok[3], strlen (tok[3]));
	for (k=0; k<nalt; k++)
		PutTypedString (&body, alts[k], strlen (alts[k]));
	/* FILTER */
	if (nfilt == 0)
		BPut8 (&body, 0);
	else
		{
		PutDesc (&body, T_INT8, nfilt);
		for (k=0; k<nfilt; k++)
			{
			int d = DictFind (filt[k]);

			if (d < 0)
				FailLine (lineNo, "undeclared FILTER");
			BPut8 (&body, (unsigned) d);
			}
		}
	/* INFO */
	for (k=0; k<ninfo; k++)
		{
		char *eq = strchr (info[k], '='), *val = NULL;
		int d, kind;

		if (eq)
			{ *eq = '\0'; val = eq + 1; }
		d = DictFind (info[k]);
		if (d < 0)
			FailLine (lineNo, "undeclared INFO key");
		kind = dict[d].kind;
		PutTypedInt (&body, d);
		if (kind == K_FLAG || val == NULL)
			BPut8 (&body, 0);
		else if (kind == K_STRING)
			PutTypedString (&body, val, strlen (val));
		else
			{
			nv = Split (val, ',', vals, 64);
			if (kind == K_INT)
				{
				long lo = 0, hi = 0, x;
				int type;

				for (v=0; v<nv; v++)
					{
					x = atol (vals[v]);
					if (x < lo) lo = x;
					if (x > hi) hi = x;
					}
				type = (lo >= -120 && hi <= 127) ? T_INT8 : (lo >= -32760 && hi <= 32767) ? T_INT16 : T_INT32;
				PutDesc (&body, type, nv);
				for (v=0; v<nv; v++)
					{
					x = atol (vals[v]);
					if (type == T_INT8) BPut8 (&body, (unsigned)(x & 0xFF));
					else if (type == T_INT16) BPut16 (&body, (unsigned)(x & 0xFFFF));
					else BPut32 (&body, (uint32_t) x);
					}
				}
			else
				{
				PutDesc (&body, T_FLOAT, nv);
				for (v=0; v<nv; v++)
					BPutFloat (&body, strtof (vals[v], NULL));
				}
			}
		}
	BPut16 (&shared, (unsigned) ninfo);
	BPut16 (&shared, (unsigned) nallele);
	/* n_sample (24 bits) and n_fmt (8 bits) are written after FORMAT is known */
	{
	size_t where = shared.n;

	BPut32 (&shared, 0);
	BPutRaw (&shared, body.d, body.n);
	free (body.d);

	/* FORMAT: keys and, for every key, the values of all the samples */
	nkeys = Split (tok[8], ':', keys, MAXKEYS);
	{
	char	***cell = (char***) malloc (nsamples * sizeof(char**));
	char	*parts[MAXKEYS];
	int		*kd = (int*) malloc (nkeys * sizeof(int));

	for (s=0; s<nsamples; s++)
		{
		int np;

		cell[s] = (char**) malloc (nkeys * sizeof(char*));
		np = Split (tok[9+s], ':', parts, MAXKEYS);
		for (k=0; k<nkeys; k++)
			cell[s][k] = (k < np) ? parts[k] : (char*) ".";
		/* the text of the parts lives in the line buffer, which stays valid until the record is done */
		}
	for (k=0; k<nkeys; k++)
		{
		int kind, len = 1;

		kd[k] = DictFind (keys[k]);
		if (kd[k] < 0)
			FailLine (lineNo, "undeclared FORMAT key");
		kind = dict[kd[k]].kind;
		PutTypedInt (&indiv, kd[k]);
		if (strcmp (keys[k], "GT") == 0)
			{
			/* genotypes are integer vectors in BCF: (allele + 1) << 1 | phased, 0 for a missing allele (htslib requires this for GT) */
			int ploidy = 1;
			char *z;

			for (z=cell[0][k]; *z; z++)
				if (*z == '|' || *z == '/')
					ploidy++;
			PutDesc (&indiv, T_INT8, ploidy);
			for (s=0; s<nsamples; s++)
				{
				int np = 0;

				z = cell[s][k];
				while (np < ploidy)
					{
					int allele = -1, phased = 0;

					if (np > 0)
						{
						if (*z != '|' && *z != '/')
							FailLine (lineNo, "inconsistent ploidy in GT");
						phased = (*z == '|');
						z++;
						}
					if (*z == '.')
						z++;
					else if (*z >= '0' && *z <= '9')
						{
						allele = 0;
						while (*z >= '0' && *z <= '9')
							allele = allele * 10 + (*z++ - '0');
						}
					else
						FailLine (lineNo, "bad GT");
					BPut8 (&indiv, (unsigned)(((allele + 1) << 1) | phased));
					np++;
					}
				}
			}
		else if (kind == K_STRING)
			{
			size_t mx = 0;

			for (s=0; s<nsamples; s++)
				if (strlen (cell[s][k]) > mx)
					mx = strlen (cell[s][k]);
			PutDesc (&indiv, T_CHAR, (long) mx);
			for (s=0; s<nsamples; s++)
				{
				size_t l = strlen (cell[s][k]);

				BPutRaw (&indiv, cell[s][k], l);
				for (; l<mx; l++)
					BPut8 (&indiv, 0);
				}
			}
		else
			{
			long	*iv;
			unsigned char *miss;
			char	tmp[8192];

			for (s=0; s<nsamples; s++)
				{
				int c = 1;
				char *z;

				for (z=cell[s][k]; *z; z++)
					if (*z == ',')
						c++;
				if (c > len)
					len = c;
				}
			iv = (long*) malloc ((size_t) nsamples * len * sizeof(long));
			miss = (unsigned char*) malloc ((size_t) nsamples * len);
			for (s=0; s<nsamples; s++)
				{
				int nn;
				char *vv[64];

				if (strlen (cell[s][k]) >= sizeof(tmp))
					FailLine (lineNo, "field too long");
				strcpy (tmp, cell[s][k]);
				nn = Split (tmp, ',', vv, 64);
				for (v=0; v<len; v++)
					{
					long *pv = &iv[(size_t) s * len + v];
					unsigned char *pm = &miss[(size_t) s * len + v];

					if (v >= nn)
						{ *pv = 0; *pm = 2; }
					else if (strcmp (vv[v], ".") == 0)
						{ *pv = 0; *pm = 1; }
					else
						{
						*pm = 0;
						if (kind == K_FLOAT)
							{
							float f = strtof (vv[v], NULL);
							memcpy (pv, &f, sizeof(float));
							}
						else
							*pv = atol (vv[v]);
						}
					}
				}
			if (kind == K_FLOAT)
				{
				PutDesc (&indiv, T_FLOAT, len);
				for (s=0; s<nsamples*len; s++)
					{
					float f;

					if (miss[s] == 1) BPut32 (&indiv, F_MISSING);
					else if (miss[s] == 2) BPut32 (&indiv, F_END);
					else { memcpy (&f, &iv[s], sizeof(float)); BPutFloat (&indiv, f); }
					}
				}
			else
				PutIntVectors (&indiv, iv, miss, nsamples, len);
			free (iv); free (miss);
			}
		}
	for (s=0; s<nsamples; s++)
		free (cell[s]);
	free (cell); free (kd);
	}
	Put32le (shared.d + where, (uint32_t) nsamples | ((uint32_t) nkeys << 24));
	}
	}
	rec->n = 0;
	BPut32 (rec, (uint32_t) shared.n);
	BPut32 (rec, (uint32_t) indiv.n);
	BPutRaw (rec, shared.d, shared.n);
	BPutRaw (rec, indiv.d, indiv.n);
	free (shared.d); free (indiv.d);
	}

/* Converts the VCF text file vcfPath into the BCF file bcfPath; returns 0 if all went well */
/* ---------- CSI index ---------- */
/* Virtual file offset (BGZF block address << 16 | offset inside the block) of the next byte written */
static uint64_t BgzfTell (Bgzf *b)
	{
	return ((uint64_t) ftello (b->fp) << 16) | (uint64_t) b->n;
	}

typedef struct { int tid; uint32_t bin; uint64_t beg, end; } CsiChunk;
typedef struct
	{
	CsiChunk *c;
	size_t	n, cap;
	int		curTid;
	int64_t	curBin;
	uint64_t curBeg;
	uint64_t firstOff[8], lastEnd[8];	/* range of the offsets of the records of each contig */
	long	nrec[8];
	} CsiIndex;

/* bin of the interval [beg,end) (0-based) in a CSI index of the given minimum shift and depth */
static int CsiReg2Bin (int64_t beg, int64_t end, int minShift, int depth)
	{
	int		l, s = minShift;
	int64_t	t = (((int64_t)1 << (3*depth)) - 1) / 7;	/* first bin of the lowest level */

	for (--end, l=depth; l>0; --l, s+=3, t -= (int64_t)1 << (3*l))
		if (beg >> s == end >> s)
			return (int) (t + (beg >> s));
	return 0;
	}

static void CsiPushChunk (CsiIndex *x, int tid, uint32_t bin, uint64_t beg, uint64_t end)
	{
	if (x->n == x->cap)
		{
		x->cap = x->cap ? 2 * x->cap : 1024;
		x->c = (CsiChunk*) realloc (x->c, x->cap * sizeof(CsiChunk));
		if (!x->c)
			{
			fprintf (stderr, "Out of memory writing the BCF index\n");
			exit (-1);
			}
		}
	x->c[x->n].tid = tid; x->c[x->n].bin = bin; x->c[x->n].beg = beg; x->c[x->n].end = end;
	x->n++;
	}

/* registers a record (starting at the virtual offset off and ending at endOff) with its contig, 0-based position and length;
   consecutive records of the same bin share a chunk */
static void CsiAdd (CsiIndex *x, int minShift, int depth, int tid, int64_t pos, int64_t rlen, uint64_t off, uint64_t endOff)
	{
	int64_t bin = CsiReg2Bin (pos, pos + (rlen > 0 ? rlen : 1), minShift, depth);
	int		first = (x->nrec[tid] == 0);
	long	total = 0;
	int		t;

	for (t=0; t<8; t++)
		total += x->nrec[t];
	if (first)
		x->firstOff[tid] = off;
	if (total > 0 && (tid != x->curTid || bin != x->curBin))
		CsiPushChunk (x, x->curTid, (uint32_t) x->curBin, x->curBeg, off);
	if (total == 0 || tid != x->curTid || bin != x->curBin)
		{
		x->curTid = tid; x->curBin = bin; x->curBeg = off;
		}
	x->lastEnd[tid] = endOff;
	x->nrec[tid]++;
	}

static int CmpChunk (const void *a, const void *b)
	{
	const CsiChunk *x = (const CsiChunk*) a, *y = (const CsiChunk*) b;

	if (x->tid != y->tid)
		return x->tid < y->tid ? -1 : 1;
	if (x->bin != y->bin)
		return x->bin < y->bin ? -1 : 1;
	return (x->beg > y->beg) - (x->beg < y->beg);
	}

static void BPut64 (Bytes *b, uint64_t v)
	{
	BPut32 (b, (uint32_t) (v & 0xFFFFFFFFu));
	BPut32 (b, (uint32_t) (v >> 32));
	}

/* Writes the CSI index (min_shift 14 like "bcftools index"; the depth grows to cover the contig) in a BGZF file.
   CellCoal records are one base long, so they fall in the bins of the lowest level; every bin is a run of consecutive
   records and has a chunk per run, and its loffset is the first record of the bin. Each contig also has the usual
   pseudo-bin with its range of offsets and its number of records. */
static int WriteCsi (CsiIndex *x, int minShift, int depth, int nref, const char *path)
	{
	FILE	*fp = fopen (path, "wb");
	Bytes	b = {0,0,0};
	Bgzf	bz;
	int		tid;
	size_t	i = 0, j, k, first;
	uint32_t metaBin = (uint32_t) ((((int64_t)1 << (3*(depth+1))) - 1) / 7 + 1);

	if (!fp)
		return -1;
	{
	long total = 0;
	int t;

	for (t=0; t<8; t++)
		total += x->nrec[t];
	if (total > 0)
		CsiPushChunk (x, x->curTid, (uint32_t) x->curBin, x->curBeg, x->lastEnd[x->curTid]);
	}
	if (x->n > 1)
		qsort (x->c, x->n, sizeof(CsiChunk), CmpChunk);

	BPutRaw (&b, "CSI\1", 4);
	BPut32 (&b, (uint32_t) minShift);
	BPut32 (&b, (uint32_t) depth);
	BPut32 (&b, 0);			/* l_aux */
	BPut32 (&b, (uint32_t) nref);
	for (tid=0; tid<nref; tid++)
		{
		size_t nbin = 0;

		first = i;
		for (j=i; j<x->n && x->c[j].tid == tid; j++)
			if (j == first || x->c[j].bin != x->c[j-1].bin)
				nbin++;
		i = j;
		BPut32 (&b, (uint32_t) (nbin + (x->nrec[tid] > 0 ? 1 : 0)));
		for (j=first; j<i; j=k)
			{
			uint32_t bin = x->c[j].bin;

			for (k=j; k<i && x->c[k].bin == bin; k++)
				;
			BPut32 (&b, bin);
			BPut64 (&b, x->c[j].beg);			/* loffset */
			BPut32 (&b, (uint32_t) (k - j));
			for (; j<k; j++)
				{
				BPut64 (&b, x->c[j].beg);
				BPut64 (&b, x->c[j].end);
				}
			}
		if (x->nrec[tid] > 0)
			{
			BPut32 (&b, metaBin);
			BPut64 (&b, 0);
			BPut32 (&b, 2);
			BPut64 (&b, x->firstOff[tid]);
			BPut64 (&b, x->lastEnd[tid]);
			BPut64 (&b, (uint64_t) x->nrec[tid]);	/* mapped */
			BPut64 (&b, 0);						/* unmapped */
			}
		}
	BPut64 (&b, 0);			/* records without coordinate */

	bz.fp = fp; bz.n = 0;
	BgzfWrite (&bz, b.d, b.n);
	BgzfClose (&bz);
	free (b.d);
	return fclose (fp);
	}

int ConvertVCFtoBCF (const char *vcfPath, const char *bcfPath)
	{
	FILE	*in = fopen (vcfPath, "r"), *out = fopen (bcfPath, "wb");
	char	*line = NULL;
	size_t	cap = 0;
	ssize_t	len;
	long	lineNo = 0;
	int		nsamples = 0, headerDone = NO, nIdx;
	Bytes	text = {0,0,0}, rec = {0,0,0};
	Bgzf	bz;
	char	magic[5] = {'B','C','F',2,2};
	CsiIndex csi;
	int		csiDepth = 5;

	if (!in || !out)
		return -1;
	CrcInit ();
	memset (&csi, 0, sizeof(csi));
	while (((int64_t)1 << (14 + 3*csiDepth)) <= (int64_t) numSites + 1)
		csiDepth++;
	ndict = 0; ncontig = 0;
	bz.fp = out; bz.n = 0;

	/* header: the "##" lines (with IDX added), the "##FILTER" for PASS, and the #CHROM line */
	nIdx = 1;
	while ((len = getline (&line, &cap, in)) > 0)
		{
		lineNo++;
		if (line[0] != '#')
			continue;	/* e.g. the replicate markers before the header */
		if (strncmp (line, "##", 2) == 0)
			{
			if (strncmp (line, "##contig=", 9) == 0)
				{
				const char *p = strstr (line, "ID=");
				int k = 0;

				if (p && ncontig < 8)
					{
					p += 3;
					while (*p && *p != ',' && *p != '>' && k < 63)
						contigName[ncontig][k++] = *p++;
					contigName[ncontig][k] = '\0';
					HeaderLineWithIdx (&text, line, ncontig);
					ncontig++;
					}
				}
			else if (strncmp (line, "##INFO=", 7) == 0 || strncmp (line, "##FORMAT=", 9) == 0 || strncmp (line, "##FILTER=", 9) == 0)
				{
				HeaderEntry (line, strncmp (line, "##FILTER=", 9) != 0);
				HeaderLineWithIdx (&text, line, nIdx++);
				}
			else
				{
				BPutRaw (&text, line, strlen (line));
				if (strncmp (line, "##fileformat", 12) == 0)
					{
					static const char pass[] = "##FILTER=<ID=PASS,Description=\"All filters passed\",IDX=0>\n";

					BPutRaw (&text, pass, sizeof(pass) - 1);
					}
				}
			}
		else	/* #CHROM line */
			{
			char *t;
			int tabs = 0;

			for (t=line; *t; t++)
				if (*t == '\t')
					tabs++;
			nsamples = tabs - 8;
			BPutRaw (&text, line, strlen (line));
			headerDone = YES;
			break;
			}
		}
	if (headerDone == NO || nsamples < 1)
		{
		fprintf (stderr, "\n!!! BCF conversion: the VCF has no header\n\n");
		return -1;
		}
	BPut8 (&text, 0);
	{
	unsigned char l4[4];

	BgzfWrite (&bz, magic, 5);
	Put32le (l4, (uint32_t) text.n);
	BgzfWrite (&bz, l4, 4);
	BgzfWrite (&bz, text.d, text.n);
	}
	while ((len = getline (&line, &cap, in)) > 0)
		{
		lineNo++;
		while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r'))
			line[--len] = '\0';
		if (len == 0 || line[0] == '#' || strchr (line, '\t') == NULL)
			continue;	/* blank lines and replicate markers */
		ConvertRecord (&rec, line, nsamples, lineNo);
		{
		uint64_t off = BgzfTell (&bz);
		uint32_t tid, pos, rlen;

		BgzfWrite (&bz, rec.d, rec.n);
		memcpy (&tid, rec.d + 8, 4); memcpy (&pos, rec.d + 12, 4); memcpy (&rlen, rec.d + 16, 4);	/* little-endian hosts */
		if (tid < 8)
			CsiAdd (&csi, 14, csiDepth, (int) tid, (int64_t) pos, (int64_t) rlen, off, BgzfTell (&bz));
		}
		}
	BgzfClose (&bz);
	free (text.d); free (rec.d); free (line);
	fclose (in);
	if (fclose (out) != 0)
		return -1;
	{
	char	csiPath[4*MAX_NAME+16];
	int		rc;

	snprintf (csiPath, sizeof(csiPath), "%s.csi", bcfPath);
	rc = WriteCsi (&csi, 14, csiDepth, ncontig, csiPath);
	free (csi.c);
	return rc;
	}
	}
