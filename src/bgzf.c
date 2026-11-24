/*
 *  bgzf.c
 *  CellCoal: BGZF writer (blocks of up to 64 kB compressed with zlib when it is available), shared by the BCF and BAM writers
 */
#include "cellcoal.h"
#include "bgzf.h"
#ifdef HAVE_ZLIB
#include <zlib.h>
#endif

/* ---------- BGZF writer ---------- */

static uint32_t crcTable[256];
void CrcInit (void)
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

uint32_t Crc32 (const unsigned char *d, size_t n)
	{
#ifdef HAVE_ZLIB
	return (uint32_t) crc32 (0L, (const Bytef*) d, (uInt) n);
#else
	uint32_t c = 0xFFFFFFFFu;
	size_t i;

	for (i=0; i<n; i++)
		c = crcTable[(c ^ d[i]) & 0xFF] ^ (c >> 8);
	return c ^ 0xFFFFFFFFu;
#endif
	}

void Put32le (unsigned char *p, uint32_t v)
	{
	p[0] = v & 0xFF; p[1] = (v >> 8) & 0xFF; p[2] = (v >> 16) & 0xFF; p[3] = (v >> 24) & 0xFF;
	}

/* writes one BGZF block (a gzip member with the BC extra field) with the n bytes of d */
void BgzfBlock (FILE *fp, const unsigned char *d, int n)
	{
	static const unsigned char head[16] = {0x1f,0x8b,8,4, 0,0,0,0, 0,0xff, 6,0, 'B','C', 2,0};
	unsigned char out[0x10000 + 64], tail[8];
	int clen = -1, total;

#ifdef HAVE_ZLIB
	if (n > 0)
		{
		/* one deflate stream for every thread, reset for each block: making and freeing the 256 kB state of zlib for every block of 64 kB
		   costs as much as a part of the compression itself, and the bytes that are written are the same */
		static __thread z_stream zs;
		static __thread int zsOk = 0, zsLevel = 0;

		if (zsOk && zsLevel != compressLevel)
			{
			deflateEnd (&zs);
			zsOk = 0;
			}
		if (!zsOk)
			{
			memset (&zs, 0, sizeof(zs));
			zsOk = (deflateInit2 (&zs, compressLevel, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) == Z_OK);
			zsLevel = compressLevel;
			}
		else
			zsOk = (deflateReset (&zs) == Z_OK);
		if (zsOk)
			{
			zs.next_in = (Bytef*) d;
			zs.avail_in = n;
			zs.next_out = out + 18;
			zs.avail_out = 0x10000 - 26;
			if (deflate (&zs, Z_FINISH) == Z_STREAM_END)
				clen = (int) zs.total_out;
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

void BgzfWrite (Bgzf *b, const void *buf, size_t len)
	{
	const unsigned char *p = (const unsigned char*) buf;

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

void BgzfClose (Bgzf *b)
	{
	if (b->n > 0)
		BgzfBlock (b->fp, b->buf, b->n);
	{
	/* the standard 28-byte BGZF end-of-file marker (an empty block) */
	static const unsigned char eof[28] = {0x1f,0x8b,0x08,0x04,0,0,0,0,0,0xff,0x06,0,0x42,0x43,0x02,0,0x1b,0,0x03,0,0,0,0,0,0,0,0,0};

	fwrite (eof, 1, sizeof(eof), b->fp);
	}
	}


/* Virtual file offset (BGZF block address << 16 | offset inside the block) of the next byte written */
uint64_t BgzfTell (Bgzf *b)
	{
	return ((uint64_t) ftello (b->fp) << 16) | (uint64_t) b->n;
	}

