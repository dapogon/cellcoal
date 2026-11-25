/* BAM writer with its index (bam.c), used by the read simulator */
#ifndef BAM_H
#define BAM_H
#include "bgzf.h"

static const char BASES[] = "ACGT";

typedef struct { long pos; size_t off; int len; int key; } RecIdx;

typedef struct
	{
	unsigned char	*buf;
	size_t			n, cap;
	RecIdx			*idx;
	long			ni, capi;
	int				key;		/* source genome and haplotype of the records being added: ties in the position are sorted by it */
	} RecSet;

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

void Grow (unsigned char **b, size_t *cap, size_t need);
void AddRecord (RecSet *rs, const char *name, int flag, long pos, long mpos, long tlen, const char *seq, int len, int q, const unsigned char *qv, const char *rg, int hap, int mol, const uint32_t *cig, int ncig, int reflen, int mapq);
void BamOpen (BamOut *bo, const char *path, const char *rg);
void BamFlush (BamOut *bo, RecSet *rs, long limit);
void BamClose (BamOut *bo);

#endif
