#ifndef bgzf_h
#define bgzf_h
#include <stdio.h>
#include <stdint.h>
#include <stddef.h>

typedef struct { FILE *fp; unsigned char buf[0xff00]; int n; } Bgzf;

void		Put32le (unsigned char *p, uint32_t v);
void		BgzfBlock (FILE *fp, const unsigned char *d, int n);
void		BgzfWrite (Bgzf *b, const void *data, size_t len);
void		BgzfClose (Bgzf *b);
void		CrcInit (void);
uint32_t	Crc32 (const unsigned char *d, size_t n);
uint64_t	BgzfTell (Bgzf *b);

#endif
