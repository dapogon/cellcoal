/*
	CellCoal - the few things that differ between POSIX systems (Linux, macOS) and Windows (MSYS2/MinGW-w64).
*/
#ifndef platform_h
#define platform_h

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#ifdef _WIN32
	#include <direct.h>
	#include <io.h>
	#define MAKE_DIR(path)		_mkdir (path)
	#define LSTAT(path, st)		stat (path, st)
	#define SHELL_QUOTE_CHAR	'"'
	#include <stdio.h>
	#include <stdlib.h>
	#include <string.h>
	#define fseeko(f, off, whence)	_fseeki64 (f, (__int64) (off), whence)
	#define ftello(f)				((long long) _ftelli64 (f))
	/* getline is not in the Windows C library */
	static inline long long PortableGetline (char **line, size_t *cap, FILE *fp)
		{
		size_t len = 0;
		int c;

		if (*line == NULL || *cap == 0)
			{
			*cap = 256;
			*line = (char*) malloc (*cap);
			if (*line == NULL)
				return -1;
			}
		while ((c = fgetc (fp)) != EOF)
			{
			if (len + 2 > *cap)
				{
				char *bigger = (char*) realloc (*line, *cap * 2);

				if (bigger == NULL)
					return -1;
				*line = bigger;
				*cap *= 2;
				}
			(*line)[len++] = (char) c;
			if (c == '\n')
				break;
			}
		if (len == 0)
			return -1;
		(*line)[len] = '\0';
		return (long long) len;
		}
	#define getline(l, c, f)		PortableGetline (l, c, f)
	#define gmtime_r(t, r)		(gmtime_s (r, t) == 0 ? (r) : NULL)
	#define localtime_r(t, r)	(localtime_s (r, t) == 0 ? (r) : NULL)
#else
	#define MAKE_DIR(path)		mkdir (path, S_IRWXU)
	#define LSTAT(path, st)		lstat (path, st)
	#define SHELL_QUOTE_CHAR	'\''
#endif

/* In-memory FILE streams (a line of the VCF/CATG is formatted in a thread and written later in order).
   POSIX: open_memstream. Windows has none: a temporary file is used and its content is read back when the stream is closed.
   The buffer and its length are valid after MemStreamClose; the caller frees the buffer. */
#ifdef _WIN32
	#include <stdio.h>
	#include <stdlib.h>
	static inline FILE *MemStreamOpen (char **buf, size_t *len)
		{
		*buf = NULL;
		*len = 0;
		return tmpfile ();
		}
	static inline void MemStreamClose (FILE *f, char **buf, size_t *len)
		{
		long n;

		fflush (f);
		n = ftell (f);
		if (n < 0)
			n = 0;
		*buf = (char*) malloc ((size_t) n + 1);
		if (*buf == NULL)
			{
			fprintf (stderr, "Out of memory\n");
			exit (-1);
			}
		rewind (f);
		*len = fread (*buf, 1, (size_t) n, f);
		(*buf)[*len] = '\0';
		fclose (f);
		}
#else
	#define MemStreamOpen(buf, len)			open_memstream (buf, len)
	#define MemStreamClose(f, buf, len)		fclose (f)
#endif

#endif
