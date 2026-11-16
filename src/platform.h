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
#else
	#define MAKE_DIR(path)		mkdir (path, S_IRWXU)
	#define LSTAT(path, st)		lstat (path, st)
	#define SHELL_QUOTE_CHAR	'\''
#endif

#endif
