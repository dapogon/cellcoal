/*
	CellCoal - path of the running executable on Windows (kept apart because windows.h defines names that clash with
	the single-letter constants of definitions.h). On other systems this file is empty.
*/
#ifdef _WIN32
#include <windows.h>
#include <stddef.h>

/* writes the path of the executable to buf (size n); returns 1 on success, 0 otherwise */
int WindowsExecutablePath (char *buf, size_t n)
	{
	DWORD len = GetModuleFileNameA (NULL, buf, (DWORD) n);

	return (len > 0 && len < n) ? 1 : 0;
	}
#endif
