/*
	CellCoal - parameter sweeps (--sweep)
	  cellcoal <options> --sweep "C=2,5,10 E=0.001,0.01" [-oSweep] [--sweep-keep]
	runs the simulation once for every combination (grid) of the listed values, with the same seed in all of them, asks each
	run for its summary table (--summary) and collects the mean and the standard error over replicates of every metric in
	<dir>/sweep.tsv (one row per combination). All the replicates of all the runs are in <dir>/sweep_replicates.tsv.
	A term is KEY=VALUES where KEY is an option letter (C, E, s, ...) or a long name (coverage, cells, cov-ado-corr ...) and
	VALUES a comma-separated list; a value can also be a range LOW..HIGH..STEP; options with several numbers
	(e.g. -A 0.05 0.01 0) are written with colons: A=0.05:0.01:0,0.1:0.02:0. The run folders are deleted unless --sweep-keep.
*/
#include "cellcoal.h"
#include <unistd.h>
#ifdef _WIN32
#include <windows.h>
#endif

#define SW_MAX_PARAMS	8
#define SW_MAX_VALUES	256
#define SW_MAX_RUNS		10000
#define SW_MAX_COLS		32

typedef struct
	{
	char	key[64];			/* option as given */
	char	opt[64];			/* "-C" or "--cov-ado-corr" */
	int		isLong;
	char	*vals[SW_MAX_VALUES];
	int		nvals;
	} SweepParam;


static void SweepFail (const char *msg, const char *what)
	{
	fprintf (stderr, "\n!!! SWEEP ERROR: %s%s%s%s\n\n", msg, what ? " (\"" : "", what ? what : "", what ? "\")" : "");
	exit (-1);
	}

/* quotes an argument for the shell: '...' on POSIX, "..." on Windows (cmd.exe) */
static char *SweepQuote (const char *s)
	{
	size_t i;
	char *q = (char*) malloc (4 * strlen (s) + 3), *p;

	if (!q)
		SweepFail ("out of memory", NULL);
	p = q;
	*p++ = SHELL_QUOTE_CHAR;
	for (i=0; s[i]; i++)
		{
#ifndef _WIN32
		if (s[i] == '\'')
			{ memcpy (p, "'\\''", 4); p += 4; continue; }
#endif
		*p++ = s[i];
		}
	*p++ = SHELL_QUOTE_CHAR;
	*p = '\0';
	return q;
	}

/* adds the values of a range LOW..HIGH..STEP or the single value to the parameter */
static void SweepAddValue (SweepParam *sp, const char *v)
	{
	const char *d = strstr (v, "..");

	if (d)
		{
		char *end, part[64], *d2;
		double lo, hi, step, x;
		int i;

		if (d == v || (size_t)(d - v) >= sizeof(part))
			SweepFail ("bad range (use LOW..HIGH..STEP)", v);
		memcpy (part, v, d - v);
		part[d - v] = '\0';
		lo = strtod (part, &end);
		if (*end != '\0')
			SweepFail ("bad range (use LOW..HIGH..STEP)", v);
		d2 = strstr (d + 2, "..");
		if (!d2 || (size_t)(d2 - (d + 2)) >= sizeof(part))
			SweepFail ("bad range (use LOW..HIGH..STEP)", v);
		memcpy (part, d + 2, d2 - (d + 2));
		part[d2 - (d + 2)] = '\0';
		hi = strtod (part, &end);
		if (*end != '\0' || part[0] == '\0')
			SweepFail ("bad range (use LOW..HIGH..STEP)", v);
		step = strtod (d2 + 2, &end);
		if (*end != '\0' || !(step > 0) || hi < lo)
			SweepFail ("bad range (use LOW..HIGH..STEP with STEP > 0 and HIGH >= LOW)", v);
		for (i=0, x=lo; x <= hi + 1e-9 * step; i++, x = lo + (i) * step)
			{
			char buf[64];

			if (sp->nvals >= SW_MAX_VALUES)
				SweepFail ("too many values in a parameter", sp->key);
			snprintf (buf, sizeof(buf), "%.10g", x);
			sp->vals[sp->nvals++] = strdup (buf);
			}
		return;
		}
	if (sp->nvals >= SW_MAX_VALUES)
		SweepFail ("too many values in a parameter", sp->key);
	sp->vals[sp->nvals++] = strdup (v);
	}

static int SweepIsFlagToken (const char *t)
	{
	if (t[0] != '-' || t[1] == '\0')
		return NO;
	if (isdigit ((unsigned char) t[1]))
		return (t[2] == '\0');	/* -1 ... -9, -0 are options; -0.5 is a number */
	return (t[1] != '.');
	}

/* statistics of a column */
typedef struct { double sum[SW_MAX_COLS], sumsq[SW_MAX_COLS]; long n; } SweepAcc;

static int SweepReadSummary (const char *path, int run, char **paramVals, int nparams, SweepAcc *acc, int *ncols, char names[][32], FILE *fpAll, int writeAllHeader)
	{
	FILE *fp = fopen (path, "r");
	char line[4096];
	int first = YES, i;

	if (!fp)
		return -1;
	while (fgets (line, sizeof(line), fp))
		{
		char *tok, *save = NULL, copy[4096];
		int c = 0;

		strcpy (copy, line);
		if (first)
			{
			first = NO;
			if (writeAllHeader == YES)
				{
				for (tok = strtok_r (copy, "\t\r\n", &save); tok; tok = strtok_r (NULL, "\t\r\n", &save))
					{
					if (c > 0 && c <= SW_MAX_COLS)
						snprintf (names[c-1], 32, "%s", tok);
					c++;
					}
				*ncols = c - 1;
				fprintf (fpAll, "run");
				for (i=0; i<nparams; i++)
					fprintf (fpAll, "\tparam%d", i+1);
				fprintf (fpAll, "\treplicate");
				for (i=0; i<*ncols; i++)
					fprintf (fpAll, "\t%s", names[i]);
				fprintf (fpAll, "\n");
				}
			continue;
			}
		fprintf (fpAll, "%d", run);
		for (i=0; i<nparams; i++)
			fprintf (fpAll, "\t%s", paramVals[i]);
		for (tok = strtok_r (copy, "\t\r\n", &save); tok; tok = strtok_r (NULL, "\t\r\n", &save))
			{
			double x = atof (tok);

			fprintf (fpAll, "\t%s", tok);
			if (c > 0 && c <= SW_MAX_COLS)
				{
				acc->sum[c-1] += x;
				acc->sumsq[c-1] += x * x;
				}
			c++;
			}
		fprintf (fpAll, "\n");
		acc->n++;
		}
	fclose (fp);
	return 0;
	}

void RunSweep (int argc, char **argv, const char *spec, int keep)
	{
	SweepParam params[SW_MAX_PARAMS];
	int nparams = 0, i, j, nruns = 1, ncols = 0, seedGiven = NO;
	char *specCopy = strdup (spec), *tok, *save = NULL, dir[1024] = "Sweep", exe[2048], names[SW_MAX_COLS][32];
	char **base = (char**) calloc (argc + 1, sizeof(char*));
	int nbase = 0;
	FILE *fpSweep, *fpAll;
	long seed;
	char path[1400];

	if (!specCopy || !base)
		SweepFail ("out of memory", NULL);
	memset (params, 0, sizeof(params));

	/* parse the specification */
	for (tok = strtok_r (specCopy, " ;\t", &save); tok; tok = strtok_r (NULL, " ;\t", &save))
		{
		char *eq = strchr (tok, '='), *vs, *vsave = NULL;
		SweepParam *sp;
		const char *lopt;

		if (!eq || eq == tok || eq[1] == '\0')
			SweepFail ("each term must be KEY=VALUE[,VALUE...]", tok);
		if (nparams >= SW_MAX_PARAMS)
			SweepFail ("too many parameters in the sweep", NULL);
		sp = &params[nparams++];
		*eq = '\0';
		snprintf (sp->key, sizeof(sp->key), "%s", tok);
		lopt = SweepOptionForName (tok);
		if (!lopt)
			SweepFail ("unknown option in the sweep", tok);
		snprintf (sp->opt, sizeof(sp->opt), "%s", lopt);
		sp->isLong = (lopt[1] == '-');
		if (sp->opt[1] == 'o' || sp->opt[1] == '#' || sp->opt[1] == 'F' || strcmp (sp->opt, "-n") == 0)
			SweepFail ("this option cannot be swept", tok);
		for (vs = strtok_r (eq + 1, ",", &vsave); vs; vs = strtok_r (NULL, ",", &vsave))
			SweepAddValue (sp, vs);
		if (sp->nvals == 0)
			SweepFail ("no values for", sp->key);
		nruns *= sp->nvals;
		if (nruns > SW_MAX_RUNS)
			SweepFail ("the grid has too many combinations (maximum 10000)", NULL);
		}
	if (nparams == 0)
		SweepFail ("--sweep needs at least one KEY=VALUES term", NULL);

	/* base command line: the program name, without --sweep and its swept and fixed (-o, -#) options */
	for (i=1; i<argc; i++)
		{
		const char *t = argv[i];
		int skip = NO;

		if (strcmp (t, "--sweep") == 0)
			{ i++; continue; }
		if (strcmp (t, "--sweep-keep") == 0 || strcmp (t, "--summary") == 0)
			continue;
		if (SweepIsFlagToken (t))
			{
			if (t[1] == 'o' && t[2] != '\0')
				{
				snprintf (dir, sizeof(dir), "%s", t + 2);
				skip = YES;
				}
			else if (t[1] == '#')
				{
				seedGiven = YES;
				seed = atol (t + 2);
				skip = YES;
				}
			else
				for (j=0; j<nparams; j++)
					if (strcmp (t, params[j].opt) == 0 || (!params[j].isLong && t[1] == params[j].opt[1] && t[1] != '-' ))
						skip = YES;
			if (skip)
				{
				while (i + 1 < argc && !SweepIsFlagToken (argv[i+1]))
					i++;
				continue;
				}
			}
		base[nbase++] = argv[i];
		}

	if (seedGiven == NO)
		seed = (long) (time (NULL) % 1000000000L) + 1;
	if (seed <= 0)
		seed = 1;

	{
#ifdef _WIN32
	DWORD len = GetModuleFileNameA (NULL, exe, (DWORD) sizeof(exe));

	if (len == 0 || len >= sizeof(exe))
		snprintf (exe, sizeof(exe), "%s", argv[0]);
#else
	ssize_t len = readlink ("/proc/self/exe", exe, sizeof(exe) - 1);

	if (len > 0)
		exe[len] = '\0';
	else
		snprintf (exe, sizeof(exe), "%s", argv[0]);	/* macOS: the path used to start the program */
#endif
	}

	MAKE_DIR (dir);
	snprintf (path, sizeof(path), "%s/sweep.tsv", dir);
	if ((fpSweep = fopen (path, "w")) == NULL)
		SweepFail ("cannot write", path);
	snprintf (path, sizeof(path), "%s/sweep_replicates.tsv", dir);
	if ((fpAll = fopen (path, "w")) == NULL)
		SweepFail ("cannot write", path);

	fprintf (stderr, "\nSweep: %d runs (", nruns);
	for (j=0; j<nparams; j++)
		fprintf (stderr, "%s%s x%d", j ? ", " : "", params[j].key, params[j].nvals);
	fprintf (stderr, "), seed %ld, folder %s\n", seed, dir);

	for (i=0; i<nruns; i++)
		{
		int idx[SW_MAX_PARAMS], r = i, status;
		char *vals[SW_MAX_PARAMS], cmd[16384], rundir[1300], logf[1300], *q;
		SweepAcc acc;
		size_t len;

		for (j=nparams-1; j>=0; j--)
			{
			idx[j] = r % params[j].nvals;
			r /= params[j].nvals;
			vals[j] = params[j].vals[idx[j]];
			}
		snprintf (rundir, sizeof(rundir), "%s/run_%04d", dir, i + 1);
		snprintf (logf, sizeof(logf), "%s/run_%04d.log", dir, i + 1);
		q = SweepQuote (exe);
		len = snprintf (cmd, sizeof(cmd), "%s", q);
		free (q);
		for (j=0; j<nbase && len < sizeof(cmd) - 512; j++)
			{
			q = SweepQuote (base[j]);
			len += snprintf (cmd + len, sizeof(cmd) - len, " %s", q);
			free (q);
			}
		for (j=0; j<nparams; j++)
			{
			char *vcopy = strdup (vals[j]), *p, *sv = NULL;
			int first = YES;

			for (p = strtok_r (vcopy, ":", &sv); p; p = strtok_r (NULL, ":", &sv))
				{
				char *arg;
				size_t al = strlen (params[j].opt) + strlen (p) + 2;

				arg = (char*) malloc (al);
				if (first == YES)
					{
					if (params[j].isLong)
						{
						q = SweepQuote (params[j].opt);
						len += snprintf (cmd + len, sizeof(cmd) - len, " %s", q);
						free (q);
						snprintf (arg, al, "%s", p);
						}
					else
						snprintf (arg, al, "%s%s", params[j].opt, p);
					}
				else
					snprintf (arg, al, "%s", p);
				first = NO;
				q = SweepQuote (arg);
				len += snprintf (cmd + len, sizeof(cmd) - len, " %s", q);
				free (q);
				free (arg);
				}
			free (vcopy);
			}
		{
		char tail[2800], *qd, *ql;

		qd = SweepQuote (rundir);
		ql = SweepQuote (logf);
		snprintf (tail, sizeof(tail), " %c-o%s%c -y0 %c-#%ld%c --summary > %s 2>&1", SHELL_QUOTE_CHAR, rundir, SHELL_QUOTE_CHAR, SHELL_QUOTE_CHAR, seed, SHELL_QUOTE_CHAR, ql);
		/* the directory is quoted inside -o */
		(void) qd;
		len += snprintf (cmd + len, sizeof(cmd) - len, "%s", tail);
		free (qd);
		free (ql);
		}
#ifdef _WIN32
		/* cmd.exe removes the first and the last quote of a command line that starts with a quote */
		{
		char wrapped[sizeof(cmd) + 4];

		snprintf (wrapped, sizeof(wrapped), "\"%s\"", cmd);
		status = system (wrapped);
		}
#else
		status = system (cmd);
#endif
		if (status != 0)
			{
			fprintf (stderr, "\n!!! SWEEP ERROR: run %d of %d failed (see %s)\n", i + 1, nruns, logf);
			exit (-1);
			}
		memset (&acc, 0, sizeof(acc));
		snprintf (path, sizeof(path), "%s/summary.tsv", rundir);
		if (SweepReadSummary (path, i + 1, vals, nparams, &acc, &ncols, names, fpAll, i == 0 ? YES : NO) != 0)
			SweepFail ("a run did not write its summary table", path);
		if (i == 0)
			{
			fprintf (fpSweep, "run");
			for (j=0; j<nparams; j++)
				fprintf (fpSweep, "\t%s", params[j].key);
			fprintf (fpSweep, "\tn");
			for (j=0; j<ncols; j++)
				fprintf (fpSweep, "\t%s\t%s_se", names[j], names[j]);
			fprintf (fpSweep, "\n");
			}
		fprintf (fpSweep, "%d", i + 1);
		for (j=0; j<nparams; j++)
			fprintf (fpSweep, "\t%s", vals[j]);
		fprintf (fpSweep, "\t%ld", acc.n);
		for (j=0; j<ncols; j++)
			{
			double m = acc.n > 0 ? acc.sum[j] / acc.n : 0.0, var = 0.0;

			if (acc.n > 1)
				{
				var = (acc.sumsq[j] - acc.n * m * m) / (acc.n - 1);
				if (var < 0)
					var = 0;
				}
			fprintf (fpSweep, "\t%.8g\t%.8g", m, acc.n > 1 ? sqrt (var / acc.n) : 0.0);
			}
		fprintf (fpSweep, "\n");
		fflush (fpSweep);
		if (keep == NO)
			{
			RemoveDir (rundir);
			remove (logf);
			}
		fprintf (stderr, "  run %d/%d done\n", i + 1, nruns);
		}
	fclose (fpSweep);
	fclose (fpAll);
	fprintf (stderr, "Sweep table: %s/sweep.tsv\n\n", dir);
	exit (0);
	}
