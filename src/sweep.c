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
#ifndef _WIN32
#include <sys/wait.h>
#endif
#ifdef _WIN32
int WindowsExecutablePath (char *buf, size_t n);
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
typedef struct { double sum[SW_MAX_COLS], sumsq[SW_MAX_COLS]; long n, nc[SW_MAX_COLS]; } SweepAcc;	/* n: replicates; nc: values that are not NA, per column */

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
			if (c > 0 && c <= SW_MAX_COLS && strcmp (tok, "NA") != 0)
				{
				acc->sum[c-1] += x;
				acc->sumsq[c-1] += x * x;
				acc->nc[c-1]++;
				}
			c++;
			}
		fprintf (fpAll, "\n");
		acc->n++;
		}
	fclose (fp);
	return 0;
	}

/* a run finished: keep its summary table (this is what --sweep-resume looks for) and remove its folder unless --sweep-keep */
static void SweepFinishRun (const char *dir, int run, int keep)
	{
	char from[1400], to[1400], rundir[1300], logf[1300];
	FILE *in, *out;
	char buf[8192];
	size_t n;

	snprintf (from, sizeof(from), "%s/run_%04d/summary.tsv", dir, run);
	snprintf (to, sizeof(to), "%s/runs/summary_%04d.tsv.tmp", dir, run);
	in = fopen (from, "rb");
	out = fopen (to, "wb");
	if (!in || !out)
		SweepFail ("a run did not write its summary table", from);
	while ((n = fread (buf, 1, sizeof(buf), in)) > 0)
		fwrite (buf, 1, n, out);
	fclose (in);
	fclose (out);
	{
	char final[1400];

	snprintf (final, sizeof(final), "%s/runs/summary_%04d.tsv", dir, run);
	remove (final);
	if (rename (to, final) != 0)
		SweepFail ("cannot keep the summary of a run", final);
	}
	if (keep == NO)
		{
		snprintf (rundir, sizeof(rundir), "%s/run_%04d", dir, run);
		snprintf (logf, sizeof(logf), "%s/run_%04d.log", dir, run);
		RemoveDir (rundir);
		remove (logf);
		}
	}

void RunSweep (int argc, char **argv, const char *spec, int keep, int jobs, int resume)
	{
	SweepParam params[SW_MAX_PARAMS];
	int nparams = 0, i, j, nruns = 1, ncols = 0, seedGiven = NO;
	char *specCopy = strdup (spec), *tok, *save = NULL, dir[1024] = "Sweep", exe[2048], names[SW_MAX_COLS][32];
	char **base = (char**) calloc (argc + 1, sizeof(char*));
	int nbase = 0;
	char **cmds;
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
		if (strcmp (t, "--sweep-jobs") == 0)
			{ i++; continue; }
		if (strcmp (t, "--sweep-keep") == 0 || strcmp (t, "--sweep-resume") == 0 || strcmp (t, "--summary") == 0)
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
	if (!WindowsExecutablePath (exe, sizeof(exe)))
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
	snprintf (path, sizeof(path), "%s/runs", dir);
	MAKE_DIR (path);

	/* the commands of all the runs */
	cmds = (char**) calloc (nruns, sizeof(char*));
	if (!cmds)
		SweepFail ("out of memory", NULL);
	for (i=0; i<nruns; i++)
		{
		int idx[SW_MAX_PARAMS], r = i;
		char *vals[SW_MAX_PARAMS], rundir[1300], logf[1300], *q, *cmd;
		size_t len, cap = 16384;

		cmd = (char*) malloc (cap);
		if (!cmd)
			SweepFail ("out of memory", NULL);
		for (j=nparams-1; j>=0; j--)
			{
			idx[j] = r % params[j].nvals;
			r /= params[j].nvals;
			vals[j] = params[j].vals[idx[j]];
			}
		snprintf (rundir, sizeof(rundir), "%s/run_%04d", dir, i + 1);
		snprintf (logf, sizeof(logf), "%s/run_%04d.log", dir, i + 1);
		q = SweepQuote (exe);
		len = snprintf (cmd, cap, "%s", q);
		free (q);
		for (j=0; j<nbase && len < cap - 1024; j++)
			{
			q = SweepQuote (base[j]);
			len += snprintf (cmd + len, cap - len, " %s", q);
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
						len += snprintf (cmd + len, cap - len, " %s", q);
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
				len += snprintf (cmd + len, cap - len, " %s", q);
				free (q);
				free (arg);
				}
			free (vcopy);
			}
		{
		char *ql = SweepQuote (logf);

		len += snprintf (cmd + len, cap - len, " %c-o%s%c -y0 %c-#%ld%c --summary > %s 2>&1", SHELL_QUOTE_CHAR, rundir, SHELL_QUOTE_CHAR, SHELL_QUOTE_CHAR, seed, SHELL_QUOTE_CHAR, ql);
		free (ql);
		}
		cmds[i] = cmd;
		}

	/* a sweep can only be resumed if it is the same one */
	{
	char specFile[1400], text[20000];
	FILE *fp;
	size_t n = 0;

	snprintf (specFile, sizeof(specFile), "%s/sweep_spec.txt", dir);
	text[0] = '\0';
	for (i=0; i<nruns && n < sizeof(text) - 20000/ (nruns > 4 ? nruns : 4) - 1; i++)
		n += snprintf (text + n, sizeof(text) - n, "%s\n", cmds[i]);
	if (resume == YES && (fp = fopen (specFile, "r")) != NULL)
		{
		char old[20000];
		size_t m = fread (old, 1, sizeof(old) - 1, fp);

		fclose (fp);
		old[m] = '\0';
		if (strcmp (old, text) != 0)
			SweepFail ("--sweep-resume: the folder holds a different sweep (other options, values or seed)", dir);
		}
	else if (resume == YES)
		fprintf (stderr, "No earlier sweep in %s: starting a new one\n", dir);
	if (resume == NO || (fp = fopen (specFile, "r")) == NULL)
		{
		if ((fp = fopen (specFile, "w")) != NULL)
			{
			fputs (text, fp);
			fclose (fp);
			}
		}
	else
		fclose (fp);
	}

	fprintf (stderr, "\nSweep: %d runs (", nruns);
	for (j=0; j<nparams; j++)
		fprintf (stderr, "%s%s x%d", j ? ", " : "", params[j].key, params[j].nvals);
	fprintf (stderr, "), seed %ld, folder %s", seed, dir);
	if (jobs > 1)
		fprintf (stderr, ", %d at a time", jobs);
	fprintf (stderr, "\n");

	/* execute the runs that are not finished (the summary of a finished run is kept in <dir>/runs) */
	{
	int next = 0, running = 0, done = 0;
#ifndef _WIN32
	int *whichRun = (int*) calloc (jobs > 0 ? jobs : 1, sizeof(int));
	pid_t *slotPid = (pid_t*) calloc (jobs > 0 ? jobs : 1, sizeof(pid_t));
	int threadsEach = 0;
	long cores = sysconf (_SC_NPROCESSORS_ONLN);

	if (!whichRun || !slotPid)
		SweepFail ("out of memory", NULL);
	if (jobs > 1 && cores > 0)
		threadsEach = (int) (cores / jobs) > 0 ? (int) (cores / jobs) : 1;
#endif

	for (;;)
		{
		/* start runs while there is a free slot */
		while (next < nruns && running < jobs)
			{
			struct stat st;
			char kept[1400];

			snprintf (kept, sizeof(kept), "%s/runs/summary_%04d.tsv", dir, next + 1);
			if (resume == YES && stat (kept, &st) == 0 && st.st_size > 0)
				{
				fprintf (stderr, "  run %d/%d already done\n", next + 1, nruns);
				next++;
				done++;
				continue;
				}
#ifdef _WIN32
			{
			char wrapped[16400];
			int status;

			/* cmd.exe removes the first and the last quote of a command line that starts with a quote */
			snprintf (wrapped, sizeof(wrapped), "\"%s\"", cmds[next]);
			status = system (wrapped);
			if (status != 0)
				{
				fprintf (stderr, "\n!!! SWEEP ERROR: run %d of %d failed (see %s/run_%04d.log)\n", next + 1, nruns, dir, next + 1);
				exit (-1);
				}
			SweepFinishRun (dir, next + 1, keep);
			next++;
			done++;
			fprintf (stderr, "  run %d/%d done\n", done, nruns);
			}
#else
			{
			pid_t pid = fork ();
			int slot;

			if (pid < 0)
				SweepFail ("cannot start a run (fork failed)", NULL);
			if (pid == 0)
				{
				if (threadsEach > 0 && getenv ("OMP_NUM_THREADS") == NULL)
					{
					char nt[16];

					snprintf (nt, sizeof(nt), "%d", threadsEach);
					setenv ("OMP_NUM_THREADS", nt, 1);
					}
				execl ("/bin/sh", "sh", "-c", cmds[next], (char*) NULL);
				_exit (127);
				}
			for (slot=0; slot<jobs; slot++)
				if (slotPid[slot] == 0)
					break;
			slotPid[slot] = pid;
			whichRun[slot] = next;
			running++;
			next++;
			}
#endif
			}
#ifndef _WIN32
		if (running == 0)
			break;
		{
		int status, slot;
		pid_t pid = wait (&status);

		if (pid < 0)
			break;
		for (slot=0; slot<jobs; slot++)
			if (slotPid[slot] == pid)
				break;
		if (slot == jobs)
			continue;
		slotPid[slot] = 0;
		running--;
		if (!WIFEXITED (status) || WEXITSTATUS (status) != 0)
			{
			fprintf (stderr, "\n!!! SWEEP ERROR: run %d of %d failed (see %s/run_%04d.log); the finished runs are kept: use --sweep-resume to continue\n",
				whichRun[slot] + 1, nruns, dir, whichRun[slot] + 1);
			/* let the other running children end by themselves */
			exit (-1);
			}
		SweepFinishRun (dir, whichRun[slot] + 1, keep);
		done++;
		fprintf (stderr, "  run %d/%d done\n", done, nruns);
		}
#else
		break;
#endif
		}
	}

	/* aggregate the summaries in the order of the grid */
	snprintf (path, sizeof(path), "%s/sweep.tsv", dir);
	if ((fpSweep = fopen (path, "w")) == NULL)
		SweepFail ("cannot write", path);
	snprintf (path, sizeof(path), "%s/sweep_replicates.tsv", dir);
	if ((fpAll = fopen (path, "w")) == NULL)
		SweepFail ("cannot write", path);
	for (i=0; i<nruns; i++)
		{
		int idx[SW_MAX_PARAMS], r = i;
		char *vals[SW_MAX_PARAMS];
		SweepAcc acc;

		for (j=nparams-1; j>=0; j--)
			{
			idx[j] = r % params[j].nvals;
			r /= params[j].nvals;
			vals[j] = params[j].vals[idx[j]];
			}
		memset (&acc, 0, sizeof(acc));
		snprintf (path, sizeof(path), "%s/runs/summary_%04d.tsv", dir, i + 1);
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
			long	n = acc.nc[j];
			double	m = n > 0 ? acc.sum[j] / n : 0.0, var = 0.0;

			if (n == 0)
				{
				fprintf (fpSweep, "\tNA\tNA");	/* the metric is not available in any replicate (e.g. no reads) */
				continue;
				}
			if (n > 1)
				{
				var = (acc.sumsq[j] - n * m * m) / (n - 1);
				if (var < 0)
					var = 0;
				}
			fprintf (fpSweep, "\t%.8g\t%.8g", m, n > 1 ? sqrt (var / n) : 0.0);
			}
		fprintf (fpSweep, "\n");
		}
	fclose (fpSweep);
	fclose (fpAll);
	fprintf (stderr, "Sweep table: %s/sweep.tsv\n\n", dir);
	exit (0);
	}
