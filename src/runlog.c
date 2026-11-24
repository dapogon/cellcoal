#include "cellcoal.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "definitions.h"

/********************* WriteRunRecord  ************************/
/* Writes Results/run.json: everything needed to repeat the run and to know how it was run (version, command line, seed,
   threads, environment variables, genome and read storage modes, times). Failing to write it never stops the run. */
void JsonString (FILE *fp, const char *s)
	{
	fputc ('"', fp);
	for (; s && *s; s++)
		{
		if (*s == '"' || *s == '\\')
			{ fputc ('\\', fp); fputc (*s, fp); }
		else if ((unsigned char) *s < 0x20)
			fprintf (fp, "\\u%04x", (unsigned char) *s);
		else
			fputc (*s, fp);
		}
	fputc ('"', fp);
	}

void WriteRunRecord (double wallSecs, double cpuSecs)
	{
	char		fileArg[2*MAX_NAME+8], when[64];
	FILE		*fp;
	time_t		now = time (NULL);
	struct tm	tmv;
	int			threads = 1, k;
	const char	*envNames[] = {"CELLCOAL_SPARSE", "CELLCOAL_DENSE_MB", "CELLCOAL_STREAM", "CELLCOAL_STREAM_MB", "CELLCOAL_GL_CACHE_MB", "OMP_NUM_THREADS"};
	int			first = YES;

	snprintf (fileArg, sizeof(fileArg), "%s/%s", resultsDir, "run.json");
	if ((fp = fopen (fileArg, "w")) == NULL)
		return;
#ifdef _OPENMP
	threads = omp_get_max_threads ();
#endif
	gmtime_r (&now, &tmv);
	strftime (when, sizeof(when), "%Y-%m-%dT%H:%M:%SZ", &tmv);
	fprintf (fp, "{\n  \"program\": \"%s\",\n  \"version\": ", PROGRAM_NAME);
	JsonString (fp, VERSION_NUMBER);
	fprintf (fp, ",\n  \"finished\": \"%s\",\n  \"command_line\": ", when);
	JsonString (fp, CommandLine ? CommandLine : "");
	{
	/* the options of the run (a valid --config file), with the seed that was used */
	char	*cmd = (char*) calloc (strlen (CommandLine ? CommandLine : "") + 64, 1), *tok, *save = NULL, *copy = strdup (CommandLine ? CommandLine : "");

	for (tok = strtok_r (copy, " ", &save); tok; tok = strtok_r (NULL, " ", &save))
		{
		if (tok[0] == '-' && tok[1] == '#')
			continue;
		strcat (cmd, tok);
		strcat (cmd, " ");
		}
	snprintf (cmd + strlen (cmd), 48, "-#%ld", originalSeed);
	fprintf (fp, ",\n  \"options\": ");
	WriteOptionsJson (fp, cmd, "  ");
	free (cmd);
	free (copy);
	}
	fprintf (fp, ",\n  \"seed\": %ld,\n  \"replicates\": %d,\n  \"cells\": %d,\n  \"sites\": %d,\n", originalSeed, numDataSets, numCells, numSites);
	fprintf (fp, "  \"threads\": %d,\n  \"results_dir\": ", threads);
	JsonString (fp, resultsDir);
	fprintf (fp, ",\n  \"reads_simulated\": %s,\n  \"genome_storage\": \"%s\",\n  \"reads_mode\": \"%s\",\n", (doNGS == YES) ? "true" : "false",
			 (sparseMode == YES) ? "sparse" : "dense", (doNGS == NO) ? "none" : (streamMode == YES) ? ((readsOnly == YES) ? "streamed, reads only" : "streamed") : "stored");
	fprintf (fp, "  \"environment\": {");
	for (k=0; k<(int)(sizeof(envNames)/sizeof(envNames[0])); k++)
		if (getenv (envNames[k]) != NULL)
			{
			fprintf (fp, "%s\n    \"%s\": ", (first == YES) ? "" : ",", envNames[k]);
			JsonString (fp, getenv (envNames[k]));
			first = NO;
			}
	fprintf (fp, "%s},\n  \"compiler\": ", (first == YES) ? "" : "\n  ");
#ifdef __VERSION__
	JsonString (fp, __VERSION__);
#else
	JsonString (fp, "unknown");
#endif
#ifdef _OPENMP
	fprintf (fp, ",\n  \"openmp\": true");
#else
	fprintf (fp, ",\n  \"openmp\": false");
#endif
	fprintf (fp, ",\n  \"wall_seconds\": %.3f,\n  \"cpu_seconds\": %.3f\n}\n", wallSecs, cpuSecs);
	fclose (fp);
	}
