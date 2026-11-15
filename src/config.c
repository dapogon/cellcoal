/*
	CellCoal - configuration files (--config)
	A run can be described in a JSON or YAML file instead of a long list of options. The keys are the option
	letters (e.g. "s": 20, "A": [0.05, 0.01, 0], "1": true) or the long names of the table below
	(e.g. "cells": 20, "amplification_error": [0.05, 0.01, 0], "print_snv_genotypes": true).
	Also run.json (Results/run.json), whose "options" object lists the options of that run, can be used.
	YAML: only "key: value", "key: [a, b]" and "- item" lines are understood (no nesting), with # comments.
	The file is turned into the command-line options that the program already understands, so
	  cellcoal --config run.yaml -n5          is the same as          cellcoal -s20 ... -n5
	(options after --config are read later, so they override the file).
*/
#include "cellcoal.h"

static const struct { const char *name; char flag; } optionNames[] =
	{
	{"parameter_file",'F'}, {"replicates",'n'}, {"cells",'s'}, {"sites",'l'}, {"population_size",'e'},
	{"growth_rate",'g'}, {"birth_rate",'K'}, {"death_rate",'L'}, {"demographic_periods",'h'},
	{"transforming_branch_length",'k'}, {"healthy_tip_branch_length",'q'}, {"lineage_rate_shape",'i'},
	{"tree_rate_switches",'N'}, {"alphabet",'b'}, {"germline_snp_rate",'c'}, {"mutation_rate",'u'},
	{"deletion_rate",'d'}, {"cnloh_rate",'H'}, {"fixed_mutations",'j'}, {"signatures",'S'},
	{"alt_model",'m'}, {"alt_model_proportion",'p'}, {"alt_model_rate",'w'}, {"base_frequencies",'f'},
	{"ti_tv",'t'}, {"site_rate_shape",'a'}, {"mutation_matrix",'r'}, {"ado",'D'}, {"ado_cells",'P'},
	{"ado_sites",'Q'}, {"allelic_imbalance",'I'}, {"amplification_error",'A'}, {"doublet_rate",'B'},
	{"genotyping_error",'G'}, {"coverage",'C'}, {"coverage_dispersion",'V'}, {"coverage_dispersion_sites",'O'},
	{"coverage_dispersion_cells",'J'}, {"haploid_coverage_reduction",'R'}, {"sequencing_error",'E'},
	{"error_matrix",'X'}, {"print_snv_genotypes",'1'}, {"print_snv_haplotypes",'2'},
	{"print_full_genotypes",'3'}, {"print_full_haplotypes",'4'}, {"print_ancestors",'5'}, {"print_trees",'6'},
	{"print_times",'7'}, {"print_catg",'8'}, {"print_pileup",'9'}, {"print_bulk",'0'}, {"no_vcf",'Z'},
	{"gzip_vcf",'z'}, {"min_vaf",'M'}, {"print_true_haplotypes",'Y'}, {"replicate_folders",'v'},
	{"iupac",'x'}, {"results_folder",'o'}, {"tree_file",'T'}, {"genome_file",'U'}, {"tumor_names",'W'},
	{"verbosity",'y'}, {"seed",'#'}, {NULL,0}
	};

#define CFG_MAX_TOKENS 4096
typedef struct { char *key; char *vals[64]; int nvals; int isBool; int boolValue; int objStart, objEnd; } CfgItem;
typedef struct { CfgItem *items; int n, cap; } CfgList;

static void CfgError (const char *path, const char *msg, const char *what)
	{
	fprintf (stderr, "\n!!! CONFIG ERROR in %s: %s%s%s\n\n", path, msg, (what && *what) ? ": " : "", (what) ? what : "");
	exit (-1);
	}

static CfgItem *CfgNew (CfgList *l, char *key)
	{
	CfgItem *it;

	if (l->n == l->cap)
		{
		l->cap = l->cap ? 2 * l->cap : 64;
		l->items = (CfgItem*) realloc (l->items, l->cap * sizeof(CfgItem));
		if (!l->items)
			{
			fprintf (stderr, "Could not allocate the configuration items\n");
			exit (-1);
			}
		}
	it = &l->items[l->n++];
	memset (it, 0, sizeof(*it));
	it->key = key;
	return it;
	}

static void CfgAddVal (CfgItem *it, char *v, const char *path)
	{
	if (it->nvals >= 64)
		CfgError (path, "too many values for the key", it->key);
	it->vals[it->nvals++] = v;
	}

/* ---------- JSON ---------- */
typedef struct { char *s; int pos; const char *path; } JP;

static void JSkip (JP *p)
	{
	while (p->s[p->pos] == ' ' || p->s[p->pos] == '\t' || p->s[p->pos] == '\n' || p->s[p->pos] == '\r')
		p->pos++;
	}

/* parses a JSON string at the current position (the opening quote), in place; returns the unescaped text */
static char *JString (JP *p)
	{
	char *start, *q;

	if (p->s[p->pos] != '"')
		CfgError (p->path, "a string was expected", p->s + p->pos);
	p->pos++;
	start = q = p->s + p->pos;
	while (p->s[p->pos] != '"')
		{
		if (p->s[p->pos] == '\0')
			CfgError (p->path, "unterminated string", NULL);
		if (p->s[p->pos] == '\\' && p->s[p->pos+1] != '\0')
			{
			char c = p->s[++p->pos];
			*q++ = (c == 'n') ? '\n' : (c == 't') ? '\t' : c;
			p->pos++;
			}
		else
			*q++ = p->s[p->pos++];
		}
	*q = '\0';
	p->pos++;
	return start;
	}

/* copy of the number or literal at the current position, up to one of the delimiters */
static char *JScalarCopy (JP *p, const char *delims, const char *what)
	{
	int s0 = p->pos;
	char *copy;

	while (p->s[p->pos] && strchr (delims, p->s[p->pos]) == NULL)
		p->pos++;
	if (p->pos == s0)
		CfgError (p->path, "a value was expected for", what);
	copy = (char*) malloc (p->pos - s0 + 1);
	memcpy (copy, p->s + s0, p->pos - s0);
	copy[p->pos - s0] = '\0';
	return copy;
	}

static void JObject (JP *p, CfgList *l);

static void JValue (JP *p, CfgItem *it)
	{
	JSkip (p);
	if (p->s[p->pos] == '"')
		CfgAddVal (it, JString (p), p->path);
	else if (p->s[p->pos] == '[')
		{
		p->pos++;
		JSkip (p);
		while (p->s[p->pos] != ']')
			{
			if (p->s[p->pos] == '"')
				CfgAddVal (it, JString (p), p->path);
			else if (p->s[p->pos] == '[' || p->s[p->pos] == '{')
				CfgError (p->path, "nested arrays or objects are not allowed in a value of", it->key);
			else
				CfgAddVal (it, JScalarCopy (p, ",] \t\r\n", it->key), p->path);
			JSkip (p);
			if (p->s[p->pos] == ',')
				{ p->pos++; JSkip (p); }
			else if (p->s[p->pos] != ']')
				CfgError (p->path, "',' or ']' was expected in", it->key);
			}
		p->pos++;
		}
	else if (p->s[p->pos] == '{')
		{
		/* remember where the object is (it is parsed again if it is the "options" object) */
		int depth = 0, inStr = 0;
		it->objStart = p->pos;
		for (;;)
			{
			char c = p->s[p->pos];
			if (c == '\0')
				CfgError (p->path, "unterminated object in", it->key);
			if (inStr)
				{ if (c == '\\') p->pos++; else if (c == '"') inStr = 0; }
			else if (c == '"')
				inStr = 1;
			else if (c == '{')
				depth++;
			else if (c == '}' && --depth == 0)
				{ p->pos++; break; }
			p->pos++;
			}
		it->objEnd = p->pos;
		}
	else
		{
		char *v = JScalarCopy (p, ",}] \t\r\n", it->key);

		if (strcmp (v, "true") == 0)
			it->isBool = YES, it->boolValue = YES;
		else if (strcmp (v, "false") == 0 || strcmp (v, "null") == 0)
			it->isBool = YES, it->boolValue = NO;
		else
			CfgAddVal (it, v, p->path);
		}
	}

static void JObject (JP *p, CfgList *l)
	{
	JSkip (p);
	if (p->s[p->pos] != '{')
		CfgError (p->path, "a JSON object was expected", NULL);
	p->pos++;
	JSkip (p);
	while (p->s[p->pos] != '}')
		{
		char *key = JString (p);
		CfgItem *it = CfgNew (l, key);

		JSkip (p);
		if (p->s[p->pos] != ':')
			CfgError (p->path, "':' was expected after", key);
		p->pos++;
		JValue (p, it);
		JSkip (p);
		if (p->s[p->pos] == ',')
			{ p->pos++; JSkip (p); }
		else if (p->s[p->pos] != '}')
			CfgError (p->path, "',' or '}' was expected after", key);
		}
	p->pos++;
	}

/* ---------- YAML (flat subset) ---------- */
static char *YTrim (char *s)
	{
	char *e;

	while (*s == ' ' || *s == '\t')
		s++;
	e = s + strlen (s);
	while (e > s && (e[-1] == ' ' || e[-1] == '\t' || e[-1] == '\r'))
		*--e = '\0';
	return s;
}

static char *YUnquote (char *v)
	{
	size_t n = strlen (v);

	if (n >= 2 && ((v[0] == '"' && v[n-1] == '"') || (v[0] == '\'' && v[n-1] == '\'')))
		{
		v[n-1] = '\0';
		return v + 1;
		}
	return v;
	}

static void YAddScalar (CfgItem *it, char *v, const char *path)
	{
	v = YUnquote (YTrim (v));
	if (strcmp (v, "true") == 0 || strcmp (v, "false") == 0 || strcmp (v, "null") == 0)
		{ it->isBool = YES; it->boolValue = (strcmp (v, "true") == 0); }
	else if (*v)
		CfgAddVal (it, v, path);
	}

static void YamlParse (char *s, CfgList *l, const char *path)
	{
	char *line, *next;
	CfgItem *last = NULL;

	for (line = s; line; line = next)
		{
		char *c, *colon, *inStr = NULL;
		int q = 0;

		next = strchr (line, '\n');
		if (next)
			*next++ = '\0';
		/* comments: a # outside quotes */
		for (c = line; *c; c++)
			{
			if (*c == '"' || *c == '\'')
				q = (q == *c) ? 0 : (q ? q : *c);
			else if (*c == '#' && q == 0 && (c == line || c[-1] == ' ' || c[-1] == '\t'))
				{ *c = '\0'; break; }
			}
		(void) inStr;
		line = YTrim (line);
		if (*line == '\0' || strcmp (line, "---") == 0)
			continue;
		if (line[0] == '-' && (line[1] == ' ' || line[1] == '\t'))
			{
			if (!last)
				CfgError (path, "a list item without a key", line);
			last->isBool = NO;
			CfgAddVal (last, YUnquote (YTrim (line + 1)), path);
			continue;
			}
		colon = strchr (line, ':');
		if (!colon)
			CfgError (path, "a line 'key: value' was expected", line);
		*colon = '\0';
		last = CfgNew (l, YUnquote (YTrim (line)));
		colon = YTrim (colon + 1);
		if (*colon == '[')
			{
			char *end = strrchr (colon, ']'), *tok, *save = NULL;

			if (!end)
				CfgError (path, "']' is missing", last->key);
			*end = '\0';
			for (tok = strtok_r (colon + 1, ",", &save); tok; tok = strtok_r (NULL, ",", &save))
				{
				tok = YUnquote (YTrim (tok));
				if (*tok)
					CfgAddVal (last, tok, path);
				}
			}
		else if (*colon)
			YAddScalar (last, colon, path);
		}
	}

/* ---------- from items to options ---------- */
static int KeyToFlag (const char *key, const char *path)
	{
	int i;

	if (key[0] != '\0' && key[1] == '\0')
		return key[0];
	for (i=0; optionNames[i].name; i++)
		if (strcmp (optionNames[i].name, key) == 0)
			return optionNames[i].flag;
	CfgError (path, "unknown option", key);
	return 0;
	}

static void ItemsToArgs (CfgList *l, const char *path, char **out, int *n, int onlyKnown)
	{
	int k, v;

	for (k=0; k<l->n; k++)
		{
		CfgItem *it = &l->items[k];
		char *tok;
		int flag;

		if (it->objStart != it->objEnd || (it->nvals == 0 && it->isBool == NO))
			{
			if (onlyKnown)
				continue;
			CfgError (path, "the key has no value", it->key);
			}
		flag = KeyToFlag (it->key, path);
		if (it->isBool == YES && it->nvals == 0)
			{
			if (it->boolValue == NO)
				continue;
			tok = (char*) malloc (3);
			snprintf (tok, 3, "-%c", flag);
			out[(*n)++] = tok;
			continue;
			}
		for (v=0; v<it->nvals; v++)
			{
			if (*n >= CFG_MAX_TOKENS - 1)
				CfgError (path, "too many options", NULL);
			if (v == 0)
				{
				size_t len = strlen (it->vals[0]) + 3;
				tok = (char*) malloc (len);
				snprintf (tok, len, "-%c%s", flag, it->vals[0]);
				}
			else
				tok = strdup (it->vals[v]);
			out[(*n)++] = tok;
			}
		}
	}

static char *ReadWholeFile (const char *path)
	{
	FILE *fp = fopen (path, "r");
	long size;
	char *buf;

	if (!fp)
		{
		fprintf (stderr, "\n!!! CONFIG ERROR: cannot open the configuration file %s\n\n", path);
		exit (-1);
		}
	fseek (fp, 0, SEEK_END);
	size = ftell (fp);
	rewind (fp);
	buf = (char*) malloc (size + 1);
	if (!buf || (long) fread (buf, 1, size, fp) != size)
		{
		fprintf (stderr, "\n!!! CONFIG ERROR: cannot read %s\n\n", path);
		exit (-1);
		}
	buf[size] = '\0';
	fclose (fp);
	return buf;
	}

/* The options of a configuration file, as command-line tokens (without the program name) */
static int ConfigTokens (const char *path, char **out)
	{
	char *text = ReadWholeFile (path), *p = text;
	CfgList l = {NULL, 0, 0};
	int n = 0, i, hasOptions = NO;

	while (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r' || (unsigned char) *p == 0xEF || (unsigned char) *p == 0xBB || (unsigned char) *p == 0xBF)
		p++;
	if (*p == '{')
		{
		JP jp;

		jp.s = text; jp.pos = (int) (p - text); jp.path = path;
		JObject (&jp, &l);
		for (i=0; i<l.n; i++)
			if (strcmp (l.items[i].key, "options") == 0 && l.items[i].objStart != l.items[i].objEnd)
				{
				CfgList sub = {NULL, 0, 0};
				JP sp;

				sp.s = text; sp.pos = l.items[i].objStart; sp.path = path;
				JObject (&sp, &sub);
				ItemsToArgs (&sub, path, out, &n, NO);
				hasOptions = YES;
				break;
				}
		if (hasOptions == NO)
			ItemsToArgs (&l, path, out, &n, NO);
		}
	else
		{
		YamlParse (p, &l, path);
		ItemsToArgs (&l, path, out, &n, NO);
		}
	return n;
	}

/* Replaces "--config FILE" in the arguments by the options of the file (the program name stays first) */
char **ExpandConfigArgs (int *argc, char **argv)
	{
	int i, n = 0, found = NO;
	char **out;

	for (i=1; i<*argc; i++)
		if (strcmp (argv[i], "--config") == 0)
			found = YES;
	if (found == NO)
		return argv;
	out = (char**) calloc (CFG_MAX_TOKENS + *argc + 1, sizeof(char*));
	if (!out)
		{
		fprintf (stderr, "Could not allocate the expanded arguments\n");
		exit (-1);
		}
	out[n++] = argv[0];
	for (i=1; i<*argc; i++)
		{
		if (strcmp (argv[i], "--config") == 0)
			{
			if (i + 1 >= *argc)
				{
				fprintf (stderr, "\n!!! CONFIG ERROR: --config needs a file name\n\n");
				exit (-1);
				}
			n += ConfigTokens (argv[++i], out + n);
			}
		else
			out[n++] = argv[i];
		}
	out[n] = NULL;
	*argc = n;
	return out;
	}

/* Writes the options of a command line (as recorded in CommandLine: " prog -s20 -A0.05 0.01 0 ...") as a JSON object whose
   keys are the option letters; it can be read back with --config. */
static void JsonValue (FILE *fp, const char *v)
	{
	int isNum = (*v != '\0' && strspn (v, "0123456789+-.eE") == strlen (v));

	if (isNum)
		fprintf (fp, "%s", v);
	else
		{
		fputc ('"', fp);
		for (; *v; v++)
			{
			if (*v == '"' || *v == '\\')
				fputc ('\\', fp);
			fputc (*v, fp);
			}
		fputc ('"', fp);
		}
	}

void WriteOptionsJson (FILE *fp, const char *commandLine, const char *indent)
	{
	char *copy = strdup (commandLine ? commandLine : ""), *tok, *save = NULL;
	char *flagTok[512], *vals[512][32];
	int nvals[512], nflags = 0, k, v;

	for (tok = strtok_r (copy, " ", &save); tok; tok = strtok_r (NULL, " ", &save))
		{
		if (tok[0] == '-' && tok[1] != '\0')
			{
			if (nflags == 512)
				break;
			flagTok[nflags] = tok;
			nvals[nflags] = 0;
			if (tok[2] != '\0')
				vals[nflags][nvals[nflags]++] = tok + 2;
			nflags++;
			}
		else if (nflags > 0 && nvals[nflags-1] < 32)
			vals[nflags-1][nvals[nflags-1]++] = tok;
		}
	fprintf (fp, "{");
	for (k=0; k<nflags; k++)
		{
		fprintf (fp, "%s\n%s  \"%c\": ", k ? "," : "", indent, flagTok[k][1]);
		if (nvals[k] == 0)
			fprintf (fp, "true");
		else if (nvals[k] == 1)
			JsonValue (fp, vals[k][0]);
		else
			{
			fprintf (fp, "[");
			for (v=0; v<nvals[k]; v++)
				{
				fprintf (fp, "%s", v ? ", " : "");
				JsonValue (fp, vals[k][v]);
				}
			fprintf (fp, "]");
			}
		}
	fprintf (fp, "\n%s}", indent);
	free (copy);
	}
