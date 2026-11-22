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

#include "readopts.h"

/* generic table of the long options of the read-level simulator (see readopts.h) */
typedef struct { const char *key, *opt; char flag; int kind; void *ptr; double def; } ExtraOpt;
#define X_STRUCT(key, opt, flag, kind, var, def) {#key, opt, flag, kind, &var, def},
static const ExtraOpt extraOpts[] = { READ_OPTS(X_STRUCT) {NULL, NULL, 0, 0, NULL, 0} };
static void CfgError (const char *path, const char *msg, const char *what);
#define X_NAMES(key, opt, flag, kind, var, def) {#key, flag},

static char extraGiven[64];

/* 1 if the long option of the table (--read-len, --dup-mean ...) was given on the command line or in the configuration file */
int ExtraOptionGiven (const char *opt)
	{
	int e;

	for (e=0; extraOpts[e].key; e++)
		if (strcmp (opt, extraOpts[e].opt) == 0)
			return extraGiven[e];
	return 0;
	}

static int IsExtraOpt (const char *a)
	{
	int e;

	for (e=0; extraOpts[e].key; e++)
		if (strcmp (a, extraOpts[e].opt) == 0)
			return 1;
	return 0;
	}

static int ExtraBoolFlag (int flag)
	{
	int e;

	for (e=0; extraOpts[e].key; e++)
		if (extraOpts[e].flag == flag && extraOpts[e].kind == 0)
			return 1;
	return 0;
	}

/* parses the option out[*i] if it is one of the table; returns 1 when it was */
/* options of the table whose value may be 0 (the range is checked later), and the one that may be negative */
static int ExtraZeroAllowed (const char *opt)
	{
	static const char *zero[] = {"--bulk-depth", "--seg-loss-rate", "--seg-gain-rate", "--qual-sd", "--strand-bias", "--dup-mean", NULL};
	int	k;

	for (k=0; zero[k]; k++)
		if (strcmp (opt, zero[k]) == 0)
			return 1;
	return 0;
	}

static int ExtraValueAllowed (const char *opt, double v)
	{
	if (strcmp (opt, "--qual-slope") == 0)
		return 1;					/* negative: the quality improves along the read (range checked later) */
	return ExtraZeroAllowed (opt) ? (v >= 0.0) : (v > 0.0);
	}

static int ParseExtra (char **out, int *i, int n)
	{
	int e;

	for (e=0; extraOpts[e].key; e++)
		if (strcmp (out[*i], extraOpts[e].opt) == 0)
			{
			const ExtraOpt *o = &extraOpts[e];

			extraGiven[e] = 1;
			if (o->kind == 0)
				*(int*) o->ptr = YES;
			else
				{
				char msg[160];

				if (*i + 1 >= n)
					{
					snprintf (msg, sizeof(msg), "%s needs a value", o->opt);
					CfgError ("the command line", msg, NULL);
					}
				(*i)++;
				if (o->kind == 1)
					*(int*) o->ptr = atoi (out[*i]);
				else if (o->kind == 2)
					*(double*) o->ptr = atof (out[*i]);
				else
					*(char**) o->ptr = out[*i];
				if ((o->kind == 1 || o->kind == 2) && !ExtraValueAllowed (o->opt, atof (out[*i])))
					{
					snprintf (msg, sizeof(msg), "%s has to be %s", o->opt, ExtraZeroAllowed (o->opt) ? ">= 0" : "> 0");
					CfgError ("the command line", msg, out[*i]);
					}
				}
			return 1;
			}
	return 0;
	}

static void JsonExtras (FILE *fp, const char *indent, int *emitted)
	{
	int e;

	for (e=0; extraOpts[e].key; e++)
		{
		const ExtraOpt *o = &extraOpts[e];
		int set = 0;

		if (o->kind == 0)
			set = *(int*) o->ptr;
		else if (o->kind == 1)
			set = (*(int*) o->ptr != (int) o->def);
		else if (o->kind == 2)
			set = (*(double*) o->ptr != o->def);
		else
			set = (strlen (*(char**) o->ptr) > 0 && strcmp (*(char**) o->ptr, "chr1") != 0);
		if (!set)
			continue;
		fprintf (fp, "%s\n%s  \"%s\": ", *emitted ? "," : "", indent, o->key);
		if (o->kind == 0)
			fprintf (fp, "true");
		else if (o->kind == 1)
			fprintf (fp, "%d", *(int*) o->ptr);
		else if (o->kind == 2)
			fprintf (fp, "%.10g", *(double*) o->ptr);
		else
			fprintf (fp, "\"%s\"", *(char**) o->ptr);
		(*emitted)++;
		}
	}

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
	{"verbosity",'y'}, {"seed",'#'}, {"fast",'\1'}, {"bcf",'\3'}, {"cov_ado_correlation",'\2'},
	{"cov_doublet_correlation",'\4'}, {"ado_doublet_correlation",'\5'}, {"summary",'\6'},
	{"sweep",'\7'}, {"sweep_keep",'\16'}, {"sweep_jobs",'\17'}, {"sweep_resume",'\20'},
	{"amplification_error_cell_variance",'\21'}, {"no_outgroup",'\22'}, {"min_alt_reads",'\23'}, {"coverage_window",'\24'}, {"min_alt_vaf",'\25'}, {"ado_block",'\26'}, {"ado_origin",'\27'}, {"ado_seed",'\30'}, {"signature_file",'\31'}, {"amp_signature",'\32'}, {"gain_sd",'\33'}, READ_OPTS(X_NAMES) {"gain_len",'\34'}, {NULL,0}
	};

#define CFG_MAX_TOKENS 4096
typedef struct { char *key; char *vals[64]; int nvals; int isBool; int boolValue; int objStart, objEnd; int line, col; } CfgItem;
typedef struct { CfgItem *items; int n, cap; } CfgList;

/* the text of the file being read (a copy that the parsers do not modify), to show the line of an error */
static char *cfgOrig = NULL;

/* prints the line of the file with a caret under the column of the error */
static void CfgShowLine (int line, int col)
	{
	const char *p = cfgOrig;
	int k;

	if (!p || line < 1)
		return;
	for (k=1; k<line && p; k++)
		{
		p = strchr (p, '\n');
		if (p)
			p++;
		}
	if (!p)
		return;
	fprintf (stderr, "  %4d | ", line);
	for (; *p && *p != '\n' && *p != '\r'; p++)
		fputc (*p, stderr);
	fprintf (stderr, "\n       | %*s^\n", col > 1 ? col - 1 : 0, "");
	}

static void CfgErrorAt (const char *path, int line, int col, const char *msg, const char *what, const char *hint)
	{
	if (line > 0)
		fprintf (stderr, "\n!!! CONFIG ERROR in %s, line %d, column %d: %s%s%s", path, line, col, msg, (what && *what) ? ": " : "", (what) ? what : "");
	else
		fprintf (stderr, "\n!!! CONFIG ERROR in %s: %s%s%s", path, msg, (what && *what) ? ": " : "", (what) ? what : "");
	if (hint && *hint)
		fprintf (stderr, "  (did you mean '%s'?)", hint);
	fprintf (stderr, "\n");
	CfgShowLine (line, col);
	fprintf (stderr, "\n");
	exit (-1);
	}

static void CfgError (const char *path, const char *msg, const char *what)
	{
	CfgErrorAt (path, 0, 0, msg, what, NULL);
	}

/* line and column (1-based) of the character at pos in text */
static void PosToLineCol (const char *text, int pos, int *line, int *col)
	{
	int k, ln = 1, start = 0;

	for (k=0; k<pos && text[k]; k++)
		if (text[k] == '\n')
			{ ln++; start = k + 1; }
	*line = ln;
	*col = pos - start + 1;
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
		CfgErrorAt (path, it->line, it->col, "too many values for the key", it->key, NULL);
	it->vals[it->nvals++] = v;
	}

/* ---------- JSON ---------- */
typedef struct { char *s; int pos; const char *path; } JP;

/* error at the current position of the JSON text (cfgOrig has the same offsets as the text being parsed) */
static void JErr (JP *p, const char *msg, const char *what)
	{
	int line, col;

	if (cfgOrig)
		{
		PosToLineCol (cfgOrig, p->pos, &line, &col);
		CfgErrorAt (p->path, line, col, msg, what, NULL);
		}
	CfgError (p->path, msg, what);
	}

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
		{
		char got[2] = {p->s[p->pos], 0};
		JErr (p, p->s[p->pos] ? "a quoted string (a key) was expected, found" : "the file ended while a quoted key was expected", got);
		}
	p->pos++;
	start = q = p->s + p->pos;
	while (p->s[p->pos] != '"')
		{
		if (p->s[p->pos] == '\0')
			JErr (p, "unterminated string", NULL);
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
		JErr (p, "a value was expected for", what);
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
				JErr (p, "nested arrays or objects are not allowed in a value of", it->key);
			else
				CfgAddVal (it, JScalarCopy (p, ",] \t\r\n", it->key), p->path);
			JSkip (p);
			if (p->s[p->pos] == ',')
				{ p->pos++; JSkip (p); }
			else if (p->s[p->pos] != ']')
				JErr (p, "',' or ']' was expected in", it->key);
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
				JErr (p, "unterminated object in", it->key);
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
		JErr (p, "a JSON object was expected", NULL);
	p->pos++;
	JSkip (p);
	while (p->s[p->pos] != '}')
		{
		int keyPos;
		char *key;
		CfgItem *it;

		keyPos = p->pos;
		key = JString (p);
		it = CfgNew (l, key);
		if (cfgOrig)
			PosToLineCol (cfgOrig, keyPos, &it->line, &it->col);

		JSkip (p);
		if (p->s[p->pos] != ':')
			JErr (p, "':' was expected after", key);
		p->pos++;
		JValue (p, it);
		JSkip (p);
		if (p->s[p->pos] == ',')
			{ p->pos++; JSkip (p); }
		else if (p->s[p->pos] != '}')
			JErr (p, "',' or '}' was expected after", key);
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
	int lineNo = 0;

	for (line = s; line; line = next)
		{
		char *c, *colon, *inStr = NULL, *raw = line;
		int q = 0, col;

		lineNo++;

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
		col = (int) (line - raw) + 1;
		if (line[0] == '-' && (line[1] == ' ' || line[1] == '\t'))
			{
			if (!last)
				CfgErrorAt (path, lineNo, col, "a list item ('-') without a key before it", NULL, NULL);
			last->isBool = NO;
			CfgAddVal (last, YUnquote (YTrim (line + 1)), path);
			continue;
			}
		colon = strchr (line, ':');
		if (!colon)
			CfgErrorAt (path, lineNo, col, "a line 'key: value' was expected", line, NULL);
		*colon = '\0';
		last = CfgNew (l, YUnquote (YTrim (line)));
		last->line = lineNo;
		last->col = col;
		colon = YTrim (colon + 1);
		if (*colon == '[')
			{
			char *end = strrchr (colon, ']'), *tok, *save = NULL;

			if (!end)
				CfgErrorAt (path, lineNo, col, "']' is missing in the list of", last->key, NULL);
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
/* edit distance between two keys (case-insensitive; '-' and ' ' count as '_') */
static int EditDistance (const char *a, const char *b)
	{
	int la = (int) strlen (a), lb = (int) strlen (b), i, j, d[64][64];

	if (la > 63 || lb > 63)
		return 99;
	for (i=0; i<=la; i++)
		d[i][0] = i;
	for (j=0; j<=lb; j++)
		d[0][j] = j;
	for (i=1; i<=la; i++)
		for (j=1; j<=lb; j++)
			{
			char x = a[i-1], y = b[j-1];
			int cost, m;

			if (x == '-' || x == ' ') x = '_';
			if (y == '-' || y == ' ') y = '_';
			cost = (tolower ((unsigned char) x) != tolower ((unsigned char) y));
			m = d[i-1][j-1] + cost;
			if (d[i-1][j] + 1 < m) m = d[i-1][j] + 1;
			if (d[i][j-1] + 1 < m) m = d[i][j-1] + 1;
			d[i][j] = m;
			}
	return d[la][lb];
	}

/* the known long name closest to a mistyped key, or NULL if none is close enough */
static const char *SimilarKey (const char *key)
	{
	int i, best = 99, dist, limit = (int) strlen (key) / 3;
	const char *name = NULL;

	if (limit < 2)
		limit = 2;
	for (i=0; optionNames[i].name; i++)
		{
		dist = EditDistance (key, optionNames[i].name);
		if (dist < best)
			{ best = dist; name = optionNames[i].name; }
		}
	return (best <= limit) ? name : NULL;
	}

static int KeyToFlag (const CfgItem *it, const char *path)
	{
	const char *key = it->key;
	int i;

	for (i=0; optionNames[i].name; i++)
		if (strcmp (optionNames[i].name, key) == 0)
			return optionNames[i].flag;
	if (key[0] != '\0' && key[1] == '\0')
		{
		for (i=0; optionNames[i].name; i++)
			if (optionNames[i].flag == key[0])
				return key[0];
		CfgErrorAt (path, it->line, it->col, "unknown option letter", key, NULL);
		}
	CfgErrorAt (path, it->line, it->col, "unknown option", key, SimilarKey (key));
	return 0;
	}

/* options with a value that are written with two dashes: --cov-ado-corr etc. (flag -> option name) */
static const char *LongValueOption (int flag)
	{
	if (flag == '\2') return "--cov-ado-corr";
	if (flag == '\4') return "--cov-doublet-corr";
	if (flag == '\5') return "--ado-doublet-corr";
	if (flag == '\7') return "--sweep";
	if (flag == '\17') return "--sweep-jobs";
	if (flag == '\21') return "--amp-cell-var";
	if (flag == '\23') return "--min-alt-reads";
	if (flag == '\24') return "--coverage";
	if (flag == '\25') return "--min-alt-vaf";
	if (flag == '\26') return "--ado-block";
	if (flag == '\27') return "--ado-origin";
	if (flag == '\30') return "--ado-seed";
	if (flag == '\31') return "--signature-file";
	if (flag == '\32') return "--amp-signature";
	if (flag == '\33') return "--gain-sd";
	if (flag == '\34') return "--gain-len";
	{
	int e;
	for (e=0; extraOpts[e].key; e++)
		if (extraOpts[e].flag == flag)
			return extraOpts[e].opt;
	}
	return NULL;
	}

/* the keys of the configuration file with the option each one stands for (end of --help-all) */
void PrintConfigKeys (FILE *fp)
	{
	int		i, e, col = 0;
	char	what[64];

	fprintf (fp, "\n\n*** KEYS OF THE CONFIGURATION FILE (--config) ***");
	fprintf (fp, "\nEvery option can be written in the file by its letter or by its name; e.g. in YAML:  cells: 20   ado: 0.1   amplification_error: [0.05, 0.01, 0]");
	for (i=0; optionNames[i].name; i++)
		{
		if (optionNames[i].flag >= 33 && optionNames[i].flag < 127)
			snprintf (what, sizeof(what), "-%c", optionNames[i].flag);
		else if (LongValueOption (optionNames[i].flag) != NULL)
			snprintf (what, sizeof(what), "%s", LongValueOption (optionNames[i].flag));
		else
			snprintf (what, sizeof(what), "--%s", optionNames[i].name);
		fprintf (fp, "%s%-34s %-18s", (col % 2 == 0) ? "\n" : "  ", optionNames[i].name, what);
		col++;
		}
	for (e=0; extraOpts[e].key; e++)
		{
		int dup = 0;
		for (i=0; optionNames[i].name; i++)		/* already listed with the options (the long options are in both tables) */
			if (strcmp (optionNames[i].name, extraOpts[e].key) == 0)
				{ dup = 1; break; }
		for (i=0; i<e && !dup; i++)
			if (strcmp (extraOpts[i].key, extraOpts[e].key) == 0)
				dup = 1;
		if (dup)
			continue;
		fprintf (fp, "%s%-34s %-18s", (col % 2 == 0) ? "\n" : "  ", extraOpts[e].key, extraOpts[e].opt);
		col++;
		}
	fprintf (fp, "\n");
	}

static int IsOptionName (const char *key)
	{
	int i;

	for (i=0; optionNames[i].name; i++)
		if (strcmp (optionNames[i].name, key) == 0)
			return 1;
	return 0;
	}

static char *Join2 (const char *a, const char *b)
	{
	char *s = (char*) malloc (strlen (a) + strlen (b) + 1);

	strcpy (s, a);
	strcat (s, b);
	return s;
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
			CfgErrorAt (path, it->line, it->col, "the key has no value", it->key, NULL);
			}
		{
		char norm[64];

		/* a group written as a string (ado: mean=0.1,block=3000) or as a list (output: [snv-genotypes, vcf-gz]) */
		if (it->nvals >= 1 && IsGroupKeyName (it->key, norm, sizeof(norm)) && (strchr (it->vals[0], '=') != NULL || !IsOptionName (it->key)))
			{
			size_t len = 1;
			char *joined;
			int w;

			for (w=0; w<it->nvals; w++)
				len += strlen (it->vals[w]) + 1;
			joined = (char*) malloc (len);
			joined[0] = '\0';
			for (w=0; w<it->nvals; w++)
				{
				if (w > 0)
					strcat (joined, ",");
				strcat (joined, it->vals[w]);
				}
			out[(*n)++] = Join2 ("--", norm);
			out[(*n)++] = joined;
			continue;
			}
		}
		flag = KeyToFlag (it, path);
		if (it->isBool == YES && it->nvals == 0)
			{
			if (it->boolValue == NO)
				continue;
			tok = (char*) malloc (16);
			if (flag == '\1')
				snprintf (tok, 16, "--fast");
			else if (flag == '\3')
				snprintf (tok, 16, "--bcf");
			else if (flag == '\6')
				snprintf (tok, 16, "--summary");
			else if (flag == '\22')
				snprintf (tok, 16, "--no-outgroup");
			else if (flag == '\16')
				snprintf (tok, 16, "--sweep-keep");
			else if (flag == '\20')
				snprintf (tok, 16, "--sweep-resume");
			else if (ExtraBoolFlag (flag))
				snprintf (tok, 16, "%s", LongValueOption (flag));
			else if (LongValueOption (flag))
				continue;
			else
				snprintf (tok, 16, "-%c", flag);
			out[(*n)++] = tok;
			continue;
			}
		for (v=0; v<it->nvals; v++)
			{
			if (*n >= CFG_MAX_TOKENS - 1)
				CfgError (path, "too many options", NULL);
			if (v == 0 && LongValueOption (flag))
				{
				out[(*n)++] = strdup (LongValueOption (flag));
				tok = strdup (it->vals[0]);
				}
			else if (v == 0)
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
	/* in text mode (Windows) fewer bytes than the file size can be read: CRLF becomes LF */
	if (!buf)
		{
		fprintf (stderr, "\n!!! CONFIG ERROR: cannot read %s\n\n", path);
		exit (-1);
		}
	size = (long) fread (buf, 1, size, fp);
	buf[size] = '\0';
	fclose (fp);
	return buf;
	}

/* The options of a configuration file, as command-line tokens (without the program name) */
static int ConfigTokens (const char *path, char **out)
	{
	char *text = ReadWholeFile (path), *p = text;

	cfgOrig = strdup (text);
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

/* The command-line option ("-C", "--cov-ado-corr") of a key of --sweep: an option letter or a long name (dashes or underscores) */
const char *SweepOptionForName (const char *name)
	{
	static char buf[80];
	char norm[64];
	int i;

	snprintf (norm, sizeof(norm), "%s", name);
	for (i=0; norm[i]; i++)
		if (norm[i] == '-')
			norm[i] = '_';
	for (i=0; optionNames[i].name; i++)
		{
		char f = optionNames[i].flag;

		if (f == '\1' || f == '\3' || f == '\6' || f == '\7' || f == '\16' || f == '\22')
			continue;	/* switches: nothing to sweep */
		if (strcmp (optionNames[i].name, norm) == 0 || (norm[1] == '\0' && norm[0] == f) ||
			(f == '\2' && strcmp (norm, "cov_ado_corr") == 0) || (f == '\4' && strcmp (norm, "cov_doublet_corr") == 0) ||
			(f == '\5' && strcmp (norm, "ado_doublet_corr") == 0) || (f == '\21' && strcmp (norm, "amp_cell_var") == 0))
			{
			const char *lo = LongValueOption (f);

			if (lo)
				return lo;
			snprintf (buf, sizeof(buf), "-%c", f);
			return buf;
			}
		}
	return NULL;
	}

/* Replaces "--config FILE" in the arguments by the options of the file (the program name stays first) */
char **ExpandConfigArgs (int *argc, char **argv)
	{
	int i, k, n = 0, found = NO;
	char **out;

	for (i=1; i<*argc; i++)
		{
		if (strcmp (argv[i], "--version") == 0)
			{
			printf ("CellCoal %s\n", VERSION_NUMBER);
			exit (0);
			}
		if (strcmp (argv[i], "--help") == 0)
			{
			PrintHelp (stdout);
			exit (0);
			}
		if (strcmp (argv[i], "--help-all") == 0)	/* the same as --help */
			{
			PrintHelp (stdout);
			exit (0);
			}
		}
	for (i=1; i<*argc; i++)
		if (strcmp (argv[i], "--check") == 0 || strcmp (argv[i], "--config") == 0 || strcmp (argv[i], "--fast") == 0 || strcmp (argv[i], "--no-outgroup") == 0 || strcmp (argv[i], "--cov-ado-corr") == 0 || strcmp (argv[i], "--bcf") == 0 ||
			strcmp (argv[i], "--cov-doublet-corr") == 0 || strcmp (argv[i], "--ado-doublet-corr") == 0 || strcmp (argv[i], "--summary") == 0 ||
			strcmp (argv[i], "--sweep") == 0 || strcmp (argv[i], "--sweep-keep") == 0 ||
			strcmp (argv[i], "--sweep-jobs") == 0 || strcmp (argv[i], "--sweep-resume") == 0 || strcmp (argv[i], "--amp-cell-var") == 0 || strcmp (argv[i], "--min-alt-reads") == 0 || strcmp (argv[i], "--coverage") == 0 || strcmp (argv[i], "--min-alt-vaf") == 0 || strcmp (argv[i], "--ado-block") == 0 || strcmp (argv[i], "--ado-origin") == 0 || strcmp (argv[i], "--ado-seed") == 0 || strcmp (argv[i], "--signature-file") == 0 || strcmp (argv[i], "--amp-signature") == 0 || strcmp (argv[i], "--gain-sd") == 0 || strcmp (argv[i], "--gain-len") == 0 || IsExtraOpt (argv[i]))
			found = YES;
	for (i=1; i<*argc; i++)
		if (argv[i][0] == '-' && argv[i][1] == '-' && IsGroupName (argv[i] + 2))
			found = YES;
	if (found == NO)
		return argv;
	out = (char**) Checked (calloc(CFG_MAX_TOKENS + *argc + 1, sizeof(char*)), "Could not allocate the expanded arguments\n");
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
	/* grouped options (--ado mean=0.1,block=3000 ...) are turned into the single options */
	{
	char **expanded = (char**) Checked (calloc(CFG_MAX_TOKENS * 4 + *argc + 1, sizeof(char*)), "Could not allocate the expanded arguments\n");

	n = ExpandGroups (out, n, expanded);
	out = expanded;
	}
	/* --sweep SPEC: the program runs itself once per combination and exits */
	{
	const char *sweepSpec = NULL;
	int sweepKeep = NO, sweepResume = NO, sweepJobs = 1;

	for (i=1; i<n; i++)
		{
		if (strcmp (out[i], "--sweep") == 0)
			{
			if (i + 1 >= n)
				CfgError ("the command line", "--sweep needs a specification such as \"C=2,5,10 E=0.001,0.01\"", NULL);
			sweepSpec = out[++i];
			}
		else if (strcmp (out[i], "--sweep-keep") == 0)
			sweepKeep = YES;
		else if (strcmp (out[i], "--sweep-resume") == 0)
			sweepResume = YES;
		else if (strcmp (out[i], "--sweep-jobs") == 0)
			{
			if (i + 1 >= n || atoi (out[i+1]) < 1 || atoi (out[i+1]) > 256)
				CfgError ("the command line", "--sweep-jobs needs a number of runs at a time (1 to 256)", NULL);
			sweepJobs = atoi (out[++i]);
			}
		}
	if (sweepSpec)
		{
		out[n] = NULL;
		RunSweep (n, out, sweepSpec, sweepKeep, sweepJobs, sweepResume);
		}
	}
	/* --fast (from the command line or a configuration file) switches on the fast mode and is not passed on */
	for (i=1, k=1; i<n; i++)
		{
		if (strcmp (out[i], "--cov-ado-corr") == 0 || strcmp (out[i], "--cov-doublet-corr") == 0 || strcmp (out[i], "--ado-doublet-corr") == 0)
			{
			const char *name = out[i];
			double *target = (strcmp (name, "--cov-ado-corr") == 0) ? &covAdoCorr : (strcmp (name, "--cov-doublet-corr") == 0) ? &covDoubletCorr : &adoDoubletCorr;
			char msg[96];

			if (i + 1 >= n)
				{
				snprintf (msg, sizeof(msg), "%s needs a value between -1 and 1", name);
				CfgError ("the command line", msg, NULL);
				}
			*target = atof (out[++i]);
			if (!(*target >= -1.0 && *target <= 1.0))
				{
				snprintf (msg, sizeof(msg), "%s has to be between -1 and 1", name);
				CfgError ("the command line", msg, out[i]);
				}
			}
		else if (strcmp (out[i], "--amp-cell-var") == 0)
			{
			if (i + 1 >= n)
				CfgError ("the command line", "--amp-cell-var needs a variance (> 0)", NULL);
			ampCellVar = atof (out[++i]);
			if (!(ampCellVar > 0.0))
				CfgError ("the command line", "--amp-cell-var has to be > 0", out[i]);
			}
		else if (strcmp (out[i], "--min-alt-vaf") == 0)
			{
			if (i + 1 >= n)
				CfgError ("the command line", "--min-alt-vaf needs a fraction (0 to 1)", NULL);
			minAltVaf = atof (out[++i]);
			if (minAltVaf < 0.0 || minAltVaf > 1.0)
				CfgError ("the command line", "--min-alt-vaf has to be between 0 and 1", out[i]);
			}
		else if (strcmp (out[i], "--ado-block") == 0)
			{
			if (i + 1 >= n)
				CfgError ("the command line", "--ado-block needs a mean segment length in sites (> 0)", NULL);
			adoBlockLen = atof (out[++i]);
			if (!(adoBlockLen > 0.0))
				CfgError ("the command line", "--ado-block has to be > 0", out[i]);
			}
		else if (strcmp (out[i], "--signature-file") == 0)
			{
			if (i + 1 >= n)
				CfgError ("the command line", "--signature-file needs a file name", NULL);
			signatureFile = out[++i];
			}
		else if (strcmp (out[i], "--amp-signature") == 0)
			{
			if (i + 1 >= n)
				CfgError ("the command line", "--amp-signature needs a file name", NULL);
			ampSignatureFile = out[++i];
			}
		else if (strcmp (out[i], "--gain-sd") == 0)
			{
			if (i + 1 >= n)
				CfgError ("the command line", "--gain-sd needs a standard deviation (> 0)", NULL);
			gainSd = atof (out[++i]);
			if (!(gainSd > 0.0))
				CfgError ("the command line", "--gain-sd has to be > 0", out[i]);
			}
		else if (strcmp (out[i], "--gain-len") == 0)
			{
			if (i + 1 >= n)
				CfgError ("the command line", "--gain-len needs a mean segment length in sites (> 0)", NULL);
			gainLen = atof (out[++i]);
			if (!(gainLen > 0.0))
				CfgError ("the command line", "--gain-len has to be > 0", out[i]);
			}
		else if (strcmp (out[i], "--ado-origin") == 0)
			{
			if (i + 1 >= n)
				CfgError ("the command line", "--ado-origin needs a position (>= 0)", NULL);
			adoOrigin = atoll (out[++i]);
			if (adoOrigin < 0)
				CfgError ("the command line", "--ado-origin has to be >= 0", out[i]);
			}
		else if (strcmp (out[i], "--ado-seed") == 0)
			{
			if (i + 1 >= n)
				CfgError ("the command line", "--ado-seed needs a positive integer", NULL);
			adoSeedOpt = atol (out[++i]);
			if (adoSeedOpt < 1)
				CfgError ("the command line", "--ado-seed has to be >= 1", out[i]);
			}
		else if (strcmp (out[i], "--coverage") == 0)
			{
			if (i + 1 >= n)
				CfgError ("the command line", "--coverage needs a window size in sites (>= 1)", NULL);
			coverageWindow = atoi (out[++i]);
			if (coverageWindow < 1)
				CfgError ("the command line", "--coverage has to be >= 1", out[i]);
			}
		else if (strcmp (out[i], "--min-alt-reads") == 0)
			{
			if (i + 1 >= n)
				CfgError ("the command line", "--min-alt-reads needs a number of reads (>= 0)", NULL);
			minAltReads = atoi (out[++i]);
			if (minAltReads < 0)
				CfgError ("the command line", "--min-alt-reads has to be >= 0", out[i]);
			}
		else if (ParseExtra (out, &i, n))
			;
		else if (strcmp (out[i], "--fast") == 0)
			fastMode = YES;
		else if (strcmp (out[i], "--no-outgroup") == 0)
			numOutgroup = 0;
		else if (strcmp (out[i], "--bcf") == 0)
			doBCF = YES;
		else if (strcmp (out[i], "--summary") == 0)
			doSummary = YES;
		else if (strcmp (out[i], "--check") == 0)
			checkOnly = YES;
		else if (strcmp (out[i], "--sweep-keep") == 0 || strcmp (out[i], "--sweep-resume") == 0)
			;
		else if (strcmp (out[i], "--sweep-jobs") == 0)
			i++;	/* only meaningful together with --sweep */
		else
			out[k++] = out[i];
		}
	n = k;
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
	/* the long options are taken out of the command line before it is recorded: write them from their variables, with the names of the configuration files */
	{
	int emitted = nflags;

#define JSON_EXTRA(cond, name, text, value) \
	if (cond) \
		{ \
		fprintf (fp, "%s\n%s  \"" name "\": " text, emitted ? "," : "", indent, value); \
		emitted++; \
		}
	JSON_EXTRA (fastMode == YES, "fast", "%s", "true") JSON_EXTRA (numOutgroup == 0, "no_outgroup", "%s", "true")
	JSON_EXTRA (doBCF == YES, "bcf", "%s", "true")
	JSON_EXTRA (doSummary == YES, "summary", "%s", "true")
	JSON_EXTRA (covAdoCorr != 0.0, "cov_ado_correlation", "%.10g", covAdoCorr)
	JSON_EXTRA (covDoubletCorr != 0.0, "cov_doublet_correlation", "%.10g", covDoubletCorr)
	JSON_EXTRA (adoDoubletCorr != 0.0, "ado_doublet_correlation", "%.10g", adoDoubletCorr)
	JSON_EXTRA (gainSd > 0.0, "gain_sd", "%.10g", gainSd) JSON_EXTRA (gainLen > 0.0, "gain_len", "%.10g", gainLen) JSON_EXTRA (ampSignatureFile != NULL, "amp_signature", "\"%s\"", ampSignatureFile) JSON_EXTRA (signatureFile != NULL, "signature_file", "\"%s\"", signatureFile) JSON_EXTRA (adoBlockLen > 0.0, "ado_block", "%.10g", adoBlockLen) JSON_EXTRA (adoOrigin > 0, "ado_origin", "%lld", adoOrigin) JSON_EXTRA (adoSeedOpt > 0, "ado_seed", "%ld", adoSeedOpt) JSON_EXTRA (minAltVaf > 0.0, "min_alt_vaf", "%.10g", minAltVaf) JSON_EXTRA (coverageWindow > 0, "coverage_window", "%d", coverageWindow) JSON_EXTRA (minAltReads > 0, "min_alt_reads", "%d", minAltReads) JSON_EXTRA (ampCellVar > 0.0, "amplification_error_cell_variance", "%.10g", ampCellVar)
#undef JSON_EXTRA
	JsonExtras (fp, indent, &emitted);
	}
	fprintf (fp, "\n%s}", indent);
	free (copy);
	}
