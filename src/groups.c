/*
	CellCoal - grouped options
	The parameters that belong together are given in one option, as a comma-separated list of keys:
	  --ado mean=0.1,block=3000        --output snv-genotypes,vcf-gz        --amplification error=0.05:0.01:0
	Values with several numbers are separated by ':'. A key without a value is a switch.
	Each key stands for one of the single options (-D, -P, --ado-block ...), and the group is turned into those
	options before the command line is read, so a group and the options it stands for give the same run.
	The table below is also what --help prints.
*/
#include "cellcoal.h"

/* kinds of key:  F switch   L single option -X (several numbers separated by ':')   T single option -X with a text value
                  S long option with one value   O several long options, one number each (separated by ':') */
typedef struct { const char *group, *key; char kind; const char *target, *what; } GroupKey;

static const struct { const char *name, *what; } groupInfo[] =
	{
	{"tree",          "population history and tree shape"},
	{"mutation",      "somatic and germline mutation"},
	{"substitution",  "base substitution model, alternative model and rate variation among sites"},
	{"ado",           "allelic dropout (all modes)"},
	{"errors",        "sequencing and genotyping errors"},
	{"amplification", "whole-genome amplification (needs coverage)"},
	{"depth",         "read depth and what varies it (needs coverage)"},
	{"cell-corr",     "correlations between the parameters of the cells (-1 to 1)"},
	{"call-filter",   "which sites are written to the VCF, CATG and PILEUP files (read counts)"},
	{"reads",         "paired-end read simulator (needs coverage; at least bam or fastq)"},
	{"chromosome",    "a real or long chromosome in blocks"},
	{"output",        "results folder and files"},
	{NULL, NULL}
	};

static const GroupKey groupKeys[] =
	{
	{"tree", "growth",        'L', "-g", "exponential growth rate of the population (e.g. growth=1e-5)"},
	{"tree", "birth",         'L', "-K", "cell birth rate (birth-death process; with death)"},
	{"tree", "death",         'L', "-L", "cell death rate"},
	{"tree", "periods",       'L', "-h", "demographic periods: number, then for each its population size at the beginning and at the end and its duration in generations (e.g. periods=2:100:100:40000:1000:2000:20000)"},
	{"tree", "transform",     'L', "-k", "transforming branch length"},
	{"tree", "healthy-tip",   'L', "-q", "healthy tip branch length"},
	{"tree", "lineage-shape", 'L', "-i", "shape of the gamma distribution for rate variation among lineages"},
	{"tree", "rate-switches", 'L', "-N", "number of tree rate switches, then the rate multipliers (e.g. rate-switches=2:10:0.2)"},
	{"tree", "file",          'T', "-T", "user tree file (Newick)"},
	{"tree", "no-outgroup",   'F', "--no-outgroup", "tumor only: no healthy outgroup cell in the outputs"},

	{"mutation", "rate",           'L', "-u", "somatic mutation rate per site per generation"},
	{"mutation", "germline-rate",  'L', "-c", "germline SNP rate"},
	{"mutation", "deletion-rate",  'L', "-d", "deletion rate"},
	{"mutation", "cnloh-rate",     'L', "-H", "CNLOH rate"},
	{"mutation", "loss-rate",      'S', "--seg-loss-rate", "copy-number losses of segments (of a random copy of the haplotype): rate per haplotype and generation"},
	{"mutation", "gain-rate",      'S', "--seg-gain-rate", "copy-number gains of segments (an extra copy of a random copy of the haplotype): rate per haplotype and generation"},
	{"mutation", "segment-length", 'S', "--seg-len", "mean length in sites of the lost and gained segments (default 1000)"},
	{"mutation", "max-copies",     'S', "--seg-max-copies", "maximum number of copies of a site, both haplotypes together (default 4)"},
	{"mutation", "loss-length",    'S', "--seg-loss-len", "mean length in sites of the lost segments (default: segment-length)"},
	{"mutation", "gain-length",    'S', "--seg-gain-len", "mean length in sites of the gained segments (default: segment-length)"},
	{"mutation", "fixed",          'L', "-j", "fixed number of mutations"},
	{"mutation", "signatures",     'L', "-S", "trinucleotide signatures: number, then the signature and its weight for each (e.g. signatures=2:3:0.8:13:0.2)"},
	{"mutation", "signature-file", 'S', "--signature-file", "96-context profiles in COSMIC format used by signatures instead of the 30 built-in ones"},
	{"mutation", "alphabet",       'L', "-b", "0 binary, 1 DNA"},
	{"mutation", "genome",         'T', "-U", "user diploid genome in FASTA (older: use chromosome fasta=...)"},

	{"substitution", "model",            'L', "-m", "alternative mutation model: 0 ISMhap, 1 Mk2, 2 finite-sites DNA (the default is ISM diploid)"},
	{"substitution", "model-proportion", 'L', "-p", "proportion of sites under the alternative model"},
	{"substitution", "model-rate",       'L', "-w", "relative mutation rate of the alternative model"},
	{"substitution", "freqs",            'L', "-f", "nucleotide frequencies A:C:G:T"},
	{"substitution", "titv",             'L', "-t", "transition/transversion ratio"},
	{"substitution", "matrix",           'L', "-r", "mutation matrix, 16 numbers (ACGT x ACGT)"},
	{"substitution", "site-shape",       'L', "-a", "shape of the gamma distribution for rate variation among sites"},

	{"ado", "mean",   'L', "-D", "allelic dropout (fixed)"},
	{"ado", "cells",  'L', "-P", "allelic dropout varying among cells: mean:variance (e.g. cells=0.1:0.001)"},
	{"ado", "sites",  'L', "-Q", "allelic dropout varying among sites: mean:variance"},
	{"ado", "block",  'S', "--ado-block", "dropout in contiguous segments of this mean length in sites (needed in reads mode)"},
	{"ado", "origin", 'S', "--ado-origin", "position of the first site on the chromosome (chromosome pieces)"},
	{"ado", "seed",   'S', "--ado-seed", "seed of the dropout segments (chromosome pieces)"},

	{"errors", "sequencing", 'L', "-E", "sequencing error"},
	{"errors", "matrix",     'L', "-X", "sequencing error matrix, 16 numbers (ACGT x ACGT)"},
	{"errors", "genotyping", 'L', "-G", "genotyping error of the genotypes mode (no coverage): mean:variance (e.g. genotyping=0.01:0.001)"},

	{"amplification", "error",     'L', "-A", "amplification error: mean:variance:two-templates, where two-templates is 0 (all four templates are possible, the default) or 1 (only two templates) (e.g. error=0.05:0.01:0)"},
	{"amplification", "cell-var",  'S', "--amp-cell-var", "variance of the mean amplification error among cells"},
	{"amplification", "signature", 'S', "--amp-signature", "96-context profile of the substitutions of the amplification error (needs error)"},
	{"amplification", "gain-sd",   'S', "--gain-sd", "log-sd of the amplification gain of every haplotype and cell (needs gain-len)"},
	{"amplification", "gain-len",  'S', "--gain-len", "mean length in sites of the segments of constant gain"},

	{"depth", "mean",              'L', "-C", "mean read depth (coverage)"},
	{"depth", "genotypes",         'L', "-V", "overdispersion within a genotype (negative binomial)"},
	{"depth", "cells",             'L', "-J", "overdispersion among cells, shared by all sites"},
	{"depth", "sites",             'L', "-O", "overdispersion among sites, shared by all cells"},
	{"depth", "haploid-reduction", 'L', "-R", "depth factor of a haplotype alone (dropout or deletion of the other)"},
	{"depth", "imbalance",         'L', "-I", "maternal allelic imbalance: mean:variance"},
	{"depth", "doublets",          'L', "-B", "doublet rate per cell: mean:variance"},
	{"depth", "window",            'S', "--coverage", "write the mean depth of every cell over the non-variable sites in windows of this many sites (coverage.tsv)"},

	{"cell-corr", "ado-depth",    'S', "--cov-ado-corr", "correlation between the dropout (ado cells) and the depth (depth cells) of the cells"},
	{"cell-corr", "ado-doublet",  'S', "--ado-doublet-corr", "correlation between the dropout and being a doublet"},
	{"cell-corr", "depth-doublet",'S', "--cov-doublet-corr", "correlation between the depth and being a doublet"},

	{"call-filter", "reads",      'S', "--min-alt-reads", "a site is written only if a cell has a non-reference call with at least this many reads of the allele"},
	{"call-filter", "vaf",        'S', "--min-alt-vaf", "... and at least this fraction of the reads of that cell"},
	{"call-filter", "pooled-vaf", 'L', "-M", "minimum pooled VAF of the sites listed in the VCF"},

	{"reads", "bam",        'F', "--bam", "write BAM files (and indexes) of the amplified cells"},
	{"reads", "fastq",      'F', "--fastq", "write FASTQ files of the amplified cells"},
	{"reads", "length",     'S', "--read-len", "read length (default 150)"},
	{"reads", "insert",     'O', "--insert-mean|--insert-sd", "insert size: mean:sd (default 400:60)"},
	{"reads", "duplicates", 'S', "--dup-mean", "mean extra PCR copies of a molecule (default 0.3)"},
	{"reads", "contig",     'S', "--contig", "name of the contig"},
	{"reads", "qual-sd",    'S', "--qual-sd", "log-sd of the base quality variation among bases (0 to 3)"},
	{"reads", "qual-slope", 'S', "--qual-slope", "decline of the base quality along the read (-6 to 6)"},
	{"reads", "damage", 'S', "--damage", "damage of the library before the PCR, FROM>TO:RATE (e.g. G>T:0.002,C>T:0.001)"},
	{"reads", "seq-context", 'S', "--seq-context", "file with the error multiplier of trinucleotide contexts (e.g. GGC 4)"},
	{"reads", "strand-bias",'S', "--strand-bias", "fraction of the amplification errors seen on one strand only (0 to 1)"},
	{"reads", "bulk-depth", 'S', "--bulk-depth", "also simulate a bulk healthy sample (normal.bam) at this depth"},

	{"chromosome", "blocks",       'S', "--blocks", "number of blocks of the chromosome (the sites option gives the sites per block)"},
	{"chromosome", "fasta",        'S', "--fasta", "reference FASTA (plain or .gz)"},
	{"chromosome", "region",       'S', "--region", "region of the FASTA, chr:start-end"},
	{"chromosome", "germline-vcf", 'S', "--germline-vcf", "germline SNVs to put in the healthy genome"},
	{"chromosome", "region-rates", 'S', "--region-rates", "relative mutation rate of regions (BED: chrom start end multiplier); needs --blocks or --fasta"},

	{"output", "dir",              'T', "-o", "results folder"},
	{"output", "snv-genotypes",    'F', "-1", "SNV genotypes"},
	{"output", "snv-haplotypes",   'F', "-2", "SNV haplotypes"},
	{"output", "full-genotypes",   'F', "-3", "genotypes of all sites"},
	{"output", "full-haplotypes",  'F', "-4", "haplotypes of all sites"},
	{"output", "ancestors",        'F', "-5", "ancestral genotypes and haplotypes"},
	{"output", "trees",            'F', "-6", "trees"},
	{"output", "times",            'F', "-7", "coalescent times"},
	{"output", "catg",             'F', "-8", "read counts in CATG format"},
	{"output", "pileup",           'F', "-9", "read counts in PILEUP format"},
	{"output", "bulk",             'F', "-0", "pooled (bulk) read counts and VAF per site"},
	{"output", "true-haplotypes",  'F', "-Y", "true SNV haplotypes"},
	{"output", "iupac",            'F', "-x", "haplotypes in IUPAC (consensus) code"},
	{"output", "no-vcf",           'F', "-Z", "do not write the VCF (needs catg, pileup or bulk)"},
	{"output", "vcf-gz",           'F', "-z", "write the VCF compressed with gzip"},
	{"output", "bcf",              'F', "--bcf", "write a BCF file instead of a VCF"},
	{"output", "compress-level",   'S', "--compress-level", "zlib level (1 fastest, 9 smallest; default 6) of the BAM, BCF and FASTQ.gz files"},
	{"output", "phase",             'F', "--phase", "show the phase (maternal|paternal copies) in the genotypes with copy-number events"},
	{"output", "summary",          'F', "--summary", "summary.tsv: one row of metrics per replicate"},
	{"output", "replicate-folders",'F', "-v", "one folder per replicate"},
	{"output", "tumor-names",      'F', "-W", "name the cells with tumor nomenclature"},
	{"output", "verbosity",        'L', "-y", "messages on the screen: 0 none, 1 summary, 2 replicates, 3 events"},
	{NULL, NULL, 0, NULL, NULL}
	};

/* the name of the group number i (0 ...), or NULL after the last one */
const char *GroupNameAt (int i)
	{
	int n;

	for (n=0; groupInfo[n].name; n++)
		;
	return (i >= 0 && i < n) ? groupInfo[i].name : NULL;
	}

int IsGroupName (const char *name)
	{
	int g;

	for (g=0; groupInfo[g].name; g++)
		if (strcmp (groupInfo[g].name, name) == 0)
			return 1;
	return 0;
	}

/* configuration-file keys use '_' where the options use '-' */
int IsGroupKeyName (const char *key, char *normalized, size_t size)
	{
	size_t k;

	for (k=0; key[k] && k+1 < size; k++)
		normalized[k] = (key[k] == '_') ? '-' : key[k];
	normalized[k] = '\0';
	return IsGroupName (normalized);
	}

static void GroupFail (const char *group, const char *item, const char *msg)
	{
	int k;

	fprintf (stderr, "\n!!! PARAMETER ERROR: --%s: %s: %s\n", group, msg, item);
	fprintf (stderr, "    keys of --%s:", group);
	for (k=0; groupKeys[k].group; k++)
		if (strcmp (groupKeys[k].group, group) == 0)
			fprintf (stderr, " %s", groupKeys[k].key);
	fprintf (stderr, "\n\n");
	exit (-1);
	}

static char *Join (const char *a, const char *b)
	{
	char *s = (char*) Checked (malloc (strlen (a) + strlen (b) + 1), "Could not allocate the expanded options\n");

	strcpy (s, a);
	strcat (s, b);
	return s;
	}

/* Replaces every "--group list" of in[1..n-1] by the options it stands for; returns the new number of arguments */
int ExpandGroups (char **in, int n, char **out)
	{
	int i, m = 0;

	out[m++] = in[0];
	for (i=1; i<n; i++)
		{
		const char *g = (in[i][0] == '-' && in[i][1] == '-') ? in[i] + 2 : NULL;
		char *list, *item, *save = NULL;

		if (g == NULL || !IsGroupName (g))
			{
			out[m++] = in[i];
			continue;
			}
		if (i + 1 >= n)
			GroupFail (g, "", "a list of keys is needed");
		list = strdup (in[++i]);
		for (item = strtok_r (list, ",", &save); item; item = strtok_r (NULL, ",", &save))
			{
			char *eq = strchr (item, '='), *val = NULL, *vals[40], *p, *sv2 = NULL;
			const GroupKey *k = NULL;
			int e, nv = 0;

			if (eq)
				{ *eq = '\0'; val = eq + 1; }
			for (e=0; groupKeys[e].group; e++)
				if (strcmp (groupKeys[e].group, g) == 0 && strcmp (groupKeys[e].key, item) == 0)
					{ k = &groupKeys[e]; break; }
			if (k == NULL)
				GroupFail (g, item, "unknown key");
			if (k->kind == 'F')
				{
				if (val != NULL && strcmp (val, "1") != 0 && strcmp (val, "true") != 0)
					GroupFail (g, item, "this key is a switch and takes no value");
				out[m++] = strdup (k->target);
				continue;
				}
			if (val == NULL || *val == '\0')
				GroupFail (g, item, "this key needs a value (key=value)");
			if (k->kind == 'T' || k->kind == 'S')
				vals[nv++] = val;
			else
				for (p = strtok_r (val, ":", &sv2); p && nv < 40; p = strtok_r (NULL, ":", &sv2))
					vals[nv++] = p;
			if (k->kind == 'L' || k->kind == 'T')
				{
				out[m++] = Join (k->target, vals[0]);	/* -D0.1 */
				for (e=1; e<nv; e++)
					out[m++] = strdup (vals[e]);
				}
			else if (k->kind == 'S')
				{
				out[m++] = strdup (k->target);
				out[m++] = strdup (vals[0]);
				}
			else	/* 'O': "--a|--b" one number each */
				{
				char *t = strdup (k->target), *tp, *sv3 = NULL;

				e = 0;
				for (tp = strtok_r (t, "|", &sv3); tp && e < nv; tp = strtok_r (NULL, "|", &sv3), e++)
					{
					out[m++] = strdup (tp);
					out[m++] = strdup (vals[e]);
					}
				}
			}
		free (list);
		}
	out[m] = NULL;
	return m;
	}

/* The grouped options, as printed by --help: every key with its meaning and the single option it stands for */
void PrintGroupHelp (FILE *fp)
	{
	int g, k;

	for (g=0; groupInfo[g].name; g++)
		{
		fprintf (fp, "\n--%s KEY[=VALUE],...: %s", groupInfo[g].name, groupInfo[g].what);
		for (k=0; groupKeys[k].group; k++)
			if (strcmp (groupKeys[k].group, groupInfo[g].name) == 0)
				fprintf (fp, "\n    %-18s %s", groupKeys[k].key, groupKeys[k].what);
		}
	}

/* The single options (letters and long ones) that each key stands for, for readers of older commands */
void PrintGroupEquivalents (FILE *fp)
	{
	int gi, k, len;

	for (gi=0; groupInfo[gi].name; gi++)
		{
		fprintf (fp, "\n--%s:", groupInfo[gi].name);
		len = (int) strlen (groupInfo[gi].name) + 3;
		for (k=0; groupKeys[k].group; k++)
			if (strcmp (groupKeys[k].group, groupInfo[gi].name) == 0)
				{
				int add = (int) (strlen (groupKeys[k].key) + strlen (groupKeys[k].target) + 3);

				if (len + add > 110)
					{
					fprintf (fp, "\n    ");
					len = 4;
					}
				fprintf (fp, " %s=%s", groupKeys[k].key, groupKeys[k].target);
				len += add;
				}
		}
	}
