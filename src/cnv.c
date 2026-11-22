/*
 *  cnv.c  --  copy-number events simulated together with the SNVs, branch by branch
 *
 *  When --seg-loss-rate or --seg-gain-rate is given, the genomes of the cells are not made in two steps (first all the SNVs, then the
 *  copy-number changes). All the events happen on the branches of the tree and are applied in the order of the tree (from the root to the tips,
 *  and by their position along the branch), so a gain can come before or after an SNV, a loss can remove an SNV that was already there, etc.
 *
 *  Events (all of them are ISM events: a site mutates once in the tree):
 *    SNV     a substitution in one copy of a haplotype of one site.
 *    DEL     a point deletion (-d): one copy of a haplotype of one site is lost.
 *    CNLOH   a point cnLOH (-H): the copies of a haplotype take the allele of a copy of the other haplotype (the copy number does not change).
 *    LOSS    a segment loss (--seg-loss-rate): one copy of a haplotype of a segment is lost.
 *    GAIN    a segment gain (--seg-gain-rate): one copy of a haplotype of a segment is duplicated, up to --seg-max-copies copies per site.
 *
 *  DNA copies. Each haplotype of a cell starts with one copy, the "original" molecule, stored in data[] as always (lineage 0; DELETION when
 *  it is lost). The copies that come from a gain are stored in lists (lineages 1, 2, ...): each one has its sites, its alleles (shared
 *  between the nodes until a mutation changes them) and a lineage number, so the phase of every SNV is known. Every copy mutates at the same
 *  rate: the SNVs of the original molecules are generated at the start, those of a new copy when it is made. The other events (deletions, cnLOH,
 *  losses and gains) have a rate per haplotype and they hit one of its copies at random.
 *
 *  Once the events are applied, the copies of every node are put in the order that the output uses: the first copy that is still there is
 *  in data[] and the others are in the lists, so a haplotype only shows a DELETION when it has no copies at all.
 */
#include "cellcoal.h"

#define MAXC	16							/* maximum number of copies of a haplotype that this code can store */

typedef struct { int s, e, id; Geno *snap; } Seg;		/* extra copy id on the sites [s,e); snap[0] is the allele of site s */
typedef struct { Seg *v; int n, cap; } SegList;
int numPlannedSNV;		/* SNV events planned on the lineages of the two haplotypes */
enum { EV_SNV, EV_DEL, EV_CNLOH, EV_LOSS, EV_GAIN };
typedef struct { int type, hap, lineage, site, len; double u; char done; } Ev;
typedef struct { Ev *v; int n, cap; } EvList;
typedef struct { int type, hap, start, len; TreeNode *branch; } Applied;	/* events with an effect, for the truth file */
typedef struct { int a, b, n, nOther; int ids[MAXC]; } Run;

int				maxCopies = 4;				/* --seg-max-copies: maximum number of copies of a site, both haplotypes together */
int				copyGainOn = 0;				/* 1 when some cell has extra copies */
static SegList	*lin = NULL;				/* extra copies by 2 * node + haplotype */
static EvList	*nodeEv = NULL;				/* pending events by node */
static int		cnNodes = 0, nextLineage = 0;
static Geno		**blocks = NULL;
static int		numBlocks = 0, capBlocks = 0;
static Applied	*applied = NULL;
static int		numApplied = 0, capApplied = 0;
static TreeNode	**order = NULL, **below = NULL;	/* non-root nodes of the tree in preorder; scratch list of a subtree */
static int		numOrder = 0;
static double	orderTotal = 0.0;
static long int	*rseed = NULL;				/* random number seed of the replicate */

/************************************* memory *************************************/
void ResetSegments (void)
	{
	int		k;

	numPlannedSNV = 0;

	if (lin != NULL)
		{
		for (k=0; k<2*cnNodes; k++)
			free (lin[k].v);
		free (lin);
		lin = NULL;
		}
	if (nodeEv != NULL)
		{
		for (k=0; k<cnNodes; k++)
			free (nodeEv[k].v);
		free (nodeEv);
		nodeEv = NULL;
		}
	for (k=0; k<numBlocks; k++)
		free (blocks[k]);
	numBlocks = 0;
	numApplied = 0;
	cnNodes = 0;
	free (order);
	free (below);
	order = below = NULL;
	copyGainOn = NO;
	numSegLoss = numSegGain = 0;
	}

static Geno *NewBlock (int n)
	{
	Geno	*g = (Geno *) Checked (malloc (n > 0 ? n : 1), "Could not allocate the alleles of a copy\n");

	if (numBlocks == capBlocks)
		{
		capBlocks = capBlocks ? 2 * capBlocks : 64;
		blocks = (Geno **) Checked (realloc (blocks, capBlocks * sizeof (Geno *)), "Could not allocate the list of copies\n");
		}
	blocks[numBlocks++] = g;
	return g;
	}

static void SegAdd (SegList *l, int s, int e, int id, Geno *snap)
	{
	if (l->n == l->cap)
		{
		l->cap = l->cap ? 2 * l->cap : 4;
		l->v = (Seg *) Checked (realloc (l->v, l->cap * sizeof (Seg)), "Could not allocate the extra copies of the segments\n");
		}
	l->v[l->n].s = s;
	l->v[l->n].e = e;
	l->v[l->n].id = id;
	l->v[l->n++].snap = snap;
	}

/* removes the sites [a,b) of the copy id (id < 0: of all the extra copies) */
static void SegClip (SegList *l, int id, int a, int b)
	{
	int		k, n = l->n;

	for (k=0; k<n; k++)
		{
		Seg	sg = l->v[k];

		if ((id >= 0 && sg.id != id) || sg.e <= a || sg.s >= b)
			continue;
		if (sg.s < a && sg.e > b)		/* a piece on each side */
			{
			l->v[k].e = a;
			SegAdd (l, b, sg.e, sg.id, sg.snap + (b - sg.s));
			}
		else if (sg.s < a)
			l->v[k].e = a;
		else if (sg.e > b)
			{
			l->v[k].snap += (b - sg.s);
			l->v[k].s = b;
			}
		else
			{
			l->v[k] = l->v[l->n-1];		/* it is all inside [a,b) */
			l->n--;
			k--;
			n--;
			}
		}
	}

static void AddEvent (int node, int type, int hap, int lineage, int site, int len, double u)
	{
	EvList	*l = &nodeEv[node];

	if (l->n == l->cap)
		{
		l->cap = l->cap ? 2 * l->cap : 8;
		l->v = (Ev *) Checked (realloc (l->v, l->cap * sizeof (Ev)), "Could not allocate the list of events\n");
		}
	l->v[l->n].type = type;
	l->v[l->n].hap = hap;
	l->v[l->n].lineage = lineage;
	l->v[l->n].site = site;
	l->v[l->n].len = len;
	l->v[l->n].u = u;
	l->v[l->n++].done = 0;
	}

static void LogApplied (int type, int hap, int start, int len, TreeNode *branch)
	{
	if (numApplied == capApplied)
		{
		capApplied = capApplied ? 2 * capApplied : 16;
		applied = (Applied *) Checked (realloc (applied, capApplied * sizeof (Applied)), "Could not allocate the list of events\n");
		}
	applied[numApplied].type = type;
	applied[numApplied].hap = hap;
	applied[numApplied].start = start;
	applied[numApplied].len = len;
	applied[numApplied++].branch = branch;
	}

/************************************* tree helpers *************************************/
static void Preorder (TreeNode *p, int *n)
	{
	if (p == NULL)
		return;
	if (p->isHealthyRoot != YES)
		order[(*n)++] = p;
	Preorder (p->left, n);
	Preorder (p->right, n);
	}

static void Subtree (TreeNode *p, int *n)
	{
	if (p == NULL)
		return;
	below[(*n)++] = p;
	Subtree (p->left, n);
	Subtree (p->right, n);
	}

/* the branch that a uniform number x in [0,1) falls on, with the branches in the order of the tree and their lengths */
static TreeNode *PickBranch (double x)
	{
	double	cum = 0.0, target = x * orderTotal;
	int		k;

	for (k=0; k<numOrder-1; k++)
		{
		cum += order[k]->branchLength;
		if (cum >= target)
			break;
		}
	return order[k];
	}

/************************************* copies of a site *************************************/
/* the copies that a haplotype of a node has at a site, in the order of their lineage: ids and alleles; returns their number */
static int Present (int node, int h, long site, int *ids, Geno *als)
	{
	const SegList	*l = &lin[2 * node + h];
	Geno			a = GetG (h, node, site);
	int				n = 0, k;

	if (a != DELETION)
		{
		ids[n] = 0;
		als[n++] = a;
		}
	for (k=0; k<l->n; k++)
		if (site >= l->v[k].s && site < l->v[k].e && n < MAXC)
			{
			ids[n] = l->v[k].id;
			als[n++] = l->v[k].snap[site - l->v[k].s];
			}
	return n;
	}

/* the allele of the copy id (0: the original one) of every node of the subtree at the site */
static void SetCopyAllele (int nb, int h, int id, int site, Geno v)
	{
	Geno	*oldSnap = NULL, *newSnap = NULL;
	int		k, j, len;

	if (id == 0)
		{
		for (k=0; k<nb; k++)
			SetG (h, below[k]->label, site, v);
		return;
		}
	for (k=0; k<nb; k++)
		{
		SegList	*l = &lin[2 * below[k]->label + h];

		for (j=0; j<l->n; j++)
			if (l->v[j].id == id && site >= l->v[j].s && site < l->v[j].e)
				{
				if (newSnap == NULL)
					{
					oldSnap = l->v[j].snap;
					len = l->v[j].e - l->v[j].s;
					newSnap = NewBlock (len);
					memcpy (newSnap, oldSnap, len);
					newSnap[site - l->v[j].s] = v;
					}
				if (l->v[j].snap == oldSnap)		/* the nodes of the subtree share the same alleles */
					{
					l->v[j].snap = newSnap;
					break;
					}
				/* a node of the subtree with its own alleles: it only happens if the subtree was not uniform, which cannot be */
				l->v[j].snap[site - l->v[j].s] = v;
				break;
				}
		}
	}

/* the runs of sites of [s,e) in which the haplotype has the same copies (and the other haplotype the same number of copies) */
static int BuildRuns (int node, int h, int s, int e, Run **out)
	{
	Run		*runs = NULL;
	int		nr = 0, cap = 0, j, n, no;
	int		ids[MAXC], idso[MAXC];
	Geno	als[MAXC];

	for (j=s; j<e; j++)
		{
		n = Present (node, h, j, ids, als);
		no = Present (node, 1 - h, j, idso, als);
		if (nr == 0 || runs[nr-1].n != n || runs[nr-1].nOther != no || memcmp (runs[nr-1].ids, ids, n * sizeof (int)) != 0)
			{
			if (nr == cap)
				{
				cap = cap ? 2 * cap : 8;
				runs = (Run *) Checked (realloc (runs, cap * sizeof (Run)), "Could not allocate the runs of sites\n");
				}
			runs[nr].a = j;
			runs[nr].n = n;
			runs[nr].nOther = no;
			memcpy (runs[nr].ids, ids, n * sizeof (int));
			nr++;
			}
		runs[nr-1].b = j + 1;
		}
	*out = runs;
	return nr;
	}

/************************************* events *************************************/
/* the SNVs that a new copy will get: a Poisson on every branch below the gain, at the same rate as any other copy */
static void NewCopyMutations (TreeNode *x, double u0, int h, int id, int a, int b, long int *seed)
	{
	int		k, nb = 0, m, i, y;
	double	eff, u;

	Subtree (x, &nb);
	for (k=0; k<nb; k++)
		{
		y = below[k]->label;
		eff = below[k]->branchLength * (below[k] == x ? 1.0 - u0 : 1.0);
		m = RandomPoisson (eff * (double) (b - a), seed);
		for (i=0; i<m; i++)
			{
			u = (below[k] == x) ? u0 + (1.0 - u0) * RandomUniform (seed) : RandomUniform (seed);
			AddEvent (y, EV_SNV, h, id, a + RandomUniformTo (b - a, seed), 1, u);
			}
		}
	}

static int ApplySNV (TreeNode *x, const Ev *ev, int nb, int ismHap)
	{
	int		ids[MAXC], n, k, j, h = ev->hap, site = ev->site;
	Geno	als[MAXC], a, v;
	double	ran;

	n = Present (x->label, h, site, ids, als);
	for (k=0; k<n && ids[k] != ev->lineage; k++)
		;
	if (k == n)
		return 0;									/* the copy is not there any more */
	if ((ismHap == NO && allSites[site].numMutations > 0) ||
		(ismHap == YES && ((h == MATERNAL && allSites[site].numMutationsMaternal > 0) || (h == PATERNAL && allSites[site].numMutationsPaternal > 0))))
		return 0;									/* the site has already mutated (ISM) */
	a = als[k];
	if (alphabet == DNA)
		{
		ran = RandomUniform (rseed) * cumMij[a][3];
		v = a;
		for (j=0; j<4; j++)
			if (ran <= cumMij[a][j])
				{
				v = j;
				break;
				}
		}
	else if (a == 0)
		v = 1;
	else if (a == 1 && allSites[site].isSNP == YES)	/* germline SNP: the mutation goes back to the reference */
		v = 0;
	else
		return 0;
	SetCopyAllele (nb, h, ev->lineage, site, v);
	if (h == MATERNAL)
		allSites[site].numMutationsMaternal++;
	else
		allSites[site].numMutationsPaternal++;
	allSites[site].numMutations++;
	numMU++;
	return 1;
	}

/* a loss of the segment [s,e) of the haplotype: every run of sites loses one of its copies, chosen at random */
static int ApplyLoss (TreeNode *x, int h, int s, int e, int nb)
	{
	Run		*runs;
	int		nr, r, k, j, id, effect = 0;

	nr = BuildRuns (x->label, h, s, e, &runs);
	for (r=0; r<nr; r++)
		{
		if (runs[r].n == 0)
			continue;
		id = runs[r].ids[RandomUniformTo (runs[r].n, rseed)];
		if (id == 0)
			for (k=0; k<nb; k++)
				for (j=runs[r].a; j<runs[r].b; j++)
					SetG (h, below[k]->label, j, DELETION);
		else
			for (k=0; k<nb; k++)
				SegClip (&lin[2 * below[k]->label + h], id, runs[r].a, runs[r].b);
		effect = 1;
		}
	free (runs);
	return effect;
	}

/* a gain of the segment [s,e): every run of sites gets a new copy of one of its copies, chosen at random */
static int ApplyGain (TreeNode *x, int h, int s, int e, double u, int nb)
	{
	Run		*runs;
	int		nr, r, k, j, i, id, n, ids[MAXC], effect = 0, newId, a, b;
	Geno	als[MAXC], *snap;

	nr = BuildRuns (x->label, h, s, e, &runs);
	for (r=0; r<nr; r++)
		{
		if (runs[r].n == 0 || runs[r].n + runs[r].nOther >= maxCopies || runs[r].n >= MAXC - 1)
			continue;								/* nothing to copy, or too many copies already */
		a = runs[r].a;
		b = runs[r].b;
		id = runs[r].ids[RandomUniformTo (runs[r].n, rseed)];
		snap = NewBlock (b - a);
		for (j=a; j<b; j++)
			{
			n = Present (x->label, h, j, ids, als);
			for (i=0; i<n && ids[i] != id; i++)
				;
			snap[j - a] = als[i];
			}
		newId = ++nextLineage;
		for (k=0; k<nb; k++)
			SegAdd (&lin[2 * below[k]->label + h], a, b, newId, snap);
		copyGainOn = YES;
		NewCopyMutations (x, u, h, newId, a, b, rseed);
		/* NewCopyMutations used the scratch list of the subtree: rebuild it */
		nb = 0;
		Subtree (x, &nb);
		effect = 1;
		}
	free (runs);
	return effect;
	}

/* a cnLOH at the site: the copies of the receiving haplotype take the allele of one copy of the other haplotype */
static int ApplyCNLOH (TreeNode *x, int h, int site, int nb)
	{
	int		idsR[MAXC], idsD[MAXC], nR, nD, k;
	Geno	alsR[MAXC], alsD[MAXC], v;

	nR = Present (x->label, h, site, idsR, alsR);
	nD = Present (x->label, 1 - h, site, idsD, alsD);
	if (nR == 0 || nD == 0)
		return 0;
	v = alsD[RandomUniformTo (nD, rseed)];
	for (k=0; k<nR; k++)
		SetCopyAllele (nb, h, idsR[k], site, v);
	return 1;
	}

/************************************* events of the replicate *************************************/
static void GenerateEvents (int ismHap)
	{
	int		h, i, m, site, trials, len, numMutations;
	double	gens = totalTreeLength / mutationRate, meanLoss, meanGain;
	char	*usedSnv[2], *usedDel[2], *usedLoh[2];

	for (h=0; h<2; h++)
		{
		usedSnv[h] = (char *) Checked (calloc (numSites, 1), "Could not allocate the used sites\n");
		usedDel[h] = (char *) Checked (calloc (numSites, 1), "Could not allocate the used sites\n");
		usedLoh[h] = (char *) Checked (calloc (numSites, 1), "Could not allocate the used sites\n");
		}
	meanLoss = segLossLen > 0 ? segLossLen : segLen;
	meanGain = segGainLen > 0 ? segGainLen : segLen;
	for (h=MATERNAL; h<=PATERNAL; h++)
		{
		/* SNVs of the original copy, as in the infinite sites model: a Poisson (or a fixed number) of mutations, each one on a branch chosen
		   according to its length and on a site that has not mutated */
		if (doSimulateFixedNumMutations == YES)
			{
			if (h == MATERNAL)
				numMutations = numSNVmaternal = RandomBinomial (0.5, numFixedMutations, rseed);
			else
				numMutations = numFixedMutations - numSNVmaternal;
			}
		else
			numMutations = RandomPoisson (totalTreeLength * numSites, rseed);
		if ((ismHap == NO && numISMmutations + numMutations > numSites) || (ismHap == YES && numISMmutations + numMutations > 2 * numSites))
			{
			fprintf (stderr, "\n\nERROR: the infinite sites model (ISM) has been violated. There will be more mutations (%d existing + %d proposed)", numISMmutations, numMutations);
			fprintf (stderr, "\nthan available sites under this model (%d).", ismHap == NO ? numSites : 2 * numSites);
			fprintf (stderr, "\nTry using a smaller mutation rate or effective population size (i.e. smaller theta).\n");
			exit (-1);
			}
		numISMmutations += numMutations;
		trials = 0;
		for (m=0; m<numMutations; m++)
			{
			site = RandomUniformTo (numSites, rseed);
			while ((ismHap == NO && (usedSnv[MATERNAL][site] || usedSnv[PATERNAL][site])) || (ismHap == YES && usedSnv[h][site]))
				{
				site = RandomUniformTo (numSites, rseed);
				if (trials++ > 100 * numSites)
					{
					fprintf (stderr, "\n\nERROR: after %d trials cannot find an unmuted site\n", 100 * numSites);
					exit (-1);
					}
				}
			usedSnv[h][site] = 1;
			{ AddEvent (PickBranch (RandomUniform (rseed))->label, EV_SNV, h, 0, site, 1, RandomUniform (rseed)); numPlannedSNV++; }
			}
		/* point cnLOH and deletions */
		for (i=0; i<2; i++)
			{
			double	rate = (i == 0) ? CNLOHrate : deletionRate;
			char	*used = (i == 0) ? usedLoh[h] : usedDel[h];

			if (rate <= 0)
				continue;
			m = RandomPoisson (numSites * gens * rate, rseed);
			if (m > numSites)
				{
				fprintf (stderr, "\n\nERROR: the haploid infinite sites model (ISM) for %s has been violated: %d events for %d sites.\nTry using a smaller rate.\n",
					i == 0 ? "CNLOH events" : "deletions", m, numSites);
				exit (-1);
				}
			while (m-- > 0)
				{
				site = RandomUniformTo (numSites, rseed);
				for (trials=0; used[site]; trials++)
					{
					site = RandomUniformTo (numSites, rseed);
					if (trials > 100 * numSites)
						{
						fprintf (stderr, "\n\nERROR: after %d trials cannot find a site without events\n", 100 * numSites);
						exit (-1);
						}
					}
				used[site] = 1;
				AddEvent (PickBranch (RandomUniform (rseed))->label, i == 0 ? EV_CNLOH : EV_DEL, h, -1, site, 1, RandomUniform (rseed));
				}
			}
		/* segment losses and gains */
		for (i=0; i<2; i++)
			{
			double	rate = (i == 0) ? lossRate : gainRate, mean = (i == 0) ? meanLoss : meanGain;

			if (rate <= 0)
				continue;
			m = RandomPoisson (gens * rate, rseed);
			while (m-- > 0)
				{
				TreeNode	*br = PickBranch (RandomUniform (rseed));

				site = RandomUniformTo (numSites, rseed);
				len = (mean <= 1.0) ? 1 : 1 + (int) (log (1.0 - RandomUniform (rseed)) / log (1.0 - 1.0 / mean));
				if (site + len > numSites)
					len = numSites - site;
				AddEvent (br->label, i == 0 ? EV_LOSS : EV_GAIN, h, -1, site, len, RandomUniform (rseed));
				if (i == 0)
					numSegLoss++;
				else
					numSegGain++;
				}
			}
		}
	for (h=0; h<2; h++)
		{
		free (usedSnv[h]);
		free (usedDel[h]);
		free (usedLoh[h]);
		}
	}

/* applies the events node by node (from the root to the tips) and, inside a node, by position along the branch */
static void ApplyEvents (int ismHap)
	{
	int		k, i, best, nb, effect;
	TreeNode *x;

	for (k=0; k<numOrder; k++)
		{
		x = order[k];
		for (;;)
			{
			EvList	*l = &nodeEv[x->label];
			Ev		ev;

			best = -1;
			for (i=0; i<l->n; i++)
				if (!l->v[i].done && (best < 0 || l->v[i].u < l->v[best].u))
					best = i;
			if (best < 0)
				break;
			l->v[best].done = 1;
			ev = l->v[best];				/* a copy: NewCopyMutations adds events and can move the list */
			nb = 0;
			Subtree (x, &nb);
			effect = 0;
			switch (ev.type)
				{
				case EV_SNV:
					ApplySNV (x, &ev, nb, ismHap);
					break;
				case EV_DEL:
					if (ApplyLoss (x, ev.hap, ev.site, ev.site + 1, nb))
						{
						numDEL++;
						allSites[ev.site].numDeletions++;
						if (ev.hap == MATERNAL)
							allSites[ev.site].numDeletionsMaternal++;
						else
							allSites[ev.site].numDeletionsPaternal++;
						}
					break;
				case EV_CNLOH:
					if (ApplyCNLOH (x, ev.hap, ev.site, nb))
						{
						numCNLOH++;
						allSites[ev.site].numCNLOH++;
						if (ev.hap == MATERNAL)
							allSites[ev.site].numCNLOHmaternal++;
						else
							allSites[ev.site].numCNLOHpaternal++;
						}
					break;
				case EV_LOSS:
					effect = ApplyLoss (x, ev.hap, ev.site, ev.site + ev.len, nb);
					break;
				case EV_GAIN:
					effect = ApplyGain (x, ev.hap, ev.site, ev.site + ev.len, ev.u, nb);
					break;
				}
			if (effect)
				LogApplied (ev.type, ev.hap, ev.site, ev.len, x);
			}
		}
	}

/* The copies of every node in the order of the output: the first one that is still there goes to data[] and the others stay in the lists */
static int CompareSegs (const void *x, const void *y)
	{
	const Seg	*a = (const Seg *) x, *b = (const Seg *) y;

	return (a->id > b->id) - (a->id < b->id);
	}

static void Normalize (void)
	{
	int		k, h, j, i;

	for (k=0; k<numOrder; k++)
		for (h=0; h<2; h++)
			{
			SegList	*l = &lin[2 * order[k]->label + h];
			int		*pid = NULL, *ps = NULL, *pe = NULL, np = 0, cap = 0, node = order[k]->label;

			if (l->n == 0)
				continue;
			qsort (l->v, l->n, sizeof (Seg), CompareSegs);
			for (i=0; i<l->n; i++)
				for (j=l->v[i].s; j<l->v[i].e; j++)
					if (GetG (h, node, j) == DELETION)
						{
						SetG (h, node, j, l->v[i].snap[j - l->v[i].s]);
						if (np > 0 && pid[np-1] == l->v[i].id && pe[np-1] == j)
							pe[np-1] = j + 1;
						else
							{
							if (np == cap)
								{
								cap = cap ? 2 * cap : 8;
								pid = (int *) Checked (realloc (pid, cap * sizeof (int)), "Could not allocate the promoted copies\n");
								ps = (int *) Checked (realloc (ps, cap * sizeof (int)), "Could not allocate the promoted copies\n");
								pe = (int *) Checked (realloc (pe, cap * sizeof (int)), "Could not allocate the promoted copies\n");
								}
							pid[np] = l->v[i].id;
							ps[np] = j;
							pe[np++] = j + 1;
							}
						}
			for (i=0; i<np; i++)
				SegClip (l, pid[i], ps[i], pe[i]);
			free (pid);
			free (ps);
			free (pe);
			qsort (l->v, l->n, sizeof (Seg), CompareSegs);
			}
	}

/* the extra copies (after the first one, which is in data[]) of the haplotype of a node at a site: their alleles in out[]; returns their number */
int ExtraCopies (int node, int genome, long site, Geno *out)
	{
	const SegList	*l;
	int				k, n = 0;

	if (lin == NULL || node >= cnNodes)
		return 0;
	l = &lin[2 * node + genome];
	for (k=0; k<l->n && n < MAXC; k++)
		if (site >= l->v[k].s && site < l->v[k].e)
			out[n++] = l->v[k].snap[site - l->v[k].s];
	return n;
	}

/* the largest number of extra copies that the haplotype of a node has at some site */
int MaxExtraCopies (int node, int genome)
	{
	const SegList	*l;
	int				mx = 0, j, k, n, s0, e0;

	if (lin == NULL || copyGainOn == NO || node >= cnNodes)
		return 0;
	l = &lin[2 * node + genome];
	for (k=0; k<l->n; k++)
		{
		s0 = l->v[k].s;
		e0 = l->v[k].e;
		for (j=s0; j<e0; j++)		/* the copies are counted site by site: segments of different copies overlap */
			{
			int		c = 0;

			for (n=0; n<l->n; n++)
				if (j >= l->v[n].s && j < l->v[n].e)
					c++;
			if (c > mx)
				mx = c;
			}
		}
	return mx;
	}

/* sites where an extra copy has an allele that differs from the root: they are variable even if the first copies are not */
void MarkExtraCopyVariants (char *isVar)
	{
	int		k, h, i, j;

	if (lin == NULL || copyGainOn == NO)
		return;
	for (k=0; k<=numCells; k++)
		for (h=0; h<2; h++)
			{
			const SegList	*l = &lin[2 * k + h];

			for (i=0; i<l->n; i++)
				for (j=l->v[i].s; j<l->v[i].e; j++)
					if (l->v[i].snap[j - l->v[i].s] != data[h][OUTGROUP_ROOT][j])
						isVar[j] = 1;
			}
	}

/* the whole genome of the replicate: SNVs, point cnLOH and deletions, segment losses and gains, branch by branch */
void SimulateCopyNumberGenomes (TreeNode *root, long int *seed)
	{
	int		ismHap, k;

	if (doGeneticSignatures == YES || (propAltModelSites > 0 && !(propAltModelSites == 1.0 && altModel == ISMhap)))
		{
		fprintf (stderr, "\n\nERROR: the copy-number events (--seg-loss-rate, --seg-gain-rate) are simulated with the infinite sites models (diploid or haploid ISM) only\n");
		fprintf (stderr, "       They cannot be used with genetic signatures, Mk2 or finite sites models.\n");
		exit (-1);
		}
	ismHap = (propAltModelSites == 1.0 && altModel == ISMhap) ? YES : NO;
	ResetSegments ();
	rseed = seed;
	cnNodes = numNodes + 1;
	lin = (SegList *) Checked (calloc (2 * cnNodes, sizeof (SegList)), "Could not allocate the extra copies of the segments\n");
	nodeEv = (EvList *) Checked (calloc (cnNodes, sizeof (EvList)), "Could not allocate the events of the nodes\n");
	order = (TreeNode **) Checked (calloc (cnNodes, sizeof (TreeNode *)), "Could not allocate the list of branches\n");
	below = (TreeNode **) Checked (calloc (cnNodes, sizeof (TreeNode *)), "Could not allocate the list of branches\n");
	numOrder = 0;
	Preorder (root, &numOrder);
	orderTotal = 0.0;
	for (k=0; k<numOrder; k++)
		orderTotal += order[k]->branchLength;
	nextLineage = 0;
	GenerateEvents (ismHap);
	ApplyEvents (ismHap);
	Normalize ();
	}

/************************************* truth file *************************************/
/* the true number of copies of every cell, by segments of sites with the same copies of both haplotypes */
void WriteCnvSegments (int replicate)
	{
	char	file[4*MAX_NAME], name[4*MAX_NAME];
	FILE	*fp;
	int		i, j, m, p, pm = -1, pp = -1, start = 0;
	Geno	ex[MAX_COPIES];

	if (doPrintSeparateReplicates == YES)
		{
		snprintf (file, sizeof(file), "%s/cnv_segments_dir", resultsDir);
		MAKE_DIR (file);
		snprintf (file, sizeof(file), "%s/cnv_segments_dir/cnv_segments.%04d", resultsDir, replicate+1);
		fp = fopen (file, "w");
		}
	else
		{
		snprintf (file, sizeof(file), "%s/cnv_segments", resultsDir);
		fp = fopen (file, replicate == 0 ? "w" : "a");
		}
	if (!fp)
		{
		fprintf (stderr, "Can't open \"%s\"\n", file);
		exit (-1);
		}
	if (doPrintSeparateReplicates == NO)
		fprintf (fp, "[#%d]\n", replicate+1);
	fprintf (fp, "cell\tstart\tend\tmaternal_copies\tpaternal_copies\ttotal_copies\n");
	for (i=0; i<numCells; i++)
		{
		if (doUserTree == NO)
			snprintf (name, sizeof (name), "%s%04d", inCellName, i + 1);
		else
			snprintf (name, sizeof (name), "%s", cellNames[i]);
		pm = pp = -1;
		start = 0;
		for (j=0; j<=numSites; j++)
			{
			if (j < numSites)
				{
				m = (GetG (MATERNAL, i, j) != DELETION) ? 1 + ExtraCopies (i, MATERNAL, j, ex) : 0;
				p = (GetG (PATERNAL, i, j) != DELETION) ? 1 + ExtraCopies (i, PATERNAL, j, ex) : 0;
				}
			else
				m = p = -2;
			if (m != pm || p != pp)
				{
				if (pm >= 0)
					fprintf (fp, "%s\t%d\t%d\t%d\t%d\t%d\n", name, start + 1, j, pm, pp, pm + pp);
				start = j;
				pm = m;
				pp = p;
				}
			}
		}
	fclose (fp);
	}

static void LeafNames (TreeNode *p, char *buf, size_t cap, int *count)
	{
	char	name[4*MAX_NAME];

	if (p == NULL)
		return;
	if (p->left == NULL && p->right == NULL && p->label >= 0 && p->label < numCells)
		{
		if (doUserTree == NO)
			snprintf (name, sizeof (name), "%s%04d", inCellName, p->label + 1);
		else
			snprintf (name, sizeof (name), "%s", cellNames[p->label]);
		if (strlen (buf) + strlen (name) + 2 < cap)
			{
			if (*count > 0)
				strcat (buf, ",");
			strcat (buf, name);
			}
		(*count)++;
		}
	LeafNames (p->left, buf, cap, count);
	LeafNames (p->right, buf, cap, count);
	}

/* Truth of the copy-number events of a replicate: file cnv_events (or cnv_events_dir/cnv_events.NNNN with -v), one row per event that had an
   effect, in the order in which they were applied, with the type, the haplotype, the segment (1-based, inclusive) and the cells that carry it. */
void WriteCnvEvents (int replicate)
	{
	char	file[4*MAX_NAME], *buf;
	FILE	*fp;
	int		k, count;
	size_t	cap = 64 * (size_t) (numCells + 1) + 64;

	if (doPrintSeparateReplicates == YES)
		{
		snprintf (file, sizeof(file), "%s/cnv_events_dir", resultsDir);
		MAKE_DIR (file);
		snprintf (file, sizeof(file), "%s/cnv_events_dir/cnv_events.%04d", resultsDir, replicate+1);
		fp = fopen (file, "w");
		}
	else
		{
		snprintf (file, sizeof(file), "%s/cnv_events", resultsDir);
		fp = fopen (file, replicate == 0 ? "w" : "a");
		}
	if (!fp)
		{
		fprintf (stderr, "Can't open \"%s\"\n", file);
		exit (-1);
		}
	if (doPrintSeparateReplicates == NO)
		fprintf (fp, "[#%d]\n", replicate+1);
	fprintf (fp, "type\thaplotype\tstart\tend\tnum_cells\tcells\n");
	buf = (char *) Checked (malloc (cap), "Could not allocate the list of cells of the events\n");
	for (k=0; k<numApplied; k++)
		{
		static const char *names[] = {"snv", "deletion", "cnloh", "loss", "gain"};

		if (applied[k].type != EV_LOSS && applied[k].type != EV_GAIN)
			continue;
		buf[0] = '\0';
		count = 0;
		LeafNames (applied[k].branch, buf, cap, &count);
		fprintf (fp, "%s\t%s\t%d\t%d\t%d\t%s\n", names[applied[k].type], applied[k].hap == MATERNAL ? "maternal" : "paternal",
			applied[k].start + 1, applied[k].start + applied[k].len, count, buf);
		}
	free (buf);
	fclose (fp);
	}
