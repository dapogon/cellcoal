//
//  cellcoal.h
//  CellCoal
//
//  Created by David Posada on 21/11/2017.
//  Copyright © 2017 David Posada. All rights reserved.
//

#ifndef cellcoal_h
#define cellcoal_h

#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <errno.h>
#include "platform.h"

#include "definitions.h"

/* Scratch matrices of a thread that simulates reads */
typedef struct
	{
	double	**ampMat, **ampPat, *sbuf, *probsT;
	}
	ReadScratch;

/* Prototypes */
void 	PrintHeader (FILE *fp);
void	PrintCheck (FILE *fp);
extern int	checkOnly;
extern void 	PrintUsage (FILE *fp);
void 	PrintDate (FILE *fp);
void 	PrintDefaults (FILE *fp);
void		ReadUntil (FILE *fv, char stopChar, char *what);
void		PrintRunInformation (FILE *fp);
void		ReadParametersFromFile (void);
void		ReadParametersFromCommandLine (int argc, char **argv);
void 	PrintCommandLine (FILE *fp, int argc,char **argv);
void		PrepareGlobalFiles (int argc, char **argv);
void		PrepareSeparateFiles (int replicate);
void		WriteADOCells (const double *rate, int replicate);
void		KeepCellADO (const double *rate);
int			ConvertVCFtoBCF (const char *vcfPath, const char *bcfPath);
extern int	doBCF;	/* --bcf: write the VCF as BCF */
extern int	adoReplicate;	/* replicate being simulated, for the ADO file */
void		ReadUserGenome (FILE *fp);
int			CountUserTreeLeaves (const char *fileName);
void		ReadUserTree (FILE *fp);
void		FreeUserTree (void);
void		CheckTree (char *treeString);
void 	RelabelUserTree (TreeNode *p);
void		PrintTree (TreeNode *treeRoot, FILE *fp);
void		WriteTree (TreeNode *p, FILE *fp);
void		PrintTimes (int listPosition, FILE *fp);
void		ListTimes (int position, FILE *fp);
void		PrintSNVGenotypes (FILE *fp);
void		PrintSNVHaplotypes (FILE *fp, int PrintTrueVariants);
void 	PrintFullGenotypes (FILE *fp);
void 	PrintFullHaplotypes (FILE *fp);
void 	PrintTrueFullHaplotypes (FILE *fp);
void		PrintSiteInfo (FILE *fp, int site);
void		PrintCATG (FILE *fp);
void		PrintVCF (FILE *fp);
void		PrintPILEUP (FILE *fp);
void		AllelicDropout (long int *seed);
void		GenotypeError (long int *seed);
void		GenerateReadCounts (long int *seed);
void 	StreamReads (long int *seed, FILE *fpVCFout, FILE *fpPILEUPout, FILE *fpBulkOut, FILE *fpCATGbody);
void 	PrintCATGHeader (FILE *fp);
extern unsigned char	*mlStore;
extern long				mlStoreCap;
void 	PrintBulk (FILE *fp);
void 	WriteRunRecord (double wallSecs, double cpuSecs);
int 	VCFSitePasses (int j);
FILE *	OpenVCF (const char *path);
void 	CloseVCF (FILE *fp);
void 	SiteReadCounts (CellSiteStr *cellstr, int cell, int site, double *probs, double **ngsEij, double **ampEijmat, double **ampEijpat, long int *seed);
void 	GenotypeLikelihoods (CellSiteStr *cellstr, int cell, int site, double *probs, double **ngsEij, double **ampEijmat, double **ampEijpat, long int *seed, double *glOut);
void	BuildAmpMatrices (CellSiteStr *c, double **ampMat, double **ampPat);
void	SiteGL (CellSiteStr *c, int i, int j, double *glOut);
const double	*GetSiteGL (int i, int j, double *buf);
void 	SetupDoublets (long int *seed, long windowSites);
void 	SetupCellCoverage (long int *seed);
extern char		WhichNuc (int nucleotide);
int		WhichNucChar (char nucleotide);
char		WhichMut (int state);
char		WhichIUPAC (int allele1, int allele2);
char 	WhichConsensusBinary (int allele1, int allele2);
void     MakeCoalescenceTree (int numCells, int N, long int *seed);
int		Index (TreeNode *p);
int 		Label (TreeNode *p);
void 	EvolveSitesOnTree (TreeNode *p, int genome, long int *seed);
void		RelabelNodes (TreeNode *p);
void     MakeTreeNonClock (TreeNode *p, long int *seed);
void		TreeRateSwitch (TreeNode *p, long int *seed);
void 	BranchRateSwitch (TreeNode *p, float multiplier, int *numTips, long int *seed);
void     InitializeGenomes (TreeNode *p, long int *seed);
void     AddGermlineVariation (long int *seed);
void     InheritRootGenome (TreeNode *p);
void     SimulateMk2 (TreeNode *p, int genome, long int *seed);
void     SimulateMk2forSite (TreeNode *p, int genome, int site, long int *seed);
void 	SimulateISM (TreeNode *p,  int genome, int doISMhaploid, long int *seed);
void 	SimulateISMforSite (TreeNode *p,  int genome, int site, int doISMhaploid, long int *seed);
void 	SimulateISMDNAforSite (TreeNode *p,  int genome, int site, int doISMhaploid, long int *seed);
void 	SimulateFiniteDNA (TreeNode *p, int genome, long int *seed);
void		SimulateFiniteDNAforSite (TreeNode *p, int genome, int site, long int *seed);
void				FreeTriNucleotides (void);
void 	SimulateSignatureISM (TreeNode *p,  int genome, long int *seed);
void		SimulateSignatureISMforSite (TreeNode *p,  int genome, int site, int newState, long int *seed);
void		SimulateTriNucFreqGenome (int cell, long int *seed);
void 	CountTriNucFrequencies (Geno *genome_array, int genome);
int		ChooseTrinucleotideSite (long int *seed, int *newState, int genome);
void		FillSubstitutionMatrix (double ch_prob[4][4], double branchLength);
void		JCmodel (double Pij[4][4], double branchLength);
void		HKYmodel (double Pij[4][4], double branchLength);
void		GTRmodel (double Pij[4][4], double branchLength);
void 	EvolveDeletionsOnTree (TreeNode *p, int genome, long int *seed);
void 	EvolveCNLOHonTree (TreeNode *p, int genome, long int *seed);
void		SimulateDeletionforSite (TreeNode *p, int genome, int site, long int *seed);
void		SimulateCNLOHforSite (TreeNode *p, int genome, int site, long int *seed);
int 		CountTrueVariants(void);
int		CountAllelesInObservedGenotypes (void);
int		CountAllelesInMLGenotypes (void);
double	SumBranches (TreeNode *p);
double 	RandomUniform (long int *seed);
int 		RandomUniformTo (int max, long int *seed);
int		RandomPoisson (double lambda, long int *seed);
double	RandomExponential (double mean, long int *seed);
int		RandomBinomial (double prob, int numTrials, long int *seed);
int		RandomNegativeBinomial (double mean, double dispersion,  long int *seed);
double	RandomBetaMeanVar (double mean, double var, long int *seed);
/* static double	RandomBeta (double alpha, double beta, long int *seed);*/
double	RandomGamma (double shape, long int *seed);
double	RandomGamma1 (double s, long int *seed);
double	RandomGamma2 (double s, long int *seed);
int		ChooseUniformState (double *freq, long int *seed);
int		RemoveDir(char *path);
int 		CheckMatrixSymmetry(double matrix[4][4]);
void 	RecordTriNucObservation(TriNucStr *trin, int site);
void 	RemoveTriNucObservation(TriNucStr *trin, int site);
void		AllocateCellStructure (void);
int		CompareGenotypes (int a1, int a2, int b1, int b2);
#define NUM_TRUTH_METRICS	12
int		TruthMetrics (double *m);	/* m[NUM_TRUTH_METRICS]: called_fraction, genotype_error, snv_precision, snv_recall, het_error, hom_error, het_to_hom, doublet_error, singlet_error, ado_simulated */
#ifdef COUNT_ML_GENOTYPE_ERRORS
	void		CountMLGenotypingErrors (void);
#endif

extern int 		EigenREV (double Root[], double Cijk[]);
extern void		PrepareGeneticSignatures(void);


/* Functions shared between modules */
void MLPairAt (int i, int q, int *m, int *p);
int MLCountsAtSite (int j);
void NewReadScratch (ReadScratch *rs);
void FreeReadScratch (ReadScratch *rs);
void SimulateSiteReads (int jj, uint64_t streamBase, double **ngsEij, ReadScratch *rs);
void SetupReadModel (double **ngsEij);

char **ExpandConfigArgs (int *argc, char **argv);
int JoinSplitDecimalArgs (int argc, char **argv);
const char *SweepOptionForName (const char *name);
void RunSweep (int argc, char **argv, const char *spec, int keep, int jobs, int resume);
void WriteOptionsJson (FILE *fp, const char *commandLine, const char *indent);

/* Global variables */
extern TreeNode		*treeNodes, *coalTreeMRCA, *healthyRoot, *healthyTip;
extern SiteStr			*allSites;
extern TriNucStr		*triNucleotideMaternal, *triNucleotidePaternal;
extern CellStr			*cell;
extern long int		userSeed, originalSeed;
extern Geno		***data;
extern int      *SNVsites, *DefaultModelSites, *AltModelSites, *variantSites;
extern int		tipLabel, intLabel;
extern int		ploidy, numCells, N, *Nbegin, *Nend,  *cumDuration, numSites, numDataSets, numPeriods, maxTreeRateSwitches;
extern int		noisy, numNodes, numAltModelSites, numDefaultModelSites, numISMmutations, altModel;
extern int		numCA, numMU, numDEL, numCNLOH, numProposedMU, numSNVs, numFixedMutations, numSNVmaternal, zeroSNVs;
extern int		numISMdeletions, numISMCNLOH, stringPrecision;
/* statistics accumulated over replicates (sum, sum of squares; mean and variance after the last replicate) */
typedef struct { double sum, sumSq, mean, var; } Accum;
typedef struct { Accum tmrca, ca, mu, snv, del, cnloh, mlErr, called; } RunStats;
extern RunStats runStats;
void AccumAdd (Accum *a, double x);
void AccumFinish (Accum *a, int n);
extern double	expNumMU, expVarNumMU;
extern double	theta, healthyTipBranchLength, healthyTipBranchLengthRatio, transformingBranchLengthRatio, totalTreeLength;
extern double	mutationRate, nonISMRelMutRate, propAltModelSites, altModelMutationRate, deletionRate, CNLOHrate;
extern char		SNVgenotypesFile[MAX_NAME], SNVhaplotypesFile[MAX_NAME], trueHaplotypesFile[MAX_NAME],MLhaplotypesFile[MAX_NAME], fullGenotypesFile[MAX_NAME], fullHaplotypesFile[MAX_NAME];
extern char		treeFile[MAX_NAME], timesFile[MAX_NAME], CATGfile[MAX_NAME], VCFfile[MAX_NAME], PILEUPfile[MAX_NAME], bulkFile[MAX_NAME], logFile[MAX_NAME], settingsFile[MAX_NAME], userTreeFile[MAX_NAME], userGenomeFile[MAX_NAME], parameterFile[MAX_NAME];
extern char		SNVgenotypesDir[MAX_NAME], SNVhaplotypesDir[MAX_NAME], trueHaplotypesDir[MAX_NAME], MLhaplotypesDir[MAX_NAME], fullGenotypesDir[MAX_NAME], fullHaplotypesDir[MAX_NAME];
extern char		treeDir[MAX_NAME], timesDir[MAX_NAME], CATGdir[MAX_NAME], VCFdir[MAX_NAME], PILEUPdir[MAX_NAME], bulkDir[MAX_NAME];
extern char		resultsDir[MAX_NAME], treeDir[MAX_NAME], timesDir[MAX_NAME], File[4*MAX_NAME], *CommandLine, *treeString, *taxonName, **cellNames, *treeRateInfo;
extern char		inCellName[MAX_NAME], outCellName[MAX_NAME], inRootCellName[MAX_NAME], outRootCellName[MAX_NAME];
extern char		GLmodel[MAX_NAME];
extern int		doPrintSNVgenotypes, doPrintSNVhaplotypes, doPrintTrueHaplotypes, doPrintFullHaplotypes, doPrintFullGenotypes, doPrintTree, doUserTree, doUserGenome, doSpecificParameterFile;
extern int		doPrintTimes, doPrintAncestors, doPrintCATG, doPrintPILEUP, doPrintBulk, doPrintVCF, doZipVCF, doPrintSeparateReplicates, doPrintIUPAChaplotypes;
extern int		doExponential, doDemographics, doOthsukiInnanCoal, doSimulateData, doSimulateFixedNumMutations, doNGS, doTumorNames, taxonNamesAreChars;
extern int		doJC, doHKY, doGTR, doGTnR, doGeneticSignatures;
extern int      rateVarAmongSites, rateVarAmongLineages, rateVarCoverage, doTreeRateSwitches, numTreeRateSwitches, equalBaseFreq, alphabet, thereIsMij, thereIsEij;
extern double	*periodGrowth, growthRate, birthRate, deathRate, *treeRateSwitchMultiplier;
extern double	meanGenotypingError, varGenotypingError, sequencingError;
extern double	fixedADOrate, doADOcell, meanADOcell, varADOcell, doADOsite, meanADOsite, varADOsite, thereisADO;
extern double	doAllelicImbalance, meanAllelicImbalance, varAllelicImbalance;
extern double	coverage, haploidCoverageReduction;
extern double	meanAmplificationError, varAmplificationError, meanDoubletRate, varDoubletRate;
extern double	TMRCA, expTMRCA, expVarTMRCA;
extern double   titv, kappa, beta, freqR, freqY, freqAG, freqCT, freq[4], cumfreq[4], Mij[4][4], cumMij[4][4], Eij[4][4], cumEij[4][4], alphaSites, alphaBranches;
extern double	Cijk[256], Root[4];
extern double	SNPrate, alphaCoverage, alphaSiteCoverage, alphaCellCoverage;
extern int		OUTGROUP_ROOT, INGROUP_ROOT;

/* Genotype storage. In dense mode, data[genome][node][site] holds the genotype of every node at every site.
   In sparse mode (very large data sets), only the genome of the root is stored densely (data[genome][OUTGROUP_ROOT]);
   a site where some node differs from the root gets a column sparseCols[site] with the genotypes of all the nodes
   (maternal nodes first, then paternal). Sites without a column are identical to the root in all the nodes. */
extern int		sparseMode;
extern Geno		**sparseCols;
extern int		numNodes;
extern int		ignoredNumCells;
extern int		numOutgroup;
#define NUM_TAXA (numCells + numOutgroup)	/* cells that appear in the outputs (the outgroup is cell number numCells) */
extern long		numSparseCols;
/* In sparse mode, ADO and genotyping errors are not written in the genomes (that would create a column for almost
   every site): they are computed on the fly for the cells (see ObservedAllele) once obsActive is set. */
extern int			obsActive;
/* Per-cell, per-site data lives in cell[i].site / cell[i].doublet. When it is streamed (see StreamReads), those arrays only
   hold a window of sites starting at siteOff; otherwise siteOff is 0 and they hold all the sites. */
extern long			siteOff;
#define CS(i,j)	(cell[i].site[(long)(j) - siteOff])
#define CD(i,j)	(cell[i].doublet[(long)(j) - siteOff])
extern double		minVAFvcf;	/* the VCF lists only sites whose pooled VAF in reads is at least this */
extern struct timespec	wallStart;	/* start of the run, for the elapsed (wall clock) time */
extern int			readsOnly;		/* only the pooled (bulk) read file is printed: no genotype likelihoods are needed */
extern int			streamMode;		/* reads are simulated and printed in blocks of sites (see StreamReads) */
extern FILE			*catgBody;	/* temporary file with the CATG site lines of a streamed replicate */
extern int			streamBlock;		/* number of sites in a block */
extern int			*doubletDonor;	/* donor cell of each doublet cell when streaming */
extern int			sparseADO;
extern int			*errCount;
void		CorrelateDoublets (long int *seed);
void		WriteCellParams (int replicate);
void		StoreFullBlock (int j0, int nb);
void		LoadFullRow (int i);
void		FullMLPair (int i, int j, int *m, int *p);
void		CloseFullStore (void);
Geno		ObservedAllele (int genome, int node, int site, Geno trueAllele);
void		SparseAllelicDropout (long int *seed);
void		SparseGenotypeError (long int *seed);
void		ResetObservationLayer (void);
Geno		*NewSparseColumn (int site);
void		DecideSparseMode (void);
void		ResetSparseColumns (void);
static inline Geno GetG (int genome, int node, int site)
	{
	if (__builtin_expect (sparseMode, 0) && node != OUTGROUP_ROOT)
		{
		Geno *c = sparseCols[site];
		Geno b = c ? c[(long)genome * numNodes + node] : data[genome][OUTGROUP_ROOT][site];
		if (__builtin_expect (obsActive, 0) && node <= numCells && (sparseADO || (errCount != NULL && errCount[site] > 0)))
			return ObservedAllele (genome, node, site, b);
		return b;
		}
	return data[genome][node][site];
	}
static inline Geno SetG (int genome, int node, int site, Geno value)
	{
	if (__builtin_expect (sparseMode, 0) && node != OUTGROUP_ROOT)
		{
		Geno *c = sparseCols[site];
		if (c == NULL)
			{
			if (value == data[genome][OUTGROUP_ROOT][site])
				return value;
			c = NewSparseColumn (site);
			}
		c[(long)genome * numNodes + node] = value;
		return value;
		}
	data[genome][node][site] = value;
	return value;
	}
extern int		readingParameterFile, simulateOnlyTwoTemplates;
extern int 		doGATK, do2T, do4T, doGATK_ADO, do2T_ADO, do4T_ADO;
extern double	**siteGL;	/* per site: genotype likelihoods of all cells (GL_NUM doubles each), kept only for SNV sites and within the memory budget; NULL otherwise (recomputed when needed) */
extern size_t	glBudget, glUsed;
extern int		TipNodeNum, IntNodeNum;
extern char		*maternalUserGenome, *paternalUserGenome;
extern int		complementBase[4];
extern int		targetTriChange[6];
extern int		*triMutationsCounter;
extern int		numUserSignatures;
extern double	*signatureWeight;
extern int		*signatureID;
extern double 	Qij[16], mr;
extern double 	***selectedSignature;
extern double	****geneticSignature;
extern double	**signatureProbs;
extern double	*triNucFreq;

/* The genotype of a cell without allelic dropout and genotyping errors (sparse mode keeps them apart from the genomes) */
static inline Geno GetTrueG (int genome, int node, int site)
	{
	if (__builtin_expect (sparseMode, 0) && node != OUTGROUP_ROOT)
		{
		Geno *c = sparseCols[site];
		return c ? c[(long)genome*numNodes+node] : data[genome][OUTGROUP_ROOT][site];
		}
	return data[genome][node][site];
	}

#ifdef CHECK_MUT_DISTRIBUTION
extern int		*MutCount, *SiteMut, dataSetsWithSNVs;
extern double	sumPos, meansumPos;
#endif

/* File pointers */
extern FILE			*fpSNVgenotypes, *fpFullGenotypes, *fpSNVhaplotypes, *fpTrueHaplotypes, *fpFullHaplotypes, *fpMLhaplotypes, *fpTrees, *fpTimes, *fpCATG, *fpVCF, *fpPILEUP, *fpBulk, *fpLog, *fpUserTree, *fpUserGenome, *fpSims;


/* per-cell, per-site data are needed for NGS reads and for the haplotype, full-genotype and ancestor outputs (see main) */
#define NEED_CELL_SITES (doNGS == YES || doPrintSNVhaplotypes == YES || doPrintTrueHaplotypes == YES || \
                         doPrintFullHaplotypes == YES || doPrintFullGenotypes == YES || \
                         doPrintCATG == YES || doPrintPILEUP == YES || doPrintBulk == YES || doPrintAncestors == YES)

#endif /* coaltumor_h */

/* Independent random streams for the read-count simulation.
   Each site gets its own stream (a 64-bit SplitMix64 generator seeded from the main stream and the site
   number), so that the sites can be simulated in parallel and the result does not depend on the number of
   threads or on the order in which they are processed. While a thread has a stream active, RandomUniform
   ignores the seed it is given and draws from the stream instead. */
extern double			*cellCovMult;	/* between-cell coverage multiplier of every cell (-J) */
extern __thread double		siteCovMult;	/* between-site coverage multiplier of the site being simulated (-O) */
extern __thread int			streamActive;
extern double				covAdoCorr;	/* --cov-ado-corr: correlation between the coverage and the ADO rate of the cells */
extern double				covDoubletCorr, adoDoubletCorr;	/* --cov-doublet-corr, --ado-doublet-corr: correlation of being a doublet with the coverage / ADO rate of the cells */
extern int					doSummary, numSNVsTrue;	/* --summary: write summary.tsv, one row of metrics per replicate */
extern FILE				*fpSummary;
extern int					minAltReads;
extern double				minAltVaf;
extern double				adoBlockLen;
extern long long			adoOrigin;
extern long				adoSeedOpt;
void		BuildAdoBlocks (long int *seed);
void		FreeAdoBlocks (void);
int			AdoBlockLost (int node, int genome, long site);
extern int			adoBlocksOn;
extern int					coverageWindow, covNumWindows;
extern long long			*covSum, *covSites;
extern FILE					*fpCoverage;
void						CoverageAddSite (int site);
void						CoverageWrite (int replicate);
int							CellAltCall (int cellIndex, int site);
void						InvariantGenotypeCounts (long long *cnt);
extern double				ampCellVar;	/* --amp-cell-var: variance among cells of the mean amplification error (-A) */
extern double				*cellAmpMean;	/* mean amplification error of every cell (with --amp-cell-var) */
extern double				*cellAdoRate;	/* ADO rate of every cell (kept only when covAdoCorr != 0) */
extern int					fastMode;	/* --fast: cheaper samplers and log10/10^x for the reads (same distributions, different random numbers) */
double	Log10Approx (double x);
double	Pow10Approx (double x);
void	Log10Batch (const double *x, double *y, int n);
void	Pow10Batch (const double *x, double *y, int n);
extern __thread uint64_t	streamState;

static inline uint64_t SplitMix64 (uint64_t *x)
	{
	uint64_t z = (*x += 0x9E3779B97F4A7C15ULL);

	z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
	z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
	return z ^ (z >> 31);
	}

/* starts the stream of a given item (site) from a base value taken from the main stream */
static inline void StartStream (uint64_t base, uint64_t item)
	{
	uint64_t x = base + item * 0xD1B54A32D192ED03ULL;

	streamState = SplitMix64 (&x);	/* mix, so that neighbouring items start far apart */
	streamActive = 1;
	}

