//
//  coaltumor.h
//  CoalTumor
//
//  Created by David Posada on 21/11/2017.
//  Copyright © 2017 David Posada. All rights reserved.
//

#ifndef coaltumor_h
#define coaltumor_h

#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <ftw.h>

#define PROGRAM_NAME            "CoalTumor"
#define PROGRAM_NAME_UPPERCASE	"COALTUMOR"
#define VERSION_NUMBER          "version 0.1"
#define	MAX_NAME                120
#define	MAX_LINE                3500
#define	NO                      0
#define	YES                     1
#define	infinity                999999
#define	MATERNAL                0
#define	PATERNAL                1
#define BINARY                  0
#define	DNA                     1

#define ISMhap					0
#define Mk						1
#define finiteDNA				2

#define A						0
#define C						1
#define G						2
#define T						3
#define ADO						9
#define DELETION				7

#define ISMhap					0
#define Mk2						1
#define finiteDNA				2

#define CHECK_MUT_DISTRIBUTION
#undef CHECK_MUT_DISTRIBUTION

#define MYDEBUG
#undef MYDEBUG

typedef struct node
    {
    struct node		*left, *right, *anc;
    int				index, label, isHealthyTip, isHealthyRoot;
    double			length, time, branchLength;
    char			*name;
    }
    TreeNode;

typedef struct site
    {
    int		isSNV, isSNP, isVariant;
	int		numMutations, numMutationsMaternal, numMutationsPaternal;
    int     numDeletions, numDeletionsMaternal, numDeletionsPaternal;
	int		hasADO;
	int		hasGenotypeError;
	int		countA, countC, countG, countT, countACGT, countCellswithData, countDropped;
	int		referenceAllele;
	int		*alternateAlleles;
	int		numAltAlleles;
	double	branchSum;
	double	deletionBranchSum;
    double  rateMultiplier;
	}
    SiteStr;

/*
typedef struct cell
    {
    int				index;
	int				**genome;
	int				**readCount;
	double			**genLike;
    char			*name;
    }
    CellStr;
*/

/* Prototypes */
static void 	PrintHeader (FILE *fp);
static void 	PrintDate (FILE *fp);
static void 	PrintUsage (void);
static void 	PrintDefaults (FILE *fp);
static void		ReadUntil (FILE *fv, char stopChar, char *what);
static void		PrintRunInformation (FILE *fp);
static void		ReadParametersFromFile (void);
static void		ReadParametersFromCommandLine (int argc, char **argv);
static void 	PrintCommandLine (FILE *fp, int argc,char **argv);
static void		PrepareGlobalFiles (int argc, char **argv);
static void		PrepareSeparateFiles (int replicate);
static void		ReadUserTree (FILE *fp);
static void		CheckTree (char *treeString);
static void 	RelabelUserTree (TreeNode *p);
static void		PrintTree (TreeNode *treeRoot, FILE *fp);
static void		WriteTree (TreeNode *p, FILE *fp);
static void		PrintTimes (int listPosition, FILE *fp);
static void		ListTimes (int position, FILE *fp);
static void		PrintSNVGenotypes (FILE *fp);
static void		PrintSNVHaplotypes (FILE *fp, int PrintTrueVariants);
static void 	PrintFullGenotypes (FILE *fp);
static void 	PrintFullHaplotypes (FILE *fp);
static void		PrintSiteInfo (FILE *fp, int site);
static void		AllelicDropout (long int *seed);
static void		GenotypeError (long int *seed);
static void		GenerateReadCounts (long int *seed);
static char		WhichNuc (int nucleotide);
static char		WhichMut (int state);
static char		WhichIUPAC (int allele1, int allele2);
static char 	WhichConsensusBinary (int allele1, int allele2);
static void     MakeCoalescenceTree (int numCells, int N, long int *seed);
static int		Index (TreeNode *p);
static int 		Label (TreeNode *p);
static void 	EvolveSitesOnTree (TreeNode *p, int genome, long int *seed);
static void		RelabelNodes (TreeNode *p);
static void     MakeTreeNonClock (TreeNode *p, long int *seed);
static void     InitializeGenomes (TreeNode *p, long int *seed);
static void     AddGermlineVariation (long int *seed);
static void     SimulateMk2 (TreeNode *p, int genome, long int *seed);
static void     SimulateMk2forSite (TreeNode *p, int genome, int site, long int *seed);
static void 	SimulateISM (TreeNode *p,  int genome, int doISMhaploid, long int *seed);
static void 	SimulateISMforSite (TreeNode *p,  int genome, int site, int doISMhaploid, long int *seed);
static void 	SimulateISMDNAforSite (TreeNode *p,  int genome, int site, int doISMhaploid, long int *seed);
static void 	SimulateFiniteDNA (TreeNode *p, int genome, long int *seed);
static void		SimulateFiniteDNAforSite (TreeNode *p, int genome, int site, long int *seed);
static void		FillSubstitutionMatrix (double ch_prob[4][4], double branchLength);
static void		JCmodel (double Pij[4][4], double branchLength);
static void		HKYmodel (double Pij[4][4], double branchLength);
static void		GTRmodel (double Pij[4][4], double branchLength);
static void 	LoadGeneticSignatures (void);
static void 	EvolveDeletionsOnTree (TreeNode *p, int genome, long int *seed);
static void		SimulateDeletionforSite (TreeNode *p, int genome, int site, long int *seed);
static int 		CountTrueVariants(void);
static int		CountAlleles (void);
static double	SumBranches (TreeNode *p);
static double 	RandomUniform (long int *seed);
static int 		RandomUniformTo (int max, long int *seed);
static int		RandomPoisson (double lambda, long int *seed);
static double	RandomExponential (double mean, long int *seed);
static int		RandomBinomial (double prob, int numTrials, long int *seed);
static int		RandomNegativeBinomial (double mean, double dispersion,  long int *seed);
static double	RandomBeta (double mean, double var, long int *seed);
static double	RandomGamma (double shape, long int *seed);
static double	RandomGamma1 (double s, long int *seed);
static double	RandomGamma2 (double s, long int *seed);
static int		ChooseUniformState (double *freq, long int *seed);
static int		Unlink_callback(const char *fpath, const struct stat *sb, int typeflag, struct FTW *ftwbuf);
static int		RemoveDir(char *path);
static int 		CheckMatrixSymmetry(double matrix[4][4]);
extern int 		EigenREV (double Root[], double Cijk[]);


/* Global variables */
TreeNode		*treeNodes, *coalTreeMRCA, *healthyRoot, *healthyTip;
SiteStr			*allSites;
//CellStr			*cell;
long int		userSeed, originalSeed;
static int		***data;
static int      *SNVsites, *DefaultModelSites, *AltModelSites, *variantSites;
static int		tipLabel, intLabel;
static int		ploidy, numCells, N, *Nbegin, *Nend,  *cumDuration, numSites, numDataSets, numPeriods;
static int		noisy, numNodes, numAltModelSites, numDefaultModelSites, numISMmutations, altModel;
static int		numCA, numMU, numDEL, numProposedMU, numSNVs, numFixedSNVs, numSNVmaternal, zeroSNVs;
static int		numISMdeletions;
static double	meanNumSNVs, meanNumCA, meanNumMU, meanNumDEL;
static double 	cumNumSNVs, cumNumCA, cumNumMU, cumNumDEL;
static double	cumNumMUSq, cumNumSNVsSq, cumNumDELSq;
static double	varNumMU, varNumSNVs, varNumDEL;
static double	expNumMU, expVarNumMU;
static double	theta, healthyTipBranchLength, transformingBranchLength;
static double	mutationRate, nonISMRelMutRate, propAltModelSites, altModelMutationRate;
static double	deletionRate;
static char		SNVgenotypesFile[MAX_NAME], SNVhaplotypesFile[MAX_NAME], SNVtrueHaplotypesFile[MAX_NAME],fullGenotypesFile[MAX_NAME], fullHaplotypesFile[MAX_NAME];
static char		treeFile[MAX_NAME], timesFile[MAX_NAME], CATGfile[MAX_NAME], VCFfile[MAX_NAME], logFile[MAX_NAME], settingsFile[MAX_NAME], userTreeFile[MAX_NAME];
static char		SNVgenotypesDir[MAX_NAME], SNVhaplotypesDir[MAX_NAME], SNVtrueHaplotypesDir[MAX_NAME], fullGenotypesDir[MAX_NAME], fullHaplotypesDir[MAX_NAME];
static char		treeDir[MAX_NAME], timesDir[MAX_NAME], CATGdir[MAX_NAME], VCFdir[MAX_NAME];
static char		resultsDir[MAX_NAME], treeDir[MAX_NAME], timesDir[MAX_NAME], File[MAX_NAME], *CommandLine, *treeString, *taxonName, **cellNames;
static int		doPrintSNVgenotypes, doPrintSNVhaplotypes, doPrintSNVtrueHaplotypes, doPrintFullHaplotypes, doPrintFullGenotypes, doPrintTree, doUserTree;
static int		doPrintTimes, doPrintAncestors, doPrintCATG, doPrintSeparateReplicates, doPrintIUPAChaplotypes;
static int		doExponential, doDemographics, doSimulateData, doSimulateFixedNumSNVs,doSimulateReadCounts, taxonNamesAreChars;
static int		doJC, doHKY, doGTR, doGTnR, doGeneticSignatures, geneticSignature;
static int      rateVarAmongSites, rateVarAmongLineages, rateVarCoverage, equalBaseFreq, alphabet, thereIsMij, thereIsEij, coverage;
static double	*periodGrowth, growthRate, sequencingError, ADOrate, singleAlleleCoverageReduction, genotypingError, meanAmplificationError, varAmplificationError;
static double	TMRCA, cumTMRCA, cumTMRCASq, meanTMRCA, expTMRCA, varTMRCA, expVarTMRCA;
static double   titv, kappa, beta, freqR, freqY, freqAG, freqCT, freq[4], cumfreq[4], Mij[4][4], cumMij[4][4], Eij[4][4], cumEij[4][4], alphaSites, alphaBranches;
static double	Rmat[6], NRmat[12], Cijk[256], Root[4];
static double	SNPrate, alphaCoverage;
static int		HEALTHY_ROOT, TUMOR_ROOT;
static int		readingParameterFile, simulateOnlyTwoTemplates;
static int		TipNodeNum, IntNodeNum;
extern double 	Qij[16], mr;


#ifdef CHECK_MUT_DISTRIBUTION
static int		*MutCount, *SiteMut, dataSetsWithSNVs;
static double	sumPos, meansumPos;
#endif

/* File pointers */
FILE			*fpSNVgenotypes, *fpFullGenotypes, *fpSNVhaplotypes, *fpSNVtrueHaplotypes, *fpFullHaplotypes, *fpTrees, *fpTimes, *fpCATG, *fpVCF, *fpLog, *fpUserTree;



#endif /* coaltumor_h */
