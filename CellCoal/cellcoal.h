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
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <ftw.h>

#include "definitions.h"


/* Prototypes */
static void 	PrintHeader (FILE *fp);
extern void 	PrintUsage (void);
static void 	PrintDate (FILE *fp);
static void 	PrintDefaults (FILE *fp);
static void		ReadUntil (FILE *fv, char stopChar, char *what);
static void		PrintRunInformation (FILE *fp);
static void		ReadParametersFromFile (void);
static void		ReadParametersFromCommandLine (int argc, char **argv);
static void 	PrintCommandLine (FILE *fp, int argc,char **argv);
static void		PrepareGlobalFiles (int argc, char **argv);
static void		PrepareSeparateFiles (int replicate);
static void		ReadUserGenome (FILE *fp);
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
extern char		WhichNuc (int nucleotide);
static int		WhichNucChar (char nucleotide);
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
static void 	SimulateSignatureISM (TreeNode *p,  int genome, long int *seed);
static void		SimulateSignatureISMforSite (TreeNode *p,  int genome, int site, int newState, long int *seed);
static void		SimulateTriNucFreqGenome (int cell, long int *seed);
static void 	CountTriNucFrequencies (int *genome_array, int genome);
static int		ChooseTrinucleotideSite (long int *seed, int *newState, int genome);
static void		FillSubstitutionMatrix (double ch_prob[4][4], double branchLength);
static void		JCmodel (double Pij[4][4], double branchLength);
static void		HKYmodel (double Pij[4][4], double branchLength);
static void		GTRmodel (double Pij[4][4], double branchLength);
static void 	EvolveDeletionsOnTree (TreeNode *p, int genome, long int *seed);
static void 	EvolveCNLOHonTree (TreeNode *p, int genome, long int *seed);
static void		SimulateDeletionforSite (TreeNode *p, int genome, int site, long int *seed);
static void		SimulateCNLOHforSite (TreeNode *p, int genome, int site, long int *seed);
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
static void 	RecordTriNucObservation(TriNucStr *trin, int site);
static void 	RemoveTriNucObservation(TriNucStr *trin, int site);

extern int 		EigenREV (double Root[], double Cijk[]);
extern void		PrepareGeneticSignatures(void);

/* Global variables */
TreeNode		*treeNodes, *coalTreeMRCA, *healthyRoot, *healthyTip;
SiteStr			*allSites;
TriNucStr		*triNucleotideMaternal, *triNucleotidePaternal;
//CellStr			*cell;
long int		userSeed, originalSeed;
static int		***data;
static int      *SNVsites, *DefaultModelSites, *AltModelSites, *variantSites;
static int		tipLabel, intLabel;
static int		ploidy, numCells, N, *Nbegin, *Nend,  *cumDuration, numSites, numDataSets, numPeriods;
static int		noisy, numNodes, numAltModelSites, numDefaultModelSites, numISMmutations, altModel;
static int		numCA, numMU, numDEL, numCNLOH, numProposedMU, numSNVs, numFixedMutations, numSNVmaternal, zeroSNVs;
static int		numISMdeletions, numISMCNLOH;
static double	meanNumSNVs, meanNumCA, meanNumMU, meanNumDEL, meanNumCNLOH;
static double 	cumNumSNVs, cumNumCA, cumNumMU, cumNumDEL, cumNumCNLOH, cumCountMLgenotypeErrors;
static double	cumNumMUSq, cumNumSNVsSq, cumNumDELSq, cumNumCNLOHSq;
static double	varNumMU, varNumSNVs, varNumDEL, varNumCNLOH;
static double	expNumMU, expVarNumMU;
static double	theta, healthyTipBranchLength, transformingBranchLength, totalTreeLength;
static double	mutationRate, nonISMRelMutRate, propAltModelSites, altModelMutationRate, deletionRate, CNLOHrate;
static char		SNVgenotypesFile[MAX_NAME], SNVhaplotypesFile[MAX_NAME], trueHaplotypesFile[MAX_NAME],fullGenotypesFile[MAX_NAME], fullHaplotypesFile[MAX_NAME];
static char		treeFile[MAX_NAME], timesFile[MAX_NAME], CATGfile[MAX_NAME], VCFfile[MAX_NAME], logFile[MAX_NAME], settingsFile[MAX_NAME], userTreeFile[MAX_NAME], userGenomeFile[MAX_NAME];
static char		SNVgenotypesDir[MAX_NAME], SNVhaplotypesDir[MAX_NAME], trueHaplotypesDir[MAX_NAME], fullGenotypesDir[MAX_NAME], fullHaplotypesDir[MAX_NAME];
static char		treeDir[MAX_NAME], timesDir[MAX_NAME], CATGdir[MAX_NAME], VCFdir[MAX_NAME];
static char		resultsDir[MAX_NAME], treeDir[MAX_NAME], timesDir[MAX_NAME], File[MAX_NAME], *CommandLine, *treeString, *taxonName, **cellNames;
static int		doPrintSNVgenotypes, doPrintSNVhaplotypes, doPrintTrueHaplotypes, doPrintFullHaplotypes, doPrintFullGenotypes, doPrintTree, doUserTree, doUserGenome;
static int		doPrintTimes, doPrintAncestors, doPrintCATG, doPrintSeparateReplicates, doPrintIUPAChaplotypes;
static int		doExponential, doDemographics, doSimulateData, doSimulateFixedNumMutations,doSimulateReadCounts, taxonNamesAreChars;
static int		doJC, doHKY, doGTR, doGTnR, doGeneticSignatures;
static int      rateVarAmongSites, rateVarAmongLineages, rateVarCoverage, equalBaseFreq, alphabet, thereIsMij, thereIsEij, coverage, countMLgenotypeErrors;
static double	*periodGrowth, growthRate, sequencingError, ADOrate, allelicImbalance, genotypingError, meanAmplificationError, varAmplificationError;
static double	TMRCA, cumTMRCA, cumTMRCASq, meanTMRCA, expTMRCA, varTMRCA, expVarTMRCA;
static double   titv, kappa, beta, freqR, freqY, freqAG, freqCT, freq[4], cumfreq[4], Mij[4][4], cumMij[4][4], Eij[4][4], cumEij[4][4], alphaSites, alphaBranches;
static double	Rmat[6], NRmat[12], Cijk[256], Root[4];
static double	SNPrate, alphaCoverage;
static int		HEALTHY_ROOT, TUMOR_ROOT;
static int		readingParameterFile, simulateOnlyTwoTemplates;
static int		TipNodeNum, IntNodeNum;
static char		*maternalUserGenome, *paternalUserGenome;
static int		complementBase[4] = {3,2,1,0};
static int		targetTriChange[6] = {0,2,3,0,1,2};
static int		*triMutationsCounter;
static int		numUserSignatures;
static double	*signatureWeight;
static int		*signatureID;
extern double 	Qij[16], mr;
extern double 	***selectedSignature;
extern double	****geneticSignature;
extern double	**signatureProbs;
extern double	*triNucFreq;

#ifdef CHECK_MUT_DISTRIBUTION
static int		*MutCount, *SiteMut, dataSetsWithSNVs;
static double	sumPos, meansumPos;
#endif

/* File pointers */
FILE			*fpSNVgenotypes, *fpFullGenotypes, *fpSNVhaplotypes, *fpTrueHaplotypes, *fpFullHaplotypes, *fpTrees, *fpTimes, *fpCATG, *fpVCF, *fpLog, *fpUserTree, *fpUserGenome;



#endif /* coaltumor_h */
