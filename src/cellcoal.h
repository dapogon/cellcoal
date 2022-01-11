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
extern void 	PrintUsage (FILE *fp);
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
static void 	PrintTrueFullHaplotypes (FILE *fp);
static void		PrintSiteInfo (FILE *fp, int site);
static void		PrintCATG (FILE *fp);
static void		PrintVCF (FILE *fp);
static void		PrintPILEUP (FILE *fp);
static void		AllelicDropout (long int *seed);
static void		GenotypeError (long int *seed);
static void		GenerateReadCounts (long int *seed);
static void 	SiteReadCounts (CellSiteStr *cellstr, int cell, int site, double *probs, double **ngsEij, double **ampEijmat, double **ampEijpat, long int *seed);
static void 	GenotypeLikelihoods (CellSiteStr *cellstr, int cell, int site, double *probs, double **ngsEij, double **ampEijmat, double **ampEijpat, long int *seed);
static void 	MakeDoublets (double *probs, double **ngsEij, double **ampEijmat, double **ampEijpat, long int *seed);
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
static void		TreeRateSwitch (TreeNode *p, long int *seed);
static void 	BranchRateSwitch (TreeNode *p, float multiplier, int *numTips, long int *seed);
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
static int		CountAllelesInObservedGenotypes (void);
static int		CountAllelesInMLGenotypes (void);
static double	SumBranches (TreeNode *p);
static double 	RandomUniform (long int *seed);
static int 		RandomUniformTo (int max, long int *seed);
static int		RandomPoisson (double lambda, long int *seed);
static double	RandomExponential (double mean, long int *seed);
static int		RandomBinomial (double prob, int numTrials, long int *seed);
static int		RandomNegativeBinomial (double mean, double dispersion,  long int *seed);
static double	RandomBetaMeanVar (double mean, double var, long int *seed);
/* static double	RandomBeta (double alpha, double beta, long int *seed);*/
static double	RandomGamma (double shape, long int *seed);
static double	RandomGamma1 (double s, long int *seed);
static double	RandomGamma2 (double s, long int *seed);
static int		ChooseUniformState (double *freq, long int *seed);
static int		Unlink_callback (const char *fpath, const struct stat *sb, int typeflag, struct FTW *ftwbuf);
static int		RemoveDir(char *path);
static int 		CheckMatrixSymmetry(double matrix[4][4]);
static void 	RecordTriNucObservation(TriNucStr *trin, int site);
static void 	RemoveTriNucObservation(TriNucStr *trin, int site);
static void		AllocateCellStructure (void);
#ifdef COUNT_ML_GENOTYPE_ERRORS
	static int		CompareGenotypes (int a1, int a2, int b1, int b2);
	static void		CountMLGenotypingErrors (void);
#endif

extern int 		EigenREV (double Root[], double Cijk[]);
extern void		PrepareGeneticSignatures(void);

/* Global variables */
TreeNode		*treeNodes, *coalTreeMRCA, *healthyRoot, *healthyTip;
SiteStr			*allSites;
TriNucStr		*triNucleotideMaternal, *triNucleotidePaternal;
CellStr			*cell;
long int		userSeed, originalSeed;
static int		***data;
static int      *SNVsites, *DefaultModelSites, *AltModelSites, *variantSites;
static int		tipLabel, intLabel;
static int		ploidy, numCells, N, *Nbegin, *Nend,  *cumDuration, numSites, numDataSets, numPeriods, maxTreeRateSwitches;
static int		noisy, numNodes, numAltModelSites, numDefaultModelSites, numISMmutations, altModel;
static int		numCA, numMU, numDEL, numCNLOH, numProposedMU, numSNVs, numFixedMutations, numSNVmaternal, zeroSNVs;
static int		numISMdeletions, numISMCNLOH, stringPrecision;
static double	meanNumSNVs, meanNumCA, meanNumMU, meanNumDEL, meanNumCNLOH;
static double 	cumNumSNVs, cumNumCA, cumNumMU, cumNumDEL, cumNumCNLOH, cumCountMLgenotypeErrors, cumCountCalledGenotypes;
static double	cumNumMUSq, cumNumSNVsSq, cumNumDELSq, cumNumCNLOHSq;
static double	varNumMU, varNumSNVs, varNumDEL, varNumCNLOH;
static double	expNumMU, expVarNumMU;
static double	theta, healthyTipBranchLength, transformingBranchLength, healthyTipBranchLengthRatio, transformingBranchLengthRatio, totalTreeLength;
static double	mutationRate, nonISMRelMutRate, propAltModelSites, altModelMutationRate, deletionRate, CNLOHrate;
static char		SNVgenotypesFile[MAX_NAME], SNVhaplotypesFile[MAX_NAME], trueHaplotypesFile[MAX_NAME],MLhaplotypesFile[MAX_NAME], fullGenotypesFile[MAX_NAME], fullHaplotypesFile[MAX_NAME];
static char		treeFile[MAX_NAME], timesFile[MAX_NAME], CATGfile[MAX_NAME], VCFfile[MAX_NAME], PILEUPfile[MAX_NAME], logFile[MAX_NAME], settingsFile[MAX_NAME], userTreeFile[MAX_NAME], userGenomeFile[MAX_NAME], parameterFile[MAX_NAME];
static char		SNVgenotypesDir[MAX_NAME], SNVhaplotypesDir[MAX_NAME], trueHaplotypesDir[MAX_NAME], MLhaplotypesDir[MAX_NAME], fullGenotypesDir[MAX_NAME], fullHaplotypesDir[MAX_NAME];
static char		treeDir[MAX_NAME], timesDir[MAX_NAME], CATGdir[MAX_NAME], VCFdir[MAX_NAME], PILEUPdir[MAX_NAME];
static char		resultsDir[MAX_NAME], treeDir[MAX_NAME], timesDir[MAX_NAME], File[MAX_NAME], *CommandLine, *treeString, *taxonName, **cellNames, *treeRateInfo;
static char		inCellName[MAX_NAME], outCellName[MAX_NAME], inRootCellName[MAX_NAME], outRootCellName[MAX_NAME];
static char		GLmodel[MAX_NAME];
static int		doPrintSNVgenotypes, doPrintSNVhaplotypes, doPrintTrueHaplotypes, doPrintFullHaplotypes, doPrintFullGenotypes, doPrintTree, doUserTree, doUserGenome, doSpecificParameterFile;
static int		doPrintTimes, doPrintAncestors, doPrintCATG, doPrintPILEUP, doPrintSeparateReplicates, doPrintIUPAChaplotypes;
static int		doExponential, doDemographics, doOthsukiInnanCoal, doSimulateData, doSimulateFixedNumMutations, doNGS, doTumorNames, taxonNamesAreChars;
static int		doJC, doHKY, doGTR, doGTnR, doGeneticSignatures;
static int      rateVarAmongSites, rateVarAmongLineages, rateVarCoverage, doTreeRateSwitches, numTreeRateSwitches, equalBaseFreq, alphabet, thereIsMij, thereIsEij;
static double	*periodGrowth, growthRate, birthRate, deathRate, *treeRateSwitchMultiplier;
static double	meanGenotypingError, varGenotypingError, sequencingError;
static double	fixedADOrate, doADOcell, meanADOcell, varADOcell, doADOsite, meanADOsite, varADOsite, thereisADO;
static double	doAllelicImbalance, meanAllelicImbalance, varAllelicImbalance;
static double	coverage, haploidCoverageReduction;
static double	meanAmplificationError, varAmplificationError, meanDoubletRate, varDoubletRate;
static double	TMRCA, cumTMRCA, cumTMRCASq, meanTMRCA, expTMRCA, varTMRCA, expVarTMRCA;
static double 	cumCountMLgenotypeErrors, cumCountMLgenotypeErrorsSq, meanMLgenotypeError, varMLgenotypeError;
static double 	cumCountCalledGenotypes, cumCountCalledGenotypesSq, meanCalledGenotypes, varCalledGenotypes;
static double   titv, kappa, beta, freqR, freqY, freqAG, freqCT, freq[4], cumfreq[4], Mij[4][4], cumMij[4][4], Eij[4][4], cumEij[4][4], alphaSites, alphaBranches;
static double	Rmat[6], NRmat[12], Cijk[256], Root[4];
static double	SNPrate, alphaCoverage;
static int		OUTGROUP_ROOT, INGROUP_ROOT;
static int		readingParameterFile, simulateOnlyTwoTemplates;
static int 		doGATK, do2T, do4T, doGATK_ADO, do2T_ADO, do4T_ADO;
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
FILE			*fpSNVgenotypes, *fpFullGenotypes, *fpSNVhaplotypes, *fpTrueHaplotypes, *fpFullHaplotypes, *fpMLhaplotypes, *fpTrees, *fpTimes, *fpCATG, *fpVCF, *fpPILEUP, *fpLog, *fpUserTree, *fpUserGenome, *fpSims;


#endif /* coaltumor_h */
