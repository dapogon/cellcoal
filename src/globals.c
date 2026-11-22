/* Global variables of CellCoal (declared extern in cellcoal.h) */
#include "cellcoal.h"

unsigned char	*mlStore = NULL;
long				mlStoreCap = 0;
TreeNode		*treeNodes, *coalTreeMRCA, *healthyRoot, *healthyTip;
SiteStr			*allSites;
TriNucStr		*triNucleotideMaternal, *triNucleotidePaternal;
CellStr			*cell;
long int		userSeed, originalSeed;
Geno		***data;
int      *SNVsites, *DefaultModelSites, *AltModelSites, *variantSites;
int		tipLabel, intLabel;
int		ploidy, numCells, N, *Nbegin, *Nend,  *cumDuration, numSites, numDataSets, numPeriods, maxTreeRateSwitches;
int		noisy, numNodes, numAltModelSites, numDefaultModelSites, numISMmutations, altModel;
int		numCA, numMU, numDEL, numCNLOH, numProposedMU, numSNVs, numFixedMutations, numSNVmaternal, zeroSNVs;
int		numISMdeletions, numISMCNLOH, stringPrecision;
RunStats	runStats;
double	expNumMU, expVarNumMU;
double	theta, healthyTipBranchLength, healthyTipBranchLengthRatio, transformingBranchLengthRatio, totalTreeLength;
double	mutationRate, nonISMRelMutRate, propAltModelSites, altModelMutationRate, deletionRate, CNLOHrate;
char		SNVgenotypesFile[MAX_NAME], SNVhaplotypesFile[MAX_NAME], trueHaplotypesFile[MAX_NAME],MLhaplotypesFile[MAX_NAME], fullGenotypesFile[MAX_NAME], fullHaplotypesFile[MAX_NAME];
char		treeFile[MAX_NAME], timesFile[MAX_NAME], CATGfile[MAX_NAME], VCFfile[MAX_NAME], PILEUPfile[MAX_NAME], bulkFile[MAX_NAME], logFile[MAX_NAME], settingsFile[MAX_NAME], userTreeFile[MAX_NAME], userGenomeFile[MAX_NAME], parameterFile[MAX_NAME];
char		SNVgenotypesDir[MAX_NAME], SNVhaplotypesDir[MAX_NAME], trueHaplotypesDir[MAX_NAME], MLhaplotypesDir[MAX_NAME], fullGenotypesDir[MAX_NAME], fullHaplotypesDir[MAX_NAME];
char		treeDir[MAX_NAME], timesDir[MAX_NAME], CATGdir[MAX_NAME], VCFdir[MAX_NAME], PILEUPdir[MAX_NAME], bulkDir[MAX_NAME];
char		resultsDir[MAX_NAME], treeDir[MAX_NAME], timesDir[MAX_NAME], File[4*MAX_NAME], *CommandLine, *paramFileArgs, *treeString, *taxonName, **cellNames, *treeRateInfo;
char		inCellName[MAX_NAME], outCellName[MAX_NAME], inRootCellName[MAX_NAME], outRootCellName[MAX_NAME];
char		GLmodel[MAX_NAME];
int		doPrintSNVgenotypes, doPrintSNVhaplotypes, doPrintTrueHaplotypes, doPrintFullHaplotypes, doPrintFullGenotypes, doPrintTree, doUserTree, doUserGenome, doSpecificParameterFile;
int		doPrintTimes, doPrintAncestors, doPrintCATG, doPrintPILEUP, doPrintBulk, doPrintVCF, doZipVCF, doPrintSeparateReplicates, doPrintIUPAChaplotypes;
int		doExponential, doDemographics, doOthsukiInnanCoal, doSimulateData, doSimulateFixedNumMutations, doNGS, doTumorNames, taxonNamesAreChars;
int		doJC, doHKY, doGTR, doGTnR, doGeneticSignatures;
int      rateVarAmongSites, rateVarAmongLineages, rateVarCoverage, doTreeRateSwitches, numTreeRateSwitches, equalBaseFreq, alphabet, thereIsMij, thereIsEij;
double	*periodGrowth, growthRate, birthRate, deathRate, *treeRateSwitchMultiplier;
double	meanGenotypingError, varGenotypingError, sequencingError;
double	fixedADOrate, doADOcell, meanADOcell, varADOcell, doADOsite, meanADOsite, varADOsite, thereisADO;
double	doAllelicImbalance, meanAllelicImbalance, varAllelicImbalance;
double	coverage, haploidCoverageReduction;
double	meanAmplificationError, varAmplificationError, meanDoubletRate, varDoubletRate;
double	TMRCA, expTMRCA, expVarTMRCA;
double   titv, kappa, beta, freqR, freqY, freqAG, freqCT, freq[4], cumfreq[4], Mij[4][4], cumMij[4][4], Eij[4][4], cumEij[4][4], ampSubRow[4][4], alphaSites, alphaBranches;
double	Cijk[256], Root[4];
double	SNPrate, alphaCoverage, alphaSiteCoverage, alphaCellCoverage;
int		OUTGROUP_ROOT, INGROUP_ROOT;
int		sparseMode = 0;
Geno		**sparseCols = NULL;
int		numNodes = 0;
int		numOutgroup = 1;	/* 1 = the healthy outgroup cell is part of the outputs; 0 with --no-outgroup */
int		ignoredNumCells = 0;	/* value of -s ignored because a user tree fixes the number of cells (0 = none) */
long		numSparseCols = 0;
int			obsActive = 0;
long			siteOff = 0;
double		minVAFvcf = 0.0;	/* the VCF lists only sites whose pooled VAF in reads is at least this */
struct timespec	wallStart;	/* start of the run, for the elapsed (wall clock) time */
int			checkOnly = 0;		/* --check: validate the options and stop without simulating */
int			readsOnly = 0;		/* only the pooled (bulk) read file is printed: no genotype likelihoods are needed */
int			streamMode = 0;		/* reads are simulated and printed in blocks of sites (see StreamReads) */
FILE			*catgBody = NULL;	/* temporary file with the CATG site lines of a streamed replicate */
int			streamBlock = 0;		/* number of sites in a block */
int			*doubletDonor = NULL;	/* donor cell of each doublet cell when streaming */
int			sparseADO = 0;
int			adoBlocksOn = 0;
int			*errCount = NULL;
int		readingParameterFile, simulateOnlyTwoTemplates;
int 		doGATK, do2T, do4T, doGATK_ADO, do2T_ADO, do4T_ADO;
double	**siteGL;	/* per site: genotype likelihoods of all cells (GL_NUM doubles each), kept only for SNV sites and within the memory budget; NULL otherwise (recomputed when needed) */
size_t	glBudget, glUsed;
int		TipNodeNum, IntNodeNum;
char		*maternalUserGenome, *paternalUserGenome;
int		complementBase[4] = {3,2,1,0};
int		targetTriChange[6] = {0,2,3,0,1,2};
int		*triMutationsCounter;
int		numUserSignatures;
double	*signatureWeight;
int		*signatureID;
int		*MutCount, *SiteMut, dataSetsWithSNVs;
double	sumPos, meansumPos;
FILE			*fpSNVgenotypes, *fpFullGenotypes, *fpSNVhaplotypes, *fpTrueHaplotypes, *fpFullHaplotypes, *fpMLhaplotypes, *fpTrees, *fpTimes, *fpCATG, *fpVCF, *fpPILEUP, *fpBulk, *fpLog, *fpUserTree, *fpUserGenome, *fpSims;
double			*cellCovMult = NULL;	/* between-cell coverage multiplier of every cell (-J) */
__thread double		siteCovMult = 1.0;	/* between-site coverage multiplier of the site being simulated (-O) */
__thread int			streamActive = 0;
int					fastMode = 0;
int					doBCF = 0;
double				covAdoCorr = 0.0;
double				*cellAdoRate = NULL;
int					coverageWindow = 0;	/* --coverage W: depth per cell of the sites that are not variable, in windows of W sites (0 = off) */
long long			*covSum = NULL;	/* [window * numCells + cell]: sum of the depth over the non-variable sites of the window */
long long			*covSites = NULL;	/* [window]: number of non-variable sites of the window */
int					covNumWindows = 0;
FILE				*fpCoverage = NULL;
int					doBam = 0, doFastq = 0, readLen = 150;	/* read-level simulator (--bam, --fastq, --read-len) */
double				insertMean = 400.0, insertSd = 60.0, bulkDepth = 0.0, dupMean = 0.3, qualSd = 0.0, qualSlope = 0.0, strandBias = 0.0;
char				*contigName = "chr1";
double				gainSd = 0.0;	/* --gain-sd: standard deviation of the log amplification gain of a haplotype (0 = off) */
double				lossRate = 0.0, gainRate = 0.0, segLen = 1000.0;	/* --seg-loss-rate, --seg-len: copy-number losses of segments (rate per genome and generation, mean length in sites) */
int					numSegLoss = 0, numSegGain = 0, snvsLostBySegments = 0;
double				gainLen = 0.0;	/* --gain-len: mean length (sites) of the segments of constant gain */
int					gainOn = 0;
char				*ampSignatureFile = NULL;	/* --amp-signature: 96-context profile of the amplification errors (COSMIC format) */
int					ampSigOn = 0;
double				ampSigMult[4][4][4];	/* [5' base][base][3' base]: rate of the amplification error relative to the mean (-A) */
double				ampSigRow[4][4][4][4];	/* [5'][base][3'][new base]: probability of each substitution given an error */
char				*signatureFile = NULL;	/* --signature-file: table with the 96-context profiles (COSMIC format) that replaces the built-in signatures */
int					numSignatureProfiles = NUM_SIGNATURES;	/* number of profiles available (30 built-in, or the columns of --signature-file) */
double				adoBlockLen = 0.0;	/* --ado-block L: ADO in contiguous segments of mean length L sites (0 = independent per site) */
long long			adoOrigin = 0;	/* --ado-origin O: position (0-based, in sites) of the first site of this data set along the chromosome */
long				adoSeedOpt = 0;	/* --ado-seed S: seed of the dropout segments and of the ADO rate of the cells (0 = derived from the main seed) */
double				minAltVaf = 0.0;	/* --min-alt-vaf: minimum fraction of the reads of the cell that support the non-reference allele */
int					minAltReads = 0;	/* --min-alt-reads: a site is variable only when a cell has a non-reference call supported by at least this many reads (0 = any non-reference call) */
double				ampCellVar = 0.0;	/* --amp-cell-var: variance among cells of the mean amplification error (-A) */
double				*cellAmpMean = NULL;	/* mean amplification error of every cell (with --amp-cell-var) */
double				covDoubletCorr = 0.0, adoDoubletCorr = 0.0;
int					doSummary = 0, numSNVsTrue = 0;
FILE				*fpSummary = NULL;

int					adoReplicate = 0;
__thread uint64_t	streamState;
