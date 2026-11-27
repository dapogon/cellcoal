/*
	CellCoal - command-line and parameter-file options
	See cellcoal.h for the shared declarations. This file was split from cellcoal.c.
*/
#include "cellcoal.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "definitions.h"

/******************** CopyName **************************/
/* Copies a name given by the user (file or folder name) checking that it fits in the buffer */
static void CopyName (char *dest, const char *src, const char *what)
	{
	if (strlen (src) >= MAX_NAME)
		{
		fprintf (stderr, "\n!!! PARAMETER ERROR: The name for %s is too long (maximum %d characters)\n\n", what, MAX_NAME-1);
		PrintUsage (stderr);
		exit (-1);
		}
	strcpy (dest, src);
	}


/******************** MissingArgument **************************/
/* Reports an option that is missing one or more of its values */
static char *MissingArgument (char flag)
	{
	fprintf (stderr, "\n!!! PARAMETER ERROR: Missing value(s) for option -%c\n\n", flag);
	PrintUsage (stderr);
	exit (-1);
	return NULL;
	}


/******************** ReadParametersFromCommandLine **************************/
/*
 USED IN ORDER
 F n s l e K L g h k q i b u d H j S m p w c f t a r G C V A E D P Q I R B M 1 2 3 4 5 6 7 8 9 Y v x o T U W y #
 
 USED
 a b c d e f g h i j k l m n o p q r s t u v w x y
 A B C D E F G H I   K L M     P Q R S T U V W X Y
 1 2 3 4 5 6 7 8 9 #
*/

/* PowerShell (Windows, and pwsh anywhere) splits a native argument like -E0.001 into
   "-E0" and ".001" before the program sees it. Join such pieces again: an option with no
   dot followed by an argument that starts with a dot and a digit (no CellCoal argument
   starts like that). Works in place on argv; returns the new argc. */
int JoinSplitDecimalArgs (int argc, char **argv)
{
	int		i, n;
	char	*joined;

	if (argc < 3)
		return argc;
	for (i=1, n=1; i<argc; i++)
		{
		if (n > 1 && argv[i][0] == '.' && isdigit ((unsigned char) argv[i][1])
			&& argv[n-1][0] == '-' && argv[n-1][1] != '-' && argv[n-1][1] != '\0'
			&& strchr (argv[n-1], '.') == NULL
			&& isdigit ((unsigned char) argv[n-1][strlen (argv[n-1]) - 1]))
			{
			joined = (char *) Checked (malloc(strlen (argv[n-1]) + strlen (argv[i]) + 1), "Could not allocate the command line arguments\n");
			strcpy (joined, argv[n-1]);
			strcat (joined, argv[i]);
			argv[n-1] = joined;	/* kept until the end of the run, like argv */
			}
		else
			argv[n++] = argv[i];
		}
	argv[n] = NULL;
	return n;
}

void ReadParametersFromCommandLine (int argc,char **argv)
{
    int		i, j;
/* Next command-line value of an option that takes several arguments: stop with an error
   instead of reading past the end of argv when a value is missing. */
#define NEXTARG(i)	(((i)+1 < argc) ? argv[++(i)] : MissingArgument (flag))
    char	flag, ch;
    double  sumPi, sum;
    float	argument;
	double  argdouble;
	
	if (argc > 0 && CommandLine == NULL)	/* a parameter file is read by this same function: keep the first command line */
		{
		CommandLine = (char *) Checked (calloc(MAX_LINE, sizeof(char)), "Could not allocate the command line string\n");
		/* bounded copy: a long command line must not overflow the buffer */
		for (i=0; i<argc; i++)
			{
			size_t used = strlen (CommandLine);
			if (used + 2 >= MAX_LINE)
				break;
			snprintf (CommandLine + used, MAX_LINE - used, " %s", argv[i]);
			}
		}
	
    for (i=1; i<argc; i++)
        {
        if (argv[i][0] != '-' || argv[i][1] == '\0')
            {
            fprintf (stderr, "\n!!! PARAMETER ERROR: Incorrect parameter: %s\n\n", argv[i]);
            PrintUsage (stderr);
            }
        argv[i]++;
        flag =* argv[i];
        argv[i]++;
        argument = -9999;
			
        switch (flag)
			{
			/* user-defined parameter file */
			case 'F':
				ch = *argv[i];
				if(!isspace(ch))
					CopyName (parameterFile, argv[i], "parameterFile");
				i = argc;
				doSpecificParameterFile = YES;
				break;
			case 'n':
                argument = atof(argv[i]);
                numDataSets = (int) argument;
                if (numDataSets <1)
                    {
                    fprintf (stderr, "\n!!! PARAMETER ERROR: Bad number of replicates (%d)\n\n", numDataSets);
                    PrintUsage(stderr);
                    }
                break;
			case 's':
                argument = atof(argv[i]);
                numCells = (int) argument;
                if (numCells < 2)
                    {
                    fprintf (stderr, "\n!!! PARAMETER ERROR: Bad sample size (%d)\n\n", numCells);
                    PrintUsage(stderr);
                    }
                break;
            case 'l':
                argument = atof(argv[i]);
                numSites = (int) argument;
                if (numSites<1)
                    {
                    fprintf (stderr, "\n!!! PARAMETER ERROR: Bad sequence length (%d)\n\n", numSites);
                    PrintUsage(stderr);
                    }
                break;
            case 'e':
                argument = atof(argv[i]);
                N = (int) argument;
                popSizeGiven = 1;
                if (N < 1)
                    {
                    fprintf (stderr, "\n!!! PARAMETER ERROR: Bad effective population size (%d)\n\n", N);
                    PrintUsage(stderr);
                    }
                break;
             case 'g':
                growthRate = atof(argv[i]);
                if (growthRate < 0)
                    {
                    fprintf (stderr, "\n!!! PARAMETER ERROR: Bad growth rate (%f)\n\n", growthRate);
                    PrintUsage(stderr);
                    }
                if (growthRate != 0)
                    {
                    doExponential = YES;
                    if (doDemographics == YES || doOthsukiInnanCoal == YES)
                        {
                        fprintf (stderr, "\n!!! PARAMETER ERROR: Cannot have both exponential (--tree growth) and other demographics (--tree periods, --tree birth, --tree death)\n\n");
                        exit (-1);
                        }
                    }
                break;
			case 'K':
				birthRate = atof(argv[i]);
				if (birthRate < 0)
					{
					fprintf (stderr, "\n!!! PARAMETER ERROR: Bad birth rate (%f)\n\n", birthRate);
					PrintUsage(stderr);
					}
				if (birthRate != 0)
					{
					doOthsukiInnanCoal = YES;
					if (doDemographics == YES || doExponential == YES)
						{
						fprintf (stderr, "\n!!! PARAMETER ERROR: Cannot have both birth rate (--tree birth) and exponential (--tree growth) or demographics (--tree periods)\n\n");
						exit (-1);
						}
					}
				break;
			case 'L':
				deathRate = atof(argv[i]);
				if (deathRate < 0)
					{
					fprintf (stderr, "\n!!! PARAMETER ERROR: Bad death rate (%f)\n\n", deathRate);
					PrintUsage(stderr);
					}
				if (deathRate > birthRate)
					{
					fprintf (stderr, "\n!!! PARAMETER ERROR: Bad death rate (%f), it should be smaller than the birth rate (%f)\n\n", deathRate, birthRate);
					PrintUsage(stderr);
					}
				if (deathRate != 0)
					{
					doOthsukiInnanCoal = YES;
					if (doDemographics == YES || doExponential == YES)
						{
						fprintf (stderr, "\n!!! PARAMETER ERROR: Cannot have both death rate (--tree death) and exponential (--tree growth) or demographics (--tree periods)\n\n");
						exit (-1);
						}
					}
				break;
			case 'h':
				argument = atof(argv[i]);
                numPeriods = (int) argument;
				if (numPeriods > 0)
					doDemographics = YES;
                if (doExponential == YES || doOthsukiInnanCoal == YES)
                    {
                    fprintf (stderr, "\n!!! PARAMETER ERROR: Cannot have both demographics periods (--tree periods) and other types of growth (--tree growth, --tree birth, --tree death)\n\n");
                    exit (-1);
                    }
                if (numPeriods <= 0)
                    {
                    fprintf (stderr, "\n!!! PARAMETER ERROR: Bad number of periods (%d)\n\n", numPeriods);
                    PrintUsage(stderr);
                    }
                Nbegin = 	(int *) calloc(numPeriods+1, sizeof(int));
                Nend =		(int *) calloc(numPeriods+1, sizeof(int));
                cumDuration =	(int *) calloc(numPeriods+1, sizeof(int));
                periodGrowth =	(double *) calloc(numPeriods+1, sizeof(double));
                if (Nbegin == NULL || Nend == NULL || cumDuration == NULL)
                    {
                    fprintf (stderr, "\n!!! PARAMETER ERROR: Could not allocate demographic vectors (%lu)\n", (unsigned long) (numPeriods * sizeof(int)));
                    exit (-1);
                    }
                for (j=1; j<=numPeriods; j++)
                    {
                    argument = atof(NEXTARG(i));
                    Nbegin[j] = (int) argument;
                    argument = atof(NEXTARG(i));
                    Nend[j] = (int) argument;
                    argument = atof(NEXTARG(i));
                    cumDuration[j] = (int) argument + cumDuration[j-1];
                    }
                break;
            case 'k':
                transformingBranchLengthRatio = atof(argv[i]);
              if (transformingBranchLengthRatio < 0)
                    {
                    fprintf (stderr, "\n!!! PARAMETER ERROR: Bad transforming branch length (%f)\n\n", transformingBranchLengthRatio);
                    PrintUsage(stderr);
                    }
                break;
			case 'q':
                healthyTipBranchLengthRatio = atof(argv[i]);
                if (healthyTipBranchLengthRatio < 0)
                    {
                    fprintf (stderr, "\n!!! PARAMETER ERROR: Bad healthy tip branch length (%f)\n\n", healthyTipBranchLengthRatio);
                    PrintUsage(stderr);
                    }
                break;
             case 'i':
                alphaBranches = atof(argv[i]);
                    if (alphaBranches <= 0)
                    {
                    fprintf (stderr, "\n!!! PARAMETER ERROR: Bad gamma alphaBranches shape (%f)\n\n", alphaBranches);
                    PrintUsage(stderr);
                    }
                rateVarAmongLineages = YES;
                break;
			case 'N':
				argument = atof(argv[i]);
				if (argument < 1 || argument > numCells)
					{
					fprintf (stderr, "\n!!! PARAMETER ERROR: Bad number of tree rate switches (%d), it has to be a single number between 1 and the number of cells (%d)\n\n", (int) argument, numCells);
					PrintUsage(stderr);
					}

				if (rateVarAmongLineages == YES)
					{
					fprintf (stderr, "\n!!! PARAMETER ERROR: Tree rate switches (--tree rate-switches) cannot be simulated together with gamma rate variation among lineages (--tree lineage-shape)");
					PrintUsage(stderr);
					}

				numTreeRateSwitches = (int) argument;
				doTreeRateSwitches = YES;
				treeRateSwitchMultiplier = (double*) Checked (calloc(numTreeRateSwitches, sizeof(double)), "Could not allocate the treeRateSwitchMultiplier array\n");

				for (j=0; j<numTreeRateSwitches; j++)
					{
					argument = atof(NEXTARG(i));
					treeRateSwitchMultiplier[j] = argument;
					if (treeRateSwitchMultiplier[j] <= 0)
						{
						fprintf (stderr, "\n!!! PARAMETER ERROR: Bad tree rate switch multiplier (%f), it has to be > 0 \n\n", treeRateSwitchMultiplier[j]);
						PrintUsage(stderr);
						}
					}
				break;
			case 'b':
                argument = atof(argv[i]);
                alphabet = (int) argument;
                   if (alphabet < 0 || alphabet > 1)
                    {
                    fprintf (stderr, "\n!!! PARAMETER ERROR: Bad alphabet (%d)\n\n", alphabet);
                    PrintUsage(stderr);
                    }
               break;
			case 'u':
                mutationRate = atof(argv[i]);
                if (mutationRate < 0)
                    {
                    fprintf (stderr, "\n!!! PARAMETER ERROR: Bad mutation rate (%f)\n\n", mutationRate);
                    PrintUsage(stderr);
                    }
                break;
            case 'd':
                deletionRate = atof(argv[i]);
                if (deletionRate < 0)
                    {
                    fprintf (stderr, "\n!!! PARAMETER ERROR: Bad deletion rate (%f)\n\n", deletionRate);
                    PrintUsage(stderr);
                    }
                break;
            case 'H':
                CNLOHrate = atof(argv[i]);
                if (CNLOHrate < 0)
                    {
                    fprintf (stderr, "\n!!! PARAMETER ERROR: Bad CNLOH rate (%f)\n\n", CNLOHrate);
                    PrintUsage(stderr);
                    }
                break;
			case 'j':
                argument = atof(argv[i]);
                numFixedMutations = (int) argument;
                if (numFixedMutations < 1)
                    {
                    fprintf (stderr, "\n!!! PARAMETER ERROR: Bad number of fixed mutations (%d)\n\n", numFixedMutations);
                    PrintUsage(stderr);
                    }
                doSimulateFixedNumMutations = YES;
                if (propAltModelSites > 0)
                    {
                    fprintf (stderr, "\n!!! PARAMETER ERROR: You cannot specify a fixed number of mutations if there is any non-ISM site. Set the proportion of non-ISM sites to zero\n\n");
					PrintUsage(stderr);
                    }
                break;
			case 'S':
                argument = atof(argv[i]);
                numUserSignatures = (int) argument;
                if (numUserSignatures < 1 || numUserSignatures > 30)
                    {
                    fprintf (stderr, "\n!!! PARAMETER ERROR: Bad  number of genetic signatures (%d), it has to be a single number between 1 and 30\n\n", numUserSignatures);
                    PrintUsage(stderr);
                    }
                doGeneticSignatures = YES;

				if (alphabet == BINARY)
					{
					fprintf (stderr, "\n!!! PARAMETER ERROR: Trinucleotide signatures cannot be simulated with the binary alphabet");
					PrintUsage(stderr);
					}

				signatureWeight = (double*) Checked (calloc(numUserSignatures, sizeof(double)), "Could not allocate the signatureWeight array\n");

				signatureID = (int*) Checked (calloc(numUserSignatures, sizeof(int)), "Could not allocate the signatureID array\n");
				sum = 0;
				for (j=0; j<numUserSignatures; j++)
					{
					argument = atof(NEXTARG(i));
					signatureID[j] = (int) argument;
					argument = atof(NEXTARG(i));
					signatureWeight[j] = argument;
					if (signatureWeight[j] <= 0)
						{
						fprintf (stderr, "\n!!! PARAMETER ERROR: Bad signatureWeight (%f)\n\n", signatureWeight[j]);
						PrintUsage(stderr);
						}
					sum += signatureWeight[j];
					}
				if (sum != 1.0)
					for (j=0; j<numUserSignatures; j++)
						signatureWeight[j] /= sum;
				    break;
			case 'm':
                argument = atof(argv[i]);
                altModel = (int) argument;
                if (altModel < 0 || altModel > 2)
                    {
                    fprintf (stderr, "\n!!! PARAMETER ERROR: Bad alternative mutation model (%d)\n\n", altModel);
                    PrintUsage(stderr);
                    }
				if (alphabet == DNA && propAltModelSites > 0)
					{
					if (altModel == Mk)
						{
						fprintf (stderr, "\n!!! PARAMETER ERROR: The DNA alphabet and the model (%d) specified are incompatible", (int) argument);
						PrintUsage(stderr);
						}
					}
				else if (alphabet == BINARY && propAltModelSites > 0)
					{
					if (altModel == finiteDNA)
						{
						fprintf (stderr, "\n!!! PARAMETER ERROR: The binary alphabet and the model (%d) specified are incompatible", (int) argument);
						PrintUsage(stderr);
						}
					}
               break;
		   case 'p':
                propAltModelSites = atof(argv[i]);
                if (propAltModelSites < 0 || propAltModelSites > 1)
                    {
                    fprintf (stderr, "\n!!! PARAMETER ERROR: Bad proportion of alternative model sites (%f). It has to be between 0 and 1\n\n", propAltModelSites);
                    PrintUsage(stderr);
                    }
				if (propAltModelSites > 0 && altModel != ISMhap && doSimulateFixedNumMutations == YES)
                    {
                    fprintf (stderr, "\n!!! PARAMETER ERROR: You cannot specify a proportion of non-ISM  sites bigger than zero if the number of mutations is fixed\n\n");
                    PrintUsage(stderr);
                    }
 				if (alphabet == DNA && propAltModelSites > 0)
					{
					if (altModel == Mk)
						{
						fprintf (stderr, "\n!!! PARAMETER ERROR: The DNA alphabet and the alt model (%d) specified are incompatible", (int) altModel);
						PrintUsage(stderr);
						}
					}
				else if (alphabet == BINARY && propAltModelSites > 0)
					{
					if (altModel == finiteDNA)
						{
						fprintf (stderr, "\n!!! PARAMETER ERROR: The binary alphabet and the alt model (%d) specified are incompatible", (int) altModel);
						PrintUsage(stderr);
						}
					}
               break;
            case 'w':
                nonISMRelMutRate = atof(argv[i]);
//               if (nonISMRelMutRate < 0 || nonISMRelMutRate > 1)
                 if (nonISMRelMutRate <= 0)
                    {
                    fprintf (stderr, "\n!!! PARAMETER ERROR: Bad relative nonISM/ISM mutation rate (%f)\n\n", nonISMRelMutRate);
                    PrintUsage(stderr);
                    }
                break;
			case 'c':
                SNPrate = atof(argv[i]);
                if (SNPrate < 0)
                    {
                    fprintf (stderr, "\n!!! PARAMETER ERROR: Bad germline SNP rate (%f)\n\n", SNPrate);
                    PrintUsage(stderr);
                    }
                break;
			case 'f':
                freq[0] = atof(argv[i]);
                freq[1] = atof(NEXTARG(i));
                freq[2] = atof(NEXTARG(i));
                freq[3] = atof(NEXTARG(i));
					
                if ((freq[0] == freq[1]) && (freq[0] == freq[2]) && (freq[0] == freq[3])) /* TODO fix with loop */
                    equalBaseFreq = YES;
                else
                    equalBaseFreq = NO;
					
                sumPi = freq[0] + freq[1] + freq[2] + freq[3];
                if (sumPi !=1.0)
                    {
                    freq[0]/=sumPi;
                    freq[1]/=sumPi;
                    freq[2]/=sumPi;
                    freq[3]/=sumPi;
                    }
                break;
			case 't':
	        	titv = atof(argv[i]);
				if (titv < 0)
					{
					fprintf (stderr, "\n!!! PARAMETER ERROR: Bad ti/tv (%f)\n\n", titv);
					PrintUsage(stderr);
					}
				else
					{
					doJC = NO;
					doHKY = YES;
					doGTR = NO;
					doGTnR = NO;
					}
				if (thereIsMij == YES)
					{
					fprintf (stderr, "\n!!! PARAMETER ERROR: Cannot specify a mutation matrix (GTR model) and a ti/tv (HKY model) at the same time\n\n");
					PrintUsage(stderr);
					}
				break;
			case 'a':
                alphaSites = atof(argv[i]);
                if (alphaSites < 0)
                    {
                    fprintf (stderr, "\n!!! PARAMETER ERROR: Bad gamma alphaSites shape (%f)\n\n", alphaSites);
                    PrintUsage(stderr);
                    }
                rateVarAmongSites = YES;
                break;
			case 'r':
				if (doHKY == YES)
					{
					fprintf (stderr, "\n!!! PARAMETER ERROR: Cannot specify a mutation matrix (GTR model) and a ti/tv (HKY model) at the same time\n\n");
					PrintUsage(stderr);
					}
					/*	AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT */
					Mij[0][0]= atof(argv[i]);
					Mij[0][1] = atof(NEXTARG(i));
					Mij[0][2] = atof(NEXTARG(i));
					Mij[0][3] = atof(NEXTARG(i));
					Mij[1][0] = atof(NEXTARG(i));
					Mij[1][1] = atof(NEXTARG(i));
					Mij[1][2] = atof(NEXTARG(i));
					Mij[1][3] = atof(NEXTARG(i));
					Mij[2][0] = atof(NEXTARG(i));
					Mij[2][1] = atof(NEXTARG(i));
					Mij[2][2] = atof(NEXTARG(i));
					Mij[2][3] = atof(NEXTARG(i));
					Mij[3][0] = atof(NEXTARG(i));
					Mij[3][1] = atof(NEXTARG(i));
					Mij[3][2] = atof(NEXTARG(i));
					Mij[3][3] = atof(NEXTARG(i));
					if (Mij[0][0] != 0  || Mij[1][1] != 0 || Mij[2][2] != 0 || Mij[3][3] != 0)
						{
						fprintf(stderr, "\n!!! PARAMETER ERROR: Bad general rate matrix: diagonals should be 0 \n\n");
						PrintUsage(stderr);
						}
					thereIsMij = YES;
					if (CheckMatrixSymmetry (Mij) == YES)
						{
						doJC = NO;
						doHKY = NO;
						doGTR = YES;
						doGTnR = NO;
						}
					else
						{
						doJC = NO;
						doHKY = NO;
						doGTR = NO;
						doGTnR = YES;
						}
				break;
   			case 'G':
				meanGenotypingError = atof(argv[i]);
				varGenotypingError = atof(NEXTARG(i));
				if (meanGenotypingError < 0 || meanGenotypingError > 1)
					{
					fprintf(stderr, "\n!!! PARAMETER ERROR: Bad mean genotyping error (%f)\n\n", meanGenotypingError);
					PrintUsage(stderr);
					}
				if (varGenotypingError <= 0 || (meanGenotypingError > 0 && varGenotypingError >= (meanGenotypingError * (1.0 - meanGenotypingError))))
					{
					fprintf(stderr, "\n!!! PARAMETER ERROR: Bad variance genotyping error (%f); it has to be > 0 and < mean*(1-mean)\n\n", varGenotypingError);
					PrintUsage(stderr);
					}
				doNGS = NO;
				break;
			case 'C':
				coverage = atof(argv[i]);
				if (!(coverage > 0))
					{
					fprintf (stderr, "\n!!! PARAMETER ERROR: Bad coverage (%f)\n\n", coverage);
					PrintUsage(stderr);
					}
				if (coverage > 0)
					doNGS = YES;
				if (meanGenotypingError > 0 && doNGS == YES)
					{
					fprintf (stderr, "\n!!! PARAMETER ERROR: Cannot specify a coverage larger than 0, which implies read count generation, and a genotyping error at the same time\n\n");
					PrintUsage(stderr);
					}
				break;
			case 'J':
				alphaCellCoverage = atof(argv[i]);
				if (alphaCellCoverage <= 0)
					{
					fprintf (stderr, "\n!!! PARAMETER ERROR: Bad between-cell coverage dispersion (%f)\n\n", alphaCellCoverage);
					PrintUsage(stderr);
					}
				break;
			case 'O':
				alphaSiteCoverage = atof(argv[i]);
				if (alphaSiteCoverage <= 0)
					{
					fprintf (stderr, "\n!!! PARAMETER ERROR: Bad between-site coverage dispersion (%f)\n\n", alphaSiteCoverage);
					PrintUsage(stderr);
					}
				break;
            case 'V':
               alphaCoverage = atof(argv[i]);
				if (alphaCoverage <= 0)
                    {
                    fprintf (stderr, "\n!!! PARAMETER ERROR: Bad coverage dispersion (%f)\n\n", alphaCoverage);
                    PrintUsage(stderr);
                    }
                rateVarCoverage = YES;
				break;
   			case 'A':
					meanAmplificationError = atof(argv[i]);
					varAmplificationError = atof(NEXTARG(i));
					simulateOnlyTwoTemplates = (int) atof(NEXTARG(i));
				if (meanAmplificationError < 0 || meanAmplificationError > 1)
					{
					fprintf(stderr, "\n!!! PARAMETER ERROR: Bad mean amplification error (%f)\n\n", meanAmplificationError);
					PrintUsage(stderr);
					}
				if (varAmplificationError <= 0 || (meanAmplificationError > 0 && varAmplificationError >= (meanAmplificationError * (1.0 - meanAmplificationError))))
					{
					fprintf(stderr, "\n!!! PARAMETER ERROR: Bad variance amplification error (%f); it has to be > 0 and < mean*(1-mean)\n\n", varAmplificationError);
					PrintUsage(stderr);
					}
				if (simulateOnlyTwoTemplates < 0 || simulateOnlyTwoTemplates > 1)
					{
					fprintf(stderr, "\n!!! PARAMETER ERROR: Bad simulateOnlyTwoTemplates error (%d); it has to be 0 (assume 2 templates) or 1 (assume 4 templates)", simulateOnlyTwoTemplates);
					PrintUsage(stderr);
					}
				break;
            case 'E':
                sequencingError = atof(argv[i]);
				if (sequencingError < 0 || sequencingError > 1)
                    {
                    fprintf (stderr, "\n!!! PARAMETER ERROR: Bad sequencing error (%f)\n\n", sequencingError);
                    PrintUsage(stderr);
                    }
                break;
			case 'D':
              fixedADOrate = atof(argv[i]);
				if (fixedADOrate < 0 || fixedADOrate > 1)
                    {
                    fprintf (stderr, "\n!!! PARAMETER ERROR: Bad fixed allelic dropout (%f)\n\n", fixedADOrate);
                    PrintUsage(stderr);
                    }
                if (meanADOcell > 0 || meanADOsite > 0)
					{
					fprintf (stderr, "\n!!! PARAMETER ERROR: Fixed ADO cannot be specified at the same time as variable ADO\n\n");
					PrintUsage(stderr);
					}
				if (fixedADOrate > 0)
					thereisADO = YES;
				break;
			case 'P':
				meanADOcell = atof(argv[i]);
				varADOcell = atof(NEXTARG(i));
				if (meanADOcell < 0 || meanADOcell == 0 || meanADOcell > 1)
					{
					fprintf(stderr, "\n!!! PARAMETER ERROR: Bad mean ADO per cell (%f)\n\n", meanADOcell);
					PrintUsage(stderr);
					}
				if (varADOcell <= 0 || (meanADOcell > 0 && varADOcell >= (meanADOcell * (1.0 - meanADOcell))))
					{
					fprintf(stderr, "\n!!! PARAMETER ERROR: Bad variance ADO per cell (%f); it has to be > 0 and < mean*(1-mean)\n\n", varADOcell);
					PrintUsage(stderr);
					}
				if (fixedADOrate > 0)
					{
					fprintf(stderr, "\n!!! PARAMETER ERROR: ADO has to be either fixed or variable per cell/site, cannot be both");
					PrintUsage(stderr);
					}
				doADOcell = YES;
				thereisADO = YES;
				break;
			case 'Q':
				meanADOsite = atof(argv[i]);
				varADOsite = atof(NEXTARG(i));
				if (meanADOsite < 0 || meanADOsite > 1)
					{
					fprintf(stderr, "\n!!! PARAMETER ERROR: Bad mean ADO per site (%f)\n\n", meanADOsite);
					PrintUsage(stderr);
					}
				if (varADOsite <= 0 || (meanADOsite > 0 && varADOsite >= (meanADOsite * (1.0 - meanADOsite))))
					{
					fprintf(stderr, "\n!!! PARAMETER ERROR: Bad variance ADO per site (%f); it has to be > 0 and < mean*(1-mean)\n\n", varADOsite);
					PrintUsage(stderr);
					}
				if (fixedADOrate > 0)
					{
					fprintf(stderr, "\n!!! PARAMETER ERROR: ADO has to be either fixed or variable per cell/site, cannot be both");
					PrintUsage(stderr);
					}
				thereisADO = YES;
				doADOsite = YES;
				break;
			case 'I':
					meanAllelicImbalance = atof(argv[i]);
					varAllelicImbalance = atof(NEXTARG(i));
				if (meanAllelicImbalance < 0 || meanAllelicImbalance > 1)
					{
					fprintf(stderr, "\n!!! PARAMETER ERROR: Bad mean allelic imbalance (%f)\n\n", meanAllelicImbalance);
					PrintUsage(stderr);
					}
				if (varAllelicImbalance <= 0 || (meanAllelicImbalance > 0 && varAllelicImbalance >= (meanAllelicImbalance * (1.0 - meanAllelicImbalance))))
					{
					fprintf(stderr, "\n!!! PARAMETER ERROR: Bad variance allelic imbalance (%f); it has to be > 0 and < mean*(1-mean)\n\n", varAllelicImbalance);
					PrintUsage(stderr);
					}
				doAllelicImbalance = YES;
				break;
			case 'R':
				haploidCoverageReduction = atof(argv[i]);
				if (haploidCoverageReduction < 0 || haploidCoverageReduction > 1)
                    {
                    fprintf (stderr, "\n!!! PARAMETER ERROR: Bad haploid coverage reduction (%f)\n\n", haploidCoverageReduction);
                    PrintUsage(stderr);
                    }
                break;
			case 'B':
					meanDoubletRate = atof(argv[i]);
					varDoubletRate = atof(NEXTARG(i));
				if (meanDoubletRate < 0 || meanDoubletRate > 1)
					{
					fprintf(stderr, "\n!!! PARAMETER ERROR: Bad mean doublet rate (%f)\n\n", meanDoubletRate);
					PrintUsage(stderr);
					}
				if (varDoubletRate <= 0 || (meanDoubletRate > 0 && varDoubletRate >= (meanDoubletRate * (1.0 - meanDoubletRate))))
					{
					fprintf(stderr, "\n!!! PARAMETER ERROR: Bad variance doublet rate (%f); it has to be > 0 and < mean*(1-mean)\n\n", varAllelicImbalance);
					PrintUsage(stderr);
					}
				break;
			case 'X':
					/*	AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT */
					Eij[0][0]= atof(argv[i]);
					Eij[0][1] = atof(NEXTARG(i));
					Eij[0][2] = atof(NEXTARG(i));
					Eij[0][3] = atof(NEXTARG(i));
					Eij[1][0] = atof(NEXTARG(i));
					Eij[1][1] = atof(NEXTARG(i));
					Eij[1][2] = atof(NEXTARG(i));
					Eij[1][3] = atof(NEXTARG(i));
					Eij[2][0] = atof(NEXTARG(i));
					Eij[2][1] = atof(NEXTARG(i));
					Eij[2][2] = atof(NEXTARG(i));
					Eij[2][3] = atof(NEXTARG(i));
					Eij[3][0] = atof(NEXTARG(i));
					Eij[3][1] = atof(NEXTARG(i));
					Eij[3][2] = atof(NEXTARG(i));
					Eij[3][3] = atof(NEXTARG(i));
					if (Eij[0][0] != 0  || Eij[1][1] != 0 || Eij[2][2] != 0 || Eij[3][3] != 0)
						{
						fprintf(stderr, "\n!!! PARAMETER ERROR: Bad error matrix: diagonals should be 0 \n\n");
						PrintUsage(stderr);
						}
				thereIsEij = YES;
				break;
            case '1':
                doPrintSNVgenotypes = YES;
                break;
            case '2':
                doPrintSNVhaplotypes = YES;
                break;
            case '3':
                doPrintFullGenotypes = YES;
                break;
            case '4':
                doPrintFullHaplotypes = YES;
                break;
            case '5':
                doPrintAncestors = YES;
                break;
            case '6':
                doPrintTree = YES;
                break;
            case '7':
                doPrintTimes = YES;
                break;
            case '8':
				doPrintCATG = YES;
   				if (doNGS == NO)
					{
					fprintf (stderr, "\n!!! PARAMETER ERROR: Cannot print the CATG format when coverage is <= 0\n\n");
					PrintUsage(stderr);
					}
               break;
			case '9':
				doPrintPILEUP = YES;
				if (doNGS == NO)
					{
					fprintf (stderr, "\n!!! PARAMETER ERROR: Cannot print the PILEUP format when coverage is <= 0\n\n");
					PrintUsage(stderr);
					}
				break;
           case 'Y':
                doPrintTrueHaplotypes = YES;
                break;
 			case 'v':
                doPrintSeparateReplicates = YES;
                break;
			case 'x':
                doPrintIUPAChaplotypes = YES;
                break;
 			case 'o':
				ch = *argv[i];
				if(!isspace(ch))
					CopyName (resultsDir, argv[i], "resultsDir");
			break;
 			case 'T':
				ch = *argv[i];
				if(!isspace(ch))
					CopyName (userTreeFile, argv[i], "userTreeFile");
				doUserTree = YES;
			break;
 			case 'U':
				ch = *argv[i];
				if(!isspace(ch))
					CopyName (userGenomeFile, argv[i], "userGenomeFile");
				doUserGenome = YES;
			break;
			case 'Z':
				doPrintVCF = NO;
				break;
			case 'z':
				doZipVCF = YES;
				break;
			case 'M':
				minVAFvcf = atof(argv[i]);
				if (minVAFvcf < 0 || minVAFvcf > 1)
					{
					fprintf (stderr, "\n!!! PARAMETER ERROR: Bad minimum VAF (%f)\n\n", minVAFvcf);
					PrintUsage(stderr);
					}
				break;
			case '0':
				doPrintBulk = YES;
				if (doNGS == NO)
					{
					fprintf (stderr, "\n!!! PARAMETER ERROR: Cannot print pooled read counts when coverage is <= 0\n\n");
					PrintUsage(stderr);
					}
				break;
			case 'W':
				doTumorNames = YES;
			break;
			case 'y':
                argument = atof(argv[i]);
                noisy = (int) argument;
                if (noisy < 0)
                    {
                    fprintf (stderr, "\n!!! PARAMETER ERROR: Bad noisy value (%d)\n\n", noisy);
                    PrintUsage(stderr);
                    }
                break;
            /*case 'z':
                argument = atof(argv[i]);
                numNodes = (int) argument;
                if (numNodes<1)
                    {
                    fprintf (stderr, "\n!!! PARAMETER ERROR: Bad number of nodes (%d)\n\n", numNodes);
                    PrintUsage(stderr);
                    }
                break;*/
            case '#':
                argdouble = atof(argv[i]);
                userSeed = (long int) argdouble;
                if (userSeed < 0)
                    {
                    fprintf (stderr, "\n!!! PARAMETER ERROR: Bad user seed for random number (%ld)\n\n", userSeed);
                    PrintUsage(stderr);
                    }
                break;
			case '?':
                PrintHelp (stdout);
                exit (0);
                break;
            default :
                if (flag == '-')		/* a long option that does not exist: say which, and the closest one */
                    {
                    const char *wrong = argv[i] - 2, *similar = SimilarLongOption (wrong);

                    fprintf (stderr, "\n!!! PARAMETER ERROR: Unknown option: %s", wrong);
                    if (similar)
                        fprintf (stderr, "  (did you mean %s?)", similar);
                    fprintf (stderr, "\n\n");
                    }
                else
                    fprintf (stderr, "\n!!! PARAMETER ERROR: Incorrect parameter: %c\n\n", flag);
                PrintUsage(stderr);
                break;
        }
    }
}


/***************************** ReadParametersFromFile *******************************/
/* A parameter file (-F, or the file "parameters") holds the options of the command line without the dash: "n10 s10 l1000 A0.1 0.01 0 ...",
   with text between brackets ignored. It is turned into the argument vector of the equivalent command line and read by the same parser,
   so a file and a command line cannot disagree. Options of the command line that take several numbers (h, N, S, f, r, G, A, P, Q, I, B, X) take
   the following words of the file (-h, -N and -S as many as their first number says). */
static int ExtraWords (const char *word)
	{
	int		n = atoi (word + 1);		/* the number of periods (-h), of switches (-N) or of signatures (-S) */

	switch (word[0])
		{
		case 'h': return 3 * n;
		case 'N': return n;
		case 'S': return 2 * n;
		case 'A': return 2;
		case 'f': return 3;
		case 'r': case 'X': return 15;
		case 'G': case 'P': case 'Q': case 'I': case 'B': return 1;
		}
	return 0;
	}

void ReadParametersFromFile (void)
	{
	size_t	cap = 1 << 16, n = 0, nw = 0, k;
	int		ch, argc = 1, i, depth = 0;
	char	*text = (char*) malloc (cap), **word = NULL, **argv, *p;

	if (feof (stdin))
		{
		fprintf (stderr, "FILE ERROR: Unable to read parameters from stdin\n");
		exit (-1);
		}
	while ((ch = fgetc (stdin)) != EOF)		/* the file without the comments in brackets */
		{
		if (n + 2 >= cap)
			text = (char*) realloc (text, cap *= 2);
		if (ch == '[')
			depth = 1;
		else if (ch == ']' && depth)
			{
			depth = 0;
			ch = ' ';
			}
		if (!depth)
			text[n++] = (char) ch;
		}
	if (depth)
		{
		fprintf (stderr, "FILE ERROR: The closing bracket is missing in the parameter file\n");
		exit (-1);
		}
	text[n] = '\0';
	for (p = strtok (text, " \t\r\n\v\f"); p != NULL; p = strtok (NULL, " \t\r\n\v\f"))
		{
		word = (char**) realloc (word, (nw + 1) * sizeof(char*));
		word[nw++] = p;
		}
	argv = (char**) calloc (nw + 2, sizeof(char*));
	argv[0] = (char*) "cellcoal";
	for (k=0; k<nw; )
		{
		int		extra = ExtraWords (word[k]);
		size_t	len = strlen (word[k]) + 2;
		char	*arg;

		arg = (char*) malloc (len);
		snprintf (arg, len, "-%s", word[k++]);
		argv[argc++] = arg;
		for (i=0; i<extra && k<nw; i++)
			argv[argc++] = word[k++];
		}
	paramFileArgs = (char*) calloc (MAX_LINE, 1);
	for (i=1; i<argc; i++)		/* the equivalent command line, for the log */
		{
		size_t used = strlen (paramFileArgs);

		if (used + strlen (argv[i]) + 2 >= MAX_LINE)
			break;
		snprintf (paramFileArgs + used, MAX_LINE - used, " %s", argv[i]);
		}
	ReadParametersFromCommandLine (argc, argv);
	}

/* default values of all the options */
void SetDefaultOptions (void)
	{
    numDataSets = 50;			/* the number of samples to simulate */
    numCells = 20;				/* number of ingroup cells in each data set */
    ploidy = 2;                 /* we assume diploid genomes */
    numSites = 1000;			/* number of sites to simulate = length of the genome */
    N = 10000;					/* effective population size */
    numPeriods = 0;				/* number of distinct demographic periods */
    doDemographics = NO;		/* whether to implement demographics */
    doExponential = NO;			/* whether to do exponential growth */
    growthRate = 1.0e-3;		/* rate for the exponential population growth */
	birthRate = 1;				/* birth rate in the Ohtuski and Innan 2017 coalescent model*/
	deathRate = 0;				/* death rate in the Ohtuski and Innan 2017 coalescent model*/
   	mutationRate = 1.0e-7;		/* nucleotide mutation rate per site per generation */
	SNPrate = 0.0;				/* germline variation rate for rooth healthy genome */
	rateVarAmongLineages = NO;	/* modify rate variation among branches (to non-clock) */
	doTreeRateSwitches = NO;	/* do rate switches along the tree */
	numTreeRateSwitches = 0;	/* number of rate switches along the tree */
	alphabet = DNA;          	/* alphabet 0/1 or DNA */
    altModel = 0;				/* by default the alternative model will be ISM haploid */
    propAltModelSites = 0;		/* proportion of sites that will mutate according to alternative model */
    nonISMRelMutRate = 1.0;		/* relative rate alternative/default model for sites */
    equalBaseFreq = YES;        /* DNA base frequencies */
	freq[0] = freq[1] = freq[2] = freq[3] = 0.25;
	titv = 0.5;					/* transition/transversion rate ratio */
	thereIsMij = NO;			/* mutation rate matrix*/
	Mij[0][0] = Mij[1][1] = Mij[2][2] = Mij[3][3] = 0;  /* mutation probabilities */
	Mij[0][1] = Mij[0][2] = Mij[0][3] = 1.0;
	Mij[1][0] = Mij[1][2] = Mij[1][3] = 1.0;
	Mij[2][0] = Mij[2][1] = Mij[2][3] = 1.0;
	Mij[3][0] = Mij[3][1] = Mij[3][2] = 1.0;
	thereIsEij = NO;			/* error rate matrix*/
	Eij[0][0] = Eij[1][1] = Eij[2][2] = Eij[3][3] = 0;  /* sequencing error probabilities */
	Eij[0][1] = Eij[0][2] = Eij[0][3] = 1.0;
	Eij[1][0] = Eij[1][2] = Eij[1][3] = 1.0;
	Eij[2][0] = Eij[2][1] = Eij[2][3] = 1.0;
	Eij[3][0] = Eij[3][1] = Eij[3][2] = 1.0;
	doJC = YES;
	doHKY = NO;
	doGTR = NO;
	doGTnR = NO;
	rateVarAmongSites = NO;     	/* rate variation among different sites along the genome */
    alphaSites = infinity;      	/* alpha shape of the gamma distribution for rate variation among sites */
    alphaBranches = infinity;   	/* alpha shape of the gamma distribution for rate variation among lineages */
    alphaCellCoverage = 0;			/* between-cell coverage dispersion (-J): 0 = none */
    alphaSiteCoverage = 0;			/* between-site coverage dispersion (-O): 0 = none */
    alphaCoverage = infinity;   	/* alpha shape of the gamma distribution for coverage */
	doSimulateFixedNumMutations = NO;	/* whether to simulate a fixed number of mutations */
	doUserTree = NO;				/* whether to assume a user tree instead od making the coalescent */
	doUserGenome = NO;				/* whether to use a user genome instead of a simulated one */
    doPrintSNVgenotypes = NO;		/* whether to print SNVs */
    doPrintSNVhaplotypes = NO;  	/* whether to print haplotypes */
    doPrintTrueHaplotypes = NO;  	/* whether to print the true haplotypes s */
    doPrintFullHaplotypes = NO;		/* whether to print sequences */
    doPrintFullGenotypes = NO;		/* whether to print all genotypes (variable + invariable) */
    doPrintTree = NO;				/* whether to print the coalescent tree */
    doPrintTimes = NO;				/* whether to print coalescent times */
	doPrintAncestors = NO;      	/* whether to print data for ancestral cells */
	doNGS = NO;  					/* produce reads by default */
	doTumorNames = NO;				/* use specific name for taxa in the tumor scenario */
	stringPrecision = 12;			/* precision for taxa names */
	doPrintCATG = NO;				/* whether to print read counts for SNVs in CATG format*/
	doPrintPILEUP = NO;				/* whether to print read counts for SNVs in PILEUP format*/
	doPrintBulk = NO;				/* whether to print pooled (bulk) read counts and VAF per site */
	doPrintVCF = YES;				/* whether to print the VCF file when reads are simulated */
	doZipVCF = NO;					/* whether the VCF is written compressed with gzip */
	minVAFvcf = 0.0;				/* minimum pooled VAF of the sites in the VCF */
	doSimulateData = YES;			/* whether to simulate any data (or just look at the expectations; useful for debugging) */
	doPrintSeparateReplicates = NO; /* whether to put every replica in its own file */
	doPrintIUPAChaplotypes = NO;	/* whether to print IUPAC halotypes */
	doGeneticSignatures = NO;		/* whether to use a genetic signature to model trinucleotide mutations */
	numUserSignatures = 0;			/* by default we do not use a genetic signature */
	transformingBranchLengthRatio = 0.5; /* ratio of the transforming branch length compare to the tumor MRCA length */
	healthyTipBranchLengthRatio = 1.0; 	/* ratio the branch leading to the healthy cell to the (tumor MRCA + transfoming branch) length*/
 	coverage = 5;					/* NGS  depth for read counts */
	rateVarCoverage = NO;			/* there is coverage dispersion */
	fixedADOrate = 0;				/* allelic dropout */
	thereisADO = NO;				/* allelic dropout flag */
	doADOcell = NO;					/* ADO does not change per cell by default */
	doADOsite = NO;					/* ADO does not change per site by default */
	meanADOcell = 0;     			/* mean of beta distribution for ADO expected for a cell */
	varADOcell = 1.0e-06;     		/* var of beta distribution for ADO expected for a cell */
	meanADOsite = 0;				/* mean of beta distribution for ADO expected for a site */
	varADOsite = 1.0e-06; ;			/* var of beta distribution for ADO expected for a site */
	sequencingError = 0;			/* NGS error rate */
	meanGenotypingError = 0;		/* mean of beta distribution for errors added directly in the genotypes */
	varGenotypingError = 1.0e-06;   /* variance of beta distribution for genotype errors */
	meanAmplificationError = 0; 	/* mean of beta distribution for WGA errors */
	varAmplificationError = 1.0e-06; /* variance of beta distribution for WGA errors */
	simulateOnlyTwoTemplates = NO;	/* whether simulating maximum of two templates after single-cell amplification, or there can be all four */
	haploidCoverageReduction = 0.5; /* proportion of reads produced when a single allele is present */
	doAllelicImbalance = NO;		/* by default there is no allelic imbalance */
	meanAllelicImbalance = 0.5;		/* beta mean proportion of maternal/ paternal reads */
	varAllelicImbalance = 1.0e-06;	/* beta var proportion of maternal/ paternal reads */
	meanDoubletRate = 0;			/* mean of beta distribution for doublet errors */
	varDoubletRate = 1.0e-06;		/* variance of beta distribution for doublet errors */
	numNodes = 1000;            	/* initial number of nodes allocated to build the coalescent trees */
	doSpecificParameterFile = NO;	/* by default do not use a custom parameter file */
	userSeed = 0;					/* seed entered by the user */
	noisy = 1;						/* level of information to be printed in the screen (see below) */
	}
