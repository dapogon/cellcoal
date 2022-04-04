//
//  macros.h
//  CellCoal
//
//  Created by David Posada on 11/04/2018.
//  Copyright © 2018 David Posada. All rights reserved.
//

#ifndef macros_h
#define macros_h

#define PROGRAM_NAME            "CellCoal"
#define PROGRAM_NAME_UPPERCASE	"CELLCOAL"
#define PROGRAM_NAME_LOWERCASE	"cellcoal"
#define VERSION_NUMBER          "version 1.3.1"

#define	MAX_NAME                120
#define	MAX_LINE                500000
#define	MAX_GENOME              1000000
#define	NO                      0
#define	YES                     1
#define	infinity                99999
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
#define MISSING					-9

#define CG_AT					0 	// C>A or G>T
#define CG_GC					1	// C>G or G>C
#define CG_TA					2	// C>T or G>A
#define TA_AT					3	// T>A or A>T
#define TA_CG					4	// T>C or A>G
#define TA_GC					5	// T>G or A>C

#define ISMhap					0
#define Mk2						1
#define finiteDNA				2

#define NUM_SIGNATURES			30

#define trinuc(i,j,k)   (i*16 + j*4 + k)
#define trimut(i,j,k)   (i*24 + j*4 + k)

#define CHECK_MUT_DISTRIBUTION
#undef CHECK_MUT_DISTRIBUTION

#define MYDEBUG
#undef MYDEBUG

#define LOAD_INTERNAL_SIGNATURES

#define PRINT_TRIMUTATIONS
#undef PRINT_TRIMUTATIONS

#define PRINT_TRINUC_GENOME
#undef PRINT_TRINUC_GENOME

#define PRINT_TRIMUTCOUNTER
#undef PRINT_TRIMUTCOUNTER

#define COUNT_ML_GENOTYPE_ERRORS
#undef COUNT_ML_GENOTYPE_ERRORS


typedef struct tnode
    {
    struct tnode	*left, *right, *anc;
	int				index, label, isHealthyTip, isHealthyRoot;
	int				numLeftDescendants, numRightDescendants;
    double			length, time, branchLength, meanDepth;
	double			rateSwitch;
	char			*name;
    }
    TreeNode;

typedef struct
    {
    int		isSNV, isSNP, isVariant;
	int		numMutations, numMutationsMaternal, numMutationsPaternal;
    int     numDeletions, numDeletionsMaternal, numDeletionsPaternal;
    int     numCNLOH, numCNLOHmaternal, numCNLOHpaternal;
	int		hasADO;
	int		hasGenotypeError;
	int		countA, countC, countG, countT, countN, countACGT; // alelle counts in genotypes
	int		readCountA, readCountC, readCountG, readCountT, readCountACGT; // alelle counts in reads
	int		countCellswithData, countDropped;
	int		referenceAllele;
	int		*alternateAlleles;
	int		numAltAlleles;
    double  rateMultiplier;
	}
    SiteStr;


typedef struct
    {
	int 	tempLength;
	int 	numAvailablePositions;
	int 	*position;
	}
    TriNucStr;


typedef struct csite
    {
	int		numReads;
	int		*readCount;
	int		trueMaternalAllele, truePaternalAllele;
	int		maternalAllele, paternalAllele;
	int		thereIsMaternalAllele, thereIsPaternalAllele;
	double  maternalSiteAmplificationError, paternalSiteAmplificationError;
	double	**logGL;
	double	**normLogGL;
	int		MLmatAllele;
	int		MLpatAllele;
	int		numReadsDoublet;
	int		*readCountDoublet;
	int		maternalAlleleDoublet, paternalAlleleDoublet;
	double	**logGLdoublet;
	double	**normLogGLdoublet;
	int		MLmatAlleleDoublet;
	int		MLpatAlleleDoublet;
	}
    CellSiteStr;


typedef struct
    {
	struct csite	*site;
	int				hasDoublet;
	}
    CellStr;

#endif /* macros_h */
