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
#define VERSION_NUMBER          "version 1.4.0"

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



/* Genotype likelihoods of a cell and site: after the phases are combined the 4x4 log10 likelihood matrix is
   symmetric, so only the 10 values with a1<=a2 are kept, plus the maximum; the normalized matrix is always
   logGL - max, so it is rebuilt with the same arithmetic when needed. A "gl" is an array of GL_NUM doubles.
   For big runs the likelihoods are not stored at all: they are a pure function of the data in CellSiteStr and
   are recomputed when they are printed (see GetSiteGL). */
#define GL_NUM		11
#define GL_MAX		10
static inline int GL_Index (int a1, int a2)
	{
	int lo = a1 < a2 ? a1 : a2, hi = a1 < a2 ? a2 : a1;
	return lo * (7 - lo) / 2 + hi;
	}
#define GLV(gl,a1,a2)	((gl)[GL_Index((a1),(a2))])
#define GLN(gl,a1,a2)	((gl)[GL_Index((a1),(a2))] - (gl)[GL_MAX])

/* 48 bytes per cell and site (alleles fit in a signed char: values are 0-3, DELETION=7, ADO=9, MISSING=-9) */
typedef struct csite
    {
	int			numReads;
	int			readCount[4];
	double  	maternalSiteAmplificationError, paternalSiteAmplificationError;
	signed char	trueMaternalAllele, truePaternalAllele;
	signed char	maternalAllele, paternalAllele;
	signed char	thereIsMaternalAllele, thereIsPaternalAllele;
	signed char	MLmatAllele, MLpatAllele;
	}
    CellSiteStr;

/* data for cells that are doublets; only allocated for those cells */
typedef struct
    {
	int		numReadsDoublet;
	int		readCountDoublet[4];
	int		maternalAlleleDoublet, paternalAlleleDoublet;
	double	logGLdoublet[4][4];
	double	normLogGLdoublet[4][4];
	int		MLmatAlleleDoublet;
	int		MLpatAlleleDoublet;
	}
    DoubletSiteStr;


typedef struct
    {
	struct csite	*site;
	DoubletSiteStr	*doublet;	/* NULL unless hasDoublet == YES */
	int				hasDoublet;
	}
    CellStr;

#endif /* macros_h */
