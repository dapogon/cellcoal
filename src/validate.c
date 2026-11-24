#include "cellcoal.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "definitions.h"

/* Checks of the options: the modes, the combinations, the read-level options and the chromosome (blocks, regions) */

/* checks that need several options at once (after reading the command line) */
void CheckOptionCombinations (void)
	{
	if (doBCF == YES)
		{
		if (doZipVCF == YES || doPrintVCF == NO)
			{
			fprintf (stderr, "\n!!! PARAMETER ERROR: --bcf replaces the VCF file: it cannot be combined with --output vcf-gz or --output no-vcf\n\n");
			exit (-1);
			}
		if (numDataSets > 1 && doPrintSeparateReplicates == NO)
			{
			fprintf (stderr, "\n!!! PARAMETER ERROR: --bcf with several replicates needs --output replicate-folders (one BCF file per replicate)\n\n");
			exit (-1);
			}
		}
	if (ampCellVar > 0.0 && (meanAmplificationError <= 0.0 || ampCellVar >= meanAmplificationError * (1.0 - meanAmplificationError)))
		{
		fprintf (stderr, "\n!!! PARAMETER ERROR: --amp-cell-var needs an amplification error (--amplification error) and a variance < mean*(1-mean) of --amplification error\n\n");
		exit (-1);
		}
	if (covAdoCorr != 0.0 && (alphaCellCoverage <= 0.0 || doADOcell == NO))
		{
		fprintf (stderr, "\n!!! PARAMETER ERROR: --cov-ado-corr needs coverage variation among cells (--depth cells) and ADO variation among cells (--ado cells)\n\n");
		exit (-1);
		}
	if (covDoubletCorr != 0.0 && (alphaCellCoverage <= 0.0 || meanDoubletRate <= 0))
		{
		fprintf (stderr, "\n!!! PARAMETER ERROR: --cov-doublet-corr needs coverage variation among cells (--depth cells) and doublets (--depth doublets)\n\n");
		exit (-1);
		}
	if (adoDoubletCorr != 0.0 && (doADOcell == NO || meanDoubletRate <= 0))
		{
		fprintf (stderr, "\n!!! PARAMETER ERROR: --ado-doublet-corr needs ADO variation among cells (--ado cells) and doublets (--depth doublets)\n\n");
		exit (-1);
		}
	}

/* The three ways CellCoal produces data: genotypes (no -C), read counts (-C) and reads (--bam, --fastq) */
int RunMode (void)
	{
	if (doBam || doFastq)
		return 2;
	return (doNGS == YES) ? 1 : 0;
	}

/* options of a mode that the other modes do not use: an error in reads mode, a warning without -C */
void CheckModeOptions (void)
	{
	int		reads = (RunMode () == 2), counts = (RunMode () == 1), genotypes = (RunMode () == 0);
	struct { const char *name; int given; int countsOnly; } o[] =
		{
		{"-I (allelic imbalance)", doAllelicImbalance == YES, 1},
		{"-V (coverage overdispersion)", rateVarCoverage == YES, 1},
		{"-O (coverage overdispersion between sites)", alphaSiteCoverage > 0.0, 1},
		{"-R (haploid coverage reduction)", haploidCoverageReduction != 0.5, 1},
		{"-M (minimum pooled VAF)", minVAFvcf > 0.0, 1},
		{"--fast", fastMode != 0, 1},
		{"--min-alt-reads", minAltReads > 0, 1},
		{"--min-alt-vaf", minAltVaf > 0.0, 1},
		{"--coverage", coverageWindow > 0, 1},
		{"-8, -9, -0, -Z, -z, --bcf (outputs of read counts)", doPrintCATG == YES || doPrintPILEUP == YES || doPrintBulk == YES || doBCF == YES || doZipVCF == YES, 1},
		{NULL, 0, 0}
		};
	struct { const char *name; int given; } both[] =
		{
		{"-J (coverage overdispersion between cells)", alphaCellCoverage > 0.0},
		{"-A (amplification error)", meanAmplificationError > 0.0},
		{"-E (sequencing error)", sequencingError > 0.0},
		{"-X (error matrix)", thereIsEij == YES},
		{"--gain-sd/--gain-len", gainSd > 0.0 || gainLen > 0.0},
		{"--amp-signature", ampSignatureFile != NULL},
		{NULL, 0}
		};
	struct { const char *name; int given; } readsOnlyArg[] =
		{
		{"--read-len", ExtraOptionGiven ("--read-len")},
		{"--insert-mean", ExtraOptionGiven ("--insert-mean")},
		{"--insert-sd", ExtraOptionGiven ("--insert-sd")},
		{"--dup-mean", ExtraOptionGiven ("--dup-mean")},
		{"--bulk-depth", ExtraOptionGiven ("--bulk-depth")},
		{"--qual-sd", ExtraOptionGiven ("--qual-sd")},
		{"--qual-slope", ExtraOptionGiven ("--qual-slope")},
		{"--strand-bias", ExtraOptionGiven ("--strand-bias")},
		{"--seq-context", ExtraOptionGiven ("--seq-context")},
		{"--damage", ExtraOptionGiven ("--damage")},
		{"--indel-rate", ExtraOptionGiven ("--indel-rate")},
		{"--indel-ins", ExtraOptionGiven ("--indel-ins")},
		{"--indel-len", ExtraOptionGiven ("--indel-len")},
		{"--indel-max", ExtraOptionGiven ("--indel-max")},
		{"--indel-tail", ExtraOptionGiven ("--indel-tail")},
		{"--indel-tail-exp", ExtraOptionGiven ("--indel-tail-exp")},
		{"--indel-repeat", ExtraOptionGiven ("--indel-repeat")},
		{"--germline-indel-rate", ExtraOptionGiven ("--germline-indel-rate")},
		{"--mappability", ExtraOptionGiven ("--mappability")},
		{"--map-drop", ExtraOptionGiven ("--map-drop")},
		{"--homopolymer-error", ExtraOptionGiven ("--homopolymer-error")},
		{NULL, 0}
		};
	int		k;

	if (reads)
		for (k=0; o[k].name; k++)
			if (o[k].given)
				{
				fprintf (stderr, "\n!!! PARAMETER ERROR: %s is an option of the read counts (-C) and is not used with --bam/--fastq (reads mode): the coverage varies with --depth cells and --gain-sd\n\n", o[k].name);
				exit (-1);
				}
	if (!reads)
		for (k=0; readsOnlyArg[k].name; k++)
			if (readsOnlyArg[k].given)
				{
				fprintf (stderr, "\n!!! PARAMETER ERROR: %s is an option of the reads mode and needs --bam and/or --fastq\n\n", readsOnlyArg[k].name);
				exit (-1);
				}
	if (genotypes)
		{
		for (k=0; o[k].name; k++)
			if (o[k].given)
				fprintf (stderr, "\nNote: %s is ignored without -C (genotypes mode: no reads are simulated)\n", o[k].name);
		for (k=0; both[k].name; k++)
			if (both[k].given)
				fprintf (stderr, "\nNote: %s is ignored without -C (genotypes mode: no reads are simulated)\n", both[k].name);
		}
	(void) counts;
	}

void CheckReadOptions (void)
	{
	if (gainSd > 0.0 || gainLen > 0.0)
		{
		if (!(gainSd > 0.0 && gainLen > 0.0))
			{
			fprintf (stderr, "\nERROR: --gain-sd and --gain-len go together\n");
			exit (-1);
			}
		if (doNGS == NO)
			{
			fprintf (stderr, "\nERROR: --gain-sd needs read counts (-C)\n");
			exit (-1);
			}
		}

	if (ampSignatureFile != NULL)
		{
		if (doNGS == NO || meanAmplificationError <= 0.0)
			{
			fprintf (stderr, "\nERROR: --amp-signature needs read counts (-C) and an amplification error (--amplification error)\n");
			exit (-1);
			}
		if (alphabet == BINARY)
			{
			fprintf (stderr, "\nERROR: --amp-signature needs the DNA alphabet\n");
			exit (-1);
			}
		SetupAmpSignature (ampSignatureFile);
		}

	if (adoBlockLen > 0.0 && (doADOsite == YES || (fixedADOrate <= 0 && doADOcell == NO)))
		{
		fprintf (stderr, "\nERROR: --ado-block needs allelic dropout (--ado mean or --ado cells) and cannot be combined with ADO among sites (--ado sites)\n");
		exit (-1);
		}

	/* read-level simulator: --bam, --fastq */
	if (doBam || doFastq)
		{
		if (doNGS == NO || alphabet == BINARY || ploidy != 2)
			{
			fprintf (stderr, "\nERROR: --bam/--fastq need the mean coverage (-C) and diploid DNA\n");
			exit (-1);
			}
		if ((fixedADOrate > 0 || doADOcell == YES || doADOsite == YES) && adoBlockLen <= 0.0)
			{
			fprintf (stderr, "\nERROR: in reads mode the allelic dropout has to come in segments (--ado-block L, e.g. 3000): with dropout independent at every site almost no fragment of %d bp survives\n", 2 * readLen);
			exit (-1);
			}
		if (readLen < 20 || readLen > 400 || insertMean < readLen || insertSd < 0.0 || dupMean < 0.0)
			{
			fprintf (stderr, "\nERROR: bad --read-len/--insert-mean/--insert-sd/--dup-mean\n");
			exit (-1);
			}
		if (qualSd < 0.0 || qualSd > 3.0 || qualSlope < -6.0 || qualSlope > 6.0 || strandBias < 0.0 || strandBias > 1.0)
			{
			fprintf (stderr, "\nERROR: bad --qual-sd (0 to 3), --qual-slope (-6 to 6) or --strand-bias (0 to 1)\n");
			exit (-1);
			}
		if (strlen (seqContextFile) > 0 && sequencingError <= 0.0)
			{
			fprintf (stderr, "\nERROR: --seq-context needs a sequencing error (-E) greater than 0\n");
			exit (-1);
			}
		if (indelRate < 0.0 || indelIns < 0.0 || indelIns > 1.0 || indelLen < 1.0 || indelMax < 1 || indelMax > 100 || indelTail < 0.0 || indelTail > 1.0 || indelTailExp < 0.0 || indelRepeat < 0.0 || germIndelRate < 0.0)
			{
			fprintf (stderr, "\nERROR: bad --indel-rate (0 or more), --indel-ins (0 to 1), --indel-len (1 or more), --indel-max (1 to 100), --indel-tail (0 to 1), --indel-tail-exp, --indel-repeat (0 or more) or --germline-indel-rate (0 or more)\n");
			exit (-1);
			}
		if (mapDrop < 0.0 || mapDrop > 1.0 || (mapDrop > 0.0 && strlen (mappabilityFile) == 0) || (mapDrop > 0.0 && doBam == NO))
			{
			fprintf (stderr, "\nERROR: --map-drop has to be from 0 to 1 and needs --mappability and --bam (the FASTQ files are the reads before the alignment)\n");
			exit (-1);
			}
		if (homoError < 0.0 || homoError > 0.1)
			{
			fprintf (stderr, "\nERROR: bad --homopolymer-error (0 to 0.1)\n");
			exit (-1);
			}
		/* the VCF of read counts and genotype likelihoods does not exist in this mode */
		doPrintVCF = NO;
		}
	if (compressLevel < 1 || compressLevel > 9)
		{
		fprintf (stderr, "\nERROR: --compress-level has to be between 1 and 9\n");
		exit (-1);
		}
	}

/* the options that need a chromosome (blocks, a real reference, regions of rates), and the setup of the chromosome */
void CheckChromosomeOptions (void)
	{
	if (coverage > 1e6)
		{
		fprintf (stderr, "\nERROR: the mean coverage -C (%g) is too large: the maximum is 1000000\n", coverage);
		exit (-1);
		}

	/* single-nucleotide deletions (-d) are a model of isolated sites: with a chromosome (blocks, real reference) or reads, deletions are small indels */
	if (deletionRate > 0.0 && (blocksN > 0 || strlen (fastaFile) > 0 || doBam == YES || doFastq == YES))
		{
		fprintf (stderr, "\nERROR: -d (single-nucleotide deletions of the sites) cannot be used with a chromosome in blocks, a reference (--fasta) or reads.\n");
		fprintf (stderr, "Use --indel-rate (with --indel-ins, --indel-len, --indel-max) to simulate small insertions and deletions along the chromosome.\n");
		exit (-1);
		}

	/* regional mutation rates need a chromosome with coordinates: a simulated one in blocks or a real one (--blocks, --fasta) */
	if (strlen (regionRatesFile) > 0)
		{
		if (blocksN == 0 && strlen (fastaFile) == 0)
			{
			fprintf (stderr, "\nERROR: --region-rates needs a chromosome with coordinates: use --blocks (simulated chromosome) or --fasta (real chromosome).\n");
			fprintf (stderr, "Without them the sites of the simulation have no position, so there are no regions.\n");
			exit (-1);
			}
		if (doGeneticSignatures == YES)
			{
			fprintf (stderr, "\nERROR: --region-rates cannot be combined with mutational signatures\n");
			exit (-1);
			}
		regionRatesOn = YES;
		rateVarAmongSites = YES;
		}

	/* chromosome in blocks (--blocks, --fasta): the reference and the number of blocks */
	if (blocksN > 0 || strlen (fastaFile) > 0)
		{
		if (doUserGenome == YES || (SNPrate > 0 && (strlen (germlineVcf) > 0 || strlen (germlinePop) > 0)))
			{
			fprintf (stderr, "\nERROR: --blocks/--fasta cannot be combined with --mutation genome, and --germline-vcf with --mutation germline-rate\n");
			exit (-1);
			}
		BlocksSetup ();
		if (strlen (fastaFile) > 0)
			{
			doUserGenome = YES;
			snprintf (userGenomeFile, sizeof(userGenomeFile), "%s", fastaFile);
			}
		if (adoSeedOpt == 0)
			adoSeedOpt = originalSeed + 1;
		}
	}
