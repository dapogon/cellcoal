#ifndef readopts_h
#define readopts_h
/* Long options of the read-level simulator (--bam, --fastq ...) and of the chromosome in blocks (--blocks, --fasta ...): X(key, option, flag, kind, variable, default)
   kind 0 = switch, 1 = integer, 2 = number, 3 = text. The key is the name used in a configuration file (dashes -> underscores). */
#define READ_OPTS(X) \
	X(bam,            "--bam",            '\200', 0, doBam,       0)    \
	X(fastq,          "--fastq",          '\201', 0, doFastq,     0)    \
	X(read_len,       "--read-len",       '\202', 1, readLen,     150)  \
	X(insert_mean,    "--insert-mean",    '\203', 2, insertMean,  400)  \
	X(insert_sd,      "--insert-sd",      '\204', 2, insertSd,    60)   \
	X(amplicon_len,   "--amplicon-len",   '\205', 2, ampliconLen, 2000) \
	X(amplicon_cv,    "--amplicon-cv",    '\206', 2, ampliconCV,  0.5)  \
	X(bulk_depth,     "--bulk-depth",     '\207', 2, bulkDepth,   0)    \
	X(contig,         "--contig",         '\210', 3, contigName,  0)    \
	X(dup_mean,       "--dup-mean",       '\211', 2, dupMean,     0.3)    \
	X(blocks,         "--blocks",         '\212', 1, blocksN,     0)      \
	X(fasta,          "--fasta",          '\213', 3, fastaFile,   0)      \
	X(region,         "--region",         '\214', 3, regionSpec,  0)      \
	X(germline_vcf,   "--germline-vcf",   '\215', 3, germlineVcf, 0)
#endif
