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
	X(bulk_depth,     "--bulk-depth",     '\207', 2, bulkDepth,   0)    \
	X(contig,         "--contig",         '\210', 3, contigName,  0)    \
	X(dup_mean,       "--dup-mean",       '\211', 2, dupMean,     0.3)    \
	X(blocks,         "--blocks",         '\212', 1, blocksN,     0)      \
	X(fasta,          "--fasta",          '\213', 3, fastaFile,   0)      \
	X(region,         "--region",         '\214', 3, regionSpec,  0)      \
	X(germline_vcf,   "--germline-vcf",   '\215', 3, germlineVcf, 0)      \
	X(qual_sd,        "--qual-sd",        '\216', 2, qualSd,     0)      \
	X(qual_slope,     "--qual-slope",     '\217', 2, qualSlope,  0)      \
	X(strand_bias,    "--strand-bias",    '\220', 2, strandBias, 0)    \
	X(loss_rate,      "--loss-rate",      '\221', 2, lossRate,   0)      \
	X(loss_len,       "--loss-len",       '\222', 2, lossLen,    1000)
#endif
