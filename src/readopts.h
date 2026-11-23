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
	X(germline_pop,   "--germline-pop",   '\234', 3, germlinePop, 0)      \
	X(qual_sd,        "--qual-sd",        '\216', 2, qualSd,     0)      \
	X(qual_slope,     "--qual-slope",     '\217', 2, qualSlope,  0)      \
	X(strand_bias,    "--strand-bias",    '\220', 2, strandBias, 0)    \
	X(seg_loss_rate,  "--seg-loss-rate",      '\221', 2, lossRate,   0)      \
	X(seg_len,        "--seg-len",       '\222', 2, segLen,    1000)  \
	X(seg_gain_rate,  "--seg-gain-rate",  '\223', 2, gainRate,   0)    \
	X(seg_loss_len,   "--seg-loss-len",   '\224', 2, segLossLen, 0)    \
	X(seg_gain_len,   "--seg-gain-len",   '\225', 2, segGainLen, 0)    \
	X(seg_max_copies, "--seg-max-copies", '\226', 1, maxCopies,  4)    \
	X(phase,          "--phase",          '\227', 0, doPhase,    0)    \
	X(compress_level, "--compress-level", '\230', 1, compressLevel, 6) \
	X(region_rates,   "--region-rates",   '\231', 3, regionRatesFile, 0) \
	X(seq_context,    "--seq-context",    '\232', 3, seqContextFile, 0) \
	X(damage,         "--damage",         '\233', 3, damageSpec, 0)
#endif
