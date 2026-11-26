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
	X(deme_n,         "--deme-n",         '\235', 3, demeN,       0)      \
	X(deme_size,      "--deme-size",      '\236', 3, demeSize,    0)      \
	X(deme_source,    "--deme-source",   '\237', 3, demeSource,  0)      \
	X(deme_time,      "--deme-time",     '\240', 3, demeTime,    0)      \
	X(deme_founders,  "--deme-founders", '\241', 3, demeFounders, 0)     \
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
	X(damage,         "--damage",         '\233', 3, damageSpec, 0)      \
	X(indel_rate,     "--indel-rate",     '\242', 2, indelRate,  0)      \
	X(indel_ins,      "--indel-ins",      '\243', 2, indelIns,   0.5)    \
	X(indel_len,      "--indel-len",      '\244', 2, indelLen,   2)      \
	X(indel_max,      "--indel-max",      '\245', 1, indelMax,   10)      \
	X(homopolymer_error, "--homopolymer-error", '\246', 2, homoError, 0) \
	X(indel_tail,     "--indel-tail",     '\247', 2, indelTail,  0)      \
	X(indel_tail_exp, "--indel-tail-exp", '\250', 2, indelTailExp, 2)    \
	X(indel_repeat,   "--indel-repeat",   '\251', 2, indelRepeat, 1)      \
	X(germline_indel_rate, "--germline-indel-rate", '\252', 2, germIndelRate, 0) \
	X(mappability,    "--mappability",    '\253', 3, mappabilityFile, 0) \
	X(map_drop,       "--map-drop",       '\254', 2, mapDrop,    0)       \
	X(targets,        "--targets",        '\255', 3, targetsFile, 0)    \
	X(off_target,     "--off-target",     '\256', 2, offTarget,  0.05) \
	X(clone_n,        "--clone-n",        '\257', 1, cloneN,     0)      \
	X(clone_time,     "--clone-time",     '\260', 3, cloneTime,  0)      \
	X(clone_parent,   "--clone-parent",   '\261', 3, cloneParent, 0)     \
	X(clone_fitness,  "--clone-fitness",  '\262', 3, cloneFitness, 0)    \
	X(clone_founders, "--clone-founders", '\263', 3, cloneFounders, 0)  \
	X(clone_cells,    "--clone-cells",    '\264', 2, cloneCells0, 0)    \
	X(clone_birth,    "--clone-birth",    '\265', 2, cloneBirth,  0)     \
	X(clone_death,    "--clone-death",    '\266', 2, cloneDeath,  0)     \
	X(clone_sizes,    "--clone-sizes",    '\267', 3, cloneSizes,  0)     \
	X(clone_freq,     "--clone-freq",     '\270', 3, cloneFreq,   0)
#endif
