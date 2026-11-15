/*
	CellCoal - streamed reads and the CATG, PILEUP, VCF and bulk outputs
	See cellcoal.h for the shared declarations. This file was split from cellcoal.c.
*/
#include "cellcoal.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "definitions.h"

/********************* PrintCATG  ************************/
/**	CATG output format with read counts
	http://nbviewer.jupyter.org/gist/dereneaton/d2fd4e70d29f5ee5d195/testing_cat.ipynb#View-the-.cat-results-files

	The first line has the number of samples and the number of sites.
	Following this is a transposed data matrix with the consensus base calls where each row is a site and each column a different sample.
	To the right of each site is a tab-separated list of counts of the four bases at that site in the order C,A,T,G.

	12 44500
	1A0	1B0	1C0	1D0	2E0	2F0	2G0	2H0	3I0
	YCCCCCCCCCCC	10,0,10,0	20,0,0,0	20,0,0,0	20,0,0, ...
	GGGGGGGGGGGG	0,0,0,20	0,0,0,20	0,0,0,20	0,0,0,20 ...
	AAAAAAAAAAAA	0,20,0,0	0,20,0,0	0,20,0,0	0,20,0,0 ...
	CCCCCCCCCCCC	20,0,0,0	20,0,0,0	20,0,0,0	20,0,0,0 ...
	AAAAAAAAAAAA	0,20,0,0	0,20,0,0	0, 20,0,0	0,20,0,0 ...
	...
*/
/* Maximum likelihood alleles of the SNV sites. When the reads are streamed the window of sites does not hold them
   any more, so they are kept in a compact store: one byte per cell and SNV (maternal allele in the low nibble,
   paternal in the high nibble, 15 = missing). */
void MLPairAt (int i, int q, int *m, int *p)
	{
	int j = SNVsites[q];

	if (streamMode == YES)
		{
		unsigned char b = mlStore[(long) q * (numCells+1) + i];

		*m = ((b & 15) == 15) ? MISSING : (int) (b & 15);
		*p = ((b >> 4) == 15) ? MISSING : (int) (b >> 4);
		}
	else if (cell[i].hasDoublet == NO)
		{
		*m = CS(i,j).MLmatAllele;
		*p = CS(i,j).MLpatAllele;
		}
	else
		{
		*m = CD(i,j).MLmatAlleleDoublet;
		*p = CD(i,j).MLpatAlleleDoublet;
		}
	}

/* ML alleles of the site j (the q-th SNV) of all the cells, taken from the window into the store */
static void StoreMLAlleles (int q, int j)
	{
	int		i, m, p;
	unsigned char *dst = mlStore + (long) q * (numCells+1);

	for (i=0; i<numCells+1; i++)
		{
		if (cell[i].hasDoublet == NO)
			{
			m = CS(i,j).MLmatAllele;
			p = CS(i,j).MLpatAllele;
			}
		else
			{
			m = CD(i,j).MLmatAlleleDoublet;
			p = CD(i,j).MLpatAlleleDoublet;
			}
		dst[i] = (unsigned char) (((m == MISSING) ? 15 : m) | (((p == MISSING) ? 15 : p) << 4));
		}
	}

void PrintCATGHeader (FILE *fp)
	{
	int	i;

	fprintf (fp,"%d %d\n",numCells+1, numSNVs);
	for (i=0; i<numCells; i++)
		fprintf (fp,"%s%04d ", inCellName,i+1);

	fprintf (fp,"%s ", outCellName);
	}

/* one line of the CATG file: starts with a newline, then the IUPAC codes and the read counts of the cells */
static void CATGSiteLine (FILE *fp, int j)
	{
	int	i;

	fprintf (fp,"\n");
	for (i=0; i<numCells+1; i++)
		fprintf (fp,"%c",WhichIUPAC(CS(i,j).MLmatAllele, CS(i,j).MLpatAllele));

		/* decide per cell (the original tested cell[numCells+1].hasDoublet, one element past the end of the array) */
	for (i=0; i<numCells+1; i++)
		{
		if (cell[i].hasDoublet == NO)
			fprintf (fp,"\t%d,%d,%d,%d", CS(i,j).readCount[C], CS(i,j).readCount[A], CS(i,j).readCount[T], CS(i,j).readCount[G]);
		else
			fprintf (fp,"\t%d,%d,%d,%d", CD(i,j).readCountDoublet[C], CD(i,j).readCountDoublet[A], CD(i,j).readCountDoublet[T], CD(i,j).readCountDoublet[G]);
		}
	}

void PrintCATG (FILE *fp)
	{
	int	snv;

	PrintCATGHeader (fp);
	for (snv=0; snv<numSNVs; snv++)
		CATGSiteLine (fp, SNVsites[snv]);
	fprintf (fp,"\n");
	}


/********************* PrintPILEUP  ************************/
/**	PILEUP output format with read counts
	http://samtools.sourceforge.net/pileup.shtml
 
 seq1 272 T 24  ,.$.....,,.,.,...,,,.,..^+. <<<+;<<<<<<<<<<<=<;<;7<&
 seq1 273 T 23  ,.....,,.,.,...,,,.,..A <<<;<<<<<<<<<3<=<<<;<<+
 seq1 274 T 23  ,.$....,,.,.,...,,,.,...    7<7;<;<<<<<<<<<=<;<;<<6
 seq1 275 A 23  ,$....,,.,.,...,,,.,...^l.  <+;9*<<<<<<<<<=<<:;<<<<
 seq1 276 G 22  ...T,,.,.,...,,,.,....  33;+<<7=7<<7<&<<1;<<6<
 seq1 277 T 22  ....,,.,.,.C.,,,.,..G.  +7<;<<<<<<<&<=<<:;<<&<
 seq1 278 G 23  ....,,.,.,...,,,.,....^k.   %38*<<;<7<<7<=<<<;<<<<<
 seq1 279 C 23  A..T,,.,.,...,,,.,..... ;75&<<<<<<<<<=<<<9<<:<<
 where each line consists of
 chromosome,		1-based coordinate, 		reference base, 		the number of reads covering the site, 		read bases  		base qualities.
 
 
 At the read base column, a dot stands for a match to the reference base on the forward strand, a comma for a match on the reverse strand, `ACGTN' for a mismatch on the forward strand and `acgtn' for a mismatch on the reverse strand. A pattern `\+[0-9]+[ACGTNacgtn]+' indicates there is an insertion between this reference position and the next reference position.
 
 The byte representing quality runs from 0x21 (lowest quality; '!' in ASCII) to 0x7e (highest quality; '~' in ASCII). Here are the quality value characters in left-to-right increasing order of quality (ASCII):

 */
static void PileupSite (FILE *fp, int j)
	{
	int		i, k;
	int		referenceAllele;

	/* PILEUP: CHROMOSOME */
	fprintf (fp,"%d", 1);
	
	/* PILEUP: POSITION */
	fprintf (fp,"\t%d", j+1);

	/* PILEUP: REFERENCE allele(s) => healthy root alleles */
	referenceAllele = allSites[j].referenceAllele;
	fprintf (fp,"\t%c", WhichNuc(referenceAllele));
	
	for (i=0; i<numCells+1; i++)
		{
		/* PILEUP: NUMBER OF READS */
		fprintf (fp, "\t%d\t", CS(i,j).numReads);
		
		/* PILEUP: READ BASES */
		if (allSites[j].countA > 0)
			{
			if (referenceAllele != A)
				for (k=0; k<CS(i,j).readCount[A]; k++)
					fprintf (fp, "A");
			else
				for (k=0; k<CS(i,j).readCount[A]; k++)
					fprintf (fp, ".");
			}
		
		if (allSites[j].countC > 0)
			{
			if (referenceAllele != C)
				for (k=0; k<CS(i,j).readCount[C]; k++)
					fprintf (fp, "C");
			else
				for (k=0; k<CS(i,j).readCount[C]; k++)
					fprintf (fp, ".");
			}

		if (allSites[j].countG > 0)
			{
			if (referenceAllele != G)
				for (k=0; k<CS(i,j).readCount[G]; k++)
					fprintf (fp, "G");
			else
				for (k=0; k<CS(i,j).readCount[G]; k++)
					fprintf (fp, ".");
			}

		if (allSites[j].countT > 0)
			{
			if (referenceAllele != T)
				for (k=0; k<CS(i,j).readCount[T]; k++)
					fprintf (fp, "T");
			else
				for (k=0; k<CS(i,j).readCount[T]; k++)
					fprintf (fp, ".");
			}

		fprintf (fp, "\t");
		
		/* PILEUP: READ QUALITIES */  /* put always highest quality */
		for (k=0; k<CS(i,j).numReads; k++)
			fprintf (fp, "~");
			}
	
	fprintf (fp,"\n");
	}

void PrintPILEUP (FILE *fp)
	{
	int		snv;

	for (snv=0; snv<numSNVs; snv++)
		PileupSite (fp, SNVsites[snv]);
	}



/********************* PrintVCF  ************************/
/**	VCF output format with read counts
*/

/********************* fast VCF formatting helpers ************************/
/* Output-identical replacements for the many tiny fprintf calls used to write each VCF genotype field */

static inline char *VCF_PutUInt (char *p, unsigned long long n)
	{
	char	t[24];
	int		k = 0;

	/* the values written are almost always small: skip the digit reversal for them */
	if (n < 10)
		{
		*p++ = (char)('0' + n);
		return p;
		}
	if (n < 100)
		{
		*p++ = (char)('0' + n / 10);
		*p++ = (char)('0' + n % 10);
		return p;
		}
	if (n < 1000)
		{
		*p++ = (char)('0' + n / 100);
		*p++ = (char)('0' + (n / 10) % 10);
		*p++ = (char)('0' + n % 10);
		return p;
		}
	do
		{
		t[k++] = (char)('0' + n % 10);
		n /= 10;
		}
	while (n);
	while (k)
		*p++ = t[--k];
	return p;
	}

static inline char *VCF_PutInt (char *p, int v)
	{
	if (v < 0)
		{
		*p++ = '-';
		return VCF_PutUInt (p, (unsigned long long)(-(long long)v));
		}
	return VCF_PutUInt (p, (unsigned long long)v);
	}

/* same text as sprintf("%.2f", x). Fast path when x*100 is not within 1e-6 of a rounding tie
   (the error of x*100 is < 1e-8 for |x| < 1e5); otherwise defer to the C library */
static char *VCF_PutF2 (char *p, double x)
	{
	double				ax = fabs (x), v, fl, fr;
	unsigned long long	n;

	if (!(ax < 1e5))
		return p + sprintf (p, "%.2f", x);
	v  = ax * 100.0;
	fl = (double)(unsigned long long) v;	/* v >= 0 and < 1e7, so this is floor (v) */
	fr = v - fl;
	if (fabs (fr - 0.5) < 1e-6)
		return p + sprintf (p, "%.2f", x);
	n = (unsigned long long)(fr > 0.5 ? fl + 1.0 : fl);
	if (signbit (x))
		*p++ = '-';
	p = VCF_PutUInt (p, n / 100);
	*p++ = '.';
	*p++ = (char)('0' + (n % 100) / 10);
	*p++ = (char)('0' + (n % 100) % 10);
	return p;
	}

/* a log10 likelihood; a value that is not finite (a likelihood of exactly zero, e.g. with no sequencing error) is written as '.', the VCF missing value */
static inline char *VCF_PutLL (char *p, double x)
	{
	if (!isfinite (x))
		{
		*p++ = '.';
		return p;
		}
	return VCF_PutF2 (p, x);
	}

#define VCF_ROOM(fp,buf,p)	do { if ((p) - (buf) > 3000) { fwrite ((buf), 1, (size_t)((p) - (buf)), (fp)); (p) = (buf); } } while (0)

/* writes the tab-started genotype fields (GT ... TG) of one cell at one SNV site */
static void VCFCellFields (FILE *fp, CellSiteStr *cs, const double *glp, DoubletSiteStr *ds, int hasDoublet, SiteStr *st, int referenceAllele, int npairs, const int *pa1, const int *pa2)
	{
	char	buf[4096], *p = buf;
	double	(*LG)[4], (*NLG)[4], LGs[4][4], NLGs[4][4];
	char	txt[4][4][400], txtN[4][4][400];	/* text of the G10 / G10N values (a %.2f of a huge double can be long) */
	int		txtLen[4][4], txtNLen[4][4];
	int		*rc, numReads, MLm, MLp, l, k, a1, a2;
	int		MLmatAllele = cs->MLmatAllele, MLpatAllele = cs->MLpatAllele;

	for (k=0; k<4; k++)
		for (l=0; l<4; l++)
			{
			LGs[k][l]  = GLV (glp, k, l);
			NLGs[k][l] = GLN (glp, k, l);
			}
	if (hasDoublet == NO)
		{
		LG = LGs;  NLG = NLGs;  rc = cs->readCount;  numReads = cs->numReads;
		MLm = cs->MLmatAllele;  MLp = cs->MLpatAllele;
		}
	else
		{
		LG = ds->logGLdoublet;  NLG = ds->normLogGLdoublet;  rc = ds->readCountDoublet;  numReads = ds->numReadsDoublet;
		MLm = ds->MLmatAlleleDoublet;  MLp = ds->MLpatAlleleDoublet;
		}

	/* GT */
	*p++ = '\t';
	if (MLmatAllele == MISSING)
		*p++ = '.';
	else if (MLmatAllele == referenceAllele)
		*p++ = '0';
	else
		for (l=0; l<st->numAltAlleles; l++)
			if (st->alternateAlleles[l] == MLmatAllele)
				{
				p = VCF_PutInt (p, l+1);
				break;
				}
	*p++ = '|';
	if (MLpatAllele == MISSING)
		*p++ = '.';
	else if (MLpatAllele == referenceAllele)
		*p++ = '0';
	else
		for (l=0; l<st->numAltAlleles; l++)
			if (st->alternateAlleles[l] == MLpatAllele)
				{
				p = VCF_PutInt (p, l+1);
				break;
				}

	/* DP, RC */
	*p++ = ':';
	p = VCF_PutInt (p, numReads);
	*p++ = ':';
	p = VCF_PutInt (p, rc[A]);  *p++ = ',';
	p = VCF_PutInt (p, rc[C]);  *p++ = ',';
	p = VCF_PutInt (p, rc[G]);  *p++ = ',';
	p = VCF_PutInt (p, rc[T]);

	/* G10: the text of each of the 10 values is kept, because GL repeats a subset of them */
	*p++ = ':';
	for (k=0; k<4; k++)
		for (l=k; l<4; l++)
			{
			if (!(k == 0 && l == 0))
				*p++ = ',';
			txtLen[k][l] = (int)(VCF_PutLL (txt[k][l], LG[k][l]) - txt[k][l]);
			memcpy (p, txt[k][l], (size_t)txtLen[k][l]);
			p += txtLen[k][l];
			VCF_ROOM (fp, buf, p);
			}
	/* G10N */
	*p++ = ':';
	for (k=0; k<4; k++)
		for (l=k; l<4; l++)
			{
			if (!(k == 0 && l == 0))
				*p++ = ',';
			txtNLen[k][l] = (int)(VCF_PutLL (txtN[k][l], NLG[k][l]) - txtN[k][l]);
			memcpy (p, txtN[k][l], (size_t)txtNLen[k][l]);
			p += txtNLen[k][l];
			VCF_ROOM (fp, buf, p);
			}
	/* GL */
	*p++ = ':';
	for (k=0; k<npairs; k++)
		{
		if (k > 0)
			*p++ = ',';
		if (pa1[k] <= pa2[k])
			{
			memcpy (p, txt[pa1[k]][pa2[k]], (size_t)txtLen[pa1[k]][pa2[k]]);
			p += txtLen[pa1[k]][pa2[k]];
			}
		else
			p = VCF_PutLL (p, LG[pa1[k]][pa2[k]]);
		VCF_ROOM (fp, buf, p);
		}
	/* GLN */
	*p++ = ':';
	for (k=0; k<npairs; k++)
		{
		if (k > 0)
			*p++ = ',';
		if (pa1[k] <= pa2[k])
			{
			memcpy (p, txtN[pa1[k]][pa2[k]], (size_t)txtNLen[pa1[k]][pa2[k]]);
			p += txtNLen[pa1[k]][pa2[k]];
			}
		else
			p = VCF_PutLL (p, NLG[pa1[k]][pa2[k]]);
		VCF_ROOM (fp, buf, p);
		}
	/* PL  (the non-doublet logGL is tested for infinity in both cases, as in the original) */
	*p++ = ':';
	for (k=0; k<npairs; k++)
		{
		a1 = pa1[k];  a2 = pa2[k];
		if (k > 0)
			*p++ = ',';
		if (!isfinite (LGs[a1][a2]) || !isfinite (LG[a1][a2]))
			*p++ = '.';	/* infinite likelihood: PL cannot be represented as an integer, so use the VCF missing value */
		else
			p = VCF_PutInt (p, (int)round (-10 * LG[a1][a2]));
		VCF_ROOM (fp, buf, p);
		}
	/* PLN */
	*p++ = ':';
	for (k=0; k<npairs; k++)
		{
		a1 = pa1[k];  a2 = pa2[k];
		if (k > 0)
			*p++ = ',';
		if (!isfinite (LGs[a1][a2]) || !isfinite (NLG[a1][a2]))
			*p++ = '.';	/* infinite likelihood: PL cannot be represented as an integer, so use the VCF missing value */
		else
			p = VCF_PutInt (p, (int)round (-10 * NLG[a1][a2]));
		VCF_ROOM (fp, buf, p);
		}
	/* ML */
	*p++ = ':';
	if (numReads > 0)
		{
		*p++ = WhichNuc (MLm);  *p++ = '/';  *p++ = WhichNuc (MLp);
		}
	else
		{
		*p++ = '.';  *p++ = '/';  *p++ = '.';
		}
	/* NG */
	*p++ = ':';
	*p++ = WhichNuc (MLmatAllele);  *p++ = '|';  *p++ = WhichNuc (MLpatAllele);
	/* DG */
	*p++ = ':';
	if (hasDoublet == NO)
		{
		*p++ = '.';  *p++ = '/';  *p++ = '.';
		}
	else
		{
		*p++ = WhichNuc (ds->maternalAlleleDoublet);  *p++ = '|';  *p++ = WhichNuc (ds->paternalAlleleDoublet);
		}
	/* TG */
	*p++ = ':';
	*p++ = WhichNuc (cs->trueMaternalAllele);  *p++ = '|';  *p++ = WhichNuc (cs->truePaternalAllele);

	fwrite (buf, 1, (size_t)(p - buf), fp);
	}

/* writes the complete VCF line of one SNV site; it only reads shared data, so lines of different sites can be built concurrently */
static void VCFSiteLine (FILE *fp, int j)
	{
	if (VCFSitePasses (j) == NO)
		return;
	int		i, k, l, referenceAllele;

		
	/* VFC: CHROMOSOME */
	fprintf (fp,"\n%d", 1);
	
	/* VFC: POSITION */
	fprintf (fp,"\t%d", j+1);
	
	/* VFC: ID */
	if (allSites[j].isSNV == NO)
		fprintf (fp,"\tinv%06d", j+1);
	else
		fprintf (fp,"\tsnv%06d", j+1);

	/* VFC: REFERENCE allele(s) => healthy root alleles */
	referenceAllele = allSites[j].referenceAllele;
	fprintf (fp,"\t%c", WhichNuc(referenceAllele));
	
	/* VFC: ALTERNATE allele(s) seen in the reads (so not necessarily true) */
	fprintf (fp,"\t");
	if (allSites[j].numAltAlleles == 0)
		fprintf (fp,".");
	else
		{
		fprintf (fp,"%c", WhichNuc(allSites[j].alternateAlleles[0]));
		for (k=1; k<allSites[j].numAltAlleles; k++)
			fprintf (fp,",%c", WhichNuc(allSites[j].alternateAlleles[k]));
		}

	/* VFC: QUALITY */  /* Qphred probability that a SNV exists at this site; I will put missing info here */
	fprintf (fp, "\t.");

	/* VFC: FILTER */
	if (allSites[j].countACGT < allSites[j].countN)  /* if less than half of the samples have missing data */
		fprintf (fp, "\ts50");
	else
		fprintf (fp, "\tPASS");

	/* VFC: INFO: AA (ancestral allele) ; ancestral allele is the reference alleles in this case */
	fprintf (fp, "\tAA=%c", WhichNuc(referenceAllele));

	/* VFC: INFO: NS (number of samples with data) */
	fprintf (fp, ";NS=%d", allSites[j].countCellswithData);

	/* VFC: INFO: AF (true alternate allele frequencies) */
	fprintf (fp, ";AF");
	{
	const char *afSep = "=";
	if (allSites[j].countA > 0 && referenceAllele != A)
		{ fprintf (fp,"%s%4.3f", afSep, (double) allSites[j].countA / allSites[j].countACGT); afSep = ","; }
	if (allSites[j].countC > 0 && referenceAllele != C)
		{ fprintf (fp,"%s%4.3f", afSep, (double) allSites[j].countC / allSites[j].countACGT); afSep = ","; }
	if (allSites[j].countG > 0 && referenceAllele != G)
		{ fprintf (fp,"%s%4.3f", afSep, (double) allSites[j].countG / allSites[j].countACGT); afSep = ","; }
	if (allSites[j].countT > 0 && referenceAllele != T)
		{ fprintf (fp,"%s%4.3f", afSep, (double) allSites[j].countT / allSites[j].countACGT); afSep = ","; }
	}

	/* VFC: INFO : SOMATIC */
	fprintf (fp, ";SOMATIC");
	
	/* VFC: FORMAT */
	fprintf (fp, "\tGT:DP:RC:G10:G10N:GL:GLN:PL:PLN:ML:NG:DG:TG");
	
	{
	int npairs = 0, pa1[10], pa2[10];
	for (k=0; k<=allSites[j].numAltAlleles; k++)
		for (l=0; l<=k; l++)
			{
			pa1[npairs] = (k == 0) ? allSites[j].referenceAllele : allSites[j].alternateAlleles[k-1];
			pa2[npairs] = (l == 0) ? allSites[j].referenceAllele : allSites[j].alternateAlleles[l-1];
			npairs++;
			}
	for (i=0; i<numCells+1; i++)
		{
		double glbuf[GL_NUM];
		const double *glp = GetSiteGL (i, j, glbuf);
		VCFCellFields (fp, &CS(i,j), glp, cell[i].doublet ? &CD(i,j) : NULL, cell[i].hasDoublet, &allSites[j], referenceAllele, npairs, pa1, pa2);
		}
	}
	}


static void PrintVCFHeader (FILE *fp)
	{
	int		i;

	/* print file header */
	fprintf (fp,"##fileformat=VCFv4.3");
	fprintf (fp,"\n##filedate=");
	PrintDate (fp);
	fprintf (fp,"##source=CellCoal SNV simulation");
	fprintf (fp,"\n##contig=<ID=1>");
	fprintf (fp,"\n##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral allele\">");
	fprintf (fp,"\n##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">");
	fprintf (fp,"\n##INFO=<ID=AF,Number=R,Type=Float,Description=\"True alternate/s allele frequency (ignores doublets)\">");
	fprintf (fp,"\n##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic SNV flag\">");
	fprintf (fp,"\n##FILTER=<ID=s50,Description=\"Less than half of samples have data\">");
	fprintf (fp,"\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Called ML genotype\">");
	fprintf (fp,"\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">");
	fprintf (fp,"\n##FORMAT=<ID=RC,Number=4,Type=Integer,Description=\"Read counts for ACGT\">");
	fprintf (fp,"\n##FORMAT=<ID=G10,Number=G,Type=Float,Description=\"Log10-scaled genotype likelihoods for all 10 genotypes (order is AA AC AG AT CC CG CT GG GT TT)\">");
	fprintf (fp,"\n##FORMAT=<ID=G10N,Number=G,Type=Float,Description=\"Normalized log10-scaled genotype likelihoods for all 10 genotypes (order is AA AC AG AT CC CG CT GG GT TT)\">");
	fprintf (fp,"\n##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Log10-scaled genotype likelihoods\">");
	fprintf (fp,"\n##FORMAT=<ID=GLN,Number=G,Type=Float,Description=\"Normalized log10-scaled genotype likelihoods\">");
	fprintf (fp,"\n##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">");
	fprintf (fp,"\n##FORMAT=<ID=PLN,Number=G,Type=Integer,Description=\"Normalized Phred-scaled genotype likelihoods\">");	fprintf (fp,"\n##FORMAT=<ID=ML,Number=1,Type=String,Description=\"Maximum likelihood genotype\">");
	fprintf (fp,"\n##FORMAT=<ID=NG,Number=1,Type=String,Description=\"Genotype in the NGS library (considers ADO/DEL; ignores doublets)\">");
	fprintf (fp,"\n##FORMAT=<ID=DG,Number=1,Type=String,Description=\"Genotype in the NGS library for the second cell if there is a doublet\">");
	fprintf (fp,"\n##FORMAT=<ID=TG,Number=1,Type=String,Description=\"True genotype (considers DEL; ignores ADO/doublets\">");
	fprintf (fp,"\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	
	if (doUserTree == NO)
		{
		for (i=0; i<numCells; i++)
			fprintf (fp,"\t%s%04d", inCellName, i+1);
		fprintf (fp,"\t%s", outCellName);
		}
	else
		{
		for (i=0; i<numCells; i++)
			fprintf (fp,"\t%s", cellNames[i]);
		fprintf (fp,"\t%s", outCellName);
		}
	}

void PrintVCF (FILE *fp)
	{
	int		snv;

	PrintVCFHeader (fp);

	/* print only apparent SNVs */
#ifdef _OPENMP
	if (omp_get_max_threads() > 1 && numSNVs > 1)
		{
		/* build the lines of a chunk of sites in memory in parallel, then write them in order */
		int		chunk = 4 * omp_get_max_threads(), q, nq;
		char	**lineBuf = (char**) calloc (chunk, sizeof(char*));
		size_t	*lineLen = (size_t*) calloc (chunk, sizeof(size_t));
		if (!lineBuf || !lineLen)
			{
			fprintf (stderr, "Could not allocate the VCF line buffers\n");
			exit (-1);
			}
		for (snv=0; snv<numSNVs; snv+=chunk)
			{
			nq = (numSNVs - snv < chunk) ? numSNVs - snv : chunk;
			#pragma omp parallel for schedule(dynamic,1)
			for (q=0; q<nq; q++)
				{
				FILE *ms = open_memstream (&lineBuf[q], &lineLen[q]);
				if (!ms)
					{
					fprintf (stderr, "Could not open a memory stream for the VCF line\n");
					exit (-1);
					}
				VCFSiteLine (ms, SNVsites[snv+q]);
				fclose (ms);
				}
			for (q=0; q<nq; q++)
				{
				fwrite (lineBuf[q], 1, lineLen[q], fp);
				free (lineBuf[q]);
				}
			}
		free (lineBuf);
		free (lineLen);
		}
	else
#endif
		for (snv=0; snv<numSNVs; snv++)
			VCFSiteLine (fp, SNVsites[snv]);
	
	fprintf (fp, "\n");
	}


/********************* BulkSite  ************************/
/* One line per site: reads of the four nucleotides pooled over the tumour cells, pooled depth, alternate reads and VAF
   (as in a bulk sample), and the true fraction of alternate alleles among the tumour cells. Only sites with alternate
   reads or with a true alternate allele are printed. */
static void BulkSite (FILE *fp, int j)
	{
	int		i, ref = allSites[j].referenceAllele;
	long	reads[4] = {0, 0, 0, 0}, depth, altReads, trueAlt = 0;
	double	vaf;

	for (i=0; i<numCells; i++)
		{
		if (cell[i].hasDoublet == YES)
			{
			reads[0] += CD(i,j).readCountDoublet[A];
			reads[1] += CD(i,j).readCountDoublet[C];
			reads[2] += CD(i,j).readCountDoublet[G];
			reads[3] += CD(i,j).readCountDoublet[T];
			}
		else
			{
			reads[0] += CS(i,j).readCount[A];
			reads[1] += CS(i,j).readCount[C];
			reads[2] += CS(i,j).readCount[G];
			reads[3] += CS(i,j).readCount[T];
			}
		trueAlt += (CS(i,j).trueMaternalAllele != ref) + (CS(i,j).truePaternalAllele != ref);
		}
	depth = reads[0] + reads[1] + reads[2] + reads[3];
	altReads = depth - reads[ref];
	if (altReads == 0 && trueAlt == 0)
		return;
	vaf = (depth > 0) ? (double) altReads / (double) depth : 0.0;
	fprintf (fp, "%d\t%c\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%.6f\t%ld\t%.6f\n", j+1, WhichNuc(ref), reads[0], reads[1], reads[2], reads[3],
			 depth, altReads, vaf, trueAlt, (double) trueAlt / (2.0 * numCells));
	}

void PrintBulk (FILE *fp)
	{
	int		j;

	for (j=0; j<numSites; j++)
		BulkSite (fp, j);
	}

/* inverse of the standard normal distribution function (Acklam's rational approximation, relative error < 1.2e-9) */
static double NormalQuantile (double p)
	{
	static const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00};
	static const double b[] = {-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01};
	static const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00};
	static const double d[] = {7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00};
	double q, r;

	if (p < 0.02425)
		{
		q = sqrt (-2.0 * log (p));
		return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
		}
	if (p > 1.0 - 0.02425)
		{
		q = sqrt (-2.0 * log (1.0 - p));
		return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
		}
	q = p - 0.5;
	r = q * q;
	return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1.0);
	}

static const double	*sortKey;
static int CmpByKey (const void *x, const void *y)
	{
	double a = sortKey[*(const int*)x], b = sortKey[*(const int*)y];

	if (a < b) return -1;
	if (a > b) return 1;
	return *(const int*)x - *(const int*)y;
	}

/* keeps the ADO rate of every cell (drawn with -P) for the correlation with the coverage (--cov-ado-corr) */
void KeepCellADO (const double *rate)
	{
	if (covAdoCorr == 0.0)
		return;
	free (cellAdoRate);
	cellAdoRate = (double*) malloc ((numCells+1) * sizeof(double));
	if (!cellAdoRate)
		{
		fprintf (stderr, "Could not allocate the cell ADO rates\n");
		exit (-1);
		}
	memcpy (cellAdoRate, rate, (numCells+1) * sizeof(double));
	}

/* Gaussian-copula style correlation between the ADO rate and the coverage multiplier of the cells (--cov-ado-corr rho).
   The coverage multipliers already drawn (gamma, -J) are reassigned to the cells so that their ranks follow
   s = rho z + sqrt(1-rho^2) e, where z are the normal scores of the ranks of the ADO rates (beta, -P) and e is standard
   normal noise. Only the assignment changes: the distributions of both quantities are exactly the requested ones.
   rho > 0: cells with more dropout have more coverage; rho < 0: cells with more dropout have less coverage. */
static void CorrelateCoverageWithADO (long int *seed)
	{
	int		i, n = numCells, *byAdo, *byScore;
	double	*z, *score, *cov, u1, u2, e;

	byAdo = (int*) malloc (n * sizeof(int));
	byScore = (int*) malloc (n * sizeof(int));
	z = (double*) malloc (n * sizeof(double));
	score = (double*) malloc (n * sizeof(double));
	cov = (double*) malloc (n * sizeof(double));
	if (!byAdo || !byScore || !z || !score || !cov)
		{
		fprintf (stderr, "Could not allocate the coverage-ADO correlation\n");
		exit (-1);
		}
	for (i=0; i<n; i++)
		byAdo[i] = byScore[i] = i;
	sortKey = cellAdoRate;
	qsort (byAdo, n, sizeof(int), CmpByKey);
	for (i=0; i<n; i++)
		z[byAdo[i]] = NormalQuantile ((i + 0.5) / n);
	for (i=0; i<n; i++)
		{
		u1 = RandomUniform (seed);
		u2 = RandomUniform (seed);
		e = sqrt (-2.0 * log (u1)) * cos (2.0 * M_PI * u2);
		score[i] = covAdoCorr * z[i] + sqrt (1.0 - covAdoCorr * covAdoCorr) * e;
		}
	sortKey = score;
	qsort (byScore, n, sizeof(int), CmpByKey);
	for (i=0; i<n; i++)
		cov[i] = cellCovMult[i];
	sortKey = cov;
	qsort (byAdo, n, sizeof(int), CmpByKey);	/* byAdo now holds the cells' coverage order (reused as scratch) */
	{
	/* k-th smallest coverage value goes to the cell with the k-th smallest score */
	double *sorted = (double*) malloc (n * sizeof(double));

	for (i=0; i<n; i++)
		sorted[i] = cov[byAdo[i]];
	for (i=0; i<n; i++)
		cellCovMult[byScore[i]] = sorted[i];
	free (sorted);
	}
	free (byAdo); free (byScore); free (z); free (score); free (cov);
	}

/********************* SetupCellCoverage  ************************/
/* Coverage varies between cells (-J): the depth expected in every cell is multiplied by a gamma variate with mean 1,
   drawn once per cell and replicate from the main random stream (the same in stored and streamed reads). */
void SetupCellCoverage (long int *seed)
	{
	int	i;

	free (cellCovMult);
	cellCovMult = (double*) malloc ((numCells+1) * sizeof(double));
	if (!cellCovMult)
		{
		fprintf (stderr, "Could not allocate the cell coverage multipliers\n");
		exit (-1);
		}
	for (i=0; i<numCells+1; i++)
		cellCovMult[i] = (alphaCellCoverage > 0.0) ? RandomGamma (alphaCellCoverage, seed) / alphaCellCoverage : 1.0;
	if (covAdoCorr != 0.0 && cellAdoRate != NULL && numCells > 1)
		CorrelateCoverageWithADO (seed);
	}

/********************* SetupDoublets  ************************/
/* Chooses which cells are doublets and their donor cells (from the main random stream), and allocates the window of sites
   of the doublet cells. The reads of the doublets are then simulated site by site, from the stream of each site. */
void SetupDoublets (long int *seed, long windowSites)
	{
	int		i;
	double	doubletRate;

	if (meanDoubletRate <= 0)
		return;
	doubletDonor = (int*) calloc (numCells, sizeof(int));
	for (i=0; i<numCells; i++)
		{
		doubletRate = RandomBetaMeanVar (meanDoubletRate, varDoubletRate, seed);
		if (RandomUniform (seed) < doubletRate)
			{
			cell[i].hasDoublet = YES;
			doubletDonor[i] = RandomUniformTo (numCells, seed);
			cell[i].doublet = (DoubletSiteStr*) calloc (windowSites, sizeof(DoubletSiteStr));
			if (!cell[i].doublet)
				{
				fprintf (stderr, "Could not allocate the doublet window of cell %d\n", i);
				exit (-1);
				}
			}
		}
	}

/********************* StreamReads  ************************/
/* Simulates the reads of very large data sets without storing them. Sites are processed in blocks: the reads and genotype
   likelihoods of all the cells are simulated (in parallel over the sites of the block, each site with its own random
   stream, as in GenerateReadCounts), the sites that look like SNVs are printed (VCF, PILEUP), and the block is discarded.
   The window of sites held in cell[].site starts at siteOff. */
void StreamReads (long int *seed, FILE *fpVCFout, FILE *fpPILEUPout, FILE *fpBulkOut, FILE *fpCATGbody)
	{
	int			i, j, j0, nb, q, nq, nSNVs = 0;
	double		**ngsEij;
	uint64_t	streamBase;
	long int	s1, s2;
	char		**lineBuf = NULL;
	size_t		*lineLen = NULL;
	int			chunk = 1;

	ngsEij = (double**) calloc (4, sizeof(double*));
	for (i=0; i<4; i++)
		ngsEij[i] = (double*) calloc (4, sizeof(double));
	SetupReadModel (ngsEij);

	SetupDoublets (seed, streamBlock);
	SetupCellCoverage (seed);

	RandomUniform (seed);
	s1 = *seed;
	RandomUniform (seed);
	s2 = *seed;
	streamBase = ((uint64_t) s1 << 32) ^ (uint64_t) s2;

	if (fpVCFout)
		PrintVCFHeader (fpVCFout);
#ifdef _OPENMP
	chunk = 4 * omp_get_max_threads();
#endif
	lineBuf = (char**) calloc (chunk, sizeof(char*));
	lineLen = (size_t*) calloc (chunk, sizeof(size_t));

	for (j0=0; j0<numSites; j0+=streamBlock)
		{
		int		*blockSNVs;

		nb = (numSites - j0 < streamBlock) ? numSites - j0 : streamBlock;
		siteOff = j0;
		glUsed = 0;

		#pragma omp parallel
			{
			ReadScratch	rs;
			int			jj;

			NewReadScratch (&rs);
			#pragma omp for schedule(dynamic,1)
			for (jj=j0; jj<j0+nb; jj++)
				SimulateSiteReads (jj, streamBase, ngsEij, &rs);
			FreeReadScratch (&rs);
			}

		blockSNVs = (int*) malloc (nb * sizeof(int));
		nq = 0;
		for (j=j0; j<j0+nb; j++)
			if (allSites[j].isSNV == YES)
				{
				SNVsites[nSNVs++] = j;
				blockSNVs[nq++] = j;
				}

		/* maximum likelihood alleles of the block's SNVs, kept for the genotype and haplotype files */
		if ((doPrintSNVgenotypes == YES || doPrintSNVhaplotypes == YES) && nq > 0)
			{
			long need = (long) nSNVs * (numCells+1);

			if (need > mlStoreCap)
				{
				long newCap = need + need / 2 + (long) (numCells+1) * 1024;
				unsigned char *t = (unsigned char*) realloc (mlStore, (size_t) newCap);

				if (!t)
					{
					fprintf (stderr, "Could not allocate the store of ML alleles\n");
					exit (-1);
					}
				mlStore = t;
				mlStoreCap = newCap;
				}
			#pragma omp parallel for schedule(dynamic,1)
			for (q=0; q<nq; q++)
				StoreMLAlleles (nSNVs - nq + q, blockSNVs[q]);
			}

		if (doPrintFullGenotypes == YES || doPrintFullHaplotypes == YES)
			StoreFullBlock (j0, nb);

		/* CATG lines of the block, built in parallel and written in order to the temporary body file */
		if (fpCATGbody)
			for (i=0; i<nq; i+=chunk)
				{
				int		m = (nq - i < chunk) ? nq - i : chunk;

				#pragma omp parallel for schedule(dynamic,1)
				for (q=0; q<m; q++)
					{
					FILE *ms = open_memstream (&lineBuf[q], &lineLen[q]);
					if (!ms)
						{
						fprintf (stderr, "Could not open a memory stream for the CATG line\n");
						exit (-1);
						}
					CATGSiteLine (ms, blockSNVs[i+q]);
					fclose (ms);
					}
				for (q=0; q<m; q++)
					{
					fwrite (lineBuf[q], 1, lineLen[q], fpCATGbody);
					free (lineBuf[q]);
					}
				}

		/* VCF lines of the block, built in parallel and written in order */
		if (fpVCFout)
			for (i=0; i<nq; i+=chunk)
				{
				int		m = (nq - i < chunk) ? nq - i : chunk;

				#pragma omp parallel for schedule(dynamic,1)
				for (q=0; q<m; q++)
					{
					FILE *ms = open_memstream (&lineBuf[q], &lineLen[q]);
					if (!ms)
						{
						fprintf (stderr, "Could not open a memory stream for the VCF line\n");
						exit (-1);
						}
					VCFSiteLine (ms, blockSNVs[i+q]);
					fclose (ms);
					}
				for (q=0; q<m; q++)
					{
					fwrite (lineBuf[q], 1, lineLen[q], fpVCFout);
					free (lineBuf[q]);
					}
				}
		if (fpPILEUPout)
			for (q=0; q<nq; q++)
				PileupSite (fpPILEUPout, blockSNVs[q]);
		if (fpBulkOut)
			for (j=j0; j<j0+nb; j++)
				BulkSite (fpBulkOut, j);

		free (blockSNVs);
		for (j=j0; j<j0+nb; j++)
			{
			free (siteGL[j]);
			siteGL[j] = NULL;
			}
		}
	if (fpVCFout)
		fprintf (fpVCFout, "\n");

	siteOff = 0;
	numSNVs = nSNVs;
	free (lineBuf);
	free (lineLen);
	free (doubletDonor);
	doubletDonor = NULL;
	for (i=0; i<4; i++)
		free (ngsEij[i]);
	free (ngsEij);
	}


/********************* VCF file (plain or gzip) ************************/
/* Opens the VCF file, compressing it with gzip when requested (-z): the data is piped to "gzip -c" */
static char		bcfTextPath[4*MAX_NAME];

FILE *OpenVCF (const char *path)
	{
	char	cmd[2*MAX_NAME+32];
	FILE	*fp;

	if (doBCF == YES)	/* the VCF is written as text next to the result and converted to BCF when it is closed */
		{
		snprintf (bcfTextPath, sizeof(bcfTextPath), "%s", path);
		return fopen (path, "w");
		}

	if (doZipVCF == NO)
		return fopen (path, "w");
	if (strchr (path, '\'') != NULL)
		{
		fprintf (stderr, "Cannot compress the VCF when the path has a quote: %s\n", path);
		exit (-1);
		}
	snprintf (cmd, sizeof(cmd), "gzip -c > '%s.gz'", path);
	fp = popen (cmd, "w");
	if (fp == NULL)
		fprintf (stderr, "Could not run gzip to write the compressed VCF\n");
	return fp;
	}

void CloseVCF (FILE *fp)
	{
	if (doBCF == YES)
		{
		char	bcfPath[4*MAX_NAME+8];

		fclose (fp);
		snprintf (bcfPath, sizeof(bcfPath), "%s.bcf", bcfTextPath);
		if (ConvertVCFtoBCF (bcfTextPath, bcfPath) != 0)
			{
			fprintf (stderr, "\n!!! Could not write the BCF file %s\n\n", bcfPath);
			exit (-1);
			}
		remove (bcfTextPath);
		return;
		}
	if (doZipVCF == YES)
		pclose (fp);
	else
		fclose (fp);
	}

/* The VCF lists the apparent SNVs whose pooled VAF in reads (alternate reads / reads) is at least minVAFvcf (-M) */
int VCFSitePasses (int j)
	{
	long	depth, alt;

	if (minVAFvcf <= 0.0)
		return YES;
	depth = allSites[j].readCountACGT;
	if (depth <= 0)
		return NO;
	alt = depth - ((allSites[j].referenceAllele == A) ? allSites[j].readCountA : (allSites[j].referenceAllele == C) ? allSites[j].readCountC :
				   (allSites[j].referenceAllele == G) ? allSites[j].readCountG : allSites[j].readCountT);
	return ((double) alt / (double) depth >= minVAFvcf) ? YES : NO;
	}



/********************* ML alleles of all the sites (-3, -4 with reads) ************************/
/* The full genotype and haplotype files list the ML alleles of every site. When the reads are streamed the window of
   sites does not keep them, so StreamReads writes them by blocks into a temporary file (one byte per cell and site, the
   maternal allele in the low nibble, the paternal in the high nibble, 15 = missing) laid out by rows (cells), and the
   files are written row by row from it. */
static FILE			*fullStore = NULL;
static unsigned char	*fullRow = NULL;

void StoreFullBlock (int j0, int nb)
	{
	int		i, j, m, p;
	unsigned char *buf = (unsigned char*) malloc (nb);

	if (!buf)
		{
		fprintf (stderr, "Could not allocate the block of ML alleles\n");
		exit (-1);
		}
	if (fullStore == NULL)
		{
		fullStore = tmpfile ();
		if (!fullStore)
			{
			fprintf (stderr, "Could not open a temporary file for the full genotypes\n");
			exit (-1);
			}
		}
	for (i=0; i<numCells+1; i++)
		{
		for (j=0; j<nb; j++)
			{
			if (cell[i].hasDoublet == NO)
				{
				m = CS(i,j0+j).MLmatAllele;
				p = CS(i,j0+j).MLpatAllele;
				}
			else
				{
				m = CD(i,j0+j).MLmatAlleleDoublet;
				p = CD(i,j0+j).MLpatAlleleDoublet;
				}
			buf[j] = (unsigned char) (((m == MISSING) ? 15 : m) | (((p == MISSING) ? 15 : p) << 4));
			}
		if (fseeko (fullStore, (off_t) i * (off_t) numSites + (off_t) j0, SEEK_SET) != 0 || fwrite (buf, 1, nb, fullStore) != (size_t) nb)
			{
			fprintf (stderr, "Could not write the temporary file of full genotypes\n");
			exit (-1);
			}
		}
	free (buf);
	}

/* when streaming, loads the ML alleles of the cell i (all the sites) */
void LoadFullRow (int i)
	{
	if (streamMode == NO || doNGS == NO)
		return;
	if (fullRow == NULL)
		fullRow = (unsigned char*) malloc (numSites);
	if (!fullRow || fseeko (fullStore, (off_t) i * (off_t) numSites, SEEK_SET) != 0 || fread (fullRow, 1, numSites, fullStore) != (size_t) numSites)
		{
		fprintf (stderr, "Could not read the temporary file of full genotypes\n");
		exit (-1);
		}
	}

void FullMLPair (int i, int j, int *m, int *p)
	{
	if (streamMode == YES)
		{
		unsigned char b = fullRow[j];

		*m = ((b & 15) == 15) ? MISSING : (int) (b & 15);
		*p = ((b >> 4) == 15) ? MISSING : (int) (b >> 4);
		}
	else if (cell[i].hasDoublet == NO)
		{
		*m = CS(i,j).MLmatAllele;
		*p = CS(i,j).MLpatAllele;
		}
	else
		{
		*m = CD(i,j).MLmatAlleleDoublet;
		*p = CD(i,j).MLpatAlleleDoublet;
		}
	}

void CloseFullStore (void)
	{
	if (fullStore)
		fclose (fullStore);
	fullStore = NULL;
	free (fullRow);
	fullRow = NULL;
	}
