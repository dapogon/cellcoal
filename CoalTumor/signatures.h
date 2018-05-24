//
//  signatures.h
//  CoalTumor
//
//  Created by David Posada on 06/04/2018.
//  Copyright © 2018 David Posada. All rights reserved.
//

#ifndef signatures_h
#define signatures_h

void 		PrepareGeneticSignatures(void);
void 		LoadSignatures(void);
void		LoadTrinucleotideFrequencies(void);

extern void PrintUsage(void);
extern char	WhichNuc (int nucleotide);
extern int	userSignature;

double 		***selectedSignature; 	// [5prime][change][3prime]
double		****geneticSignature; 	// [signature][5prime][change][3prime]
double 		*signatureProbs;
double 		*triNucFreq;


#endif /* signatures_h */
