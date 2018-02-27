//
//  eigein.h
//  CoalTumor
//
//  Created by David Posada on 21/11/2017.
//  Everything below is shamelessly taken from Yang's Paml package */

//

#ifndef eigein_h
#define eigein_h

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

FILE			*fpout, *fperr;

/* Everything below is shamelessly taken from Yang's Paml package */
int abyx (double a, double x[], int n);
int xtoy (double x[], double y[], int n);
int matinv( double x[], int n, int m, double space[]);
int eigen(int job, double AA[], int n, double rr[], double ri[], double vr[], double vi[], double w[]);
void balance(double mat[], int n, int *low, int *hi, double scale[]);
void unbalance(int n, double vr[], double vi[], int low, int hi, double scale[]);
int realeig(int job, double mat[], int n,int low, int hi, double valr[], double vali[], double vr[], double vi[]);
void elemhess(int job, double mat[], int n, int low, int hi, double vr[], double vi[], int work[]);


typedef struct { double re, im; } complex;
#define csize(a) (fabs(a.re)+fabs(a.im))

complex compl (double re,double im);
complex conjj (complex a);
complex cplus (complex a, complex b);
complex cminus (complex a, complex b);
complex cby (complex a, complex b);
complex cdiv (complex a,complex b);
complex cexpp (complex a);
complex cfactor (complex x, double a);
int cxtoy (complex x[], complex y[], int n);
int cmatby (complex a[], complex b[], complex c[], int n,int m,int k);
int cmatout (FILE * fout, complex x[], int n, int m);
int cmatinv( complex x[], int n, int m, double space[]);

#endif /* eigein_h */
