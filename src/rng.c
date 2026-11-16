/*
	CellCoal - random number generators and small helpers
	See cellcoal.h for the shared declarations. This file was split from cellcoal.c.
*/
#include "cellcoal.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "definitions.h"

/***************** Index ***************/
/* Returns index for a given node */

int Index (TreeNode *p)
    {
    return (p == NULL) ? -1 : p->index+1;
    }


/***************** Label ***************/
/* Returns label for a given node */

int Label (TreeNode *p)
    {
    return (p->anc == NULL && p->left == NULL && p->right == NULL) ? -1 : p->label+1;
    }


/************************************* ChooseTrinucMutation **********************************************/
/* Chooses uniformy a random trinucleotide mutation (out of 96) for a given signature 3D matrix with 96 types */

int ChooseTrinucMutation (double ***prob, long int *seed)
	{
	int		i, j = 0, k = 0;
	double	ran, cumProb;
	
	cumProb = 0;
	ran = RandomUniform(seed);
	for (i=0; i<4; i++)
		for (j=0; j<6; j++)
			for (k=0; k<4; k++)
				{
				cumProb += prob[i][j][k];
				if (ran <= cumProb)
					break;
				}
	return trimut(i,j,k);
	}



/************************************* ChooseUniformState **********************************************/
/* Chooses uniformy a random state according to a vector of state probabilities */

int ChooseUniformState (double *prob, long int *seed)
	{
	int			chosenState;
	double		ran, cumProb;
	
	chosenState = 0;
	cumProb = prob[chosenState];
	ran = RandomUniform(seed);

	while (ran > cumProb)
		cumProb += prob[++chosenState];
	
	return chosenState;
	}


/***************************** RandomUniform **********************************/
/*  It returns a random uniform variate in range 0..1. It is described in
 Park, S. K. and K. W. Miller.  1988.  Random number generators: good
 ones are hard to find.  Communications of the ACM, 31(10):1192-1201.
 */

/* One step of the minimal standard generator s <- 16807*s mod (2^31-1), computed with a multiplication and a
   Mersenne reduction instead of Schrage's divisions. For a valid state (1 <= s <= 2^31-2) it gives exactly the
   same value as the original code; any other state goes through the original code. */
#define LCG_M		2147483647UL
#define LCG_A1		16807UL
#define LCG_A2		282475249UL
#define LCG_A3		1622650073UL
#define LCG_A4		984943658UL

static inline int64_t LCG_MulMod (int64_t s, uint64_t a)
	{
	uint64_t	x = a * (uint64_t) s;
	uint64_t	r = (x & LCG_M) + (x >> 31);
	return (int64_t) (r >= LCG_M ? r - LCG_M : r);
	}

double RandomUniform (long int *seed)
	{
	int64_t	lo, hi, test;

	if (streamActive)
		return (double)((SplitMix64 (&streamState) >> 11) + 0.5) * (1.0 / 9007199254740992.0);	/* in (0,1) */
	if (*seed >= 1 && *seed < 2147483647L)
		*seed = LCG_MulMod (*seed, LCG_A1);
	else
		{
		hi = (*seed) / 127773;
		lo = (*seed) % 127773;
		test = 16807 * lo - 2836 * hi;
		if (test > 0)
			*seed = test;
		else
			*seed = test + 2147483647;
		}
	return (double)(*seed) / (double)2147483647;
	}



/********************** RandomUniformTo ****************************/
/* it returns random uniform in range 0...max-1          */

int RandomUniformTo (int max, long int *seed)
	{
    double	rd;
    rd = RandomUniform (seed);
    return (floor(rd*max));
	}




/************************ fast mode (--fast) ******************************/
/* Used only for the reads (per-site streams) when --fast is given. The distributions are the same as in the exact
   samplers, but the random numbers are consumed differently, so the data sets differ for a given seed. */

/* log10 with a absolute error below 1e-9: x = m 2^e with m in [sqrt(1/2), sqrt(2)), log m = 2 atanh ((m-1)/(m+1)) */
double Log10Approx (double x)
	{
	union { double d; uint64_t u; } v;
	double	m, f, f2, r;
	int		e;

	if (!(x > 0.0) || !(x < INFINITY))
		return log10 (x);
	v.d = x;
	if ((v.u >> 52) == 0)	/* subnormal */
		return log10 (x);
	e = (int)(v.u >> 52) - 1023;
	v.u = (v.u & 0x000FFFFFFFFFFFFFULL) | 0x3FF0000000000000ULL;
	m = v.d;				/* [1,2) */
	if (m > 1.4142135623730951)
		{
		m *= 0.5;
		e++;
		}
	f = (m - 1.0) / (m + 1.0);
	f2 = f * f;
	r = 2.0 * f * (1.0 + f2 * (1.0/3 + f2 * (1.0/5 + f2 * (1.0/7 + f2 * (1.0/9)))));
	return (r + e * 0.6931471805599453) * 0.4342944819032518;
	}

/* 10^x with a relative error below 1e-8 (0 below -400, library for huge arguments) */
double Pow10Approx (double x)
	{
	union { double d; uint64_t u; } v;
	double	y, n, r, p;
	int		ni;

	if (x == 0.0)
		return 1.0;
	if (x < -300.0 || x > 300.0 || x != x)
		return x < -400.0 ? 0.0 : pow (10.0, x);
	y = x * 3.321928094887362;		/* log2 (10) */
	ni = (int)(y + (y >= 0 ? 0.5 : -0.5));
	n = ni;
	r = (y - n) * 0.6931471805599453;	/* |r| <= 0.3466 */
	p = 1.0 + r * (1.0 + r * (1.0/2 + r * (1.0/6 + r * (1.0/24 + r * (1.0/120 + r * (1.0/720 + r * (1.0/5040 + r * (1.0/40320))))))));
	v.u = (uint64_t)(ni + 1023) << 52;	/* 2^ni, normal range since |x| <= 300 */
	return p * v.d;
	}


/* Batch versions of Log10Approx / Pow10Approx for up to 16 values (--fast). The loops have no calls and no data-dependent
   branches in the main pass, so the compiler can vectorize them; the rare special values (x <= 0, not finite or
   subnormal for the logarithm; |x| > 300 or NaN for the power) are repaired afterwards with the scalar functions. */
void Log10Batch (const double *x, double *y, int n)
	{
	int		i, special = 0;

	for (i=0; i<n; i++)
		{
		union { double d; uint64_t u; } v;
		double	m, f, f2, r;
		int		e;
		uint64_t	big;

		v.d = x[i];
		e = (int)((v.u >> 52) & 0x7FF) - 1023;
		special |= !(x[i] > 1e-300 && x[i] < 1e300);
		v.u = (v.u & 0x000FFFFFFFFFFFFFULL) | 0x3FF0000000000000ULL;
		m = v.d;
		big = m > 1.4142135623730951;
		m = big ? m * 0.5 : m;
		e += (int) big;
		f = (m - 1.0) / (m + 1.0);
		f2 = f * f;
		r = 2.0 * f * (1.0 + f2 * (1.0/3 + f2 * (1.0/5 + f2 * (1.0/7 + f2 * (1.0/9)))));
		y[i] = (r + e * 0.6931471805599453) * 0.4342944819032518;
		}
	if (special)
		for (i=0; i<n; i++)
			if (!(x[i] > 1e-300 && x[i] < 1e300))
				y[i] = Log10Approx (x[i]);
	}

void Pow10Batch (const double *x, double *y, int n)
	{
	int		i, special = 0;

	for (i=0; i<n; i++)
		{
		union { double d; uint64_t u; } v;
		double	t, r, p;
		int		ni;

		special |= !(x[i] >= -300.0 && x[i] <= 300.0);
		t = x[i] * 3.321928094887362;
		t = t < -1000.0 ? -1000.0 : (t > 1000.0 ? 1000.0 : t);
		ni = (int)(t + (t >= 0 ? 0.5 : -0.5));
		r = (t - ni) * 0.6931471805599453;
		p = 1.0 + r * (1.0 + r * (1.0/2 + r * (1.0/6 + r * (1.0/24 + r * (1.0/120 + r * (1.0/720 + r * (1.0/5040 + r * (1.0/40320))))))));
		v.u = (uint64_t)(ni + 1023) << 52;
		y[i] = p * v.d;
		}
	if (special)
		for (i=0; i<n; i++)
			if (!(x[i] >= -300.0 && x[i] <= 300.0))
				y[i] = Pow10Approx (x[i]);
	}

/* Poisson by inversion with a single uniform (lambda < 30) or the PTRS transformed rejection of Hormann (1993) */
static int PoissonFast (double lambda, long int *seed)
	{
	if (lambda < 30.0)
		{
		double	p = exp (-lambda), cdf, u;
		int		k = 0;

		u = RandomUniform (seed);
		cdf = p;
		while (u > cdf && k < 1000)
			{
			k++;
			p *= lambda / k;
			cdf += p;
			}
		return k;
		}
	else
		{
		double	slam = sqrt (lambda), loglam = log (lambda);
		double	b = 0.931 + 2.53 * slam, a = -0.059 + 0.02483 * b;
		double	invalpha = 1.1239 + 1.1328 / (b - 3.4), vr = 0.9277 - 3.6224 / (b - 2.0);
		double	U, V, us;
		int64_t	k;

		for (;;)
			{
			U = RandomUniform (seed) - 0.5;
			V = RandomUniform (seed);
			us = 0.5 - fabs (U);
			k = (int64_t) floor ((2.0 * a / us + b) * U + lambda + 0.43);
			if (us >= 0.07 && V <= vr)
				return (int) k;
			if (k < 0 || (us < 0.013 && V > us))
				continue;
			if (log (V) + log (invalpha) - log (a / (us * us) + b) <= -lambda + k * loglam - lgamma (k + 1.0))
				return (int) k;
			}
		}
	}

/* binomial by geometric jumps between successes (exact; one log per success instead of one uniform per trial) */
static int BinomialFast (double prob, int numTrials, long int *seed)
	{
	double	lq;
	int		pos = 0, count = 0, flip = 0;

	if (prob > 0.5)
		{
		prob = 1.0 - prob;
		flip = 1;
		}
	if (prob < 0.25 && numTrials > 8)
		{
		lq = log1p (-prob);
		for (;;)
			{
			double jump = 1.0 + floor (log (RandomUniform (seed)) / lq);

			if (jump > (double)(numTrials - pos))
				break;
			pos += (int) jump;
			count++;
			}
		}
	else
		{
		uint64_t	thr53 = (uint64_t)(prob * 9007199254740992.0);
		int			i;

		for (i=0; i<numTrials; i++)
			count += (SplitMix64 (&streamState) >> 11) < thr53;
		}
	return flip ? numTrials - count : count;
	}

/********************* RandomPoisson ********************/
/* Generates a random number from a Poisson distibution with
 mean lambda.
 */

int RandomPoisson (double lambda, long int *seed)
	{
	int		poissonNumber;
	double	sum;

	/* Fast path. The original algorithm adds Exponential(lambda) draws until the sum exceeds 1,
	   which costs one log() per draw. The decision "sum > 1" is the same as "prod(U) < exp(-lambda)",
	   so we multiply the uniforms and only fall back to the exact sum of logs (using the stored
	   uniforms) when the product is within 1e-9 (relative) of the threshold. The number of uniforms
	   consumed and the returned value are identical to the original algorithm. */
	if (fastMode && streamActive && lambda > 0.0)
		return PoissonFast (lambda, seed);
	if (lambda > 0.0 && lambda < 500.0)
		{
		enum { CAP = 4096 };
		double	U[CAP], u, P, thr, lo, hi;
		int		k, i;

		thr = exp (-lambda);
		lo = thr * (1.0 - 1e-9);
		hi = thr * (1.0 + 1e-9);
		P  = 1.0;
		k  = 0;
		for (;;)
			{
			do
				u = RandomUniform (seed);
			while (u == 0);
			if (k < CAP)
				U[k] = u;
			k++;
			P *= u;
			if (P < lo)
				return (k - 1);
			if (P > hi && k < CAP)
				continue;
			/* ambiguous (or too many draws to store): exact computation as in the original */
			if (k <= CAP)
				{
				sum = 0;
				for (i=0; i<k; i++)
					sum += -log (U[i]) / lambda;
				if (sum > 1.0)
					return (k - 1);
				}
			else if (P < hi)
				return (k - 1);
			}
		}

	/* a mean of billions is absurd here (the callers reject such numbers of mutations); the loop below would need
	   that many draws and its counter would overflow */
	if (lambda > 2.0e9)
		return 2000000000;
	sum = 0;
	poissonNumber = -1;
	while (sum <= 1.0)
		{
		sum += RandomExponential (lambda, seed);
		poissonNumber++;
		}
	return poissonNumber;
	}


/********************* RandomExponential ********************/
/* Generates a random number from a Exponential distibution with
 mean lambda.
 */

double RandomExponential (double lambda, long int *seed)
	{
		
    double 	exponentialNumber, U;
		
    do
		U = RandomUniform (seed);
	while (U == 0);
		
    exponentialNumber = -log (U) / lambda;
		
    return exponentialNumber;
	}


/********************* RandomBinomial ********************/
/* Generates a random number from a Binomial distibution with sucess probabilty p
	and n trials using the direct method (sum of Bernoulli variables).
 */

int RandomBinomial (double prob, int numTrials, long int *seed)
	{
	/* Same stream and result as summing Bernoulli(RandomUniform(seed) < prob), but the comparison
	   (double)s/2147483647 < prob is turned into the exactly equivalent integer comparison s < t,
	   where t is the smallest integer with (double)t/2147483647 >= prob (the division is monotonic).
	   This avoids one floating-point division per trial. */
	const int64_t	m = 2147483647LL;
	const double M = 2147483647.0;
	int64_t		s, t, hi, lo, test;
	int			i, sum;

	if (streamActive)	/* the fast path below works on the state of the minimal standard generator */
		{
		/* Bernoulli trial: the top 53 bits of the draw are below prob * 2^53 (probability prob to within 2^-53) */
		uint64_t	thr53;

		if (!(prob > 0.0))
			return 0;
		if (prob >= 1.0)
			return numTrials > 0 ? numTrials : 0;
		if (fastMode)
			return BinomialFast (prob, numTrials, seed);
		thr53 = (uint64_t)(prob * 9007199254740992.0);
		for (sum=0, i=0; i<numTrials; i++)
			sum += (SplitMix64 (&streamState) >> 11) < thr53;
		return sum;
		}

	if (!(prob > 0.0))
		t = 0;
	else if (prob >= 1.0)
		t = m;
	else
		{
		t = (int64_t) (prob * M);
		if (t < 0) t = 0;
		if (t > m) t = m;
		while (t > 0 && (double)(t - 1) / M >= prob)
			t--;
		while (t < m && (double)t / M < prob)
			t++;
		}

	s = *seed;
	sum = 0;
	i = 0;
	if (s >= 1 && s < m)
		{
		/* valid state: the next four states are independent products of the current one (a, a^2, a^3, a^4 mod m),
		   so the four steps can be computed at the same time; same states, same order as one step at a time */
		for (; i+4 <= numTrials; i+=4)
			{
			int64_t s1 = LCG_MulMod (s, LCG_A1);
			int64_t s2 = LCG_MulMod (s, LCG_A2);
			int64_t s3 = LCG_MulMod (s, LCG_A3);
			int64_t s4 = LCG_MulMod (s, LCG_A4);
			sum += (s1 < t) + (s2 < t) + (s3 < t) + (s4 < t);
			s = s4;
			}
		for (; i<numTrials; i++)
			{
			s = LCG_MulMod (s, LCG_A1);
			if (s < t)
				sum++;
			}
		}
	for (; i<numTrials; i++)
		{
		hi = s / 127773;
		lo = s % 127773;
		test = 16807 * lo - 2836 * hi;
		s = (test > 0) ? test : test + m;
		if (s < t)
			sum++;
		}
	*seed = s;
	return sum;
	}

/********************* RandomNegativeBinomial ********************/
/*
 *    Random variates from the negative binomial distribution.
 *
 *  NOTES
 *
 *    x = the number of failures before the n-th success
 *
 *  REFERENCE
 *
 *    Devroye, L. (1986).
 *    Non-Uniform Random Variate Generation.
 *    New York:Springer-Verlag.  Pages 488 and 543.
 *
 *  METHOD
 *
 *    Generate lambda as gamma with shape parameter "dispersion" (aka size) and scale
 *    parameter "mean/dispersion".  Return a Poisson deviate with mean lambda.

**** NOTE: Extracted from rnbinom.c R code:
	rpois(rgamma(size, (1 - prob) / prob));
	rpois(rgamma(size, mu / size));

The negative binomial distribution with dispersion = n and prob = p has density

p(x) = Gamma(x+n)/(Gamma(n) x!) p^n (1-p)^x

for x = 0, 1, 2, ..., n > 0 and 0 < p <= 1.

A negative binomial distribution can arise as a mixture of Poisson distributions with mean distributed as a Γ (pgamma) distribution with scale parameter (1 - prob)/prob and shape parameter dispersion. In this model prob = scale/(1+scale), and the mean is dispersion * (1 - prob)/prob. The variance in this parametrization is n (1-p)/p^2.

The alternative parameterization, often used in ecology, and the one used here, is by the mean mu, and the dispersion parameter, where prob = dispersion/(dispersion+mu). The variance is mu + mu^2/dispersion in this parametrization.
*/

int RandomNegativeBinomial (double mean, double dispersion, long int *seed)
	{
	int		poissonRand;
	double	gammaRand;
	
	/* the RandomGamma function here has mean 1, so we need to scale it ourselves */
	gammaRand = mean/dispersion * RandomGamma (dispersion, seed);
	poissonRand = RandomPoisson (gammaRand, seed);
		
	return poissonRand;
	}


/**************************** RandomBeta *************************/
/*	Generates a beta random number given shape1 and shape2

An algorithm for generating beta variates B(α,β) is to generate X/(X + Y), where X is a gamma variate with parameters (α, 1) and Y is an independent gamma variate with parameters (β, 1).[53]
From Numerical Recipes 3rd Edition: The Art of Scientific Computing. 2007. Press et al. ISBN-13: 978-0521880688
*/
/*
double RandomBeta (double alpha, double beta, long int *seed)
	{
	double gamma1, gamma2, randBeta;
	
	gamma1 = RandomGamma(alpha, seed);
	gamma2 = RandomGamma(beta, seed);
	randBeta = gamma1 / (gamma1 + gamma2);

	return randBeta;
	}
*/

/**************************** RandomBetaMeanVar *************************/
/*	Generates a beta random number given mean and variance

An algorithm for generating beta variates B(α,β) is to generate X/(X + Y), where X is a gamma variate with parameters (α, 1) and Y is an independent gamma variate with parameters (β, 1).[53]
From Numerical Recipes 3rd Edition: The Art of Scientific Computing. 2007. Press et al. ISBN-13: 978-0521880688

note that the variance has to be > 0
*/

double RandomBetaMeanVar (double mean, double var, long int *seed)
	{
	double sample_size, shape1, shape2, gamma1, gamma2, randBeta;
	
	/* assuming variance < mean (1-mean) */
	sample_size = (mean * (1.0 - mean) / var) - 1.0;
	
	shape1 = mean * sample_size;
	shape2 = (1.0 - mean) * sample_size;
	
	gamma1 = RandomGamma(shape1, seed);
	gamma2 = RandomGamma(shape2, seed);
	randBeta = gamma1 / (gamma1 + gamma2);
	
	return randBeta;
	}



/* Gamma functions are taken from Ziheng Yang's PAML package.
 http://abacus.gene.ucl.ac.uk/*/


/**************************** RandomGamma *************************/
/*	Generates a gamma number using routines in Ziheng's
	Yang tools.h in PAML
 
	Random standard gamma (Mean=Var=s,  with shape par=s, scale par=1)
	r^(s-1)*exp(-r)
 
	J. Dagpunar (1988) Principles of random variate generation,
	Clarendon Press, Oxford
 
	Calling rndgamma1() if s<1 or rndgamma2() if s>1 or exponential if s=1
 */

double	RandomGamma (double shape, long int *seed)
{
    double gammaNumber = 0;
	
    if (shape <= 0)
        fprintf (stderr, "!!! ERROR: problems with gamma variable generation, shape < 0");
    else if (shape < 1)
        gammaNumber = RandomGamma1 (shape, seed);
    else if (shape > 1)
        gammaNumber = RandomGamma2 (shape, seed);
    else
        gammaNumber = -log (RandomUniform(seed));
	
    return (gammaNumber);
	
	
}

/*************** RandomGamma1 ***************/
double RandomGamma1 (double s, long int *seed)
{
    /* Random standard gamma for s<1
     switching method
     */
    double			r, x=0.0, small=1e-37, w;
    static __thread double   a, p, uf, ss=10.0, d;
	
    if (s!=ss)
    {
        a  = 1.0-s;
        p  = a/(a+s*exp(-a));
        uf = p*pow(small/a,s);
        d  = a*log(a);
        ss = s;
    }
    for (;;)
    {
        r = RandomUniform(seed);
        if (r > p)
			{
            x = a-log((1.0-r)/(1.0-p));
            w=a*log(x)-d;  /* this was with comma in line above before 270917*/
			}
        else if (r>uf)
            {
            x = a*pow(r/p,1/s);
			w=x; /* this was with comma in line above before 270917*/
			}
        else
            return (0.0);
        r = RandomUniform(seed);
        if (1.0-r <= w && r > 0.0)
            if (r*(w+1.0) >= 1.0 || -log(r) <= w)
                continue;
        break;
    }
    return (x);
}


/*************** RandomGamma2 ***************/
double RandomGamma2 (double s, long int *seed)
{
    /* Random standard gamma for s>1
     Best's (1978) t distribution method
     */
    double			r ,d, f, g, x;
    static __thread double	b, h, ss=0;

    /* for an enormous shape (above 1e7) the gamma is practically normal (mean s, variance s), and Best's rejection
       loop would not terminate because b+f loses all precision */
    if (s > 1e7)
        {
        do {
            r = RandomUniform(seed);
            d = RandomUniform(seed);
            x = s + sqrt(s) * sqrt(-2.0*log(r > 0.0 ? r : 1e-300)) * cos(6.283185307179586*d);
            } while (x <= 0.0);
        return x;
        }
    if (s!=ss)
    {
        b  = s-1.0;
        h  = sqrt(3.0*s-0.75);
        ss = s;
    }
    for (;;)
    {
        r = RandomUniform(seed);
        g = r-r*r;
        f = (r-0.5)*h/sqrt(g);
        x = b+f;
        if (x <= 0.0)
            continue;
        r = RandomUniform(seed);
        d = 64*r*r*g*g*g;
        if (d*x < x-2.0*f*f || log(d) < 2*(b*log(x/b)-f))
            break;
    }
    return (x);
}



/************************ ReadUntil **************************/
/* Utility to skip comments between brackets. Borrowed from Andy
 Rambaut and Nick Grassly */

void ReadUntil(FILE *fv, char stopChar, char *what)
{
    char ch;
	
    ch=fgetc(fv);
    while (!feof(fv) && ch!=stopChar)
        ch=fgetc(fv);
	
    if (feof(fv) || ch!=stopChar) {
        fprintf (stderr, "%s missing", what);
        exit(0);
    }
}



/************************ RemoveDir **************************/
/* Removes a directory and its contents (portable: POSIX and Windows); returns 0 when everything was removed or the path did not exist */
int RemoveDir (char *path)
	{
	struct stat	st;
	DIR			*d;
	struct dirent	*e;
	int			rc = 0;

	if (LSTAT (path, &st) != 0)
		return 0;
	if (!S_ISDIR (st.st_mode))
		return unlink (path);
	if ((d = opendir (path)) == NULL)
		return -1;
	while ((e = readdir (d)) != NULL)
		{
		char	*sub;
		size_t	len;

		if (strcmp (e->d_name, ".") == 0 || strcmp (e->d_name, "..") == 0)
			continue;
		len = strlen (path) + strlen (e->d_name) + 2;
		sub = (char*) malloc (len);
		if (!sub)
			{
			rc = -1;
			break;
			}
		snprintf (sub, len, "%s/%s", path, e->d_name);
		if (RemoveDir (sub) != 0)
			rc = -1;
		free (sub);
		}
	closedir (d);
	if (rmdir (path) != 0)
		rc = -1;
	return rc;
	}


/************************ CheckMatrixSymmetry **************************/
/* Checks whether a given matrix is symmetric */

int CheckMatrixSymmetry(double matrix[4][4])
	{
	int i,j;
	
	for(i=0; i<4; i++)
		for(j=0; j<4; j++)
			if(matrix[i][j] != matrix[j][i])
				return NO;
	return YES;
	}


