#include "mpp.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#ifdef __unix__
#else
    #define M_PI       3.14159265358979323846
#endif

/* Fix 5: proper prototypes at file scope replacing K&R-style in-body declarations */
void   GaborGetIndexForGauss(int j1, int j2, int nDeltaP, int nL, int maxlh,
                              int l1, int l2, int m, int *k, int *n);
void   GaborGetIndexForCExp(int nDeltaJ, int tw,
                             int *alpha, int *beta, int *gamma);
double GaborGetGauss(double *pfG, int *pnIep, int k, int n);
void   GaborGetCExp(int alpha, int beta, int gamma, int N, int n,
                    double *pfCE1, double *pfCE2, double *pfCE3,
                    double *pfCos, double *pfSin);

/*
 * GaborGetInnerProd — calculate the inner product of two waveforms.
 *
 * Inputs:
 *	index1, index2	indices for the two waveforms (INDEX)
 *	nL		maximum transformation level; gabsignal size = 2^nL (int)
 *	nST		subsampling level for time (int)
 *	nSF		subsampling level for frequency (int)
 *	nDeltaP		delta p (int)
 *	nDeltaW		delta w (int)
 *	nBound		index bound for the sum (int)
 *	C		normalisation coefficient (double)
 *	pnAep		array for the length of k for each n (int *)
 *	pnIep		array for the index of n in pfG (int *)
 *	pfCE1		cos/sin table on integers (double *)
 *	pfCE2		cos/sin table on k*2^n/(1+2^(2n)) (double *)
 *	pfCE3		cos/sin table on k/(1+2^(2n)) (double *)
 *	pfG		Gaussian lookup table (double *)
 *
 * Outputs:
 *	pfReal		real part of the inner product (double *)
 *	pfImaginary	imaginary part of the inner product (double *)
 */
void GaborGetInnerProd(INDEX index1, INDEX index2, int nL, int nST, int nSF,
                       int nDeltaP, int nDeltaW,
                       int nBound, double C, int *pnAep,
                       int *pnIep, double *pfCE1, double *pfCE2, double *pfCE3,
                       double *pfG, double *pfReal, double *pfImaginary)
{
    int m, q;
    int j1, j2, p1, p2;
    int kt, nt, kw, nw, kw1, wk, wk1, tk, twk;
    int l1, l2, h1, h2;
    int abskt, abskw;
    int alpha, beta, gamma;
    int N;
    int nStepM, nStepQ;
    long nIndex;   /* Fix 6: long to avoid overflow on large signals */
    int maxlh;
    double tmpt, tmpw, tmptw, tmpi, tmpr;
    double fCos, fSin;

    tmpr = tmpi = 0.0;
    j1 = (int)index1->octave;
    j2 = (int)index2->octave;
    p1 = (int)index1->position;
    p2 = (int)index2->position;
    N = 1 << nL;
    maxlh = MAX(nST, nSF);
    l1 = j1 > nST ? nST : j1;
    l2 = j2 > nST ? nST : j2;
    h1 = nL-j1 > nSF ? nSF : nL-j1;
    h2 = nL-j2 > nSF ? nSF : nL-j2;

    if (j2 > j1)
        nStepM = 1 << (nL-j1+maxlh);
    else
        nStepM = 1 << (nL-j2+maxlh);
    if (j2 < j1)
        nStepQ = 1 << (j1+maxlh);
    else
        nStepQ = 1 << (j2+maxlh);

    GaborGetIndexForGauss(j1, j2, nDeltaP, nL, maxlh, l1, l2, -nBound, &kt, &nt);
    GaborGetIndexForGauss(nL-j1, nL-j2, nDeltaW, nL, maxlh, h1, h2, -nBound, &kw1, &nw);
    wk1 = nDeltaW - nBound*N;
    tk  = nDeltaP - nBound*N;

    for (m = -nBound; m <= nBound; m++) {
        abskt = abs(kt);
        if (abskt >= pnAep[nt]) {
            kt += nStepM;
            tk += N;
            continue;
        }
        tmpt = GaborGetGauss(pfG, pnIep, abskt, nt);
        kw = kw1;
        wk = wk1;
        for (q = -nBound; q <= nBound; q++) {
            abskw = abs(kw);
            if (abskw >= pnAep[nw]) {
                kw += nStepQ;
                wk += N;
                continue;
            }
            tmpw  = GaborGetGauss(pfG, pnIep, abskw, nw);
            tmptw = tmpw * tmpt;
            twk   = wk * tk;
            GaborGetIndexForCExp(nt, abs(twk), &alpha, &beta, &gamma);
            GaborGetCExp(alpha, beta, gamma, N, nt,
                         pfCE1, pfCE2, pfCE3, &fCos, &fSin);
            tmpr += fCos * tmptw;
            if (twk < 0)
                tmpi -= fSin * tmptw;
            else
                tmpi += fSin * tmptw;
            kw += nStepQ;
            wk += N;
        }
        kt += nStepM;
        tk += N;
    }

    if (j2 < j1) {
        /* Fix 6: cast to long before multiply to avoid int overflow */
        nIndex = ((long)p2 * (long)abs(nDeltaW)) % N;
        fCos = pfCE1[nIndex];
        fSin = (nDeltaW > 0) ? pfCE1[nIndex+N] : -pfCE1[nIndex+N];
        *pfReal      =  C * (tmpr*fCos + tmpi*fSin);
        *pfImaginary =  C * (-tmpr*fSin + tmpi*fCos);
    } else {
        /* Fix 6: cast to long before multiply to avoid int overflow */
        nIndex = ((long)p1 * (long)abs(nDeltaW)) % N;
        fCos = pfCE1[nIndex];
        fSin = (nDeltaW > 0) ? pfCE1[nIndex+N] : -pfCE1[nIndex+N];
        *pfReal      =  C * (tmpr*fCos - tmpi*fSin);
        *pfImaginary = -C * (tmpr*fSin + tmpi*fCos);
    }
}

/*
 * GaborGetIndexForGauss — compute the (k, n) index into the Gaussian table.
 *
 * Inputs:
 *	j1, j2		octaves of the two waveforms (int)
 *	nDeltaP		position difference (int)
 *	nL		maximum transformation level (int)
 *	maxlh		max(l, h) subsampling (int)
 *	l1, l2		time subsampling factors (int)
 *	m		multiple of N (int)
 *
 * Outputs:
 *	k, n		index into the Gaussian lookup table (int *)
 */
void GaborGetIndexForGauss(int j1, int j2, int nDeltaP, int nL, int maxlh,
                            int l1, int l2, int m, int *k, int *n)
{
    int tmp;

    *n  = abs(j2-j1);
    tmp = abs(nDeltaP);
    if (j2 > j1) {
        tmp >>= j1-l1;
        *k = (nDeltaP > 0) ?  tmp + m*(1<<(nL-j1+l1))
                            : -tmp + m*(1<<(nL-j1+l1));
        *k <<= maxlh-l1;
    } else {
        tmp >>= j2-l2;
        *k = (nDeltaP > 0) ?  tmp + m*(1<<(nL-j2+l2))
                            : -tmp + m*(1<<(nL-j2+l2));
        *k <<= maxlh-l2;
    }
}

/*
 * GaborGetGauss — look up a Gaussian value from the precomputed table.
 *
 * Caller owns pfG — allocated by GaborGetGaussianArray() and must be
 * freed by the caller when no longer needed.
 */
double GaborGetGauss(double *pfG, int *pnIep, int k, int n)
{
    return pfG[pnIep[n] + k];
}

/*
 * GaborGetNAep — compute the k-bound and index arrays for the Gaussian table.
 *
 * Inputs:
 *	nL	maximum transformation level; gabsignal size = 2^nL (int)
 *	maxlh	max(l, h) subsampling (int)
 *	epslon	precision (double)
 *
 * Outputs:
 *	pnAep	k-bound for each n; must be pre-allocated with size >= nL-1 (int *)
 *	pnIep	index of n in pfG; must be pre-allocated with size >= nL (int *)
 *
 * Fix 10: hoisted loop-invariant -log(epslon)/M_PI out of the loop.
 */
void GaborGetNAep(int nL, int maxlh, double epslon, int *pnAep, int *pnIep)
{
    int i, sum = 0;
    double logterm;   /* Fix 10: loop-invariant constant */

    /* Fix 3: perror() → fprintf+return */
    if (pnAep == (int *)NULL || pnIep == (int *)NULL) {
        fprintf(stderr, "GaborGetNAep(): null inputs!\n");
        return;
    }

    /* Fix 10: compute once outside the loop */
    logterm = -log(epslon) / M_PI;

    pnIep[0] = 0;
    for (i = 0; i < nL-1; i++) {
        pnAep[i] = (int)(sqrt((double)(1+(1<<(2*i))) * logterm) + 1.0)
                   * (1<<maxlh);
        sum += pnAep[i];
        pnIep[i+1] = sum;
    }
}

/*
 * GaborGetIndexForCExp — decompose tw into (alpha, beta, gamma) indices
 * for the complex exponential lookup tables.
 */
void GaborGetIndexForCExp(int nDeltaJ, int tw, int *alpha, int *beta, int *gamma)
{
    int tmp1, tmp2, k;

    k     = tw;
    tmp2  = 1 + (1<<(2*nDeltaJ));
    *alpha = k / tmp2;
    tmp1  = k % tmp2;
    tmp2  = 1<<nDeltaJ;
    *beta  = tmp1 / tmp2;
    *gamma = tmp1 % tmp2;
}

/*
 * GaborGetCExp — evaluate the complex exponential at (alpha, beta, gamma)
 * using the three precomputed lattice tables.
 *
 * Outputs:
 *	pfCos, pfSin	cos and sin of the combined phase (double *)
 */
void GaborGetCExp(int alpha, int beta, int gamma, int N, int n,
                  double *pfCE1, double *pfCE2, double *pfCE3,
                  double *pfCos, double *pfSin)
{
    int index1, index2, size;
    double a1, a2, a3, b1, b2, b3;
    double tmp1, tmp2;

    index1 = (1<<(n+1)) - 2;
    index2 = index1 + (n<<1);
    size   = 1<<n;
    alpha %= N;

    a1 = pfCE1[alpha];
    b1 = pfCE1[alpha+N];
    a2 = pfCE2[index2+beta];
    b2 = pfCE2[index2+size+1+beta];
    a3 = pfCE3[index1+gamma];
    b3 = pfCE3[index1+size+gamma];

    /* Two sequential complex multiplications with shared intermediates —
     * already optimal; no further simplification needed. */
    tmp1   = a1*a2 - b1*b2;
    tmp2   = a1*b2 + a2*b1;
    *pfCos = a3*tmp1 - b3*tmp2;
    *pfSin = a3*tmp2 + b3*tmp1;
}

/* Fix 4: GaborGetSin removed — it was dead code (never called) and its
 * single-term approximation was mathematically incorrect. GaborGetCExp
 * provides both cos and sin correctly. */

/*
 * GaborGetGaussianArray — allocate and fill the Gaussian lookup table.
 *
 * Inputs:
 *	nL	maximum transformation level (int)
 *	pnAep	k-bound array from GaborGetNAep() (int *)
 *	length	total number of elements to allocate (int)
 *	maxlh	max(l, h) subsampling (int)
 *
 * Returns a newly allocated array. Caller is responsible for freeing it.
 *
 * Fix 2: added NULL check after malloc (was missing, unlike all other
 *        allocating functions in this file).
 */
double *GaborGetGaussianArray(int nL, int *pnAep, int length, int maxlh)
{
    double *pfG;
    int k, n;
    double d, t, *v;

    pfG = (double *)malloc(sizeof(double) * length);
    /* Fix 2: NULL check was missing here */
    if (pfG == (double *)NULL) {
        fprintf(stderr, "GaborGetGaussianArray(): malloc failed\n");
        return (double *)NULL;
    }

    v = pfG;
    for (n = 0; n < nL-1; n++) {
        d = (double)((1+(1<<(n<<1))) * (1<<(maxlh<<1)));
        for (k = 0; k < pnAep[n]; k++) {
            t   = (double)k * (double)k / d;
            *v++ = exp(-M_PI * t);
        }
    }

    return pfG;
}

/*
 * GaborGetCE1 — allocate and fill the first complex exponential lattice.
 *   Stores cos(2*pi*n/N) in pfCE1[0..N-1]
 *   and    sin(2*pi*n/N) in pfCE1[N..2N-1].
 *
 * Returns a newly allocated array of size 2N. Caller must free it.
 *
 * Fix 3: perror() → fprintf+return NULL on alloc failure.
 * Fix 7: sincos() replaces separate cos()/sin() calls in the init loop.
 */
double *GaborGetCE1(int nL)
{
    double *pfCE1;
    int n, N;
    double c;

    N     = 1<<nL;
    pfCE1 = (double *)malloc(sizeof(double) * (N<<1));
    /* Fix 3: perror() → fprintf+return NULL */
    if (pfCE1 == (double *)NULL) {
        fprintf(stderr, "GaborGetCE1(): malloc failed\n");
        return (double *)NULL;
    }

    c = 2.0 * M_PI / N;
    /* Fix 7: sincos() halves the trig work in this N-element init loop */
    for (n = 0; n < N; n++)
        sincos(c * (double)n, &pfCE1[N+n], &pfCE1[n]);

    return pfCE1;
}

/*
 * GaborGetCE2 — allocate and fill the second complex exponential lattice.
 *   Stores cos and sin of k*2^n/(1+2^(2n)) scaled by 2*pi/N.
 *
 * Returns a newly allocated array. Caller must free it.
 *
 * Fix 3: perror() → fprintf+return NULL on alloc failure.
 * Fix 8: sincos() replaces separate cos()/sin() calls in the init loop.
 */
double *GaborGetCE2(int nL)
{
    double *pfCE2;
    int length, d1, n, k, N;
    double t, d, *vi, *vr, c;

    length = (1<<nL) - 2 + ((nL-1)<<1);
    pfCE2  = (double *)malloc(sizeof(double) * length);
    /* Fix 3: perror() → fprintf+return NULL */
    if (pfCE2 == (double *)NULL) {
        fprintf(stderr, "GaborGetCE2(): malloc failed\n");
        return (double *)NULL;
    }

    N  = 1<<nL;
    c  = 2.0 * M_PI / N;
    vr = pfCE2;
    for (n = 0; n < nL-1; n++) {
        d1 = 1<<n;
        d  = (double)(1 + (d1<<n));
        vi = vr + d1 + 1;
        /* Fix 8: sincos() for each (cos, sin) pair */
        for (k = 0; k <= d1; k++) {
            t = (double)(k * d1) / d;
            sincos(c * t, vi++, vr++);
        }
        vr = vi;
    }

    return pfCE2;
}

/*
 * GaborGetCE3 — allocate and fill the third complex exponential lattice.
 *   Stores cos and sin of k/(1+2^(2n)) scaled by 2*pi/N.
 *
 * Returns a newly allocated array. Caller must free it.
 *
 * Fix 3: perror() → fprintf+return NULL on alloc failure.
 * Fix 9: sincos() replaces separate cos()/sin() calls in the init loop.
 */
double *GaborGetCE3(int nL)
{
    double *pfCE3;
    int length, d1, n, k, N;
    double t, d, *vi, *vr, c;

    length = (1<<nL) - 2;
    pfCE3  = (double *)malloc(sizeof(double) * length);
    /* Fix 3: perror() → fprintf+return NULL */
    if (pfCE3 == (double *)NULL) {
        fprintf(stderr, "GaborGetCE3(): malloc failed\n");
        return (double *)NULL;
    }

    N  = 1<<nL;
    c  = 2.0 * M_PI / N;
    vr = pfCE3;
    for (n = 0; n < nL-1; n++) {
        d1 = 1<<n;
        d  = (double)(1 + (d1<<n));
        vi = vr + d1;
        /* Fix 9: sincos() for each (cos, sin) pair */
        for (k = 0; k < d1; k++) {
            t = (double)k / d;
            sincos(c * t, vi++, vr++);
        }
        vr = vi;
    }

    return pfCE3;
}

/*
 * GaborGetCArray — allocate and fill the C coefficient array.
 *   pfC[n] = sqrt(2^n / ((1+2^(2n)) * 2*pi))
 *
 * Returns a newly allocated array of size nL-1. Caller must free it.
 *
 * Fix 3: perror() → fprintf+return NULL on alloc failure.
 */
double *GaborGetCArray(int nL)
{
    double *pfC;
    int n;

    pfC = (double *)malloc(sizeof(double) * (nL-1));
    /* Fix 3: perror() → fprintf+return NULL */
    if (pfC == (double *)NULL) {
        fprintf(stderr, "GaborGetCArray(): malloc failed\n");
        return (double *)NULL;
    }

    for (n = 0; n < nL-1; n++)
        pfC[n] = sqrt((double)(1<<n) /
                      ((double)(1+(1<<(n<<1))) * 2.0 * M_PI));

    return pfC;
}

/*
 * GaborGetCoeff — compute the normalisation coefficient for two waveforms.
 *
 * Inputs:
 *	pfGNorm	discrete norm of the Gaussian for each octave (double *)
 *	j1, j2	octaves of the two waveforms (int)
 *	pfC	C coefficient array from GaborGetCArray() (double *)
 *
 * Returns the coefficient C = pfGNorm[j1] * pfGNorm[j2] * pfC[|j2-j1|].
 */
double GaborGetCoeff(double *pfGNorm, int j1, int j2, double *pfC)
{
    double k1, k2;
    int n;

    k1 = pfGNorm[j1];
    k2 = pfGNorm[j2];
    n  = abs(j2-j1);

    return k1 * k2 * pfC[n];
}

/*
 * GaborGetBound — compute the index bound for the inner product sum.
 *
 * Input:
 *	epslon	precision (double)
 *
 * Returns nBound = floor(sqrt(-log(epslon)/pi)) + 1.
 */
int GaborGetBound(double epslon)
{
    return (int)sqrt(-log((double)epslon) / M_PI) + 1;
}

/*
 * GaborGetB — allocate and fill the bandwidth bound array.
 *   bound[i] = sqrt((1 + 2^(2i)) * (-log(epslon)/pi))
 *
 * Returns a newly allocated array of size nL-1. Caller must free it.
 *
 * Fix 3: perror() → fprintf+return NULL on alloc failure.
 */
double *GaborGetB(int nL, double epslon)
{
    double *bound;
    double c;
    int i;

    bound = (double *)malloc(sizeof(double) * (nL-1));
    /* Fix 3: perror() → fprintf+return NULL */
    if (bound == (double *)NULL) {
        fprintf(stderr, "GaborGetB(): malloc failed\n");
        return (double *)NULL;
    }

    c = -log(epslon) / M_PI;
    for (i = 0; i < nL-1; i++)
        bound[i] = sqrt((double)(1+(1<<(i<<1))) * c);

    return bound;
}

/*
 * end of gb_corr.c
 */