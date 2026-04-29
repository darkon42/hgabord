/*..........................................................................*/
/*                                                                          */
/*      ------------------------------------------------------------        */
/*      (C) 1993 Copyright New York University, All Right Reserved.         */
/*                                                                          */
/*..........................................................................*/
/****************************************************************************/
/*                                                                          */
/*  math.c           Some useful mathematical functions                     */
/*                                                                          */
/*  Removed (all obsolete as of C99):                                       */
/*    log2()   — now standard in <math.h>; was conflicting with stdlib      */
/*    iexp2()  — trivial wrapper around (1<<j); never called                */
/*    fexp2()  — replaced by exp2() from <math.h>; never called            */
/*    rint()   — now standard in <math.h>; was wrong (truncated not round)  */
/*    nint()   — use (int)round() instead; never called                     */
/*    random() — now standard POSIX; SCO/i386 platform dead; never called   */
/*                                                                          */
/****************************************************************************/
//#include "mpp.h"
#include <math.h>

/* M_2_SQRTPI = 2/sqrt(pi) and M_SQRT1_2 = 1/sqrt(2) are both defined in
 * <math.h> on any C99/POSIX platform. The local defines are kept only as
 * a fallback for non-conforming environments. */
#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.12837916709551257390
#endif
#ifndef M_SQRT1_2
#define M_SQRT1_2  0.707106781186547524401
#endif

/*
 * gaussian — g(x) = exp(-(a*(x-b))^2)
 */
double gaussian(double x, double a, double b)
{
    double d = a * (x - b);
    return exp(-(d * d));
}

/*
 * gaussianL2 — L2-normalised Gaussian:
 *   g(t) = (1 / sqrt(2*pi*sigma)) * exp(-t^2 / (2*sigma^2))
 *
 * M_2_SQRTPI = 2/sqrt(pi), M_SQRT1_2 = 1/sqrt(2).
 * The commented-out duplicate inside the original is removed.
 */
double gaussianL2(double t, double sigma)
{
    double tmp = sqrt(sigma);
    return exp(-t*t / (2.0*sigma*sigma))
           * (M_2_SQRTPI * M_SQRT1_2) / (2.0 * tmp);
}

/*
 * singauss — g(x) = (exp(-(a*(x+b))^2) - exp(-(a*(x-b))^2)) / 2
 */
double singauss(double x, double a, double b)
{
    double y;
    y  = exp(-(a*a*(x+b)*(x+b)));
    y -= exp(-(a*a*(x-b)*(x-b)));
    return y / 2.0;
}

/*
 * find2power — return smallest m such that 2^m >= n.
 *
 * Fix 8: original used long for m2 which overflows for large n on 32-bit
 * platforms. Using unsigned long long makes the shift safe up to 2^63.
 * For n <= 0 returns 0 by convention.
 */
int find2power(int n)
{
    int m;
    unsigned long long m2;

    if (n <= 1) return 0;
    m  = 0;
    m2 = 1ULL;
    while ((long long)m2 - n < 0) {
        m++;
        m2 <<= 1;
    }
    return m;
}

/*
 * end of math.c
 */
 