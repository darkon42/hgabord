/*..........................................................................*/
/*                                                                          */
/*      ------------------------------------------------------------*/
/*      (C) 1993 Copyright New York University, All Right Reserved.         */
/*                                                                          */
/*..........................................................................*/
/****************************************************************************/
/*                                                                          */
/*  math.c           Some useful mathematical functions                     */
/*                                                                          */
/****************************************************************************/
#include "mpp.h"

//#define M_LN2      0.693147180559945309417
#define M_2_SQRTPI 1.12837916709551257390
//#define M_SQRT1_2  0.707106781186547524401

#ifndef sun4
/* log in base 2 */
double log2(double x)
{
  double y;
  y= ((double)log(x))/ M_LN2;
  return(y);
} /* this function is defined (as double) in math.h on the sun4 */
#endif


/* 2^j  (j > 0) */
int iexp2(int j)
{
  return( 1 << j);
}

/* 2^j  (j can be >= 0 or < 0 ) */
double fexp2(int j)
{
  int k;
  double s;
  s = 1.0;
  if (j >= 0) {
    return( (double)(1 << j));
  }
  else {
    for (k = j; k < 0 ; k++)
      s /= 2.0;
    return(s);
  }
}



/*
 * gauss function g(x) = exp(-(a*(x-b))**2)
  */
double gaussian(double x,double a,double b)
{
   double y;

   y=(double)exp(-(double)(a*a*(x-b)*(x-b)));
   return(y);
}
/*
 *
 * compute the L2 normalized gaussian:
 *	 g(t)=1/sqrt(2*pi*sigma) exp(-t*t/(2*sigma*sigma)
 *
 * M_2_SQRTPI is 2/sqrt(pi) defined in math.h
 * M_SQRT1_2 is 1/sqrt(2) defined in math.h
 *
 *
 */
double gaussianL2(double t,double sigma)
{
    double y, tmp;

    tmp = sqrt(sigma);
    y = exp(-t*t/(2.0*sigma*sigma))*
		(double)(M_2_SQRTPI*M_SQRT1_2)/(2.0*tmp);
/*
    return((double)(exp(-t*t/(2.0*sigma*sigma))*
		(double)(M_2_SQRTPI*M_SQRT1_2)/(2.0*sqrt(sigma))));
*/
    return((double)y);
}
/*
 * singauss function g(x) = [exp(-(a*(x+b))**2) -
                             exp(-(a*(x - b))**2)] / 2.0
 */
double singauss(double x,double a,double b)
{
double y;

	y=(double)exp(-(double)(a*a*(x+b)*(x+b)));
	y -= (double)exp(-(double)(a*a*(x-b)*(x-b)));
	y /= 2.0;
return(y);
}



int find2power(int n)
{
   long m, m2;

   m = 0;
   m2 = 1<<m; /* 2 to the power of m */
   while (m2-n < 0) {
	m++;
	m2 <<= 1; /* m2 = m2*2 */
   }
   return(m);
}
#ifdef hp
int rint(r)
double r;
{
    return((int)r);
}
#endif
#ifndef sun4
int nint(double r)
{
    return((int)r);
}
#endif
#ifdef i386sco
int random()
{
    return((int)Irand48());
}
#endif
/*
 * end of math.c
 */
