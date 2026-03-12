/* ************************************************************************ */
/*		Operations on complex gabsignals				    */
/* 									    */
/* 	The structure ofcomplex gabsignal is GABSIGNAL, as defined in 	    */
/*	file gabsignals.h . 						    */
/* 	The size of a complex gabsignal is 2N, with N number of complex points */
/*									    */
/* 	The real part is stored in s->values[0], ..., s->values[N-1]	    */
/*	The imaginary part is  in  s->values[N], ..., s->values[2N-1]	    */
/*									    */
/* ************************************************************************ */
#include <stdio.h>
#include <math.h>
#include "mpp.h"

extern GABSIGNAL temporary;

/*
 *	functions defined in this file:
 *
 *	ComplexConjugate();		complex conjugate 	
 *	ComplexMultiplication();	Multiplication of 2 complex gabsignals  
 *	ComplexMulCoeff();		Multiplication of a gabsignal by a scalar 
 *	ComplexSubstraction();		Substraction of 2 complex gabsignals
 *	ComplexAddition();		Addition	" 	"
 *	ComplexLinComb();		Linear Combination of 2 complex gabsignals
 *
 *	CreateExponential();		Create function e^{-ix omega}
 *
 *	Real2Complex();   change a real gabsignal into a complex one
 *	Complex2Real();   Extract the real part of the gabsignal
 *
 */
/* ************************************************************************ */
/* 		Complex Conjugate of a gabsignal				    */
/* 									    */
/*	in: GABSIGNAL, *GABSIGNAL						    */
/*	out: 								    */
/* ************************************************************************ */
void ComplexConjugate(GABSIGNAL sigin, GABSIGNAL *psigout)
{
int index;
int HalfSize;

void change_gabsignal(GABSIGNAL gabsignal, int size);

if ( sigin->size %2 != 0 )
	perror("ComplexConjugate: size of complex gabsignal should be even");

HalfSize = sigin->size/2;

/* sigout allocated ? */
if ( *psigout == NULL )
	*psigout = new_struct_gabsignal(); /* allocation of sigout */

/* sigout = sigin ? */
if ( *psigout != sigin )
	change_gabsignal(*psigout, sigin->size);

for (index = 0; index <=  HalfSize-1; index ++)
	{
	(*psigout)->values[index] = sigin->values[index];
	(*psigout)->values[HalfSize + index] = -1.*sigin->values[HalfSize + index];
	}
return;
} /* end of ComplexConjugate */
/* ************************************************************************ */
/* 		Complex multiplication  of two gabsignals			    */
/* 									    */
/*	in: GABSIGNAL, GABSIGNAL, *GABSIGNAL					    */
/*	out: 								    */
/* ************************************************************************ */
void ComplexMultiplication(GABSIGNAL sig1, GABSIGNAL sig2, GABSIGNAL *psigout)
{
int index;
int HalfSize;
double temp; /* very important temporary variable */

void change_gabsignal(GABSIGNAL gabsignal, int size);

if ( (sig1 == NULL) || (sig2 == NULL) )
	perror("ComplexMultiplication: input gabsignal not allocated");

if ( sig1->size != sig2->size )
	perror("ComplexMultiplication: size the gabsignals must be equal");

if (( sig1->size %2 != 0 ) || ( sig2->size %2 != 0 ))
	perror("ComplexMultiplication: size of complex gabsignal should be even");

/* sigout allocated ? */
if ( *psigout == NULL )
	*psigout = new_struct_gabsignal(); /* allocation of sigout */

/* sigout = sig1 or sig2 ? */
if ( (*psigout != sig1) && (*psigout != sig2) )
	change_gabsignal(*psigout, sig1->size);

HalfSize = sig1->size/2;

for (index = 0; index <=  HalfSize-1; index ++)
	{
	/* Real part ( saved in temp ) 			*/
	/* Does not affect (*psigout)->values[index] 		*/
	/* Which is used later					*/
	/* Avoids a big problem when *psigout = sig1 or sig2 !! */

	temp = sig1->values[index] * sig2->values[index];
	temp -= sig1->values[HalfSize + index] * sig2->values[HalfSize + index];
	
	/* Imaginary part */
	(*psigout)->values[HalfSize + index] = sig1->values[index] * sig2->values[HalfSize + index] + 	sig1->values[HalfSize + index] *sig2->values[index];
	(*psigout)->values[index] = temp;
	}
return;
} /* end of ComplexMultiplication */
/* ************************************************************************ */
/* 		Complex Multiplication by a scalar of a gabsignal		    */
/* 									    */
/*	in: GABSIGNAL, double, *GABSIGNAL					    */
/*	out: 							    */
/* ************************************************************************ */
void ComplexMulCoeff(GABSIGNAL sigin, double coeff, GABSIGNAL *psigout)
{
int index;

void change_gabsignal(GABSIGNAL gabsignal, int size);

if (sigin == NULL)
	perror("ComplexMulCoeff: input gabsignal not allocated");

if ( sigin->size %2 != 0 )
	perror("ComplexMulCoeff: size of complex gabsignal should be even");

/* sigout allocated ? */
if ( *psigout == NULL )
	*psigout = new_struct_gabsignal(); /* allocation of sigout */

if(*psigout != sigin) change_gabsignal(*psigout, sigin->size);


for (index = 0; index <=  sigin->size - 1; index ++)
	(*psigout)->values[index] = coeff * sigin->values[index];


return;
} /* end of ComplexMulCoeff */
/* ************************************************************************ */
/* 		Complex substraction  of two gabsignals			    */
/* 									    */
/*	in: GABSIGNAL, GABSIGNAL, *GABSIGNAL					    */
/*	out: 								    */
/* ************************************************************************ */
//GABSIGNAL ComplexSubstraction(GABSIGNAL sig1, GABSIGNAL sig2, GABSIGNAL *psigout)
void ComplexSubstraction(GABSIGNAL sig1, GABSIGNAL sig2, GABSIGNAL *psigout)
{
int index;

void change_gabsignal(GABSIGNAL gabsignal, int size);

if ( (sig1 == NULL) || (sig2 == NULL) )
	perror("ComplexMultiplication: input gabsignal not allocated");

if ( sig1->size != sig2->size )
	perror("ComplexSubstraction: size the gabsignals must be equal");

if (( sig1->size %2 != 0 ) || ( sig2->size %2 != 0 ))
	perror("ComplexSubstraction: size of complex gabsignal should be even");

/* sigout allocated ? */
if ( *psigout == NULL )
	*psigout = new_struct_gabsignal(); /* allocation of sigout */

/* sigout = sig1 or sig2 ? */
if ( (*psigout != sig1) && (*psigout != sig2) )
	change_gabsignal(*psigout, sig1->size);


for (index = 0; index <=  sig1->size-1; index ++)
	(*psigout)->values[index] = sig1->values[index] - sig2->values[index];

return;
} /* end of ComplexSubstraction */
/* ************************************************************************ */
/* 		Complex addition  of two gabsignals			    */
/* 									    */
/*	in: GABSIGNAL, GABSIGNAL						    */
/*	out: GABSIGNAL							    */
/* ************************************************************************ */
//GABSIGNAL ComplexAddition(GABSIGNAL sig1, GABSIGNAL sig2, GABSIGNAL *psigout)
void ComplexAddition(GABSIGNAL sig1, GABSIGNAL sig2, GABSIGNAL *psigout)
{
int index;

void change_gabsignal(GABSIGNAL gabsignal, int size);

if ( (sig1 == NULL) || (sig2 == NULL) )
	perror("ComplexMultiplication: input gabsignal not allocated");

if ( sig1->size != sig2->size )
	perror("ComplexAddition: size the gabsignals must be equal");

if (( sig1->size %2 != 0 ) || ( sig2->size %2 != 0 ))
	perror("ComplexAddition: size of complex gabsignal should be even");

/* sigout allocated ? */
if ( *psigout == NULL )
	*psigout = new_struct_gabsignal(); /* allocation of sigout */

/* sigout = sig1 or sig2 ? */
if ( (*psigout != sig1) && (*psigout != sig2) )
	change_gabsignal(*psigout, sig1->size);


for (index = 0; index <=  sig1->size-1; index ++)
	(*psigout)->values[index] = sig1->values[index] + sig2->values[index];

return;
} /* end of ComplexAddition */
/* ************************************************************************ */
/* 		Linear combination  of two gabsignals			    */
/* 									    */
/*	in: double, GABSIGNAL, double, GABSIGNAL, *GABSIGNAL			    */
/*	out: 								    */
/* 									    */
/*	computes lambda * sig1 + mu * sig2				    */
/* ************************************************************************ */
void ComplexLinComb(double lambda, GABSIGNAL sig1, double mu, GABSIGNAL sig2, GABSIGNAL *psigout)
{
int index;
void change_gabsignal(GABSIGNAL gabsignal, int size);

if ( (sig1 == NULL) || (sig2 == NULL) )
	perror("ComplexMultiplication: input gabsignal not allocated");

if ( sig1->size != sig2->size )
	perror("ComplexLinComb: size the gabsignals must be equal");

if (( sig1->size %2 != 0 ) || ( sig2->size %2 != 0 ))
	perror("ComplexLinComb: size of complex gabsignal should be even");

/* sigout allocated ? */
if ( *psigout == NULL )
	*psigout = new_struct_gabsignal(); /* allocation of sigout */

/* sigout = sig1 or sig2 ? */
if ( (*psigout != sig1) && (*psigout != sig2) )
	change_gabsignal(*psigout, sig1->size);


for (index = 0; index <=  sig1->size-1; index ++)
	(*psigout)->values[index] = (lambda * sig1->values[index]) + (mu * sig2->values[index]);

} /* end of ComplexLinComb */
/* ************************************************************************ */
/* 			Create Exponential function			    */
/* 									    */
/*	in: double, double, double, int, *GABSIGNAL				    */
/*	out: 								    */
/* 									    */
/*	computes exp(i*position*omega)					    */
/* 	between min and max - scale					    */
/*	where scale = (max - min)/size					    */
/*	The length of the gabsignal is 2*size, because it is complex	    */
/* ************************************************************************ */
void CreateExponential(double position, double min, double max, int size, GABSIGNAL *psigout)
//int size; /* size of the real gabsignal */
{
int index;
double x, scale;
void change_gabsignal(GABSIGNAL gabsignal, int size);

/* sigout allocated ? */
if ( *psigout == NULL ){
	*psigout = new_struct_gabsignal(); /* allocation of sigout */
}
change_gabsignal((*psigout), 2*size);
(*psigout)->scale = (double) (max-min)/size;
scale = (max - min)/(double)(size);


for ( index = 0; index <= size -1 ; index++ )
	{
	x = scale * (double)index + min;
	/* real part */
	(*psigout)->values[index] = (double) cos((double) (position * x) );

	/* imaginary part */
	(*psigout)->values[size + index] = (double) sin((double) (position * x) );

	}

return;
} /* end of CreateExponential */
/* ************************************************************************ */
/* 	Transformation of a real gabsignal into a complex one		    */
/* 									    */
/*	in: GABSIGNAL, *GABSIGNAL						    */
/*	out: 								    */
/* ************************************************************************ */
void Real2Complex(GABSIGNAL sigin, GABSIGNAL *psigout)
{
int i;
void change_gabsignal(GABSIGNAL gabsignal, int size);
int sig_copy(GABSIGNAL input,GABSIGNAL output);

/* sigout allocated ? */
if ( *psigout == NULL )
	*psigout = new_struct_gabsignal(); /* allocation of sigout */

if ( *psigout != sigin)
	{
	change_gabsignal(*psigout, 2 * sigin->size);

	/* real part */
	for ( i = 0; i <= sigin->size - 1; i ++ )
		(*psigout)->values[i] = sigin->values[i];
	/* imaginary part */
	for ( i = 0; i <= sigin->size - 1; i ++ )
		(*psigout)->values[i + sigin->size] = 0;
	}


if ( *psigout == sigin )
	{
	change_gabsignal(temporary, 2 * sigin->size);

	/* real part */
	for ( i = 0; i <= sigin->size - 1; i ++ )
		temporary->values[i] = sigin->values[i];
	/* imaginary part */
	for ( i = 0; i <= sigin->size - 1; i ++ )
		temporary->values[i + sigin->size] = 0;

	sig_copy(temporary, *psigout);
	}
} /* end of Real2Complex */
/* ************************************************************************ */
/* 			Real part of a complex gabsignal			    */
/* 									    */
/*	in: GABSIGNAL, *GABSIGNAL						    */
/*	out: 								    */
/* ************************************************************************ */
void Complex2Real(GABSIGNAL sigin, GABSIGNAL *psigout)
{
int i;
int Size;
void change_gabsignal(GABSIGNAL gabsignal, int size);
int sig_copy(GABSIGNAL input,GABSIGNAL output);

/* sigout allocated ? */
if ( *psigout == NULL )
	*psigout = new_struct_gabsignal(); /* allocation of sigout */

Size = sigin->size / 2;

if ( *psigout != sigin)
	{
		change_gabsignal(*psigout, Size);

	/* real part */
	for ( i = 0; i <= Size - 1; i ++ )
		(*psigout)->values[i] = sigin->values[i];
	}

if ( *psigout == sigin )
	{
	change_gabsignal(temporary,  Size);

	/* real part */
	for ( i = 0; i <= Size - 1; i ++ )
		temporary->values[i] = sigin->values[i];

	sig_copy(temporary, *psigout);
	}
} /* end of Complex2Real */
/*
 * take the imaginary part of the complex gabsignal
 */
void Complex2Imaginary(GABSIGNAL sigin, GABSIGNAL *psigout)
{
int i;
int Size;
void change_gabsignal(GABSIGNAL gabsignal, int size);
int sig_copy(GABSIGNAL input,GABSIGNAL output);

/* sigout allocated ? */
if ( *psigout == NULL )
	*psigout = new_struct_gabsignal(); /* allocation of sigout */

Size = sigin->size / 2;

if ( *psigout != sigin)
	{
		change_gabsignal(*psigout, Size);

	/* real part */
	for ( i = Size; i <= sigin->size - 1; i ++ )
		(*psigout)->values[i-Size] = sigin->values[i];
	}

if ( *psigout == sigin )
	{
	change_gabsignal(temporary,  Size);

	/* real part */
	for ( i = Size; i <= sigin->size - 1; i ++ )
		temporary->values[i-Size] = sigin->values[i];

	sig_copy(temporary, *psigout);
	}
} /* end of Complex2Real */
/*
 * taking complex modulation from a complex
 * gabsignal
 */
void ComplexModulation(GABSIGNAL sigin,GABSIGNAL sigout)
{
    int i, size;
	void change_gabsignal(GABSIGNAL gabsignal, int size);
	
    if (sigin == (GABSIGNAL)NULL || sigout == (GABSIGNAL)NULL)
	perror("ComplexModulation(): null argument!");

    size = sigin->size/2;

    if (sigout->size != size)
	change_gabsignal(sigout,size);

    for (i=0;i<size;i++)
	sigout->values[i] = sqrt((double)(sigin->values[i]*sigin->values[i]+
		sigin->values[i+size]*sigin->values[i+size]));
}
/*
 * multiply a complex gabsignal by a complex number
 */
void ComplexMulCN(GABSIGNAL gabsignal,double v_real,double v_imag)
{
    int i, size;
    double temp;

    if (gabsignal == (GABSIGNAL)NULL)
	perror("ComplexMulCN(): null input!");

    size = gabsignal->size/2;
    for (i=0;i<size;i++)
	{
	temp = v_real*gabsignal->values[i] -
			v_imag*gabsignal->values[i+size];
	gabsignal->values[i+size] = v_imag*gabsignal->values[i] +
			v_real*gabsignal->values[i+size];
	gabsignal->values[i] = temp;
	}
}
/*
 * searching the maximum modula from a complex array
 *
 * The first half of the array contains the real part and the second
 * half of the array contains the imaginary part
 *
 * Inputs:
 *	value		the complex array (double *)
 *	size		size of the complex array (int)
 *	size_op		size for the operation (int)
 *
 * Outputs:
 *	modula		the maximum modula square found (double *)
 *	v_r		the real part of the maximum modula found (double *)
 *	v_i		the imaginary part of the maximum modula 
 *			found (double *)
 *	index		the index of the maximum modula (int)
 *
 */
void complex_array_max(double *value,int size,int size_op,double *modula,double *v_r,double *v_i,int *index)
{
    int i;
    double m, *p_r, *p_i;

    *modula = 0.0;
    *v_r = 0.0;
    *v_i = 0.0;
    *index = 0;
    p_r = value;
    p_i = value+size;
    for (i=0;i<size_op;i++)
	{
	m = (*p_r)*(*p_r)+(*p_i)*(*p_i);
	if (m>*modula)
	   {
	   *modula = m;
	   *v_r = *p_r;
	   *v_i = *p_i;
	   *index = i;
	   }
	p_r++;
	p_i++;
	}
}
/*
 * end of complex_op.c
 */
