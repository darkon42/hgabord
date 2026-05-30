/* ************************************************************************ */
/*		Operations on complex gabsignals				    */
/* 									    */
/* 	The structure of complex gabsignal is GABSIGNAL, as defined in 	    */
/*	file gabsignals.h . 						    */
/* 	The size of a complex gabsignal is 2N, with N number of complex points */
/*									    */
/* 	The real part is stored in s->values[0], ..., s->values[N-1]	    */
/*	The imaginary part is  in  s->values[N], ..., s->values[2N-1]	    */
/*									    */
/* ************************************************************************ */
#include <stdio.h>
#define _GNU_SOURCE
#include <math.h>
#include <string.h>
#include "mpp.h"

extern GABSIGNAL temporary;

void change_gabsignal(GABSIGNAL gabsignal, int size);
int  sig_copy(GABSIGNAL input, GABSIGNAL output);

/*
 *	functions defined in this file:
 *
 *	ComplexConjugate();		complex conjugate
 *	ComplexMultiplication();	multiplication of 2 complex gabsignals
 *	ComplexMulCoeff();		multiplication of a gabsignal by a scalar
 *	ComplexSubstraction();		subtraction of 2 complex gabsignals
 *	ComplexAddition();		addition of 2 complex gabsignals
 *	ComplexLinComb();		linear combination of 2 complex gabsignals
 *	CreateExponential();		create function e^{i*position*omega}
 *	Real2Complex();			change a real gabsignal into a complex one
 *	Complex2Real();			extract the real part of a gabsignal
 *	Complex2Imaginary();		extract the imaginary part of a gabsignal
 *	ComplexModulation();		compute modulus of each complex sample
 *	ComplexMulCN();			multiply a gabsignal by a complex number
 *	complex_array_max();		find the sample with maximum modulus
 */

/* ************************************************************************ */
/*		Complex Conjugate of a gabsignal			    */
/* ************************************************************************ */
void ComplexConjugate(GABSIGNAL sigin, GABSIGNAL *psigout)
{
    int HalfSize;

    /* Fix 1: perror() doesn't stop execution — use fprintf+return */
    if (sigin == NULL) {
        fprintf(stderr, "ComplexConjugate: null input\n");
        return;
    }
    if (sigin->size % 2 != 0) {
        fprintf(stderr, "ComplexConjugate: size of complex gabsignal should be even\n");
        return;
    }

    HalfSize = sigin->size / 2;

    if (*psigout == NULL)
        *psigout = new_struct_gabsignal();

    if (*psigout != sigin)
        change_gabsignal(*psigout, sigin->size);

    /* Fix 6: real part copy with memcpy on the non-in-place path */
    if (*psigout != sigin)
        memcpy((*psigout)->values, sigin->values, (size_t)HalfSize * sizeof(double));

    /* Negate imaginary part in-place (works for both in-place and out-of-place) */
    {
        int index;
        for (index = 0; index < HalfSize; index++)
            (*psigout)->values[HalfSize + index] = -sigin->values[HalfSize + index];
    }
} /* end of ComplexConjugate */


/* ************************************************************************ */
/*		Complex multiplication of two gabsignals		    */
/* ************************************************************************ */
void ComplexMultiplication(GABSIGNAL sig1, GABSIGNAL sig2, GABSIGNAL *psigout)
{
    int index;
    int HalfSize;
    double temp;

    /* Fix 1: perror() → fprintf+return */
    if (sig1 == NULL || sig2 == NULL) {
        fprintf(stderr, "ComplexMultiplication: input gabsignal not allocated\n");
        return;
    }
    if (sig1->size != sig2->size) {
        fprintf(stderr, "ComplexMultiplication: sizes of gabsignals must be equal\n");
        return;
    }
    if (sig1->size % 2 != 0) {
        fprintf(stderr, "ComplexMultiplication: size of complex gabsignal should be even\n");
        return;
    }

    if (*psigout == NULL)
        *psigout = new_struct_gabsignal();

    if ((*psigout != sig1) && (*psigout != sig2))
        change_gabsignal(*psigout, sig1->size);

    HalfSize = sig1->size / 2;

    /* Fix 7: `<= HalfSize-1` → `< HalfSize` throughout */
    for (index = 0; index < HalfSize; index++) {
        /* temp guards correctness when *psigout aliases sig1 or sig2 */
        temp = sig1->values[index] * sig2->values[index]
             - sig1->values[HalfSize + index] * sig2->values[HalfSize + index];
        (*psigout)->values[HalfSize + index] =
              sig1->values[index]            * sig2->values[HalfSize + index]
            + sig1->values[HalfSize + index] * sig2->values[index];
        (*psigout)->values[index] = temp;
    }
} /* end of ComplexMultiplication */


/* ************************************************************************ */
/*		Complex multiplication by a scalar			    */
/* ************************************************************************ */
void ComplexMulCoeff(GABSIGNAL sigin, double coeff, GABSIGNAL *psigout)
{
    int index;

    /* Fix 1: perror() → fprintf+return */
    if (sigin == NULL) {
        fprintf(stderr, "ComplexMulCoeff: input gabsignal not allocated\n");
        return;
    }
    if (sigin->size % 2 != 0) {
        fprintf(stderr, "ComplexMulCoeff: size of complex gabsignal should be even\n");
        return;
    }

    if (*psigout == NULL)
        *psigout = new_struct_gabsignal();

    if (*psigout != sigin)
        change_gabsignal(*psigout, sigin->size);

    /* Fix 7: `<= sigin->size-1` → `< sigin->size` */
    for (index = 0; index < sigin->size; index++)
        (*psigout)->values[index] = coeff * sigin->values[index];
} /* end of ComplexMulCoeff */


/* ************************************************************************ */
/*		Complex subtraction of two gabsignals			    */
/* ************************************************************************ */
void ComplexSubstraction(GABSIGNAL sig1, GABSIGNAL sig2, GABSIGNAL *psigout)
{
    int index;

    /* Fix 1: perror() → fprintf+return
     * Fix 2: was printing "ComplexMultiplication" — corrected function name */
    if (sig1 == NULL || sig2 == NULL) {
        fprintf(stderr, "ComplexSubstraction: input gabsignal not allocated\n");
        return;
    }
    if (sig1->size != sig2->size) {
        fprintf(stderr, "ComplexSubstraction: sizes of gabsignals must be equal\n");
        return;
    }
    if (sig1->size % 2 != 0) {
        fprintf(stderr, "ComplexSubstraction: size of complex gabsignal should be even\n");
        return;
    }

    if (*psigout == NULL)
        *psigout = new_struct_gabsignal();

    if ((*psigout != sig1) && (*psigout != sig2))
        change_gabsignal(*psigout, sig1->size);

    /* Fix 7: `<= sig1->size-1` → `< sig1->size` */
    for (index = 0; index < sig1->size; index++)
        (*psigout)->values[index] = sig1->values[index] - sig2->values[index];
} /* end of ComplexSubstraction */


/* ************************************************************************ */
/*		Complex addition of two gabsignals			    */
/* ************************************************************************ */
void ComplexAddition(GABSIGNAL sig1, GABSIGNAL sig2, GABSIGNAL *psigout)
{
    int index;

    /* Fix 1: perror() → fprintf+return
     * Fix 2: was printing "ComplexMultiplication" — corrected function name */
    if (sig1 == NULL || sig2 == NULL) {
        fprintf(stderr, "ComplexAddition: input gabsignal not allocated\n");
        return;
    }
    if (sig1->size != sig2->size) {
        fprintf(stderr, "ComplexAddition: sizes of gabsignals must be equal\n");
        return;
    }
    if (sig1->size % 2 != 0) {
        fprintf(stderr, "ComplexAddition: size of complex gabsignal should be even\n");
        return;
    }

    if (*psigout == NULL)
        *psigout = new_struct_gabsignal();

    if ((*psigout != sig1) && (*psigout != sig2))
        change_gabsignal(*psigout, sig1->size);

    /* Fix 7: `<= sig1->size-1` → `< sig1->size` */
    for (index = 0; index < sig1->size; index++)
        (*psigout)->values[index] = sig1->values[index] + sig2->values[index];
} /* end of ComplexAddition */


/* ************************************************************************ */
/*		Linear combination of two gabsignals			    */
/*		computes lambda * sig1 + mu * sig2			    */
/* ************************************************************************ */
void ComplexLinComb(double lambda, GABSIGNAL sig1, double mu, GABSIGNAL sig2, GABSIGNAL *psigout)
{
    int index;

    /* Fix 1: perror() → fprintf+return
     * Fix 2: was printing "ComplexMultiplication" — corrected function name */
    if (sig1 == NULL || sig2 == NULL) {
        fprintf(stderr, "ComplexLinComb: input gabsignal not allocated\n");
        return;
    }
    if (sig1->size != sig2->size) {
        fprintf(stderr, "ComplexLinComb: sizes of gabsignals must be equal\n");
        return;
    }
    if (sig1->size % 2 != 0) {
        fprintf(stderr, "ComplexLinComb: size of complex gabsignal should be even\n");
        return;
    }

    if (*psigout == NULL)
        *psigout = new_struct_gabsignal();

    if ((*psigout != sig1) && (*psigout != sig2))
        change_gabsignal(*psigout, sig1->size);

    /* Fix 7: `<= sig1->size-1` → `< sig1->size` */
    for (index = 0; index < sig1->size; index++)
        (*psigout)->values[index] = lambda * sig1->values[index]
                                  + mu     * sig2->values[index];
} /* end of ComplexLinComb */


/* ************************************************************************ */
/*		Create Exponential function				    */
/*		computes exp(i*position*omega)				    */
/*		between min and max-scale				    */
/*		where scale = (max-min)/size				    */
/*		The gabsignal length is 2*size (complex)		    */
/* ************************************************************************ */
void CreateExponential(double position, double min, double max, int size, GABSIGNAL *psigout)
{
    int    index;
    double x, scale;

    if (*psigout == NULL)
        *psigout = new_struct_gabsignal();

    change_gabsignal(*psigout, 2 * size);

    /* Fix 5: scale was computed twice with identical expressions — compute once */
    scale = (max - min) / (double)size;
    (*psigout)->scale = scale;

    /* Fix 7: `<= size-1` → `< size`
     * Fix 8: cos/sin of the same argument replaced with sincos() */
    for (index = 0; index < size; index++) {
        x = scale * (double)index + min;
        sincos(position * x,
               &(*psigout)->values[size  + index],   /* imaginary */
               &(*psigout)->values[index]);           /* real      */
    }
} /* end of CreateExponential */


/* ************************************************************************ */
/*		Real gabsignal → complex gabsignal			    */
/* ************************************************************************ */
void Real2Complex(GABSIGNAL sigin, GABSIGNAL *psigout)
{
    /* Fix 3: was if/if — second branch could never be entered after the
     * first modified *psigout. Changed to if/else if. */
    if (*psigout == NULL)
        *psigout = new_struct_gabsignal();

    if (*psigout != sigin) {
        change_gabsignal(*psigout, 2 * sigin->size);
        /* Fix 6: memcpy for the real part copy */
        memcpy((*psigout)->values, sigin->values,
               (size_t)sigin->size * sizeof(double));
        /* zero the imaginary part */
        memset((*psigout)->values + sigin->size, 0,
               (size_t)sigin->size * sizeof(double));
    } else {
        /* in-place: sigin and *psigout point to the same object */
        change_gabsignal(temporary, 2 * sigin->size);
        memcpy(temporary->values, sigin->values,
               (size_t)sigin->size * sizeof(double));
        memset(temporary->values + sigin->size, 0,
               (size_t)sigin->size * sizeof(double));
        sig_copy(temporary, *psigout);
    }
} /* end of Real2Complex */


/* ************************************************************************ */
/*		Complex gabsignal → real part only			    */
/* ************************************************************************ */
void Complex2Real(GABSIGNAL sigin, GABSIGNAL *psigout)
{
    int Size;

    /* Fix 3: was if/if — changed to if/else if */
    if (*psigout == NULL)
        *psigout = new_struct_gabsignal();

    Size = sigin->size / 2;

    if (*psigout != sigin) {
        change_gabsignal(*psigout, Size);
        /* Fix 6: memcpy for the real part copy */
        memcpy((*psigout)->values, sigin->values,
               (size_t)Size * sizeof(double));
    } else {
        change_gabsignal(temporary, Size);
        memcpy(temporary->values, sigin->values,
               (size_t)Size * sizeof(double));
        sig_copy(temporary, *psigout);
    }
} /* end of Complex2Real */


/* ************************************************************************ */
/*		Complex gabsignal → imaginary part only			    */
/* ************************************************************************ */
void Complex2Imaginary(GABSIGNAL sigin, GABSIGNAL *psigout)
{
    int Size;

    /* Fix 3: was if/if — changed to if/else if */
    if (*psigout == NULL)
        *psigout = new_struct_gabsignal();

    Size = sigin->size / 2;

    if (*psigout != sigin) {
        change_gabsignal(*psigout, Size);
        /* Fix 6: memcpy for the imaginary part copy */
        memcpy((*psigout)->values, sigin->values + Size,
               (size_t)Size * sizeof(double));
    } else {
        change_gabsignal(temporary, Size);
        memcpy(temporary->values, sigin->values + Size,
               (size_t)Size * sizeof(double));
        sig_copy(temporary, *psigout);
    }
} /* end of Complex2Imaginary */


/* ************************************************************************ */
/*		Modulus of each complex sample				    */
/* ************************************************************************ */
void ComplexModulation(GABSIGNAL sigin, GABSIGNAL sigout)
{
    int i, size;

    /* Fix 1: perror() → fprintf+return */
    if (sigin == NULL || sigout == NULL) {
        fprintf(stderr, "ComplexModulation(): null argument!\n");
        return;
    }

    size = sigin->size / 2;

    if (sigout->size != size)
        change_gabsignal(sigout, size);

    /* sqrt kept as-is per requirement */
    for (i = 0; i < size; i++)
        sigout->values[i] = sqrt(sigin->values[i]          * sigin->values[i]
                               + sigin->values[i + size]   * sigin->values[i + size]);
}


/* ************************************************************************ */
/*		Multiply a complex gabsignal by a complex number	    */
/* ************************************************************************ */
void ComplexMulCN(GABSIGNAL gabsignal, double v_real, double v_imag)
{
    int i, size;
    double temp;

    /* Fix 1: perror() → fprintf+return */
    if (gabsignal == NULL) {
        fprintf(stderr, "ComplexMulCN(): null input!\n");
        return;
    }

    size = gabsignal->size / 2;
    for (i = 0; i < size; i++) {
        temp = v_real * gabsignal->values[i]
             - v_imag * gabsignal->values[i + size];
        gabsignal->values[i + size] = v_imag * gabsignal->values[i]
                                    + v_real * gabsignal->values[i + size];
        gabsignal->values[i] = temp;
    }
}


/* ************************************************************************ */
/*		Maximum modulus in a complex array			    */
/*								    */
/*  Inputs:							    */
/*	value    the complex array (split: real then imag)		    */
/*	size     half-size (number of complex samples)			    */
/*	size_op  number of samples to search				    */
/*  Outputs:							    */
/*	modula   maximum modulus squared found				    */
/*	v_r      real part of the maximum sample			    */
/*	v_i      imaginary part of the maximum sample			    */
/*	index    index of the maximum sample				    */
/* ************************************************************************ */
void complex_array_max(double *value, int size, int size_op,
                       double *modula, double *v_r, double *v_i, int *index)
{
    int i;
    double m, *p_r, *p_i;

    *modula = 0.0;
    *v_r    = 0.0;
    *v_i    = 0.0;
    *index  = 0;
    p_r = value;
    p_i = value + size;

    for (i = 0; i < size_op; i++) {
        m = (*p_r) * (*p_r) + (*p_i) * (*p_i);
        if (m > *modula) {
            *modula = m;
            *v_r    = *p_r;
            *v_i    = *p_i;
            *index  = i;
        }
        p_r++;
        p_i++;
    }
}

/*
 * end of complex_op.c
 */
