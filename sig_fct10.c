/*..........................................................................*/
/*                                                                          */
/*      ------------------------------------------------------------        */
/*      (C) 1993 Copyright New York University, All Right Reserved.         */
/*      Modify by Zhifeng Zhang and Mike Orszag, 1992                       */
/*                                                                          */
/*..........................................................................*/
/****************************************************************************/
/*                                                                          */
/*  sig_fct10.c   Miscellaneous useful functions on gabsignals:             */
/*                    Input : 1 gabsignal                                   */
/*                    Output: 0 gabsignal                                   */
/*                                                                          */
/****************************************************************************/
#include "mpp.h"
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>

/* Fix 7: proper prototypes at file scope replacing K&R in-body declarations */
int    sig_copy(GABSIGNAL input, GABSIGNAL output);
void   delete_gabsignal(GABSIGNAL gabsignal);
GABSIGNAL new_struct_gabsignal(void);

/* Fix 8: removed redundant 'extern' — sig_mean and sig_variance are defined
 * in this file, not externally. Forward declarations suffice. */
double sig_mean(GABSIGNAL input);
double sig_variance(GABSIGNAL input);

/*
 * sig_zero — set all values in a gabsignal to zero.
 *
 * Fix 9: manual loop replaced with memset. IEEE 754 guarantees all-zero
 *        bits represents 0.0, and memset uses SIMD internally.
 */
int sig_zero(GABSIGNAL gabsignal)
{
    memset(gabsignal->values, 0,
           (size_t)gabsignal->size * sizeof(double));
    return 0;
}

/*
 * sig_abs — replace each value with its absolute value.
 */
int sig_abs(GABSIGNAL gabsignal)
{
    int j;
    for (j = 0; j < gabsignal->size; j++)
        gabsignal->values[j] = fabs(gabsignal->values[j]);
    return 0;
}

/*
 * sig_log — replace each value with log(|value|).
 */
int sig_log(GABSIGNAL gabsignal)
{
    int j;
    for (j = 0; j < gabsignal->size; j++)
        gabsignal->values[j] = log(fabs(gabsignal->values[j]));
    return 0;
}

/*
 * sig_exp — replace each value with exp(value).
 *
 * Fix 3/4: perror() → fprintf+return. The second guard was also unsafe —
 *          if gabsignal was NULL the second check would crash immediately.
 *          Now both guards return cleanly before any dereference.
 */
int sig_exp(GABSIGNAL gabsignal)
{
    int j;

    if (gabsignal == (GABSIGNAL)NULL) {
        fprintf(stderr, "sig_exp(): null input!\n");
        return -1;
    }
    if (gabsignal->values == (double *)NULL) {
        fprintf(stderr, "sig_exp(): null values pointer!\n");
        return -1;
    }
    for (j = 0; j < gabsignal->size; j++)
        gabsignal->values[j] = exp(gabsignal->values[j]);
    return 0;
}

/*
 * sig_square — replace each value with its square.
 */
int sig_square(GABSIGNAL gabsignal)
{
    int j;
    for (j = 0; j < gabsignal->size; j++) {
        double v = gabsignal->values[j];
        gabsignal->values[j] = v * v;
    }
    return 0;
}

/*
 * sig_add_num — add a scalar to every value.
 */
void sig_add_num(GABSIGNAL gabsignal, double num)
{
    int j;
    for (j = 0; j < gabsignal->size; j++)
        gabsignal->values[j] += num;
}

/*
 * sig_sub_num — subtract a scalar from every value.
 * Fix 6: comment was wrong ("Add a number") — corrected.
 */
int sig_sub_num(GABSIGNAL gabsignal, double num)
{
    int j;
    for (j = 0; j < gabsignal->size; j++)
        gabsignal->values[j] -= num;
    return 0;
}

/*
 * sig_mult_num — multiply every value by a scalar.
 * Fix 6: comment was wrong ("Add a number") — corrected.
 */
int sig_mult_num(GABSIGNAL gabsignal, double num)
{
    int j;
    for (j = 0; j < gabsignal->size; j++)
        gabsignal->values[j] *= num;
    return 0;
}

/*
 * sig_div_num — divide every value by a scalar.
 * Fix 6: comment was wrong ("Add a number") — corrected.
 */
void sig_div_num(GABSIGNAL gabsignal, double num)
{
    int j;
    if (!num)
        return;
    for (j = 0; j < gabsignal->size; j++)
        gabsignal->values[j] /= num;
}

/*
 * sig_pos_min_max — find the positions of the min and max values.
 *
 * Outputs:
 *   pmin   position (index) of the minimum value (double *)
 *   pmax   position (index) of the maximum value (double *)
 *
 * Fix 5: magic numbers 999999/-999999 replaced with DBL_MAX/-DBL_MAX
 *        from <float.h>. The old values would give wrong results for
 *        signals whose values fall entirely outside ±999999.
 */
int sig_pos_min_max(GABSIGNAL gabsignal, double *pmin, double *pmax)
{
    int j;
    double vmin, vmax;

    vmin = DBL_MAX;
    vmax = -DBL_MAX;
    for (j = 0; j < gabsignal->size; j++) {
        if (vmax < gabsignal->values[j]) {
            vmax  = gabsignal->values[j];
            *pmax = (double)j;
        }
        if (vmin > gabsignal->values[j]) {
            vmin  = gabsignal->values[j];
            *pmin = (double)j;
        }
    }
    return 0;
}

/*
 * sig_min_max — find the minimum and maximum values.
 *
 * Outputs:
 *   vmin   minimum value found (double *)
 *   vmax   maximum value found (double *)
 */
int sig_min_max(GABSIGNAL gabsignal, double *vmin, double *vmax)
{
    int j;

    *vmin = MAX_VALUE;
    *vmax = MIN_VALUE;
    for (j = 0; j < gabsignal->size; j++) {
        if (*vmax < gabsignal->values[j]) *vmax = gabsignal->values[j];
        if (*vmin > gabsignal->values[j]) *vmin = gabsignal->values[j];
    }
    return 0;
}

/*
 * sig_mean — compute the arithmetic mean of a gabsignal.
 */
double sig_mean(GABSIGNAL input)
{
    int j;
    double sum = 0.0;
    for (j = 0; j < input->size; j++)
        sum += input->values[j];
    return sum / (double)input->size;
}

/*
 * sig_variance — compute the variance of a gabsignal.
 *
 * Fix 10: eliminated the allocation+copy+square pass. A single-pass
 *         algorithm computes var = E[X^2] - E[X]^2 inline with no
 *         temporary GABSIGNAL, which is also numerically more stable.
 * Fix 1:  the original new_struct_gabsignal() was not NULL-checked,
 *         and the allocation is now gone entirely.
 */
double sig_variance(GABSIGNAL input)
{
    int j;
    double mean, var = 0.0;

    mean = sig_mean(input);
    for (j = 0; j < input->size; j++) {
        double d = input->values[j] - mean;
        var += d * d;
    }
    return var / (double)input->size;
}

/*
 * integral — compute the sum of all values (rectangular integration).
 */
double integral(GABSIGNAL gabsignal)
{
    int j;
    double sum = 0.0;
    for (j = 0; j < gabsignal->size; j++)
        sum += gabsignal->values[j];
    return sum;
}

/*
 * sig_L2_norm_sq — compute the squared L2 norm of a gabsignal.
 *
 * Fix 10: eliminated the allocation+copy+square pass. Computed directly
 *         in a single pass with no temporary GABSIGNAL.
 * Fix 2:  the original new_struct_gabsignal() was not NULL-checked,
 *         and the allocation is now gone entirely.
 */
double sig_L2_norm_sq(GABSIGNAL input)
{
    int j;
    double norm = 0.0;
    for (j = 0; j < input->size; j++)
        norm += input->values[j] * input->values[j];
    return norm;
}

/*
 * farray_L2_sq_norm — compute the squared L2 norm of a raw double array.
 */
double farray_L2_sq_norm(double *value, int size)
{
    int i;
    double norm = 0.0;
    double *v = value;
    for (i = 0; i < size; i++, v++)
        norm += (*v) * (*v);
    return norm;
}

/*
 * farray_L2_norm — compute the L2 norm of a raw double array.
 *
 * Fix 11: now calls farray_L2_sq_norm() to eliminate the duplicate loop.
 */
double farray_L2_norm(double *value, int size)
{
    return sqrt(farray_L2_sq_norm(value, size));
}

/*
 * SigAbsMx — find the maximum absolute value in a gabsignal and its position.
 *
 * Outputs:
 *   MxValue    the sample with the largest absolute value (double *)
 *   position   its index (double *)
 *
 * Fix 3: perror() → fprintf+return so execution stops on NULL input
 *        instead of crashing on gabsignal->size below.
 */
void SigAbsMx(GABSIGNAL gabsignal, double *MxValue, double *position)
{
    int i;

    /* Fix 3: perror() → fprintf+return */
    if (gabsignal == (GABSIGNAL)NULL) {
        fprintf(stderr, "SigAbsMx(): null argument!\n");
        return;
    }

    *MxValue  = 0.0;
    *position = 0.0;
    for (i = 0; i < gabsignal->size; i++) {
        if (fabs(gabsignal->values[i]) > fabs(*MxValue)) {
            *MxValue  = gabsignal->values[i];
            *position = (double)i;
        }
    }
}

/*
 * end of sig_fct10.c
 */