#include "mpp.h"
#include "cwt1d.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define BOUNDN(j,epslon,sigma) ((double)(1<<(j))*(sigma)*sqrt(-2.0*log((epslon))))
#define INDEXF(j,L,n) (abs(n)<<((L)-(j)-1))

#define M_PI       3.14159265358979323846

extern double epslonG;

/* Fix 6: single extern declaration at file scope is sufficient —
 * the redundant local declaration inside GaborFreeFilter is removed. */
extern void clear_gabsignal(GABSIGNAL gabsignal);

/* Fix 5: proper prototypes replacing K&R-style in-body declarations */
GABSIGNAL   new_gabsignal(int size);
GABSIGNAL   FreeSignal(GABSIGNAL gabsignal);
FILTER     *GaborAllocFilter(int MaxOctave);
void        change_gabsignal(GABSIGNAL gabsignal, int size);
void        CreateExponential(double position, double min, double max,
                               int size, GABSIGNAL *psigout);
GABSIGNAL   Gaussian2Signal(double min_t, double max_t, double size_g,
                              int size_s, double sigma);

/*
 * GaborAllocFilter — allocate an array of MaxOctave+1 FILTER slots,
 * all initialised to NULL.
 *
 * Inputs:
 *	MaxOctave	number of octaves (int)
 *
 * Returns a newly allocated FILTER array, or NULL on failure.
 * Caller is responsible for freeing via GaborFreeFilter().
 *
 * Fix 4: perror() replaced with fprintf+return NULL so the caller
 *        receives NULL instead of a crash on the next dereference.
 */
FILTER *GaborAllocFilter(int MaxOctave)
{
    FILTER *pfilter;
    int j;

    pfilter = (GABSIGNAL *)malloc((unsigned)((MaxOctave+1) *
                                              sizeof(struct gabsignal)));
    /* Fix 4: perror() → fprintf+return NULL */
    if (pfilter == (GABSIGNAL *)NULL) {
        fprintf(stderr,
                "GaborAllocFilter(): malloc failed for %d octaves\n",
                MaxOctave);
        return (FILTER *)NULL;
    }

    for (j = 0; j <= MaxOctave; j++)
        pfilter[j] = (FILTER)NULL;

    return pfilter;
}

/*
 * GaborFreeFilter — free all filter slots and the filter array itself.
 *
 * Inputs:
 *	pfilter		filter array to free (FILTER *)
 *	NumFilter	number of slots to free (int)
 *
 * Returns NULL so the caller can write:
 *   pfilter = GaborFreeFilter(pfilter, n);
 * to both free and NULL the caller's pointer in one step.
 */
FILTER *GaborFreeFilter(FILTER *pfilter, int NumFilter)
{
    int j;
    /* Fix 6: redundant local declaration of clear_gabsignal removed —
     * the extern at file scope is sufficient. */

    if (pfilter == (FILTER *)NULL)
        return (FILTER *)NULL;

    if (NumFilter > 0) {
        for (j = 0; j < NumFilter; j++) {
            if (pfilter[j] != (GABSIGNAL)NULL) {
                clear_gabsignal(pfilter[j]);
                free((char *)pfilter[j]);
            }
        }
    }
    free((char *)pfilter);

    /* Fix 2: note — this sets the LOCAL pointer to NULL. The caller must
     * use the return value to NULL their own pointer:
     *   pfilter = GaborFreeFilter(pfilter, n); */
    return (FILTER *)NULL;
}

/*
 * GaborBuildFilter — build the Gabor filter bank for a given decomposition.
 *
 * Inputs:
 *	MaxOctave	number of octaves (int)
 *	SigSize		gabsignal size — must be a power of 2 (int)
 *	sigma		deviation of the Gaussian window (double)
 *	fNorm		receives a newly allocated normalisation array (double **)
 *
 * Returns the built filter bank, or NULL on failure.
 * Caller is responsible for freeing pfilter via GaborFreeFilter() and
 * freeing *fNorm with free().
 *
 * Fix 1: filterL is freed on all error paths so it is never leaked.
 * Fix 4: perror() → fprintf+return NULL throughout.
 * Fix 7: *fNorm allocated with MaxOctave elements (indices 1..MaxOctave-1
 *        are used; index 0 is unused but kept for consistent indexing).
 */
FILTER *GaborBuildFilter(int MaxOctave, int SigSize, double sigma,
                          double **fNorm)
{
    FILTER *pfilter = (FILTER *)NULL;
    FILTER filterL  = (FILTER)NULL;   /* Fix 1: initialise to NULL for safe cleanup */
    double boundN, boundR, boundL, B1;
    int boundI, multN, index;
    int i, j, n;
    double *value;
    double norm, factor;

    /* Fix 5: all prototypes moved to file scope — no local declarations here */

    /* Allocate filter array */
    pfilter = GaborAllocFilter(MaxOctave);
    if (pfilter == (FILTER *)NULL) {
        fprintf(stderr, "GaborBuildFilter(): GaborAllocFilter failed\n");
        return (FILTER *)NULL;
    }

    /* Fix 7: allocate MaxOctave elements (index 0 unused, 1..MaxOctave-1 written) */
    *fNorm = (double *)malloc(sizeof(double) * MaxOctave);
    /* Fix 4: perror() → fprintf+goto cleanup */
    if (*fNorm == (double *)NULL) {
        fprintf(stderr, "GaborBuildFilter(): malloc for fNorm failed\n");
        goto cleanup;
    }

    /*
     * Calculate the non-periodic Gaussian over the finest grid.
     */
    boundN  = BOUNDN(MaxOctave-1, epslonG, sigma);
    filterL = Gaussian2Signal(0.0, boundN, boundN,
                               (int)boundN+1,
                               (double)(1<<(MaxOctave-1)) * sigma);
    if (filterL == (FILTER)NULL) {
        fprintf(stderr, "GaborBuildFilter(): Gaussian2Signal failed\n");
        goto cleanup;
    }

    /*
     * Calculate the basic Gabor functions and store in the filter bank.
     */
    for (j = 1; j < MaxOctave; j++) {
        boundN = BOUNDN(j, epslonG, sigma);
        boundI = (int)(boundN / (double)SigSize + 1.5);
        B1     = MIN(boundN, (double)(1<<(MaxOctave-1)));

        pfilter[j] = new_gabsignal((int)B1+1);
        if (pfilter[j] == (GABSIGNAL)NULL) {
            fprintf(stderr,
                    "GaborBuildFilter(): new_gabsignal failed at octave %d\n",
                    j);
            /* Fix 1: filterL is freed via cleanup label — not leaked */
            goto cleanup;
        }

        value  = pfilter[j]->values;
        factor = sqrt((double)(1<<(MaxOctave-j-1)));

        for (i = 0; i <= boundI; i++) {
            multN  = i<<MaxOctave;
            boundR = MIN(boundN - (double)multN, B1);
            for (n = 0; n <= (int)boundR; n++) {
                index = INDEXF(j, MaxOctave, n+multN);
                /* Fix 8: indented consistently with surrounding code */
                if (index >= filterL->size)
                    continue;
                value[n] += filterL->values[index] * factor;
            }
            if (i == 0)
                continue;
            boundR = MIN(boundN + (double)multN, B1);
            boundL = MAX(0, (double)multN - boundN);
            for (n = (int)boundL; n <= (int)boundR; n++) {
                index = INDEXF(j, MaxOctave, n-multN);
                /* Fix 8: indented consistently */
                if (index >= filterL->size)
                    continue;
                value[n] += filterL->values[index] * factor;
            }
        }

        /* Normalise this octave's filter */
        B1   = pfilter[j]->size - 1;
        norm = 0.0;
        for (n = 0; n <= (int)B1; n++)
            norm += value[n] * value[n];
        norm = 2.0*norm - value[0]*value[0];
        if ((int)B1 == (1<<(MaxOctave-1)))
            norm -= value[(int)B1] * value[(int)B1];
        norm = sqrt(norm);
        (*fNorm)[j] = 1.0 / norm;
        for (n = 0; n <= (int)B1; n++)
            value[n] /= norm;
    }

    /*
     * Store exp(i*2*pi*k/N) in the last filter slot.
     */
    CreateExponential((double)(PI2/SigSize), 0.0, (double)SigSize,
                      SigSize, &pfilter[MaxOctave]);

    /* Free the temporary fine-grid Gaussian */
    filterL = FreeSignal(filterL);
    return pfilter;

cleanup:
    /* Fix 1: free filterL on any error path so it is never leaked */
    if (filterL != (FILTER)NULL)
        filterL = FreeSignal(filterL);
    if (*fNorm != (double *)NULL) {
        free(*fNorm);
        *fNorm = (double *)NULL;
    }
    GaborFreeFilter(pfilter, MaxOctave);
    return (FILTER *)NULL;
}

/*
 * end of gb_filter.c
 */