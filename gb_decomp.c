#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mpp.h"
#include "cwt1d.h"

/*
 * decomposition for the gabor library
 *
 * Migrated from custom srfft (fft.c) to FFTW3.
 *
 * Key layout note:
 *   GABSIGNAL->values is a split-complex array of total size DoubleSigSize:
 *     values[0 .. SigSize-1]            = real part
 *     values[SigSize .. DoubleSigSize-1] = imaginary part
 *
 *   FFTW operates on interleaved fftw_complex (re,im pairs), so we
 *   pack/unpack around every transform call using the helpers below.
 *   Plans are created once per unique size and reused for performance.
 *
 * Inputs:
 *   trans                Stores the result of the transformation (GABSIGNAL *)
 *   gabsignal            Input (complex) gabsignal to be decomposed (GABSIGNAL)
 *   filter               Stored basic gabor functions (FILTER *)
 *   SubsampleOctaveTime  Octave number at which to begin subsampling in time
 *   SubsampleOctaveFreq  Octave number at which to begin subsampling in freq
 *   MinOctave            Minimum octave index
 *   MaxOctave            Number of octaves
 *   ShiftOctave          Octave at which to switch to the second formula
 *
 * Outputs:
 *   trans                Updated in-place and returned
 *
 * Bugs (inherited):
 *   1. Handles complex gabsignal only.
 *   2. Does not validate SubsampleOctaveTime/Freq or ShiftOctave ranges.
 *   3. Requires gabsignal->size to be a power of 2.
 */

/* -----------------------------------------------------------------------
 * Plan cache: one forward + one inverse plan per signal size.
 * FFTW plans are expensive to create; reuse them across calls.
 * ----------------------------------------------------------------------- */
typedef struct {
    int          size;       /* SigSize this plan was built for */
    fftw_plan    forward;
    fftw_plan    inverse;
    fftw_complex *buf;       /* scratch buffer owned by this entry */
} PlanEntry;

#define MAX_PLAN_CACHE 32
static PlanEntry plan_cache[MAX_PLAN_CACHE];
static int       plan_cache_count = 0;

/* Call this at program exit (or via atexit()) to release all FFTW resources. */
void GaborFFTCleanup(void)
{
    int i;
    for (i = 0; i < plan_cache_count; i++) {
        fftw_destroy_plan(plan_cache[i].forward);
        fftw_destroy_plan(plan_cache[i].inverse);
        fftw_free(plan_cache[i].buf);
    }
    plan_cache_count = 0;
    fftw_cleanup();
}

/*
 * Look up or create FFTW plans for a given SigSize.
 * Returns a pointer to the PlanEntry (never NULL — exits on alloc failure).
 */
static PlanEntry *get_plans(int SigSize)
{
    int i;
    PlanEntry *e;

	// At the top of get_plans(), first time only:
	static int cleanup_registered = 0;
	if (!cleanup_registered) {
		atexit(GaborFFTCleanup);
		cleanup_registered = 1;
	}

    for (i = 0; i < plan_cache_count; i++)
        if (plan_cache[i].size == SigSize)
            return &plan_cache[i];

    if (plan_cache_count >= MAX_PLAN_CACHE) {
        fprintf(stderr, "GaborDecomp: plan cache full (MAX_PLAN_CACHE=%d)\n",
                MAX_PLAN_CACHE);
        exit(1);
    }

    e = &plan_cache[plan_cache_count];
    e->size = SigSize;
    e->buf  = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * SigSize);
    if (!e->buf) { perror("fftw_malloc"); exit(1); }

    /*
     * FFTW_MEASURE lets FFTW benchmark a few strategies at plan-creation time
     * and pick the fastest for this machine. Use FFTW_ESTIMATE if you want
     * faster startup at the cost of slightly slower transforms.
     */
    e->forward = fftw_plan_dft_1d(SigSize, e->buf, e->buf,
                                  FFTW_FORWARD,  FFTW_MEASURE);
    e->inverse = fftw_plan_dft_1d(SigSize, e->buf, e->buf,
                                  FFTW_BACKWARD, FFTW_MEASURE);
    //if (!e->forward || !e->inverse) {
    //    fprintf(stderr, "GaborDecomp: fftw_plan_dft_1d failed for N=%d\n", SigSize);
    //    exit(1);
    //}
	
	if (!e->forward || !e->inverse) {
		fprintf(stderr, "GaborDecomp: fftw_plan_dft_1d failed for N=%d\n", SigSize);
		if (e->forward) fftw_destroy_plan(e->forward);
		if (e->inverse) fftw_destroy_plan(e->inverse);
		fftw_free(e->buf);
    exit(1);
}



    plan_cache_count++;
    return e;
}

/* -----------------------------------------------------------------------
 * Pack split-complex GABSIGNAL layout into interleaved fftw_complex buffer.
 *   xr[0..N-1] = real part
 *   xi[0..N-1] = imaginary part  (xi = xr + N in the gabsignal)
 * ----------------------------------------------------------------------- */
static void pack_interleaved(const double *xr, const double *xi,
                             fftw_complex *out, int N)
{
    int i;
    for (i = 0; i < N; i++) {
        out[i][0] = xr[i];
        out[i][1] = xi[i];
    }
}

/* Unpack interleaved fftw_complex back to split-complex layout. */
static void unpack_interleaved(const fftw_complex *in,
                               double *xr, double *xi, int N)
{
    int i;
    for (i = 0; i < N; i++) {
        xr[i] = in[i][0];
        xi[i] = in[i][1];
    }
}

/* -----------------------------------------------------------------------
 * gabor_fft: forward DFT of a GABSIGNAL (split-complex, in-place).
 *   The transform is unnormalised, matching the old FFT() behaviour.
 * ----------------------------------------------------------------------- */
static void gabor_fft(GABSIGNAL sig)
{
    int SigSize = sig->size >> 1;
    double *xr   = sig->values;
    double *xi   = sig->values + SigSize;
    PlanEntry *e = get_plans(SigSize);

    pack_interleaved(xr, xi, e->buf, SigSize);
    fftw_execute(e->forward);
    unpack_interleaved(e->buf, xr, xi, SigSize);
}

/* -----------------------------------------------------------------------
 * gabor_ifft: inverse DFT of a GABSIGNAL (split-complex, in-place).
 *   FFTW's backward transform is unnormalised (output scaled by N), so we
 *   divide by N afterwards — identical to the old srifft() behaviour.
 * ----------------------------------------------------------------------- */
static void gabor_ifft(GABSIGNAL sig)
{
    int i;
    int SigSize = sig->size >> 1;
    double *xr   = sig->values;
    double *xi   = sig->values + SigSize;
    double fac   = 1.0 / SigSize;
    PlanEntry *e = get_plans(SigSize);

    pack_interleaved(xr, xi, e->buf, SigSize);
    fftw_execute(e->inverse);
    unpack_interleaved(e->buf, xr, xi, SigSize);

    for (i = 0; i < SigSize; i++) {
        xr[i] *= fac;
        xi[i] *= fac;
    }
}

/* -----------------------------------------------------------------------
 * GaborSubsample — unchanged from original
 * ----------------------------------------------------------------------- */
void GaborSubsample(GABSIGNAL gabsignal)
{
    int newsize, new_sig_size, i;
    double *value1, *value2, *value3, *value4;

    if (gabsignal == (GABSIGNAL)NULL)
        perror("GaborSubsample(): null input!");
		
	if (gabsignal == (GABSIGNAL)NULL) {
		fprintf(stderr, "GaborSubsample(): null input\n");
		return;
	}
    //if (gabsignal->size <= 0) return;

    newsize      = gabsignal->size >> 1;
    new_sig_size = newsize >> 1;

    gabsignal->size = newsize;

    value1 = gabsignal->values;
    value2 = gabsignal->values + new_sig_size;
    value3 = gabsignal->values + newsize;
    value4 = value3 + new_sig_size;

    for (i = 0; i < new_sig_size; i++) {
        *value1 += *value2;
        *value1++ /= 2.0;
        *value2++ = ((*value3++) + (*value4++)) / 2.0;
    }
}

/* -----------------------------------------------------------------------
 * GaborDecomp — same algorithm as before, FFT/IFFT calls replaced.
 * ----------------------------------------------------------------------- */
GABSIGNAL *GaborDecomp(GABSIGNAL *trans, GABSIGNAL gabsignal, FILTER *filter,
                       int SubsampleOctaveTime,
                       int SubsampleOctaveFreq,
                       int MinOctave,
                       int MaxOctave,
                       int ShiftOctave)
{
    int i, j, k, k_length, SigSize, DoubleSigSize, shift, size;
    int SampleRate_t = 1;
    int HalfSigSize, index;
    int LnSigSize;
    int n, n_length, n_2_length, SampleRate_f;
    GABSIGNAL sigtmp = (GABSIGNAL)NULL;
    double *value_r, *value_i, *value1, *value2, *value3, *value4, *value5;
    double *sig_value_r = NULL;
    double *sig_value_i = NULL;
    double SqrtSigSize;

    /* Forward-declare functions defined elsewhere */
    double *darray_malloc();
    GABSIGNAL new_gabsignal(), *GaborAllocFilter();
    int delete_gabsignal();
    void change_gabsignal(), sig_range_copy(),
         farray_translate(double *, double *, int, int);

    if (trans == (GABSIGNAL *)NULL)
        trans = GaborAllocFilter(MaxOctave - MinOctave + 3);
    
	//if (gabsignal == (GABSIGNAL)NULL || gabsignal->size == 0)
     //   perror("Empty gabsignal!");
		
	if (gabsignal == (GABSIGNAL)NULL || gabsignal->size == 0) {
		fprintf(stderr, "GaborDecomp: empty gabsignal\n");
		return NULL; }

    DoubleSigSize = gabsignal->size;
    SigSize       = DoubleSigSize >> 1;
    HalfSigSize   = SigSize >> 1;
    SqrtSigSize   = (double)sqrt((double)SigSize);
    LnSigSize     = (int)log2((double)SigSize);

    /* Store original signal in trans[0] */
    trans[0] = gabsignal;

    /*
     * Compute the Fourier transform of the original signal and store it
     * in trans[MaxOctave - MinOctave + 2].
     */
    if (trans[MaxOctave - MinOctave + 2] == (GABSIGNAL)NULL)
        trans[MaxOctave - MinOctave + 2] = new_gabsignal(DoubleSigSize);
    else if (trans[MaxOctave - MinOctave + 2]->size != DoubleSigSize)
        change_gabsignal(trans[MaxOctave - MinOctave + 1], DoubleSigSize);

    value1 = trans[MaxOctave - MinOctave + 2]->values;
    value2 = trans[MaxOctave - MinOctave + 2]->values + SigSize;
    value3 = trans[0]->values;
    value4 = trans[0]->values + SigSize;
    for (i = 0; i < SigSize; i++) {
        *value1++ = *value3++;
        *value2++ = *value4++;
    }

    /* --- replaced: FFT(trans[...], &trans[...]) --- */
    gabor_fft(trans[MaxOctave - MinOctave + 2]);

    /*
     * Calculate <f, g_j,k,n> using the first formula (j < ShiftOctave).
     */
    if (sigtmp == (GABSIGNAL)NULL)
        sigtmp = new_gabsignal(DoubleSigSize);
    sig_value_r = darray_malloc(SigSize);
    sig_value_i = darray_malloc(SigSize);

    for (j = MinOctave; j < ShiftOctave; j++) {

        /* k dimension */
        if (j > (LnSigSize - SubsampleOctaveFreq)) {
            k_length    = 1 << (LnSigSize - 1);
            SampleRate_f = 1;
        } else {
            k_length    = 1 << (SubsampleOctaveFreq + j - 2);
            SampleRate_f = 1 << (LnSigSize - SubsampleOctaveFreq - j + 1);
        }

        /* n dimension */
        if (j < SubsampleOctaveTime) {
            SampleRate_t = 1;
            n_length = SigSize;
        } else {
            SampleRate_t = 1 << (j - SubsampleOctaveTime + 1);
            n_length = SigSize >> (j - SubsampleOctaveTime + 1);
        }
        n_2_length = n_length << 1;

        size = k_length * n_2_length + n_2_length + 2;   /* PJF */

        if (trans[j - MinOctave + 1] == (GABSIGNAL)NULL)
            trans[j - MinOctave + 1] = new_gabsignal(size);
        else if (trans[j - MinOctave + 1]->size != size)
            change_gabsignal(trans[j - MinOctave + 1], size);

        value_r = trans[j - MinOctave + 1]->values;

        for (k = 0; k <= k_length; k++) {

            value1 = trans[MaxOctave - MinOctave + 2]->values;
            value2 = trans[MaxOctave - MinOctave + 2]->values + SigSize;

            /* Shift: f^(w + 2^(L-L1-j+1)*k) */
            shift = -SampleRate_f * k;
            farray_translate(value1, sig_value_r, SigSize, shift);
            farray_translate(value2, sig_value_i, SigSize, shift);

            /* Multiply by filter in frequency domain */
            value1 = sigtmp->values;
            value2 = sigtmp->values + SigSize;
            value3 = filter[LnSigSize - j]->values;
            value4 = sig_value_r;
            value5 = sig_value_i;

            value1[0] = value3[0] * value4[0] * SqrtSigSize;
            value2[0] = value3[0] * value5[0] * SqrtSigSize;

            for (i = 1; i < HalfSigSize; i++) {
                if (i < filter[LnSigSize - j]->size) {
                    value1[i]          = value3[i] * value4[i]          * SqrtSigSize;
                    value2[i]          = value3[i] * value5[i]          * SqrtSigSize;
                    value1[SigSize-i]  = value3[i] * value4[SigSize-i]  * SqrtSigSize;
                    value2[SigSize-i]  = value3[i] * value5[SigSize-i]  * SqrtSigSize;
                } else {
                    value1[i] = value2[i] = 0.0;
                    value1[SigSize-i] = value2[SigSize-i] = 0.0;
                }
            }

            if (filter[LnSigSize - j]->size > HalfSigSize) {
                value1[HalfSigSize] = value3[HalfSigSize] * value4[HalfSigSize] * SqrtSigSize;
                value2[HalfSigSize] = value3[HalfSigSize] * value5[HalfSigSize] * SqrtSigSize;
            } else {
                value1[HalfSigSize] = value2[HalfSigSize] = 0.0;
            }

            /* Subsample in time domain if needed */
            if (j >= SubsampleOctaveTime)
                for (i = SubsampleOctaveTime; i <= j; i++)
                    GaborSubsample(sigtmp);

            /* --- replaced: IFFT(sigtmp, &sigtmp) --- */
            gabor_ifft(sigtmp);

            value1  = sigtmp->values;
            value2  = sigtmp->values + n_length;
            value_i = value_r + n_length;
            for (i = 0; i < n_length; i++) {
                *value_r++ = *value1++;
                *value_i++ = *value2++;
            }

            /* Restore size modified by GaborSubsample */
            sigtmp->size = DoubleSigSize;
            value_r = value_i;

        } /* end k loop */
    } /* end first j loop */

    /*
     * Calculate <f, g_j,k,n> using the second formula (j >= ShiftOctave).
     */
    for (j = ShiftOctave; j <= MaxOctave; j++) {

        /* k dimension */
        if (j > (LnSigSize - SubsampleOctaveFreq)) {
            k_length    = (1 << LnSigSize) - 1;
            SampleRate_f = 1;
        } else {
            k_length    = 1 << (SubsampleOctaveFreq + j - 2);
            SampleRate_f = 1 << (LnSigSize - SubsampleOctaveFreq - j + 1);
        }

        /* n dimension */
        if (j < SubsampleOctaveTime) {
            SampleRate_t = 1;
            n_length = SigSize;
        } else {
            SampleRate_t = 1 << (j - SubsampleOctaveTime + 1);
            n_length = SigSize >> (j - SubsampleOctaveTime + 1);
        }
        n_2_length = n_length << 1;

        size = k_length * n_2_length + n_2_length + 2;  /* PJF */

        if (trans[j - MinOctave + 1] == (GABSIGNAL)NULL)
            trans[j - MinOctave + 1] = new_gabsignal(size);
        else if (trans[j - MinOctave + 1]->size != size)
            change_gabsignal(trans[j - MinOctave + 1], size);

        for (n = 0; n < n_length; n++) {

            shift = n * SampleRate_t;

            /*
             * Multiply signal by shifted filter window:
             *   sigtmp[i] = signal[i] * filter_j[(i - shift) mod SigSize]
             * The filter is real and symmetric so we index into its stored
             * half via the same mirror logic as the original.
             */
            value1 = sigtmp->values;
            value2 = value1 + SigSize;
            value3 = trans[0]->values;
            value4 = value3 + SigSize;
            value5 = filter[j]->values;

            for (i = 0; i < SigSize; i++) {
                index = (SigSize + i - shift) % SigSize;
                if (index < filter[j]->size && index != 0
                        && index != HalfSigSize) {
                    value1[i] = value3[i] * value5[index];
                    value2[i] = value4[i] * value5[index];
                } else if (SigSize - index < filter[j]->size
                        && index != 0 && index != HalfSigSize) {
                    value1[i] = value3[i] * value5[SigSize - index];
                    value2[i] = value4[i] * value5[SigSize - index];
                } else if (index == 0) {
                    value1[i] = value3[i] * value5[0];
                    value2[i] = value4[i] * value5[0];
                } else if (index == HalfSigSize) {
                    if (filter[j]->size > HalfSigSize) {
                        value1[i] = value3[i] * value5[HalfSigSize];
                        value2[i] = value4[i] * value5[HalfSigSize];
                    } else {
                        value1[i] = value2[i] = 0.0;
                    }
                } else {
                    value1[i] = value2[i] = 0.0;
                }
            }

            /* Subsample in frequency domain */
            for (i = j; i <= LnSigSize - SubsampleOctaveFreq; i++)
                GaborSubsample(sigtmp);

            /* --- replaced: FFT(sigtmp, &sigtmp) --- */
            gabor_fft(sigtmp);

            value1  = sigtmp->values;
            value2  = value1 + (sigtmp->size >> 1);
            value_r = trans[j - MinOctave + 1]->values + n;

            for (k = 0; k <= k_length; k++) {
                value_i  = value_r + n_length;
                *value_r = (*value1++) * (double)SampleRate_f;
                value_r += n_2_length;
                *value_i = (*value2++) * (double)SampleRate_f;
            }

            sigtmp->size = DoubleSigSize;

        } /* end n loop */
    } /* end second j loop */

    /* Normalise the stored Fourier transform of the original signal */
    value1 = trans[MaxOctave - MinOctave + 2]->values;
    value2 = value1 + SigSize;
    for (i = 0; i < SigSize; i++) {
        *value1++ /= SqrtSigSize;
        *value2++ /= SqrtSigSize;
    }

    /* Free local allocations */
    free((char *)sig_value_r);
    free((char *)sig_value_i);
    delete_gabsignal(sigtmp);

	/* Caller is responsible for freeing trans[1..MaxOctave-MinOctave+2]
	* and trans[MaxOctave-MinOctave+2] via delete_gabsignal(). */
	
    return trans;
} /* end of GaborDecomp() */