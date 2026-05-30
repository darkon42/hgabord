/*
 * gabord.h — public API for the Gabor Matching Pursuit library
 *
 * Include this header in any project that links against libgabord.a.
 * Set the stopping-criteria globals before each call to gabord().
 */

#ifndef GABORD_H
#define GABORD_H

/*
 * gabord() — Gabor Matching Pursuit decomposition.
 *
 * data:     input signal (double[], SigSize samples)
 * SigSize:  number of samples; rounded up to the next power of 2 internally
 * booksize: OUTPUT — number of atoms decomposed
 *
 * Returns a heap-allocated double[5 * *booksize].  Each atom occupies five
 * consecutive doubles: octave, freq_id, position, coefficient, phase.
 * The caller must free() the returned pointer.  Returns NULL on error.
 *
 * Stopping criteria globals must be set before calling gabord().
 * If none are set, the library defaults apply (citer=1, max_num_iter=1000).
 */
double *gabord(double *data, int SigSize, int *booksize);

/* --- Stopping criteria (OR-ed; decomposition stops on the first met) --- */

extern int    cpct;         /* 1 = stop when book energy >= th_pct% of signal */
extern double th_pct;       /* energy percentage threshold                     */

extern int    citer;        /* 1 = stop after max_num_iter atoms               */
extern int    max_num_iter; /* maximum atom count (default 500)                */

extern int    cth;          /* 1 = stop when last atom energy < thatomnrj      */
extern double thatomnrj;    /* atom energy threshold                           */

extern int    ccoh;         /* 1 = coherence criterion (SOT==SOF==1 only)      */

/* Release all library-internal heap allocations.
 * Call once at program exit after the last gabord() call.
 * FFTW plans are freed automatically via atexit() and need not be mentioned. */
void gabord_cleanup(void);

/* Explicitly zero all per-thread state (threadprivate variables).
 * Call once per worker thread inside a parallel region, BEFORE the first
 * gabord() call, when running with --jobs N > 1.  Without this, new OpenMP
 * threads have undefined initial values for threadprivate variables. */
void gabord_reset(void);

#endif /* GABORD_H */
