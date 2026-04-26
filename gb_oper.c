#include "mpp.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/* real part of complex multiplication */
#define CMREAL(a1,b1,a2,b2)      ((a1)*(a2)-(b1)*(b2))
/* imaginary part of complex multiplication */
#define CMIMAGINARY(a1,b1,a2,b2) ((a1)*(b2)+(a2)*(b1))
/* calculate the index for the transform given by k2 and p2 */
#define TINDXREAL(k2,p2,Ljl2)    (((k2)<<((Ljl2)+1))+(p2))
#define TINDXIMAG(k2,p2,Ljl2)    (((k2)<<((Ljl2)+1))+(1<<(Ljl2))+(p2))
#define SSFACTOR(j,l)            (((j)>(l))?(l):(j))
#define STEP(j,l)                (1<<((j)-(l)))

#ifdef __unix__
#else
    #define M_SQRT2    1.41421356237309504880
#endif

/* global variables */
extern double *cur_norm;
extern double *pfG, *pfCE1;
extern double *pfCE2, *pfCE3, *pfC;
extern double *pfB;
extern int *pnAep, *pnIep;
extern int cur_l, cur_h;   /* defined in ng_cmd.c */
extern double epslonG;     /* defined in ng_cmd.c */
extern int nBoundG;        /* defined in ng_cmd.c */

/* Fix 2: global INDEX objects allocated once in GaborGetGaborNorm and
 * reused across calls. Call GaborOperCleanup() at program exit to free. */
INDEX indexG1 = (INDEX)NULL;
INDEX indexG2 = (INDEX)NULL;

/* Fix 5: proper prototypes at file scope replacing all K&R-style
 * in-body declarations. */
WORD      AllocWord(void);
INDEX     AllocIndex(void);
GABSIGNAL new_gabsignal(int size);
void      SigAbsMx(GABSIGNAL sig, double *maxVal, double *maxPos);
void      complex_array_max(double *value, int size, int size_op,
                             double *modula, double *v_r, double *v_i,
                             int *index);
double    GaborGetCoeff(double *pfGNorm, int j1, int j2, double *pfC);
void      GaborGetInnerProd(INDEX index1, INDEX index2, int nL,
                             int nST, int nSF, int nDeltaP, int nDeltaW,
                             int nBound, double C, int *pnAep, int *pnIep,
                             double *pfCE1, double *pfCE2, double *pfCE3,
                             double *pfG, double *pfReal, double *pfImaginary);

/* Fix 4: all functions converted from K&R to ANSI prototype style —
 * forward declarations for file-local functions */
static void   GaborArrayMax(FILTER *filter, double *value, int size,
                             long SigSize, int Log2SigSize, int octave,
                             int freq, int SampleRate_n,
                             double *modula, double *v_r, double *v_i,
                             int *index);
static double GaborGetGaborNorm(double modula, double c_r, double c_i,
                                 long SigSize, int num_octave, int octave,
                                 long freq, long translation);
static int    getNbhd(double f[3][3], double *vector, int j,
                      long *k, long *p, int l2, int h2, int L, int flag);
static double interpolation(double f[3][3], long *k, long *p,
                             double dk, double dp, double dkF, double dpF,
                             long N, int flag);
static void   innSigWvForm(GABSIGNAL gabsignal, FILTER *filter,
                            int j, long k, long p, int L,
                            double *vr, double *vi);
static double findMaxInNbdhd(double *vector, long *k2, long *p2,
                              int Ljl2, long L2);

/*
 * GaborOperCleanup — free the persistent global INDEX objects.
 * Call at program exit (or register with atexit()).
 * Fix 2: indexG1/indexG2 were allocated but never freed.
 */
void GaborOperCleanup(void)
{
    if (indexG1 != (INDEX)NULL) { free(indexG1); indexG1 = (INDEX)NULL; }
    if (indexG2 != (INDEX)NULL) { free(indexG2); indexG2 = (INDEX)NULL; }
}

/*
 * GaborGetMaxFrmTrans — find the maximum modulus coefficient in the transform.
 *
 * Inputs:
 *	trans			gabor transform array (GABSIGNAL *)
 *	filter			basic gabor functions (FILTER *)
 *	MinOctave		minimum octave index (int)
 *	MaxOctave		maximum octave index (int)
 *	num_octave		total number of octaves = log2(SigSize) (int)
 *	SubsampleOctaveTime	octave at which time subsampling begins (int)
 *	SubsampleOctaveFreq	octave at which freq subsampling begins (int)
 *
 * Returns an allocated WORD with the selected atom, or NULL on error.
 * Caller owns the returned WORD.
 *
 * Bugs:
 *	1. requires gabsignal size to be a power of 2
 *	2. does not validate MinOctave, MaxOctave, SubsampleOctave* ranges
 */
WORD GaborGetMaxFrmTrans(GABSIGNAL *trans, FILTER *filter,
                          int MinOctave, int MaxOctave, int num_octave,
                          int SubsampleOctaveTime, int SubsampleOctaveFreq)
{
    WORD word = (WORD)NULL;
    int k, j, SampleRate_f, SampleRate_n;
    long SigSize, HalfSigSize, DoubleSigSize;
    int freq, bndBadK;
    int n_length, n_2_length, kend, index;
    double modula, value_r, value_i;
    double MaxModula = 0.0, MaxOct = 0.0, MaxFreq = 0.0, MaxTranslation = 0.0;
    double MaxReal = 0.0, MaxImg = 0.0, MaxNorm = 1.0;
    double *pointer;
    int LnSigSize;

    /* Fix 3: perror() → fprintf+return NULL */
    if (trans == (GABSIGNAL *)NULL || filter == (FILTER *)NULL) {
        fprintf(stderr, "GaborGetMaxFrmTrans(): null argument!\n");
        return (WORD)NULL;
    }
    if (trans[0] == (GABSIGNAL)NULL) {
        fprintf(stderr, "GaborGetMaxFrmTrans(): trans[0] is null!\n");
        return (WORD)NULL;
    }

    DoubleSigSize = trans[0]->size;
    SigSize       = DoubleSigSize >> 1;
    HalfSigSize   = SigSize >> 1;
    LnSigSize     = (int)log2((double)SigSize);

    /* The Dirac basis: j=0 */
    if ((MinOctave == 1) && (MaxOctave == LnSigSize - 1)) {
        SigAbsMx(trans[0], &MaxReal, &MaxTranslation);
        MaxOct    = 0.0;
        MaxFreq   = 0.0;
        MaxModula = fabs(MaxReal * MaxReal) / 2.0;
        MaxImg    = 0.0;
        MaxNorm   = 1.0;
    }

    /* The Gabor basis: j = MinOctave..MaxOctave */
    for (j = MinOctave; j <= MaxOctave; j++) {
        /* Fix 3: perror() → fprintf+return NULL */
        if (trans[j - MinOctave + 1] == (GABSIGNAL)NULL) {
            fprintf(stderr,
                    "GaborGetMaxFrmTrans(): trans[%d] is null!\n",
                    j - MinOctave + 1);
            return (WORD)NULL;
        }
        if (j > LnSigSize - SubsampleOctaveFreq) {
            kend        = 1 << (LnSigSize-1);
            bndBadK     = 1 << (LnSigSize-j);
            SampleRate_f = 1;
        } else {
            kend        = 1 << (SubsampleOctaveFreq+j-2);
            bndBadK     = 1 << (SubsampleOctaveFreq-1);
            SampleRate_f = 1 << (LnSigSize-SubsampleOctaveFreq-j+1);
        }
        if (j >= SubsampleOctaveTime) {
            n_length    = SigSize >> (j-SubsampleOctaveTime+1);
            n_2_length  = n_length << 1;
            SampleRate_n = 1 << (j-SubsampleOctaveTime+1);
        } else {
            n_length    = SigSize;
            n_2_length  = DoubleSigSize;
            SampleRate_n = 1;
        }
        pointer = trans[j - MinOctave + 1]->values;
        for (k = 0; k <= kend; k++) {
            if ((k < bndBadK && k > 0) || (k > (kend-bndBadK) && k < kend)) {
                pointer += n_2_length;
                continue;
            }
            freq = k * SampleRate_f;
            GaborArrayMax(filter, pointer, n_length, SigSize, LnSigSize,
                          j, freq, SampleRate_n,
                          &modula, &value_r, &value_i, &index);
            if (freq == 0 || freq == HalfSigSize)
                modula /= 2.0;
            if (modula > MaxModula) {
                MaxOct         = (double)j;
                MaxFreq        = (double)freq;
                MaxTranslation = (double)(SampleRate_n * index);
                MaxModula      = modula;
                MaxReal        = value_r;
                MaxImg         = value_i;
            }
            pointer += n_2_length;
        }
    }

    /* The Fourier basis: j = LnSigSize */
    if ((MinOctave == 1) && (MaxOctave == LnSigSize - 1)) {
        if (trans[MaxOctave-MinOctave+2] == (GABSIGNAL)NULL) {
            fprintf(stderr,
                    "GaborGetMaxFrmTrans(): Fourier trans slot is null!\n");
            return (WORD)NULL;
        }
        pointer = trans[MaxOctave-MinOctave+2]->values;
        complex_array_max(pointer, SigSize, (SigSize>>1)+1,
                          &modula, &value_r, &value_i, &index);
        if (index == 0 || index == HalfSigSize)
            modula /= 2.0;
        if (modula > MaxModula) {
            MaxModula      = modula;
            MaxOct         = (double)LnSigSize;
            MaxFreq        = (double)index;
            MaxTranslation = 0.0;
            MaxReal        = value_r;
            MaxImg         = value_i;
        }
    }

    /* Install maximum into word */
    word = AllocWord();
    if (MaxOct == (double)LnSigSize) {
        MaxNorm = (MaxFreq == 0.0 || MaxFreq == (double)HalfSigSize)
                  ? 1.0 / sqrt((double)SigSize)
                  : M_SQRT2 / sqrt((double)SigSize);
        word->value = MaxNorm;
    } else if (MaxOct == 0.0) {
        word->value = MaxNorm;
    } else {
        MaxNorm = GaborGetGaborNorm(MaxModula, MaxReal, MaxImg, SigSize,
                                    LnSigSize, (int)MaxOct,
                                    (long)MaxFreq, (long)MaxTranslation);
        word->value = MaxNorm;
    }

    if (MaxOct == 0.0 || MaxFreq == 0.0 || MaxFreq == (double)HalfSigSize) {
        MaxModula = sqrt(2.0 * MaxModula);
    } else {
        if (MaxOct != (double)LnSigSize) {
            MaxModula = sqrt(MaxModula) * MaxNorm;
            MaxReal  *= MaxNorm;
        } else {
            MaxModula = M_SQRT2 * sqrt(MaxModula);
            MaxReal  *= M_SQRT2;
        }
    }

    word->coeff          = MaxModula;
    word->index->phase   = (MaxImg > 0)
                           ?  acos(MaxReal / MaxModula)
                           : -acos(MaxReal / MaxModula);
    word->index->id       = MaxFreq;
    word->index->octave   = MaxOct;
    word->index->position = MaxTranslation;

    return word;
}

/*
 * GaborGetResidue — subtract the selected atom from the signal residue.
 *
 * Inputs:
 *	trans		gabor transform array (GABSIGNAL *)
 *	filter		basic gabor functions (FILTER *)
 *	word		the selected atom (WORD)
 *	num_octave	log2(SigSize) (int)
 *
 * Output:
 *	trans[0]	updated residue signal (GABSIGNAL)
 *
 * Fix 3: perror() → fprintf+return.
 * Fix 6: sincos() replaces separate cos()/sin() calls.
 * Fix 9: modulo-in-loop replaced with running index counters.
 */
void GaborGetResidue(GABSIGNAL *trans, FILTER *filter,
                     WORD word, int num_octave)
{
    int i, shift, index, index1, SigSize, octave, freq;
    double cos_phi, sin_phi, *value, tmp;

    /* Fix 3: perror() → fprintf+return */
    if (trans == (GABSIGNAL *)NULL || word == (WORD)NULL) {
        fprintf(stderr, "GaborGetResidue(): null input!\n");
        return;
    }
    if (trans[0] == (GABSIGNAL)NULL) {
        fprintf(stderr, "GaborGetResidue(): trans[0] is null!\n");
        return;
    }

    SigSize = trans[0]->size >> 1;

    /* Fix 6: sincos() replaces cos()/sin() on the same argument */
    sincos((double)word->index->phase, &sin_phi, &cos_phi);

    shift  = (int)word->index->position;
    octave = (int)word->index->octave;
    freq   = (int)word->index->id;

    /* Fix 3: perror() → fprintf+return */
    if (shift > SigSize || shift < 0) {
        fprintf(stderr, "GaborGetResidue(): illegal translation %d\n", shift);
        return;
    }
    shift = SigSize - shift;

    if (octave == 0) {
        /* Dirac basis */
        word->value = 1.0;
        trans[0]->values[(int)word->index->position] = 0.0;
        return;
    }

    value = trans[0]->values;
    tmp   = word->coeff * word->value;

    if (octave == num_octave) {
        /* Fourier basis
         * Fix 9: running index replaces (freq*i)%SigSize per iteration */
        int freqStep = freq % SigSize;
        index = 0;
        for (i = 0; i < SigSize; i++) {
            *value++ -= tmp *
                (filter[num_octave]->values[index]          * cos_phi -
                 filter[num_octave]->values[index + SigSize] * sin_phi);
            index += freqStep;
            if (index >= SigSize) index -= SigSize;
        }
        return;
    }

    /* Gabor basis
     * Fix 9: running counters replace (freq*i)%SigSize and (shift+i)%SigSize */
    {
        int freqStep = freq % SigSize;
        index  = 0;   /* tracks (freq*i) % SigSize */
        index1 = shift % SigSize;  /* tracks (shift+i) % SigSize */
        for (i = 0; i < SigSize; i++) {
            double filt;
            if (index1 < filter[octave]->size)
                filt = filter[octave]->values[index1];
            else if (SigSize - index1 < filter[octave]->size)
                filt = filter[octave]->values[SigSize - index1];
            else
                filt = 0.0;

            if (filt != 0.0)
                *value -= tmp * filt *
                    (filter[num_octave]->values[index]          * cos_phi -
                     filter[num_octave]->values[SigSize + index] * sin_phi);
            value++;

            index  += freqStep;  if (index  >= SigSize) index  -= SigSize;
            index1 += 1;         if (index1 >= SigSize) index1 -= SigSize;
        }
    }
}

/*
 * GaborArrayMax — find the maximum modulus in a transform slice.
 * Fix 4: converted from K&R to ANSI prototype style.
 */
static void GaborArrayMax(FILTER *filter, double *value, int size,
                           long SigSize, int Log2SigSize, int octave,
                           int freq, int SampleRate_n,
                           double *modula, double *v_r, double *v_i,
                           int *index)
{
    int i;
    double m, *p_r, *p_i;

    /* suppress unused-parameter warnings for parameters used by caller
     * for context but not needed inside this function */
    (void)filter; (void)SigSize; (void)Log2SigSize;
    (void)octave; (void)freq; (void)SampleRate_n;

    *modula = 0.0;
    *v_r    = 0.0;
    *v_i    = 0.0;
    *index  = 0;
    p_r = value;
    p_i = value + size;
    for (i = 0; i < size; i++) {
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
 * GaborGetGaborNorm — compute the normalisation factor for a Gabor atom.
 * Fix 4: converted from K&R to ANSI prototype style.
 */
static double GaborGetGaborNorm(double modula, double c_r, double c_i,
                                 long SigSize, int num_octave,
                                 int octave, long freq, long translation)
{
    double norm, innR, innI, C;
    double cos2Phi, sin2Phi;

    if (freq == 0 || freq == (SigSize >> 1))
        return 1.0;

    /* Allocate persistent global INDEX objects (Fix 2: freed in GaborOperCleanup) */
    if (indexG1 == (INDEX)NULL) indexG1 = AllocIndex();
    if (indexG2 == (INDEX)NULL) indexG2 = AllocIndex();

    indexG1->octave    = indexG2->octave    = (double)octave;
    indexG1->position  = indexG2->position  = (double)translation;
    indexG1->id        = (double)freq;
    indexG2->id        = (double)(SigSize - freq);

    C = GaborGetCoeff(cur_norm, octave, octave, pfC);

    cos2Phi = (c_r*c_r - c_i*c_i) / modula;
    sin2Phi = 2.0 * c_r * c_i    / modula;

    GaborGetInnerProd(indexG1, indexG2, num_octave, cur_l, cur_h,
                      0, -2*freq, nBoundG, C,
                      pnAep, pnIep,
                      pfCE1, pfCE2, pfCE3, pfG,
                      &innR, &innI);

    norm = 2.0 / (1.0 + CMREAL(cos2Phi, sin2Phi, innR, innI));
    return sqrt(norm);
}

/*
 * GaborGetSizeTime — return the number of time samples at a given octave.
 * Fix 4: converted from K&R to ANSI prototype style.
 */
int GaborGetSizeTime(int octave, int nL, int nS)
{
    if (octave < nS)
        return 1 << nL;
    return 1 << (nL - octave + nS - 1);
}

/*
 * GaborGetSizeFreq — return the number of frequency samples at a given octave.
 * Fix 4: converted from K&R to ANSI prototype style.
 */
int GaborGetSizeFreq(int octave, int nL, int nS)
{
    if (octave > nL - nS)
        return 1 << (nL-1);
    return 1 << (octave + nS - 2);
}

/*
 * GaborGetWaveForm — reconstruct the waveform for a given atom.
 *
 * Inputs:
 *	word		atom descriptor (WORD)
 *	pfFilter	filter bank (FILTER *)
 *	N		signal size (int)
 *	L		log2(N) (int)
 *
 * Returns a newly allocated GABSIGNAL. Caller must free it.
 *
 * Fix 1: NULL check added after new_gabsignal().
 * Fix 4: converted from K&R to ANSI prototype style.
 * Fix 7: sincos() replaces separate cos()/sin() calls.
 */
GABSIGNAL GaborGetWaveForm(WORD word, FILTER *pfFilter, int N, int L)
{
    int i, octave, position, freq, index, index1;
    double cos_phi, sin_phi, *v, norm;
    GABSIGNAL gabsignal;

    octave   = (int)word->index->octave;
    position = (int)word->index->position;
    freq     = (int)word->index->id;
    norm     = word->value;

    /* Fix 7: sincos() replaces cos()/sin() on the same argument */
    sincos((double)word->index->phase, &sin_phi, &cos_phi);
    cos_phi *= norm;
    sin_phi *= norm;

    gabsignal = new_gabsignal(N);
    /* Fix 1: NULL check was missing */
    if (gabsignal == (GABSIGNAL)NULL) {
        fprintf(stderr, "GaborGetWaveForm(): new_gabsignal(%d) failed\n", N);
        return (GABSIGNAL)NULL;
    }
    v = gabsignal->values;

    /* Dirac basis */
    if (octave == 0) {
        v[position] = (fabs((double)word->index->phase) < 1.0e-5)
                      ? 1.0 : -1.0;
        return gabsignal;
    }

    /* Fourier basis */
    if (octave == L) {
        for (i = 0; i < N; i++) {
            index = (freq * i) % N;
            *v++  = pfFilter[L]->values[index]   * cos_phi
                  - pfFilter[L]->values[index+N]  * sin_phi;
        }
        return gabsignal;
    }

    /* Gabor basis */
    for (i = 0; i < N; i++) {
        index  = (freq * i) % N;
        index1 = (N - position + i) % N;
        if (index1 < pfFilter[octave]->size)
            *v++ = pfFilter[octave]->values[index1] *
                   (pfFilter[L]->values[index]  * cos_phi -
                    pfFilter[L]->values[index+N] * sin_phi);
        else if (N - index1 < pfFilter[octave]->size)
            *v++ = pfFilter[octave]->values[N-index1] *
                   (pfFilter[L]->values[index]  * cos_phi -
                    pfFilter[L]->values[index+N] * sin_phi);
        else
            *v++ = 0.0;
    }

    return gabsignal;
}

/*
 * getMaxFrmNewton — refine the selected atom using Newton's method
 * on the finer grid.
 *
 * Fix 3: perror() → fprintf (no return needed — default case only).
 * Fix 4: converted from K&R to ANSI prototype style.
 */
void getMaxFrmNewton(WORD word, GABSIGNAL *trans, FILTER *filter,
                     int MinOctave, int MaxOctave, int L,
                     int l, int h, int lc, int hc)
{
    long km, pm, kc, pc;
    long k, p, N;
    int jm, jc, j;
    int i, dkC, dpC, dkF, dpF, lc2, hc2, lf2, hf2;
    long badK, N2;
    double vr, vi, modula;
    double f[3][3], alpha, fm, coeff, nrm;

    j = (int)word->index->octave;
    if (j == 0 || j == L)
        return;

    k  = (long)word->index->id;
    p  = (long)word->index->position;
    N  = 1 << L;
    N2 = N >> 1;

    fm = 0.0;
    jm = j;
    km = k;
    pm = p;

    for (i = -1; i < 2; i++) {
        jc = j + i;
        if (jc == 0 || jc == L)
            continue;
        lc2 = SSFACTOR(jc, lc);
        hc2 = SSFACTOR(L-jc, hc);
        lf2 = SSFACTOR(jc, l);
        hf2 = SSFACTOR(L-jc, h);
        dpC = STEP(jc, lc2);
        dkC = STEP(L-jc, hc2);
        dpF = STEP(jc, lf2);
        dkF = STEP(L-jc, hf2);
        kc = k;
        pc = p;
        switch (getNbhd(f, trans[jc - MinOctave + 1]->values,
                        jc, &kc, &pc, lc2, hc2, L, i)) {
            case -2:
                continue;
            case -1:
                alpha = interpolation(f, &kc, &pc,
                                      (double)dkC, (double)dpC,
                                      (double)dkF, (double)dpF, N, -1);
                break;
            case 1:
                alpha = interpolation(f, &kc, &pc,
                                      (double)dkC, (double)dpC,
                                      (double)dkF, (double)dpF, N, 0);
                break;
            default:
                /* Fix 3: perror() → fprintf; can't return here since we're
                 * in a switch inside a loop — just report and continue */
                fprintf(stderr,
                        "getMaxFrmNewton(): internal error in getNbhd\n");
                continue;
        }
        badK = 1 << (L-jc);
        if ((kc > 0 && kc < badK) || (kc > (N2-badK) && kc < N2))
            continue;
        if (alpha > fm) {
            jm = jc;
            km = kc;
            pm = pc;
            fm = alpha;
        }
    }

    if (jm == j && km == k && pm == p)
        return;

    innSigWvForm(trans[0], filter, jm, km, pm, L, &vr, &vi);
    modula = vr*vr + vi*vi;
    nrm    = GaborGetGaborNorm(modula, vr, vi, N, L, jm, km, pm);
    modula = sqrt(modula);
    coeff  = modula * nrm;
    if (coeff < word->coeff)
        return;

    word->coeff          = coeff;
    word->value          = nrm;
    word->index->phase   = (vi > 0.0)
                           ?  acos(vr / modula)
                           : -acos(vr / modula);
    word->index->octave   = (double)jm;
    word->index->id       = (double)km;
    word->index->position = (double)pm;
}

/*
 * getNbhd — collect the 3x3 neighbourhood of modulus values around (k2, p2).
 * Fix 4: converted from K&R to ANSI prototype style.
 */
static int getNbhd(double f[3][3], double *vector, int j,
                   long *k, long *p, int l2, int h2, int L, int flag)
{
    long indxReal, indxImag;
    int Ljl2;
    long k2, p2, pw, kw, pBnd, L2, badK;
    int i, m;
    double vr, vi;

    L2   = 1 << (j+h2-1);
    Ljl2 = L - j + l2;
    pBnd = 1 << Ljl2;
    p2   = (*p) >> (j-l2);
    k2   = (*k) >> (L-j-h2);
    badK = 1 << h2;

    if ((k2 < badK && k2 > 0) || (k2 > (L2-badK) && k2 < L2))
        return -2;

    if (flag == 0) {
        indxReal    = TINDXREAL(k2, p2, Ljl2);
        indxImag    = TINDXIMAG(k2, p2, Ljl2);
        vr          = vector[indxReal];
        vi          = vector[indxImag];
        f[1][1]     = vr*vr + vi*vi;
    } else {
        f[1][1] = findMaxInNbdhd(vector, &k2, &p2, Ljl2, L2);
        if ((k2 < badK && k2 > 0) || (k2 >= (L2-badK) && k2 < L2))
            return -2;
    }

    if (k2 == badK || k2 == 0 || k2 == (L2-badK) || k2 == L2) {
        pw = p2 - 1;
        if (pw < 0)   pw += pBnd;
        else if (pw > pBnd) pw -= pBnd;
        indxReal = TINDXREAL(k2, pw, Ljl2);
        indxImag = TINDXIMAG(k2, pw, Ljl2);
        vr = vector[indxReal]; vi = vector[indxImag];
        f[1][0] = vr*vr + vi*vi;

        pw = p2 + 1;
        if (pw < 0)   pw += pBnd;
        else if (pw > pBnd) pw -= pBnd;
        indxReal = TINDXREAL(k2, pw, Ljl2);
        indxImag = TINDXIMAG(k2, pw, Ljl2);
        vr = vector[indxReal]; vi = vector[indxImag];
        f[1][2] = vr*vr + vi*vi;

        if (k2 == 0 || k2 == L2) {
            f[1][0] /= 2.0;
            f[1][2] /= 2.0;
        }
        return -1;
    }

    for (i = -1; i < 2; i++) {
        kw = k2 + i;
        for (m = -1; m < 2; m++) {
            if (i == 0 && m == 0)
                continue;
            pw = p2 + m;
            if (pw < 0)   pw += pBnd;
            else if (pw > pBnd) pw -= pBnd;
            indxReal = TINDXREAL(kw, pw, Ljl2);
            indxImag = TINDXIMAG(kw, pw, Ljl2);
            vr = vector[indxReal]; vi = vector[indxImag];
            f[i+1][m+1] = vr*vr + vi*vi;
        }
    }
    *k = k2 << (L-j-h2);
    *p = p2 << (j-l2);
    return 1;
}

/*
 * interpolation — Newton step to refine the maximum location.
 *
 * Fix 4: converted from K&R to ANSI prototype style.
 * Fix 8: while(ktmp<0)/while(ptmp<0) replaced with single if — the Newton
 *        correction is at most one period negative so a single addition suffices.
 */
static double interpolation(double f[3][3], long *k, long *p,
                             double dk, double dp, double dkF, double dpF,
                             long N, int flag)
{
    double fkk, fpp, fkp, fk, fp;
    double det, ktmp, ptmp;
    double fe;
    double dks, dps;

    dks = dk * dk;
    dps = dp * dp;

    if (flag < 0) {
        fkk = 1.0;
        fkp = fk = 0.0;
    } else {
        fkk = (f[2][1] - 2.0*f[1][1] + f[0][1]) / dks;
        fkp = (f[2][2] - f[2][0] - f[0][2] + f[0][0]) / (4.0*dk*dp);
        fk  = (f[2][1] - f[0][1]) / (2.0*dk);
    }
    fpp = (f[1][2] - 2.0*f[1][1] + f[1][0]) / dps;
    fp  = (f[1][2] - f[1][0]) / (2.0*dp);
    det = fkk*fpp - fkp*fkp;

    ktmp = (double)(*k) - (fpp*fk - fkp*fp) / det;
    ptmp = (double)(*p) - (fkk*fp - fkp*fk) / det;

    dk  = ktmp - (double)(*k);
    dks = dk * dk;
    dp  = ptmp - (double)(*p);
    dps = dp * dp;

    fe = f[1][1] + fk*dk + fp*dp + (fkk*dks + fpp*dps)/2.0 + fkp*dk*dp;

    /* Fix 8: single if replaces while loop — Newton correction is at most
     * one period away from zero */
    if (ktmp < 0.0) ktmp += (double)N;
    if (ptmp < 0.0) ptmp += (double)N;

    *k  = (long)(ktmp / dkF + 0.5) * (long)dkF;
    *p  = (long)(ptmp / dpF + 0.5) * (long)dpF;
    *k %= N;
    if (*k > (N >> 1)) *k = N - (*k);
    *p %= N;

    return fe;
}

/*
 * innSigWvForm — compute the inner product of a signal with a Gabor atom.
 *
 * Fix 4: converted from K&R to ANSI prototype style.
 * Fix 10: (k*(p+t))%N and (k*(p-t))%N replaced with running accumulators.
 */
static void innSigWvForm(GABSIGNAL gabsignal, FILTER *filter,
                          int j, long k, long p, int L,
                          double *vr, double *vi)
{
    double *vs, *vf, *vCos, *vSin, tmp;
    long t, N, indxE;

    N    = 1 << L;
    vs   = gabsignal->values;
    vf   = filter[j]->values;
    vCos = filter[L]->values;
    vSin = vCos + N;
    *vr  = 0.0;
    *vi  = 0.0;

    /* Fix 10: indxE = (k*(p+t))%N — k is added each step, wrap when >= N */
    indxE = (k * p) % N;   /* starting value at t=0 */
    for (t = 0; t <= N-p-1; t++) {
        if (t < filter[j]->size)
            tmp = vf[t];
        else if (N-t < filter[j]->size)
            tmp = vf[N-t];
        else
            tmp = 0.0;
        if (tmp != 0.0) {
            tmp  *= vs[t+p];
            *vr  += tmp * vCos[indxE];
            *vi  -= tmp * vSin[indxE];
        }
        indxE += k;
        if (indxE >= N) indxE -= N;
    }

    /* Fix 10: indxE = (k*(p-t))%N — k is subtracted each step */
    indxE = (k * (p-1)) % N;   /* starting value at t=1 */
    if (indxE < 0) indxE += N;
    for (t = 1; t <= p; t++) {
        if (N-t < filter[j]->size)
            tmp = vf[N-t];
        else if (t < filter[j]->size)
            tmp = vf[t];
        else
            tmp = 0.0;
        if (tmp != 0.0) {
            tmp  *= vs[p-t];
            *vr  += tmp * vCos[indxE];
            *vi  -= tmp * vSin[indxE];
        }
        indxE -= k;
        if (indxE < 0) indxE += N;
    }
}

/*
 * findMaxInNbdhd — find the maximum modulus in the 3x3 neighbourhood
 * of (k2, p2) and update (k2, p2) to point at it.
 * Fix 4: converted from K&R to ANSI prototype style.
 */
static double findMaxInNbdhd(double *vector, long *k2, long *p2,
                              int Ljl2, long L2)
{
    int k, p;
    long kw, pw, pBnd, pm, km;
    long indxReal, indxImag;
    double fm, f, vr, vi;

    fm   = 0.0;
    km   = *k2;
    pm   = *p2;
    pBnd = 1 << Ljl2;

    for (k = -1; k < 2; k++) {
        kw = (*k2) + k;
        if (kw < 0 || kw > L2)
            continue;
        for (p = -1; p < 2; p++) {
            pw       = ((*p2) + p + pBnd) % pBnd;
            indxReal = TINDXREAL(kw, pw, Ljl2);
            indxImag = TINDXIMAG(kw, pw, Ljl2);
            vr = vector[indxReal];
            vi = vector[indxImag];
            f  = vr*vr + vi*vi;
            if (kw == L2) f /= 2.0;
            if (f > fm) {
                pm = pw;
                km = kw;
                fm = f;
            }
        }
    }
    *k2 = km;
    *p2 = pm;
    return fm;
}

/*
 * end of gb_oper.c
 */
 