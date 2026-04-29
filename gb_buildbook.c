#include <stdio.h>
#include <math.h>
#include "mpp.h"
#include "cwt1d.h"

#ifdef __unix__
    #include <sys/types.h>
    #include <sys/times.h>
    #include <sys/param.h>
#else
    #define M_PI_2     1.57079632679489661923
#endif

/* real part of complex multiplication */
#define CMREAL(a1,b1,a2,b2)        ((a1)*(a2)-(b1)*(b2))
/* imaginary part of complex multiplication */
#define CMIMAGINARY(a1,b1,a2,b2)   ((a1)*(b2)+(a2)*(b1))

/* Fix 4: inline complex multiply — avoids recomputing the same four
 * products when CMREAL and CMIMAGINARY are called back-to-back on
 * identical arguments inside the hot k2 loop. */
static inline void cmul(double ar, double ai, double br, double bi,
                        double *rr, double *ri)
{
    *rr = ar*br - ai*bi;
    *ri = ar*bi + br*ai;
}

/* Fix 3: proper prototypes replacing K&R-style unprototyped declarations */
WORD       GaborGetMaxFrmTrans(GABSIGNAL *trans, FILTER *filter,
                                int MinOctave, int MaxOctave, int LnSize,
                                int SubsampleOctaveTime, int SubsampleOctaveFreq);
void       GaborGetResidue(GABSIGNAL *trans, FILTER *filter,
                            WORD word, int LnSize);
void       BookAppend(BOOK book, WORD word);
void       UpdateGabor(GABSIGNAL *trans, WORD word,
                        int MinOctave, int MaxOctave, int L,
                        int l, int h, int lf, int hf,
                        double *pfG, double *pfC, double *pfFilterNorm,
                        FILTER *pfFilter, double *pfB,
                        double *pfCE1, int *pnIep);
void       GaborUpdateFourier(GABSIGNAL *trans, FILTER *filter, WORD word,
                               int MinOctave, int MaxOctave, int L);
void       getMaxFrmNewton(WORD word, GABSIGNAL *trans, FILTER *filter,
                            int MinOctave, int MaxOctave, int LnSize,
                            int l, int h,
                            int SubsampleOctaveTime, int SubsampleOctaveFreq);

/*
 * Build book from a gabor transform
 *
 * Inputs:
 *	trans			gabor transform (GABSIGNAL *)
 *	filter			the basic gabor functions, used for updating (FILTER *)
 *	book			book that stores the results (BOOK)
 *	MinOctave		minimum octave of the decomposition (int)
 *	MaxOctave		maximum octave of the decomposition (int)
 *	LnSize			log2 of the gabsignal size (int)
 *	SubsampleOctaveTime	octave at which to begin subsampling in time (int)
 *	SubsampleOctaveFreq	octave at which to begin subsampling in frequency (int)
 *	l			subsampling rate for time (int)
 *	h			subsampling rate for frequency (int)
 *	pnIep			index array for pfG (int *)
 *	pfC			normalization factor array (double *)
 *	pfB			bandwidth array (double *)
 *	pfG			Gaussian array (double *)
 *	pfCE1			complex exponential table (double *)
 *	pfFilterNorm		filter normalization array (double *)
 *
 * Outputs:
 *	book			updated with the selected word appended
 *	trans			updated in-place with the residue transform
 */
void GaborBuildBook(GABSIGNAL *trans, FILTER *filter, BOOK book,
                    int MinOctave, int MaxOctave, int LnSize,
                    int SubsampleOctaveTime, int SubsampleOctaveFreq,
                    int l, int h, int *pnIep, double *pfC, double *pfB,
                    double *pfG, double *pfCE1, double *pfFilterNorm)
{
    WORD word;

    /* Fix 1: perror() doesn't stop execution — use fprintf+return */
    if (trans == (GABSIGNAL *)NULL || book == (BOOK)NULL) {
        fprintf(stderr, "GaborBuildBook(): null arguments!\n");
        return;
    }

    /* Get the maximum coefficient from the transform */
    word = GaborGetMaxFrmTrans(trans, filter, MinOctave, MaxOctave, LnSize,
                               SubsampleOctaveTime, SubsampleOctaveFreq);

    /* Fix 2: guard against GaborGetMaxFrmTrans returning NULL */
    if (word == (WORD)NULL) {
        fprintf(stderr, "GaborBuildBook(): GaborGetMaxFrmTrans returned NULL\n");
        return;
    }

    if (word->index->octave != 0.0
        && word->index->octave != (double)LnSize
        && (l > SubsampleOctaveTime || h > SubsampleOctaveFreq))
        getMaxFrmNewton(word, trans, filter, MinOctave, MaxOctave, LnSize,
                        l-1, h-1, SubsampleOctaveTime-1, SubsampleOctaveFreq-1);

    /* Subtract the selected atom from the transform (compute residue) */
    GaborGetResidue(trans, filter, word, LnSize);

    UpdateGabor(trans, word, MinOctave, MaxOctave, LnSize,
                SubsampleOctaveTime-1, SubsampleOctaveFreq-1, l-1, h-1,
                pfG, pfC, pfFilterNorm, filter, pfB, pfCE1, pnIep);

    GaborUpdateFourier(trans, filter, word, MinOctave, MaxOctave, LnSize);

    /* Append the selected word to the book */
    BookAppend(book, word);
}


/*
 * GaborUpdateFourier — update the stored Fourier transform of the signal
 * by subtracting the contribution of the selected word.
 */
void GaborUpdateFourier(GABSIGNAL *trans, FILTER *filter, WORD word,
                         int MinOctave, int MaxOctave, int L)
{
    double *v1r, *v1i, *v2r, *v2i, *vg;
    double coeff, tmp, phi1, cosPhi, sinPhi;
    long j1, k1, p1;
    long index1, index2;
    long k2, N;

    N    = (trans[0]->size) >> 1;
    j1   = (long)  word->index->octave;
    k1   = (long)  word->index->id;
    p1   = (long)  word->index->position;
    phi1 =         word->index->phase;
    coeff = word->coeff * word->value;
    v1r  = trans[MaxOctave - MinOctave + 2]->values;
    v1i  = v1r + N;

    if (j1 == 0)
        {
        /* the case of dirac */
        double sqrtN = sqrt((double)N); /* Fix 6: compute once */
        v2r = filter[L]->values;
        v2i = v2r + N;
        if (fabs(phi1) > M_PI_2)
            coeff = -coeff / sqrtN;
        else
            coeff /= sqrtN;
			
		double *v1r0 = v1r, *v1i0 = v1i;                                                                                                                                                             
		//#pragma omp parallel for private(index1)                  
		  for (k2 = 0; k2 < N; k2++)              
			  {                                                                                                                                                                                        
			  index1 = (k2 * p1) % N;
			  v1r0[k2] -= coeff * v2r[index1];                                                                                                                                                         
			  v1i0[k2] += coeff * v2i[index1];                      
			  }    
        }
    else if (j1 == L)
        {
        /* the case of Fourier */
        v1r[k1] = 0.0;
        v1i[k1] = 0.0;
        if (k1 != 0)
            {
            v1r[N - k1] = 0.0;
            v1i[N - k1] = 0.0;
            }
        }
    else
        {
        /* the case of gabor */
        /* Fix 5: sincos() replaces separate cos()/sin() calls */
        sincos(phi1, &sinPhi, &cosPhi);
        coeff /= 2.0;
        v2r = filter[L]->values;
        v2i = v2r + N;
        vg  = filter[L - j1]->values;

		double *v1r0 = v1r, *v1i0 = v1i;
		//#pragma omp parallel for private(index1, index2, tmp)

        for (k2 = 0; k2 < N; k2++)
            {
            double cr, ci;

            /* First term: k2 - k1 */
            index1 = k2 - k1 + N;
            index2 = (index1 * p1) % N;  /* wraps independently of index1 */
            index1 %= N;
            if (index1 < filter[L-j1]->size)
                tmp = vg[index1];
            else if ((N - index1) < filter[L-j1]->size)
                tmp = vg[N - index1];
            else
                tmp = 0.0;
            if (tmp != 0.0)
                {
                tmp *= coeff;
                /* Fix 4: cmul() replaces back-to-back CMREAL/CMIMAGINARY */
                cmul(v2r[index2], -v2i[index2], cosPhi, sinPhi, &cr, &ci);
                //*v1r -= tmp * cr;
                //*v1i -= tmp * ci;
				v1r0[k2] -= tmp * cr;
				v1i0[k2] -= tmp * ci;  
                }

            /* Second term: k2 + k1 */
            index1 = k2 + k1;
            index2 = (index1 * p1) % N;
            index1 %= N;
            if (index1 < filter[L-j1]->size)
                tmp = vg[index1];
            else if ((N - index1) < filter[L-j1]->size)
                tmp = vg[N - index1];
            else
                tmp = 0.0;
            if (tmp != 0.0)
                {
                tmp *= coeff;
                /* Fix 4: cmul() — note conjugate sinPhi for second term */
                cmul(v2r[index2], -v2i[index2], cosPhi, -sinPhi, &cr, &ci);
				
				v1r0[k2] -= tmp * cr;
				v1i0[k2] -= tmp * ci;
                //*v1r -= tmp * cr;
                //*v1i -= tmp * ci;
                }

            //v1r++;
            //v1i++;
            }
        }
}

/*
 * end of gb_buildbook.c
 */