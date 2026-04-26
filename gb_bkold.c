#include <stdio.h>
#include "mpp.h"
#include "cwt1d.h"

/*
 * Build book from a gabor transform
 *
 * Inputs:
 *	trans			gabor transform (GABSIGNAL *)
 *	filter			the basic gabor functions, used for updating (FILTER *)
 *	book			book that stores the results (BOOK)
 *	num_octave		number of octaves in the transform (int)
 *	ShiftOctave		the octave that begins to use the second formula (int)
 *	SubsampleOctaveTime	the octave that begins to subsample in translation (int)
 *	SubsampleOctaveFreq	the octave that begins to subsample in frequency (int)
 *	LamdaNoise		noise threshold — NOTE: currently unused, reserved for
 *				future use as a stopping criterion (double)
 *
 * Outputs:
 *	book			updated with the selected word appended
 *	trans			updated in-place with the residue transform
 */

/* Fix 5: proper prototypes instead of unprototyped K&R-style declarations */
WORD      GaborGetMaxFrmTrans(GABSIGNAL *trans, FILTER *filter,
                               int min_octave, int max_octave, int num_octave,
                               int SubsampleOctaveTime, int SubsampleOctaveFreq);
void      GaborGetResidue(GABSIGNAL *trans, FILTER *filter, WORD word, int num_octave);
void      BookAppend(BOOK book, WORD word);
GABSIGNAL *GaborDecomp(GABSIGNAL *trans, GABSIGNAL gabsignal, FILTER *filter,
                        int SubsampleOctaveTime, int SubsampleOctaveFreq,
                        int MinOctave, int MaxOctave, int ShiftOctave);

void GaborBuildBookOld(GABSIGNAL *trans, FILTER *filter, BOOK book,
                       int num_octave, int ShiftOctave,
                       int SubsampleOctaveTime, int SubsampleOctaveFreq,
                       double LamdaNoise)
{
    WORD word;

    /* Fix 3: perror() doesn't stop execution — use fprintf+return */
    if (trans == (GABSIGNAL *)NULL || book == (BOOK)NULL) {
        fprintf(stderr, "GaborBuildBookOld(): null arguments!\n");
        return;
    }

    /* Fix 4: LamdaNoise is accepted but not yet used; documented above.
     * Suppress the unused-parameter warning explicitly. */
    (void)LamdaNoise;

    /* Get the maximum coefficient from the transform */
    word = GaborGetMaxFrmTrans(trans, filter, 1, num_octave - 1, num_octave,
                               SubsampleOctaveTime, SubsampleOctaveFreq);

    /* Subtract the selected atom from the transform (compute residue) */
    GaborGetResidue(trans, filter, word, num_octave);

    /* Append the selected word to the book */
    BookAppend(book, word);

    /* Recompute the Gabor transform of the residue signal.
     * trans contents are modified in-place; the return value is the same
     * pointer (GaborDecomp only returns a different pointer when trans is
     * NULL on entry, which is guarded above). */
    GaborDecomp(trans, trans[0], filter,
                SubsampleOctaveTime, SubsampleOctaveFreq,
                1, num_octave, ShiftOctave);
}

/*
 * end of gb_bkold.c
 */