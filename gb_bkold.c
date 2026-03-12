#include "mpp.h"
#include "cwt1d.h"
/*
 *
 * Build book from a gabor transform
 *
 * Inputs:
 *	trans		gabor transform (GABSIGNAL *)
 *	filter		the basic gabor fuctions, used for updating
 *			(FILTER *)
 *	threshold	precision that stop the loop (double)
 *	num_octave	number of octaves in the transform
 *	sig_norm	the L2 norm of the original gabsignal (double)
 *	num_iter	current number of iterations (int)
 *	max_num_iter	number of iterations allowed in the loop
 *	ShiftOctave	the octave that begins to use the second formula
 *			compute the projection
 *	SubsampleOctaveTime	the octave that begins to subsample in
 *				translation
 *	SubsampleOctaveFreq	the octave that begins to subsample in
 *				frequency
 *
 * Outputs:
 *	book		book that stores the results (BOOK)
 *
 */
void GaborBuildBookOld(GABSIGNAL *trans, FILTER *filter, BOOK book, int num_octave,int ShiftOctave, int SubsampleOctaveTime, int SubsampleOctaveFreq,double LamdaNoise)
{
    WORD word;
    WORD GaborGetMaxFrmTrans();
    GABSIGNAL *GaborDecomp();
    void BookAppend(), GaborGetResidue();

    if (trans == (GABSIGNAL *)NULL || book == (BOOK)NULL)
	perror("GaborBuildBook(): null arguments!");
/*
 * get the maximum from the trans and put it into word
 */
    word = GaborGetMaxFrmTrans(trans,filter, 1, num_octave - 1, num_octave,SubsampleOctaveTime,
		SubsampleOctaveFreq);
/*
 * get the residue
 */
    GaborGetResidue(trans,filter,word,num_octave);
/*
 * put the word into book
 */
    BookAppend(book,word);
/*
 * calculate the residue transform
 */
   trans = GaborDecomp(trans,trans[0],filter,
		SubsampleOctaveTime,
		SubsampleOctaveFreq,
		num_octave,ShiftOctave);
    return;
}
/*
 * end of ng_buildbook.c
 */
