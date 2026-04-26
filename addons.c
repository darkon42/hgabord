#include "mpp.h"
#include <stdio.h>
#include <string.h>

extern void change_gabsignal(GABSIGNAL gabsignal, int size);

/*
 * calculate a gabsignal created by the L2 normalized gaussian
 *
 * Inputs:
 *	min_t	minimum value for t (double)
 *	max_t	maximum value for t (double)
 *	size_g	number of points in t need to be calculated (int)
 *	size_s	size of the gabsignal need to be created (int)
 *	sigma	the deviation for the gaussian
 *
 * Remark:
 *	if size_s > size_g, then the gabsignal is filled with zero for
 *	those indices greater than size_g.
 *
 * Returns NULL if size_s < size_g.
 *
 * Fix 2: loop was `i <= size_g` (off-by-one, wrote size_g+1 values into
 *        a size_s-element buffer when size_s == size_g). Changed to `<`.
 */
GABSIGNAL Gaussian2Signal(double min_t,
                           double max_t,
                           double size_g,
                           int    size_s,
                           double sigma)
{
    GABSIGNAL gabsignal;
    double t, scale;
    int i;
    double gaussianL2(double t, double sigma);

    if (size_s < (int)size_g)
        return (GABSIGNAL)NULL;

    gabsignal = new_gabsignal(size_s);
    if (gabsignal == (GABSIGNAL)NULL) {
        fprintf(stderr, "Gaussian2Signal: new_gabsignal(%d) failed\n", size_s);
        return (GABSIGNAL)NULL;
    }

    scale = (max_t - min_t) / size_g;
    gabsignal->scale = scale;
    t = min_t;

    /* Fix 2: was `i <= (int)size_g` — off-by-one heap overwrite when
     * size_s == size_g. The loop now writes exactly size_g values. */
    for (i = 0; i < (int)size_g; i++) {
        gabsignal->values[i] = gaussianL2(t, sigma);
        t += scale;
    }

    return gabsignal;
}


/*--------------------------------------------------------------------------*/
/*
 * farray_copy — copy an array of doubles.
 *
 * Fix 5: replaced the hand-written pointer loop with memcpy, which the
 * C library implements with SIMD on every modern platform.
 */
/*--------------------------------------------------------------------------*/
void farray_copy(double *input, int sigsize, double *output)
{
    memcpy(output, input, (size_t)sigsize * sizeof(double));
}


/*--------------------------------------------------------------------------*/
/*
 * farray_translate — circular-shift an array of doubles by `shift` positions.
 *
 * Fix 1: perror() was used for the null check but doesn't stop execution;
 *        replaced with fprintf+return.
 *
 * Fix 6: the original loop used `i % size` on every iteration and two
 *        while-loops to normalise shift. Since i runs over exactly one
 *        full period, the modulo wraps exactly once. We replace the whole
 *        thing with two memcpy calls — no division, no branching.
 *
 *        Shift normalisation: ((shift % size) + size) % size maps any
 *        integer shift into [0, size) in constant time.
 */
/*--------------------------------------------------------------------------*/
void farray_translate(double *f_in, double *f_out, int size, int shift)
{
    int tail;

    if (f_in == NULL || f_out == NULL) {
        fprintf(stderr, "farray_translate(): null input!\n");
        return;
    }

    /* Fix 6: branchless normalisation replaces two while loops */
    shift = ((shift % size) + size) % size;

    if (shift == 0) {
        memcpy(f_out, f_in, (size_t)size * sizeof(double));
        return;
    }

    /* The translated array is: f_in[size-shift .. size-1] ++ f_in[0 .. size-shift-1]
     * i.e. two contiguous chunks — no per-element modulo needed. */
    tail = size - shift;
    memcpy(f_out,        f_in + tail, (size_t)shift * sizeof(double));
    memcpy(f_out + shift, f_in,       (size_t)tail  * sizeof(double));
}


/*--------------------------------------------------------------------------*/
/*
 * sig_copy — copy a GABSIGNAL into another.
 *
 * Fix 1: perror() doesn't stop execution; replaced with fprintf+return.
 * Fix 7: element-by-element copy replaced with memcpy.
 */
/*--------------------------------------------------------------------------*/
int sig_copy(GABSIGNAL input, GABSIGNAL output)
{
    if (input == (GABSIGNAL)NULL || output == (GABSIGNAL)NULL) {
        fprintf(stderr, "sig_copy(): null input!\n");
        return -1;
    }

    change_gabsignal(output, input->size);

    /* Fix 7: was a manual element loop — memcpy is SIMD-optimised */
    memcpy(output->values, input->values, (size_t)input->size * sizeof(double));

    output->scale  = input->scale;
    output->shift  = input->shift;
    output->firstp = input->firstp;
    output->lastp  = input->lastp;
    output->param  = input->param;

    return 0;
}


/*--------------------------------------------------------------------------*/
/*
 * sig_range_copy — copy a sub-range [min, max) of a GABSIGNAL.
 *
 * Fix 1: perror() doesn't stop execution; replaced with fprintf+return.
 * Fix 3: the bounds check used input->size after a null check that didn't
 *        return, so a NULL input would have crashed on the second check.
 *        Now there is a single guarded block: return before any dereference.
 * Fix 8: element-by-element copy replaced with memcpy.
 */
/*--------------------------------------------------------------------------*/
int sig_range_copy(GABSIGNAL input, GABSIGNAL output, int min, int max)
{
    int size;

    /* Fix 3: all guard checks before any dereference */
    if (input == (GABSIGNAL)NULL || output == (GABSIGNAL)NULL) {
        fprintf(stderr, "sig_range_copy(): null input!\n");
        return -1;
    }
    if (min < 0 || max < min || max > input->size) {
        fprintf(stderr, "sig_range_copy(): illegal argument"
                        " (min=%d max=%d size=%d)\n", min, max, input->size);
        return -1;
    }

    size = max - min;
    if (output->size != size)
        change_gabsignal(output, size);

    /* Fix 8: was a manual element loop — memcpy is SIMD-optimised */
    memcpy(output->values, input->values + min, (size_t)size * sizeof(double));

    output->scale  = input->scale;
    output->shift  = input->shift;
    output->firstp = min;
    output->lastp  = max - 1;
    output->param  = input->param;

    return 0;
}


/*--------------------------------------------------------------------------*/
/*
 * BookAppend — append a WORD to a BOOK.
 *
 * Fix 4: perror() doesn't stop execution; a NULL book or word would have
 *        caused a segfault or silent corruption on the very next line.
 *        Now returns immediately on bad input.
 */
/*--------------------------------------------------------------------------*/
void BookAppend(BOOK book, WORD word)
{
    if (book == NULL || word == NULL) {
        fprintf(stderr, "BookAppend: NULL input\n");
        return;
    }

    if (book->last == NULL) {
        book->first = word;
        book->last  = word;
    } else {
        book->last->next = word;
        book->last       = word;
    }

    book->size   += 1;
    book->energy += word->coeff * word->coeff;
}


/*--------------------------------------------------------------------------*/
/* Utility wrappers                                                          */
/*--------------------------------------------------------------------------*/

void error(char *str)
{
    fprintf(stderr, "%s\n", str);
}

void error_option(char *str)
{
    fprintf(stderr, "Option error: %s\n", str);
}

void warning(char *str)
{
    fprintf(stderr, "Warning: %s\n", str);
}