/*******************************************************************/
/*  MPP gabsignal processing program.                              */
/* (C) 1993 Copyright New York University, All Right Reserved.     */
/*  Francois Bergeaud, Mike Orszag.                                */
/*******************************************************************/

/****************************************************************************/
/*  sig_alloc.c   Functions which deal with the dynamical allocation        */
/*                of memory for GABSIGNALs                                  */
/****************************************************************************/

#include "mpp.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * darray_malloc — allocate and zero-initialise an array of doubles.
 *
 * Returns a pointer to the allocated array, or NULL on failure.
 * Caller is responsible for freeing the returned pointer.
 *
 * Fix 4: perror((char*)EXIT_FAILURE) was undefined behaviour (cast of
 *        integer 1 to char* passed as string). Replaced with perror(NULL)
 *        since fprintf already printed the descriptive message.
 * Fix 5: added return NULL so the caller receives NULL instead of a
 *        broken pointer after the error message.
 */
double *darray_malloc(int size)
{
    double *ptr;

    ptr = (double *)calloc((size_t)size, sizeof(double));
    if (ptr == (double *)NULL) {
        fprintf(stderr,
                "darray_malloc: can't allocate %d doubles\n", size);
        /* Fix 4: was perror((char*)EXIT_FAILURE) — undefined behaviour */
        return (double *)NULL;
    }
    return ptr;
}

/*
 * init_gabsignal — zero-initialise all fields of an allocated GABSIGNAL.
 *
 * Fix 5: perror() → fprintf+return so execution stops on NULL input
 *        instead of crashing on gabsignal->values below.
 */
void init_gabsignal(GABSIGNAL gabsignal)
{
    /* Fix 5: perror() → fprintf+return */
    if (gabsignal == (GABSIGNAL)NULL) {
        fprintf(stderr, "init_gabsignal(): null argument!\n");
        return;
    }

    gabsignal->values      = (double *)NULL;
    gabsignal->size_alloca = 0;
    gabsignal->size        = 0;
    gabsignal->name[0]     = '\0';
    gabsignal->scale       = 1.0;
    gabsignal->shift       = 0.0;
    gabsignal->firstp      = 0;
    gabsignal->lastp       = 0;
    gabsignal->param       = 1.0;
}

/*
 * new_struct_gabsignal — allocate and initialise a bare GABSIGNAL struct
 * (without a values array).
 *
 * Returns the allocated struct, or NULL on failure.
 * Caller is responsible for eventually calling delete_gabsignal().
 *
 * Fix 1/5: perror() → fprintf+return NULL so the caller receives NULL
 *          instead of crashing on init_gabsignal(NULL) below.
 */
GABSIGNAL new_struct_gabsignal(void)
{
    GABSIGNAL gabsignal;

    gabsignal = (GABSIGNAL)malloc(sizeof(struct gabsignal));
    /* Fix 1/5: perror() → fprintf+return NULL */
    if (gabsignal == (GABSIGNAL)NULL) {
        fprintf(stderr,
                "new_struct_gabsignal: malloc failed\n");
        return (GABSIGNAL)NULL;
    }

    init_gabsignal(gabsignal);
    return gabsignal;
}

/*
 * new_gabsignal — allocate a GABSIGNAL struct with a zero-initialised
 * values array of 'size' doubles.
 *
 * Returns the allocated GABSIGNAL, or NULL on failure.
 * Caller is responsible for calling delete_gabsignal() when done.
 *
 * Fix 2/3: added NULL check after new_struct_gabsignal(), and free the
 *          struct if darray_malloc fails so no memory is leaked.
 * Fix 8: removed redundant zero-init loop — darray_malloc uses calloc
 *        which already zero-initialises the allocation.
 */
GABSIGNAL new_gabsignal(int size)
{
    GABSIGNAL gabsignal;

    gabsignal = new_struct_gabsignal();
    /* Fix 3: NULL check was missing */
    if (gabsignal == (GABSIGNAL)NULL)
        return (GABSIGNAL)NULL;

    gabsignal->values = darray_malloc(size);
    /* Fix 2: free the struct if values allocation fails */
    if (gabsignal->values == (double *)NULL) {
        free(gabsignal);
        return (GABSIGNAL)NULL;
    }

    gabsignal->size        = size;
    gabsignal->size_alloca = size;
    gabsignal->lastp       = size - 1;
    gabsignal->param       = 1.0;
    /* Fix 8: calloc in darray_malloc already zeroed the array */

    return gabsignal;
}

/*
 * delete_gabsignal — free the values array and the GABSIGNAL struct itself.
 */
void delete_gabsignal(GABSIGNAL gabsignal)
{
    if (gabsignal != (GABSIGNAL)NULL) {
        if (gabsignal->values != (double *)NULL)
            free(gabsignal->values);
        free(gabsignal);
    }
}

/*
 * clear_gabsignal — free the values array and re-initialise all fields,
 * but keep the struct itself allocated.
 *
 * Fix 6: the NULL check was commented out — restored. A NULL gabsignal
 *        would crash on gabsignal->values below.
 */
void clear_gabsignal(GABSIGNAL gabsignal)
{
    /* Fix 6: NULL check was commented out */
    if (gabsignal == (GABSIGNAL)NULL) {
        fprintf(stderr, "clear_gabsignal(): null argument!\n");
        return;
    }

    if (gabsignal->values != (double *)NULL) {
        free(gabsignal->values);
        init_gabsignal(gabsignal);
    }
}

/*
 * is_empty — return YES if gabsignal->values is NULL, NO otherwise.
 */
int is_empty(GABSIGNAL gabsignal)
{
    return (gabsignal->values != (double *)NULL) ? NO : YES;
}

/*
 * change_gabsignal — resize a GABSIGNAL to hold 'size' doubles.
 *
 * If the existing allocation is large enough it is reused (and zeroed).
 * Otherwise the old array is freed and a new one allocated.
 *
 * Fix 7: error message was going to stdout via printf — changed to
 *        fprintf(stderr,...) consistent with the rest of the file.
 * Fix 9: zero-init loop only runs when reusing the existing buffer.
 *        When reallocating, darray_malloc's calloc already zeroes it.
 */
void change_gabsignal(GABSIGNAL gabsignal, int size)
{
    if (gabsignal->size_alloca < size) {
        /* Need a larger buffer — free the old one and allocate fresh */
        clear_gabsignal(gabsignal);
        gabsignal->values = darray_malloc(size);
        if (gabsignal->values == (double *)NULL) {
            /* Fix 7: was printf() to stdout */
            fprintf(stderr,
                    "change_gabsignal: darray_malloc(%d) failed!\n", size);
            return;
        }
        gabsignal->size_alloca = size;
        /* Fix 9: calloc in darray_malloc already zeroed the new buffer */
    } else {
        /* Reuse existing buffer — must zero it explicitly */
        memset(gabsignal->values, 0, (size_t)size * sizeof(double));
    }

    gabsignal->size   = size;
    gabsignal->scale  = 1.0;
    gabsignal->shift  = 0.0;
    gabsignal->firstp = 0;
    gabsignal->lastp  = size - 1;
    gabsignal->param  = 1.0;
}

/*
 * FreeSignal — free a GABSIGNAL and return NULL so the caller can write:
 *   sig = FreeSignal(sig);
 * to both free and NULL their pointer in one step.
 */
GABSIGNAL FreeSignal(GABSIGNAL gabsignal)
{
    delete_gabsignal(gabsignal);
    return (GABSIGNAL)NULL;
}

/*
 * end of sig_alloc.c
 */