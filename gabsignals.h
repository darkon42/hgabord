/*..........................................................................*/
/*                                                                          */
/*   ------------------------Author Z. Zhang---------------------           */
/*         -------- (C) 1993 Copyright, All Right Reserved.--------         */
/*                                                                          */
/*..........................................................................*/
/****************************************************************************/
/*                                                                          */
/*  gabsignals.h    Definition of the GABSIGNAL structure.                  */
/*                  Included by mpp.h — do not include directly.            */
/*                                                                          */
/****************************************************************************/

/* Fix 11: include guard added — was missing, causing redefinition errors
 * on multiple inclusions. */
#ifndef GABSIGNALS_H
#define GABSIGNALS_H

#define SIG_SIZE 32768  /* maximum size of a gabsignal when read from an
                           ascii file (see gabsignal_io.c) */

/* Fix 9: STRING_SIZE is defined in mpp.h before this file is included.
 * Guard the definition here so it doesn't redefine and trigger a warning. */
#ifndef STRING_SIZE
#define STRING_SIZE 200
#endif

/*
 * Signal structure
 */
typedef struct gabsignal {
    int     size_alloca;        /* size of the 'values' allocation          */
    int     size;               /* number of active elements                */
    double  shift;              /* shift of the signal with respect to zero */
    double *values;             /* signal sample array                      */
    double  scale;              /* signal scale                             */
    int     firstp;             /* first index unaffected by left edge      */
    int     lastp;              /* last index unaffected by right edge      */
    double  param;              /* distance between uncorrelated samples    */
    char    name[STRING_SIZE];  /* signal name (reserved, not yet used)     */
} *GABSIGNAL;

/*
 * Fix 10: corrected function name (was farray_malloc — does not exist;
 *         the actual function is darray_malloc in sig_alloc.c).
 * Fix 10: proper prototypes replacing K&R-style unprototyped declarations.
 */
extern double    *darray_malloc(int size);
extern GABSIGNAL  new_struct_gabsignal(void);

#endif /* GABSIGNALS_H */