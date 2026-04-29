/*..........................................................................*/
/*                                                                          */
/*       ------------------------Author Z. Zhang---------------------       */
/*         -------- (C) 1993 Copyright, All Right Reserved.--------         */
/*                                                                          */
/*..........................................................................*/
/****************************************************************************/
/*                                                                          */
/*  mpp.h     Basic include file for all the mpp files.                     */
/*            Has to be included in any new mpp file.                       */
/*                                                                          */
/****************************************************************************/

#ifndef MPP_H
#define MPP_H

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>   /* Fix 1: replaced obsolete <malloc.h> with <stdlib.h> */
                      /* Fix 2: removed dead #ifdef sparc / <alloca.h> block  */

/*
 * global constants
 */
#define NO              0
#define YES             1
#define TRUE            1
#define FALSE           0

#define STRING_SIZE     200
#define FILTERNAME_SIZE 10
#define FILTER_SIZE     128     /* max size of a filter */
#define MAX_NUM_GABSIGNAL  2    /* number of gabsignals */
#define MAX_NUM_SB      1       /* number of structure books */
#define MAX_NUM_FILTER  2       /* size of filter bank vector */
#define MAX_VALUE       1.0e38  /* max value */
#define MIN_VALUE       -1.0e38 /* min value */
#define MAX_NUM_ITERATION 100   /* number of iterations when building the book */
#define WAVELET_PATH    1
#define WHOLE_PATH      0
#define KEPT            0
#define DISCARD         1
#define OCTAVEMARK      2
#define FREQMARK        4
#define SHIFTMARK       8
#define PHASEMARK       16
#define CHIRP           32

/*
 * include file for gabsignal — must define STRING_SIZE first
 */
#include "gabsignals.h"

/*
 * constants for flags
 */
#define A_FLAG          0x1
#define B_FLAG          0x2
#define C_FLAG          0x4
#define D_FLAG          0x8
#define E_FLAG          0x10
#define F_FLAG          0x20
#define G_FLAG          0x40
#define H_FLAG          0x80
#define I_FLAG          0x100
#define J_FLAG          0x200
#define K_FLAG          0x400
#define L_FLAG          0x800
#define M_FLAG          0x1000
#define N_FLAG          0x2000
#define O_FLAG          0x4000
#define P_FLAG          0x8000
#define Q_FLAG          0x10000
#define R_FLAG          0x20000
#define S_FLAG          0x40000
#define T_FLAG          0x80000
#define U_FLAG          0x100000
#define V_FLAG          0x200000
#define W_FLAG          0x400000
#define X_FLAG          0x800000
#define Y_FLAG          0x1000000
#define Z_FLAG          0x2000000
#define a_FLAG          0x4000000
#define b_FLAG          0x8000000
#define c_FLAG          0x10000000
#define d_FLAG          0x20000000
#define e_FLAG          0x40000000
#define f_FLAG          0x80000000

/* Fix 5: PI2 macro was missing parentheses — 1.0/PI2 would expand to
 * 1.0/2.0*M_PI = M_PI/2.0 without them. */
#define PI2             (2.0 * M_PI)

/*
 * Useful macros
 */
/* returns YES iff inf <= x <= sup */
#define INRANGE(inf,x,sup)  ((inf) <= (x) && (x) <= (sup))
/* square of a number */
#define SQUARE(x)           ((x)*(x))
/* max of two numbers */
#define MAX(x,y)            ((x) > (y) ? (x) : (y))
/* min of two numbers */
#define MIN(x,y)            ((x) < (y) ? (x) : (y))
/* sign of a number */
#define SIGN(x)             ((x) < 0 ? (-1) : (x) > 0 ? 1 : 0)
/* masking equal operator */
#define MEQ(flag,mask)      (((flag) & (mask)) == (mask))
/* masking not operator */
#define MNOT(flag,mask)     ((((flag) & (mask)) == (mask)) ? (flag)^(mask) : (flag))

/*
 * Filter structure (alias for GABSIGNAL)
 */
typedef GABSIGNAL FILTER;

/*
 * Index structure
 */
typedef struct index {
    double id;       /* frequency of structure */
    double octave;   /* scale of structure     */
    double position; /* position of structure  */
    double phase;    /* phase of the structure */
} *INDEX;

/*
 * Word structure
 */
typedef struct word {
    double coeff;       /* coefficient value (real part) */
    double value;       /* value of the word */
    double coeff1;
    double value1;
    int    status;      /* either KEPT or DISCARD */
    INDEX  index;       /* frequency, position, and scale */
    struct word *next;  /* next node in book linked list */
} *WORD;

/*
 * Book structure
 */
typedef struct book {
    WORD   first;     /* first element in the list */
    WORD   last;      /* last element in the list */
    int    size;      /* number of words in book */
    double energy;    /* energy of book representation */
    double sigen;     /* energy of original signal */
    int    id;        /* book volume number */
    int    type;      /* book type */
    double smax;
    double smin;      /* max and min parameters of decomposition */
    int    sig_size;  /* size of the signal */
} *BOOK;

#define WAVELET       0
#define GABOR         1
#define QMF           2
#define DELTAFOURIER  3
#define NEWGABOR      4
#define FOURIER       5
#define DIRAC         6
#define WINDOWFOURIER 7

/*
 * global variables
 */
extern GABSIGNAL gabsignals[MAX_NUM_GABSIGNAL];
extern BOOK      library[MAX_NUM_SB];
extern GABSIGNAL cur_gabsignal;
extern int       cur_sig_size;
extern BOOK      old_cur_book;
extern GABSIGNAL *old_cur_filter;
extern int       Current_Book, Old_Book;
extern int       filter_type[MAX_NUM_SB];
extern GABSIGNAL *filter[MAX_NUM_SB];
extern int       cur_shift_octave;
extern int       cur_SOT;
extern int       cur_SOF;
extern int       old_cur_num_filter;
extern int       num_filter[MAX_NUM_SB];
extern int       plot_var;
extern FILE      *foutput;

#define cur_book         (library[Current_Book])
#define cur_filter       (filter[Current_Book])
#define cur_filter_type  (filter_type[Current_Book])
#define cur_num_filter   (num_filter[Current_Book])

#define BW 0
#define CP 1

/*
 * Fix 4: proper prototypes replacing K&R-style unprototyped extern declarations.
 * Fix 3: log2() extern declarations removed — log2() is standard in <math.h>
 *        since C99 and our local definition in math.c has been deleted.
 * Fix 6: dead #ifdef PROTOTYPES block removed (PROTOTYPES was never defined).
 * Fix 7: commented-out error() declaration removed.
 */
extern char      *check_format(const char *fmt);
extern GABSIGNAL *AllocFilter(int MaxOctave);
extern INDEX      AllocIndex(void);
extern WORD       AllocWord(void);
extern BOOK       AllocBook(void);
extern GABSIGNAL  new_gabsignal(int size);
extern void       warning(char *str);
extern void 	  clear_book(BOOK book);

/*
 * Fix 8: typo "traslation" corrected to "translation" in comment below.
 *
 * shift octave, subsample octave in translation and subsample octave in frequency
 */
extern int cur_shift_octave;
extern int cur_SOT;
extern int cur_SOF;

/*
 * end of mpp.h
 */

#endif /* MPP_H */