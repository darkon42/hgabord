
/*..........................................................................*/
/*                                                                          */
/*   ------------------------Author Z. Zhang--------------------- */
/*         -------- (C) 1993 Copyright, All Right Reserved.--------         */
/*                                                                          */
/*..........................................................................*/



/****************************************************************************/
/*                                                                          */
/*  gabsignals.h        Definition of the GABSIGNAL structure                     */
/*                   Is included in mpp.h                            */
/*                                                                          */
/****************************************************************************/


#define SIG_SIZE 32768 /* maximum size of a gabsignal when  (bylo 22000)
			  read from an ascii file (see file
			  gabsignal_io.c) */
			  
#define STRING_SIZE     200

/***************************/
/*  Signal structure       */
/***************************/

typedef struct gabsignal{
  int size_alloca;         /* size of the allocation of the 'values' field*/
  int  size;               /* size of gabsignal */
  double shift;             /* shifting of the gabsignal with respect to zero */
  double *values;           /* gabsignal values */
  double scale;             /* gabsignal scale */
  int firstp;              /* index of the first point not 
			      affected by left side effect */
  int lastp;               /* index of the last point not 
			      affected by right side effect */
  double param;             /* distance between two successive 
			      uncorrelated points */
  char name[STRING_SIZE];  /* name of the gabsignal (not used yet) */
} *GABSIGNAL;


/* Functions in gabsignal_alloc.c */

extern double *farray_malloc(); /* allocate an array of double */
extern GABSIGNAL new_struct_gabsignal(); /* allocate a gabsignal structure */








