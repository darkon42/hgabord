/*******************************************************************/
/*  MPP gabsignal processing program.                                 */
/* (C) 1993 Copyright New York University, All Right Reserved.     */
/*  Francois Bergeaud, Mike Orszag.                                */
/*******************************************************************/

/****************************************************************************/
/*  gabsignal_alloc.c   Functions which deal with the dynamical                */
/*                   allocation of memory for GABSIGNAL's                      */
/****************************************************************************/

#include "mpp.h"
#include <stdlib.h>

/************************************/
/* allocation of an array of doubles */
/************************************/
double *darray_malloc(int size)
{
  double *ptr;
  if ((ptr = (double*)calloc(size,sizeof(double)))) return(ptr);    /*/ PJF*/
  else {
    fprintf(stderr,"Can't find enough memory for %d doubles\n",size);
    perror((char *)EXIT_FAILURE);
  }
  return(ptr);
}
  
/* there is also a macro using "alloca" instead of "malloc" defined
   in mpp.h and called farray_alloca but it does NOT test wether
   the allocation succeeded or not */

/*************************************************/
/* Create a new gabsignal structure and returns it  */
/*************************************************/
GABSIGNAL new_struct_gabsignal()
{
  GABSIGNAL gabsignal;
  void init_gabsignal(GABSIGNAL);

  if(!(gabsignal = (GABSIGNAL) (malloc(sizeof(struct gabsignal)))))
    perror("Mem. alloc for GABSIGNAL failed\n");

  init_gabsignal(gabsignal);

  return (gabsignal);
}

void init_gabsignal(GABSIGNAL gabsignal)
{
  if (gabsignal == (GABSIGNAL)NULL)
	perror("init_gabsignal(): null argument!");

  gabsignal->values = (double *)NULL;
  gabsignal->size_alloca = 0;
  gabsignal->size = 0;
  gabsignal->name[0] = '\0';
  gabsignal->scale = 1.;
  gabsignal->shift = 0.;
  gabsignal->firstp = 0;
  gabsignal->lastp = 0;
  gabsignal->param = 1.;
}


/*************************************************/
/* Desallocate the whole 'gabsignal' structure      */
/*************************************************/
void delete_gabsignal(GABSIGNAL gabsignal)
{
  if (gabsignal)
    {
    if (gabsignal->values) free((char *)gabsignal->values);
	free((char *)gabsignal);
    }
}


/*************************************************/
/* Create a gabsignal structure and an array of     */
/* double of size 'size'. The array is put in the */
/* 'values' field of the gabsignal.                 */
/* It returns the gabsignal                         */
/*************************************************/
GABSIGNAL new_gabsignal(int size)
{
  int i;

  GABSIGNAL gabsignal = new_struct_gabsignal();
  gabsignal->values = darray_malloc(size);
  gabsignal->size = size;
  gabsignal->size_alloca = size;
  gabsignal->lastp = size-1;
  gabsignal->param = 1.;
  for (i=0;i<size;i++)
	gabsignal->values[i] = 0.0;
  return(gabsignal);
}


/*************************************************/
/* Initialization of 'gabsignal' and desallocation  */
/* of the array of double gabsignal->values          */
/*************************************************/
void clear_gabsignal(GABSIGNAL gabsignal)
{
   void init_gabsignal(GABSIGNAL);

/*   if (gabsignal == (GABSIGNAL)NULL)
	error("clear_gabsignal(): null argument!"); */

   if (gabsignal->values)
	{
	free((char *)gabsignal->values);
	init_gabsignal(gabsignal);
	}
}
    

/*************************************************/
/* Return YES if gabsignal->values is not NULL      */
/*************************************************/
int is_empty(GABSIGNAL gabsignal)
{
  if (gabsignal->values) {return(NO);} else return(YES);
}


/*************************************************/
/* Very important procedure which has to be      */
/* called each time one needs to store doubles    */
/* in a gabsignal.                                  */
/* If 'gabsignal' has already a size > 'size'       */
/* it doesn't do any memory allocation.          */
/* If not then it first desallocates             */
/* gabsignal->values and replaces it by an array    */
/* of size 'size.                                */
/* In all cases, it initializes the fields       */
/* of 'gabsignal'.                                  */
/*************************************************/
void change_gabsignal(GABSIGNAL gabsignal, int size)
{
  int i;

  if (gabsignal->size_alloca < size)
    {
      clear_gabsignal(gabsignal);
      gabsignal->values = darray_malloc(size);
      gabsignal->size_alloca = size;
      if (gabsignal->values == NULL)
	printf ("malloc(%d) failed!\n",size);
    }

  gabsignal->size = size;
  gabsignal->scale = 1.;
  gabsignal->shift = 0.;
  gabsignal->firstp = 0;
  gabsignal->lastp = size-1;
  gabsignal->param = 1.;
  for(i = 0; i < size; i++)
    gabsignal->values[i] = 0.0;
}
/*
 * free GABSIGNAL structure and return NULL
 */
GABSIGNAL FreeSignal(GABSIGNAL gabsignal)
{
    delete_gabsignal(gabsignal);
    return((GABSIGNAL)NULL);
}
