/*..........................................................................*/
/*                                                                          */
/*      ------------------------------------------------------------*/
/*      (C) 1993 Copyright New York University, All Right Reserved.         */
/*	Modify by Zhifeng Zhang and Mike Orszag, 1992                       */
/*                                                                          */
/*..........................................................................*/
/****************************************************************************/
/*                                                                          */
/*  gabsignal_functions10.c                                                    */
/*                   Miscellaneous useful functions on gabsignals:             */
/*                       Input : 1 gabsignal                                   */
/*                       Output: 0 gabsignal                                   */
/*                                                                          */
/****************************************************************************/
#include "mpp.h"

extern double sig_mean(GABSIGNAL input); /* computes the mean of a gabsignal */
extern double sig_variance(GABSIGNAL input); /* computes the variance of a gabsignal */

/*********************************/
/* Put a gabsignal to zero          */
/*********************************/
int sig_zero(GABSIGNAL gabsignal)
{
  int j;

  for (j = 0; j < gabsignal->size; j++)
    gabsignal->values[j] = 0.0;
	
   return(0);
}

/*********************************/
/* Compute the absolute value of */
/* a gabsignal.                     */
/*********************************/
int sig_abs(GABSIGNAL gabsignal)
{
  int j;

  for (j = 0; j < gabsignal->size; j++)
    gabsignal->values[j] = (double) fabs(gabsignal->values[j]);

  return(0);
}
/*********************************/
/* Compute the log of a gabsignal   */
/*********************************/
int sig_log(GABSIGNAL gabsignal)
{
	int j;
	
	for (j=0; j < gabsignal->size; j++)
	  gabsignal->values[j] = (double) log(fabs(gabsignal->values[j]));

return(0);
}

/*
 * compute the exp of a gabsignal
 */
int sig_exp(GABSIGNAL gabsignal)
{
    int j;

    if (gabsignal==(GABSIGNAL)NULL)
	perror("sig_exp(): null input!");
    if (gabsignal->values==(double *)NULL)
	perror("sig_exp(): null point!");

    for (j=0; j < gabsignal->size; j++)
	gabsignal->values[j] = (double)exp(gabsignal->values[j]);
return(0);
}
/*********************************/
/* Compute the square of a gabsignal*/
/*********************************/
int sig_square(GABSIGNAL gabsignal)
{
  int j;

  for (j = 0; j < gabsignal->size; j++)
    gabsignal->values[j] = SQUARE(fabs(gabsignal->values[j]));

return(0);
}



/*********************************/
/* Add a number to the gabsignal    */
/*********************************/
void sig_add_num(GABSIGNAL gabsignal, double num)
{
  int j;

  for (j = 0; j < gabsignal->size; j++)
    gabsignal->values[j] += num;
}

/*********************************/
/* Add a number to the gabsignal    */
/*********************************/
int sig_sub_num(GABSIGNAL gabsignal, double num)
{
  int j;

  for (j = 0; j < gabsignal->size; j++)
    gabsignal->values[j] -= num;

return(0);
}


/*********************************/
/* Add a number to the gabsignal    */
/*********************************/
int sig_mult_num(GABSIGNAL gabsignal, double num)
{
  int j;

  for (j = 0; j < gabsignal->size; j++)
    gabsignal->values[j] *= num;

return(0);
}


/*********************************/
/* Add a number to the gabsignal    */
/*********************************/
void sig_div_num(GABSIGNAL gabsignal, double num)
{
  int j;

  if (!num)
    return;
  for (j = 0; j < gabsignal->size; j++)
    gabsignal->values[j] /= num;
}

/*********************************/
/* Compute the position of the   */
/* the min and the max of a      */
/* gabsignal                        */
/*********************************/
int sig_pos_min_max(GABSIGNAL gabsignal,double *pmin,double *pmax)
{
  int j;
  double vmin, vmax;
  
  vmin = 999999;
  vmax = -999999;
  for(j=0;j<gabsignal->size;j++)
    {
      if (vmax < gabsignal->values[j]) {
	vmax = gabsignal->values[j];
	*pmax = (double)j;
      }
      if (vmin > gabsignal->values[j]) {
	vmin = gabsignal->values[j];
	*pmin = (double)j;
      }
    }
return(0);
}
/*********************************/
/* Computes the min and the max  */
/* of a gabsignal                   */
/*********************************/
int sig_min_max(GABSIGNAL gabsignal,double *vmin,double *vmax)
{
  int j;
  
  (*vmin) = MAX_VALUE;
  (*vmax) = MIN_VALUE;
  for(j=0;j<gabsignal->size;j++)
    {
      if ((*vmax) < gabsignal->values[j]) (*vmax) = gabsignal->values[j];
      if ((*vmin) > gabsignal->values[j]) (*vmin) = gabsignal->values[j];
    }
return(0);
}
/*********************************/
/* Compute the mean of a gabsignal  */
/*********************************/
double sig_mean(GABSIGNAL input)
{
  int j;
  double sum;

  sum = 0;
  for (j = 0; j < input->size; j++)
    sum += input->values[j];
  return(sum/input->size);
}
/*********************************/
/* Compute the variance of a     */
/* gabsignal                        */
/*********************************/
double sig_variance(GABSIGNAL input)
{
  GABSIGNAL input_sq = new_struct_gabsignal();
  double var;
  int sig_copy(GABSIGNAL input,GABSIGNAL output);
  void delete_gabsignal(GABSIGNAL gabsignal);
  
  var = - SQUARE(sig_mean(input));
  sig_copy(input,input_sq);
  sig_square(input_sq);
  var += sig_mean(input_sq);
  delete_gabsignal(input_sq);
  return(var);
}
/*********************************/
/* Compute the integral of a     */
/* gabsignal                        */
/*********************************/
double integral(GABSIGNAL gabsignal)
{
  int j;
  double sum;

  sum = 0;
  for (j = 0; j < gabsignal->size; j++)
    sum += gabsignal->values[j] ;
  return(sum);
}
/*********************************/
/* Compute the L2 norm of a      */
/* gabsignal                        */
/*********************************/
double sig_L2_norm_sq(GABSIGNAL input)
{
  GABSIGNAL input_sq = new_struct_gabsignal();
  double norm;
  double integral(GABSIGNAL gabsignal);
  int sig_copy(GABSIGNAL input,GABSIGNAL output);
  void delete_gabsignal(GABSIGNAL gabsignal);
  
  sig_copy(input,input_sq);
  sig_square(input_sq);
  norm = integral(input_sq);
  delete_gabsignal(input_sq);
  return(norm);
}
/*
 * compute the L2 square norm of a real array
 */
double farray_L2_sq_norm(double *value,int size)
{
   int i;
   double *v, norm=0.0;

   v = value;
   for (i=0;i<size;i++)
	{
	norm += (*v)*(*v);
	v++;
	}

   return(norm);
}
/*
 * compute the L2 norm of a real array
 */
double farray_L2_norm(double *value,int size)
{
   int i;
   double *v, norm=0.0;

   v = value;
   for (i=0;i<size;i++)
	{
	norm += (*v)*(*v);
	v++;
	}
   norm = (double)sqrt((double)norm);

   return(norm);
}
/*
 * find the max abs value in the gabsignal and its position
 */
void SigAbsMx(GABSIGNAL gabsignal,double *MxValue,double *position)
{
    int i;

    if (gabsignal == (GABSIGNAL)NULL)
	perror("SigAbsMx(): null argument!");
    *MxValue = 0.0;
    *position = 0.0;

    for (i=0;i<gabsignal->size;i++)
	if (fabs(gabsignal->values[i])>fabs(*MxValue))
	    {
	    *MxValue = gabsignal->values[i];
	    *position = (double)i;
	    }
}
