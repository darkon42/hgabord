#include "mpp.h"

//extern void change_gabsignal();
extern void change_gabsignal(GABSIGNAL gabsignal, int size);

/*
 * calculate a gabsignal created by the L2 normalized gaussian
 *
 * Inputs:
 * 	min_t	minimum value for t (double)
 *	max_t	maximum value for t (double)
 *	size_g	number of points in t need to be calculated (int)
 *	size_s  size of the gabsignal need to be created (int),
 *	sigma	the deviation for the gaussian
 * 
 * Remark:
 *	if size_s > size_g, then the gabsignal filled with zero for
 *	those points that index number greater than size_g
 *
 * Bugs:
 *	size_s must be greater than size_g, if not return null pointer
 *
 */
GABSIGNAL Gaussian2Signal(double min_t,
					   double max_t,
					   double size_g,
					   int size_s,
					   double sigma)
{
    GABSIGNAL gabsignal=(GABSIGNAL)NULL;
    double t, scale;
    int i;
	double gaussianL2(double t,double sigma);
//    GABSIGNAL new_gabsignal();

    if (size_s < size_g)
	return((GABSIGNAL)NULL);

    gabsignal = new_gabsignal(size_s);
    scale = (double)(max_t-min_t)/(double)size_g;
    gabsignal->scale = (double)scale;

    t = (double)min_t;

    for (i=0;i<=(int)size_g;i++)
	{
	gabsignal->values[i] = gaussianL2(t,(double)sigma);
	t += scale;
	}

    return(gabsignal);
}


/**************************************/
/* Copy a gabsignal in another           */
/**************************************/

/* copy an array of double ('input') of 
   size 'sigsize' in another ('output') */
/*--------------------------------------------------------------------------*/
void farray_copy(double *input, int sigsize, double *output)
{
  double *to, *from;

  for (to=output,from = input; from < input+sigsize; to++, from++)
		*to = *from;
}

/*--------------------------------------------------------------------------*/
/*
 * tranlate an array of double
 */
/*--------------------------------------------------------------------------*/
void farray_translate(double *f_in,double *f_out,int size,int shift)
{
    int i, i0, i1;
    double *value;

    if (f_in == (double *)NULL || f_out == (double *)NULL)
	perror("farray_translate(): null input!");

    while (shift < 0)
	shift += size;
    while (shift-size>0)
	shift -= size;
    i0 = size-shift;
    i1 = i0+size;
    value = f_out;
    for (i=i0;i<i1;i++)
	*value++ = f_in[i%size];
}

/*--------------------------------------------------------------------------*/
/* Signal Copy */
/*--------------------------------------------------------------------------*/
int sig_copy(GABSIGNAL input,GABSIGNAL output)
{
  int i;

  if (input == (GABSIGNAL)NULL || output == (GABSIGNAL)NULL)
	perror("sig_copy(): null input!");

  change_gabsignal(output,input->size);

  for(i=0;i<input->size;i++)
    output->values[i] = input->values[i];
  output->scale = input->scale;
  output->shift = input->shift;
  output->firstp = input->firstp;
  output->lastp = input->lastp;
  output->param = input->param;
  
  return(0);
}

/*--------------------------------------------------------------------------*/
/*
 * gabsignal copy from range min to max
 */
/*--------------------------------------------------------------------------*/
int sig_range_copy(GABSIGNAL input,GABSIGNAL output,int min,int max)
{
  int i, size;
  double *value;

  if (input == (GABSIGNAL)NULL || output == (GABSIGNAL)NULL)
	perror("sig_range_copy(): null input!");
  if (min<0 || max<min || max>input->size)
	perror("sig_range_copy(): illegal argument!");

  size = max-min;
  if (output->size != size)
	change_gabsignal(output,size);

  value = output->values;
  for (i=min;i<max;i++)
	*value++ = input->values[i];
  output->scale = input->scale;
  output->shift = input->shift;
  output->firstp = min;
  output->lastp = max-1;
  output->param = input->param;
  
  return(0);
}



/*--------------------------------------------------------------------------*/
/*
 * append a word into book
 */
/*--------------------------------------------------------------------------*/
void BookAppend(BOOK book, WORD word)
{
/* checking the inputs */

if ( book == NULL || word == NULL )
	perror("BookAppend: NULL input");	


  if(book->last == NULL) {
    book->first = word;
    book->last = book->first;
  }
  else {
    book->last->next = word;
    book->last = book->last->next;
  }
  book->size = book->size + 1;
  book->energy += word->coeff*word->coeff;
}

void error(char *str)
{
perror(str);
}

void error_option(char *str)
{
perror("Option error\n");
}

void warning(char *str)
{
   printf("%s", str);
}

