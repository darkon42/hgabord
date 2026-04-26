#define structh

#include "mpp.h"
#include "structpath.h"
#include "cwt1d.h"

#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef __unix__
    #include <sys/types.h>
    #include <sys/times.h>
    //#include <sys/param.h>
    #include <unistd.h>
    #include <strings.h>
#endif


#define DATAARRAY_SIZE 2048

/* The MAIN function */
int main()
{
	int bsize;
	
	double *bookB;
	int n, i;

	
	//int m=1024;
	double *data = NULL;
	
	extern void delete_gabsignal(GABSIGNAL gabsignal);
	extern FILTER *GaborFreeFilter(FILTER*, int);
	extern double *gabord(double *data, int SigSize, int *booksize);

	data = (double*) malloc(DATAARRAY_SIZE*sizeof(double));
	if (data==NULL) {
		fprintf(stderr, "alloc error\n");
		return(1);
	}

    int arr[DATAARRAY_SIZE];
    srand(time(0));
	 for (i = 0; i < DATAARRAY_SIZE; i++) data[i] = rand();
	

	for (i = 0; i < 100; i++)
	bookB=gabord(data, DATAARRAY_SIZE, &bsize);

	printf("Computation completed!\n");

	free(bookB);
	free(data);

	// releasing pointers
	/*free(pnIep);
	//free(pnAep);

	for (i=0; i<MAX_NUM_SB; i++)
	{
		free(library[i]);
	}
	
	free(pfB);
	free(pfC);
	free(pfG);

	for (i=0; i<MAX_NUM_GABSIGNAL; i++)
		delete_gabsignal(gabsignals[i]);
		

	for (i=0; i<MAX_NUM_SB; i++)
		if (transform[i] != (GABSIGNAL *)NULL)
			GaborFreeFilter(transform[i], TransAlloc[i]);

	for (i=0; i<MAX_NUM_SB; i++)
		if (filter[i] != (FILTER *)NULL)
			GaborFreeFilter(filter[i], num_filter[i]);
	
	delete_gabsignal(temporary);
	delete_gabsignal(gabsignal);

	free(pfCE1);
	free(pfCE2);
	free(pfCE3);
*/

	return(0);

}


