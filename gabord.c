/*************************************************************************/
/*  GABOR gabsignal processing program.                                     */
/* (C) 2001-2026 Copyright Johns Hopkins University, All Right Reserved.    */
/*                                                                       */
/*  Christophe JOUNY: Based on Mallat/Zhang sources for Matching Pursuit */
/*                    Adapted for continuous processing of EEGs          */
/*                    Added new energy threshold criteria                */
/*  Christophe JOUNY: mex-file Matlab Interface (2010)                   */
/*                                                                       */
/*  gabord.c   Matching Pursuit decomposition adapted for processing     */
/*             multiple gabsignals and long datasets                        */
/*************************************************************************/

#define structh

#include "mpp.h"
#include "cwt1d.h"

#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#include <fftw3.h>

#ifdef __unix__
    #include <sys/types.h>
    #include <sys/times.h>
    //#include <sys/param.h>
    #include <unistd.h>
    #include <strings.h>
#endif



//#define M_PI 3.14159265358979323846264338327


double *values=(double *)NULL;          /*values of gabsignal (temp var.) */

GABSIGNAL gabsignal=(GABSIGNAL)NULL;             /*gabsignal*/
// in global scope at top of gabord.c, alongside gabsignal:
GABSIGNAL sigtmp = (GABSIGNAL)NULL;

double *cur_norm;         /* = (double *)NULL;*/

int cur_MaxOctave;         /* Max octave used in the decomposition over gabor functions */
int cur_MinOctave;         /* Min octave used in the decomposition over gabor functions */

double *pfG=(double *)NULL, *pfCE1=(double *)NULL;
double *pfCE2=(double *)NULL, *pfCE3=(double *)NULL, *pfC=(double *)NULL;
double *pfB = (double *)NULL;
int *pnAep=(int *)NULL, *pnIep=(int *)NULL;

int nBoundG=4;         /* bound for the Gaussian */

int cur_l, cur_h;         /* oversubsampling octaves for the fine grid */
int sig_size = 32768;

/* input/output */
FILE *foutput;

/* library of books */
BOOK library[MAX_NUM_SB];
/* current book */
/*BOOK cur_book;*/
/* previous book */
BOOK old_cur_book;
/* array of gabsignals */
GABSIGNAL gabsignals[MAX_NUM_GABSIGNAL];
/* current gabsignal */
GABSIGNAL cur_gabsignal;
/* current gabsignal size */
int cur_sig_size=0;
/* filters */
GABSIGNAL *filter[MAX_NUM_SB];
/* current filter type */
int filter_type[MAX_NUM_SB];
/* shift octave and subsample octave */
int cur_shift_octave=0;
int cur_SOT=1;
int cur_SOF=1;
/* previous book */
GABSIGNAL *old_cur_filter;
/* number of filter */
int num_filter[MAX_NUM_SB];
/*  global variable for transformation */
GABSIGNAL *(transform[MAX_NUM_SB]);           /* GD  7/11/93 */

int Current_Book = 0;
/* number of transformation */
int TransAlloc[MAX_NUM_SB];
/*  temporary Signal */
GABSIGNAL temporary;

double epslonG=1.0e-15;         /* threshold for the Gaussian fuctions */

int cpct=0, citer=0, cth=0, ccoh=0; /* all criteria disabled*/
double th_pct;                      /* criterion for terminate analysis based on prct of NRJ */
int max_num_iter=500;                   /*  idem by number of atoms*/
double thatomnrj, last_atom_nrj;    /* threshold of atoms energy and last atoms added nrj */

/* Declared extern in mpp.h; defined here so clients need not supply them. */
int Old_Book          = 0;
int plot_var          = 0;
int old_cur_num_filter = 0;

/* gabord_initialised was a static local inside gabord(); moved to file scope
 * so it can be declared threadprivate for multi-channel parallel execution. */
int gabord_initialised = 0;

/* Per-thread private state for gabord() — each OpenMP thread gets its own
 * independent copy of all algorithm state so gabord() can be called
 * concurrently from different threads without data races.
 * Threadprivate variables for cross-file globals are declared in mpp.h. */
#ifdef _OPENMP
/* gabord_initialised and the per-call signal/filter state are private per
 * thread.  Stopping criteria (cpct, citer, …) are shared but written to the
 * same values by all threads, so the data race is benign in practice.
 * last_atom_nrj is written each iteration; make it private so threads don't
 * corrupt each other's stopping decisions. */
#pragma omp threadprivate(\
    gabord_initialised, gabsignal, sigtmp, values,\
    cur_MaxOctave, cur_MinOctave, last_atom_nrj)
#endif

/***************************************************************************************************************************/
/***************************************************************************************************************************/

double *gabord(double *data, int SigSize, int *booksize)
{

	int ShiftOctave=0;
	int LnSigSize;                 /* power of 2 of gabsignal size */
	int SubsampleOctaveTime=2;
	int SubsampleOctaveFreq = 2;
	int SOT_build;
	int SOF_build;
	double sigma;
	double epslon;
	int old_cur_MinOctave = 0, old_cur_MaxOctave = 0; /* previous values of the levels of decomposition cur_MinOctave and cur_MaxOctave */

	int lndata, Lnlndata;                                                           /* length of data required by user, half of that and power of two*/
	int sigN;                                                                               /* gabsignal size*/

	BOOK book;                                                                              /*temp book for output*/
	WORD word;                                                                              /*temp word from the book*/
	INDEX indx;                                                                             /*temp index from the word*/

	//GABSIGNAL sigtmp=(GABSIGNAL)NULL;             /*gabsignal temp for decomp.*/

	/*DECLARATIONS FROM GBUILDBOOK*/

	/* the information for the fine grid */
	int curMaxLH; /* max of cur_l and cur_h */

	double SigEng=-1.0;
	double LamdaNoise; /*factor=1.0, */
	long i;
	int iter;
	static int nL=7, maxlh=0;
#ifdef _OPENMP
	/* nL/maxlh cache which size the threadprivate pfG/pfB/pfC arrays were
	 * built for; they must be per-thread to stay consistent with those
	 * (threadprivate) buffers when gabord() runs concurrently. */
#pragma omp threadprivate(nL, maxlh)
#endif
	double res_n, res_n1, orgN=-1.0;
	int l, h;
	unsigned long flag=N_FLAG;
	int sb_index;

	/* EXTERNAL FUNCTIONS*/

	extern FILTER *GaborBuildFilter(int, int, double, double**), *GaborFreeFilter(FILTER*, int);
	//extern double log2();
	extern GABSIGNAL *GaborDecomp(GABSIGNAL*,GABSIGNAL, FILTER*,int,int, int, int, int);
	extern void change_gabsignal(GABSIGNAL, int);
	extern void sig_add_num(GABSIGNAL, double);
	extern double sig_mean(GABSIGNAL);

	/* building of the book */
	extern double *GaborGetGaussianArray(int, int*, int, int);
	extern double *GaborGetCE1(int), *GaborGetCE2(int), *GaborGetCE3(int);
	extern double *GaborGetCArray(int), *GaborGetB(int, double);
	extern int GaborGetBound();
	extern void GaborGetNAep(int, int, double, int*, int*), GaborUpdateTrans();
	extern double farray_L2_sq_norm(double*, int);
	extern BOOK AllocBook();
	extern BOOK init_book(BOOK);
	extern WORD WordListFree(WORD);
	extern void GaborBuildBook(GABSIGNAL*,FILTER*,BOOK,int,int,int,int, int, int, int, int*, double*, double*,double*, double*, double*);
	extern void GaborBuildBookOld();
	
	extern void delete_gabsignal(GABSIGNAL gabsignal);

	extern int find2power(int);
	extern void farray_copy(double*, int, double*);

	int nib=0;
	double *Fbook;


/* DEFAULTS SETTINGS*/

	lndata=(int) SigSize;

	if ( (cpct==0) & (ccoh==0) & (cth==0) & (citer==0)) {citer=1; max_num_iter=1000; }

/* initializations     from int_loop.c  */

 /* Fix F: only initialise globals on first call (per thread when parallel).
  * Worker threads in an OpenMP parallel region have UNDEFINED initial values
  * for threadprivate variables unless copyin is used.  Explicitly zero all
  * index/pointer globals before any read so we never dereference garbage. */
	if (!gabord_initialised) {
		/* Safe defaults for all per-thread index variables */
		Current_Book      = 0;
		cur_sig_size      = 0;
		cur_shift_octave  = 0;
		cur_SOT = 1; cur_SOF = 1;
		cur_l   = 0; cur_h   = 0;
		old_cur_num_filter = 0;
		old_cur_book   = (BOOK)NULL;
		old_cur_filter = (GABSIGNAL *)NULL;
		cur_gabsignal  = (GABSIGNAL)NULL;

		for (sb_index = 0; sb_index < MAX_NUM_SB; sb_index++) {
			filter[sb_index]      = (GABSIGNAL *)NULL;
			num_filter[sb_index]  = 0;
			TransAlloc[sb_index]  = 0;
			filter_type[sb_index] = -1;
			transform[sb_index]   = (GABSIGNAL *)NULL;
		}
		/* allocation of book */
		for (i = 0; i < MAX_NUM_SB; i++) {
			library[i]     = AllocBook();
			library[i]->id = i;
		}
		/* allocation of gabsignals  */
		for (i = 0; i < MAX_NUM_GABSIGNAL; i++)
			gabsignals[i] = new_struct_gabsignal();

		/* Now safe to read cur_book = library[Current_Book=0] */
		old_cur_book   = cur_book;
		cur_gabsignal  = gabsignals[0];
		temporary      = new_struct_gabsignal();
		old_cur_filter = filter[0];
		Current_Book   = 0;  /* redundant but explicit */

		gabord_initialised = 1;
	} else {
		/* On repeat calls: clear the book from the previous run so it's
		 * ready to receive new atoms, but keep filter/transform intact. */
		if (cur_book != (BOOK)NULL) {
			clear_book(cur_book);
			init_book(cur_book);
		}
	}

	//Added to prevent re-init w/o free on recall
	if (gabsignal != (GABSIGNAL)NULL)
		delete_gabsignal(gabsignal);
	
	/* Fix: round up to next power of 2 BEFORE allocating so size_alloca
	 * matches the actual buffer; sig_mean/sig_add_num iterate to
	 * gabsignal->size and would overflow the original (smaller) allocation.
	 * Indices SigSize..lndata-1 are zeroed by calloc inside darray_malloc. */
	Lnlndata = find2power(lndata);
	lndata = 1<<Lnlndata;
	gabsignal = new_gabsignal(lndata);
	/* was:
	gabsignal = new_gabsignal((int)SigSize);
	Lnlndata = find2power(lndata);
	lndata = 1<<Lnlndata;
	gabsignal->size_alloca = (int) SigSize;
	gabsignal->size = lndata; */
	gabsignal->shift = 0;
	gabsignal->scale = 1;
	gabsignal->firstp = 0;
	gabsignal->lastp= gabsignal->size - 1;
	gabsignal->param = 1;

	/*Read Channel Data*/
	/*nbc=get_double_chan(in+nbinp, gabsignal->values, n, ichan);*/
	/*INPUT GABSIGNAL*/
	for (i=0; i<SigSize; i++) gabsignal->values[i]=data[i];

	sig_add_num(gabsignal, -1*sig_mean(gabsignal));                       /* substract mean to center around zero*/
																			/* should try the median to alleviate artifacts effects */
	orgN = (double) gabsignal->size;
	sigN = gabsignal->size;
	LnSigSize = find2power(gabsignal->size);
	cur_MaxOctave = LnSigSize - 1;
	cur_SOT=1; cur_SOF=1;
	ShiftOctave = 0;

	/* previous values for the decomposition */

	old_cur_MinOctave = cur_MinOctave;
	old_cur_MaxOctave = cur_MaxOctave;
	cur_MinOctave = 1;
	SOT_build = cur_SOT + 2;
	SOF_build = cur_SOF + 2;

	SigSize = ((long long)(1))<<LnSigSize;
	if (SigSize != sigN)                                                   /* gabsignal resize to a power of 2*/
	{
		values = (double *)malloc(sizeof(double) * sigN);
		if (values == (double *)NULL) {
			fprintf(stderr, "gabord(): malloc failed for temp buffer\n");
			delete_gabsignal(gabsignal);
			gabsignal = (GABSIGNAL)NULL;
			return (double *)NULL;
		}
		farray_copy(gabsignal->values, sigN, values);
		change_gabsignal(gabsignal, (int)SigSize);
		for (i = 0; i < sigN; i++)
			gabsignal->values[i] = values[i];
		for (i = (long)sigN; i < (long)SigSize; i++)
			gabsignal->values[i] = 0.0;
		free(values);
		values = (double *)NULL;
		sigN = (int)SigSize;
	}

	if (SubsampleOctaveTime > LnSigSize)
		SubsampleOctaveTime = LnSigSize;
	if (SubsampleOctaveFreq> LnSigSize)
		SubsampleOctaveFreq = LnSigSize;
	if (ShiftOctave == 0)
		ShiftOctave = (LnSigSize+SubsampleOctaveTime-SubsampleOctaveFreq+1)/2;

	sigma = (double)(1.0/sqrt((double)(2.0*M_PI)));                 /* 1/sqrt(2*pi) */

	if (cur_num_filter!=LnSigSize+1 || cur_filter_type != NEWGABOR)
	{                                       /* free the gabor filter */
		if (cur_filter != (FILTER *)NULL)
			cur_filter = GaborFreeFilter(cur_filter,cur_num_filter);
		if (cur_transform != (GABSIGNAL *)NULL)
			cur_transform = GaborFreeFilter(cur_transform, cur_TransAlloc);
	}

	if ( ( old_cur_MinOctave != cur_MinOctave ) || ( old_cur_MaxOctave != cur_MaxOctave ) )         /* octave have changed */
	{
		if (cur_transform != (GABSIGNAL *)NULL)
			cur_transform = GaborFreeFilter(cur_transform, cur_TransAlloc);
	}

	/* Build the (threadprivate) filter bank once per thread.
	 * GaborBuildFilter writes only threadprivate/local state and creates no
	 * FFTW plans, so it needs no lock.  FFTW plan creation happens later in
	 * GaborDecomp -> get_plans(), which serialises the planner itself.
	 *
	 * The old fftw_wisdom.dat import/export (and the global critical that
	 * wrapped it) were removed: get_plans() now uses FFTW_ESTIMATE, so there
	 * is no wisdom file to read, write, or race on — the previous code read
	 * the file from disk on EVERY call and serialised all threads through one
	 * lock, which destroyed scaling under parallel use. */
	if (cur_filter == (FILTER *)NULL)
	{
		cur_filter = GaborBuildFilter(LnSigSize, SigSize, sigma, &cur_norm);
		cur_num_filter = LnSigSize + 1;
		cur_filter_type = NEWGABOR;
	}

	/* Do the decomposition */

	if (sigtmp == (GABSIGNAL)NULL)
		sigtmp = new_gabsignal(SigSize<<1);
	else if (sigtmp->size != SigSize)
		change_gabsignal(sigtmp,SigSize<<1);
		
	//for (i=0; i<SigSize; i++) sigtmp->values[i] = gabsignal->values[i]; % being replace by memcpy
	memcpy(sigtmp->values, gabsignal->values, (size_t)SigSize * sizeof(double));

	cur_transform = GaborDecomp(cur_transform,sigtmp,cur_filter,
	                            SubsampleOctaveTime,SubsampleOctaveFreq,
	                            cur_MinOctave, cur_MaxOctave,ShiftOctave);

/*
 * set the last ShiftOctave and SubsampleOctaveTime and SigSize
 * to avoid recomputation of the filters
 */
	cur_TransAlloc = cur_MaxOctave - cur_MinOctave +3;
	cur_shift_octave = ShiftOctave;
	cur_SOT = SubsampleOctaveTime;
	cur_SOF = SubsampleOctaveFreq;
	cur_sig_size = (int)SigSize;
	cur_gabsignal = gabsignal;

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

	epslon = epslonG;
	l=cur_SOT+2;
	h=cur_SOF+2;

	cur_l = l-1;
	cur_h = h-1;

	if (cur_transform == (GABSIGNAL *)NULL || cur_filter_type != NEWGABOR)
		perror("Run gdecomp first!");

	if (cur_book==(BOOK)NULL)
		cur_book = AllocBook();
	if (cur_book->first != (WORD)NULL)
		cur_book->first = WordListFree(cur_book->first);
	init_book(cur_book);

	if (cur_transform[0] == (GABSIGNAL)NULL)
		perror("Internal Error!");

	SigEng = farray_L2_sq_norm(cur_transform[0]->values, cur_transform[0]->size>>1);

	cur_book->sigen = res_n = SigEng;
	cur_book->sig_size = sigN;
	cur_book->type = NEWGABOR;


	if (find2power(cur_sig_size) != nL)
	{
		nL = find2power(cur_sig_size);
		if (pnAep != (int *)NULL) {
			free((char *)pnAep);
			pnAep = (int *)NULL;
		}
		if (pnIep != (int *)NULL) {
			free((char *)pnIep);
			pnIep = (int *)NULL;
		}
		if (pfG != (double *)NULL) {
			free((char *)pfG);
			pfG = (double *)NULL;
		}
		if (pfCE1 != (double *)NULL) {
			free((char *)pfCE1);
			pfCE1 = (double *)NULL;
		}
		if (pfCE2 != (double *)NULL) {       /* Fix 6: was missing */
			free((char *)pfCE2);
			pfCE2 = (double *)NULL;
		}
		if (pfCE3 != (double *)NULL) {       /* Fix 6: was missing */
			free((char *)pfCE3);
			pfCE3 = (double *)NULL;
		}
		if (pfC != (double *)NULL) {
			free((char *)pfC);
			pfC = (double *)NULL;
		}
		if (pfB != (double *)NULL) {
			free((char *)pfB);
			pfB = (double *)NULL;
		}
	}

	if ((flag&N_FLAG)==N_FLAG)
	{
		if (maxlh!=MAX(MAX(cur_SOT-1,cur_SOF-1),MAX(l-1,h-1)))
		{
			curMaxLH = maxlh = MAX(MAX(cur_SOT-1,cur_SOF-1),MAX(l-1,h-1));
			if (pfG != (double *)NULL)
			{
				free((char *)pfG);
				pfG = (double *)NULL;
			}
		}

		if (pfB == (double *)NULL)
			pfB = GaborGetB(nL,epslon);

		if (pfC==(double *)NULL)
			pfC = GaborGetCArray(nL);

		if (pfG==(double *)NULL)
		{
			if (pnIep==(int *)NULL)
				pnIep = (int *)malloc(sizeof(int)*nL);
			if (pnAep==(int *)NULL)
				pnAep = (int *)malloc(sizeof(int)*nL);
			if (pnAep==(int *)NULL || pnIep==(int *)NULL) {
				perror("mem. alloc. failed!");
				return NULL;
				}
			GaborGetNAep(nL,maxlh,epslon,pnAep,pnIep);
			pfG = GaborGetGaussianArray(nL,pnAep,pnIep[nL-1],maxlh);
		}

		if (pfCE1==(double *)NULL)
			pfCE1 = GaborGetCE1(nL);
		if (pfCE2==(double *)NULL)
			pfCE2 = GaborGetCE2(nL);
		if (pfCE3==(double *)NULL)
			pfCE3 = GaborGetCE3(nL);
	}

/*STOP CRITERION*/

/* compute the lamda square for white noise */

	if (SigEng < 0.0)
		perror("GaborBuildBook(): empty gabsignal!");
	if (cur_l != 3 || cur_h != 3)
		LamdaNoise = 0.0;                 /* this lambda is not available yet */
	else
		LamdaNoise = (2.076747-0.089091*log(orgN))*sqrt(log(orgN)/orgN);


	res_n1 = 0.0;
	//printf("NRJ of the signal: %f\n", SigEng);
	last_atom_nrj = SigEng;                             /* to pass the first test : add xtof*/

	//printf("Max Iter. = %d\n", max_num_iter);

	for (iter=0; iter<max_num_iter; iter++)                                         /* Iteration criteria*/
	{

		if (cpct)
			if (cur_book->energy>=SigEng*th_pct/100.0)                      /* Threshold criteria*/
				break;

		if (ccoh)
			if (sqrt(1.0-res_n1/res_n)<LamdaNoise)                          /* Coherence criteria*/
				break;

		if (cth)
		{
			if (last_atom_nrj < thatomnrj)                                        /* Atoms NRJ criteria*/
				break;
		}

		{
			int prev_size = cur_book->size;
			GaborBuildBook(cur_transform,cur_filter,cur_book,
			               cur_MinOctave,
			               cur_MaxOctave,
			               nL,
			               cur_SOT,cur_SOF,l,h,
			               pnIep,pfC,pfB,pfG,pfCE1,
			               cur_norm);
			/* If no atom was appended (e.g. illegal translation),
			 * the signal is depleted — further iterations will spin. */
			if (cur_book->size == prev_size)
				break;
		}

		last_atom_nrj = cur_book->last->coeff*cur_book->last->coeff;
		
		//printf("NRJ: %f\n", last_atom_nrj);

		if (iter>0) res_n = res_n1;
		res_n1 = res_n - last_atom_nrj;

	}
	//printf("Iteration: %d - NRJ: %f\n", iter, last_atom_nrj);
	
	book = cur_book;
	*booksize=book->size;
	Fbook=(double *)malloc(5*book->size*sizeof(double));
	word=book->first;

	while (word!=NULL)
	{
		indx=word->index;
		Fbook[nib*5]=(double)indx->octave;
		Fbook[nib*5+1]=(double)indx->id;
		Fbook[nib*5+2]=(double)indx->position;
		Fbook[nib*5+3]=(double)word->coeff;
		Fbook[nib*5+4]=(double)indx->phase;

		nib++;
		word=word->next;
	}
	//double GAD=(double)(2.0*nib/SigSize);
	//printf("GAD=%15.10f\n", GAD);
	
	//delete_gabsignal(sigtmp);

/*OUPUT BOOK*/
	return Fbook;
}
/* END of Gabor*/
/*
 * gabord_cleanup — release all heap memory allocated by the library.
 * Call once at program exit after the last gabord() call.
 * FFTW plans are freed automatically via atexit(); this handles the rest.
 */
void gabord_cleanup(void)
{
    int i;
    extern void   delete_gabsignal(GABSIGNAL gabsignal);
    extern FILTER *GaborFreeFilter(FILTER *, int);
    extern void   GaborOperCleanup(void);

    free(pnIep);  pnIep = NULL;
    free(pnAep);  pnAep = NULL;

    for (i = 0; i < MAX_NUM_SB; i++) {
        clear_book(library[i]);
        free(library[i]);
    }

    free(pfB);   pfB   = NULL;
    free(pfC);   pfC   = NULL;
    free(pfG);   pfG   = NULL;
    free(pfCE1); pfCE1 = NULL;
    free(pfCE2); pfCE2 = NULL;
    free(pfCE3); pfCE3 = NULL;
    free(cur_norm); cur_norm = NULL;

    for (i = 0; i < MAX_NUM_GABSIGNAL; i++)
        delete_gabsignal(gabsignals[i]);

    for (i = 0; i < MAX_NUM_SB; i++)
        if (transform[i] != (GABSIGNAL *)NULL)
            GaborFreeFilter(transform[i], TransAlloc[i]);

    for (i = 0; i < MAX_NUM_SB; i++)
        if (filter[i] != (FILTER *)NULL)
            GaborFreeFilter(filter[i], num_filter[i]);

    delete_gabsignal(temporary);
    delete_gabsignal(gabsignal);

    GaborOperCleanup();
}

/*
 * gabord_reset — force re-initialisation of all per-thread state.
 * Call this once per worker thread (inside a parallel region) before
 * the first gabord() call when using --jobs N > 1.  Without it, new
 * OpenMP threads have undefined initial values for threadprivate variables.
 */
void gabord_reset(void)
{
    int i;
    Current_Book      = 0;
    cur_sig_size      = 0;
    cur_shift_octave  = 0;
    cur_SOT           = 1;
    cur_SOF           = 1;
    cur_l             = 0;
    cur_h             = 0;
    old_cur_num_filter = 0;
    old_cur_book      = (BOOK)NULL;
    old_cur_filter    = (GABSIGNAL *)NULL;
    cur_gabsignal     = (GABSIGNAL)NULL;
    gabsignal         = (GABSIGNAL)NULL;
    sigtmp            = (GABSIGNAL)NULL;
    temporary         = (GABSIGNAL)NULL;
    cur_norm          = (double *)NULL;
    pfG = pfCE1 = pfCE2 = pfCE3 = pfC = pfB = (double *)NULL;
    pnAep = pnIep = (int *)NULL;
    cur_MaxOctave     = 0;
    cur_MinOctave     = 0;
    for (i = 0; i < MAX_NUM_SB; i++) {
        filter[i]      = (GABSIGNAL *)NULL;
        transform[i]   = (GABSIGNAL *)NULL;
        library[i]     = (BOOK)NULL;
        num_filter[i]  = 0;
        TransAlloc[i]  = 0;
        filter_type[i] = -1;
    }
    for (i = 0; i < MAX_NUM_GABSIGNAL; i++)
        gabsignals[i] = (GABSIGNAL)NULL;
    gabord_initialised = 0;
}
