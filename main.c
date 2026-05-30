/*
 * main.c — Gabor Matching Pursuit test driver
 *
 * Reads raw double samples from data.bin in the current directory,
 * runs gabord() 100 times and prints the first 10 atoms of the first run.
 */

#include "gabord.h"

#include <stdio.h>
#include <stdlib.h>

#define DATAARRAY_SIZE 2048

int main(void)
{
    int    bsize, i, j;
    double threshold = 10.0;
    double *bookB;
    double *data;
    FILE   *fp;
    long   filesize;
    size_t nbdata, nread;

    data = (double *)malloc(DATAARRAY_SIZE * sizeof(double));
    if (!data) { fprintf(stderr, "alloc error\n"); return 1; }

    fp = fopen("data.bin", "rb");
    if (!fp) { perror("fopen"); free(data); return 1; }

    fseek(fp, 0, SEEK_END);
    filesize = ftell(fp);
    rewind(fp);
    nbdata = (size_t)filesize / sizeof(double);

    nread = fread(data, sizeof(double), nbdata, fp);
    fclose(fp);

    if (nread != nbdata) {
        fprintf(stderr, "Read error\n");
        free(data); return 1;
    }

    cth       = 0;
    ccoh      = 0;
    cpct      = 0;
    citer     = 500;
    thatomnrj = threshold;

    for (j = 0; j < 100; j++) {
        bookB = gabord(data, DATAARRAY_SIZE, &bsize);
        if (j == 0) {
            for (i = 1; i <= 10; i++)
                printf("Atom %d: Oct=%2.0f, ID=%2.0f, Pos.=%4.0f, "
                       "Freq.=%f, Phase=%1.3f\n",
                       i,
                       bookB[5*(i-1)+0], bookB[5*(i-1)+1],
                       bookB[5*(i-1)+2], bookB[5*(i-1)+3],
                       bookB[5*(i-1)+4]);
        }
        free(bookB);
    }

    printf("\nComputation completed!\n");
    free(data);
    gabord_cleanup();
    return 0;
}
