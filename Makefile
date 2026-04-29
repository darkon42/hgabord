VER=1.0;

CC=gcc

SRC=addons.c complex_op.c gb_bkold.c gb_buildbook.c gb_corr.c gb_decomp.c gb_filter.c gb_oper.c math.c sig_alloc.c sig_fct10.c struct_alloc.c update.c

gabord:
	${CC} -O3 -march=native -std=gnu89 -ffast-math -fopenmp -o gabord gabord.c ${SRC} -lm -lfftw3

clean:
	-rm gabord
	


