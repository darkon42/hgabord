VER = 1.0

CC     = gcc
CFLAGS = -O3 -march=native -std=gnu89 -ffast-math -fopenmp -I.

LIB_SRCS = gabord.c addons.c complex_op.c gb_bkold.c gb_buildbook.c \
           gb_corr.c gb_decomp.c gb_filter.c gb_oper.c math.c sig_alloc.c \
           sig_fct10.c struct_alloc.c update.c
LIB_OBJS = $(LIB_SRCS:.c=.o)

.PHONY: all clean

all: gabord libgabord.a

gabord: main.o $(LIB_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm -lfftw3

libgabord.a: $(LIB_OBJS)
	ar rcs $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	-rm -f gabord libgabord.a *.o
