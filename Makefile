#compiler: gcc for C programs, g++ for C++ programs
CC = g++

#compiler flags
CFLAGS = -g
GSLFLAGS = -lgsl -lgslcblas

all: clean grinch

matf:
	gcc  -c -o modules/random_svd/matrix_funcs.o modules/random_svd/matrix_vector_functions_gsl.c

grinch:
	$(CC) grinch.C modules/*.C modules/random_svd/*.o -o grinch $(GSLFLAGS) $(CFLAGS)

clean:
	rm grinch
