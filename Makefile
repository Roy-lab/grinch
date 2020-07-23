#CHAGE PATHS AS NEEDED:
INCLUDE_PATH = ${CONDA_PREFIX}/include
LIBRARY_PATH = ${CONDA_PREFIX}/lib

#compiler: gcc for C programs, g++ for C++ programs
XX = g++
CC = gcc

#compiler flags
CFLAGS = -g
GSLFLAGS = -lgsl -lgslcblas

all: matf rsvd grinch

matf:
	$(CC) -c -o modules/random_svd/matrix_funcs.o modules/random_svd/matrix_vector_functions_gsl.c -I${INCLUDE_PATH}

rsvd:
	$(CC) -c -o modules/random_svd/rsvd.o modules/random_svd/low_rank_svd_algorithms_gsl.c -I${INCLUDE_PATH}

grinch:
	$(XX) grinch.C modules/*.C modules/random_svd/*.o -o grinch $(CFLAGS) -L${LIBRARY_PATH} ${GSLFLAGS} -I${INCLUDE_PATH} 

clean:
	rm grinch
	rm modules/random_svd/*.o
