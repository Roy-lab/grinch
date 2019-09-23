#include <iostream>
#include <fstream>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <list>
#include <math.h>
#include <queue>
#include <string>
#include <string.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include "initialization.H"

extern "C" {
	#include "random_svd/low_rank_svd_algorithms_gsl.h"
}

int init::nndsvd(gsl_matrix* A, gsl_matrix* W, gsl_matrix* H, int seed, double a) {
	int m = A->size1;
	int n = A->size2;
	int k = W->size2;

	gsl_matrix *U, *S, *V;
	
	randomized_low_rank_svd3(A, k, 3, 1, &U, &S, &V, seed);
	
	gsl_vector_view x = gsl_matrix_column(U, 0);
	gsl_vector_view y = gsl_matrix_column(V, 0);
	double sigma = sqrt(gsl_matrix_get(S, 0, 0));
	gsl_vector_scale(&x.vector, sigma);
	gsl_vector_scale(&y.vector, sigma);
	gsl_matrix_set_col(W, 0, &x.vector);
	gsl_matrix_set_col(H, 0, &y.vector);

	for (int i = 1; i < k; i++) {
		gsl_vector* xp = gsl_vector_calloc(m);
		gsl_vector* xn = gsl_vector_calloc(m);
		gsl_vector* yp = gsl_vector_calloc(n);
		gsl_vector* yn = gsl_vector_calloc(n);

		for (int j = 0; j < m; j++) {
			double val = gsl_matrix_get(U, j, i);
			if (val >= 0) {
				gsl_vector_set(xp, j, val);
			} else {	
				gsl_vector_set(xn, j, val*-1);
			}
		}		

		for (int j = 0; j < n; j++) {
			double val = gsl_matrix_get(V, j, i);
			if (val >= 0) {
				gsl_vector_set(yp, j, val);
			} else {	
				gsl_vector_set(yn, j, val*-1);
			}
		}		
	
		double xpNorm = gsl_blas_dnrm2(xp);
		double xnNorm = gsl_blas_dnrm2(xn);
		double ypNorm = gsl_blas_dnrm2(yp);
		double ynNorm = gsl_blas_dnrm2(yn);

		double mp = xpNorm * ypNorm;
		double mn = xnNorm * ynNorm;


		if (mp > mn) {
			gsl_vector_scale(xp, 1/xpNorm * sqrt(gsl_matrix_get(S,i,i)*mp));
			gsl_vector_scale(yp, 1/ypNorm * sqrt(gsl_matrix_get(S,i,i)*mp));	
			gsl_matrix_set_col(W, i, xp);
			gsl_matrix_set_col(H, i, yp);	
		} else {
			gsl_vector_scale(xn, 1/xnNorm * sqrt(gsl_matrix_get(S,i,i)*mn));
			gsl_vector_scale(yn, 1/ynNorm * sqrt(gsl_matrix_get(S,i,i)*mn));
			gsl_matrix_set_col(W, i, xn);
			gsl_matrix_set_col(H, i, yn);
		}

		gsl_vector_free(xp);
		gsl_vector_free(xn);
		gsl_vector_free(yp);
		gsl_vector_free(yn);
	}
	
	gsl_matrix_free(U);
	gsl_matrix_free(V);
	gsl_matrix_free(S);	
	
	for (int i=0; i < m; i++) {
		for (int j=0; j < k; j++) {
			double val = gsl_matrix_get(W, i, j);
			if (val <= 0) {
				gsl_matrix_set(W, i, j, a);
			}			
		}
	}

	for (int i=0; i < n; i++) {
		for (int j=0; j < k; j++) {
			double val = gsl_matrix_get(H, i, j);
			if (val <= 0) {
				gsl_matrix_set(H, i, j, a);
			}			
		}
	}

	return 0;
}

int init::random(gsl_matrix* W, gsl_matrix* H, int seed) {
	int rowCnt = W->size1;
	int k = W->size2;
	int colCnt = H->size1;

	const gsl_rng_type* T;
	gsl_rng* ri;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	ri = gsl_rng_alloc(T);
	gsl_rng_set(ri, seed);

	for (int i = 0; i < rowCnt; i++) {
		for (int j = 0; j < k; j++) {
			gsl_matrix_set(W, i, j, gsl_rng_uniform(ri));
		}
	}

	for (int i = 0; i < colCnt; i++) {
		for (int j = 0; j < k; j++) {
			gsl_matrix_set(H, i, j, gsl_rng_uniform(ri));
		}
	}

	gsl_rng_free(ri);
	return 0;
}

