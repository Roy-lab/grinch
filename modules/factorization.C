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
#include <sys/resource.h>
#include <sys/time.h>
#include "factorization.H"

int nmf::multiplicativeUpdate(gsl_matrix* X, gsl_matrix* U, gsl_matrix* V) {

	int rowCnt = X->size1; // # of rows in the matrix
	int colCnt = X->size2; // # of columns in the matrix, should be same as rowCnt
	int k = U->size2;	
	// account for all-zero rows in real data; offset by tiny bit 
	gsl_vector_view firstRowView = gsl_matrix_row(X, 0); 
	gsl_vector_view firstColView = gsl_matrix_column(X, 0);
	gsl_vector_add_constant(&firstRowView.vector, 0.1);
	gsl_vector_add_constant(&firstColView.vector, 0.1);

	int cycle = 0;
	double prevError = 10.0;
	int patience = 5;
	int patienceCnt = 0;

	// Matrix fatorization!!
	gsl_matrix* denom = gsl_matrix_alloc(rowCnt, k);
	gsl_matrix* numer_ = gsl_matrix_alloc(rowCnt, colCnt);
	gsl_matrix* numer  = gsl_matrix_alloc(rowCnt, k);

	gsl_matrix* denom2 = gsl_matrix_alloc(colCnt, k);
	gsl_matrix* numer_2 = gsl_matrix_alloc(colCnt, rowCnt);
	gsl_matrix* numer2 = gsl_matrix_alloc(colCnt, k);

	gsl_matrix* X_guess = gsl_matrix_alloc(rowCnt, colCnt);
	gsl_matrix* residuals = gsl_matrix_alloc(rowCnt, colCnt);

	while (patienceCnt < patience) {

		gsl_matrix_set_zero(denom);
		gsl_matrix_set_zero(numer_);
		gsl_matrix_set_zero(numer);

		gsl_matrix_set_zero(denom2);
		gsl_matrix_set_zero(numer_2);
		gsl_matrix_set_zero(numer2);

		gsl_matrix_set_zero(X_guess);
		gsl_matrix_set_zero(residuals);

		//U_new = U(XV)(UV'V)^(-1)
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, X, V, 0.0, denom); // XV
			
		gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, U, V, 0.0, numer_); // UV'
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, numer_, V, 0.0, numer); // UV'V
		
		gsl_matrix_div_elements(denom, numer);
		gsl_matrix_mul_elements(U, denom);

		// V_new = V(X'V)(VU_new'U_new)^(-1)
		gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, X, U, 0.0, denom2);// X'U
		
		gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, V, U, 0.0, numer_2); // VU'
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, numer_2, U, 0.0, numer2); // VU'U

		gsl_matrix_div_elements(denom2, numer2);
		gsl_matrix_mul_elements(V, denom2);
		
		// X_guess = UV'
		gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, U, V, 0.0, X_guess);
		
		// X - X_guess & mean squared error
		gsl_matrix_sub(X_guess, X);
		gsl_matrix_memcpy(residuals, X_guess);
		gsl_matrix_mul_elements(residuals, X_guess);
		double squared = 0.0;
		gsl_vector* temp = gsl_vector_alloc(colCnt);
		for (int i = 0; i < rowCnt; i++) {
			gsl_matrix_get_row(temp, residuals, i);
			squared = squared + gsl_blas_dasum(temp);	
		}
		gsl_vector_free(temp);
		squared = squared / (rowCnt * colCnt);
		cout << "\tCycle " << cycle << ": mean squared error =  " << squared << endl;
		if ((prevError - squared) / prevError < 0.001) {
			patienceCnt ++;
		}  

		prevError = squared;
		cycle ++;
		if (cycle > 100) {
			break;
		}
	}
	gsl_matrix_free(denom); 
	gsl_matrix_free(numer_);	 
	gsl_matrix_free(numer);
	
	gsl_matrix_free(denom2);
	gsl_matrix_free(numer_2);
	gsl_matrix_free(numer2);

	gsl_matrix_free(X_guess);
	gsl_matrix_free(residuals);

	gsl_vector_add_constant(&firstRowView.vector, -0.1);
	gsl_vector_add_constant(&firstColView.vector, -0.1);

	return 0;
}

int nmf::matrixCompletion(gsl_matrix* X_approximate, gsl_matrix* U, gsl_matrix* V) {
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, U, V, 0.0, X_approximate);
	return 0;
}
