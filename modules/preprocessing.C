#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <vector>
#include "preprocessing.H"
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

std::vector<int> prep::removeEmptyRows(gsl_matrix* X, int threshold) {
	int n = X->size1;
	int cnt[n];
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (j == 0) {
				cnt[i] = 0;
			}
			if (gsl_matrix_get(X,i,j) > 0) {
				(cnt[i])++;
			} 
		}
	}
	
	int sparse[n];
	int prev = 0; 
	for (int i = 0; i <n; i++) {
		if ((cnt[i]) < threshold) {
			sparse[i] = prev + 1;			
		} else {
			sparse[i] = prev;
		}
		prev = sparse[i];
	}

	std::vector<int> shifted;
	for (int i = 0; i <n; i++) {
		if (cnt[i] >= threshold) {
			if (sparse[i] == 0) {
				shifted.push_back(0);
			} else if (sparse[i] > 0) {
				int shift = sparse[i];
				gsl_vector_view curr = gsl_matrix_column(X, i);
				gsl_matrix_set_col(X, i - shift, &curr.vector);
				shifted.push_back(shift);
			}
		}
	}

	for (int i = 0; i <n; i++) {
		if (cnt[i] >= threshold && sparse[i] > 0) {
			gsl_vector_view curr = gsl_matrix_row(X, i);
			gsl_matrix_set_row(X, i - sparse[i], &curr.vector);
		}
	}

	return shifted;
}
