#include <iostream>
#include <fstream>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <list>
#include <vector>
#include <math.h>
#include <queue>
#include <string>
#include <string.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include "cluster.H"
#include "kmedoids.H"

int cluster::assignMaxArgIndex(gsl_matrix* U, vector<int>* clusterAssignment) {
	int n = U->size1;
	int k = U->size2;
	int numAmbiguousAssignments = 0;
	for (int i = 0; i < n; i++) {
		gsl_vector_view region = gsl_matrix_row(U, i);
		// row-normalize
		double sum = gsl_blas_dasum(&region.vector);
		gsl_vector_scale(&region.vector, 1/sum);
		// argmax		
		int c = static_cast<int>(gsl_vector_max_index(&region.vector));
		clusterAssignment->push_back(c);
	}
	return numAmbiguousAssignments;
}

int cluster::assignClusterConstrainedKMedoids(gsl_matrix* U, vector<int>* clusterAssignment){
        int n = U->size1;
        int k = U->size2;
        Kmedoids clustering = Kmedoids(k, U, 0);
        clustering.cluster();
        gsl_vector* clusters = clustering.getClusters();
	for (int i = 0; i <n; i++) {
		clusterAssignment->push_back((int)gsl_vector_get(clusters,i));
	}
	return 0;
}
int cluster::findTads(vector<int>* clusters, list<Tad>* tl) {

	int n = clusters->size();
	int currCluster = clusters->front();
	int currOffset = 0;
	int i = 0;
	vector<int>::const_iterator itr;
	for (itr = clusters->begin(); itr != clusters->end(); ++itr) {
		int c = *itr;
		if (c != currCluster) {
			Tad t;
			t.start = currOffset;
			t.end = i-1;
			tl->push_back(t);
			currCluster = c;
			currOffset = i;
		}
		i = i+1;
	}
	
	Tad t;
	t.start= currOffset;
	t.end = n-1;
	tl->push_back(t);
	return 0;
}

