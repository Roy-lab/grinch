#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort_vector.h>
#include <queue>
#include <math.h>
#include <iostream>
#include <float.h>
#include "kmedoids.H"
#include "io.H"

typedef pair<double, int> point;

Kmedoids::Kmedoids(int numClusters, gsl_matrix* data, int seed) {
	k = numClusters;
	U = data;
	seed = seed;
	maxIter = 100;
        rss = DBL_MAX;

	n = U->size1;
	medoidIndices = gsl_vector_alloc(k);
	clusterAssignment = gsl_vector_alloc(n);
	minDistance = gsl_vector_alloc(n);
	processed = gsl_vector_alloc(n);

	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

}

Kmedoids::~Kmedoids() {
	gsl_vector_free(medoidIndices);
	gsl_vector_free(clusterAssignment);
	gsl_vector_free(minDistance);
	gsl_vector_free(processed);
	gsl_rng_free(r);
}

int Kmedoids::setMaxInterations(int i) {
	maxIter = i;
	return 0;
}

gsl_vector* Kmedoids::getClusters(){
	return clusterAssignment;
}

double Kmedoids::getRss() {
	return rss;
}

double Kmedoids::getEuclideanDistance(int i, int j) {
	gsl_vector_view vecj = gsl_matrix_row(U, j);
	gsl_vector* diff = gsl_vector_alloc(U->size2);
	gsl_matrix_get_row(diff, U, i);
	gsl_vector_sub(diff, &vecj.vector);
	double dist = gsl_blas_dnrm2(diff);
	gsl_vector_free(diff);
	return dist;
}

int Kmedoids::initializeMedoidsMax() {
	for (int i = 0; i < k; i++) {
		gsl_vector_view col = gsl_matrix_column(U, i);
		int maxIdx = gsl_vector_max_index(&col.vector);
		gsl_vector_set(medoidIndices, i, maxIdx);
	}
	return 0;
}

int Kmedoids::initializeMedoidsRandom(int seed) {

	double nums[n], picks[k];
	for (int i = 0; i < n; i++) {
		nums[i] = (double) i;
	}
	gsl_ran_choose(r, picks, k, nums, n, sizeof(double)); 

	for (int i = 0; i < k; i++) {
		int medoid = picks[i];
		gsl_vector_set(medoidIndices,i,medoid);
	}
	gsl_sort_vector(medoidIndices);
	return 0;
}
 
int Kmedoids::getClosestReachableMedoids(int idx, int &upstream, int &downstream) {
	upstream = -1;
	int c = 1;
	while (upstream <0 && (idx - c) >=0) {
		upstream = gsl_vector_get(clusterAssignment, idx - c);
		c += 1;
	}
	downstream = -1; //gsl_vector_get(clusterAssignment, idx);
	c = 1;
	while (downstream <0 && (idx + c) < n) {
		downstream = gsl_vector_get(clusterAssignment, idx +c);
		c += 1;
	}
	return 0;
}

int Kmedoids::updateMedoids() {
	int clusterId = gsl_vector_get(clusterAssignment, 0);
	int clusterStart = 0;
	int clusterCnt = 0;
	for (int i = 1; i < n; i++) {
		int newClusterId = (int)gsl_vector_get(clusterAssignment, i);
		if (clusterId != newClusterId) {
			int medoid = getMedoid(clusterStart, i);
			gsl_vector_set(medoidIndices, clusterCnt, medoid);			
			clusterCnt += 1;
			clusterId = newClusterId;
			clusterStart = i;	
		}
	}
	int medoid = getMedoid(clusterStart, n);
	gsl_vector_set(medoidIndices, clusterCnt, medoid);
	return 0;
}
int Kmedoids::getMedoid(int clusterStart, int i){
	gsl_matrix* distance = gsl_matrix_calloc(i-clusterStart, i-clusterStart);
	for (int a = clusterStart; a < i; a++) {
		for (int b = a+1; b <i; b++) {
			double dist = getEuclideanDistance(a, b);
			gsl_matrix_set(distance, a-clusterStart, b-clusterStart, dist);
			gsl_matrix_set(distance, b-clusterStart, a-clusterStart, dist);
		}
	}
	gsl_vector* rowsum = gsl_vector_alloc(i-clusterStart);
	for (int a = 0; a< (i-clusterStart); a++) {
		gsl_vector_view row = gsl_matrix_row(distance, a);
		double rsum = gsl_blas_dasum(&row.vector);
		gsl_vector_set(rowsum, a, rsum);		
	}
	int medoid = gsl_vector_min_index(rowsum);
	medoid = clusterStart+medoid;
	gsl_vector_free(rowsum);
	gsl_matrix_free(distance);
	return medoid;
}
int Kmedoids::cluster(){
	//rss = DBL_MAX;
	initializeMedoidsMax();
	//initializeMedoidsRandom(seed);
	for (int numIter = 0; numIter < maxIter; numIter++) {
		gsl_vector_set_all(clusterAssignment, -1);
		gsl_vector_set_zero(processed);
		priority_queue<point, vector<point>, greater<point> > pq;
		//priority_queue<pair<double,int>, vector< pair<double,int> >, std::greater< pair<double, int> > pq;
		for (int i=0; i <k; i++) {
			int idx = (int) gsl_vector_get(medoidIndices, i);
			pq.push(make_pair(0, idx));
			gsl_vector_set(clusterAssignment, idx, idx);
			gsl_vector_set(minDistance, idx, 0);
			gsl_vector_set(processed, idx, 1);
			/*
			for (int i = 0; i < n; i++) {
				int c = (int) gsl_vector_get(clusterAssignment, i);
				if (c < 0) {
					cout<< '+';
				} else {
					cout<< c;
				}
			}
			cout << endl;
			*/
		}	
		while (!pq.empty()) {
			pair<double, int> v = pq.top();
			pq.pop();
			double priority = v.first;	
			int idx = v.second;
			if (gsl_vector_get(clusterAssignment, idx) < 0) {
				int upstream, downstream;
				getClosestReachableMedoids(idx, upstream, downstream);
				int myMedoid;
				double myMinDist;
				if (upstream == -1) {
					myMedoid = downstream;
					myMinDist = getEuclideanDistance(idx, downstream);
				} else {
					myMedoid = upstream;
					myMinDist = getEuclideanDistance(idx, upstream);
					if (downstream > -1) {
						double downDist = getEuclideanDistance(idx, downstream);
						if (downDist < myMinDist) {
							myMedoid = downstream;
							myMinDist = downDist;
						}
					}
				}
				gsl_vector_set(clusterAssignment, idx, myMedoid);
				gsl_vector_set(minDistance, idx, myMinDist);
				/*
				for (int i = 0; i < n; i++) {
					int c = (int) gsl_vector_get(clusterAssignment, i);
					if (c < 0) {
						cout<< '+';
					} else {
						cout<< c;
					}
				}
				cout << endl;
				*/
			}
			if (((idx-1) >=0) && gsl_vector_get(processed, idx-1) == 0)  {
				pq.push(make_pair(priority+getEuclideanDistance(idx,idx-1), idx-1));
				gsl_vector_set(processed,idx-1,1);			
			}
			if (((idx+1) < n) && gsl_vector_get(processed, idx+1) == 0) {
				pq.push(make_pair(priority+getEuclideanDistance(idx,idx+1), idx+1));
				gsl_vector_set(processed,idx+1,1);			
			}
		
				
		}		
		double newRss = gsl_blas_dasum(minDistance);
		cout << "\tRSE = " << newRss << endl;
		if (rss <= newRss) {
			break;
		} else {
			rss = newRss;
			updateMedoids();
		}	
	}
	return 0;
}
