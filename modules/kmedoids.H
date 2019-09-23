#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#ifndef _kmedoids_
#define _kmedoids_
using namespace std;

class Kmedoids {
	public: 
		Kmedoids(int, gsl_matrix*,int);
		~Kmedoids();

		int setMaxInterations(int);
		int cluster();
		double getRss();
		gsl_vector* getClusters();

	private:
		gsl_matrix* U;
		gsl_vector* medoidIndices;
		gsl_vector* clusterAssignment;
		gsl_vector* minDistance;
		gsl_vector* processed;
		int k;
		int n;
		int seed;
		double rss;
		int maxIter;
		const gsl_rng_type * T;
		gsl_rng * r;

		int initializeMedoidsMax();
		int initializeMedoidsRandom(int);
		int getClosestReachableMedoids(int, int&, int&);
		int updateMedoids();
		int getMedoid(int, int);
		double getEuclideanDistance(int i, int j);
};
#endif
