#include <gsl/gsl_matrix.h>
#include <list>
#include <vector>
#ifndef _cluster_
#define _cluster_
using namespace std;

typedef struct {
	int start;
	int end;
} Tad;

namespace cluster
{
	int findTads(vector<int>*, list<Tad>*);
	int assignMaxArgIndex(gsl_matrix*, vector<int>*);
	int assignClusterConstrainedKMedoids(gsl_matrix*, vector<int>*);
};
#endif
