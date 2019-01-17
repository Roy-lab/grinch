#include <gsl/gsl_matrix.h>
#include <list>
#include <vector>
#ifndef _gr_
#define _gr_
using namespace std;

namespace gr
{
	int nmf(gsl_matrix*, gsl_matrix*, gsl_matrix*, gsl_matrix*, double);
	int getWeightMatrixBinaryNeighbors(gsl_matrix*, int, vector<int>);
	int getColumnSumDiagonalMatrix(gsl_matrix*, gsl_matrix*);
};
#endif
