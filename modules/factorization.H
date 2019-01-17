#include <gsl/gsl_matrix.h>
#ifndef _nmf_
#define _nmf_
using namespace std;

namespace nmf
{
	int multiplicativeUpdate(gsl_matrix*, gsl_matrix*, gsl_matrix*);
	int matrixCompletion(gsl_matrix*, gsl_matrix*, gsl_matrix*);
};
#endif
