#include <gsl/gsl_matrix.h>
#ifndef _init_
#define _init_
using namespace std;

namespace init{
	int nndsvd(gsl_matrix*, gsl_matrix*, gsl_matrix*, int, double);
	int random(gsl_matrix*, gsl_matrix*, int);
};
#endif
