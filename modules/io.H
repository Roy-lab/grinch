#include <gsl/gsl_matrix.h>
#include <list>
#include <vector>
#include "cluster.H"
#ifndef _io_
#define _io_
using namespace std;
namespace io
{
	int printUsageToStdOut(string);
	int getMinAndMaxCoord(string, int &, int &);
	int readHiCFile(string, gsl_matrix*, int, int, double &);
	int outputTads(string, list<Tad>*, int, int); 
	int outputFactorSparseFormat(string, gsl_matrix*, int, int, vector<int>);
	int outputMatrixSparseFormat(string, gsl_matrix*, int, int, int, vector<int>);
	int outputMatrixDenseFormat(string, gsl_matrix*);
};
#endif
