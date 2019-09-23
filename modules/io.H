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
	int readMetaData(string, string &, int &, int &, int &, int &);
	int readHiCFile(string, gsl_matrix*, double &, double &);
	int writeTads(string, list<Tad>*, int, int); 
	int writeClusters(string, vector<int>, int, int); 
	int writeMatrixSparseFormat(string, gsl_matrix*);
	int writeMatrixDenseFormat(string, gsl_matrix*);
};
#endif
