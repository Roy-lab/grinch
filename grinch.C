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
#include <string>
#include <string.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include "modules/io.H"
#include "modules/initialization.H"
#include "modules/regularization.H"
#include "modules/cluster.H"
#include "modules/factorization.H"
#include "modules/preprocessing.H"
#include <sys/time.h>
#include <sys/resource.h>

int main(int argc, char **argv)
{
	struct timeval beginTime;
	gettimeofday(&beginTime,NULL);

	struct rusage bUsage;
	getrusage(RUSAGE_SELF,&bUsage);

	const char* inputFile;
	int binSize;
	
	int start = -1;
	int end = -1;
	
	int k = -1;
	int radius = 4;
	double lambda = 1;
	int threshold = 50;
	
	string outputPrefix = string("output");
	bool outputXSmoothed = false;
	bool outputFactors = false;
	bool outputGraph = false;
		
	int seed = 0; 
	string usage = string("usage.txt");

	int c;
	while((c = getopt(argc, argv, "o:b:e:k:r:l:t:mfsgh")) != -1)
		switch (c) {
			case 'o':
				outputPrefix = string(optarg);
				break;
			case 'b':
				start = atoi(optarg);
				break;
			case 'e':
				end = atoi(optarg);
				break;
			case 'k':
				k = atoi(optarg);
				break;
			case 'r':
				radius = atoi(optarg);
				break;
			case 'l':
				lambda = atof(optarg);
				break;
			case 't':
				threshold = atoi(optarg);
				break;
			case 's':
				outputXSmoothed = true;
				break;
			case 'f':
				outputFactors = true;
				break;
			case 'g':
				outputGraph =true;
				break;
			case 'h':
				io::printUsageToStdOut(usage);
				return 0;
			case '?':
				io::printUsageToStdOut(usage);
				return 1;
			default:
				io::printUsageToStdOut(usage);
				return 1;
		}	

	if ((argc - optind) < 2) {
		io::printUsageToStdOut(usage);
		return 1;
	} else {
		inputFile = argv[optind];
		binSize = atoi(argv[optind+1]);
	}

	ofstream logFile;
        logFile.open((outputPrefix+".log").c_str());
	logFile << "Hi-C resolution = " << binSize << " (basepairs per bin)" << endl;

	// Read input file		
	string inputFileName(inputFile);
	cout << "Reading interaction-counts file: " << inputFileName << endl;
	logFile << "Input file = " << inputFileName << endl;
	if (start < 0 || end < 0) {
		int b = -1;
		int e = -1;
		io::getMinAndMaxCoord(inputFileName, b, e);
		if (start < 0) {
			start = b;
		}
		if (end < 0) {
			end = e;
		}
	}
	logFile << "Starting coordinate = " << start << endl;
	logFile << "Ending coordinate = " << end << endl;

	int length = (end - start)/binSize + 1;
	if (k < 0) {
		k = length / (1000000/binSize);
	}
	logFile << "Lower dimension k = " << k << endl;
	gsl_matrix* X = gsl_matrix_calloc(length, length);
	double avg = 0.0;
	io::readHiCFile(inputFileName, X, start, binSize, avg);

	// remove sparse rows/columns
	vector<int> shifted = prep::removeEmptyRows(X, threshold);
	int n = shifted.size();
	logFile << "Sparsity threshold t = " << threshold << " (dropped " << length - n << " rows/columns with less than than " << threshold << " counts)" << endl;
	gsl_matrix_view hic = gsl_matrix_submatrix(X, 0, 0, n, n);

	// Initialize U & V
	cout << "Initializing factor matrices with NNDSVDa..." << endl;
	gsl_matrix* U = gsl_matrix_alloc(n, k);
	gsl_matrix* V = gsl_matrix_alloc(n, k);
	init::nndsvd(&hic.matrix, U, V, seed, avg);

	// Create weight matrix
	cout << "Creating neighborhood adjacency matrix for graph regularization..." << endl;
	gsl_matrix* W = gsl_matrix_alloc(n, n);
	gr::getWeightMatrixBinaryNeighbors(W, radius, shifted);
	logFile << "Strength of regularization lambda = " << lambda << endl;
	logFile << "Neighborhood radius r = " << radius << endl;

	// Graph-regularized NMF
	cout << "Factorizing input Hi-C matrix..." << endl;
	gr::nmf(&hic.matrix, U, V, W, lambda);
	gsl_matrix_free(X);

	// Optional output: smoothed count matrix before row-normalizing U
	if(outputXSmoothed) {
		gsl_matrix* X_smoothed = gsl_matrix_alloc(n,n);
		nmf::matrixCompletion(X_smoothed, U, V);	
		string xOutputFile = outputPrefix + ".smooth";
		cout << "Output smoothed matrix to file:\n\t " << xOutputFile << endl;
		io::outputMatrixSparseFormat(xOutputFile, X_smoothed, start, binSize, n, shifted); 
		gsl_matrix_free(X_smoothed);
	}

	// Clustering factors
	cout << "Clustering factor rows..." << endl;
	vector<int> clusters;
	cluster::assignMaxArgIndex(U, &clusters);

	// Shifting cluster assignments to account for filtered rows
	vector<int> l(length, -20);
	for (int i =0; i <n; i++) {
		int shift = shifted[i];
		l[i+shift] = clusters[i]; 
	}	

	// Find contiguous chunks of clusters (TADs)
	list<Tad> tl;
	cluster::findTads(&l, &tl);
	
	// Output contiguous chunks of clusters (TADs)
	string tadOutputFile = outputPrefix + ".tad";
	cout << "Output contingous cluster (putative TAD) coordinates file:\n\t" << tadOutputFile << endl;
	io::outputTads(tadOutputFile, &tl, start, binSize);

	// Optional output: weight matrix
	if (outputGraph) {
		string wOutputFile = outputPrefix + ".graph";
		cout << "Output weight matrix used in graph regularization to file:\n\t" << wOutputFile << endl;
		io::outputMatrixDenseFormat(wOutputFile, W); 
	}
	gsl_matrix_free(W);

	// Optional output: factor matrices
	if (outputFactors) {
		cout << "Output factors to file:\n\t" << outputPrefix << ".factor" << endl;	
		io::outputFactorSparseFormat(outputPrefix + ".factor", U, start, binSize, shifted);
	}

	gsl_matrix_free(U);
	gsl_matrix_free(V);

	struct timeval endTime;
	gettimeofday(&endTime,NULL);

	struct rusage eUsage;
	getrusage(RUSAGE_SELF,&eUsage);

	unsigned long int bt = beginTime.tv_sec;
	unsigned long int et = endTime.tv_sec;

	logFile << "Total time elapsed: " << et - bt << " seconds" << endl;

	unsigned long int bu = bUsage.ru_maxrss;
	unsigned long int eu = eUsage.ru_maxrss;
	
	logFile << "Memory usage: " << (eu - bu)/1000 << "MB" << endl;
	
	logFile.close();	
	return 0;
}
