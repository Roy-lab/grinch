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
#include <getopt.h>
#include "modules/io.H"
#include "modules/initialization.H"
#include "modules/regularization.H"
#include "modules/cluster.H"
#include "modules/factorization.H"
#include <sys/time.h>
#include <sys/resource.h>

int main(int argc, char **argv)
{
	struct timeval beginTime;
	gettimeofday(&beginTime,NULL);

	struct rusage bUsage;
	getrusage(RUSAGE_SELF,&bUsage);

	const char* matFile;
	const char* bedFile;
	int binSize;
	
	int start = -1;
	int end = -1;
	
	int k = -1;
	int expectedClusterSize = 1000000;
	int radiusBinNum = 4;
	int neighborhoodSize = 100000;
	double lambda = 1;
	
	string outputPrefix = string("output");
	bool outputXSmoothed = false;
	bool outputFactors = false;
	bool outputGraph = false;
	
	int seed = 0; 
	string usage = string("usage.txt");

	int c;
	while((c = getopt(argc, argv, "o:e:k:n:r:l:d:fsgh")) != -1)
		switch (c) {
			case 'o':
				outputPrefix = string(optarg);
				break;
			case 'e':
				expectedClusterSize = atoi(optarg);
				break;
			case 'k':
				k = atoi(optarg);
				break;
			case 'n':
				neighborhoodSize = atoi(optarg);
				break;
			case 'r':
				radiusBinNum = atoi(optarg);
				break;
			case 'l':
				lambda = atof(optarg);
				break;
			case 's':
				outputXSmoothed = true;
				break;
			case 'f':
				outputFactors = true;
				break;
			case 'g':
				outputGraph = true;
				break;
			case 'd':
				seed = atoi(optarg);
				break;
			case 'h':
				io::printUsageToStdOut(usage);
				return 0;
			case '?':
				io::printUsageToStdOut(usage);
				return 0;
			default:
				io::printUsageToStdOut(usage);
				return 1;
		}	

	if ((argc - optind) < 2) {
		io::printUsageToStdOut(usage);
		return 1;
	} else {
		matFile = argv[optind];
		bedFile = argv[optind+1];
	}

	ofstream logFile;
        logFile.open((outputPrefix+".log").c_str());

	// Read input files		
	cout << "Reading input files..." << endl;
	string bedFileName(bedFile);
	string matFileName(matFile);
	logFile << "HiC file = " << matFileName << endl;
	logFile << "Bed file = " << bedFileName << endl;
	string chro;
	int n;
	io::readMetaData(bedFileName, chro, binSize, n, start, end); 
	logFile << "Chromosome = " << chro << endl;
	logFile << "HiC resolution = " << binSize << " (basepairs per bin)" << endl;
	logFile << "Starting coordinate = " << start << endl;
	logFile << "Ending coordinate = " << end << endl;

	gsl_matrix* X = gsl_matrix_calloc(n, n);
	double avg, nonZeroCnt;
	io::readHiCFile(matFileName, X, avg, nonZeroCnt);
	
	// Set k
	k = n / (expectedClusterSize/binSize);
	logFile << "Lower dimension k = " << k << endl;
	
	// Initialize U & V
	cout << "Initializing factor matrices with NNDSVDa..." << endl;
	gsl_matrix* U = gsl_matrix_calloc(n, k);
	gsl_matrix* V = gsl_matrix_calloc(n, k);
	init::nndsvd(X, U, V, seed, avg);
	//init::random(U, V, seed);

	// Create weight matrix
	cout << "Creating neighborhood adjacency matrix for graph regularization..." << endl;
	gsl_matrix* W = gsl_matrix_calloc(n, n);
	radiusBinNum = neighborhoodSize / binSize;
	logFile << "Neighborhood radius = " << radiusBinNum << endl;
	logFile << "Strength of regularization lambda = " << lambda << endl;
	gr::getWeightMatrixBinaryNeighbors(W, radiusBinNum); 
	
	// Graph-regularized NMF
	cout << "Factorizing input Hi-C matrix..." << endl;
	gr::nmf(X, U, V, W, lambda);
	gsl_matrix_free(X);

	if(outputGraph) {
		io::writeMatrixDenseFormat(outputPrefix + ".graph", W);
	}
	gsl_matrix_free(W);
	
	cout << "Clustering..." << endl;	
	vector<int> clusters;
	cluster::assignClusterConstrainedKMedoids(U,&clusters);
	//io::writeClusters(outputPrefix + ".cluster", clusters, start, binSize);

	list<Tad> tads;
	cluster::findTads(&clusters, &tads);
	io::writeTads(outputPrefix + ".tads", &tads, start, binSize);
 
	// Optional output: smoothed count matrix
	if(outputXSmoothed) {
		gsl_matrix* X_smoothed = gsl_matrix_alloc(n,n);
		nmf::matrixCompletion(X_smoothed, U, V);	
		cout << "Output smoothed matrix to file:\n\t " << outputPrefix << ".smoothed" << endl;
		io::writeMatrixSparseFormat(outputPrefix + ".smoothed", X_smoothed); 
		gsl_matrix_free(X_smoothed);
	}
	
	// Optional output: factor matrices
	if (outputFactors) {
		cout << "Writing factors to file:\n\t" << outputPrefix << ".U/V" << endl;	
		io::writeMatrixDenseFormat(outputPrefix + ".U", U);
		io::writeMatrixDenseFormat(outputPrefix + ".V", V);
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
