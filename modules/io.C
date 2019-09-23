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
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <climits>
#include "io.H"
#include "cluster.H"

int io::printUsageToStdOut(string inputFile) {
	ifstream f(inputFile.c_str());
	string line;
	while(getline(f, line)) {
		cout << line << endl;
	}
	f.close();
	return 0;
}

int io::readMetaData(string inputFile, string &chro,int &binsize, int &n, int &start, int &end) {
	ifstream f(inputFile.c_str());
	int x,y,idx;
	string s;
	while(f >> s >> x >> y >> idx) {
		if (idx == 0) {
			chro = s;
			binsize = y - x;
			start = x;
		}
		n = idx;
		end = x;		
	}
	n = n + 1; // zero-indexed
	f.close();
	return 0;
}

int io::readHiCFile(string inputFile, gsl_matrix* X, double &avg, double &nonZeroCnt) { //, vector<int> &tally) {
	int n = X->size1; // assuming symmetric matrix
	int m = X->size2;
	ifstream inFile(inputFile.c_str());
	avg = 0.0;
	double diagonalCnt = 0;
	double offDiagonalCnt = 0;
	int x, y;
	double cnt;
	while(inFile >> x >> y >> cnt) {
		gsl_matrix_set(X, x, y, cnt);
		gsl_matrix_set(X, y, x, cnt);
		avg = avg + cnt;
		if (cnt > 0) {
			if (x == y) {
				diagonalCnt += 1;
			} else {
				offDiagonalCnt += 1;
			}
		}
	}
	inFile.close();
	avg = avg/(n * m * 1.0);
	nonZeroCnt = diagonalCnt + 2*offDiagonalCnt;
	return 0;
}

int io::writeTads(string outputFile, list<Tad>* tl, int start, int binSize) {

	ofstream ot;
	ot.open(outputFile.c_str());

	list<Tad>::iterator itr;
	for (itr = tl->begin(); itr != tl->end(); ++itr) {
		int b = itr->start;
		int e = itr->end;
		ot << b * binSize + start << "\t" << e * binSize + start  << endl;
	}

	ot.close();
	return 0;
}

int io::writeClusters(string outputFile, vector<int> clusters, int start, int binSize) {

	ofstream ot;
	ot.open(outputFile.c_str());

	int n = clusters.size();
	for (int i = 0; i < n; i++) {
		ot << i * binSize + start << "\t" << clusters[i]  << endl;
	}

	ot.close();
	return 0;
}

int io::writeMatrixSparseFormat(string outputFile, gsl_matrix* X) {
	ofstream ot;
	ot.open(outputFile.c_str());
	
	int rowNum = X->size1;
	int colNum = X->size2;
	for (int i = 0; i < rowNum; i++) {
		for (int j = 0; j < colNum; j++) {
			double val = gsl_matrix_get(X, i, j);
			if (val > 0) {
				ot << i << "\t" << j << "\t" << val << endl;
			}
		}
	}
	ot.close();
	return 0;
}

int io::writeMatrixDenseFormat(string outputFile, gsl_matrix* X) {
	ofstream ot;
	ot.open(outputFile.c_str());
	int rowNum = X-> size1;
	int colNum = X-> size2;
	for (int i = 0; i < rowNum; i++) {
		for (int j = 0; j < colNum; j++) {
			ot << gsl_matrix_get(X, i, j) << "\t";
		}
		ot << endl;
	}
	ot.close();
	return 0;
}
