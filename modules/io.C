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

int io::getMinAndMaxCoord(string inputFile, int &start, int &end) {
	ifstream f(inputFile.c_str());
	int min = INT_MAX;
	int max = 0;
	while(f.good()) {
		int x,y;
		string _;
		f >> _ >> x >> _ >> _ >> y >> _ >> _;
		if (x < min) {
			min = x;
		}
		if (x > max) {
			max = x;
		}
		if (y < min) {
			min = y;
		}
		if (y > max) {
			max = y;
		}
	}
	f.close();
	start = min;
	end = max;
	return 0;
}

int io::readHiCFile(string inputFile, gsl_matrix* X, int start, int binSize, double &avg) {
	int dim = X->size1; // assuming symmetric matrix
	ifstream inFile(inputFile.c_str());
	avg = 0.0;
	while(inFile.good()) {
		int x, y;
		double cnt;
		string _;
		inFile >> _ >> x >> _ >> _ >> y >> _ >> cnt;
		x = (x-start)/binSize;
		y = (y-start)/binSize;
		if (x >= 0 && x < dim && y >= 0 && y < dim) {
			if (cnt < 0) {
				cnt = 0;
			}
			gsl_matrix_set(X, x, y, cnt);
			gsl_matrix_set(X, y, x, cnt);
			avg = avg + cnt;
		}
	}
	inFile.close();
	avg = avg/(dim * dim * 1.0);
	return 0;
}

int io::outputTads(string outputFile, list<Tad>* tl, int start, int binSize) {

	ofstream ot;
	ot.open(outputFile.c_str());
	ot << "First_Bin\tLast_Bin" << endl;

	list<Tad>::iterator itr;
	for (itr = tl->begin(); itr != tl->end(); ++itr) {
		int c = itr->cluster;
		if (c >= 0) {
			int b = itr->offset;
			int e = b + (itr->len) - 1;
			ot << b * binSize + start << "\t" << e * binSize + start  << endl;
		}
	}

	ot.close();
	return 0;
}

int io::outputFactorSparseFormat(string outputFile, gsl_matrix* F, int start, int binSize, vector<int> shifted) {
	ofstream ot;
	ot.open(outputFile.c_str());
	int rowNum = F->size1;
	int colNum = F->size2;
	for (int i = 0; i <rowNum; i++) {	
		int loc = i + shifted[i];
		ot << start + loc*binSize;
		for (int j = 0; j < colNum; j++) {
			ot << "\t" << gsl_matrix_get(F, i, j);		
		}
		ot << endl;
	}
	ot.close();
	return 0;
}

int io::outputMatrixSparseFormat(string outputFile, gsl_matrix* X, int start, int binSize, int maxDistance, vector<int> shifted) {
	ofstream ot;
	ot.open(outputFile.c_str());
	ot << "Bin1_Location\tBin2_Location\tInteraction_Count" << endl;
	
	int rowNum = X->size1;
	for (int i = 0; i < rowNum; i++) {
		for (int j = i; j < maxDistance; j++) {
			int b = i + shifted[i];
			int e = j + shifted[j];
			ot << start + b*binSize << "\t" << start + e*binSize << "\t" << gsl_matrix_get(X, i, j) << endl;
		}
	}
	ot.close();
	return 0;
}

int io::outputMatrixDenseFormat(string outputFile, gsl_matrix* X) {
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
