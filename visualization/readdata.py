"""This script contains various functions to read in different
input files containing the data to visualize.
"""
import numpy as np
import pyBigWig as pbw
import pandas as pd

""" Determine the basepair coordinate that corresponds
to the first bin of a HiC matrix
@param bedFileNameList: list of bed files corresponding to
		HiC matrices of interest
"""
def readOffset(bedFileNameList):
	firstBins = []
	for bedFile in bedFileNameList:
		with open(bedFile) as f:
			firstBin = next(f).split('\t')[1]
			firstBin = int(firstBin)
			firstBins.append(firstBin)
	return firstBins

""" Read a HiC matrix file in sparse matrix
format into a symmetric numpy 2D array.
@param fname: name of the file containing the HiC data
@param b: starting index of the submatrix to visualize
@param e: ending index  of the submatrix to visualize (inclusive)
@return a symmetric numpy 2D array with the data filled in
"""
def readSymmetricSparseMatrix(fname,b,e):
	n=e-b+1
	x = np.zeros((n,n))
	f = open(fname)
	for line in f:
		row = line.split('\t')
		i = int(row[0])
		j = int(row[1])
		if i >= b and i <=e and j >=b and j <= e:
			cnt = float(row[2])
			x[i-b,j-b] = cnt
			x[j-b,i-b] = cnt	
	f.close()
	return x

""" Read a bigwig file into a 1D numpy array. Bins signals every
<binsize> basepairs by taking the max signal within the bin.
@param path: path to the file or an URL
@param chro: chromosome of interest to visualize
@param start: starting basepair coordinate of the region of interest
@param end: ending basepair coordinate of the region of interest
@param binsize: resolution of the HiC matrix to help bin signals
@return 1D numpy array of binned signals
"""
def readBigWig(path,chro,start,end,binsize):
	bw = pbw.open(path)
	binned = bw.stats(chro,start,end,nBins=int((end-start)/binsize),type="max")
	binned = np.array(binned,dtype=float)
	bw.close()
	return binned

""" Read a list of TAD files into dataframes. The start and end coordinates
of the TADs are converted to bin indices along the HiC submatrix that
will be visualized
@param listofTadFiles: list of file paths
@param start: starting basepair coordinate of the region to visualize
@param end: ending coordinate of the region to visualize
@param binsize: resolution of the HiC matrix
@return list of dataframes containing TAD start and end index
""" 
def getTads(listofTadFiles,start,end,binsize):
	tadlist = []
	for fname in listofTadFiles:
		tads = pd.read_csv(fname,sep='\t',header=None,names=['start','end'])
		tads = tads.loc[((tads['start']>=start) & (tads['start']<=end)) | ((tads['end']>=start) & (tads['end']<=end))].copy()
		tads['start'] = ((tads['start']-start)/binsize).astype(int)
		tads['end'] = ((tads['end']-start)/binsize).astype(int)
		tadlist.append(tads)
	return tadlist

""" Read a list of genes and its coordinates into a data frame.
The start and end coordinates of the genes are converted to bin indices along the
HiC submatrix that will be visualized.
@param fname: name of the file containing list of genes
@param start: starting basepair coordinate of the region to visualize
@param binsize: resolution of the HiC matrix, helps assign gene to specific bins
@return a dataframe containing gene start and end bin indices
"""
def processGeneFile(fname,start,binsize):
	genes = pd.read_csv(fname,sep='\t',usecols=[1,2,3],skiprows=1,header=None,names=['start','end','Gene_Symbol'])
	genes['start'] = ((genes['start']-start)/binsize).astype(int)
	genes['end'] = (((genes['end']-start)/binsize) + 1).astype(int)
	return genes

""" Read a file containing gene locations and their expression level 
into a data frame. The start and end coordinates of the genes are converted
to bin indices along the HiC submatrix that will be visualized
@param fname: name of the file containing list of genes
@param start: starting basepair coordinate of the region to visualize
@param n: number of bins in the HiC submatrix to be visualized (i.e. row dimension)
@param binsize: resolution of the HiC matrix, helps assign gene to specific bins
@return a vector containing total gene expresison level in each bin, and 
				a dataframe containing gene start and end bin indices
"""
def processGeneExpressionFile(fname,start,n,binsize):
	vec = np.zeros(n)
	genes = pd.read_csv(fname,sep='\t',usecols=['start','end','Gene_Symbol','tpm_mean'])
	genes['start'] = ((genes['start']-start)/binsize).astype(int)
	genes['end'] = (((genes['end']-start)/binsize)+1).astype(int)
	for idx, row in genes.iterrows():
		b = max(0,row['start'])
		e = min(n,row['end'])
		for i in range(b,e):
			vec[i] += row['tpm_mean']
	return vec,genes
