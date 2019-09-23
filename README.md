### GRiNCH: Graph-Regularized NMF and Clustering for Hi-C Data
GRiNCH applies non-negative matrix factorization (NMF) with graph regularization to discover clusters of highly interacting genomic regions from high-throughput chromosome conformation capture (Hi-C) data.

![alt text](http://pages.discovery.wisc.edu/~elee1/grinch_git/K562_chr12_5kb.png "GRINCH clusters discovered from Hi-C data, with CTCF binding enriched in cluster boundaries")

#### Usage
./grinch input.txt input.bed \[-o output_prefix\] 
```
>> ./grinch -h
>> ./grinch Huvec_chr12_25kb.txt Huvec_chr12_25kb.bed -o output/Huvec_chr12_25kb
```

#### Parameters
| Position | Parameter           | Description | Default Value/Behavior | 
| -     | -------------          | ----------------- | ---- |
| 1     | input matrix file      | Input matrix file path and name. Tab-delimited, in sparse matrix format, no header, e.g. 0 10 1201.78 | N/A | 
| 2     | input bed file         | Input bed file mapping each index in the input matrix file to a chromosomal coordinate. Tab-delimited, no header, e.g. chr1 50000  75000 0. Note: Bin size/resolution of the Hi-C data is assumed to be the same across all bins. Also, only cis-interactions are handled, i.e. the chromosome for all bins are assumed to be the same. | N/A | 
| optional | -o <output_file_prefix>    | Ouput file path and prefix. Note: will NOT create a directory if the specified directory does not exist. | 'output' | 
| optional | -k <number_of_clusters>  | Number of clusters, an integer value. |  n/(1000000/bin size) where n is the dimension of the symmetric input Hi-C matrix and bin size is the resolution in basepairs, i.e. k is set such hat the expected size of a cluster is 1Mb. | 
| optional | -e <expected_size_of_cluster>  | A different way to specify the number of clusters by the expected size of a cluster, i.e. if -e 500000, k = n/(500KB/bin size), where n is the number of bins in the input matrix. Note: -e will override -k. | 1000000, i.e. k = n / (1Mb/bin size) |  
| optional | -n <neighborhood_radius_size>  | Neighborhood radius used in regularization graph, in base pairs, and in multiples of resolution. -n 100000 would make neighborhood radius of 4 bins in 25kb resolution, and 4 adjacent regions on either side of a given regions will be used to regularize or 'smooth' the matrix factors. Increase for lower-depth or sparser data, to use more neighbors for smoothing. | 100000 | 
| optional | -l <lambda>  | Strength of regularization. | 1 | 
| optional | -s | Print the smoothed matrix to file,  in tab-delimited sparse matrix format (e.g. 25000 50000 500.3). This file can be large. | Do NOT output smoothed matrix. |
| optional | -f | Print to file the factor U and V. File may be large, since the matrix is written in a dense format, especially for higher-resolution input. | Do NOT output factor matrix. |
| optional | -g | Print to file the graph used in regularization. File may be large, since the matrix is written in a dense format.| Do NOT output graph. |

#### Input file format
* Matrix file, tab-delimited, sparse-matrix format, no header 
```
0	0	1000.2
0	1	1201.78
10	1	200.7
...
```
* Bed file, tab-delimited, no header 
```
chr1	50000	75000	0
chr1	75000	100000	1
chr1	100000	125000	2
...
```

#### Output files
* File suffixed '.tads' returns a list of putative TADs, each line with the first and last bin of each TAD (inclusive of last bin).
* File suffixed '.log' returns a plain-text file with list of parameter values used and time/memory consumption.
* Optional output file suffixed '.smoothed' returns the smoothed matrix in a tab-delimited sparse matrix format. This file may be large.
* Optional output file suffixed '.U' and '.V' returns the factors U and V respectively, in dense matrix format. File may be large, especially for higher-resolution input. Note that since the input cis-interaction Hi-C matrix is symmetric, U and V are equivalent up to some scaling factor and numerical error.
* Optional output file suffixed '.graph' returns the graph used in regularization. File may be large, since the matrix is written in a dense format.

#### Dependencies
[GSL (GNU Scientific Library)](https://www.gnu.org/software/gsl/doc/html/index.html) is used to handle matrix- and vector-related operations. The exact version used to compile the binary file (grinch) is included in this repo as a tarball (gsl.tgz). If the binary file fails to execute due to library dependency issues, try:
1. Download gsl.tgz from this repo and untar it (`tar -xzf gsl.tgz`).
2. Add the location of the untarred library to LD_LIBRARY_PATH. Assuming the library is in your current directory:
```
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:./gsl-sparse/lib/
```
3. Try executing the binary file again.

In order to implement NNDSVD initialization of factors, a fast randomized SVD algorithm from [RSVDPACK](https://github.com/sergeyvoronin/LowRankMatrixDecompositionCodes) was used. A minor modification to allow random seed specification was made to the original code from [RSVDPACK](https://github.com/sergeyvoronin/LowRankMatrixDecompositionCodes/tree/master/single_core_gsl_code). The code was compiled and included in this repo under modules/random_svd directory. For this dependency, no additional installation is required.

#### Installation
A make file is provided if you should want/need to recomplie the binary file grinch: `make grinch`. Compilers gcc and g++ are required.
