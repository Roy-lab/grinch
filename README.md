### GRiNCH: Graph-Regularized NMF and Clustering for Hi-C Data
GRiNCH applies graph-regularized non-negative matrix factorization to discover structural units of chromosomes from high-throughput chromosome conformation capture (Hi-C) data.

![alt text](http://pages.discovery.wisc.edu/~elee1/grinch_git/K562_chr12_5kb.png "GRINCH clusters discovered from Hi-C data, with CTCF binding enriched in cluster boundaries")

#### Usage
./grinch input_file resolution \[-o output_prefix\] \[-b begin_coordinate\] \[-e end_coordinate\] \[-k lower_dimension\] \[-l regularization_strength\] \[-r neighborhood_radius\] \[-t sparsity_threshold\] \[-fsg\]
```
>> ./grinch Huvec_chr12_25kb.txt 25000 -o Huvec_chr12
```

#### Parameters
| Order | Parameter              | Description | Default Value/Behavior | 
| -     | -------------          | ----------------- | ---- |
| 1     | input_file             | Tab-delimited input file path and name. Each line contains chromosome and start/end coordinate for loci 1 and loci 2, followed by the interaction count, e.g.,chr1 0 25000 chr1 25000 50000 8019.7 | N/A | 
| 2     | resolution             | Specify an integer value. The number of base-pairs in a single bin in the input Hi-C data, e.g., 5000, 10000, 25000. | N/A | 
| optional | -o output_prefix    | Can be any useful or descriptive label, used to prefix the output file names. | 'output' | 
| optional | -b begin_coordinate | Specify an integer value. If you want to find TADs in a subregion of the input file, specify the starting bin coordinate in base-pairs, e.g., -b 91000000. | smallest coordinate in input file | 
| optional | -e end_coordinate   | Specify an integer value. If you want to find TADs in a subregion of the input file, specify ending bin coordinate in base-pairs, e.g., -e 94000000. Inclusive of then bin starting with this last coordinate. | largest coodidnate in input file | 
| optional | -k lower_dimension  | Specify an integer value. The lower dimension of the factors; similar role to k in K-means clustering. You may get more than k number of TADs in the output, since a new TAD starts whenever the cluster assignment changes. |  n/(1000000/binSize) where n is the dimension of the symmetric input Hi-C matrix and bin_size is the resolution in basepairs | 
| optional | -l regularization_strenth  | Specify an integer or floating point value. Sets the strength of regularization. Increase for lower-depth or sparse data, to increase the influence of neighbors for smoothing. | 1 | 
| optional | -r neighborhood_radius    | Specify an integer value. If set to 2, 2 adjacent regions on either side of the given regions will be used to regularize or 'smooth' the matrix factors. Increase for lower-depth or sparser data, to use more neighbors for smoothing. | 4 | 
| optional | -t sparsity_threshold     | Specity an integer value. Minimum number of elements required in a row/column so that the row/column is included in the factorization. 0 means no sparse rows/columns are filtered out. 1 means completely empty rows/columns are filtered out. 50 means rows with fewer than 50 counts will be filtered out. | 50 |
| optional | -s | Print the smoothed matrix to file,  in tab-delimited sparse format (e.g. 25000 50000 500.3). This file can be large. | Do NOT output smoothed matrix. |
| optional | -f | Print to file the factor U, after row-normalization. First column in file is the genomic region, and the following k columns are the latent feature values from the factor, e.g., 50000 0.3 0.1 ..., tab-delimited. File may be large, especially for higher-resolution input. | Do NOT output factor matrix. |
| optional | -g | Print to file the graph used in regularization. File may be large, since the matrix is written in a dense format.| Do NOT output graph. |

#### Input file format
* Tab-delimited, no header: {region_1:chro} {region_1:start} {region_1:end} {region_2:chro} {region_2:start} {region_2:end}  {interaction_counts}
```
chr18	0	25000	chr18	0	25000	13472.4971021595
chr18	0	25000	chr18	25000	50000	277.3133210807
chr18	25000	50000	chr18	25000	50000	1056.3373179750
```

#### Output files
* File suffixed '.tad' returns a list of putative TADs, each line with the first and last bin of each TAD (inclusive of last bin).
* File suffixed '.log' returns a plain-text file with list of parameter values used and time/memory consumption.
* Optional output file suffixed '.smooth' returns the smoothed matrix in a tab-delimited sparse format, e.g., 25000 50000 500.3). This file may be large.
* Optional output file suffixed '.factor' returns the factor U, after row-normalization. First column in file is the genomic region, and the following k columns are the latent feature values from the factor, e.g., 50000 0.3 0.1 ..., tab-delimited. File may be large, especially for higher-resolution input.
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
