### Visualization Tutorial

#### 1. Install dependencies
The scripts were tested in Python 3.7 Anaconda environment. Here are the dependencies, all of which can be installed with conda.
* [pyBigWig](https://github.com/deeptools/pyBigWig)
* matplotlib
* seaborn
* numpy
* pandas

#### 2. Download example input data
1. Download the tarball of the dataset used to generate the example figures: <http://pages.discovery.wisc.edu/~elee1/grinch_git/data.tgz>. 
2. Place the `data/` directory at the same level as the rest of the visualization scripts (then you can simply run the scripts with the commands noted below). 

#### 3. Check out the input data format
1. See the configuration files (`pyramid.ini`,`diamond.ini`) to find which input files are used to generate each figure. The comments in the configuration files (hopefully) will explain which file is used to generated which portion of the figure.
2. Check out the README in the downloaded dataset; it includes all the data sources. The bit about downloading a list of genes and their coordinates from UCSC table browser might be particularly useful if you want to visualize your own genome/chromosomal stretch.

#### 4.1 Run example 1: Hi-C heatmaps + GRiNCH clusters + genes + epigenetic signal tracks

`python pyramid.py pyramid.ini` will output `zfp608.png`, shown below. The figure contains HiC heatmaps from a 4Mb region around zfp608 gene in chr18 of mouse embryonic stem cells (first column), neural precursor cells (second column), and cortical neurons (third column), along with GRiNCH clusters, and tracks for gene locations, H3k27ac ChIP-seq signals, and CTCF ChIP-seq signals. Note: for visualization purposes, the interaction counts in the HiC matrices have been log2-transformed.

![alt text](http://pages.discovery.wisc.edu/~elee1/grinch_git/zfp608.png "HiC heatmaps from a 4Mb region around zfp608 gene in chr18 of mouse embryonic stem cells, neural precursor cells, and cortical neurons, along with GRiNCH clusters, and tracks for gene locations, H3k27ac ChIP-seq signals, and CTCF ChIP-seq signals.")

#### 4.2 Run example 2: pair of HiC heatmaps + GRiNCH clusters
*  `python diamond.py diamond.ini` will output `es_vs_cn.png`, shown below. The figure contains HiC heatmaps from a 4Mb region around zfp608 gene in chr18 of mouse embryonic stem cells (top) and cortical neurons (bottom), along with the GRiNCH clusters. Note: for visualization purposes, the interaction counts in the HiC matrices have been log2-transformed.

![alt text](http://pages.discovery.wisc.edu/~elee1/grinch_git/es_vs_cn.png "HiC heatmaps from a 4Mb region around zfp608 gene in chr18 of mouse embryonic stem cells and cortical neurons, along with GRiNCH clusters.")
