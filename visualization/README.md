### Scripts for visualization

#### Requirements
The scripts were tested in Python 3.7 Anaconda environment. Here are the dependencies:
* [pyBigWig](https://github.com/deeptools/pyBigWig)
* matplotlib
* seaborn
* numpy
* pandas

#### Data
* The tarball of the data used to generate the figures below (also part of the repository as `*.png` files) can be downloaded from: <http://pages.discovery.wisc.edu/~elee1/grinch_git/data.tgz>. 
* Once untarred, if you place the data/ directory at the same level as the rest of the scripts, you can simply run the scripts with the commands noted below. 
* See the configuration files (`pyramid.ini`,`diamond.ini`) to find which input files are used to generate each figure.
* The directory also includes a README with data sources.

#### HiC heatmaps + GRiNCH clusters + genes + epigenetic signal tracks

The following image contains HiC heatmaps from a 4Mb region around zfp608 gene in chr18 of mouse embryonic stem cells (first column), neural precursor cells (second column), and cortical neurons (third column), along with GRiNCH clusters, and tracks for gene locations, H3k27ac ChIP-seq signals, and CTCF ChIP-seq signals. Note: for visualization purposes, the interaction counts in the HiC matrices have been log2-transformed.

![alt text](http://pages.discovery.wisc.edu/~elee1/grinch_git/zfp608.png "HiC heatmaps from a 4Mb region around zfp608 gene in chr18 of mouse embryonic stem cells, neural precursor cells, and cortical neurons, along with GRiNCH clusters, and tracks for gene locations, H3k27ac ChIP-seq signals, and CTCF ChIP-seq signals.")
To generate this image:
```
>> python pyramid.py pyramid.ini
```
and see the output file, `zfp608.png`. See `pyramid.ini` for list of input files used.

#### Pair of HiC heatmaps + GRiNCH clusters

The following image contains HiC heatmaps from a 4Mb region around zfp608 gene in chr18 of mouse embryonic stem cells (top) and cortical neurons (bottom), along with the GRiNCH clusters. Note: for visualization purposes, the interaction counts in the HiC matrices have been log2-transformed.

![alt text](http://pages.discovery.wisc.edu/~elee1/grinch_git/es_vs_cn.png "HiC heatmaps from a 4Mb region around zfp608 gene in chr18 of mouse embryonic stem cells and cortical neurons, along with GRiNCH clusters.")

To generate this image:
```
>> python diamond.py diamond.ini
```
and see the output file, `es_vs_cn.png`. See `diamond.ini` for list of input files used.
