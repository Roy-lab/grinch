### Scripts for visualization

#### Requirements

#### Data
The tarball of the data used to generate the figures below (also part of the repository as *.png files) can be downloaded from: <http://pages.discovery.wisc.edu/~elee1/grinch_git/data.tgz>. Once 

#### HiC heatmaps + GRiNCH clusters + genes + epigenetic signal tracks

The following image contains HiC heatmaps from a 4Mb region around zfp608 gene in chr18 of mouse embryonic stem cells, neural precursor cells, and cortical neurons, along with GRiNCH clusters, and tracks for gene locations, H3k27ac ChIP-seq signals, and CTCF ChIP-seq signals.
![alt text](http://pages.discovery.wisc.edu/~elee1/grinch_git/zfp608.png "HiC heatmaps from a 4Mb region around zfp608 gene in chr18 of mouse embryonic stem cells, neural precursor cells, and cortical neurons, along with GRiNCH clusters, and tracks for gene locations, H3k27ac ChIP-seq signals, and CTCF ChIP-seq signals.")
To generate this image:
```
>> python pyramid.py pyramid.ini
```
and see the output file, zfp608.ini

