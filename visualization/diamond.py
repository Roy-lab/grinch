""" This script takes in a config file containing 
parameter and input file list and draws a figure with
a pair of HiC heatmaps and the corresponding TADs/clusters
from the same region. Note for visualiation purposes the
counts in the HiC files are log2-transformed.
usage: python diamond.py [config].ini
dependencies: heatmap.py, readdata.py
"""
import matplotlib.patches as patches
import sys
import csv
import numpy as np
from configparser import ConfigParser
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import gridspec
from matplotlib.colors import Normalize
import heatmap
import readdata as rd
import seaborn as sns

def main(config_file):
	# Read in parameters from config file:
	config = ConfigParser()
	config.read(config_file)
	chro = config.get('coordinates','chro')
	start = config.getint('coordinates','start')
	end = config.getint('coordinates', 'end')
	binsize = config.getint('coordinates','binsize')
	nRows = 2
	# hic file names
	hicfiles = []
	bedfiles = []
	for i in range(nRows):
		fname = config.get('hicfiles','file'+str(i))
		hicfiles.append(fname)
		bname = config.get('bedfiles','file'+str(i))
		bedfiles.append(bname)
	# tad file names
	tadfiles = []
	for i in range(nRows):
		fname = config.get('tadfiles','file'+str(i))
		tadfiles.append(fname)
	# get list of TADs	
	tadlist = rd.getTads(tadfiles,start,end,binsize)
	# get offset coordinates (i.e. starting basepair 
	# coordinate corresponding to index 0 of the HiC matrices)
	offsets = rd.readOffset(bedfiles)
	# output file:
	output=config.get('files','output')
	
	#Drawing a grid to fit all plots
	fig = plt.figure(figsize = (5,nRows*3))
	gs = gridspec.GridSpec(figure=fig,ncols=1,nrows=nRows)

	#Draw subplots to fill in
	heatmaps = []
	for i in range(nRows):
		heatax = fig.add_subplot(gs[i,0],frame_on=False)
		heatmaps.append(heatax)

	# draw HiC matrices
	for i in range(nRows):
		b = int((start-offsets[i])/binsize) # starting index of the region to visualize
		e = int((end-offsets[i])/binsize) # ending index
		C = rd.readSymmetricSparseMatrix(hicfiles[i],b,e)
		C = np.log2(C+1)
		ax = heatmaps[i]
		print('Drawing heatmap {:d}...'.format(i+1))
		if i == 0:
			heatmap.drawRotatedHalfHeatmapUp(fig, ax,C)
			heatmap.drawTads(ax,tadlist[i],e-b+1,'Set3',True)
			# X-axis formatting
			ax.tick_params(labelbottom=True)
			tick_labels = range(start, end, int((end-start)/4))
			tick_labels = ['{0:,}'.format(x) for x in tick_labels]
			ax.set_xticklabels(tick_labels, ha = 'left',fontsize='small')
		else:
			heatmap.drawRotatedHalfHeatmapDown(fig, ax,C)
			heatmap.drawTads(ax,tadlist[i],e-b+1,'Set3',False)
			# X-axis formatting
			ax.set_xticklabels([])
		
	gs.tight_layout(fig)
	plt.savefig(output +'.png')

if __name__ == '__main__':
	main(sys.argv[1])
