""" This script takes in a config file containing 
parameter and input file list and draws a figure with
a HiC matrix heatmap, TADs, and a set of signal tracks
usage: python pyrmaid.py [config].ini
dependencies: heatmap.py, bar.py, readdata.py
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
import bar
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
	contexts = config.get('content','contexts').split(',')
	signals = config.get('content','signals').split(',')
	nCols = len(contexts)
	nSignals = len(signals)
	# hic file names
	hicfiles = []
	bedfiles = []
	for i in range(nCols):
		fname = config.get('hicfiles','file'+str(i))
		hicfiles.append(fname)
		bname = config.get('bedfiles','file'+str(i))
		bedfiles.append(bname)
	# get offset coordinates (i.e. starting basepair 
	# coordinate corresponding to index 0 of the HiC matrices)
	offsets = rd.readOffset(bedfiles)
	# track file names
	trackfiles = []
	for i in range(nSignals):
		row = []
		for j in range(nCols):
			fname = config.get(signals[i],'file'+str(j))
			row.append(fname)
		trackfiles.append(row)
	# tad file names
	tadfiles = []
	for i in range(nCols):
		fname = config.get('tadfiles','file'+str(i))
		tadfiles.append(fname)
	# get list of TADs	
	tadlist = rd.getTads(tadfiles,start,end,binsize)
	# gene list
	geneFile = config.get('files','genes')
	genes = rd.processGeneFile(geneFile,start,binsize)	
	# output file:
	output=config.get('files','output')
	
	#Drawing a grid to fit all plots
	fig = plt.figure(figsize = (nCols*4,2+nSignals))
	height_ratios = [3.5,0.1]
	for i in range(nSignals):
		height_ratios.append(0.5)
	gs = gridspec.GridSpec(figure=fig,ncols=nCols,nrows=nSignals+2,
		height_ratios=height_ratios)

	#Draw subplots to fill in
	heatmaps = []
	genetracks = []
	tracks = []
	for i in range(nCols):
		heatax = fig.add_subplot(gs[0,i],frame_on=False)
		heatmaps.append(heatax)
		geneax = fig.add_subplot(gs[1,i],frame_on=False)
		genetracks.append(geneax)
	for i in range(nSignals):
		row = []
		for j in range(nCols):
			trackax = fig.add_subplot(gs[i+2,j],frame_on=False)
			row.append(trackax)
		tracks.append(row)

	# draw HiC matrices
	for i in range(nCols):
		b = int((start-offsets[i])/binsize) # starting index of the region to visualize
		e = int((end-offsets[i])/binsize) # ending index
		C = rd.readSymmetricSparseMatrix(hicfiles[i],b,e)
		C = np.log2(C+1)
		ax = heatmaps[i]
		print('Drawing heatmap {:d}...'.format(i+1))
		heatmap.drawRotatedHalfHeatmapUp(fig,ax,C)
		heatmap.drawTads(ax,tadlist[i],e-b+1,'Set3',True)

	# draw gene location
	for i in range(nCols):
		print('Drawing gene track {:d}...'.format(i+1))
		bar.drawGenes(genetracks[i],genes,e-b+1)

	# draw rest of tracks
	colors=sns.husl_palette(nSignals+1,l=0.55)
	for i in range(nSignals):
		for j in range(nCols):
			print('Drawing {:} track {:d}...'.format(signals[i], j+1))
			vector = rd.readBigWig(trackfiles[i][j],chro,start,end,5000)
			bar.drawTrack(tracks[i][j],vector,colors[i+1])

	#X-axis 
	for i in range(nCols):
		ax = tracks[nSignals-1][i]
		ax.tick_params(labelbottom=True)
		tick_labels = range(start, end, int((end-start)/4))
		tick_labels = ['{0:,}'.format(x) for x in tick_labels]
		ax.set_xticklabels(tick_labels, ha = 'left',fontsize='medium')

	gs.tight_layout(fig)
	plt.savefig(output +'.png')

if __name__ == '__main__':
	main(sys.argv[1])
