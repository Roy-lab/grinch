""" This script contains the functions to draw heatmaps of HiC matrices
and also colorblocks of TADs.
"""
import seaborn as sns
import matplotlib.patches as patches
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

""" Ancillary function to help draw rotated half heatmap
@param n: row dimension of a symmetric matrix to be drawn as a heatmap
"""
def rotateCountMatrixUp(n):
	t = np.array([[1, 0.5],[-1, 0.5]])
	A = np.array([(y, x) for x in range(n, -1, -1) for y in range(n + 1)])
	A = np.dot(A, t)
	return A

""" Function to draw the upper triangle of a symmetric matrix
as a heatmap rotated 45 degrees as a upward pointing triangle
@param fig: overall figure object, needed to draw colorbar
@param ax: axis/subplot to draw the heatmap in
@param C: symmetric matrix to draw the heatmap for
"""
def drawRotatedHalfHeatmapUp(fig,ax, C):
	n = C.shape[0]
	A = rotateCountMatrixUp(n)
	ax.set_xlim([0,n])
	ax.set_ylim([-n/10,n])
	ax.xaxis.set_ticks(np.arange(0, n, n/4))
	ax.tick_params(which='both',direction='out',
				right=False,left=False,top=False,
				labelbottom=False,labelleft=False,
				length=4)
	cmap = plt.cm.get_cmap('Reds')
	pcm = ax.pcolormesh(A[:,1].reshape(n+1,n+1), A[:,0].reshape(n+1,n+1),
				np.flipud(C), 
				cmap=cmap)
	axins1 = inset_axes(ax, width="30%", height="5%",loc="upper right")
	cb = fig.colorbar(pcm, cax=axins1, orientation="horizontal")
	axins1.xaxis.set_ticks_position("bottom")
	cb.outline.set_visible(False)

""" Ancillary function to help draw rotated half heatmap pointing
downwards
@param n: row dimension of a symmetric matrix to be drawn as a heatmap
"""
def rotateCountMatrixDown(n):
	t = np.array([[-1, 0.5],[1, 0.5]])
	A = np.array([(y, x) for x in range(n, -1, -1) for y in range(n + 1)])
	A = np.dot(A, t)
	return A

""" Function to draw the lower triangle of a symmetric matrix
as a heatmap rotated 45 degrees as a downward pointing triangle
@param fig: overall figure object, needed to draw colorbar
@param ax: axis/subplot to draw the heatmap in
@param C: symmetric matrix to draw the heatmap for
"""
def drawRotatedHalfHeatmapDown(fig, ax, C):
	n = C.shape[0]
	A = rotateCountMatrixDown(n)
	ax.set_xlim([0,n])
	ax.set_ylim([-n,n/10])
	ax.xaxis.set_ticks(np.arange(0, n, n/4))
	ax.tick_params(which='both',direction='out',
				right=False,left=False,bottom=False,top=True,
				labelbottom=False,labelleft=False,labeltop=False,
				length=4)
	cmap = plt.cm.get_cmap('Reds') #'gist_heat_r')
	pcm = ax.pcolormesh(A[:,1].reshape(n+1,n+1), A[:,0].reshape(n+1,n+1),
				np.flipud(C), 
				cmap=cmap)
	axins1 = inset_axes(ax, width="30%", height="5%",loc="lower right")
	cb = fig.colorbar(pcm, cax=axins1, orientation="horizontal")
	axins1.xaxis.set_ticks_position("bottom")
	cb.outline.set_visible(False)

""" Function to draw the colorscale bar for a previously 
draw heatmap
@param ax: axis to draw the heatmap in
"""			
def drawHeatmapColorGradient(ax):
	ax.tick_params(which='both',
		right=False,left=False,top=False,bottom=False,
		labelbottom=False,labelleft=False, labelsize=6)
	cb = plt.colorbar()
	#cb.set_ticks([0.0, 1.8, 3.6, 5.4])

""" Function to draw the TADs as distinctly colored boxes
@param ax: axis to draw the TAD boxes in 
@param tads: list of TAD start and end indices, returend from readdata.getTads()
@param n: row dimension of the corresponding symmetric HiC matrix, i.e. number of bins
@param colormap: qualitative colormap to automatically assign colors to TADs
@param isPointingUpward: boolean; whether the associated heatmap is a upward-pointing triangle
"""
def drawTads(ax,tads,n,colormap, isPointingUpward):
	lowerLeftCorner = 0
	if isPointingUpward:
		lowerLeftCorner = -n/10
	tads.reset_index(inplace=True)
	norm = Normalize(vmin=0, vmax=tads.index.shape[0]-1)
	m = cm.ScalarMappable(norm = norm, cmap = colormap)
	for idx,row in tads.iterrows():
		start = max(row['start'],0)
		end = min(row['end'],n-1)
		length = end-start+1
		rect = patches.Rectangle((start,lowerLeftCorner),length,n/10,facecolor=m.to_rgba(idx),edgecolor=None)
		ax.add_patch(rect)

""" Function to draw a heatmap of a 1D array/vector,
which contains some signal or intensity level in the
location corresponding to each vector element/bin.
@param ax: axis to draw the heatmap in
@param vector: vector containing signal/intensity level
@param colormap: colormap to use for the heatmap
@param vmax: max value on the color scale
"""
def drawVectorHeatmap(ax,vector,colormap,vmax):
	n = vector.shape[0]
	vec = vector.reshape((1,n))
	ax.set_xlim((0,n))
	ax.tick_params(which='both',direction='out',
		right=False,left=False,top=False,
		labelbottom=False,labelleft=False,
		length=4)
	sns.heatmap(vec,vmax = vmax,cmap=colormap,ax=ax,cbar=False)
	ax.xaxis.set_ticks(np.arange(0, n, n/4))
