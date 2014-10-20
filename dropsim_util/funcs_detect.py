#!/usr/bin/env python

import numpy as np
import os
from pygoods import *
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import interpolate

class detfunc():
	"""
	One can calculate the detection completeness in any kind of bins (or at least this is 
	the goal)!
	"""
	def __init__(self, catalog):
		"""
		Initialize the catalog and color criteria.
		"""
		self.c = Ftable(catalog)
		if catalog.split('.')[-1] == 'fits':
			# if it's a FITS table
			self.cat_form = 'fits'
		else:
			self.cat_form = 'sext'
		self.bins1 = None; self.bins2 = None
		self.col1 = None; self.col2 = None
		
	def set_cols(self, colname1, array1, colname2="", array2=[]):
		"""
		Set the columns, and values in those columns, that one uses to bin the data.
		set_cols(self, colname1, array1, colname2="", array2=[])
		"""
		self.col1 = colname1
		self.col2 = colname2
		self.array1 = array1
		self.array2 = array2
		if len(array2):
			self.ndim = 2
		else: self.ndim = 1
		
	def calc_detection1d(self, bins, othercrit=None, axis=1):
		"""
		calc_detection1d(self, bins, othercrit=None, axis=1)
		Calculate the number of input and detection within each 1-D bin.
		axis asks for which array to use for binning (1 or 2).
		"""
		self.detcrit = (self.c.detect==True)
		if othercrit == None: othercrit = ones(len(self.c.d),'bool')
		if axis == 1:
			array_input = self.array1[othercrit]
		elif axis == 2:
			array_input = self.array2[othercrit]
		dx = bins[1] - bins[0]
		bins = np.concatenate((bins, [bins[-1]+dx]))
		ninput = np.histogram(array_input, bins)[0]
		if axis == 1:
			array_det = self.array1[othercrit & self.detcrit]
		elif axis == 2:
			array_det = self.array2[othercrit & self.detcrit]
		ndetect = np.histogram(array_det, bins)[0]
		if axis == 1:
			self.bins1 = bins[:-1]
		elif axis == 2:
			self.bins2 = bins[:-1]
		self.ninput1d = ninput
		self.ndetect1d = ndetect
		self.detcomp1d = ndetect.astype('float') / np.maximum(1.0, ninput.astype('float'))
		self.mask1d = (ninput <= 0)
		
	def calc_detection2d(self, bins1, bins2, othercrit=None):
		"""
		calc_detection(self, col1, col2, bins1, bins2)
		Calculate the number of input and detection within each 2-D bin.
		The bins are specified by bins1 and bins2. bins1, bins2 could be arrays denoting 
		the bin boundaries. bins1, bins2 should contain the end points of the last bin.
		self.c must have a column named "detect".
		"""
		self.detcrit = (self.c.detect == True)
		if othercrit == None: othercrit = ones(len(self.c.d), 'bool')
		#dx1 = bins1[1] - bins1[0]
		#dx2 = bins2[1] - bins2[0]
		#bins1 = np.concatenate((bins1, [bins1[-1]+dx1]))
		#bins2 = np.concatenate((bins2, [bins2[-1]+dx2]))
		ninput = np.histogram2d(self.array1[othercrit], self.array2[othercrit], 
			bins=[bins1, bins2])[0]
		array1_det = np.compress((self.c.detect==True) & othercrit, self.array1)
		array2_det = np.compress((self.c.detect==True) & othercrit, self.array2)
		ndetect = np.histogram2d(array1_det, array2_det, bins=[bins1, bins2])[0]
		self.bins1 = bins1
		self.bins2 = bins2
		self.bw1 = bins1[1] - bins1[0]
		self.bw2 = bins2[1] - bins2[0]
		self.ninput2d = ninput
		self.ndetect2d = ndetect
		self.detcomp2d = ndetect.astype('float') / np.maximum(1.0, ninput.astype('float'))
	
	def ninput_bin1d(self, limit1, othercrit=None):
		"""
		ninput_bin1d(self, limit1, othercrit=None)
		calculate the input numbers in a 1-D bin, pre-filtered by othercrit
		if othercrit == None, no pre-filtering of sources is done.
		"""
		if othercrit == None: othercrit = ones(len(self.c.d),'bool')
		bincrit = (self.array1>=limit1[0]) & (self.array1<limit1[1]) & othercrit
		n = sum(bincrit)
		return n
		
	def ndetect_bin1d(self, limit1, othercrit=None):
		"""
		ndetect_bin1d(self, limit1, othercrit=None)
		calculate the input numbers in a 1-D bin, pre-filtered by othercrit
		if othercrit == None, no pre-filtering of sources is done.
		"""
		if othercrit == None: othercrit = ones(len(self.c.d),'bool')
		bincrit = (self.array1>=limit1[0]) & (self.array1<limit1[1]) & othercrit
		n = sum(bincrit & (self.c.detect==True))
		return n
		
	def ninput_bin2d(self, limit1, limit2):
		"""
		ninput_bin(self, limit1, limit2)
		calculate the number of input points in a specific bin with limits limit1 and
		limit2.
		""" 
		bincrit = (self.array1>=limit1[0]) & (self.array1<limit1[1])
		bincrit = bincrit & (self.array2>=limit2[0]) & (self.array2<limit2[1])
		n = sum(bincrit)
		return n
		
	def ndetect_bin2d(self, limit1, limit2):
		"""
		ndetect_bin(self, limit1, limit2)
		calculate the number of detected points in a specific bin with limits limit1 and
		limit2.
		"""
		bincrit = (self.array1>=limit1[0]) & (self.array1<limit1[1])
		bincrit = bincrit & (self.array2>=limit2[0]) & (self.array2<limit2[1])
		n = sum(bincrit & (self.c.detect==True))
		return n
		
	def show_detcomp2d(self, dxtick=5, dytick=5):
		"""
		show_detcomp(self)
		Show the 2-D detection completeness map.
		"""
		plt.imshow(self.detcomp2d.swapaxes(0,1), origin='lower', aspect='auto',
			interpolation='none')
		plt.colorbar()
		# calculate aspect ratio of the bins
		plt.xticks(arange(len(self.bins1))[::dxtick]-0.5, around(self.bins1[::dxtick],1))
		plt.yticks(arange(len(self.bins2))[::dytick]-0.5, around(self.bins2[::dytick],1))
		plt.xlabel(self.col1)
		plt.ylabel(self.col2)
		
	def calc_detcomp2d_contour(self, xsub=20, ysub=20, s=0.5):
		"""
		completeness_curve(self, complevel=0.8, dcomp=0.05, xsub=5, ysub=5, s=0)
		Interpolate the detection completenesses to get a smooth completeness curve at
		the level specified by complevel (0 < complevel < 1). xsub and ysub are the number
		of subdivisions between bins in the x (bins1) and y (bins2) directions for
		interpolation purposes.
		dcomp -- the tolerance in completeness within which one finds a match with complevel
		s -- the smoothing parameter used in cubic spline interpolation
		"""
		#if othercrit == None:
		#	othercrit = ones(len(self.c.d), 'bool')
		#if ax == None:
		#	fig = plt.figure()
		#	ax = fig.add_subplot(111)
		dx = self.bins1[1] - self.bins1[0]
		dy = self.bins2[1] - self.bins2[0]
		x, y = np.mgrid[self.bins1[0]:self.bins1[-2]:complex(0,len(self.bins1)-1),
			self.bins2[0]:self.bins2[-2]:complex(0,len(self.bins2)-1)]
		tck = interpolate.bisplrep(x+dx/2., y+dy/2., self.detcomp2d, s=s)
		# now define the new grid of points to interpolate on
		xnew, ynew = np.mgrid[self.bins1[0]:self.bins1[-2]:complex(0,len(self.bins1)*xsub),
			self.bins2[0]:self.bins2[-2]:complex(0,len(self.bins2)*ysub)]
		znew = interpolate.bisplev(xnew[:,0]+dx/2., ynew[0,:]+dy/2., tck)
		#ax.contour(xnew[:,0]+dx/2., ynew[0,:]+dy/2., znew.swapaxes(0,1), complevels,
		#	**ctrkwargs)
		self.mag_new = xnew[:,0]+dx/2.
		self.logr_new = ynew[0,:]+dy/2.
		self.detcomp2d_new = znew
		#return ax
		