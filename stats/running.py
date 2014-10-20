#!/usr/bin/env python

from numpy import *
import robust

"""
Calculate the running average, median, and standard deviation of an array x
indexed by array y.
"""

class running(robust.robust):
	def __init__(self,x,y,nsamp=50,interval=50):
		self.yval = y   # the array on which we calculate the running averages, etc.
		self.xval = x   # the array that we index the array y
		self.indices = argsort(x)
		self.xval = sort(x)   # sort x
		self.yval = self.yval.take(self.indices)   # shuffle y based on sorted x
		if nsamp==0:
			nsamp = int(len(self.yval)/50)
		if interval==0:
			interval=1
		self.nsamp = nsamp   # number of data points to calculate the running average each time
		self.interval = interval
		self.npts = len(self.xval)  # total number of data points

		# The running statistics at each y value is calculated using nsamp points around that point
		# So the running statistics don't always sample the same interval in y
		if self.nsamp % 2 == 0:  # if self.nsamp is even
			self.x_bins = self.xval[nsamp/2-1:-(nsamp/2):self.interval]
			self.i_bins = arange(nsamp/2-1,self.npts-(nsamp/2),self.interval)
		else:
			self.x_bins = self.xval[(nsamp-1)/2:-(nsamp-1)/2:self.interval]
			self.i_bins = arange((nsamp-1)/2,self.npts-(nsamp-1)/2,self.interval)
		self.nbins = len(self.x_bins)
		print len(self.x_bins),len(self.i_bins)
		self.mean_bins = zeros(self.nbins)
		self.median_bins = zeros(self.nbins)
		self.std_bins = zeros(self.nbins)
		self.loc_robust_bins = zeros(self.nbins)
		self.scale_robust_bins = zeros(self.nbins)

	def sample(self, i):
		"""
		Sample nsamp points around the index i.
		"""
		if self.nsamp % 2 == 0:
			i0 = i-(self.nsamp/2-1)
			i1 = i+(self.nsamp/2)
		else:
			i0 = i-((self.nsamp-1)/2.)
			i1 = i+((self.nsamp-1)/2.)
		return self.yval[i0:i1]

	def mean(self):
		for j in range(self.nbins):
			i = self.i_bins[j]
			self.mean_bins[j] = mean(self.sample(i))

	def median(self):
		for j in range(self.nbins):
			i = self.i_bins[j]
			self.median_bins[j] = median(self.sample(i))

	def std(self):
		for j in range(self.nbins):
			i = self.i_bins[j]
			self.std_bins[j] = std(self.sample(i))

	def calc_robust(self):
		for j in range(self.nbins):
			i = self.i_bins[j]
			robust.robust.__init__(self,self.sample(i))
			self.loc_robust_bins[j],self.scale_robust_bins[j],w = self.loc_scale()