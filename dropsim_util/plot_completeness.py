#!/usr/bin/env python

from numpy import *
from pygoods import sextractor, Ftable
import matplotlib as mpl
import matplotlib.pyplot as plt

# A quick script to make completeness plots for SExtractor simulations

class plt_comp(Ftable):
	def __init__(self,catalog,bands=['H']):
		Ftable.__init__(self,catalog)  # read the FITS table
		self.bands = bands
		self.catalog = catalog

	def __len__(self):
		return len(self.d)

	def calc_completeness_map(self,mbins,logrbins,band,show=True):
		"""
		calculate completenesses in magnitude bins and logR bins. Note that the bins should
		contain the upper bound of the last bin.
		"""
		mag_in = getattr(self,'%s_mag_in'%band.lower())
		re_in = getattr(self,'re_in')
		logre_in = log10(re_in)
		n_input = histogram2d(mag_in,logre_in,bins=[mbins,logrbins])[0]
		n_detect = histogram2d(mag_in[self.detect==True],logre_in[self.detect==True],bins=[mbins,logrbins])[0]
		self.compmap = n_detect.astype('float') / maximum(n_input,1.0)
		self.mbins = mbins
		self.logrbins = logrbins
		if show:
			plt.figure()
			plt.imshow(self.compmap.swapaxes(0,1),origin='lower')
			plt.xticks(arange(len(mbins))[::2],mbins[::2])
			plt.yticks(arange(len(logrbins))[::2],logrbins[::2])
			plt.xlabel('input magnitude',size=14)
			plt.ylabel('input logRe', size=14)
			plt.title('Completeness map in %s-band' % band)
			plt.colorbar()

	def bincrit(self,band,mag_in_lim=[23.0,23.5],logre_in_lim=[0.2,0.4]):
		"""
		Construct the boolean array of which objects lay within a certain limits of input magnitude and Re.
		"""
		b = band.lower()
		mag_in = getattr(self,'%s_mag_in'%b)
		logre_in = log10(self.re_in)
		bincrit = (mag_in>=mag_in_lim[0]) & (mag_in<mag_in_lim[1])
		bincrit = bincrit & (logre_in>=logre_in_lim[0]) & (logre_in<logre_in_lim[1])
		return bincrit

	def plot_comp_fixed_re(self,band,ax=None,label="0.2<log(Re)<0.4",logre_in_lim=[0.2,0.4],
		mag_in_array=arange(21.,29.5,0.5),othercrit=None,pltkw={}):
		"""
		Plot the completeness curve at a given input log(Re) bin, as a function of input magnitude.
		mag_in_array should contain both edges of the limit.
		"""
		if ax==None:
			fig = plt.figure()
			ax = fig.add_subplot(111)
		comp_array = zeros(len(mag_in_array)-1)
		if othercrit==None:
			othercrit=ones(len(self.d),'bool')
		# calculate the completeness in each magnitude bin
		for i in range(len(comp_array)):
			mag_lim = [mag_in_array[i],mag_in_array[i+1]]
			bincrit_i = self.bincrit(band,mag_in_lim=mag_lim,logre_in_lim=logre_in_lim)&othercrit
			bincrit_d = bincrit_i & (self.detect==True)
			comp_array[i] = float(sum(bincrit_d)) / maximum(float(sum(bincrit_i)),1.0)
		# Plot
		dmag = mag_in_array[1]-mag_in_array[0]
		line = ax.plot(mag_in_array[:-1],comp_array,label=label,**pltkw)[0]
		if ax.get_xlabel()=="":
			plt.xlabel('input magnitude in %s' % band, size=14)
		if ax.get_ylabel()=="":
			plt.ylabel('completeness', size=14)
		if ax.get_title()=="":
			plt.title('')
		return line



