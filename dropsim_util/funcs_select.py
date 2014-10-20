#!/usr/bin/env python

import numpy as np
import os
from pygoods import *
from lbg_colorcrit import cc_library, colorcrit
from funcs_detect import detfunc
import matplotlib as mpl
import matplotlib.pyplot as plt
import lognormal, gauss
import draw_from_pdf as dfp

"""
# A script for calculating LBG selection functions using simulation catalogs.
# Goals:
# 1. To calculate P(M, R, z) (or P(z) in 2D bins of some other combinations)
# 2. To calculate P(M, R) 
# Will also calculate the detection completenesses.
"""

class selfunc(detfunc):
	"""
	A class for the selection functions from a dropout simulation.
	One can calculate the selection completeness in any kind of bins (or at least this is 
	the goal)!
	"""
	def __init__(self, catalog):
		"""
		Initialize the catalog and color criteria.
		"""
		detfunc.__init__(self, catalog)
		self.colorcrit = colorcrit()
		self.bins1 = None; self.bins2 = None
		self.col1 = None; self.col2 = None
		self.zarr = arange(4., 8., 0.1)
		self.ndim = 1
		self.last_Pz = zeros(len(self.zarr))
		
			
	def use_lbg(self, drop, key):
		"""
		Use a specific set of color criteria for LBGs.
		"""
		self.selcrit = self.colorcrit(drop, key)
		self.drop = drop
		self.key = key
		
		
	def calc_selection(self, drop, key, color1, color2):
		"""
		calc_selection1d(self, drop, key)
		Calculate the color-selection function P(col1, col2, z). Does NOT calculate object detection
		completenesses! Also assuming that object detection completenesses don't depend
		on input redshift and only on self.col1 and self.col2.
		""" 
		self.use_lbg(drop, key)
		self.colorname1 = self.colorcrit.color1
		self.colorname2 = self.colorcrit.color2
		self.colorcrit.select(color1, color2)  # calls lbg_colorcrit.select
		self.selcrit = self.colorcrit.crit  # selection array
		if self.ndim == 1:
			self.array1_sel = np.compress(self.colorcrit.crit, self.array1)
		elif self.ndim == 2:
			self.array1_sel = np.compress(self.colorcrit.crit, self.array1)
			self.array2_sel = np.compress(self.colorcrit.crit, self.array2)
		
	
	def P_z_1d(self, limit1, detcomp=True):
		"""
		P_z(self, limit1, zarr=self.last_zarr)
		Calculate P(z) in (self.col1) bins
		"""
		self.last_limit1 = limit1
		bincrit = (self.array1>=limit1[0]) & (self.array1<limit1[1])
		self.last_Pz = zeros(len(self.zarr))
		dz = self.zarr[1] - self.zarr[0]
		for i in range(len(self.zarr)):
			zcrit = (self.c.zin>=self.zarr[i]) & (self.c.zin<(self.zarr[i]+dz))
			if detcomp:
				nsel = sum(bincrit & zcrit & self.selcrit & self.detcrit)
				ndet = sum(bincrit & zcrit & self.detcrit)
			else:
				nsel = sum(bincrit & zcrit & self.selcrit)
				ndet = sum(bincrit & zcrit)
			self.last_Pz[i] = nsel.astype('float') / np.maximum(1.0, ndet).astype('float')

		
	def P_z_2d(self, limit1, limit2, detcomp=True):
		"""
		P_z(self, limit1, limit2, detcomp=True)
		Calculate P(z) in (self.col1, self.col2) bins
		if detcomp==True, P(z) will be multiplied by the detection completeness
		"""
		self.last_limit1 = limit1
		self.last_limit2 = limit2
		bincrit = (self.array1>=limit1[0]) & (self.array1<limit1[1])
		bincrit = bincrit & (self.array2>=limit2[0]) & (self.array2<limit2[1])
		self.last_Pz = zeros(len(self.zarr))
		dz = self.zarr[1] - self.zarr[0]
		for i in range(len(self.zarr)):
			zcrit = (self.c.z_in>=self.zarr[i]) & (self.c.z_in<(self.zarr[i]+dz))
			if detcomp:
				nsel = sum(bincrit & zcrit & self.selcrit & self.detcrit)
				ndet = sum(bincrit & zcrit & self.detcrit)
			else:
				nsel = sum(bincrit & zcrit & self.selcrit)
				ndet = sum(bincrit & zcrit)
			self.last_Pz[i] = nsel.astype('float') / np.maximum(1.0, ndet).astype('float')
			
	def P_z_1d_weighted(self, limit1, bins2, weights2, detcomp=True):
		"""
		Calculate the weighted PDF in the bin specified by limit1.
		The PDFs used to calculate the weighted mean are the ones in the 2-dimensional bins
		specified by limit1 and successive bins in bins2. Then each PDF in bins2 will be
		multiplied by the appropriate weights given by weights2 to produce the weighted
		PDF in limit1.
		bins2 should contain the **lower edge** of the bins used, the bin width will be 
		figured out automatically.
		"""
		P_z_list = []  # to store the p(z) in 2D bins
		bw = bins2[1] - bins2[0]
		for i in range(len(bins2)):  # now loop through all the bins
			lim2 = [bins2[i], bins2[i]+bw]
			self.P_z_2d(limit1, lim2, detcomp=detcomp)
			P_z_list += [self.last_Pz.copy()]
		# calculate the weighted average
		wPz = zeros(len(self.zarr))
		for i in range(len(bins2)):
			wPz = wPz + weights2[i] * P_z_list[i]
		wPz = wPz / sum(weights2)
		return wPz
		
	def lognormal_weights(self, bins, logr0, sigma):
		"""
		Calculate the weights of a lognormal distribution for the bins specified by bins2.
		The lognormal distribution is specified by the mean r0 and width sigma.
		bins only contain the lower edges of each bin.
		"""
		bw = bins[1] - bins[0]
		# draw a large number of points from a lognormal distribution, then rebin to calculate
		# the weights
		def pdf(x):
			return gauss.gauss(x, logr0, sigma)
		drawn = dfp.draw_from_pdf_func(10000, pdf, [bins[0], bins[-1]+bw])
		bins = concatenate((bins, [bins[-1]+bw]))
		weights = np.histogram(drawn, bins)[0]
		weights = weights / float(sum(weights))
		return weights
		
	def P_z_1d_weighted_lognormal(self, limit1, bins2, logr0, sigma, detcomp=True):
		"""
		Calculate the weighted P(z) for a given range of input M1500; the weights are
		lognormal size distribution at a given luminosity. Should also put in beta as well?
		"""
		weights = self.lognormal_weights(bins2, logr0, sigma)
		Pz_weighted = self.P_z_1d_weighted(limit1, bins2, weights, detcomp=detcomp)
		return Pz_weighted
			
		
class Haojing_z6(selfunc):
	"""
	A class that calculates the selection function of z~6 dropouts using Haojing's criteria
	and CANDELS data.
	"""
	def __init__(self, catalog, bins1=np.arange(-25.0,-14.5,0.5), 
		bins2=np.arange(-1.0,2.2,0.2)):
		selfunc.__init__(self, catalog)
		self.use_lbg('i','haojing')
		self.zp_f814w = 26.49113
		self.zp_f125w = 25.94333
		self.zp_f160w = 26.25
		self.bins1 = bins1
		self.bins2 = bins2
		self.zarr = arange(4.,9.6,0.1)
		
	def select(self, colname1, array1, colname2="", array2=[]):
		self.set_cols(colname1, array1, colname2, array2)
		self.calc_detection2d(self.bins1, self.bins2)
		i_magerr_iso = where(self.c.i_magerr_iso>0, self.c.i_magerr_iso, 99.0)
		i_fluxerr_iso = where(self.c.i_fluxerr_iso>0, self.c.i_fluxerr_iso, 1.e9)
		self.imag = where(1.0857/i_magerr_iso>=2.0, self.c.i_mag_iso, 
			self.zp_f814w-2.5*log10(2.*i_fluxerr_iso))
		# get the color criteria
		self.imj = self.imag - self.c.j_mag_iso
		self.jmh = self.c.j_mag_iso - self.c.h_mag_iso
		self.calc_selection('i','haojing', self.imj, self.jmh)
		# other criteria
		othercrit = (1.0857/self.c.j_magerr_iso>=5.0) & (1.0857/self.c.h_magerr_iso>=5.0)
		self.selcrit = self.selcrit & othercrit
		# calculate completeness
		
	def plot_Pz2d(self, limit1, limit2, detcomp=True, newfig=True, **pltargs):
		self.P_z_2d(limit1, limit2, detcomp=detcomp)
		if newfig:
			plt.figure()
		plt.plot(self.zarr, self.last_Pz, linestyle='-', **pltargs)
		plt.xlabel('Input Redshift')
		plt.ylabel('Completeness')
		plt.title('%s in [%.1f, %.1f]; %s in [%.1f, %.1f]' % (self.col1,limit1[0],limit1[1],
			self.col2,limit2[0],limit2[1]))
			
	