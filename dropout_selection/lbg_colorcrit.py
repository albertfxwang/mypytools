#!/usr/bin/env python

from numpy import *
import matplotlib.pyplot as plt
import numexpr as ne
import yaml
import os, sys

lcc_dir = '/Users/khuang/Dropbox/codes/mypytools/dropout_selection'
cc_library = yaml.load(open(os.path.join(lcc_dir, 'LBG_criteria.data.yml')))

class colorcrit(object):
	""" define a base class for dropout color criteria """
	def __init__(self):
		self.lib = cc_library
		self.droplist = cc_library.keys()

	def __call__( self, drop, label, lib=None):
		"""returns a colorcrit class instance
		usage: bdrops = colorcrit( drop, label )"""
		self.drop = drop
		self.label = label
		# ukey = drop+'drops_'+key
		dropkey = '%s_drop' % drop
		dropdic = self.lib[dropkey][label]
		for k in dropdic.keys():
			setattr( self, k, self.lib[dropkey][label][k] )
		# Now defines the color cuts using coefficients
		# self.cuts = []
		self.colorname1 = "%s-%s" % (self.bands[0], self.bands[1])
		if len(self.bands) == 4:
			self.colorname2 = "%s-%s" % (self.bands[2], self.bands[3])
		else:
			self.colorname2 = "%s-%s" % (self.bands[1], self.bands[2])
		return self

	def borders(self, xarray):
		self.y_low = []
		# self.y_hi = ones(len(xarray)) * yhi
		if len(self.coeffs) == 1:
			self.y_low = self.coeffs[0][0]
		else:
			for i in range(len(self.coeffs)):
				if len(self.coeffs[i]) == 2:
					# in the form of y >= c0 + c1 * x
					print self.coeffs[i][0], self.coeffs[i][1]
					self.y_low += [self.coeffs[i][0] + self.coeffs[i][1] * xarray]
				elif len(self.coeffs[i]) == 4:
					# in the form of (y >= c0 + c1 * x) | (y >= c2 + c3 * x)
					y1 = self.coeffs[i][0] + self.coeffs[i][1]*xarray
					y2 = self.coeffs[i][2] + self.coeffs[i][3]*xarray
					self.y_low += [minimum(y1, y2)]	
			self.y_low = maximum(*self.y_low)		
			
	# def plotcrit( self, xlo, yhi, ax=None, ls='-', color='black', lw=1.0):
	# 	"""Plot the color-selection boundary as lines.
	# 	usage: bdrops.colorcrit( xlo, yhi, ax=None, ls='-', color='black' )"""
	# 	if ax == None:
	# 		fig = plt.figure()
	# 		ax = fig.add_subplot( 111 )
	# 	xvertices = [xlo] + self.xvertices + [self.xvertices[-1]]
	# 	yvertices = [self.yvertices[0]] + self.yvertices + [yhi]
	# 	for i in range( 0, len( xvertices )-1 ):
	# 		ax.plot( [xvertices[i], xvertices[i+1]], [yvertices[i], yvertices[i+1]], 
	# 			ls=ls, color=color, lw=lw )

	def plotcrit_fill(self, xlo, yhi, ax=None, fc='blue', ec='none', alpha=0.2,
	                  hatch='', dx=0.01, xmax_plot=3.0, ymin_plot=0.0):
		"""
		Plot a filled polygon filled with color & hatch of choice.
		"""
		if ax == None:
			fig = plt.figure()
			ax = fig.add_subplot(1,1,1)
		# xvertices = [xlo] + self.xvertices + [self.xvertices[-1]] + [xlo]
		# yvertices = [self.yvertices[0]] + self.yvertices + [yhi] + [yhi]
		xarray = arange(xlo, self.color2lim+dx, dx)
		print "len(xarray)", len(xarray)
		self.y_hi = ones(len(xarray)) * yhi
		self.borders(xarray)
		ax.fill_between(xarray, self.y_hi, self.y_low, alpha=alpha, facecolor=fc, 
		                edgecolor=ec, hatch=hatch)
		ax.set_xlim(xmax=xmax_plot)
		ax.set_ylim(ymin=ymin_plot)

		# ax.fill(xvertices, yvertices, fc=fc, ec=ec, hatch=hatch)
		return ax
				
	def select( self, color1, color2 ):
		"""Do the color selection given the two colors
		usage: bdrops.select( color1, color2 )
		Any S/N threshold is applied elsewhere... except for i-dropouts?
		The selection boolean array is stored as self.crit"""
		self.color1 = color1
		self.color2 = color2
		self.crit = eval(self.eval)
		#self.crit = ne.evaluate(self.eval)
		self.color1_selected = color1[self.crit]
		self.color2_selected = color2[self.crit]
		
	def plot_selected( self, xarray=None, yarray=None, ax=None, pltkwarg={}):
		"""Plot any quantity for the selected sample
		usage: bdrops.plot_selected( xarray, yarray, ax=None, pltkwarg={} )
		xarray, yarray should have the same length as self.crit"""
		if ax == None:
			fig = plt.figure()
			ax = fig.add_subplot(111)
		if xarray==None: xarray=self.color2
		if yarray==None: yarray=self.color1
		ax.plot( xarray[self.crit], yarray[self.crit], **pltkwarg)
		return ax
		
		

