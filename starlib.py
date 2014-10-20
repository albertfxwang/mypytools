#!/usr/bin/env python

import numpy as N
import pysynphot as S
import glob
from sedtools.arrayspec import ABmag
from pygoods import sextractor
import matplotlib as mpl
import matplotlib.pyplot as plt
from sedtools.bands import filters as default_filters

"""
A script to read Pickles 1998 stellar spectra library.
"""
default_filters = {'f435w':S.ObsBandpass('acs,wfc2,f435w'),
	'f606w':S.ObsBandpass('acs,wfc2,f606w'),
	'f625w':S.ObsBandpass('acs,wfc2,f625w'),
	'f775w':S.ObsBandpass('acs,wfc2,f775w'),
	'f814w':S.ObsBandpass('acs,wfc2,f814w'),
	'f850lp':S.ObsBandpass('acs,wfc2,f850lp'),
	'f098m':S.ObsBandpass('wfc3,ir,f098m'),
	'f105w':S.ObsBandpass('wfc3,ir,f105w'),
	'f110w':S.ObsBandpass('wfc3,ir,f110w'),
	'f125w':S.ObsBandpass('wfc3,ir,f125w'),
	'f140w':S.ObsBandpass('wfc3,ir,f140w'),
	'f160w':S.ObsBandpass('wfc3,ir,f160w'),
	'uvimos':S.FileBandpass('/Users/khuang/Dropbox/codes/mypytools/sedtools/bandpass/U_vimos.bp'),
	'uctio':S.FileBandpass('/Users/khuang/Dropbox/codes/mypytools/sedtools/bandpass/CANDELS_FILTERS/CTIO/U_ctio.bp')
	}


class starlib(object):
	"""
	A class that reads in a catalog of (normalized) magnitudes calculated from stellar
	spectra library. Then calculates the colors in the specified bands (must be present
	in the catalog) and make plots.
	"""
	def __init__(self,catalog):
		"""
		Set the attributes of the class.
		"""
		c = sextractor(catalog)
		bands = []
		for col in c._colnames:
			setattr(self,col,c.__getattribute__(col))
			b = col[:-4]  # get the name of the band
			bands += [b]
		self.bands = bands
	
	def getmag(self,b):
		"""
		Get the magnitude corresponding to band b.
		"""
		return getattr(self,"%s_mag"%b)
				
	def colors(self,b1,b2):
		"""
		Calculate the colors between bands b1 and b2.
		"""
		if b1 not in self.bands:
			raise KeyError, "%s not in filter list." % b1
		if b2 not in self.bands:
			raise KeyError, "%s not in filter list." % b2
		mag1 = self.getmag(b1)
		mag2 = self.getmag(b2)
		color12 = mag1 - mag2
		setattr(self,"%s_%s"%(b1,b2),color12)
		return color12
		
	def plot_colorcolor(self,b1,b2,b3,b4,ax=None,pltkw={'marker':'*','linestyle':'none'},
		label="Pickles"):
		"""
		Plot stellar loci on the color-color diagram, where the Y-axis is b1-b2, and the 
		X-axis is b3-b4.
		"""
		try:
			color12 = getattr(self,"%s_%s"%(b1,b2))
		except:
			print "%s - %s color not calculated... do it now." % (b1,b2)
			self.colors(b1,b2)
			color12 = getattr(self,"%s_%s"%(b1,b2))
		try:
			color34 = getattr(self,"%s_%s"%(b3,b4))
		except:
			print "%s - %s color not calculated... do it now." % (b3,b4)
			self.colors(b3,b4)
			color34 = getattr(self,"%s_%s"%(b3,b4))
		if ax==None:
			fig = plt.figure()
			ax = fig.add_subplot(111)
		ax.plot(color34, color12, label=label, **pltkw)
		ax.set_xlabel("%s - %s" % (b3,b4))
		ax.set_ylabel("%s - %s" % (b1,b2))
		

def calcmag(outcat,filters=default_filters,normmag=10.0,normband='f606w'):
	"""
	Reads the spectra for all *.dat and compute the magnitudes in the filters supplied
	by the dictionary filters. All magnitudes are in AB mag, and AB mag is normalized to
	normmag (default is 10) mag in V band (f606w).
	"""
	datafiles = glob.glob('/Users/khuang/Dropbox/research/stellar_library/dat/uk*.dat')
	# only use spectra with coverage extending to mid-IR (UVK library)
	mags = {}
	
	bands = filters.keys()
	for i in range(len(bands)):  # make sure the first band is normband
		if bands[i] == normband:
			bands[i] = bands[0]
			bands[0] = normband
			break
	spectype_list = []
	for df in datafiles:
		spectype = df.split('/')[-1]
		spectype = spectype.split('.')[0]  # read the spectral type name
		mags[spectype] = N.zeros(len(bands),'float')
		spectype_list += [spectype]
		print spectype
		# read the data file
		c = sextractor(df)
		wave = c._1; flux = c._2
		sp = S.ArraySpectrum(wave,flux,'angstrom','flam')  # initialize a pysynphot spectrum
		for i in range(len(bands)):  # now iterate through all bands
			b = bands[i]
			m = ABmag(sp, filters[b])
			if b == normband:
				dmag = normmag - m
				mags[spectype][i] = normmag  # set f606w magnitude to normmag
			else:
				mags[spectype][i] = m + dmag  # normalize the magnitudes in other bands
		
	# now write to an output catalog... use a SExtractor table
	header = ""
	header += "# 1 SPECTYPE\n"
	ni = 1
	for i in range(len(bands)):
		header += "# %d %s_mag\n" % ((ni+i+1),bands[i])
	
	# now write the output
	f = open(outcat,'w')
	f.write(header)
	for i in range(len(spectype_list)):
		f.write('%10s ' % spectype_list[i])
		sptype = spectype_list[i]
		for j in range(len(bands)):
			f.write('%10.2f ' % mags[sptype][j])
		f.write('\n')
	f.close()