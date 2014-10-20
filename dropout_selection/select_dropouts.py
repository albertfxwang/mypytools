#!/usr/bin/env python
# A parent class for all dropout selections

from numpy import *
from pygoods import *
import pyfits
from dropout_selection import lbg_colorcrit as lcc
import matplotlib as mpl
import matplotlib.pyplot as plt
from sedtools import plot_sed as ps
import os, glob

default_catalog = "/Users/khuang/Dropbox/Research/bivariate/udrops/gs_all_tf_h_130213a_multi_PHOTOZ_v1.0.fits"
# This is the CANDELS GOODS-S MW catalog
# The filter names should adhere to the naming convention of GOODS-S catalog columns: [instrument]_[filter]
default_bands = N.array(['ctio_u','vimos_u','acs_f435w','acs_f606w','acs_f775w',
					'acs_f814w','acs_f850lp','wfc3_f098m','wfc3_f105w','wfc3_f125w',
					'wfc3_f160w','isaac_ks','hawki_ks','irac_ch1','irac_ch2',
					'irac_ch3','irac_ch4'])
default_waves = N.array([0.3734,0.3722,0.4317,0.5918,0.7693,
					0.8047,0.9055,0.9851,1.0550,1.2486,
					1.5370,2.1605,2.1463,3.5508,4.4960,
					5.7245,7.8840])
waves_dict = {}
for i in range(len(default_bands)):
	waves_dict[default_bands[i]] = default_waves[i]


def uJy_to_mag(flux):
	mag = where(flux > 0., 23.9 - 2.5 * log10(flux), 99.0)
	return mag
	
def calc_snratio(flux, fluxerr):
	sn = where(fluxerr>0., flux/fluxerr, -1.0)
	return sn

class dropouts(lcc.colorcrit):
	"""
	A parent class for dropout selection.
	To select dropouts, the catalog needs to have the following columns:
	-- [band]_mag: magnitudes in each band to calculate colors (usually MAG_ISO)
	-- 
	"""
	def __init__(self, drop, key, c, bands=default_bands, cc_lib=lcc.cc_library, 
		outcatalog="",sn_lolim={'wfc3_f160w':5.0},sn_hilim={}):
		"""
		drop  --- the dropout band 
		key   --- the keyword of the dropout selection criteria
		c     --- the Ftable instance of the input catalog
		bands --- the filters included in the catalog; should be one of the filter names in default_bands,
		          but will change later to add custom filters
		"""
		#c = Ftable(catalog)  # use a FITS table
		self.c = c
		self.bands = bands
		self.wave = N.zeros(len(bands))
		for i in range(len(self.bands)):
			self.wave[i] = waves_dict[self.bands[i]] 
		# effective wavelengths in microns
		# Calculate S/N ratios... make sure to use positive flux errors as 1-sigma upper limit
		# otherwise set S/N=-1.0
		for b in self.bands:
			# Using ISO aperture to calculate S/N:
			setattr(self.c,'SN_%s'%b,calc_snratio(getattr(c,'%s_flux'%b),getattr(c,'%s_fluxerr'%b)))
		#self.SN_uctio = calc_snratio(c.ctio_u_flux, c.ctio_u_fluxerr)
		#self.SN_uvimos = calc_snratio(c.vimos_u_flux, c.vimos_u_fluxerr)
		#self.SN_f435w = calc_snratio(c.acs_f435w_flux, c.acs_f435w_fluxerr)
		#self.SN_f606w = calc_snratio(c.acs_f606w_flux, c.acs_f606w_fluxerr)
		#self.SN_f775w = calc_snratio(c.acs_f775w_flux, c.acs_f775w_fluxerr)
		#self.SN_f814w = calc_snratio(c.acs_f814w_flux, c.acs_f814w_fluxerr)
		#self.SN_f850lp = calc_snratio(c.acs_f850lp_flux, c.acs_f850lp_fluxerr)
		#self.SN_f098m = calc_snratio(c.wfc3_f098m_flux, c.wfc3_f098m_fluxerr)
		#self.SN_f105w = calc_snratio(c.wfc3_f105w_flux, c.wfc3_f105w_fluxerr)
		#self.SN_f125w = calc_snratio(c.wfc3_f125w_flux, c.wfc3_f125w_fluxerr)
		#self.SN_f160w = calc_snratio(c.wfc3_f160w_flux, c.wfc3_f160w_fluxerr)
		
		# if S/N < 1, replace with 1sigma upper limit
		# NOTE that some sources have flux errors either equal to or smaller than 0!
		# These sources will have mag=99.0. Should exclude those sources from selection.
		for b in self.bands:
			# Calculate limiting magnitudes from S/N
			#ston = getattr(self.c,'SN_%s'%b)
			flux = getattr(self.c,'%s_flux'%b)
			#fluxerr = getattr(self.c,'%s_fluxerr'%b)
			#setattr(self.c,'%s_mag'%b,uJy_to_mag(where(ston > 1.,flux,fluxerr)))
			# Use limiting magnitude information
			mag_fromflux = where(flux>0, uJy_to_mag(flux), 99.0)
			limmag = getattr(self.c, 'limiting_magnitude_%s' % b)
			setattr(self.c, '%s_mag'%b, minimum(mag_fromflux, limmag))
			
		# use S/N ratio limit for additional criteria
		othercrit = ones(len(c.d),'bool')
		for b in sn_lolim.keys():
			print "S/N in %s >= %.1f" % (b,sn_lolim[b])
			sncrit = (getattr(self.c,'SN_%s'%b)>=sn_lolim[b])
			othercrit = othercrit & sncrit
		for b in sn_hilim.keys():
			print "S/N in %s < %.1f" % (b,sn_hilim[b])
			sncrit = (getattr(self.c,'SN_%s'%b)<sn_hilim[b])
			othercrit = othercrit & sncrit
		# enforce that flux error be > 0
		if hasattr(c,'agn_flag'):
			othercrit = othercrit & (c.agn_flag==0)
		self.crit = othercrit  # for now... not yet performed color selection

		# Now do color selection
		lcc.colorcrit.__init__(self, cc_lib)
		lcc.colorcrit.__call__(self, drop, key)
		#udrops = cc('u','h13', lib=cc_lib)  # my custom criteria
		mag1 = getattr(self.c,"%s_mag"%self.band1)
		mag2 = getattr(self.c,"%s_mag"%self.band2)
		mag3 = getattr(self.c,"%s_mag"%self.band3)
		if hasattr(self,'band4'):  # if has a 4th band:
			mag4 = getattr(self.c,"%s_mag"%self.band4)
			self.color1 = mag1-mag2
			self.color2 = mag3-mag4
		else:
			#fecrit = ((mag1<99.0)&(mag2<99.0)&(mag3<99.0))
			self.color1 = mag1-mag2
			self.color2 = mag2-mag3
		self.select(self.color1,self.color2)  # color selection
		print sum(self.crit), sum(othercrit)
		#self.crit = ( self.crit & othercrit & fecrit )  # fold in other criteria
		self.crit = (self.crit & othercrit)
		if hasattr(self,'band4'):
			self.cname1 = "%s-%s" % (self.band1,self.band2)
			self.cname2 = "%s-%s" % (self.band3,self.band4)
		else:
			self.cname1 = "%s-%s" % (self.band1,self.band2)
			self.cname2 = "%s-%s" % (self.band2,self.band3)
		self.color1_selected = self.color1[self.crit==True]
		self.color2_selected = self.color2[self.crit==True]
		if hasattr(c,'id_1'):
			self.id_selected = c.id_1[self.crit==True]
			self.ra_selected = c.ra_1[self.crit==True]
			self.dec_selected = c.dec_1[self.crit==True]
		#print "%d sources selected." % sum(self.crit)
		
		# make a new catalog
		if len(outcatalog) > 0:
			cols = []
			for i in range(len(c.Columns)):
				cname = c.Columns[i]
				carray = c.__getitem__(cname)[self.crit]
				cformat = c.d.formats[i]
				cols += [pyfits.Column(name=cname, format=cformat, array=carray)]
	
			coldefs = pyfits.ColDefs(cols)
			tbhdu = pyfits.new_table(coldefs)
			if os.path.exists(outcatalog):
				os.system('rm %s' % outcatalog)
			tbhdu.writeto(outcatalog)
	
	def reset_crit(self):
		"""
		Only reset color1, color2, and redshift attributes of self.
		"""
		self.color1_selected = self.color1[self.crit==True]
		self.color2_selected = self.color2[self.crit==True]
		if hasattr(self.c,'id_1'):
			self.id_selected = self.c.id_1[self.crit==True]
			self.ra_selected = self.c.ra_1[self.crit==True]
			self.dec_selected = self.c.dec_1[self.crit==True]
			self.get_selected_z()
		
	def get_selected_attr(self, attr):
		"""
		Set an attribute of self that belong to the subset of sources selected as dropouts,
		using the attribute attr.
		"""
		attr_array = getattr(self, attr)
		attr_selected = attr_array[self.crit==True]
		setattr(self,"%s_selected"%attr,attr_selected)
		return attr_selected
	
	def get_selected_z(self):
		"""
		Return arrays of the length = number of dropouts selected.
		"""
		sel = (self.crit==True)
		self.zqcrit_selected = compress(self.crit, self.c.spec_z_dq)
		self.zqcrit_selected = in1d(self.zqcrit_selected, [1,2])
		self.specz_selected = compress(self.crit,self.c.spec_z)
		self.specz_selected = where(self.zqcrit_selected, self.specz_selected, -9.0) # only spec-z's with high qualities
		self.specz_dq_selected = compress(self.crit,self.c.spec_z_dq)		
		self.photz_selected = compress(self.crit,self.c.photo_z)
		self.photz_median_selected = compress(self.crit,self.c.photo_z_median)
		# flag out photo-z of objects already have good spec-z
		#self.photz_selected = where(self.specz_selected>0., -9., self.photz_selected)
		#self.photz_median_selected = where(self.specz_selected>0.,-9.,self.photz_median_selected)
		self.photz_upper68_selected = compress(self.crit,getattr(self.c,'photo_z_upper_68%'))
		self.photz_lower68_selected = compress(self.crit,getattr(self.c,'photo_z_lower_68%'))
	
	def plot_zhist(self, zbin=arange(0.,8.,0.2), hist_kw={}):
		"""
		Makes a two-panel plot, one with spec-z, the other with photo-z.
		"""
		fig = plt.figure(figsize=(12,6))
		ax1 = fig.add_subplot(121)
		ax2 = fig.add_subplot(122)
		# First plot spec-z
		#zqcrit = in1d(self.c.spec_z_dq, [1,2])
		ax1.hist(self.specz_selected, zbin, **hist_kw)
		#specz_sel = compress(self.crit & zqcrit,self.c.spec_z)
		#frac_zst3 = sum(specz_sel < 2.8) / float(len(specz_sel))
		#print "Fraction of selected sources with specz<2.8:", frac_zst3
		ax1.set_xlabel('Spec-z')
		# Then plot photo-z
		#photzcrit = ((self.c.spec_z<0.) | (self.c.spec_z_dq==3)) & (self.c.photo_z>0.)
		ax2.hist(self.photz_selected, zbin, label='Best', **hist_kw)
		ax2.hist(self.photz_median_selected, zbin, label='Median', **hist_kw)
		ax2.legend(loc=1)
		ax2.set_xlabel('Photo-z')
		#photz_sel = compress(photzcrit,self.c.photo_z)
		#frac_zst3 = sum(photz_sel < 2.8) / float(len(photz_sel))
		#self.specz_sel = compress(self.crit, self.c.spec_z)
		#self.photoz_sel = compress(self.crit, self.c.photo_z)
		#self.specz_quality = compress(self.crit, self.c.spec_z_dq)
		#print "Fraction of selected sources with photz<2.8:", frac_zst3
	

	def get_fluxes(self, objectid):
		"""
		Returns the fluxes and flux errors of a given object.
		"""
		i = N.arange(len(self.c.d))[self.c.id_1==objectid][0]
		fluxes = zeros(len(self.bands))
		fluxerrs = zeros(len(self.bands))
		for j in range(len(self.bands)):
			fluxes[j] = getattr(self,'%s_flux'%b)[i]
			fluxerrs[j] = getattr(self,'%s_fluxerr'%b)[i]
		#fluxes = [self.c.ctio_u_flux[i],self.c.vimos_u_flux[i],self.c.acs_f435w_flux[i],
		#			self.c.acs_f606w_flux[i],self.c.acs_f775w_flux[i],self.c.acs_f814w_flux[i],
		#			self.c.acs_f850lp_flux[i],self.c.wfc3_f098m_flux[i],self.c.wfc3_f105w_flux[i],
		#			self.c.wfc3_f125w_flux[i],self.c.wfc3_f160w_flux[i],self.c.isaac_ks_flux[i],
		#			self.c.hawki_ks_flux[i],self.c.irac_ch1_flux[i],self.c.irac_ch2_flux[i],
		#			self.c.irac_ch3_flux[i],self.c.irac_ch4_flux[i]]
		#fluxerrs = [self.c.ctio_u_fluxerr[i],self.c.vimos_u_fluxerr[i],self.c.acs_f435w_fluxerr[i],
		#			self.c.acs_f606w_fluxerr[i],self.c.acs_f775w_fluxerr[i],self.c.acs_f814w_fluxerr[i],
		#			self.c.acs_f850lp_fluxerr[i],self.c.wfc3_f098m_fluxerr[i],self.c.wfc3_f105w_fluxerr[i],
		#			self.c.wfc3_f125w_fluxerr[i],self.c.wfc3_f160w_fluxerr[i],self.c.isaac_ks_fluxerr[i],
		#			self.c.hawki_ks_fluxerr[i],self.c.irac_ch1_fluxerr[i],self.c.irac_ch2_fluxerr[i],
		#			self.c.irac_ch3_fluxerr[i],self.c.irac_ch4_fluxerr[i]]
		#fluxes = N.array(fluxes)
		#fluxerrs = N.array(fluxerrs)
		return fluxes, fluxerrs

	def plot_photometry(self, objectid, ax=None, label=None, errorbarkw={}, title=""):
		"""
		Plot the photometric points of a given object (specified by ID number)
		"""
		i = N.arange(len(self.c.d))[self.c.id_1==objectid][0]
		fluxes,fluxerrs = self.get_fluxes(objectid)
		if ax == None:
			fig = plt.figure()
			ax = fig.add_subplot(111)
		ps.plot_photometry(self.wave,fluxes,fluxerrs,waveunit='micron',fluxunit='mJy',
			ax=ax,label=label,sigma=1.0,errorbarkw=errorbarkw,title=title)
		xticks = [0.3,0.5,1.0,3.0,5.0,8.0]
		ax.set_xticks(xticks)
		ax.set_xticklabels(xticks)
		ax.set_xlim(0.2,10.0)
		#plt.ylim(1.e-3,1.e3)

	def write_catalog(self, catalog_name, bands_include=None):
		if bands_include == None:
			bands_include = self.bands
		# Write to a catalog
		if os.path.exists(catalog_name):
			os.system('rm %s' % catalog_name)
		cols = []
		#catalog_bands = ['uvimos','uctio','f435w','f606w','f775w','f850lp','f105w','f125w','f160w']
		# first get the magnitudes in important bands
		for b in bands_include:
			cols += [pyfits.Column(name='%s_mag'%b,array=getattr(self,'%s_mag'%b)[self.crit],format='D')]
		# then get all the columns from the master catalog
		for i in range(len(self.c.Columns)):
			col = self.c.Columns[i]
			format = self.c.d.formats[i]
			newarray = getattr(self.c,col)[self.crit]
			cols += [pyfits.Column(name=col,array=newarray,format=format)]
		tbhdu = pyfits.new_table(pyfits.ColDefs(cols))
		tbhdu.writeto(catalog_name)

	
def select_bdrops_goodss(c, filename, key='g04'):
	""" Select B-dropouts
	usage: bdrops = select_bdrops_goodss(c, filename)
	where filename is the name of the output catalog.
	"""
	
	SN_b = c.acs_f435w_flux / c.acs_f435w_fluxerr
	SN_v = c.acs_f606w_flux / c.acs_f606w_fluxerr
	SN_z = c.acs_f850lp_flux / c.acs_f850lp_fluxerr
	
	bmag = where((SN_b > 1.)&(c.acs_f435w_flux>0), uJy_to_mag(c.acs_f435w_flux), uJy_to_mag(c.acs_f435w_fluxerr))
	vmag = where((SN_v > 1.)&(c.acs_f606w_flux>0), uJy_to_mag(c.acs_f606w_flux), uJy_to_mag(c.acs_f606w_fluxerr))
	zmag = where((SN_z > 1.)&(c.acs_f850lp_flux>0), uJy_to_mag(c.acs_f850lp_flux), uJy_to_mag(c.acs_f850lp_fluxerr))
	
	bmv = bmag - vmag
	vmz = vmag - zmag
	othercrit = ((SN_z >= 5.0) & (SN_v >= 3.0) & ((c.class_star <= 0.8)|(zmag >= 26.0)))
	
	cc = lcc.colorcrit()
	bdrops = cc('b', key)
	bdrops.select(bmv, vmz)
	bdrops.crit = (bdrops.crit & othercrit)
	bdrops.color1_selected = bdrops.color1[bdrops.crit]
	bdrops.color2_selected = bdrops.color2[bdrops.crit]
	print sum(bdrops.crit)
	
	cols = []
	for i in range(len(c.Columns)):
		cname = c.Columns[i]
		carray = c.__getitem__(cname)[bdrops.crit]
		cformat = c.d.formats[i]
		cols += [pyfits.Column(name=cname, format=cformat, array=carray)]
	
	coldefs = pyfits.ColDefs(cols)
	tbhdu = pyfits.new_table(coldefs)
	if os.path.exists(filename):
		os.system('rm %s' % filename)
	tbhdu.writeto(filename)
	return bdrops
	
def select_idrops_goodss(c, filename, key='g11', veto=True):
	"""Select i-dropouts
	usage: idrops = select_idrops_goodss(c, filename, key)
	where filename is the name of the output catalog.
	"""
	SN_b = c.acs_f435w_flux / c.acs_f435w_fluxerr
	SN_v = c.acs_f606w_flux / c.acs_f606w_fluxerr
	SN_i = c.acs_f775w_flux / c.acs_f775w_fluxerr
	SN_z = c.acs_f850lp_flux / c.acs_f850lp_fluxerr
	SN_J = c.wfc3_f125w_flux / c.wfc3_f125w_fluxerr
	SN_H = c.wfc3_f160w_flux / c.wfc3_f160w_fluxerr
	
	imag = where((SN_i > 1.)&(c.acs_f775w_flux>0), uJy_to_mag(c.acs_f775w_flux), uJy_to_mag(c.acs_f775w_fluxerr))
	zmag = where((SN_z > 1.)&(c.acs_f850lp_flux>0), uJy_to_mag(c.acs_f850lp_flux), uJy_to_mag(c.acs_f850lp_fluxerr))
	jmag = where((SN_J > 1.)&(c.wfc3_f125w_flux>0), uJy_to_mag(c.wfc3_f125w_flux), uJy_to_mag(c.wfc3_f125w_fluxerr))
	hmag = where((SN_H > 1.)&(c.wfc3_f160w_flux>0), uJy_to_mag(c.wfc3_f160w_flux), uJy_to_mag(c.wfc3_f160w_fluxerr))
	
	imz = imag - zmag
	zmj = zmag - jmag
	if veto:  # use B-band and V-band as veto bands
		othercrit = ((SN_H >= 5.) & (SN_z >= 3.) & (SN_J >= 3.) & (SN_b<=2.0) & (SN_v<=2.0)) #& (c.acs_f775w_fluxerr > 0.)
	else:
		othercrit = ((SN_H >= 5.) & (SN_z >= 3.) & (SN_J >= 3.)) #& (c.acs_f775w_fluxerr > 0.))
	
	cc = lcc.colorcrit()
	idrops = cc('i', key)
	idrops.select(imz, zmj)
	idrops.crit = (idrops.crit & othercrit)
	idrops.color1_selected = idrops.color1[idrops.crit]
	idrops.color2_selected = idrops.color2[idrops.crit]
	print sum(idrops.crit)
	
	cols = []
	for i in range(len(c.Columns)):
		cname = c.Columns[i]
		carray = c.__getitem__(cname)[idrops.crit]
		cformat = c.d.formats[i]
		cols += [pyfits.Column(name=cname, format=cformat, array=carray)]
	
	coldefs = pyfits.ColDefs(cols)
	tbhdu = pyfits.new_table(coldefs)
	if os.path.exists(filename):
		os.system('rm %s' % filename)
	tbhdu.writeto(filename)
	return idrops