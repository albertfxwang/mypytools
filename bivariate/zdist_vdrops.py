#!/usr/bin/env python

import numpy as N
from pygoods import sextractor, Ftable
import cPickle
import zdist
import cosmoclass
import matplotlib.pyplot as plt
from dropout_selection import lbg_colorcrit as lcc

"""
Defines a special class for the P(z) grids of V-dropouts.
"""

filter_dic = {'acs_f435w':'b','acs_f606w':'v','acs_f775w':'i',
	'acs_f850lp':'z',
	'wfc3_f125w':'j','wfc3_f160w':'h'}
zeropoints = {'wfc3_f160w':25.960,'acs_f435w':25.673,'acs_f606w':26.505,
	'acs_f775w':25.678,'acs_f850lp':24.867,
	'wfc3_f125w':26.250}

# Area in steradians
goods_area = 2.6146e-5
udf_area = 9.308e-7

class zdgrid_vdrops(zdist.zdgrid,object):
	"""
	A special class for V-dropouts.
	Selected with V-i v.s. i-z colors.
	"""
	def __init__(self, M0, M1, dM, logR0, logR1, dlogR, z0, z1, dz, area):
		zdist.zdgrid.__init__(self, M0, M1, dM, logR0, logR1, dlogR, z0, z1, dz, 
			'b', area)

	def __call__(self,simcatalog,field,plot_diag=False,ktest=None,
		bands=['acs_f435w','acs_f606w','acs_f775w','acs_f850lp'],
		sn_lolim={'acs_f850lp':5.0},
		sn_hilim={'acs_f435w':5.0},
		interpolate=False,dznew=0.1):
		"""
		simcatalog --- the simulation catalog as a FITS table.
		"""
		# Initialize i-dropout color criteria
		self.lcc = lcc.colorcrit(lib=lcc.cc_library)
		self.lcc = self.lcc('v','g04')  
		# Use Giavalisco et al. 2004 criteria for V-dropouts
		c = Ftable(simcatalog)
		# Now initialize the attributes
		for b in bands:
			bs = filter_dic[b]  # short name of filter b
			flux_iso = getattr(c,'%s_flux_iso'%bs)
			fluxerr_iso = getattr(c,'%s_fluxerr_iso'%bs)
			mag_iso = getattr(c,'%s_mag_iso'%bs)
			# calculate 1-sigma upper limits on sources with S/N < 1
			mag_iso = N.where(flux_iso/fluxerr_iso>1.0,mag_iso,
			                  zeropoints[b]-2.5*N.log10(fluxerr_iso))
			setattr(self,b+'_mag',mag_iso)  
			# MAG_ISO --> for color calculation
			setattr(self,b+'_flux',getattr(c,'%s_flux_auto'%bs))  
			# FLUX_AUTO
			setattr(self,b+'_fluxerr',getattr(c,'%s_fluxerr_auto'%bs))  
			# FLUXERR_AUTO
			setattr(self,b+'_sn',
				getattr(self,'%s_flux'%b)/getattr(self,'%s_fluxerr'%b))  
			# S/N calculated using FLUX_AUTO
		# Now construct S/N criteria
		self.sncrit = N.ones(len(c.d),'bool')  
		for b in sn_lolim.keys():  # enforce S/N lower limits
			self.sncrit = self.sncrit & (getattr(self,b+'_sn')>=sn_lolim[b])
		for b in sn_hilim.keys():  # enforce S/N upper limits (veto bands)
			self.sncrit = self.sncrit & (getattr(self,b+'_sn')<sn_hilim[b])
		print "Total number of objects satisfying the S/N criteria:", \
			sum(self.sncrit)
		print "Do selections..."
		self.color1 = self.acs_f606w_mag-self.acs_f775w_mag
		self.color2 = self.acs_f775w_mag-self.acs_f850lp_mag
		self.lcc.select(self.color1,self.color2)  # do color selection!
		self.lcc.crit = self.lcc.crit & self.sncrit  # enforce S/N criteria
		self.dropcrit = self.lcc.crit & (c.detect==True)  # just in case
		self.detect = c.detect.copy()
		print "Selection done."
		print "Total number of objects in the catalog: %d" % len(c.d)
		print "Total number selected as i-dropouts: %d" % (sum(self.dropcrit))

		# Now invoke zdist.zdgrid.__call__ to calculate P(z)
		zdist.zdgrid.__call__(self,c,interpolate=interpolate,dznew=dznew,
			plot_diag=plot_diag,ktest=ktest)

	def write(self,outname):
		zdist.write_zdgrid(self,outname)

def make_zdgrid_vdrops(simcatalog, picklename, field, z0=3.5, z1=6.0,
                       interpolate=True, dznew=0.02):
	if field=='goods':
		area = goods_area
	else:
		area = udf_area
	zdgrid = zdgrid_vdrops(-25., -15., 0.5, -1.0, 2.0, 0.2, z0, z1, 0.1, area)
	zdgrid(simcatalog, field, interpolate=interpolate, dznew=dznew)
	zdgrid.write(picklename)
