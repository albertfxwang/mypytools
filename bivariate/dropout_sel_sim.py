#!/usr/bin/env python

from numpy import *
from pygoods import *
import cPickle
import dropout_colorsel as dcs
import bdrops_sel as bs
import vdrops_sel as vs
import bdrops_sel_udf as bsu
import vdrops_sel_udf as vsu
from dropout_selection import select_dropouts as sd
import select_udrops_ubvy as suu

# Note that in each M-bin,
# completeness = (# of objects selected by color criteria) / 
#                (# of objects in the target z range)
# Note that all numbers here are the numbers detected and 
# measured by SE, so the completeness here does NOT include
# SE detection completeness

zpt_goods = {'b':25.65288, 'v':26.49341, 'i':25.64053, 'z':24.84315}
zpt_udf = {'b':25.673, 'v':26.486, 'i':25.654, 'z':24.862}

class catalog(object):
	def __init__(self, c, band):
		# initialize a catalog instance with the attributes being different columns of catalog c
		self.number = c.number
		self.x_input = c.x_input
		self.y_input = c.y_input
		self.type_input = c.type_input
		self.re_input = c.re_input
		self.ell_input = c.ell_input
		self.theta_input = c.theta_input
		self.z_input = c.z_input
		self.ebmv_input = c.ebmv_input
		self.m1500_input = c.m1500_input
		self.detect = c.detect
		self.class_star = c.class_star
		self.mag_auto = c.__getitem__('%s_mag_auto'%band)
		self.magerr_auto = c.__getitem__('%s_magerr_auto'%band)
		self.flux_auto = c.__getitem__('%s_flux_auto'%band)
		self.fluxerr_auto = c.__getitem__('%s_fluxerr_auto'%band)
		self.mag_iso = c.__getitem__('%s_mag_iso'%band)
		self.magerr_iso = c.__getitem__('%s_magerr_iso'%band)
		self.flux_iso = c.__getitem__('%s_flux_iso'%band)
		self.fluxerr_iso = c.__getitem__('%s_fluxerr_iso'%band)
		self.magin = c.__getitem__('%s_magin'%band)
		self.flux_radius_1 = c.__getitem__('%s_flux_radius_1'%band)
		self.imaflags_iso = c.__getitem__('%s_imaflags_iso'%band)


def bdrops_sel_sim(c,cs_thresh=0.95,field='goods',imaflags_thresh=4):
   """Does B-dropout selection from SExtractor simulation
      c = sextractor('lbgsim_out.cat')
   """
   # class_star threshold down to 0.5 mag brighter than z-band magnitude limit
   cb = catalog(c, 'b')
   cv = catalog(c, 'v')
   ci = catalog(c, 'i')
   cz = catalog(c, 'z')
   if field == 'goods':
      zmagcut = 26.5
      selfunc = bs.dropsel_magiso
   elif field == 'udf':
      zmagcut = 28.5
      selfunc = bsu.dropsel_magiso
   
   dropcrit = selfunc(cb,cv,ci,cz,cs_thresh=cs_thresh,zmagcut=zmagcut,imaflags_thresh=imaflags_thresh)[0]
   
   #cs_crit = (cz.class_star<=cs_thresh)|(cz.mag_auto>=(zmagcut-0.5))
   # image flags
   #flag_crit = (cv.imaflags_iso<=imaflags_thresh)&(cz.imaflags_iso<=imaflags_thresh)
   # z magnitude cut
   #zmag_crit = (cz.mag_auto <= zmagcut)
   # S/N cut - clarification on how to calculate S/N?
   #zston_crit = ((1.0857/cz.magerr_auto)>=5.)
   # color selection
   #bmags = c.b_mag_iso.copy()
   #vmags = c.v_mag_iso.copy()
   #zmags = c.z_mag_iso.copy()
   # Use 1-sigma as upper limits in B-band if S/N <= 1.0
   #bmags = dcs.flux_uplim(c.b_flux_iso, c.b_fluxerr_iso, acszpt['b'], bmags)
   #vmags = dcs.flux_uplim(c.v_flux_iso, c.v_fluxerr_iso, acszpt['v'], vmags)
   #for i in range(len(c)):
   #   if c.b_flux_iso[i]/c.b_fluxerr_iso[i] <= 1.0:
   #      if c.b_fluxerr_iso[i] > 0.:
   #         bmags[i] = flux2mag(c.b_fluxerr_iso[i],25.65288)
   #      else: bmags[i] = 99.0
   #   if c.v_flux_iso[i]/c.v_fluxerr_iso[i] <= 1.0:
   #      if c.v_fluxerr_iso[i] > 0.:
   #         vmags[i] = flux2mag(c.v_fluxerr_iso[i],26.49341)
   #      else: vmags[i] = 99.0
   #color_crit = dcs.bdrops_colorsel_g04(bmags, vmags, zmags)
  
   #dropcrit = cs_crit & color_crit & zston_crit & flag_crit & zmag_crit

   return dropcrit, cb.mag_iso.copy(), cv.mag_iso.copy(), cz.mag_iso.copy()


def write_bdrops(c,bdropscat):
   crit = bdrops_sel(c)
   f = open(bdropscat,'w')
   f.write(c._header)
   for i in range(len(c)):
      if crit[i]:
         for j in range(len(c._d)):
            f.write('%s ' % c._colentries[i][j])
         f.write('\n')
   f.close()
   return 0 


def vdrops_sel_sim(c,cs_thresh=0.95, field='goods', imaflags_thresh=4, bstonlim=2.0):
   """Does V-dropout selection from SExtractor simulation
      c = sextractor('lbgsim_out.cat')
   """
   cb = catalog(c, 'b')
   cv = catalog(c, 'v')
   ci = catalog(c, 'i')
   cz = catalog(c, 'z')
   if field == 'goods':
      zmagcut = 26.5
      selfunc = vs.dropsel_magiso
      zpt_dic = zpt_goods
   elif field == 'udf':
      zmagcut = 28.5
      selfunc = vsu.dropsel_magiso
      zpt_dic = zpt_udf
      
   dropcrit = selfunc(cb,cv,ci,cz,cs_thresh=cs_thresh,zmagcut=zmagcut,imaflags_thresh=imaflags_thresh,
   	bstonlim=bstonlim)[0]
   # class_star threshold down to zmag = 26.0 only
   #cs_crit = (c.class_star<=cs_thresh)|(c.z_mag_auto>=(zmagcut-0.5))
   # z-band magnitude cut
   #zmag_crit = (c.z_mag_auto<=zmagcut)
   # image flags
   #flag_crit = (c.v_imaflags_iso<=imaflags_thresh) & (c.i_imaflags_iso<=imaflags_thresh)\
   #   & (c.z_imaflags_iso<=imaflags_thresh)
   # S/N cut - clarification on how to calculate S/N?
   #zston_crit = ((c.z_flux_auto/c.z_fluxerr_auto) >= 5.)
   #bston_crit = ((c.b_flux_auto/c.b_fluxerr_auto) <= bstonlim)
   # completeness is sensitive to this threshold
   # color selection
   #vmags = c.v_mag_iso.copy()
   #imags = c.i_mag_iso.copy()
   #zmags = c.z_mag_iso.copy()
   # Use 1-sigma as upper limits in V-band and i-band if S/N <= 1.0
   #vmags = dcs.flux_uplim(c.v_flux_iso, c.v_fluxerr_iso, acszpt['v'], vmags)
   #imags = dcs.flux_uplim(c.i_flux_iso, c.i_fluxerr_iso, acszpt['i'], imags)
   #for i in range(len(c)):
   #   if c.v_flux_iso[i]/c.v_fluxerr_iso[i] <= 1.0:
   #      if c.v_fluxerr_iso[i] > 0.:
   #         vmags[i] = flux2mag(c.v_fluxerr_iso[i],26.49341)
   #      else: vmags[i] = 99.0
   #   if c.i_flux_iso[i]/c.i_fluxerr_iso[i] <= 1.0:
   #      if c.i_fluxerr_iso[i] > 0.:
   #         imags[i] = flux2mag(c.i_fluxerr_iso[i],25.64053)
   #      else: imags[i] = 99.0

   #color_crit = dcs.vdrops_colorsel_g04(vmags, imags, zmags)
  
   #dropcrit = cs_crit & zmag_crit & color_crit & zston_crit & bston_crit & flag_crit

   return dropcrit, cv.mag_iso.copy(), ci.mag_iso.copy(), cz.mag_iso.copy()

def udrops_sel_sim(catalog, cs_thresh=0.95, field='deep', imaflags_thresh=4,
   sn_lolim={'f160w':5.0,'f435w':3.0,'f105w':5.0}):
   ### Selects U-dropouts from simulation catalogs.
   # First, read the simulation catalog --- from a big multi-wavelength catalog
   # IMPORTANT: assume that the U-band flux has already been drawn and are included in the catalog
   ### TO DO: write scripts to draw U-band fluxes and then come back!!
   if field=='ers':
      # Select using F098M instead of F105W
      udrops = suu.udrops_ubvy098(catalog=catalog,sn_lolim=sn_lolim)
   else:
      udrops = suu.udrops_ubvy105(catalog=catalog,sn_lolim=sn_lolim)


   pass
