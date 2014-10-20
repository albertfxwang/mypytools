#!/usr/bin/env python

import numpy as np
from transfunc import kernelgrid, kernel
from dropout_kernels import KernelFactory, DropoutKernelGrid
from dropout_selection import lbg_colorcrit as lcc
import os, sys
import cPickle
from pygoods import Ftable

simcat_dir = '/Users/khuang/Dropbox/Research/bivariate/bivariate_fit/dropsim_catalogs'
kgrid_dir = '/Users/khuang/Dropbox/Research/bivariate/bivariate_fit/dropsim_kernels'
arcmin_str = 2.9089e-4
area = {}
# area so far only includes GOODS-S
# will add new keys for GOODS-N? Like deep_N, wide_N, etc.
area['udf'] = 11.0 * arcmin_str**2  # = 0.3079e-7
area['goods'] = (320. - 11.) * arcmin_str**2  
area_unit = 'arcmin2'
catalogs = {}
catalogs['goods'] = os.path.join(simcat_dir, 'bvdropsim_run8_nolya_111011.fits')
catalogs['udf'] = os.path.join(simcat_dir, 'bvdropsim_udf_nolya_out.fits')
pixscale = 0.03
zeropoints_goodsv2 = {
   'acs_f435w':25.65288,
   'acs_f606w':26.49341,
   'acs_f775w':25.64053,
   'acs_f850lp':24.84315}

class VdropsKernelFactory(KernelFactory):
   """
   For V-dropouts in GOODS + UDF.
   """
   def __init__(self, catalog, field, interpolate=True, interp_dz=0.02, 
                expand=[0.0, 0.0], 
                SN_lolim={'acs_f850lp':5.0}, SN_hilim={'acs_f435w':5.0},
                mag_in='m1500_in', re_in='re_in_arcsec'):
      z0 = 5.0
      KernelFactory.__init__(self, catalog, z0, 3.5, 6.5, 
                             interpolate=interpolate, expand=expand, 
                             mag_in=mag_in, re_in=re_in,
                             zeropoints=zeropoints_goodsv2)
      self.field = field
      self.lcc = lcc.colorcrit()
      self.lcc = self.lcc('f606w_drop', 'Gia04')
      self.SN_lolim = SN_lolim
      self.SN_hilim = SN_hilim
      self.bands=['acs_f435w','acs_f606w','acs_f775w','acs_f850lp']
      self.area = area[field]
      self.area_unit = area_unit
      self.pixscale = pixscale

   def ColorSelection(self):
      # define color selection
      # Also enforce S/N limits!
      # Require S/N >= sigma_low in detected bands; also require detection
      self.SNcrit = np.ones(len(self.c.d), 'bool')
      if len(self.SN_lolim.keys()) > 0:
         print "sn_lolim", self.SN_lolim
         for b in self.SN_lolim.keys():
            print "Check S/N in %s..." % b
            assert (b in self.bands)
            ston = getattr(self, b + '_sn')
            SN_band = (ston >= self.SN_lolim[b])
            # also guard against S/N == inf or -inf
            SN_band = np.logical_and(SN_band, (ston != -np.inf))
            SN_band = np.logical_and(SN_band, (ston != np.inf))
            self.SNcrit = np.logical_and(self.SNcrit, SN_band)
         for b in self.SN_hilim.keys():
            assert (b in self.bands)
            ston = getattr(self, b + '_sn')
            SN_band = (ston < self.SN_hilim[b])
            # also guard against S/N == inf or -inf
            SN_band = np.logical_and(SN_band, (ston != -np.inf))
            SN_band = np.logical_and(SN_band, (ston != np.inf))
            self.SNcrit = np.logical_and(self.SNcrit, SN_band)
      ## Second: perform color selection!
      print "Do selections..."
      color1 = self.acs_f606w_mag - self.acs_f775w_mag
      color2 = self.acs_f775w_mag - self.acs_f850lp_mag
      self.lcc.select(color1, color2)
      self.lcc.crit = self.lcc.crit & self.SNcrit  # Fold in S/N criteria
      self.lcc.crit = self.lcc.crit & (self.c.detect)  # just in case...
      # guard against weird S/N in F435W
      self.lcc.crit = self.lcc.crit & (self.acs_f435w_sn != -np.inf)
      self.lcc.crit = self.lcc.crit & (self.acs_f435w_sn != np.inf)
      print "Selections done."
      print "Total number of objects in the catalog: %d" % len(self.c.d)
      print "Total number selected as V-dropouts: %d" % (np.sum(self.lcc.crit))
      return self.lcc.crit

   def SelectDropout(self, sn_mode='auto'):
      """
      Main driver method that runs the whole procedure.
      """
      for b in self.bands:
         self.setMagnitudes(b, sn_mode=sn_mode)     
      isDropout = self.ColorSelection()
      return isDropout

class VdropoutKernelGrid(DropoutKernelGrid):
   """
   The top-level class that a user should call.
   """
   def __init__(self, field, kwargs_factory={}, kwargs_kgrid={}):
      Factory = VdropsKernelFactory(catalogs[field], field, **kwargs_factory)
      super(VdropoutKernelGrid, self).__init__(Factory, **kwargs_kgrid)

def make_vdrops_kernel_grid(field, filename, kwargs_factory={'re_in':'re_in_arcsec'}, kwargs_kgrid={'ylimits':[-2.4, 1.0]}):
   """
   A convient function to do everything.
   The keyword arguments **kwargs are passed to VdropoutKernelGrid.__init__().
   """
   Grid = VdropoutKernelGrid(field, kwargs_factory, kwargs_kgrid)
   f = open(filename, 'wb')
   cPickle.dump(Grid, f, 2)
   f.close()
