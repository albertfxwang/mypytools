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
area['udf'] = 5.135 * arcmin_str**2  # = 4.345e-07
area['deep'] = 69.33 * arcmin_str**2  # = 5.866e-06
area['ers'] = 46.38 * arcmin_str**2  # = 3.925e-06
area['wide'] = 48.19 * arcmin_str**2  # = 4.078e-06
area_unit = 'arcmin2'
catalogs = {}
catalogs['deep'] = os.path.join(simcat_dir, 'idropsim_goodss_deep_130529.fits')
catalogs['udf'] = os.path.join(simcat_dir, 'idropsim_goodss_udf_130519.fits')
catalogs['ers'] = catalogs['deep']
catalogs['wide'] = os.path.join(simcat_dir, 'idropsim_goodss_wide_130602.fits')
pixscale = 0.06  # the pixel scale of the images used in the simulations

class idropsKernelFactory(KernelFactory):
   """
   For i-dropouts in GOODS.
   """
   def __init__(self, catalog, field, interpolate=True, interp_dz=0.02, 
               expand=[0.0, 0.0],
               SN_lolim={'wfc3_f160w':5.0,'wfc3_f125w':3.5,'acs_f850lp':2.0},
               SN_hilim={'acs_f435w':3.0},
               mag_in='m1500_in', re_in='re_in_arcsec'):
      z0 = 6.0
      KernelFactory.__init__(self, catalog, z0, 4.5, 7.5, 
                             interpolate=interpolate,
                             expand=expand, mag_in=mag_in, re_in=re_in)
      self.field = field
      self.lcc = lcc.colorcrit()
      self.SN_lolim = SN_lolim
      self.SN_hilim = SN_hilim
      self.lcc = self.lcc('f775w_drop', 'Hua13')
      self.bands = ['acs_f435w', 'acs_f606w', 'acs_f775w', 'acs_f850lp', 'wfc3_f125w', 'wfc3_f160w']
      for b in self.bands:
         self.setMagnitudes(b)
      self.area = area[field]
      self.area_unit = area_unit
      self.pixscale = pixscale

   def ColorSelection(self):
      # define color selection
      # Also enforce S/N limits!
      # Require S/N >= sigma_low in detected bands; also require detection
      self.SNcrit = np.logical_and(np.ones(len(self.c.d), 'bool'), self.c.detect)
      if len(self.SN_lolim.keys()) > 0:
         print "sn_lolim", self.SN_lolim
         for b in self.SN_lolim.keys():
            assert (b in self.bands)
            SN_band = (getattr(self, b + '_sn') >= self.SN_lolim[b])
            self.SNcrit = np.logical_and(self.SNcrit, SN_band)
      if len(self.SN_hilim.keys()) > 0:
         print "sn_hilim", self.SN_hilim
         for b in self.SN_hilim.keys():
            assert (b in self.bands)
            SN_band = (getattr(self, b + '_sn') < self.SN_hilim[b])
            self.SNcrit = np.logical_and(self.SNcrit, SN_band)
      ## Second: perform color selection!
      print "Do selections..."
      # first, the V-i color OR V-band S/N limit
      # Note that I use different apertures for V-band S/N
      flux_f606w = self.c.v_flux_aper_2
      fluxerr_f606w = self.c.v_fluxerr_aper_2
      self.acs_f606w_sn = np.where(flux_f606w > 0, flux_f606w/fluxerr_f606w, 0.)
      colorcrit1 = np.logical_or((self.acs_f606w_mag-self.acs_f775w_mag)>2.8, 
                                    self.acs_f606w_sn<5.0)
      colorcrit1 = np.logical_or(colorcrit1, self.wfc3_f125w_mag<25.5)
      jmagcrit = (self.wfc3_f125w_mag<25.5)
      # colorcrit2 = np.logical_or((self.wfc3_f125w_mag<=26.), self.acs_f606w_sn<5.0)
      color1 = self.acs_f775w_mag - self.acs_f850lp_mag
      color2 = self.acs_f850lp_mag - self.wfc3_f125w_mag
      self.lcc.select(color1, color2)
      # colorcrit = self.lcc.crit.copy()
      self.lcc.crit = np.logical_and(self.lcc.crit, self.SNcrit)  # Fold in S/N criteria
      # self.lcc.crit = np.logical_and(self.lcc.crit, colorcrit1)  # Fold in V-i criterion
      # self.lcc.crit = np.logical_and(self.lcc.crit, colorcrit2)  
      # self.lcc.crit = np.logical_or(self.lcc.crit, jmagcrit)
      # Fold in F125W mag cut or S/N cut in F606W
      
      # isDropout = self.lcc.crit
      print "Selections done."
      print "Total number of objects in the catalog: %d" % len(self.c.d)
      print "Total number selected as i-dropouts: %d" % (np.sum(self.lcc.crit))
      return self.lcc.crit

   def SelectDropout(self):
      """
      Main driver method that runs the whole procedure.
      """
      isDropout = self.ColorSelection()
      return isDropout

class idropoutKernelGrid(DropoutKernelGrid):
   """
   The top-level class that a user should call.
   """
   def __init__(self, field, kwargs_factory={}, kwargs_kgrid={}):
      Factory = idropsKernelFactory(catalogs[field], field, **kwargs_factory)
      super(idropoutKernelGrid, self).__init__(Factory, **kwargs_kgrid)

def make_idrops_kernel_grid(field, filename, kwargs_factory={'re_in':'re_in_arcsec'}, 
                            kwargs_kgrid={'ylimits':[-2.4, 1.0]}):
   """
   A convient function to do everything.
   The keyword arguments **kwargs are passed to UdropoutKernelGrid.__init__().
   """
   # For ERS, need to also supply a custom SN_lolim argument.
   Grid = idropoutKernelGrid(field, kwargs_factory, kwargs_kgrid)
   Grid.filename = filename
   f = open(filename, 'wb')
   cPickle.dump(Grid, f, 2)
   f.close()