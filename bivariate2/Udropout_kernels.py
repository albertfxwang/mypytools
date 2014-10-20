#!/usr/bin/env python

import numpy as np
from transfunc import kernelgrid, kernel
from dropout_kernels import KernelFactory, DropoutKernelGrid
from dropout_selection import lbg_colorcrit as lcc
import os, sys
import udropsim_uflux as uu
import cPickle
from pygoods import Ftable

simcat_dir = '/Users/khuang/Dropbox/Research/bivariate/bivariate_fit/dropsim_catalogs'
kgrid_dir = '/Users/khuang/Dropbox/Research/bivariate/bivariate_fit/dropsim_kernels'
tfitsim_cat = {'udf': os.path.join(simcat_dir, 'run1_udf_130625.fits'),
               'deep': os.path.join(simcat_dir, 'run2_tfitsim_130322.fits')}
tfitsim_cat['ers'] = tfitsim_cat['deep']
tfitsim_cat['wide'] = tfitsim_cat['deep']
arcmin_str = 2.9089e-4
area = {}
# area so far only includes GOODS-S
# will add new keys for GOODS-N? Like deep_N, wide_N, etc.
area['udf'] = 4.497 * arcmin_str**2  # = 3.805226134736999e-07
area['deep'] = 67.318 * arcmin_str**2  # = 5.696246674187799e-06
area['ers'] = 46.359 * arcmin_str**2  # = 3.924536093598e-06
area['wide'] = 50.867 * arcmin_str**2  # = 4.304212537150699e-06
area_unit = 'arcmin2'
catalogs = {}
## --------------- BELOW ARE PSF-MATCHED SIMULATION CATALOGS --------------- ##
catalogs['deep'] = os.path.join(simcat_dir, 'udropssim_run4m_deep_140704.fits')
catalogs['udf'] = os.path.join(simcat_dir, 'udropssim_run4m_udf_140801.fits')
catalogs['ers'] = os.path.join(simcat_dir, 'udropssim_run4m_ers_140721.fits')
catalogs['wide'] = os.path.join(simcat_dir, 'udropssim_run4m_wide_140721.fits')
catalogs_ubz = {}
catalogs_ubz['udf'] = '/Users/khuang/CANDELS/goodss/udrops_sim/udf/udrops_run5m_udf_140801.fits'
## ------------------------------------------------------------------------- ##
pixscale = 0.06  # the pixel scale of the images used in the simulations
umag_kpdf = {}
umag_kpdf['deep'] = os.path.join(kgrid_dir, 'uvimos_gds_mag_pdf_140310.p')
umag_kpdf['udf'] = os.path.join(kgrid_dir, 'uvimos_hudf_mag_pdf_140310.p')
umag_kpdf['ers'] = umag_kpdf['deep']
umag_kpdf['wide'] = umag_kpdf['deep']
umag_1sig_kpdf = {}
umag_1sig_kpdf['deep'] = os.path.join(kgrid_dir, 'uvimos_gds_mag_1sig_pdf_140310.p')
umag_1sig_kpdf['udf'] = os.path.join(kgrid_dir, 'uvimos_hudf_mag_1sig_pdf_140310.p')
umag_1sig_kpdf['ers'] = umag_1sig_kpdf['deep']
umag_1sig_kpdf['wide'] = umag_1sig_kpdf['deep']
# color-correction because I didn't run the sims with PSF-matched images
# however, the corrections are small... not sure how much difference it makes
# ISO_COLOR + color_corr ~= INPUT_COLOR
# the corrections are calculated from simulation catalogs and applies to HST
# colors only.
vmy_corr = {}
# vmy_corr['udf'] = 0.025
# vmy_corr['deep'] = 0.050
# vmy_corr['ers'] = 0.040
# vmy_corr['wide'] = 0.052
vmy_corr['udf'] = 0.
vmy_corr['deep'] = 0.
vmy_corr['ers'] = 0.
vmy_corr['wide'] = 0.

class UdropsKernelFactory(KernelFactory):
   """
   For U-dropouts in GOODS-S (VIMOS U-band).
   """
   # This is for UBVY selection. Write another class for UBY selection if 
   # necessary.
   def __init__(self, catalog, field, droplabel="UBVY105", interpolate=True, interp_dz=0.02, n_repeat=5, expand=[0.0, 0.0], mag0=22.0, SN_lolim={'wfc3_f160w':5.0,'wfc3_f105w':5.0,'acs_f435w':3.0,'acs_f606w':5.0}, mag_in='m1500_in', re_in='re_in_arcsec', hstBands=['acs_f435w','acs_f606w','wfc3_f105w','wfc3_f160w'],uband='vimos_u'):
      z0 = 3.5
      KernelFactory.__init__(self, catalog, z0, 2.5, 4.5, 
                             interpolate=interpolate, n_repeat=n_repeat,
                             expand=expand, mag_in=mag_in, re_in=re_in)
      self.field = field
      self.lcc = lcc.colorcrit()
      self.SN_lolim = SN_lolim
      print "droplabel: ", droplabel
      self.lcc = self.lcc('uvimos', droplabel)
      print "Eval_string:", self.lcc.eval
      # if self.field.lower() == 'ers':
      #    if droplabel == None:
      #       droplabel = 'UBVY098'
      #    self.lcc = self.lcc('uvimos', droplabel)
      #    # self.bands=['acs_f435w','acs_f606w','wfc3_f098m','wfc3_f160w'] 
      # else:
      #    if droplabel == None:
      #       droplabel = 'UBVY105'
      #    self.lcc = self.lcc('uvimos', droplabel)
         # self.bands=['acs_f435w','acs_f606w','wfc3_f105w','wfc3_f160w']
      self.bands = hstBands
      self.uband = uband
      self.area = area[field]
      self.area_unit = area_unit
      self.pixscale = pixscale
      self.umag_kpdf = umag_kpdf[field]
      self.umag_1sig_kpdf = umag_1sig_kpdf[field]
      self._drawn_umag = False
      self.mag0 = mag0  
      # the lower limit of input U-band magnitude from simulations

   def testdrawMag(self):
      """
      For tests.
      """
      self.vimos_u_mag = self.c.u_mag_in
      self.vimos_u_mag_1sig = self.c.u_mag_in

   def drawUMag(self):
      # Draw U-band magnitude from probability distribution functions calculated
      # using TFIT simulation catalogs in VIMOS U-band.
      # In princple, we draw an output U-band magnitude given U-band input 
      # magnitudes to simulate photometric errors. We do it this way because
      # running enough simulations with TFIT is prohibitively expensive.
      cu = Ftable(tfitsim_cat[self.field])
      # mag0 = 22.0
      umag_in = cu.uvimos_mag_in[cu.uvimos_fitquerr > 0.]
      print "Repeat the drawing process %d times" % self.n_repeat
      vimos_u_mag = uu.draw_umag(np.tile(self.c.u_mag_in, self.n_repeat),
                    self.umag_kpdf, mag0=self.mag0)  
      vimos_u_mag_1sig = uu.draw_ulimmag(self.n_repeat*len(self.c.u_mag_in),
                         self.umag_1sig_kpdf)
      # If the drawn U-band 1-sigma limit mag is brighter than the drawn 
      # U-band magnitude, use the 1-sigma limit magnitude as the measured 
      # U-band magnitude.
      self._drawn_umag = True
      return np.minimum(vimos_u_mag, vimos_u_mag_1sig)

   def ColorSelection(self, vmy_corr=0.0, full_output=False):
      # define color selection
      # Also enforce S/N limits!
      # Require S/N >= sigma_low in detected bands; also require detection
      self.SNcrit = np.ones(len(self.c.d) * self.n_repeat, 'bool') & np.tile(self.c.detect, self.n_repeat)
      if len(self.SN_lolim.keys()) > 0:
         print "sn_lolim", self.SN_lolim
         for b in self.SN_lolim.keys():
            assert (b in self.bands)
            SN_band = (getattr(self, b + '_sn') >= self.SN_lolim[b])
            SN_band = np.tile(SN_band, self.n_repeat)
            self.SNcrit = np.logical_and(self.SNcrit, SN_band)
      ## Second: perform color selection!
      print "Do selections..."
      umag = getattr(self, '%s_mag' % self.uband)
      m1 = getattr(self, '%s_mag' % self.bands[0])
      color1 = umag - m1
      m2 = getattr(self, '%s_mag' % self.bands[1])
      # m3 = getattr(self, '%s_mag' % self.bands[2])
      if len(self.lcc.bands) == 3:
         # color1 = m1 - m2
         # color2 = m2 - m3
         color2 = m1 - m2
      else:
         m3 = getattr(self, '%s_mag' % self.bands[2])
         color2 = m2 - m3
         # m4 = getattr(self, '%s_mag' % self.bands[3])
         # color1 = self.vimos_u_mag-self.acs_f435w_mag
         # color1 = m1 - m2
         # color2 = m3 - m4
      # if self.field.lower() == 'ers':
      #    color2 = self.acs_f606w_mag-self.wfc3_f098m_mag
      # else:
      #    color2 = self.acs_f606w_mag-self.wfc3_f105w_mag
      # applies color-correction
      color2 = color2 + vmy_corr
      self.lcc.select(color1, color2)
      colorcrit = self.lcc.crit.copy()
      self.lcc.crit = self.lcc.crit & self.SNcrit  # Fold in S/N criteria
      # isDropout = self.lcc.crit
      print "Selections done."
      print "Total number of objects in the catalog: %d" % len(self.c.d)
      print "Total number of objects satisfying color cuts: %d" % (np.sum(colorcrit)/self.n_repeat)
      print "Total number of objects satisfying S/N cuts: %d" % (np.sum(self.SNcrit)/self.n_repeat)
      print "Total number selected as U-dropouts: %d" % (np.sum(self.lcc.crit)/(self.n_repeat))
      if full_output:
         return self.lcc.crit, color1, color2
      else:
         return self.lcc.crit
      
   def SelectDropout(self, full_output=False):
      """
      Main driver method that runs the whole procedure.
      """
      if not self._drawn_umag:
         for b in self.bands:
            self.setMagnitudes(b)     
         # self.vimos_u_mag = self.drawUMag()
         setattr(self, '%s_mag' % self.uband, self.drawUMag())
      isDropout = self.ColorSelection(full_output=full_output)
      return isDropout

class UdropsBoxcarKernelFactory(UdropsKernelFactory):
   def SelectDropout(self):
      """
      Return a boxcar P(z).
      """
      z0 = 3.5
      zlo = z0 - 0.5
      zhi = z0 + 0.5
      isDropout = np.where((self.c.redshift>=zlo) & (self.c.redshift<zhi),
                           True, False)
      isDropout = isDropout & (self.c.detect)
      return isDropout


class UdropoutKernelGrid(DropoutKernelGrid):
   """
   The top-level class that a user should call.
   """
   def __init__(self, field, catalogs=catalogs, droplabel="UBVY105", hstBands=['acs_f435w','acs_f606w','wfc3_f105w','wfc3_f160w'], kwargs_factory={}, kwargs_kgrid={}):
      if field == 'udf':
         kwargs_factory['mag0'] = 25.0
      else:
         kwargs_factory['mag0'] = 23.0
      Factory = UdropsBoxcarKernelFactory(catalogs[field], field, 
                                          hstBands=hstBands)
      super(UdropoutKernelGrid, self).__init__(Factory, **kwargs_kgrid)

class UdropoutBoxcarKernelGrid(DropoutKernelGrid):
   """
   A boxcar dropout kernel grid (P(z) = 1 between z=3 and z=4) for testing.
   """
   def __init__(self, field, n_repeat=1, hstBands=['acs_f435w','acs_f606w','wfc3_f105w','wfc3_f160w'], kwargs_factory={}, kwargs_kgrid={}):
      # if field == 'udf':
      #    kwargs_factory['mag0'] = 25.0
      # else:
      #    kwargs_factory['mag0'] = 23.0
      Factory = UdropsBoxcarKernelFactory(catalogs[field], field, n_repeat=1,
                                          droplabel='UBVY105', 
                                          interpolate=True,
                                          hstBands=hstBands, **kwargs_factory)
      print "Color criteria: ", Factory.lcc.criteria
      # super(UdropoutKernelGrid, self).__init__(Factory, **kwargs_kgrid)
      DropoutKernelGrid.__init__(self, Factory, **kwargs_kgrid)


def make_udrops_ubvy_kernel_grid(field, filename, droplabel="UBVY105",kwargs_factory={'n_repeat':20, 're_in':'re_in_arcsec'}, kwargs_kgrid={'ylimits':[-2.4, 1.0]}):
   """
   A convient function to do everything.
   The keyword arguments **kwargs are passed to UdropoutKernelGrid.__init__().
   """
   if field == 'ers':
      hstBands = ['acs_f435w','acs_f606w','wfc3_f098m','wfc3_f160w']
      kwargs_factory['SN_lolim'] = {'wfc3_f160w':5.0,'wfc3_f098m':5.0,'acs_f435w':3.0,'acs_f606w':3.0}
   else:
      hstBands = ['acs_f435w','acs_f606w','wfc3_f105w','wfc3_f160w']
      kwargs_factory['SN_lolim'] = {'wfc3_f160w':5.0,'wfc3_f105w':5.0,'acs_f606w':3.0,'acs_f435w':3.0}
   # For ERS, need to also supply a custom SN_lolim argument.
   Grid = UdropoutKernelGrid(field, droplabel=droplabel, 
                             catalogs=catalogs, hstBands=hstBands,
                             kwargs_factory=kwargs_factory, 
                             kwargs_kgrid=kwargs_kgrid)
   f = open(filename, 'wb')
   cPickle.dump(Grid, f, 2)
   f.close()

def make_udrops_ubz_kernel_grid(field, filename, droplabel='W14',kwargs_factory={'n_repeat':20, 're_in':'re_in_arcsec'}, kwargs_kgrid={'ylimits':[-2.5,1.0]}):
   hstBands = ['acs_f435w','acs_f850lp','wfc3_f160w']
   kwargs_factory['SN_lolim'] = {'wfc3_f160w':5.0,'acs_f435w':3.0,'acs_f850lp':5.0}
   Grid = UdropoutKernelGrid(field, droplabel=droplabel, 
                             catalogs=catalogs_ubz, hstBands=hstBands,
                             kwargs_factory=kwargs_factory, 
                             kwargs_kgrid=kwargs_kgrid)
   f = open(filename, 'wb')
   cPickle.dump(Grid, f, 2)
   f.close()
