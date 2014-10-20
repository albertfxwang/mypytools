#!/usr/bin/env python

import numpy as np
from BivDist import BivariateDistribution
from pygoods import Ftable, sextractor
import KPDFadaptnumpy as KPDF


class InterloperRLDist(BivariateDistribution):
   """
   A class for the size-magnitude distribution of LBG interlopers. Note that 
   instead of giving the *number density* in the distribution, it gives the 
   **fraction** of number density that are interlopers. Therefore, to get the
   total number density of real LBGs + interlopers, one needs to MULTIPLY the 
   distribution describing the real LBGs with an InterloperRLDist instance.
   """
   def __init__(self, maglimits, logrlimits, dmag, dlogr, field,
                mag_lolim=23.0, yunit='arcsec'):
      """
      model0: a bivmodel instance of the uncorrected RL distribution.
      sdfile: a pickled file containing the sizedist class instance.
      mag_lolim: the bright limit over which do not include interlopers,
                 because we can identify interlopers at the bright end.
      """
      super(InterloperRLDist, self).__init__(maglimits, logrlimits, 'magnitude',
                                           'logR', 'mag', 'arcsec', dmag, dlogr)
      self.field = field
      self.mag_lolim = mag_lolim
      # can use self.field to make sure that the interloper fraction is consistent
      # with the RL distribution being multiplied.

   def compute(self, intfrac, IntSizeDist):
      """
      IntSizeDist: a size-mag distribution (uniform in the mag direction) of the
                   size distribution for interlopers.
      intfrac: a FITS table with interloper fractions.
      """
      assert self.field == IntSizeDist.field, "Make sure you are using the distributions for the same field!"
      dm_intfrac = intfrac.mag0[1] - intfrac.mag0[0]
      # iterate through all magnitude bins for interloper fractions
      for i in range(len(intfrac.mag0)):
         m0 = intfrac.mag0[i]
         m1 = intfrac.mag0[i] + dm_intfrac
         if not np.logical_and(m0>=self.xlimits[0], m0<self.xlimits[1]):
            continue
         if (m1) < self.mag_lolim: 
            continue
         ix0 = np.searchsorted(self.xedges(), m0) 
         ix0 = np.minimum(np.maximum(ix0, 0), self.Nx)  # enforce range of ix0
         ix1 = np.searchsorted(self.xedges(), m1) 
         ix1 = np.minimum(np.maximum(ix1, 0), self.Nx)  # enforce range of ix1
         # print ix0, ix1
         assert (ix1 > ix0)
         self.value[ix0:ix1,:] = intfrac.interloper_frac_drops[i]
      # Then, multiply by the size distribution
      self.value = self.value * IntSizeDist.value
      return self.value

class InterloperSizeDist(BivariateDistribution): 
   """
   A class to store the size distribution (from SExtractor) for interlopers,
   which is in fact calculated using ALL objects in the field.
   Directly copied from the first incarnation of interlopers.py.
   NOTE: convert sizes to arcsec (so it eliminates the dependence on pixel scale).
   """
   def __init__(self, maglimits, logrlimits, dmag, dlogr, dmag_bin, dlogr_bin,
                field, yunit='arcsec'):
      super(InterloperSizeDist, self).__init__(maglimits, logrlimits, 'magnitude',
                                           'logR', 'mag', 'arcsec', dmag, dlogr)
      self.dx_bin = dmag_bin
      self.dy_bin = dlogr_bin
      self.Ny_pix_bin = int(round(dlogr_bin / dlogr))
      Nx_bins = (self.xlimits[1] - self.xlimits[0]) / self.dx_bin
      Nx_bins = int(round(Nx_bins))
      Ny_bins = (self.ylimits[1] - self.ylimits[0]) / self.dy_bin
      Ny_bins = int(round(Ny_bins))
      self.mbins = np.linspace(self.xlimits[0], self.xlimits[1], num=Nx_bins+1)
      self.logrbins = np.linspace(self.ylimits[0], self.ylimits[1], num=Ny_bins+1)
      self.field = field

   def compute(self, mag_array, size_array, h=0):
      # mag_array: array of all magnitudes in the catalog cs
      # size_array: array of all sizes in the catalog cs (in pixels)
      # NOTE: size_array should be in units of arcsec
      print "Please make sure that size_array has the unit of arcsec."
      size_array = size_array[size_array>10.**(self.ylimits[0])]
      mag_array = mag_array[size_array>10.**(self.ylimits[0])]
      logr_array = np.log10(size_array)
      # logr_bins is the edges of logR *bins*, not the pixel coordinates
      logr_bins = np.concatenate([self.logrbins, [self.logrbins[-1]+self.dy_bin]])
      logr_pix = self.yedges()[:-1] + self.dy / 2. # the pixel centers
      for i in range(len(self.mbins) - 1):
         # Now calculate the size distribution of possible interlopers in 
         # every magnitude bin. We use the size distribution of ALL sources in 
         # the field as an approximation for the size distribution of possible
         # low-z interlopers, because we do not know exactly which objects are
         # low-z interlopers (and if we knew, we would not have needed to do 
         # this in the first place.)
         # After we have the array of galaxy sizes in this magnitude bin, we 
         # use kernel density estimator to calculate the size distribution. 
         # We also enforce the size distribution to be >= 0, and normalize the 
         # size distribution so that sum(sizedist) * dy = 1.0.
         bincrit = np.logical_and(mag_array>=self.mbins[i], mag_array<(self.mbins[i]+self.dx_bin))
         print "Number of sources between [%.1f, %.1f]: %d" % (self.mbins[i], self.mbins[i]+self.dx_bin, bincrit.sum())
         logr_thisbin = logr_array[bincrit==True]
         ix0 = np.searchsorted(self.xedges(), self.mbins[i])
         ix0 = np.maximum(ix0, 0)
         ix1 = np.searchsorted(self.xedges(), self.mbins[i]+self.dx_bin)
         ix1 = np.minimum(ix1, self.Nx)
         if h==None:
            # No kernel density estimation, just a straight histogram
            sizedist = np.histogram(logr_thisbin, logr_bins)[0]
            sizedist = np.repeat(sizedist, int(round(self.dy_bin/self.dy)))
            # Now normalize the size distribution
            sizedist = sizedist.astype('float') / (sizedist.sum() * self.dy)
            # should also factor in the pixel size in logR
            print "len(logr_thisbin)", len(logr_thisbin)
         else:
            # use Kernel density estimator 
            if h <= 0:
               h = KPDF.UPDFOptimumBandwidth(logr_thisbin)
            # evalulate the PDF at the center of each pixel
            p = KPDF.UPDFEpanechnikov(logr_thisbin, logr_pix, h)
            sizedist = p / (p.sum() * self.dy)
            # sizedist = sizedist.repeat(self.Ny_pix_bin)
         self.value[ix0:ix1,:] = sizedist
      return self.value

   def get_sizedist(self, mag):
      # return a size distribution at the input magnitude
      ix = np.searchsorted(self.xedges(), mag) - 1
      ix = np.maximum(ix, 0)
      ix = np.minimum(ix, self.Nx)
      return self.value[ix,:]

def make_interloper_sizedist(catalog, maglimits, logrlimits, field):
   # from a catalog with columns of object magnitudes and sizes (in arcsec), 
   # calculate the size distribution. Assume that the column names are MAG and 
   # RE_ARCSEC.
   c = sextractor(catalog)
   SD = InterloperSizeDist(maglimits, logrlimits, 0.02, 0.02, 0.5, 0.2, field)
   x = SD.compute(c.mag, c.re_arcsec)
   return SD