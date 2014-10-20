#!/usr/bin/env python

import numpy as np
from pygoods import Ftable, sextractor
import KPDFadaptnumpy as KPDF
import cPickle

# for i-dropouts... the simulation catalog needs to have at least F105W
# or F098M fluxes
simcatdir = '/Users/khuang/Dropbox/Research/bivariate/bivariate_fit/simcatalogs/udrops'
limits = np.array([27.0, 30.0])
m0_arr = np.arange(27.0, 30.0, 0.5)
pixdx = 0.02  # pixel scale of the model
xmax = 2.56   # the extent in magnitude of each kernel

class kernel1d(object):
   def __init__(self, x0, x1, pdf, xmax=xmax, pixdx=pixdx, **kwargs):
      self.x0 = x0
      self.x1 = x1
      self.pixdx = pixdx
      for k in kwargs.keys():
         setattr(self, k, kwargs[k])
      self.npix = int(round((x1-x0)/pixdx))
      self.kernel = pdf

class kernelgrid_1d(object):
   def __init__(self, limits, binwidth, pixdx=pixdx, xmax=xmax, **kwargs):
      self.limits = limits
      self.binwidth = binwidth
      self.pixdx = pixdx
      self.xmax = xmax
      for k in kwargs.keys():
         setattr(self, k, kwargs[k])
      self.kernels = {}

def calc_sex_transfunc_1d(simcat, kgridname, limits=limits, field='udf', 
                          logre_lims=[0.0,1.5], dm=0.5):
   """
   Calculate the magnitude transfer functions for each apparent magnitude bin
   for the faint-end of i-dropouts.
   logre_lims are the limits of log(Re) (in pixels) that are used from 
   the simulation catalog. This is to match the observed range of log(Re)
   of i-dropouts at the faint-end in UDF.
   """
   # Read simulation catalog
   c = Ftable(simcatdir+'/'+simcat)
   y_mag_auto = c.y_mag_auto
   y_mag_in = c.y_mag_in
   kw = {'logre_lims':logre_lims, 'm0_arr':m0_arr}
   kgrid = kernelgrid_1d(limits, dm, pixdx=pixdx, xmax=xmax, **kw)
   kgrid.grid1d = np.arange(-xmax/2., xmax/2.+pixdx, pixdx)
   logre_crit = (np.log10(c.y_flux_radius_1)>=logre_lims[0]) & \
                (np.log10(c.y_flux_radius_1)<logre_lims[1]) & \
                (c.detect==True)

   for m in m0_arr:
      print "m=", m, m+kgrid.binwidth
      bincrit = (y_mag_in>=m) & (y_mag_in<(m+kgrid.binwidth))
      y_mag_in_bin = y_mag_in[(bincrit==True)&(logre_crit==True)]
      y_mag_auto_bin = y_mag_auto[(bincrit==True)&(logre_crit==True)]
      dmag = y_mag_auto_bin - y_mag_in_bin
      h = KPDF.UPDFOptimumBandwidth(dmag)
      pdf = KPDF.UPDFEpanechnikov(dmag, kgrid.grid1d, h)
      pdf = pdf / sum(pdf)  # normalize to 1.0
      kgrid.kernels['%.1f'%m] = kernel1d(m, m+kgrid.binwidth, pdf, 
                                         xmax=xmax, pixdx=pixdx)

   # Write kgrid into a pickled file
   f = open(kgridname, 'wb')
   cPickle.dump(kgrid, f, 2)
   f.close()

