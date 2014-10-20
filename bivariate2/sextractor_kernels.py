#!/usr/bin/env python
import numpy as np
from transfunc import kernelgrid, kernel
import os, sys
from pygoods import Ftable
import cPickle
import KPDFadaptnumpy as KPDF
import scipy
from scipy import linalg, stats

colnames_sex_default = {'mag_in_col': 'mag_in',
   'mag_out_col': 'mag_auto',
   'gtype_col': 'gtype',
   'detect': 'detect'}
simcat_dir = '/Users/khuang/Dropbox/Research/bivariate/galfit_transfer/simcatalogs'

class SExtractorKernelGrid(kernelgrid.KernelGrid1D):
   # a separate class for SExtractor magnitude kernel grid; these kernels will
   # contain no information on galaxy sizes at all. We use these kernels at 
   # the magnitude ranges where S/N is too low for GALFIT to behave well. We 
   # then just take the magnitude information from SExtractor and try to fit
   # the number counts of faint galaxies.
   def __init__(self, kernelFactory, filtername, xlimits=[21.0, 30.0], xname='MAG_AUTO', xunit='mag', dx=0.5):
      super(SExtractorKernelGrid, self).__init__(xlimits, xname, xunit, dx)
      # calculate kernels
      for i in range(self.Nx):
         print "Making kernel %d" % (i)
         xlo, xhi = self.get_kernel_coverage(i)
         print "Magnitude range: [%.2f, %.2f]" % (xlo, xhi)
         k = kernelFactory.run(xlo, xhi)
         self.kernels[i] = k


class SExKernelFactory(object):
   # A class that makes 1D SExtractor magnitude kernels. The structure of this
   # class closely follow that of the GALFIT kernel factory, just without the 
   # size dimension and without checking "good fits".
   def __init__(self, catalog, column_names=colnames_sex_default, disk=0, diskfrac=0.7, expand=0.0, dx=0.02, xmax=2.56, xunit='dmag', xname=r'$\Delta$ mag', Ntot=None):
      assert catalog.endswith('.fits'), "Please provide a FITS table as the catalog."
      self.c = Ftable(catalog)
      self.catalog = catalog
      self.column_names = column_names
      assert disk in [0, 1], "disk should be either 0 or 1."
      self.disk = disk
      if disk == 0:
         self.devauc = 1
      else:
         self.devauc = 0
      self.diskfrac = diskfrac
      self.expand = expand
      # self.Ntot will be the total number of galaxies used from the simulation
      # catalog, after drawing a specified fraction of disks
      if Ntot==None:
         self.Ntot = len(self.c.d)
      else:
         self.Ntot = Ntot
      self.outindex = self.draw_diskfrac()
      # Below are some basic kernel properties
      self.dx = dx  # pixel size for kernels, NOT bin width
      self.xlimits = [-xmax/2., xmax/2.]  # limtis of each kernel, centered around zero
      self.xunit = xunit
      self.xname = xname  # include LaTex syntax for use in plots

   def check_column(self, colname):
      assert hasattr(self.c, colname), "Column %s not in catalog." % colname

   def getcol(self, colname):
      self.check_column(colname)
      return getattr(self.c, colname)

   def draw_diskfrac(self):
      # gtype = [self.disk, self.devauc]
      prob = np.where(self.getcol(self.column_names['gtype_col'])==self.disk, 
                      self.diskfrac, 1.-self.diskfrac)
      indices = np.random.choice(np.arange(len(self.c.d)), size=self.Ntot, 
                                 replace=True, p=prob/prob.sum())
      return indices

   def withinbin(self, maglo, maghi):
      mag_in = self.getcol(self.column_names['mag_in_col'])
      w = (mag_in>=(maglo-self.expand)) & (mag_in<(maghi+self.expand))
      return w.take(self.outindex)

   def run(self, maglo, maghi, bw_method='scott'):
      """
      The guts of calculating kernels using kernel density estimator.
      """
      # expand the cell according to self.expand
      w = self.withinbin(maglo-self.expand, maghi+self.expand)
      recover = self.getcol('detect').take(self.outindex)
      source = np.logical_and(recover, w)
      N_input = np.sum(w)
      N_source = np.sum(source)
      K = kernel.Kernel1D(self.xlimits, self.xname, self.xunit, self.dx)
      K.source = source
      K.maglo = maglo
      K.maghi = maghi
      K.N_input = N_input
      K.N_source = N_source
      if N_source <= 4:
         print "This bin has 4 or fewer detected galaxies... skip."
         K.completeness = 0.
      else:
         K.completeness = float(N_source) / float(N_input)
         print "This bin has completeness %.2f%%" % (K.completeness*100.)
         print "N_source = %d" % N_source
         mag_in_bin = self.getcol(self.column_names['mag_in_col']).take(self.outindex)[source]
         mag_out_bin = self.getcol(self.column_names['mag_out_col']).take(self.outindex)[source]
         rr = mag_out_bin - mag_in_bin
         print "min(rr), max(rr)", rr.min(), rr.max()
         gdata = K.xcenters()
         kpdf = stats.gaussian_kde(rr, bw_method=bw_method)
         pdf1d = kpdf(gdata)
         # replace any NaN by 0.
         if np.isnan(pdf1d).any():
            print "Warning: this kernel has NaN!"
            pdf1d = np.nan_to_num(pdf1d)
         # normalize the kernel to sum to its completeness
         if pdf1d.ravel().sum() > 0:
            # pdf1d = pdf1d * (K.completeness / pdf1d.sum())
            pdf1d = pdf1d / pdf1d.sum()
         K.kernel = pdf1d
         K.bandwidth = kpdf.factor
         K.outindex = self.outindex
      return K

class GDS_SExKernelGrid(SExtractorKernelGrid):
   def __init__(self, catalog, filtername, field, kwargs_factory={}, kwargs_kgrid={}):
      if field == 'udf':
         xlimits = [23., 31.]
      else:
         xlimits = [21., 29.]
      catalog = os.path.join(simcat_dir, catalog)
      kfactory = SExKernelFactory(catalog, **kwargs_factory)
      super(GDS_SExKernelGrid, self).__init__(kfactory, filtername, 
            xlimits=xlimits, **kwargs_kgrid)

