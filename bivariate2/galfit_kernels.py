#!/usr/bin/env python

import numpy as np
from transfunc import kernelgrid, kernel
import os, sys
from pygoods import Ftable
import cPickle
import KPDFadaptnumpy as KPDF
import scipy
from scipy import linalg, stats

KPDF_func = {'gaussian': KPDF.MPDFGaussian,
             'epanechnikov': KPDF.MPDFEpanechnikov}
re_err_lim = 0.6  # upper-limit for Re_err / Re returned by GALFIT
colnames_default = {'mag_in_col': 'mag_in',
   'mag_out_col': 'mag_out',
   're_in_col': 're_in_arcsec',
   're_out_col': 're_out_arcsec',
   're_out_err_col': 're_out_err_arcsec',
   'recovered_galfit': 'recovered_galfit',
   'gtype_col': 'gtype',
   'chi2nu_col': 'chi2nu'}
colnames_old = {'mag_in_col': 'mag_in',
   'mag_out_col': 'mag_out',
   're_in_col': 're_in',
   're_out_col': 're_out',
   're_out_err_col': 're_out_err',
   'recovered_galfit': 'recovered_galfit',
   'gtype_col': 'galaxy_type',
   'chi2nu_col': 'chi2nu'}

simcat_dir = '/Users/khuang/Dropbox/Research/bivariate/galfit_transfer/simcatalogs'
kgrid_dir = '/Users/khuang/Dropbox/Research/bivariate/galfit_transfer/galfit_kernels'

class GALFITKernelGrid(kernelgrid.KernelGrid2D):
   def __init__(self, kernelFactory, filtername, xlimits=[21.0, 30.0], ylimits=[-2.4, 1.0], xname='magnitude', yname='logR', pixscale=0.03, xunit='mag', yunit='arcsec', dx=0.5, dy=0.2):
      kernelgrid.KernelGrid2D.__init__(self, xlimits=xlimits, ylimits=ylimits,
                                       xname=xname, yname=yname, xunit=xunit,
                                       yunit=yunit, dx=dx, dy=dy)
      # All the additional attributes come from kernelFactory
      self.pixscale = pixscale
      self.filter = filtername
      self.catalog = kernelFactory.catalog
      self.dxpix = kernelFactory.dx
      self.dypix = kernelFactory.dy
      # calculate kernels
      for i in range(self.Nx):
         for j in range(self.Ny):
            print "-----------------------------"
            print "Making kernel for (%d, %d)..." % (i, j)
            print "Magnitude range: [%.2f, %.2f]" % tuple(self.xedges()[i:i+2])
            print "logR range: [%.2f, %.2f]" % tuple(self.yedges()[j:j+2])
            xlo, xhi, ylo, yhi = self.get_kernel_coverage(i, j)
            k = kernelFactory.run(xlo, xhi, ylo, yhi)
            self.kernels[(i, j)] = k

   def completeness(self):
      c = np.zeros((self.Nx, self.Ny))
      for i in range(self.Nx):
         for j in range(self.Ny):
            c[i,j] = self.kernels[(i,j)].completeness
      return c

   def get_N_input(self):
      N = np.zeros((self.Nx, self.Ny))
      for i in range(self.Nx):
         for j in range(self.Ny):
            c[i, j] = self.kernels[(i,j)].N_input

   def get_N_source(self):
      N = np.zeros((self.Nx, self.Ny))
      for i in range(self.Nx):
         for j in range(self.Ny):
            c[i, j] = self.kernels[(i,j)].N_source

class GFKernelFactory(object):
   """
   The factory that makes 2D GALFIT kernels. Calculate kernels for each bin 
   and then assign to the grid. At the end returns the kernel grid.
   A superclass for all other GALFIT kernel factories.
   """
   def __init__(self, catalog, column_names=colnames_default, disk=0, 
                diskfrac=0.7, re_err_lim=re_err_lim, chi2nu_lim=5.0, 
                shape='gaussian', expand=[0.2, 0.1], dx=0.02, dy=0.02, 
                xmax=2.56, ymax=2.56, xunit='dmag', yunit='dlogR', 
                xname=r'$\Delta$ mag', yname=r'$\Delta\log R$', Ntot=None):
      # More of a template; don't expect to call this constructor directly
      # dx, dy: pixel size of the kernels (in kernel units)
      # xmax, ymax: extent of kernel
      assert catalog.endswith('.fits'), "Please provide a FITS table as the catalog."
      self.c = Ftable(catalog)
      self.catalog = catalog
      self.column_names = column_names
      self.re_err_lim = re_err_lim
      self.chi2nu_lim = chi2nu_lim 
      self.shape = shape
      assert shape in KPDF_func.keys(), "Kernel shape should be either 'gaussian' or 'epanechnikov'."
      self.expand = expand
      print "self.expand = ", self.expand
      assert disk in [0, 1], "disk should be either 0 or 1."
      self.disk = disk
      if disk == 0:
         self.devauc = 1
      else:
         self.devauc = 0
      self.diskfrac = diskfrac
      # self.Ntot will be the total number of galaxies used from the simulation
      # catalog, after drawing a specified fraction of disks
      if Ntot==None:
         self.Ntot = len(self.c.d)
         self.outindex = np.arange(len(self.c.d))
      else:
         self.Ntot = Ntot
         self.outindex = self.draw_diskfrac()
      print "self.outindex = ", self.outindex
      # Below are some basic kernel properties
      self.dx = dx  # pixel size for kernels, NOT bin width
      self.dy = dy
      self.xlimits = [-xmax/2., xmax/2.]  # limtis of each kernel, centered around zero
      self.ylimits = [-ymax/2., ymax/2.]
      self.xunit = xunit
      self.yunit = yunit
      self.xname = xname  # include LaTex syntax for use in plots
      self.yname = yname
      self._check_goodfit = False
      self.goodfit = self.check_goodfit()

   def draw_diskfrac(self):
      # gtype = [self.disk, self.devauc]
      prob = np.where(self.getcol(self.column_names['gtype_col'])==self.disk, 
                      self.diskfrac, 1.-self.diskfrac)
      indices = np.random.choice(np.arange(len(self.c.d)), size=self.Ntot, 
                                 replace=True, p=prob/prob.sum())
      return indices

   def check_column(self, colname):
      assert hasattr(self.c, colname), "Column %s not in catalog." % colname

   def getcol(self, colname):
      self.check_column(colname)
      return getattr(self.c, colname)

   def check_goodfit(self):
      re_out_err = self.getcol(self.column_names['re_out_err_col']).take(self.outindex)
      re_out = self.getcol(self.column_names['re_out_col']).take(self.outindex)
      chi2nu = self.getcol(self.column_names['chi2nu_col']).take(self.outindex)
      recover = self.getcol(self.column_names['recovered_galfit']).take(self.outindex)
      check_re = (re_out_err/re_out<=self.re_err_lim)
      check_chi2nu = (chi2nu<=self.chi2nu_lim)
      check_recovered = (recover==True) & (self.c.magflag==0) & (self.c.reflag==0) & (self.c.nflag==0)
      self._check_goodfit = True
      # return (check_re & check_chi2nu & check_recovered)
      goodfit = reduce(np.logical_and, [check_re, check_chi2nu, check_recovered])
      assert re_out[goodfit].min() > 0, \
         "Some output Re are negative (%.2f)... something wrong here?" % (re_out[goodfit].min())
      return goodfit

   def withinbin(self, maglo, maghi, logrlo, logrhi):
      mag_in = self.getcol(self.column_names['mag_in_col'])
      w = (mag_in>=(maglo-self.expand[0])) & (mag_in<(maghi+self.expand[0]))
      logre_in = np.log10(self.getcol(self.column_names['re_in_col']))
      w = w & (logre_in>=(logrlo-self.expand[1])) & (logre_in<(logrhi+self.expand[1]))
      return w.take(self.outindex)

   def run(self, maglo, maghi, logrlo, logrhi):
      """
      The guts of calculating kernels using kernel density estimator.
      """
      if not self._check_goodfit:
         self.goodfit = self.check_goodfit(re_err_lim, chi2nu_lim)
      # expand the cell according to self.expand... it is already done in 
      # self.withinbin!
      w = self.withinbin(maglo, maghi, logrlo, logrhi)
      recover = np.logical_and(w, 
                self.getcol(self.column_names['recovered_galfit']).take(self.outindex))
      source = np.logical_and(recover, self.goodfit)
      N_input = np.sum(w)
      N_recover = np.sum(recover)
      N_source = np.sum(source)
      K = kernel.Kernel2D(self.xlimits, self.ylimits, self.xname, self.yname,
                          self.xunit, self.yunit, self.dx, self.dy)
      K.source = source
      K.maglo = maglo
      K.maghi = maghi
      K.logrlo = logrlo
      K.logrhi = logrhi
      K.N_input = N_input
      K.N_recover = N_recover
      K.N_source = N_source
      # K.w = w
      last_invcovariance = np.ones((2, 2))
      if N_source <= 4:
         print "This bin has 4 or fewer good GALFIT points... skip."
         K.completeness = 0.
         K.covariance = 0.
         K.lam = 0.
         K.bandwidth = 0.
      else:
         K.completeness = float(N_source) / N_recover
         print "This bin has completeness %.2f%%" % (K.completeness*100.)
         print "N_source = %d" % N_source
         mag_in = self.getcol(self.column_names['mag_in_col']).take(self.outindex)[source]
         mag_out = self.getcol(self.column_names['mag_out_col']).take(self.outindex)[source]
         logre_in = np.log10(self.getcol(self.column_names['re_in_col']).take(self.outindex)[source])
         logre_out = np.log10(self.getcol(self.column_names['re_out_col']).take(self.outindex)[source])
         print "min(mag_in), max(mag_in) = ", np.min(mag_in), np.max(mag_in)
         print "min(logR_in), max(logR_in) = ", np.min(logre_in), np.max(logre_in)
         input_bin = np.array([mag_in, logre_in])
         output_bin = np.array([mag_out, logre_out])
         deviation_bin = output_bin - input_bin
         gdata = KPDF.MPDF2DGrid2Array(K.xcenters(), K.ycenters())
         # print "np.shape(gdata)", gdata.shape
         covariance = np.cov(deviation_bin, rowvar=1)
         rr = deviation_bin.T
         try:
            invcovariance = np.linalg.inv(covariance)
         except: 
            print "*** Warning, singular inverse convariance matrix, trying Moore-Penrose inverse ***" 
            try:
               invcovariance = np.linalg.pinv(covariance, rcond=1.e-10)
            except:
               # This will bomb if it happens on the first pass
               print "*** No joy... using last good invcovariance ***" 
               invcovariance = last_invcovariance 
         last_invcovariance = invcovariance
         renorm_constant = np.sqrt(scipy.linalg.det(covariance))
         bandwidth = KPDF.MPDFOptimumBandwidth(rr)
         lam = np.ones((gdata.shape[0]), "float")
         func = KPDF_func[self.shape]   # function to estimate PDF
         pdf2d = func(rr, rr, bandwidth, lam, invcovariance, renorm_constant)
         lam = KPDF.MPDFAdapt(pdf2d)
         pdf2d = func(rr, gdata, bandwidth, lam, invcovariance, renorm_constant)
         pdf2d = np.reshape(pdf2d, (K.Nx, K.Ny))
         # replace any NaN by 0.
         if np.isnan(pdf2d).any():
            print "Warning: this kernel has NaN!"
            # pdf2d = np.where(np.isnan(pdf2d), 0.)
            pdf2d = np.nan_to_num(pdf2d)
         # normalize the kernel to sum to its completeness
         if pdf2d.ravel().sum() > 0:
            pdf2d = pdf2d * (K.completeness / pdf2d.ravel().sum())
         K.kernel = pdf2d
         K.covariance = covariance
         K.lam = lam
         K.bandwidth = bandwidth
      return K


class GDS_GFKernelGrid(GALFITKernelGrid):
   def __init__(self, catalog, filtername, field, pixscale=0.03, kwargs_factory={}, kwargs_kgrid=dict(ylimits=[-2.4,1.0],yunit='arcsec')):
      if field == 'udf':
         xlimits = [23., 30.]
      else:
         xlimits = [21., 28.]
      catalog = os.path.join(simcat_dir, catalog)
      kfactory = GFKernelFactory(catalog, **kwargs_factory)
      super(GDS_GFKernelGrid, self).__init__(kfactory, filtername, 
            pixscale=pixscale, xlimits=xlimits, **kwargs_kgrid)

class GDS_GFKernelGrid_60mas(GALFITKernelGrid):
   def __init__(self, catalog, filtername, field, pixscale=0.06, kwargs_factory={}, kwargs_kgrid=dict(ylimits=[-1.,2.4],yunit='pixel',pixscale=0.06)):
      if field == 'udf':
         xlimits = [23., 30.]
      else:
         xlimits = [21., 28.]
      catalog = os.path.join(simcat_dir, catalog)
      kfactory = GFKernelFactory(catalog, column_names=colnames_old, **kwargs_factory)
      super(GDS_GFKernelGrid_60mas, self).__init__(kfactory, filtername, xlimits=xlimits, **kwargs_kgrid)
