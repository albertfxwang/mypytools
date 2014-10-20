#!/usr/bin/env python

import numpy as np
from transfunc import distribution
from BivDist import BivariateDistribution
from pygoods import sextractor
from dropout_kernels import DropoutKernelGrid
from hconvolve import hconvolve
from multiprocessing import Manager, Queue, Process
import copy
# import multiprocessing

# Defines the bivariate size-luminosity distribution class & factory

# mc_f606w = MagConvert('/Users/khuang/Dropbox/Research/bivariate/bivariate_fit/kcorr/M1500_to_f606w_omega_m_0.3.txt')
# print isinstance(mc_f606w, MagConvert)
# print type(mc_f606w)

class MagConvert(object):
   # defines magnitude conversion class (from M_1500 to m)
   def __init__(self, mconvfile):
      c = sextractor(mconvfile)
      self.zarr = c.z
      self.dmarr = c.dm
      self.__doc__ = """
   The first column of mconvfile needs to be redshift, and second column needs to be 
   dm (where M_1500 + dm = m_app).
   """
   def __call__(self, z):
      # does interpolation of magnitude conversion factor dm
      i = np.searchsorted(self.zarr, z)
      if i == 0:
         dz = z - self.zarr[i]  # dz <= 0
         ddm = (self.dmarr[1] - self.dmarr[0]) / (self.zarr[1] - self.zarr[0])*dz
         dm = self.dmarr[0] - ddm
      elif i == len(self.zarr):
         dz = z - self.zarr[i-1]
         ddm = (self.dmarr[-1] - self.dmarr[-2]) / (self.zarr[-1] - self.zarr[-2])*dz
         dm = self.dmarr[-1] + ddm
      else:    
         dz = z - self.zarr[i-1]
         ddm = (self.dmarr[i] - self.dmarr[i-1]) / (self.zarr[i] - self.zarr[i-1])*dz
         dm = self.dmarr[i-1] + ddm
      return dm

   def revert2z(self, dm):
      """
      Given dm, find the corresponding redshift z. Only works when dm is also 
      a monotonically increasing function of z (fortunately it is).
      """
      i = np.searchsorted(self.dmarr, dm)
      if i == 0:
         ddm = dm - self.dmarr[0]
         dz = (self.zarr[1] - self.zarr[0]) / (self.dmarr[1] - self.dmarr[0]) * ddm
         z = self.zarr[0] - dz
      elif i == len(self.dmarr):
         ddm = dm - self.dmarr[-1]
         dz = (self.zarr[-1] - self.zarr[-2]) / (self.dmarr[-1] - self.dmarr[-2]) * ddm
         z = self.zarr[-1] + dz
      else:
         ddm = dm - self.dmarr[i]
         dz = (self.zarr[i] - self.zarr[i-1]) / (self.dmarr[i] - self.dmarr[i-1]) * ddm
         z = self.zarr[i] + dz
      return z

class RLDistribution(BivariateDistribution):
   def __init__(self, maglimits, logrlimits, dmag, dlogr, pixscale=0.03,
                M0=-21.0, MDR=7.0):
      super(RLDistribution, self).__init__(maglimits, logrlimits, 'magnitude',
                                           'logR', 'mag', 'arcsec', dmag, dlogr)
      self.pixscale = pixscale
      self.alpha = 0.
      self.mstar = 0.
      self.logr0 = 0.
      self.sigma = 0.
      self.beta = 0.
      self.phistar = 0.
      self.z = 0.
      self.M0 = M0
      self.MDR = MDR

   def compute(self, alpha, mstar, logr0, sigma, beta, magshift=0.):
      """
      Calculate the parameterized size-luminosity distribution.
      """
      self.value_last = self.value.copy()
      self.alpha = alpha
      self.mstar = mstar
      self.mstar2 = mstar + magshift
      M0 = self.M0 + magshift
      if np.logical_or((self.mstar2-self.xlimits[1])>self.MDR, (self.xlimits[0]-self.mstar2)>self.MDR):
         print "Warning: mstar is too far away from the distribution limits; mstar2=%.2f, mag limits=%s" % (self.mstar2, self.xlimits)
      self.logr0 = logr0
      self.sigma = sigma
      self.beta = beta
      npix = [self.Nx, self.Ny]
      M_array = self.xedges()[:-1]
      logr_array = self.yedges()[:-1]
      # Compute the distribution functions at the grid points and insert 
      # into array
      result = np.zeros(npix, "float")
      # directly computing Schechter function without importing schechter.py
      delta_M = M_array - self.mstar2
      lf = 0.4 * np.log(10.) * 10.**(-0.4*(self.alpha + 1.) * delta_M) * np.exp(-10.**(-0.4 * delta_M))
      # this is the part of LF without the normalization constant
      for i in range(len(M_array)):
         expo = -1. * (np.log(10.))**2
         expo = expo * ((logr_array - self.logr0 + 0.4 * self.beta * (M_array[i] - M0)))**2
         expo = expo / (2. * self.sigma**2)
         sd = np.log(10.) / (np.sqrt(2. * np.pi) * self.sigma) *np.exp(expo)
         result[i] = lf[i]*sd
      self.value = result
      # enforce dynamical range in M so that LF does not extend to infinity
      M_uplim = self.mstar2 + self.MDR
      M_lolim = self.mstar2 - self.MDR
      self.value = self.window_function(M_lolim, M_uplim, self.ylimits[0], self.ylimits[1])
      if M_lolim >= self.xlimits[1]:
         self.value = np.zeros(npix)
      if M_uplim <= self.xlimits[0]:
         self.value = np.zeros(npix)
      return self.value

   def overlap(self, maglo, maghi, logrlo, logrhi):
      """
      Determine whether the window defined by [maglo:maghi,logrlo:loghi] 
      overlaps with the limits of the distribution.
      """
      x_overlap = (maglo < self.xlimits[1]) | (maghi > self.xlimits[0])
      y_overlap = (logrlo < self.ylimits[1]) | (logrhi > self.ylimits[0])
      return (x_overlap & y_overlap)

   def sum(self):
      """
      Calculate the sum of the probability density.
      Returns np.sum(self.value)*self.dx*self.dy
      """
      return self.value.sum() * self.dx * self.dy

      
class RLDistributionFactory(object):
   # processes the size-luminosity distribution
   # subclasses will define methods with more specific elements for each 
   # dropout sample
   # Default pixel size is 0.02 in both dimensions... should I change it later?
   def __init__(self, maglimits, logrlimits, mag_convert, pixdx=0.02, pixscale=0.03):
      self.maglimits = maglimits
      self.logrlimits = logrlimits
      self.pixscale = pixscale
      self.RLDist = RLDistribution(self.maglimits, self.logrlimits, pixdx, 
                                   pixdx, pixscale=pixscale)
      self.Pijk_dVdz = 0.
      self._calc_Pijk_dVdz = False
      self.mag_convert = mag_convert

   def reset_RLDist(self):
      self.RLDist.value = np.zeros(self.RLDist.value.shape)

   def calc_Pijk_dVdz(self, DK):
      """
      Calculate P(M, R, z)*dV/dz for all redshift bins using dropout kernel 
      grid DK.
      Returns a 3D array Pk, which contains a 2D model in the (m, logR) space for 
      every redshift bin zk. At each redshift zk, the 2D model has the value 
      P_ijk(z)*dV/dz; since P_ijk is divided into many (m, logR) cells, 
      Pk[zk] is a 2D distribution that has constant values in each (m,logR) cell 
      (with the values being P_ijk(z)). Pk[zk] is multiplied (pixel by pixel) by 
      the bivariate distribution phi(M,logR) later on and integrate over M to 
      get the corrected bivariate distribution phi'(m,logR). The conversion 
      between absolute mag M and apparent mag m is done when "pasting" P_ijk(z)
      onto the appropriate (m,logR) cell.
      """
      print "Calculate P(m, logR, z)*dV/dz..."
      z_array = DK.zcenters() # the central redshift of each bin
      Pijk = np.zeros((len(z_array), self.RLDist.Nx, self.RLDist.Ny))  
      if not DK._calc_dVdz:
         DK.calc_dVdz()
      self.dVdz = DK.dVdz
      self.z_array = z_array
      # a 3D array, with 2D arrays the same shape as the (m, logR) model in 
      # each redshift slice
      # calculate Pk in the (m, logR) plane
      for k in range(DK.Nz):
         zk = z_array[k]
         magshift = self.mag_convert(zk)
         for key in DK.kernels.keys():
            kernel = DK.kernels[key]
            maglo, maghi = np.array(DK.get_kernel_coverage(*key)[:2]) + magshift
            logrlo, logrhi = DK.get_kernel_coverage(*key)[2:]
            if self.RLDist.overlap(maglo, maghi, logrlo, logrhi):
               ix0, iy0 = self.RLDist.get_index(maglo, logrlo)
               ix1, iy1 = self.RLDist.get_index(maghi, logrhi)
               Pijk[k, ix0:ix1, iy0:iy1] = kernel(zk)
         Pijk[k] = Pijk[k] * DK.dVdz[k]
      self.Pijk_dVdz = Pijk
      self._calc_Pijk_dVdz = True

   def RLdist_boxcar(self, parameters, z0, dz=1.0):
      """
      Calculate the size-luminosity distribution if P(i,j,k) is a boxcar within
      [z0-dz/2, z0+dz/2].
      """
      zlo = z0 - dz / 2.
      zhi = z0 + dz / 2.
      print zlo, zhi
      # z_array = self.z_array[(self.z_array>=zlo)&(self.z_array<zhi)]
      # dVdz = self.dVdz[(self.z_array>=zlo)&(self.z_array<zhi)]
      dist_array = np.zeros(self.RLDist.value.shape)
      for i in range(len(self.z_array)):
         zi = self.z_array[i]
         if (zi >= zlo) & (zi < zhi):
            dist = self.RLdist_from_z(parameters, zi, None)
            dist_array = dist_array + dist.value * self.dVdz[i]
      return dist_array

   def RLdist_from_z(self, parameters, z, DK):
      # calculate the (m, logr) distribution slice at a redshift given free 
      # parameters and a class that calculates magnitude shift from redshift.
      # Also use dropout-selection kernels P(M, logR, z) to multiply by 
      # the effective volume.
      magshift = self.mag_convert(z)
      if DK == None:
         dist = self.RLDist.compute(*parameters, magshift=magshift)
         return self.RLDist
      # assert isinstance(DK, DropoutKernelGrid)
      if not self._calc_Pijk_dVdz:
         self.calc_Pijk_dVdz(DK)
      k = DK.get_zindex(z)
      distribution_z = self.RLDist.compute(*parameters, magshift=magshift) * self.Pijk_dVdz[k]
      return distribution_z

   def RLdist_dropout_selection(self, parameters, DK):
      # calculate each slice of (m, logR) distribution from z, and then 
      # sum them up
      ## Consider using multiprocessing to speed up...
      distribution_list = map(lambda x: self.RLdist_from_z(parameters, x, DK), 
                              DK.zcenters())
      self.RLDist.value_last = self.RLDist.value.copy()
      self.RLDist.value = reduce(lambda x, y: x+y, distribution_list)
      return self.RLDist

   def apply_transfer_2d(self, kgrid, maglim=99.0, nproc_model=4):
      """
      Convolve with GALFIT transfer function kernels in each bin. Can also 
      just apply kernels brighter than a certain mag limit (maglim).
      """
      # print "gfmaglim: ", maglim
      KG = kgrid
      self.RLDist.value_last = self.RLDist.value.copy()
      assert abs(self.RLDist.dx - KG.dxpix)<=1.e-4, "Pixel dimensions in x between the distribution and the kernels do not agree."
      assert abs(self.RLDist.dy - KG.dypix)<=1.e-4, "Pixel dimensions in y between the distribution and the kernels do not agree."
      convolved = np.zeros(self.RLDist.value.shape, "float")
      # Figure out which section in model to be convolved with which kernel
      chunksize = len(KG.kernels.keys()) / nproc_model
      k_indices = np.linspace(0, len(KG.kernels.keys()), num=nproc_model+1)
      k_indices = np.around(k_indices).astype('int')
      manager = Manager()
      q_out = manager.Queue()
      self.processes = []
      def worker(keys, q):
         for key in keys:
            k = KG.kernels[key[0],key[1]]
            if k.completeness > 0:
               maglo, maghi = np.array(KG.get_kernel_coverage(*key)[:2])
               logrlo, logrhi = KG.get_kernel_coverage(*key)[2:]
               # if maglo >= maglim:
               if not self.RLDist.overlap(maglo, maghi, logrlo, logrhi):
                  # this kernel is too faint... continue
                  continue
               # elif self.RLDist.overlap(maglo, maghi, logrlo, logrhi):
               else:
                  section = self.RLDist.window_function(maglo, maghi, logrlo, logrhi)
                  # If we already parallel the kernels, then DO NOT use
                  # parallelization in hconvolve again! Python doesn't like that.
                  section = hconvolve(section, k.kernel, threads=1)
                  q.put(section)
      for i in range(nproc_model):
         i0 = k_indices[i]
         i1 = k_indices[i+1]
         klist = KG.kernels.keys()[i0:i1]
         p = Process(target=worker, args=(klist, q_out))
         self.processes += [p]
         p.start()
      # join the processes
      for i in range(nproc_model):
         self.processes[i].join()
      # collect output
      while not q_out.empty():
         convolved = convolved + q_out.get()
      self.RLDist.value = convolved
      return self.RLDist 

   def apply_transfer_1d(self, kgrid, maglim=26.5):
      """
      Apply a 1D SExtractor magnitude kernels at the FAINT END. This only 
      happens for magnitudes FAINTER than maglim.
      """
      convolved = np.zeros(self.RLDist.value.shape, "float")
      # temporary store the GALFIT kernel-corrected distribution
      value_temp = self.RLDist.value.copy()
      # self.RLDist.value_last = value_temp
      # bring back the distribution before GALFIT kernel corrections
      self.RLDist.value = self.RLDist.value_last.copy()
      for key in kgrid.kernels.keys():
         k = kgrid.kernels[key]
         if k.completeness > 0:
            maglo, maghi = kgrid.get_kernel_coverage(key)
            if (maglo <= maglim) & (maghi > maglim):
               # apply 1-D kernel
               section = self.RLDist.window_function(maglo, maghi, 
                        self.RLDist.ylimits[0], self.RLDist.ylimits[1])
               k2d = np.zeros((len(k.kernel), 5))
               k2d[:,2] = k.kernel
               section = hconvolve(section, k2d, threads=1)
               convolved = convolved + section
      # now combine this with self.RLDist.value
      self.RLDist.value = value_temp + convolved
      return self.RLDist

   def calculate(self, parameters, DK=None, GK=None, SK=None, gfmaglim=26.5, semaglim=26.5, nproc_model=4, return_copy=False, verbose=0):
      # Take everything into account
      self.reset_RLDist()
      if DK == None:
         m = self.RLdist_from_z(parameters, self.z0, None)
      else:
         if verbose: print "apply dropout selection..."
         m = self.RLdist_dropout_selection(parameters, DK)
      if GK != None:
         if verbose: print "apply_transfer_2d..."
         m = self.apply_transfer_2d(GK, nproc_model=nproc_model, maglim=gfmaglim)
      if SK != None:
         if verbose: print "apply_transfer_1d..."
         m = self.apply_transfer_1d(SK, maglim=semaglim)
      if return_copy:
         return copy.deepcopy(m)
      else:
         return m

