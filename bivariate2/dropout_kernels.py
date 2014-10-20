#!/usr/bin/env python

import numpy as np
from transfunc import kernelgrid, kernel
import cosmoclass
import os, sys
from scipy.interpolate import interp1d, RectBivariateSpline, UnivariateSpline
from pygoods import Ftable
import cPickle
import copy

cosmo_default = cosmoclass.cosmoclass(70.0, 0.3, 0.7)
filter_dic = {
   'vimos_u':'u',
   'acs_f435w':'b',
   'acs_f606w':'v',
   'acs_f775w':'i',
   'acs_f850lp':'z',
   'wfc3_f125w':'j',
   'wfc3_f160w':'h',
   'wfc3_f105w':'y',
   'wfc3_f098m':'y'}
zeropoints = {
   'wfc3_f160w':25.960,
   'acs_f435w':25.673,
   'acs_f606w':26.505,
   'acs_f775w':25.678,
   'acs_f850lp':24.867,
   'wfc3_f098m':26.270,  
   # I used this as the zero-point in the simulations, although the true zeropoint should be 25.681
   'wfc3_f125w':26.250,
   'vimos_u':26.158,
   'wfc3_f105w':26.270}
zeropoints_goodsv2 = {
   'acs_f435w':25.65288,
   'acs_f606w':26.49341,
   'acs_f775w':25.64053,
   'acs_f850lp':24.84315}

class DropoutKernelGrid(kernelgrid.KernelGrid2D):
   def __init__(self, kernelFactory, xlimits=[-25.0, -15.0], ylimits=[-2.4, 1.0], 
                xname='M1500', yname='logR', xunit='mag', yunit='arcsec', 
                dx=0.5, dy=0.2, run=True):
      kernelgrid.KernelGrid2D.__init__(self, xlimits=xlimits, ylimits=ylimits,
                                       xname=xname, yname=yname, xunit=xunit,
                                       yunit=yunit, dx=dx, dy=dy)
      # All the additional attributes come from kernelFactory
      self.area = kernelFactory.area  # in arcmin^2
      self.area_unit = kernelFactory.area_unit
      self.z0 = kernelFactory.z0  # central (nominal) redshift of the LBG sample
      self.dz = kernelFactory.dz
      print "dz = ", self.dz
      # self.zlimits = [self.z0-1.5, self.z0+1.5]
      self.zlimits = kernelFactory.zlimits
      print "self.zlimits", self.zlimits
      self.Nz = kernelFactory.Nz
      self.cosmo = kernelFactory.cosmo
      self.pixscale = kernelFactory.pixscale
      self.Ndetect = np.zeros((self.Nx, self.Ny), 'int')
      self.Ninput = np.zeros((self.Nx, self.Ny), 'int')
      self.Nselect = np.zeros((self.Nx, self.Ny), 'int')
      self.dVdz = 0.
      self._calc_dVdz = False
      self.interpolated = kernelFactory.interpolate
      self.filename = ""
      # calculate kernels
      if kernelFactory.interpolate:
         self.zlimits_old = copy.deepcopy(self.zlimits)
         self.dz_old = copy.deepcopy(self.dz)
         self.Nz_old = copy.deepcopy(self.Nz)
      if run:
         for i in range(self.Nx):
            for j in range(self.Ny):
               print "Making kernel for (%d, %d)" % (i, j)
               xlo, xhi, ylo, yhi = self.get_kernel_coverage(i, j)
               k = kernelFactory.run(xlo, xhi, ylo, yhi)
               self.kernels[(i, j)] = k
               self.Ndetect[i, j] = k.Ndetect.sum()
               self.Ninput[i, j] = k.Ninput.sum()
               self.Nselect[i, j] = k.Nselect.sum()
               if kernelFactory.interpolate:
                  self.zlimits = k.xlimits
                  self.dz = k.dx
                  self.Nz = k.Nx
               

   def zedges(self):
      """
      Returns the redshift bins over which the kernel function is defined.
      """
      return np.linspace(self.zlimits[0], self.zlimits[1], num=self.Nz+1, 
                         endpoint=True)

   def zcenters(self):
      """
      Return the centers of each redshift bin.
      """
      return self.zedges()[:-1] + self.dz / 2.0

   def get_zindex(self, z):
      """ 
      Return the index of z within self.zedges().
      """
      return np.maximum(0, np.searchsorted(self.zedges(), z) - 1)

   def comoving_volume(self, z):
      return self.cosmo.comoving_volume(z, unit='Mpc3')

   def calc_dVdz(self):
      dVdz = np.zeros(self.Nz)
      comoving_V = np.array(map(self.comoving_volume, self.zedges()))
      self.dVdz = (comoving_V[1:] - comoving_V[:-1]) * self.area / (4.*np.pi) 
      self._calc_dVdz = True     

   def pickle(self, filename):
      # Will this work?
      self.filename = filename
      f = open(filename, 'wb')
      cPickle.dump(self, f, 2)
      f.close()


# Below are "factories" for dropout kernel grids. Write different factories
# for each dropout sample. Within each dropout sample, the factory generates
# kernel grid for each field (UDF, GOODS-S-Deep, GOODS-S-Wide, etc.)
class KernelFactory(object):
   """
   The factory that makes kernel grids. Calculate kernels for each bin and then
   assign to the grid. At the end returns the kernel grid.
   """
   def __init__(self, catalog, z0, zlo, zhi, dz=0.1, interpolate=True, 
                interp_dz=0.02, n_repeat=1, expand=[0., 0.],
                cosmo=cosmo_default, mag_in='m1500_in', re_in='re_in_arcsec',
                zeropoints=zeropoints):
      # More of a template; don't expect to call this constructor directly
      assert catalog.endswith('.fits'), "Please provide a FITS table as the catalog."
      self.c = Ftable(catalog)
      self.mag_in_col = mag_in
      self.re_in_col = re_in
      self.zeropoints = zeropoints
      assert hasattr(self.c, self.mag_in_col), "Column M1500_in not in catalog."
      assert hasattr(self.c, self.re_in_col), "Column %s not in catalog." % self.re_in_col
      self.z0 = z0  # the nominal redshift of the LBG sample
      self.zlimits = [zlo, zhi]
      self.dz = dz
      self.Nz = int(round((zhi - zlo) / dz))
      self.interpolate = interpolate
      self.interp_dz = interp_dz
      self.n_repeat = n_repeat
      self.expand = expand
      self.cosmo = cosmo
      self.bands = []  # will be defined by subclass
      # self.isDropout will be constructed by subclasses

   def getMag(self, flux, fluxerr, magzero, nsigma=1.0):
      """
      Calculate magnitudes from flux and flux error. If S/N < nsigma, use
      flux error to calculate 1-sigma upper limit.
      """
      # First, set mag=99.0 if non detected (fluxerr == -1.0)
      detection = (fluxerr > 0.)
      mag = np.where(detection, magzero - 2.5 * np.log10(flux), 99.0)
      # Next, calculate magnitude limits if S/N < nsigma
      uplim = np.logical_and(flux / fluxerr < nsigma, detection)
      mag = np.where(uplim, magzero - 2.5 * np.log10(fluxerr), mag)
      # whether it is the n-sigma upper limit magnitude or not
      return mag 

   def setMagnitudes(self, band, sn_mode='iso'):
      bs = filter_dic[band]  # short name of filter b
      flux_iso = getattr(self.c, '%s_flux_iso' % bs)
      fluxerr_iso = getattr(self.c, '%s_fluxerr_iso' % bs)
      # calculate 1-sigma upper limits on sources with S/N < 1
      mag_iso = self.getMag(flux_iso, fluxerr_iso, self.zeropoints[band])
      # MAG_ISO --> for color calculation
      # flux_aper = getattr(self.c, '%s_flux_aper' % bs)
      # fluxerr_aper = getattr(self.c, '%s_fluxerr_aper' % bs)
      if sn_mode=='iso':
         setattr(self, band + '_sn', flux_iso / fluxerr_iso)  
      elif sn_mode=='auto':
         flux_auto = getattr(self.c, '%s_flux_auto' % bs)
         fluxerr_auto = getattr(self.c, '%s_fluxerr_auto' % bs)
         setattr(self, band + '_sn', flux_auto / fluxerr_auto)
      # guard against S/N==-inf objects...
      mag_iso = np.where(getattr(self,band+'_sn')==-np.inf, 99.0, mag_iso)
      setattr(self, band + '_mag', np.tile(mag_iso, self.n_repeat))  
      # setattr(self, band + '_sn', flux_aper/fluxerr_aper)
      # S/N calculated using FLUX_ISO

   def SelectDropout(self):
      # Here I implement a method that returns a boxcar P(z) between z0-0.5 
      # and z0+0.5
      zlo = self.z0 - 0.5
      zhi = self.z0 + 0.5
      isDropout = np.where((self.c.redshift>=zlo) & (self.c.redshift<zhi),
                           True, False)
      return isDropout

   def run(self, maglo, maghi, logrlo, logrhi, nfloor=4):       
      """
      Main engine for calculating kernel function.
      """
      self.isDropout = self.SelectDropout()
      assert hasattr(self, 'isDropout'), "Need to perform color selection on the simulated catalog first!"     
      mag_in = getattr(self.c, 'm1500_in')
      logre_in = np.log10(getattr(self.c, self.re_in_col))
      # performs color selection, and stores the result as a boolean array
      # self.isDropout = np.tile(ColorSelection(c), n_repeat)   # needs update
      # Now iterate through all bins
      print "This kernel covers M1500 in [%.1f, %.1f] and logR in [%.1f, %.1f]" % (maglo, maghi, logrlo, logrhi)
      maglo  = maglo - self.expand[0]
      maghi  = maghi + self.expand[0]
      logrlo = logrlo - self.expand[1]
      logrhi = logrhi + self.expand[1]
      print "This kernel uses points selected from M1500 in [%.1f, %.1f] and logR in [%.1f, %.1f]" % (maglo, maghi, logrlo, logrhi)
      bincrit = (mag_in>=maglo)&(mag_in<maghi)&(logre_in>=logrlo)&(logre_in<logrhi)
      print "Total number of points:", np.sum(bincrit)
      if self.n_repeat>1:
         bincrit = np.tile(bincrit, self.n_repeat)
      assert len(bincrit)==len(self.isDropout), "isDropout and bincrit do not have the same length (%d and %d)." % (len(self.isDropout), len(bincrit))
      # Now calculate kernel
      # assert hasattr(self.c, 'z_in'), "Column z_in not in catalog."
      assert hasattr(self.c, 'redshift'), "Column redshift not in catalog."
      assert hasattr(self.c, 'detect'), "Column detect not in catalog."
      K = kernel.Kernel1D(self.zlimits, 'redshift', '', self.dz)
      # z_in = np.tile(self.c.z_in, self.n_repeat)
      z_in = np.tile(self.c.redshift, self.n_repeat)
      detect = np.tile(self.c.detect, self.n_repeat)
      z_array = np.linspace(self.zlimits[0], self.zlimits[1], self.Nz+1)
      # z_array include both endpoints
      selcrit = np.logical_and(bincrit, self.isDropout)
      detcrit = np.logical_and(bincrit, np.tile(self.c.detect, self.n_repeat))
      Pz = np.zeros(self.Nz)
      zarr_bin = z_in[bincrit==True]
      zarr_det = z_in[detcrit]
      zarr_sel = z_in[selcrit]
      K.Ninput = np.histogram(zarr_bin, z_array)[0]  
      # number counts in each redshift interval for all objects in 
      # bins (M,logR)-bin
      K.Ndetect = np.histogram(zarr_det, z_array)[0]
      K.Nselect = np.histogram(zarr_sel, z_array)[0]  
      # number counts in each redshift interval for all objects in 
      # this (M,logR)-bin 
      # selected as dropout
      Pz = K.Nselect.astype('float') / np.maximum(K.Ninput.astype('float'),1.0)
      # Enforce Pz=0 if ninput < nfloor
      Pz = np.where(K.Ninput < nfloor, 0., Pz)
      K.kernel = Pz
      if self.interpolate:
         z_array2 = z_array[:-1] + self.dz / 2.  
         # z_array2 now is the centers of each redshift bin
         fnew = interp1d(z_array2, Pz, kind='cubic')  # cubic interpolation
         z_array2 = np.around(z_array2, 2)
         Nz_interp = int(round((self.zlimits[1] - self.zlimits[0]) / self.interp_dz))
         z_array_interp = np.linspace(self.zlimits[0], self.zlimits[1], 
                                      num=Nz_interp+1)
         z_array_interp = z_array_interp + self.interp_dz / 2.0
         z_array_interp = np.compress((z_array_interp>=z_array2.min()) & (z_array_interp<=z_array2.max()), z_array_interp)
         # Now z_array_interp are the centers of the new redshift bins;
         # the redshift coordinate at which to calculate the interpolated P(z).
         # Also enforce z_array_interp to be within the boundary of z_array2
         Pz = fnew(z_array_interp)
         Pz = np.maximum(Pz, 0.)  # enforce P(z) to be positive or zero
         Pz = np.minimum(Pz, 1.)  # enforce P(z) to be <= 1
         K.xlimits = [z_array_interp[0]-self.interp_dz/2., z_array_interp[-1]+self.interp_dz/2.]
         K.dx = self.interp_dz
         K.Nx = int(round((K.xlimits[1] - K.xlimits[0]) / K.dx))
         K.kernel = Pz
      assert len(K.kernel) == K.Nx, "len(K.kernel)=%d; K.Nx=%d" % (len(K.kernel), K.Nx)
      K.value = K.kernel
      return K



