#!/usr/bin/env python

from numpy import *
import numpy as np
from pygoods import *
import pysynphot as S
from sedtools import magredshift as MR
from pysynphot import Extinction
import myutil
from sedtools import arrayspec as A
from sedtools.bands import filters
from cosmoclass import cosmoclass
import os
from scipy import integrate

## Last updated: 14/07/21

Lsun = 3.826e33  # erg / sec
pc_cm = 3.0857e18  # 1 parsec in cm
uni1500 = S.Box(1500.,100.)
iband = filters['f775w']
zband = filters['f850lp']
tenpc = 4.*pi*(10.*pc_cm)**2  # 10 parsec in centimeter
sedtools = os.getenv('mypytools') + '/sedtools'
cosmo_def = cosmoclass(70., 0.3, 0.7)  # default cosmology

class KCorrect(object):
   """
   Following Hogg et al. 2002, arXiv:0210394, calculcate the K-correction term
   K_QR, which is independent of cosmology.
   Will use pysynphot source spectrum and filter instances, but will *not* use
   its (somewhat opaque) integration methods.
   """
   def __init__(self, source, rest_filter, obs_filter, mode='AB'):
      """
      Calculates the K-correction term of a source spectrum in rest_filter into
      the observed magnitude in obs_filter at redshift z.
      source: source spectrum 
      rest_filter: filter in rest-frame
      obs_filter: filter in the observed frame
      mode: either 'AB' or 'Vega', which determines what the reference spectrum
            is.
      """
      self.f = source
      self.f.convert('flam')
      if mode == 'AB':
         defwaveset = S.refs._default_waveset
         defflux = np.ones(len(defwaveset)) * 10.**(48.6/-2.5)
         self.g = S.ArraySpectrum(defwaveset, defflux, fluxunits='fnu')
         self.g.convert('flam')
      else:
         self.g = S.FileSpectrum(sedtools + '/alpha_lyr_stis_003.fits')
         self.g.convert('flam')
      self.Q = rest_filter 
      self.R = obs_filter

   def __call__(self, z):
      # Use equation 13 in Hogg et al. 2002
      # The equation is K_QR = -2.5*log10((1+z)*(term1*term2)/(term3*term4))
      # And we calculate term1-term4 here.
      lam_e_src = self.f.wave.copy()
      lam_o_src = lam_e_src * (1. + z)
      lam_e_ref = self.g.wave.copy()
      lam_o_ref = lam_e_ref * (1. + z)
      # term1 = integrate(dlam_o*lam_o*flam(lam_o)*R(lam_o))
      y1 = lam_o_src * self.f.sample(lam_o_src/(1+z)) * self.R(lam_o_src)
      term1 = integrate.cumtrapz(y1, lam_o_src)[-1]
      # term2 = integrate(dlam_e*lam_e*g_lam(lam_e)*Q(lam_e))
      y2 = lam_e_ref * self.g.sample(lam_e_ref) * self.Q(lam_e_ref)
      term2 = integrate.cumtrapz(y2, lam_e_ref)[-1]
      # term3 = integrate(dlam_o*lam_o*g_lam(lam_o)*R(lam_0))
      y3 = lam_o_ref * self.g.sample(lam_o_ref) * self.R(lam_o_ref)
      term3 = integrate.cumtrapz(y3, lam_o_ref)[-1]
      # term4 = integrate(dlam_e*lam_e*f_lam((1+z)*lam_e)*Q(lam_e))
      y4 = lam_e_src * self.f.sample(lam_e_src) * self.Q(lam_e_src)
      term4 = integrate.cumtrapz(y4, lam_e_src)[-1]
      # print "term 1-4:", term1, term2, term3, term4
      K_QR = -2.5 * np.log10((1. / (1. + z)) * (term1 * term2) / (term3 * term4))
      return K_QR

   def calcMagnitudes(self, z, cosmo=cosmo_def):
      return cosmo.distmod(z) + self.__call__(z)

   def calcMagConvert(self, mcfile, zlo, zhi, dz, cosmo=cosmo_def):
      # added: 08/15/2014
      # To calculate the magnitude conversion from M1500 to the observed filter
      # that will be read by bivariate2.RLDist.MagConvert
      # dm = m - M = DM(z) + Kcorr(z)
      # Assume that the rest-frame bandpass is unaffected by IGM attenuation!!
      zsteps = round((zhi - zlo) / dz) + 1
      zarr = np.linspace(zlo, zhi, zsteps)
      dm = [cosmo.distmod(z)+self.__call__(z) for z in zarr]
      f = open(mcfile,'w')
      f.write('## Cosmology: H0 = %.2f, omega_m=%.2f, omega_l=%.2f\n' % (cosmo.H0, cosmo.omega_m, cosmo.omega_l))
      f.write('# 1 z\n')
      f.write('# 2 dm\n')
      for i in range(len(dm)):
         f.write('  %10f %15f\n' % (zarr[i], dm[i]))
      f.close()
      return 0 

def kcorr_2obs(sp, absmag_rest, restband, obsband, z, 
          H0=70.0, omega_m=0.3, omega_l=0.7,
          magform='abmag'):
   # defined as m_1 = M_2 + kcorr
   # m_1: apparent magnitude in observed band
   # M_2: absolute magnitude in rest-frame band
   # assumes BC03 model: model flux unit is in L_solar/angstrom
   #sp = S.FileSpectrum(spec)
   # APPLY DUST EXTINCTION ELSEWHERE!
   omega_tot = omega_m + omega_l
   
   absmag = A.ABmag(sp,restband)
   spn = sp * 10.**((absmag_rest-absmag)/-2.5)  # normalize sp to absmag_rest
   spn = spn * tenpc  # get luminosity
   
   # Now redshift spn to z
   obsmag,spec_at,spec_r = MR.mag_redshift(spn, z, obsband, magform, H0=H0,\
      omega_tot=omega_tot, omega_m=omega_m, omega_l=omega_l)
   kc = obsmag - absmag_rest
   return obsmag, kc

def kcorr_2ABmag(sp, obsmag, obsband, restband, z, H0=70., omega_m=0.3,
               omega_l=0.7):
   """
   Given an SED at redshift z, normalize the SED to the observed magnitude 
   obsmag in the filter obsband, and then calculate the rest-frame absolute
   magnitude in the rest-frame filter restband. The input SED should have
   units of luminosity per Angstrom.
   Procedure:
   1. Shift sp to z, calculate the observed magnitude obsmag0 (include IGM 
      extinction).
   2. Calculate the normalization factor f0 that normalizes sp to the observed
      magnitude obsmag: f0 = 10.**(-0.4*(obsmag-obsmag0)).
   3. Multiply sp by f0, and then calculate the absolute magnitude of sp in
      rest frame filter restband (place sp at 10 pc away).
   """
   omega_tot = omega_m + omega_l
   obsmag0, spec_at, spec_r = MR.mag_redshift(sp, z, obsband, 'abmag', 
                              H0=H0, omega_tot=omega_tot, omega_m=omega_m,
                              omega_l=omega_l)
   # Calculate the normalization factor
   f0 = 10. ** (-0.4 * (obsmag - obsmag0))
   spn = sp * f0 * (1./tenpc)
   restmag = A.ABmag(spn, restband)
   return restmag

def mconvertfile(sp, zlo, zhi, dz, mcfile, restband=uni1500, obsband=iband,
   H0=70.0, omega_m=0.3, omega_l=0.7):
   # DEPRECATED: 14/07/21
   zarr = arange(zlo, zhi+dz, dz); print len(zarr)
   dm = zeros(len(zarr))
   for i in range(len(zarr)):
      obsmag, kc = kcorr_2obs(sp, -21.0, restband, obsband, zarr[i], H0=H0,\
         omega_m=omega_m, omega_l=omega_l)
      dm[i] = kc
   f = open(mcfile,'w')
   #f.write('# %10s %15s\n' % ('z', 'dm'))
   f.write('# 1 z\n')
   f.write('# 2 dm\n')
   for i in range(len(dm)):
      f.write('  %10f %15f\n' % (zarr[i], dm[i]))
   f.close()
   return 0 


def make_mconvert_1500_i(sp, zlo, zhi, dz, mcfile, ebmv=0.15, H0=70., 
                         omega_m=0.3, omega_l=0.7):
   calz = S.Extinction(ebmv, 'xgal')  # dust obscuration
   sp_ext = sp10 * calz  # dust obscured SED
   mconvertfile(sp_ext, zlo, zhi, dz, mcfile, restband=uni1500, obsband=iband)
   return 0


def make_mconvert_1500_z(zlo, zhi, dz, mcfile, ebmv=0.15, H0=70., 
                         omega_m=0.3, omega_l=0.7):
   calz = S.Extinction(ebmv, 'xgal')  # dust obscuration
   sp_ext = sp10 * calz  # dust obscured SED
   mconvertfile(sp_ext, zlo, zhi, dz, mcfile, restband=uni1500, obsband=zband,
      H0=H0, omega_m=omega_m, omega_l=omega_l)
   return 0 


