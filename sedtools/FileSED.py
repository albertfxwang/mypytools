#!/usr/bin/env python

"""
Calculate magnitudes given a model SED and filters using pysynphot.
Makes the following scripts obsolete:
- arrayspec.py 
- magredshift.py 
Last updated: 6/15/2014
NOTE: in pysynphot, NEVER use sp.flux for any calculation!!
"""

import numpy as np
import pysynphot as S
from pysynphot import spectrum
import os, copy, sys
import cosmoclass
import igmtrans
import kcorr

defwaveset = S.refs._default_waveset
pytools = os.getenv('mypytools')
vega = S.FileSpectrum(os.path.join(pytools, 'sedtools/alpha_lyr_stis_003.fits'))
cosmo_def = cosmoclass.cosmoclass(H0=70., omega_m=0.3, omega_l=0.7)
pc_cm = 3.0857e18
area_tenpc = 4.*np.pi*(10.*pc_cm)**2  # 4*pi*R**2 with R = 10 parsec in centimeter


class FileSED(spectrum.FileSourceSpectrum):
   ## Made obsolete now with the new-found method in pysynphot?? (2014/11/17)
   def __init__(self, filename, fluxname=None, keepneg=False):
      # just use the definition from FileSourceSpectrum
      super(FileSED, self).__init__(filename, fluxname=fluxname, keepneg=keepneg)
      # by default, flux unit = Flam and wavelength unit = angstrom

   def ABmag(self, bandpass):
      if bandpass == None:
         bandpass = S.ArrayBandpass(defwaveset, np.ones(len(defwaveset)))
      stflux = 10.**(48.6/-2.5)
      # abu = S.ArraySpectrum(defwaveset, 
      #                       np.ones(len(defwaveset))*stflux,
      #                       fluxunits='fnu')
      if bandpass.wave[-1] > self.wave[-1]: 
         # if bandpass wavelenth range is longer than spectrum    
         n = np.searchsorted(bandpass.wave, self.wave[-1])
         merge = np.concatenate((self.wave, bandpass.wave[n:]))
         merge = np.sort(np.unique(merge))
         flux_standard = S.ArraySpectrum(merge, np.ones(len(merge))*stflux,
                                         fluxunits='fnu')
      else:
         flux_standard = S.ArraySpectrum(self.wave, np.ones(len(self.wave))*stflux,
                                         fluxunits='fnu') 
      numerator = (self * bandpass).integrate('fnu')
      denominator = (flux_standard * bandpass).integrate('fnu')
      ratio = numerator / denominator
      if ratio <= 0.:
         abmag = 999.
      else:
         abmag = -2.5 * np.log10(ratio)
      return abmag

   def ABmag_lambda(self, wavelength, w=100.):
      """
      Calculate the AB magnitude at a given rest-frame wavelength lambda.
      Use a boxcar filter centered at wavelength of width w.
      """
      boxcar = S.Box(wavelength, w)
      return self.ABmag(boxcar)


   def Vegamag(self, band):
      vegaflux = (vega * band).integrate()
      # vegaflux corresponds to vegamag 0 in any band
      flux = (self * band).integrate()
      vegamag = -2.5 * np.log10(flux / vegaflux)
      return vegamag

class GalaxySED(FileSED):
   def __init__(self, filename, fluxname=None, keepneg=None, sedproperties={}):
      super(GalaxySED, self).__init__(filename, fluxname=fluxname, keepneg=keepneg)
      self.ebmv = 0.
      # now record all other physical properties of the SED
      self.properties = sedproperties

   def add_dust(self, ebmv, law='xgal'):
      dust = S.Extinction(ebmv, law)
      sp2 = self * dust
      self._fluxtable = self.fluxunits.Convert(self._wavetable, 
                                               sp2.sample(self._wavetable), 
                                               'photlam')

   def add_lya(self, equiv_length):
      # look at simkcorr.py for how to do this
      raise NotImplementedError

class GalaxySEDFactory(object):
   """
   A factory that creates GalaxySED instances corresponding to a given rest-frame
   model SED.
   """
   def __init__(self, filename, cosmo=cosmo_def, sedproperties={}, normmag=-21., normband=1500.):
      self.sp = GalaxySED(filename, sedproperties=sedproperties)
      self.cosmo = cosmo
      # manipulate self.copy, leave alone self.sp
      self.copy = copy.deepcopy(self.sp)  
      if normmag < 99.0:
         # Now, normalize the SED to AB = normmag at rest-frame wavelength normwave
         self.normalize_abmag(normmag, normband)
      self.normmag = normmag
      self.normband = normband
      self.ebmv = 0.
      self.z = 0.

   def reset(self):
      self.copy = copy.deepcopy(self.sp)
      self.normalize_abmag(self.normmag, self.normband)
      self.ebmv = 0.
      self.z = 0.

   def add_dust(self, ebmv, law='xgal'):
      self.copy.add_dust(ebmv, law)
      self.ebmv = ebmv

   def add_lya(self, equiv_length):
      raise NotImplementedError

   def normalize_abmag(self, normmag, normband):
      if type(normband) in [type(0), type(0.0)]:
         # abmag = self.copy.ABmag_lambda(normband)
         filt = S.Box(normband, 100.)
         obs = S.Observation(self.copy, filt)
         abmag = obs.effstim('abmag')
      else:
         # abmag = self.copy.ABmag(normband)
         obs = S.Observation(self.copy, normband)
         abmag = obs.effstim('abmag')
      dmag = normmag - abmag
      self.copy._fluxtable = self.copy._fluxtable * 10.**(-0.4 * dmag)
      self.normband = normband
      self.normmag = normmag
      # self.sp = self.sp * 10.**(-0.4 * dmag) * area_tenpc
      # multiply the flux by 10 pic

   def redshift_no_igm(self, z):
      """
      Including Cosmology in redshifting spectrum
      Note that the input spectrum should have 'flux unit' of L_nu,
      since it will be divided by luminosity distance later.
      """
      if z > 0:
         self.copy._wavetable = self.copy._wavetable * (1. + z)
         # lambda interval is longer at higher z
         if self.cosmo.H0 > 0:
            distance = self.cosmo.lumdist(z, unit='cm') # in cm
            redshifted_flux = self.copy._fluxtable * area_tenpc / (4 * np.pi * distance**2)
            # redshifted_flux = self.copy.flux * area_tenpc / (4 * np.pi * distance**2)
         self.copy._fluxtable = redshifted_flux
         # self.copy._fluxtable = self.copy.fluxunits.Convert(self.copy._wavetable, 
         #                        redshifted_flux, 'photlam')   
      self.z = z
      return self.copy

   def redshift(self, z, igmroutine=igmtrans.meiksin):
      """
      Redshift a model SED, but include the effects of IGM opacity.
      """
      # spec: input spectrum in luminosity/lambda in REST-FRAME
      # any dust extinction should be applied outside this function
      # z: redshift
      # Initialize the IGM routine
      if igmroutine:
         getigm = igmroutine()
      igm = getigm(z)   # apply IGM opacity in the REST-FRAME
      sp_igm = self.copy * igm   # attenuated by IGM **before** redshift
      lam, flam = sp_igm.getArrays()
      # copy the wavelengths and fluxes to self.copy
      self.copy._wavetable = lam
      self.copy._fluxtable = self.copy.fluxunits.Convert(lam, 
                             flam, 'photlam')

      # now redshift the spectrum using the function defined above 
      self.redshift_no_igm(z)
      return self.copy

class ZColorFactory(GalaxySEDFactory):
   """
   Calculates the galaxy color as a function of redshift.
   """
   def __init__(self, filename, filters, cosmo=cosmo_def, sedproperties={},
                normmag=-21., normband=1500.):
      super(ZColorFactory, self).__init__(filename, cosmo=cosmo,
                                 sedproperties=sedproperties, normmag=normmag,
                                 normband=normband)
      self.filters = filters  # a list of pysynphot.SpectralElement instances

   def colors_with_z(self, z0, z1, dz, ebmv=0., extlaw='xgal', igmroutine=igmtrans.meiksin):
      self.num_z = int(round((z1 - z0) / dz))
      self.zarray = np.linspace(z0, z1, num=self.num_z+1)
      self.mags = {}
      for f in self.filters:
         fname = f.name.split(',')[-1]
         self.mags[fname] = np.zeros(len(self.zarray), 'float')
      for i in range(len(self.zarray)):
         z = self.zarray[i]
         self.add_dust(ebmv, law=extlaw)
         sys.stdout.write("z = %.2f\r" % z) 
         sys.stdout.flush()
         sp_z = self.redshift(z, igmroutine=igmroutine)
         for f in self.filters:
            fname = f.name.split(',')[-1]
            self.mags[fname][i] = sp_z.ABmag(f)
         self.reset()

class MConvertFile(GalaxySEDFactory):
   """
   Calculate the conversion from rest-frame 1500 absolute mag to observed-frame
   magnitude as a function of redshift.
   Use the galaxy SED template appropriate for the application, and include
   IGM attenuation effects if redshift is significant (i.e., almost always...)
   Also should include dust attenuation!!!
   """
   def __init__(self, sedfile, restband, cosmo=cosmo_def, normmag=-21.0, ebmv=0.15, extlaw='xgal'):
      # self.factory = GalaxySEDFactory(sedfile, cosmo=cosmo, normmag=normmag, normband=restband)
      self.Q = restband   # rest-frame filter
      self.normmag = normmag
      self.sedfile = sedfile
      self.ebmv = ebmv
      self.extlaw = extlaw
      # self.factory.sp.add_dust(ebmv)
      # self.factory.reset()
      self.cosmo = cosmo

   def __call__(self, z_range, obsband, outputfile):
      # Because m(z) = M + DM(z) + K_QR(z) := M + dm,
      # so dm = DM + K_QR.
      dm = np.zeros(len(z_range))
      sp = S.FileSpectrum(self.sedfile)
      dust = S.Extinction(self.ebmv, self.extlaw)
      sp_ext = sp * dust
      K = kcorr.KCorrect(sp_ext, self.Q, obsband, mode='AB')
      for i in range(len(z_range)):
         z = z_range[i]
         ## Use my own (clumsy) way to compute dm
         # m = self.factory.redshift(z_range[i]).ABmag(obsband)
         # dm[i] = m - self.normmag
         # self.factory.reset()
         ## Use my implementation of Hogg et al. K-correction
         dm[i] = self.cosmo.distmod(z) + K(z)
         if i % 20 == 0:
            sys.stdout.write('z = %.3f  \r' % z_range[i])
            sys.stdout.flush()
      f = open(outputfile, 'wb')
      f.write('# 1 z \n')
      f.write('# 2 dm \n')
      for i in range(len(z_range)):
         f.write('%.4f   %.6f  \n' % (z_range[i], dm[i]))
         
      f.close()
      print "Done!"

   def __str__(self):
      return "SEDFILE = %s;\nREST-BAND = %s" % (self.sedfile, str(self.Q))


