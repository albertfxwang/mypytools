#!/usr/bin/env python

# a script to fit the rest-frame UV slope from a model SED

from numpy import *
import numpy as np
import pysynphot as S
import scipy
from scipy import optimize
import matplotlib.pyplot as plt
import FileSED
from bands import filters
from PhotomTools import photomUtils as pu


class UVSlope(object):
   """
   Defines a spectrum with power-law shape:
   f_lam ~ lam ^ beta
   as is usually used for approximating the starburst rest-frame UV continuum.
   """
   def __init__(self, refwave, beta, redshift, A=1.0, verbose=False):
      rspec = S.PowerLaw(refwave, beta, waveunits='angstrom', 
                              fluxunits='flam')
      self.factory = FileSED.GalaxySEDFactory(spec=rspec)
      # Redshift the spectrum & add IGM attenuation
      self.rspec = self.factory.redshift(redshift)  # the reference spectrum
      mag = S.Observation(self.rspec, filters['f160w']).effstim('abmag')
      # normalize reference spectrum to 25.0 mag in F160W to minimize 
      # round-off errors later in fitting
      dmag = 25.0 - mag
      self.rspec = 10. ** (-0.4 * dmag) * self.rspec
      A = np.abs(A)
      self.rspec = A * self.rspec
      self.beta = beta
      self.A = A
      self.verbose = verbose

   def renorm(self, factor):
      self.spec = self.rspec * factor

   def chi2(self, bands, mags, magerrs, fit='mag'):
      """
      Calculate the total chi-square for self.spec given a list of bands 
      (pysynphot SpectralElement instances) and magnitudes and errors.
      """
      if fit == 'mag':
         mags = np.array(mags)
         magerrs = np.array(magerrs)
         mags_spec = [S.Observation(self.rspec, b).effstim('abmag') for b in bands]
         mags_spec = np.array(mags_spec)
         # chi2 = ((mags_spec - mags)**2 / mags_spec)
         chi2 = (((mags_spec - mags) / magerrs) ** 2).sum()
         if self.verbose:
            print mags_spec, self.beta, self.A, chi2
      elif fit == 'flux':
         fluxes = np.array([pu.ABmag2uJy(m) for m in mags])  # in micro-Jansky
         sn = [pu.magerr2sn(m) for m in magerrs]
         fluxerrs = fluxes / sn  # to get flux errors from S/N
         # Now sample the reference spectrum at the central wavelengths
         self.rspec.convert('muJy')  # convert to milli-Jansky 
         flux_spec = [self.rspec.sample(b.pivot()) for b in bands]
         flux_spec = np.array(flux_spec)  # convert to micro-Jansky
         chi2 = (((flux_spec - fluxes) / fluxerrs) ** 2).sum()
         self.rspec.convert('flam')  # go back to flam
      return chi2


def fit_UVslope(filternames, mags, magerrs, z, beta0=-2.0, normband=filters['f850lp'], normmag=25.0, verbose=False, fit='mag'):
   """
   Fit a UV slope (beta) to the observed rest-frame UV magnitudes.
   bands - a list of strings as filter names (e.g., f160w)
   mags, magerrs - lists of magnitude and magnitude errors
   """
   x = UVSlope(1500., beta0, z)
   mag = S.Observation(x.rspec, normband).effstim('abmag')
   A0 = 10. ** (normmag - mag)
   bands = [filters[b] for b in filternames]
   def func((beta, A)):
      y = UVSlope(1500., beta, z, A, verbose=verbose)
      return y.chi2(bands, mags, magerrs, fit=fit)
   # print beta0, A0, func((beta0, A0))
   output = scipy.optimize.fmin(func, np.array([beta0, A0]), maxiter=1000)
   return output

def show_uvslope_fit(filternames, mags, magerrs, z, beta, A):
   x = UVSlope(1500., beta, z, A=A)
   plt.figure()
   bands = [filters[b] for b in filternames]
   waves = np.array([b.pivot() for b in bands])
   plt.errorbar(waves, mags, yerr=magerrs, ls='none', fmt='x', mew=2, ms=12,
                capsize=10, elinewidth=2)
   x.rspec.convert('fnu')  # in erg cm^-2 s^-1 Hz^-1
   rspec = x.rspec * (1. / 1.e-29)  # in micro-Jansky
   plt.plot(rspec.wave, 23.9-2.5*np.log10(rspec.flux), lw=1.5, 
            label=r'$\beta=%.2f$' % beta)
   plt.xscale('log')
   plt.xlim(0.9*waves.min(), 1.1*waves.max())
   plt.ylim(np.min(mags)-1.0, np.max(mags)+1.0)
   plt.xticks(waves, filternames)
   plt.legend(loc=0)
   plt.gca().invert_yaxis()
   plt.show()


## ================================================
## Below are old (deprecated or recycled) functions
## ================================================

def flam_uv(beta, lam):
   # F_lam ~ lam ^ beta
	# the F_lambda in rest-frame UV between rest-frame lam1 and lam2
	# conventionally the UV-slope beta is determined between 1220 A
	# and 3200 A (Calzetti et al. 1999)
	# caller should make sure that lam is between appropriate wavelengths
	# the normalization is arbitrary -- i.e. flam is NOT normalized
	# F_lam \propto lam**beta
	# beta = -2.0 for a dust-free SED of a typical LBG...?
	return lam ** beta
	
def dflam_uv(beta, wave, sedflux):
	# calculate the difference between SED flux and power-law F_lam
	# first normalize at 2000 A (or the one closest to it)
	flam_powerlaw = flam_uv(beta, wave)
	if 2000. in around(wave):
		flam_powerlaw = flam_powerlaw / flam_powerlaw[wave==2000.] * sedflux[around(wave)==2000.]
	else:
		iclose = argsort(abs(wave-2000.))[0]
		flam_powerlaw = flam_powerlaw / flam_powerlaw[iclose] * sedflux[iclose]
	dflam = flam_powerlaw - sedflux
	return dflam
	
	
def fit_beta(wave, sedflux, beta0=-2.0):
	# use linear least-square method to determine the best-fit beta
	x = optimize.leastsq(dflam_uv, beta0, args=(wave,sedflux), full_output=False)
	return x
	
def fit_beta_sed(sedfile, beta0, ebmv=0.):
	# a wrapper that takes an SED, read it with pysynphot, and find the best-fit beta
	sp = S.FileSpectrum(sedfile)
	if ebmv > 0.:
		sp = sp * S.Extinction(ebmv, 'xgal')
	wave = sp.wave[(sp.wave>=1220.)&(sp.wave<=3200.)]   # the rest-frame wavelength range
	sedflux = sp.sample(wave)   # the flux from sp
	x = fit_beta(wave, sedflux, beta0=beta0)
	return x[0][0]
	
def plot_beta_sed(sp, beta):
	wave = sp.wave[(sp.wave>=1000.)&(sp.wave<=4000.)]
	sedflux = sp.sample(wave)
	iclose = argsort(abs(wave-2000.))[0]
	flam = flam_uv(beta, wave)
	flam = flam / flam[iclose] * sedflux[iclose]
	plt.plot(wave, sedflux, c='blue')
	plt.plot(wave, flam, c='green')
	return 0
	