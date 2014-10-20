#!/usr/bin/env python

# a script to fit the rest-frame UV slope from a model SED

from numpy import *
import pysynphot as S
import scipy
from scipy import optimize
import matplotlib.pyplot as plt

def flam_uv(beta, lam):
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
	