#!/usr/bin/env python

import numpy as np
import scipy
from scipy import optimize, interpolate
from stats import gauss
from scipy.integrate import cumtrapz
import pysynphot as S

def spectral_line_model(wave, a, b, l_peak, l_center, l_width, wgrid=1):
   """
   Return a model spectrum over the wavelength grid wave. The model has a
   continuum component: f_cont = a * wave + b, and a Gaussian line component:
   f_line = l_peak * gauss(l_center, l_width), where gauss() is supposed to 
   be a Gaussian function with peak 1.
   NOTE: the wavelength grid should be relative narrow, just a bit wider than 
   the line itself. Otherwise the sparser wavelength sampling might causs
   problems for interpolation in the continuum.
   """
   waveRange = wave[-1] - wave[0]
   # Need to convert wave to a regular grid before we can use gauss.gauss
   xgrid = np.arange(wave[0]-0.05*waveRange, wave[-1]+0.05*waveRange, wgrid)
   f_cont = a * xgrid + b
   f_line = gauss.gauss(xgrid, l_center, l_width, xstep=wgrid)  
   # f_cont = a * wave + b
   # f_line = gauss.gauss(wave, l_center, l_width)
   # the peak is not normalized yet
   f_line = f_line * l_peak / f_line.max()  # now it's normalized
   f_total = f_cont + f_line
   f = interpolate.interp1d(xgrid, f_total, kind='linear')
   # Then sample the model at the wavelength grid given
   return f(wave)
   # return f_total

def fit_gauss_line(wave, flux, guess, window=None, rms=-1):
   """
   Given the spectrum (preferrably within a narrow window around the line),
   fit a continuum + Gaussian profile to the spectrum. If the noise of the 
   spectrum is known at each pixel, supply it with the rms argument, and 
   the fit will be to minimize chi-square. Otherwise, use the least-squared
   regression if rms = -1.
   Initial guesses of (a, b, l_peak, l_center, l_width) are given as guess.
   The guess should be pretty close for a good fit!
   One can specify a wavelength window as [wave_lo, wave_high] within which 
   to fit the line.
   Returns the parameters a, b, l_peak, l_center, l_width.
   """
   if window != None:
      wave = wave[(wave>=window[0]) & (wave<window[1])]
      flux = flux[(wave>=window[0]) & (wave<window[1])]
   if rms == -1:
      # No RMS per pixel is specified; use least-squared linear regression
      def cost(params):
         a, b, l_peak, l_center, l_width = params
         model = spectral_line_model(wave, a, b, l_peak, l_center, l_width)
         return np.log10(np.sum((model - flux)**2))
   else:
      raise NotImplementedError
   result = optimize.fmin(cost, np.array(guess))

   return result

def calc_equiv_width(xgrid, flux, line_fit):
   """
   Calculate the equivalent width from the best-fit continuum+line model.
   Define the model over a reguar xgrid; xgrid needs to be monotonically 
   increasing with regular spacing.
   """
   dx = xgrid[1] - xgrid[0]
   f_model = spectral_line_model(xgrid, *line_fit)
   a, b = line_fit[:2]
   f_cont = interpolate.interp1d(xgrid, a * xgrid + b)
   f_line = f_model - f_cont(xgrid)  # line flux above the continuum
   flux_line = cumtrapz(f_line, xgrid)[-1]  # the integrated line flux
   EW = flux_line / f_cont(line_fit[3])  # line_fit[3] is the central wavelength
   return EW

def calc_line_contribution(EW, rest_lam, band, z):
   """
   Calculate the scaling factor between f_nu,obs and f_nu,continuum
   I.e., calculate the true stellar continuum flux given the observed flux 
   and equivalent width. This function return the factor X_EW from Smit+2014:

   f_nu,obs = f_nu,continuum * X_EW

   If EW is a list, then calculate the total contribution from all the lines
   in a given filter. Note that EW is the rest-frame equivalent width.
   rest_lam should be in the same unit as filter.wave, which is usually 
   angstrom.
   """
   if type(EW) == type(1.0) or type(EW) == type(0):
      EW = np.array([float(EW)])
      rest_lam = np.array([float(rest_lam)])
   else:
      EW = np.array(EW)
      rest_lam = np.array(rest_lam)
   obs_lam = (1. + z) * rest_lam
   print "Observed wavelengths (in angstroms):", obs_lam
   X_EW = 0.
   bandInt = cumtrapz(band.throughput/band.wave, band.wave)[-1]
   for i in range(len(EW)):
      x = (EW[i] * (1.+z) * band.sample(obs_lam[i]))
      x = x / (obs_lam[i] * bandInt)
      X_EW += x
   X_EW = 1 + X_EW
   return X_EW



