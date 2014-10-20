#!/usr/bin/env python
import numpy as np
from scipy import optimize, stats, special
from robust import MADN

def gauss(x,m,sigma,normed=False):
   sigma2 = sigma*sigma
   y = np.exp(-(x-m)**2/(2.*sigma2))/np.sqrt(2.*np.pi*sigma2)
   if normed:
     return y / np.sum(y)
   else:
     return y

def normgauss(x,m,sigma):
   g = gauss(x,m,sigma)
   return g/np.sum(g)

def fitgauss(x, clip=False, clipsigma=3.0, nmultiple=3, disp=1):
   """
   Fit a Guassian function to a distribution of values x, using maximum
   likelihood.
   """
   Ntot = len(x)
   # First get rid of crazy outliers...
   x_median = np.median(x)
   x_madn = MADN(x)
   x = np.compress(np.abs(x - x_median) <= (100. * x_madn), x)
   xmin = np.min(x); xmax = np.max(x)
   xmin = xmin - 2. * x_madn
   xmax = xmax + 2. * x_madn
   # construct a coordinate array with 3 times more grid points than data points
   xarr = np.linspace(xmin, xmax, Ntot * nmultiple)

   def gauss_likelihood((mu, sig), y):
      logl = -1. * np.log10(gauss(y, mu, sig, normed=False))
      logl = logl.sum()
      return logl
   # clip the sample using sigma if clip==True
   # clip values around the peak of the distribution... use kernel density
   # estimation to estimate the peak
   if clip:
      # runs KDE first to determine the peak of PDF
      # This might fail if crazy outliers exist...
      # So filter out the crazy outliers first
      kernel = stats.gaussian_kde(x, bw_method='scott')
      p = kernel(xarr)
      x_peak = xarr[np.argsort(p)[-1]]
      _endclip = False
      niter = 0
      while not _endclip:
         niter += 1
         x_median = np.median(x)
         x_madn = np.median(np.abs(x - x_median))
         # x = np.compress(np.abs(x - x_median) <= (clipsigma*x_madn), x)
         # now clip around the peak
         if np.sum(np.abs(x-x_peak)>(clipsigma*x_madn)) <= 0.05 * len(x):
            _endclip = True
         x = np.compress(np.abs(x - x_peak) <= (clipsigma * x_madn), x)
         if disp:
            print "iteration %d" % niter
   if disp:
      print "Final sample size in fitting:", len(x)
   # Use downhill simplex algorithm... try other algorithms later?
   out = optimize.fmin(gauss_likelihood, [np.median(x), np.std(x)], args=(x,),
                       disp=disp)
   return out

def gaussCumProb(x):
   """
   Calculates the total cumulative probability for a variable smaller than 
   mean + x_sigma if the variable follows a Gaussian distribution.
   Follows the Error Function page on Wikipedia:
   http://en.wikipedia.org/wiki/Error_function
   """
   cumprob = 0.5 + 0.5 * special.erf(x / np.sqrt(2))
   return cumprob