#!/usr/bin/env python

import numpy as np
import distributions
from scipy import stats
from scipy.stats import mstats

"""
Some simple statistical functions.
"""

def consistent_sigma(values, sigmas):
   """
   Given the values and erros for both values, test if they are 
   consistent within the errors. The quantities should be 
   val +/- sig.
   """
   assert np.min(sigmas) > 0
   imin = np.argsort(values)[0]
   minval = values[imin]
   minsig = sigmas[imin]
   imax = np.argsort(values)[-1]
   maxval = values[imax]
   maxsig = sigmas[imax]
   return (minval+minsig >= maxval-maxsig)
      
def MonteCarlo_dist(values, error):
   """
   Given an array of values, perturb each value according to the error
   assuming Gaussian distribution. Error could be a single number that is the 
   same for all values, or an array (with the same length as values) that 
   is the error for each member of values.
   """
   if type(error) in [type(1.0), type(1)]:
      error = np.ones(len(values), 'float') * error
   newValues = np.random.normal(values, error)
   return newValues

def confidence_interval(data, bestValue, p=0.68, scale=1, verbose=True):
   """
   Given an array of data and a best value, calculate the empirical confidence
   interval around bestValue that bracket a probability p.
   """
   p0 = (1. - p) / 2.
   p1 = 1. - p0
   x0, x1 = mstats.mquantiles(data, prob=[p0, p1])
   if not (bestValue > x0) and (bestValue < x1):
      print "Warning: bestValue is not within the confidence interval! Do something about it."
   interval = np.array([bestValue - x0, x1 - bestValue]) * scale
   bestValue = bestValue * scale
   if verbose:
      print "Confidence interval: %.4f (+ %.4f) (- %.4f)" % (bestValue, interval[1], interval[0])
   return interval

def MonteCarlo_average(distributions, bestValue, nsamp=1000, **cf_kwargs):
   """
   Get the confidence interval of the average from a number of MC distributions.
   distributions should be a list with nobj distributions.
   """
   nobj = len(distributions)
   output = np.zeros(nsamp)
   for i in range(nsamp):
      x = [np.random.choice(distributions[j]) for j in range(nobj)]
      output[i] = np.average(x)
   # Now get the confidence intervals
   cf = confidence_interval(output, bestValue, **cf_kwargs)
   return cf, output
