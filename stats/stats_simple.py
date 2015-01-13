#!/usr/bin/env python

import numpy as np
import distributions

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
