#!/usr/bin/env python

import numpy as np
import distributions

"""
Some simple statistical functions.
"""

def consistent_sigma(val1, sig1, val2, sig2):
   """
   Given the values and erros for both values, test if they are 
   consistent within the errors. The quantities should be 
   val +/- sig.
   """
   assert (sig1 > 0) & (sig2 > 0)
   if val1 >= val2:
      return (val2 + sig2) >= (val1 - sig1)
   else:
      return (val1 + sig1) >= (val2 - sig2)
      
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

