#!/usr/bin/env python

import numpy as np
import scipy
from scipy import optimize

"""
A routine that performs chi2 fitting when there are errors on both coordinates.
"""

def linear_chi2_xy((a,b), x_array, y_array, xerr_array, yerr_array):
   """
   Calculate the chi**2 for a straight line y = a*x + b, when both x and y 
   have errors given by xerr_array and yerr_array.
   Follow Numerical Recipe, Chapter 15.3.
   """
   chi2 = ((y_array - a - b*x_array)**2 / (yerr_array**2 + b**2 * xerr_array**2))
   chi2 = chi2.sum()
   return chi2

def linear_fit_xy(x_array, y_array, xerr_array, yerr_array, guess,
                  full_output=0):
   a0, b0 = guess
   result = optimize.fmin(linear_chi2_xy, guess, 
                          args=(x_array,y_array,xerr_array,yerr_array),
                          full_output=full_output)
   print "Result:"
   print "a=", result[0]
   print "b=", result[1]