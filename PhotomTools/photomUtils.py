#!/usr/bin/env python

"""
Some photometry-related utility functions.
"""
import numpy as np

def ABmag2uJy(mag):
   # Convert from AB magnitude to micro-Jansky
   if type(mag)==type(1.0) or type(mag)==type(1):
      if mag >= 99.0:
         return 1e-10
      else:
         return 10.**((23.9 - mag) / 2.5)
   else:
      mag = np.array(mag)
      return np.where(mag >= 99.0, 1.e-10, 10.**((23.9 - mag) / 2.5))   

def uJy2ABmag(flux):
   if type(flux)==type(1.0) or type(flux)==type(1):
      if flux > 0:
         return 23.9 - 2.5 * np.log10(flux)
      else:
         return 99.0
   else:
      flux = np.array(flux)
      return np.where(flux > 0, 23.9 - 2.5 * np.log10(flux), 99.0)
      
def magerr2sn(magerr):
   if type(magerr)==type(1.0) or type(magerr)==type(1):
      if magerr > 0:
         return 1. / (10. ** (0.4 * magerr) - 1.)
      else:
         return -1.0
   else:
      magerr = np.array(magerr)
      return np.where(magerr > 0, 1. / (10. ** (0.4 * magerr) - 1.), -1.0)

def sn2magerr(sn):
   if type(sn)==type(1.0) or type(sn)==type(1):
      if sn > 0:
         return 2.5 * np.log10(1. + 1. / sn)
      else:
         return -1.0
   else:
      sn = np.array(sn)
      return np.where(sn > 0, 2.5 * np.log10(1. + 1. / sn), -1.0)

def calcNsigMag(mag, magerr, N=1.0):
   # Given magnitude & magnitude error, calculate the expected 1-sigma 
   # magnitude. An input magnitude > 90 means it's undetected and it will 
   # just return the magnitude error
   if mag > 90:
      return magerr
   else:
      sn = magerr2sn(magerr)
      return mag + 2.5 * np.log10(sn / N)

def calcColor(mag1, magerr1, mag2, magerr2):
   """
   A simple calculation of the color mag1-mag2, and add their errors in 
   quadrature. If mag1 is undetected (>90), use magerr1 as the 1-sigma upper 
   limit and calculate the lower limit of color. On the other hand, if mag2 is
   undetected, use magerr2 to calculate the upper limit of color.
   """
   if mag1 > 90:
      color = magerr1 - mag2
      print "color >= %.3f" % color
      return color
   elif mag2 > 90:
      color = mag1 - mag2err
      print "color <= %.3f" % color
      return color
   else:
      color = mag1 - mag2
      colorerr = np.sqrt(magerr1**2 + magerr2**2)
      print "color = %.3f +/- %.3f" % (color, colorerr)
      return color, colorerr
      