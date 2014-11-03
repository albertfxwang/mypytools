#!/usr/bin/env python

import numpy as np
import scipy
from scipy import stats

class Distribution1D(object):
   """
   Given a 1-D list (array) of values, estimate the PDF with Gaussian KDE and 
   CDF at a grid of points specified by the user.
   """
   def __init__(self, values, bw_method='silverman'):
      """
      Values: a list (array) of data points based on which we estimate the PDF.
      """
      self.values = np.array(values).astype('float')
      self.kernel = stats.gaussian_kde(self.values, bw_method=bw_method)

   def PDF(self, xgrid):
      """
      Use Gaussian KDE to estimate the PDF.
      xgrid: the grid of x values where we estimate the PDF.
      """
      return [xgrid, self.kernel.evaluate(xgrid)]

   def CDF(self, xgrid):
      """
      Use Gaussian KDE to estimate the CDF.
      xgrid: the grid of x values where we estimate the CDF.
      """
      pdf = self.kernel.evaluate(xgrid)
      cdf = np.cumsum(pdf)
      return [xgrid, cdf]

   def minmax(self):
      """
      Return the (min, max) of self.values.
      """
      return self.values.min(), self.values.max()

   def conf_interval(self, xgrid=100, p=0.68, x0=None):
      """
      Estimate the confidence interval for the estimated PDF, given the xgrid
      that the user provides to estimate the CDF.
      p is the total probability enclosed by the interval.
      Optionally, the user can provide the best value x0 around which to 
      construct the confidence interval; the code doesn't really use x0 to 
      compute the confidence interval, but checks if x0 lies within the 
      estimated interval. If not, it will print a warning message and does 
      nothing.
      If xgrid is an integer, it is the number of steps in xgrid. xgrid then
      is between 0.9 * self.values.min() and 1.1 * self.values.max()
      """
      ptail = (1. - p) / 2.  # the probilities in both tails
      pLow = ptail
      pHigh = 1. - ptail
      if type(xgrid) == type(10):
         vmin, vmax = self.minmax()
         xgrid = np.linspace(vmin, vmax, xgrid+1)
      cdf = self.CDF(xgrid)[1]
      cdf = cdf / cdf.max()  # normalize max(CDF) to 1
      xinterval = np.sort(xgrid[(cdf >= pLow) & (cdf <= pHigh)])
      if x0 != None:
         if (x0 < xinterval[0]) or (x0 > xinterval[-1]):
            print "Warning: x0 is NOT within the interval [%.3f, %.3f]!!" % (xinterval[0], xinterval[-1])
         print "Confidence interval: %.4f (+ %.4f) (- %.4f)" % (x0, (xinterval[-1]-x0), (x0-xinterval[0]))
      return xinterval[0], xinterval[-1]
