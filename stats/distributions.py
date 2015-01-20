#!/usr/bin/env python

import numpy as np
import scipy
from scipy import stats

class Distribution1D(stats.rv_discrete):
   """
   Given a 1-D list (array) of values, estimate the PDF with Gaussian KDE and 
   CDF at a grid of points specified by the user.
   """
   def __init__(self, values, weights=-1, bw_method='silverman'):
      """
      Values: a list (array) of data points based on which we estimate the PDF.
      """
      # self.values = np.array(values).astype('float')
      if weights < 0:
         weights = np.ones(len(values)) / len(values) # equal weights
      assert len(weights) == len(values)
      weights = np.array(weights) / np.sum(weights)
      stats.rv_discrete.__init__(self, a=0, b=1, values=(values, weights))
      self.kernel = stats.gaussian_kde(values, bw_method=bw_method)

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
      return self.cdf(xgrid)

   def minmax(self):
      """
      Return the (min, max) of self.values.
      """
      return self.xk[0], self.xk[-1]

   def conf_interval(self, xgrid=None, p=0.68, x0=None, print_it=False, scale=1):
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
      if xgrid == None:
         xgrid = self.xk
      elif type(xgrid) == type(10):
         vmin, vmax = self.minmax()
         xgrid = np.linspace(vmin, vmax, xgrid+1)
      cdf = self.CDF(xgrid)
      # print xgrid
      # cdf = cdf / cdf.max()  # normalize max(CDF) to 1
      if (cdf[0] > pHigh) or (cdf[1] < pLow):
         # The distribution is strongly non-Gaussian, with most of the 
         # probabilities piled up at the low end... will quote the maximum
         # value instead
         print "*********************************"
         print "Distribution is strongly skewed!!"
         print "*********************************"
         if print_it: 
            print "Confidence interval: %.4f (+ %.4f) (- %.4f)" % (scale*x0, scale*(xgrid[-1]-x0), scale*(x0-xgrid[0]))
         return xgrid[0], xgrid[-1]
      xinterval = np.sort(xgrid[(cdf >= pLow) & (cdf <= pHigh)])
      xinterval = scale * xinterval
      # print len(xinterval)
      if x0 != None:
         x0 = scale * x0
         if (x0 < xinterval[0]) or (x0 > xinterval[-1]):
            print "Warning: x0 (%.3f) is NOT within the default interval [%.3f, %.3f]!!" % (x0, xinterval[0], xinterval[-1])
            if (x0 < xinterval[0]):
               print "Set CDF(x0) to pLow..."
               xlow = np.max(xgrid[xgrid<=(x0/scale)])
               pLow = self.CDF(xlow)
               pHigh = pLow + p
            else:
               print "Set CDF(x0) to pHigh..."
               xhigh = np.min(xgrid[xgrid>=(x0/scale)])
               pHigh = self.CDF(xhigh)
               pLow = pHigh - p
            print "New values for pLow, pHigh:", pLow, pHigh
            xinterval = np.sort(xgrid[(cdf >= pLow) & (cdf <= pHigh)]) * scale
         if print_it: print "Confidence interval: %.4f (+ %.4f) (- %.4f)" % (x0, (xinterval[-1]-x0), (x0-xinterval[0]))
      return xinterval[0], xinterval[-1]
