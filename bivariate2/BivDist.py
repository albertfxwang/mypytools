#!/usr/bin/env python

import numpy as np
from transfunc import distribution


class BivariateDistribution(distribution.Distribution2D):
   """
   An superclass for bivariate distributions (with variables being physically-
   measureable quantities).
   """
   def shift_x(self, delta_x):
      """
      Shifts the coordinates by delta_x. Pixel size remains the same.
      """
      self.xlimits = self.xlimits + delta_x
      
   def shift_y(self, delta_y):
      """
      Shifts the coordinates by delta_y. Pixel size remains the same.
      """
      self.ylimits = self.ylimits + delta_y

   def resample(self, xgrid, ygrid):
      """
      Resample the value of distribution onto a new pixel grid. The arguments
      xgrid and ygrid need to be increasing arrays of equal steps.
      """
      raise NotImplementedError
      
   def __add__(self, dist2):
      raise NotImplementedError
      