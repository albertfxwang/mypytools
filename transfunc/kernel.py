#!/usr/bin/env python

import numpy as np
from distribution import *

class Kernel1D(Distribution1D):
   def __init__(self, xlimits, xname, xunit, dx):
      """
      A super-class for 1D kernel functions.
      """
      super(Kernel1D, self).__init__(xlimits, xname, xunit, dx)
      self.kernel = self.value.copy()
      # delattr(self, 'value')

   def __call__(self, x):
      """
      Returns the value at position x.
      """
      assert (x>=self.xlimits[0]) & (x<self.xlimits[1]), "x is out of bounds."
      ix = self.get_index(x)
      return self.kernel[ix]

   def calcKernel(self):
      raise NotImplementedError


class Kernel2D(Distribution2D):
   def __init__(self, xlimits, ylimits, xname, yname, xunit, yunit, dx, dy):
      """
      A super-class for 2D kernel functions.
      """
      super(Kernel2D, self).__init__(xlimits, ylimits, xname, yname, xunit, yunit, dx, dy)
      self.kernel = self.value.copy()
      # delattr(self, 'value')

   def calcKernel(self):
      raise NotImplementedError


