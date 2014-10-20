#!/usr/bin/env python

import numpy as np

# Define a base class for 1D/2D distributions.

class Distribution2D(object):
   """
   A super-class for 2-D distributions.
   """
   def __init__(self, xlimits, ylimits, xname, yname, xunit, yunit, dx, dy):
      self.xlimits = xlimits
      self.ylimits = ylimits
      self.dx = dx
      self.dy = dy
      # Nx, Ny are the number of bins onto which the distribution is sampled.
      self.Nx = int(round((xlimits[1] - xlimits[0]) / dx))
      self.Ny = int(round((ylimits[1] - ylimits[0]) / dy))
      self.xname = xname
      self.yname = yname
      self.xunit = xunit
      self.yunit = yunit
      self.value = np.zeros((self.Nx, self.Ny))

   def xedges(self):
      return np.linspace(self.xlimits[0], self.xlimits[1], num=self.Nx+1)

   def yedges(self):
      return np.linspace(self.ylimits[0], self.ylimits[1], num=self.Ny+1)

   def xcenters(self):
      return self.xedges()[:-1] + self.dx / 2.0

   def ycenters(self):
      return self.yedges()[:-1] + self.dy / 2.0

   def shape(self):
      return (self.Nx, self.Ny)

   def get_index(self, x, y):
      """
      Return the pixel index given the position (x, y).
      """
      ix = np.searchsorted(self.xedges(), x) - 1 
      iy = np.searchsorted(self.yedges(), y) - 1
      ix = np.maximum(ix, 0)
      iy = np.maximum(iy, 0)
      return (ix, iy)

   def get_xindex(self, x):
      ix = np.searchsorted(self.xedges(), x) - 1 
      ix = np.maximum(ix, 0)
      ix = np.minimum(ix, self.Nx-1)
      return ix

   def get_yindex(self, y):
      iy = np.searchsorted(self.yedges(), y) - 1
      iy = np.maximum(iy, 0)
      iy = np.minimum(iy, self.Ny-1)
      return iy

   def get_xcoord(self, x):
      """
      For a given real value of x, calculate the x-coordinate value in the 
      pixel coordinate.
      """
      return (x - self.xlimits[0]) / self.dx

   def get_ycoord(self, y):
      """
      For a given real value of y, calculate the y-coordinate value in the 
      pixel coordinate.
      """
      return (y - self.ylimits[0]) / self.dy

   def normalize(self, normvalue=1.0):
      """
      Normalize the integral (discrete sum) of distribution to normvalue.
      """
      assert self.value.sum() > 0, "Distribution is zero."
      total = self.value.sum() * self.dx * self.dy
      self.value = self.value * (normvalue / total)

   def within_bounds_x(self, x):
      """
      Run an assertion statement to make sure that x is within xlimits.
      """
      assert (x>=self.xlimits[0]) & (x<=self.xlimits[1]), "x is out of bounds."

   def within_bounds_y(self, y):
      assert (y>=self.ylimits[0]) & (y<=self.ylimits[1]), "y is out of bounds."

   def window_function(self, xmin, xmax, ymin, ymax):
      """
      Returns a 2D array representing self.value that is zeroed-out everywhere 
      except within a certain window.
      """
      ixmin, iymin = self.get_index(xmin, ymin)
      ixmax, iymax = self.get_index(xmax, ymax)
      new_distribution = np.zeros(self.value.shape)
      new_distribution[ixmin:ixmax+1,iymin:iymax+1] = self.value[ixmin:ixmax+1,iymin:iymax+1]
      return new_distribution

   def update_value(self, value):
      self.value = value

   def __call__(self, x, y, value=None):
      """
      Returns the probability at position x, y.
      """
      self.within_bounds_x(x)
      self.within_bounds_y(y)
      ix, iy = self.get_index(x, y)
      if value == None:
         value = self.value
      assert np.shape(value) == self.value.shape
      val = value[ix, iy] * self.dx * self.dy
      if val < 0.:
         print "Warning: probability is negative! x=%.2f, y=%.2f, P=%.2e" % (x, y, val)
      if val == 0:
         print "Warning: zero probability... x=%.2f, y=%.2f" % (x, y)
      return val

class Distribution1D(object):
   def __init__(self, xlimits, xname, xunit, dx):
      """
      A super-class for 1D distributions.
      """
      self.xlimits = xlimits
      self.xname = xname
      self.xunit = xunit
      self.dx = dx
      self.Nx = int(round((self.xlimits[1] - self.xlimits[0]) / self.dx))
      self.value = np.zeros(self.Nx)

   def xedges(self):
      return np.linspace(self.xlimits[0], self.xlimits[1], num=self.Nx+1)

   def xcenters(self):
      return self.xedges()[:-1] + self.dx / 2.0

   def get_index(self, x):
      self.within_bounds_x(x)
      return np.maximum(0, np.searchsorted(self.xedges(), x) - 1)

   def within_bounds_x(self, x):
      """
      Run an assertion statement to make sure that x is within xlimits.
      """
      assert (x>=self.xlimits[0]) & (x<=self.xlimits[1]), "x is out of bounds."

   def __call__(self, x):
      """
      Returns the value at position x.
      """
      assert (x>=self.xlimits[0]) & (x<self.xlimits[1]), "x is out of bounds."
      ix = self.get_index(x)
      return self.value[ix]
