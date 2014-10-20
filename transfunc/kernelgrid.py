#!/usr/bin/env python

import numpy as np

class KernelGrid1D(object):
   def __init__(self, xlimits, xname, xunit, dx=0.5):
      """
      A super-class of 1D kernel grid. Each bin will have a kernel; the kernel
      could have any dimentionality.
      """
      assert (xlimits[1] > xlimits[0]), "xlimits must be in ascending order."
      self.xlimits = xlimits
      self.xname = xname
      self.xunit = xunit
      self.dx = dx
      self.Nx = int(round((self.xlimits[1] - self.xlimits[0]) / self.dx))
      self.kernels = {}  # keys are integers

   def xedges(self):
      """
      Return the lower-edge value of each bin in x.
      """
      return np.linspace(self.xlimits[0], self.xlimits[1], num=self.Nx+1,
                         endpoint=True)

   def find_bin(self, x, y):
      """
      Given a value for x, return the index of the corresponding kernel.
      """
      assert (x>=self.xlimits[0]) & (x<=self.xlimits[1]), "x is out of bounds."
      ix = np.searchsorted(self.xedges(), x) - 1
      return ix

   def get_kernel(self, ix):
      return self.kernels[ix]

   def find_kernel(self, x):
      """
      Given a value for x, return the corresponding kernel.
      """
      ix = self.find_bin(x)
      assert (ix in self.kernels.keys()), "The kernel %d does not exist." % ix
      print "Index = ", ix
      return self.kernels[ix]

   def get_kernel_coverage(self, ix):
      """
      Given the kernel indices, return the limits in x and y within which this
      kernel applies.
      """
      xlo, xhi = self.xedges()[ix:ix+2]
      return xlo, xhi


class KernelGrid2D(object):
   def __init__(self, xlimits, ylimits, xname, yname, xunit, yunit, 
                dx=0.5, dy=0.2):
      """
      A super-class of 2D kernel grid. Each bin will have a kernel; the kernel
      could have any dimentionality.
      """
      assert (xlimits[1] > xlimits[0]), "xlimits must be in ascending order."
      assert (ylimits[1] > ylimits[0]), "ylimits must be in ascending order."
      self.xlimits = xlimits
      self.ylimits = ylimits
      self.xname = xname
      self.yname = yname
      self.xunit = xunit
      self.yunit = yunit
      self.dx = dx
      self.dy = dy
      self.Nx = int(round((self.xlimits[1] - self.xlimits[0]) / self.dx))
      self.Ny = int(round((self.ylimits[1] - self.ylimits[0]) / self.dy))
      self.kernels = {}
      # self.indices = self.kernels.keys()

   def xedges(self):
      """
      Return the lower-edge value of each bin in x.
      """
      # return np.arange(self.xlimits[0], self.xlimits[1], self.dx)
      return np.linspace(self.xlimits[0], self.xlimits[1], num=self.Nx+1,
                         endpoint=True)

   def yedges(self):
      """
      Return the lower-edge value of each bin in y.
      """
      # return np.arange(self.ylimits[0], self.ylimits[1], self.dy)
      return np.linspace(self.ylimits[0], self.ylimits[1], num=self.Ny+1,
                         endpoint=True)

   def find_bin(self, x, y):
      """
      Given a pair of (x, y), return the index of the corresponding kernel.
      """
      assert (x>=self.xlimits[0]) & (x<=self.xlimits[1]), "x is out of bounds."
      assert (y>=self.ylimits[0]) & (y<=self.ylimits[1]), "y is out of bounds."
      ix = np.searchsorted(self.xedges(), x) - 1
      iy = np.searchsorted(self.yedges(), y) - 1
      return (ix, iy)

   def get_kernel(self, ix, iy):
      return self.kernels[(ix, iy)]

   def find_kernel(self, x, y):
      """
      Given a pair of (x, y), return the corresponding kernel.
      """
      (ix, iy) = self.find_bin(x, y)
      assert ((ix, iy) in self.kernels.keys()), "The kernel (%d, %d) does not exist." % (ix, iy)
      return self.kernels[(ix, iy)]

   def get_kernel_coverage(self, ix, iy):
      """
      Given the kernel indices, return the limits in x and y within which this
      kernel applies.
      """
      xlo, xhi = self.xedges()[ix:ix+2]
      ylo, yhi = self.yedges()[iy:iy+2]
      return xlo, xhi, ylo, yhi

   def find_kernel_coverage(self, x, y):
      """
      Given the (x, y) values, find the coverage of the kernel that also covers
      the input location.
      """
      (ix, iy) = self.find_bin(x, y)
      return self.get_kernel_coverage(ix, iy)

   def transfer_grid_2d(self):
      raise NotImplementedError

