#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pickle_util as pu


class Plot2DKernelGrid(object):
   """
   A base class for making plots of kernel grids.
   """
   def __init__(self, kgrid):
      self.kgrid = kgrid
      k0 = kgrid.kernels[(0,0)]  # a reference kernel
      self.xwidth_kernel = k0.Nx
      self.ywidth_kernel = k0.Ny

   def pixelgrid(self):
      """
      Returns a pixel grid to display all kernels. Individual kernels can be 
      pasted onto the grid.
      """
      k0 = self.kgrid.kernels[(0,0)]
      npix_x = k0.Nx
      npix_y = k0.Ny
      pixgrid = np.zeros((npix_x * self.kgrid.Nx, npix_y * self.kgrid.Ny))
      return pixgrid

   def get_pixelgrid_edge(self, ix, iy):
      """
      Returns the edges of each kernel on the large pixel grid.
      """
      assert (ix, iy) in self.kgrid.kernels, "Kernel (%d,%d) does not exist." % (ix, iy)
      ix0 = ix * self.xwidth_kernel
      ix1 = (ix+1) * self.xwidth_kernel
      iy0 = iy * self.ywidth_kernel
      iy1 = (iy+1) * self.ywidth_kernel
      return ix0, ix1, iy0, iy1

   def get_xcenters(self):
      pixgrid = self.pixelgrid()
      x0 = np.arange(0, pixgrid.shape[0], self.xwidth_kernel)
      return x0 + self.xwidth_kernel / 2.

   def get_ycenters(self):
      pixgrid = self.pixelgrid()
      y0 = np.arange(0, pixgrid.shape[1], self.ywidth_kernel)
      return y0 + self.ywidth_kernel / 2.


class PlotKernel(object):
   """
   A base class for making plots of single kernels.
   """
   def __init__(self, kernel):
      self.kernel = kernel

