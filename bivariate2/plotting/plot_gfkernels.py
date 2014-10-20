#!/usr/bin/env python

## IN BIVARIATE2

__version__ = "1.0"

import numpy as np
import matplotlib.pyplot as plt
import plot_kernels as pk
import os


class PlotGFKernelGrid(pk.Plot2DKernelGrid):
   def __init__(self, kgrid):
      super(PlotGFKernelGrid, self).__init__(kgrid)
      
   def show_all_kernels(self, normalize=True, vmax_factor=1., ms=30):
      """
      Show all GALFIT kernels on a large pixel grid.
      vmax_factor: set vmax to np.max(pixgrid) / vmax_factor for the entire 
                   pixel grid in plt.imshow().
      """
      pixgrid = self.pixelgrid()
      Nx, Ny = pixgrid.shape
      xlimits = self.kgrid.xlimits
      ylimits = self.kgrid.ylimits
      for key in self.kgrid.kernels:
         ix0, ix1, iy0, iy1 = self.get_pixelgrid_edge(*key)
         kvalue = self.kgrid.kernels[key].kernel
         if normalize:
            # normalize the maximum to 1.0
            if np.max(kvalue) > 0:
               kvalue = kvalue / np.max(kvalue)
            else:
               kvalue = kvalue * 0.
         assert not np.isnan(kvalue).any(), "nan is in kernel (%d,%d)" % (key[0],key[1])
         pixgrid[ix0:ix1,iy0:iy1] = kvalue
      fig = plt.figure()
      ax = fig.add_subplot(111)
      ax.imshow(pixgrid.T, vmin=0., vmax=np.max(pixgrid)/vmax_factor,
                extent=(0,Nx,0,Ny))
      ax.set_xticks(range(Nx+1)[::2*self.xwidth_kernel])
      ax.set_xticklabels(np.linspace(xlimits[0], xlimits[1], self.kgrid.Nx+1)[::2])
      ax.set_yticks(range(Ny+1)[::2*self.ywidth_kernel])
      ax.set_yticklabels(np.linspace(ylimits[0], ylimits[1], self.kgrid.Ny+1)[::2])
      ax.set_title(os.path.split(self.kgrid.catalog)[-1])
      self.pixgrid = pixgrid
      self.xcenters = self.get_xcenters()
      self.ycenters = self.get_ycenters()
      xc, yc = np.meshgrid(self.xcenters, self.ycenters)
      ax.scatter(xc, yc, marker="+", s=ms, color="white")
      ax.set_xlim(0, pixgrid.shape[0])
      ax.set_ylim(0, pixgrid.shape[1])
      return ax
