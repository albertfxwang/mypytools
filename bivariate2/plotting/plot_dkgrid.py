#!/usr/bin/env python

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import pickle_util as pu
from pygoods import *
from my_mplfonts import Helvetica
from bivariate2 import dropout_kernels as dk
from scipy.interpolate import interp1d, RectBivariateSpline, UnivariateSpline
import pickle_util as pu
import os
import pmscolors

class PlotDKGrid(object):
   """
   A class for making plots of the dropout-selection kernel grid.
   """
   def __init__(self, dkgrid, field):
      self.dkgrid = dkgrid
      self.field = field
      # self.grid are for showing things like selection/detection completeness
      # in each (M, logR) bin
      self.grid = np.zeros((dkgrid.Nx, dkgrid.Ny))
      if dkgrid.interpolated:
         self.zlimits = dkgrid.zlimits_old
         self.dz = dkgrid.dz_old
         self.Nz = dkgrid.Nz_old
      else:
         self.zlimits = dkgrid.zlimits
         self.dz = dkgrid.dz
         self.Nz = dkgrid.Nz
      self.grid_zedges = np.linspace(self.zlimits[0], self.zlimits[1], 
                                     num=self.Nz+1)
      self.grid_shape = self.grid.shape
      self.vmin = 0.
      self.vmax = 1.
      self.xname = dkgrid.xname
      self.yname = dkgrid.yname
      self.xticks = dkgrid.xedges()[::2]
      self.yticks = dkgrid.yedges()[::2]
      self.title = ""

   def combine_kernels(self, m1500_lim, logr_lim):
      """
      Combine the statistics (Ninput, Ndetect, Nselect) from all kernels
      within the specified range.
      """
      assert (m1500_lim[1] > m1500_lim[0])
      assert (logr_lim[1] > logr_lim[0])
      # Find the range of kernel indices
      i0 = self.dkgrid.find_bin(m1500_lim[0], logr_lim[0])
      i1 = self.dkgrid.find_bin(m1500_lim[1], logr_lim[1])
      # Now add all kernels within the range to klist
      klist = []
      for ix in range(i0[0], i1[0]+1):
         for iy in range(i0[1], i1[1]+1):
            klist += [self.dkgrid.kernels[(ix, iy)]]
      # Figure out what the redshift range is
      # Assume that self.dkgrid.interpolated == True
      zrange = np.linspace(*self.zlimits, num=self.Nz)
      zcenters = zrange + self.dz / 2.
      n_zbins = len(zcenters)
      # Add the statistics from all kernels
      Ninput = np.zeros(n_zbins, 'int')
      Ndetect = np.zeros(n_zbins, 'int')
      Nselect = np.zeros(n_zbins, 'int')
      for k in klist:
         Ninput = Ninput + k.Ninput
         Ndetect = Ndetect + k.Ndetect
         Nselect = Nselect + k.Nselect
      return zcenters, Ninput, Ndetect, Nselect

   def plot_Pz_single(self, m1500_lim, logr_lim, ax=None, **plot_kwargs):
      """
      Plot P(z) for a combination of kernels within the specified limits.
      """
      zcenters0, N_in, N_det, N_sel = self.combine_kernels(m1500_lim, logr_lim)
      # interpolate
      Pz0 = N_sel / np.maximum(N_in, 1).astype('float')
      f = interp1d(zcenters0, Pz0, kind='cubic')
      # evaluate at the new redshift values
      Pz_new = f(self.dkgrid.zcenters())
      Pz_new = np.maximum(Pz_new, 0.)
      # Now plot
      if ax == None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
      ax.plot(self.dkgrid.zcenters(), Pz_new, **plot_kwargs)
      ax.set_ylim(0., 1.)
      return ax

   def show_grid(self, ax=None, vmin=-1., vmax=-1.):
      """
      Show self.grid
      """
      if ax == None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
      if vmin < 0: vmin = self.vmin
      if vmax < 0: vmax = self.vmax
      m = ax.imshow(self.grid.T, 
                    extent=(0,self.grid_shape[0],0,self.grid_shape[1]),
                    vmin=vmin, vmax=vmax)
      ax.set_title(self.title, size=20)
      ax.set_xticks(np.arange(self.dkgrid.Nx)[::2])
      ax.set_xticklabels(map(lambda x:'%.1f'%x, self.xticks))
      ax.set_yticks(np.arange(self.dkgrid.Ny)[::2])
      ax.set_yticklabels(map(lambda y:'%.1f'%y, self.yticks))
      ax.set_xlabel(self.xname, size=16)
      ax.set_ylabel(self.yname, size=16)
      plt.colorbar(m)


   def selection_completeness(self, z0=None, z1=None, detect_only=False, show=True, ax=None, vmax=-1):
      if z0 == None:
         z0 = self.zlimits[0]
      if z1 == None:
         z1 = self.zlimits[1]
      assert (z1 - z0) > self.dz, "z1 must be at least z0 + self.dz"
      iz0 = np.searchsorted(self.grid_zedges, z0) - 1
      iz0 = np.maximum(iz0, 0)
      iz0 = np.minimum(iz0, self.Nz-1)
      iz1 = np.searchsorted(self.grid_zedges, z1) - 1
      iz1 = np.maximum(iz1, 0)
      iz1 = np.minimum(iz1, self.Nz-1)
      for k in self.dkgrid.kernels:
         kern = self.dkgrid.kernels[k]
         if detect_only:
            N_input = kern.Ndetect[iz0:iz1+1].sum()
         else:
            N_input = kern.Ninput[iz0:iz1+1].sum()
         N_select = kern.Nselect[iz0:iz1+1].sum()
         comp = float(N_select) / float(np.maximum(N_input, 1))
         self.grid[k[0], k[1]] = comp
      self.title = "%s\nSelection Completeness z=[%.1f,%.1f] in %s" % (self.dkgrid.filename, z0, z1, self.field.upper())
      self.vmin = 0.
      if vmax < 0:
         self.vmax = self.grid.max()
      else:
         self.vmax = vmax
      if show:
         self.show_grid(ax=ax)

## Define arguments for plotting P(z) grid
## The panel grid should have this arrangement:
## -------------------
##     1    |    2   |
## -------------------
##     3    |    4   |
## -------------------
dkgrid_dir = '/Users/khuang/Dropbox/Research/bivariate/bivariate_fit/dropsim_kernels'
dkgrid3 = pu.load_pickle(os.path.join(dkgrid_dir, 'udrops_gds_deep_kgrid_m2h_140815.p'))
dkgrid4 = pu.load_pickle(os.path.join(dkgrid_dir, 'bdrops_goods_kgrid_140317.p'))
dkgrid5 = pu.load_pickle(os.path.join(dkgrid_dir, 'vdrops_goods_kgrid_140402.p'))
dkgrid6 = pu.load_pickle(os.path.join(dkgrid_dir, 'idrops_gds_deep_kgrid_140228.p'))
dkgrid_list = [dkgrid3, dkgrid4, dkgrid5, dkgrid6]
pms = pmscolors.pmscolors()
colors = map(pms, ['Blue Purples', 'Periwinkle', 'Olive Green', 'Bright Red'])
m1500_lims = np.array([[-23.,-21.], [-21.,-19.], [-23.,-21.], [-21.,-19]])
logr_lims = np.array([[-0.5,0.5], [-0.5,0.5], [-1.5,0.5], [-1.5,0.5]])

def plot_dkgrid_multiple(dkgrids=dkgrid_list, colors=colors, m1500_lims=m1500_lims, logr_lims=logr_lims, nrows=2, zrange=[2.,7.], axes_pad=0.1):
   """
   Plot multiple P(z) kernels on a grid of plots. All P(z) kernels are supposed
   to be for different redshift ranges but in the same field.
   dkgrids     -- a list of P(z) kernel grids
   colors      -- colors for each curve (same across all panels)
   ## len(dkgrids) == len(colors) is the number of kernel grids shown in each panel.
   m1500_lims  -- the M_1500 limits in each panel
   logr_lims   -- the log(R) limits (in arcsec) in each panel
   ## len(m1500_lims) == len(logr_lims) is the number of panels.
   nrows       -- the number of rows for the grid of plots
   zrange      -- the redshift range shown in each panel
   """
   fig = plt.figure(figsize=(11,9))
   ngrids = len(m1500_lims)
   ncols = ngrids / nrows
   print ngrids % nrows
   if ngrids % nrows > 0:
      ncols += 1
   # First initialize each kernel grids
   plots = []
   for i in range(len(dkgrids)):
      plots += [PlotDKGrid(dkgrids[i], 'DEEP')]
   grid = AxesGrid(fig, 111, 
                   nrows_ncols=(nrows, ncols), 
                   axes_pad=axes_pad,
                   share_all=True, 
                   label_mode="L",
                   aspect=False)
   for j in range(len(m1500_lims)):
      ax_j = grid[j]
      for i in range(len(dkgrids)):
         plots[i].plot_Pz_single(m1500_lims[j], logr_lims[j], ax=ax_j,
                                 lw=2, color=colors[i])
      bintext = r"$%.1f \leq M_{1500} \leq %.1f$" % tuple(m1500_lims[j])
      bintext = bintext + '\n'
      bintext = bintext + r"$%.1f \leq \log R_e \leq %.1f$" % tuple(logr_lims[j])
      ax_j.text(0.95, 0.95, bintext, ha='right', va='top', size='large',
                transform=ax_j.transAxes, 
                bbox=dict(boxstyle='round', facecolor='none'))
      ax_j.set_xlim(zrange)
      ax_j.set_ylim(0., 1.1)
      ax_j.set_xlabel('Redshift')
      ax_j.set_xticks(np.linspace(3., 6., 4))
      ax_j.set_ylabel('P(z)')
      ax_j.set_yticks(np.linspace(0.25, 1.0, 4))
   plt.draw()
   return grid

