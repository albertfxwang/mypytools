#!/usr/bin/env python

from numpy import *
import numpy as np
#import mlutil
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_grid import Grid, ImageGrid
import os, glob
from my_mplfonts import Helvetica
import plot_comp


M0grid = array([22.,24.,26.,22.,24.,26.,22.,24.,26.])
r0grid = array([1.2,1.2,1.2,0.6,0.6,0.6,0.,0.,0.])
x0grid_def = array([M0grid, r0grid])
M0grid2 = array([22.,26.,22.,26.])
r0grid2 = array([1.2,1.2,0.,0.])
x0grid_def2 = array([M0grid2, r0grid2])
M0grid_udf = tile([23.,25.,27.], 3)
r0grid_udf = repeat([1.2,0.6,0.], 3)
x0grid_udf = array([M0grid_udf, r0grid_udf])

def kernel_x0(kgrid,mag0,logre0,dx=0.01):
   # search for kernel by matching bin lower bounds
   mk = None
   for k in kgrid.kernels:
      kern = kgrid.kernels[k]
      if (abs(kern.x0_cell[0]-mag0)<=dx) & (abs(kern.x0_cell[1]-logre0)<=dx):
         mk = kern
   #print mk, mag0, logre0
   return mk

def show_kernel(kgrid, mag0, logr0):
   # show the kernel image given kgrid.kernel instance
   fig = plt.figure()
   ax = plt.subplot(111)
   k = kernel_x0(kgrid, mag0, logr0)
   ax.imshow(k.kernel.swapaxes(0,1), origin='lower', cmap=mpl.cm.gnuplot)
   ax.set_xlim(0, 128); ax.set_ylim(0,128)
   ax.set_xticks([14,64,114])
   ax.set_xticklabels([-1,0,1],size=14,color='black')
   ax.set_yticks([14,64,114])
   ax.set_yticklabels([-1,0,1],size=14,color='black')
   ax.set_xlabel(r'$\Delta$mag', size=14)
   ax.set_ylabel(r'$\Delta\log_{10} R_e$', size=14)
   ax.plot([0,128],[64,64],':',c='white')
   ax.plot([64,64],[0,128],':',c='white')
   return ax, fig

def kgrid_allbins(kgrid, vmax_factor=3.):
   """
   Show the GALFIT transfer function kernels within each (mag, logRe) bin, 
   with the same coordinates as the RL distribution model.
   """
   fig = plt.figure()
   bw = kgrid.cell_widths
   bw_pix = around(kgrid.cell_widths / kgrid.pixdx).astype('int')
   # the extent of each kernel
   ksize = kgrid.xmax / kgrid.pixdx
   kcenter = around(ksize / 2).astype('int')
   print "Limits of kernel grid:", kgrid.limits
   npix = (kgrid.limits[:,1] - kgrid.limits[:,0]) / kgrid.pixdx
   npix = around(npix).astype('int')
   # kmodel is the model that shows all the kernels
   kmodel = zeros(npix, 'float')
   vmax = 0.
   for i in range(kgrid.ncells[0]):
      for j in range(kgrid.ncells[1]):
         x0 = i * bw_pix[0]
         x1 = (i+1) * bw_pix[0]
         y0 = j * bw_pix[1]
         y1 = (j+1) * bw_pix[1]
         xcenter = kcenter[0]
         ycenter = kcenter[1]
         k = kgrid.kernels[(i,j)].kernel
         kmax = max(k.ravel())
         if kmax > vmax:
            vmax = kmax
         kmodel[x0:x1,y0:y1] = k[xcenter-bw_pix[0]/2:xcenter+(bw_pix[0]+1)/2,\
                                 ycenter-bw_pix[1]/2:ycenter+(bw_pix[1]+1)/2]
   cax = plt.imshow(kmodel.swapaxes(0,1), vmin=0., vmax=vmax/vmax_factor,
                    aspect=bw[0]/bw[1], extent=(0,npix[0],0,npix[1]))
   plt.xticks(arange(0, npix[0], bw_pix[0])[::2], 
              arange(kgrid.limits[0][0], kgrid.limits[0][1], bw[0])[::2])
   plt.yticks(arange(0, npix[1], bw_pix[1])[::2],
              arange(kgrid.limits[1][0], kgrid.limits[1][1], bw[1])[::2])
   plt.xlabel('magnitude',size=14)
   plt.ylabel('log10(Re)',size=14)
   plt.title(kgrid.filename, size=20)
   plt.colorbar(cax)
   return kmodel



def tfkernels_bins(kgrid, x0grid=x0grid_def, title="",
                   colors=['red','blue','green'], pixscale=0.03):
   # Show the kernels on a 3 by 3 subplot grid
   # the lower bounds of each kernel shown is specified by x0grid[0], x0grid[1]
   # the order is in the order of subplot (x0grid[0][0] is upper left most, x0grid[0][-1]
   # is bottom right most
   if len(x0grid[0]) != 9: raise ValueError, "please give x0 for a 3 by 3 subplot grid"
   # x0grid: from top to bottom, from left to right
   # 1  2  3
   # 4  5  6
   # 7  8  9
   fig = plt.figure(figsize=(10,9))
   imgrid = ImageGrid(fig, 111, nrows_ncols=(3,3), axes_pad=0.03, share_all=True)
   vmin = 1.e-6
   vmax = 2.e-3
   #prop = {'cmap':'hot','vmin':vmin,'vmax':vmax}
   prop = {'cmap':'hot','origin':'lower'}
   for i in range(9):
      #ax = plt.subplot(3,3,i+1)
      ax = imgrid[i]
      k = kernel_x0(kgrid,x0grid[0][i],x0grid[1][i])
      print "k", k.index
      kshape = shape(k.kernel)
      #ax.imshow(k.kernel.swapaxes(0,1),**prop)
      kmax = max(k.kernel.ravel())
      kseq = [0.9*kmax, 0.5*kmax, 0.1*kmax]
      ax.contour(k.kernel.swapaxes(0,1), kseq, colors=colors, linewidths=2.0)
      ax.text(0.1, 0.1, '[%.1f, %.1f]' % (k.x0_cell[0], k.x1_cell[0]), 
              transform=ax.transAxes,
              font_properties=Helvetica(16))
      r0bin_as = pixscale * (10.**k.x0_cell[1])  # r0 of this bin in arcsec
      r1bin_as = pixscale * (10.**k.x1_cell[1])
      ax.text(0.1, 0.9, '[%.2f\", %.2f\"]' % (r0bin_as, r1bin_as), 
              transform=ax.transAxes,
              font_properties=Helvetica(16))
      #if i == 3:
      #   plt.ylabel('log10(Re_out/Re_in) [pixels]', size=18)
      #if i == 7:
      #   plt.xlabel('mag_out - mag_in', size=18)
      ax.plot([kshape[1]/2],[kshape[0]/2],'+', ms=20, mew=1.0, c='black')
      ax.set_xlim(0, kshape[1]); ax.set_ylim(0,kshape[0])
      
      #plt.title('$%.1f<m_{in}<%.1f$' % (k.x0_cell[0],k.x1_cell[0]),size=11)
      #plt.ylabel('$%.1f<\log_{10}(R_e)<%.1f$' % (k.x0_cell[1],k.x1_cell[1]),size=11)
   
   for i in range(9):
      imgrid[i].set_yticks([14,64,114])
      #imgrid[i].tick_params(axis='y', color='white')
      if i in [0,3,6]:
         imgrid[i].set_yticklabels([-1,0,1],size=14,color='black')     
         imgrid[i].axis["left"].label.set_text(r"$\delta \log R$")
         imgrid[i].axis["left"].label.set_size(20)
      if i in [2,5,8]:
         imgrid[i].axis["right"].label.set_text(r"$\delta \log R$")     
         imgrid[i].axis["right"].label.set_size(20)   
      #   ax2 = imgrid[i].twinx()
      #   yd1, yd2 = imgrid[i].yaxis.get_data_interval()
      #   ax2.yaxis.set_data_interval(yd1,yd2)
      #   yl1, yl2 = imgrid[i].get_ylim()
      #   ax2.set_ylim(yl1, yl2)
      #   #ax2.set_yticks([14,64,114])
      #   #ax2.set_yticklabels([-1,0,1],size=12,color='black')
      #   #ax2.tick_params(axis='y', color='white')
      #   ax2.yaxis.tick_right()
      #   ax2.set_ylabel(r"$\delta \log R$")
      #imgrid[i].yaxis.label.set_color('white')
      #plt.xlabel('$\delta m$', size=8)
      #plt.ylabel('$\delta \log R$', size=8)
      imgrid[i].set_xticks([14,64,114])
      imgrid[i].set_xticklabels([-1,0,1],size=16,color='black')
      #imgrid[i].tick_params(axis='x', color='white')
      imgrid[i].axis["bottom"].label.set_text(r"$\delta m$")
      imgrid[i].axis["bottom"].label.set_size(20)
      
   if title == "":
      title = kgrid.filename
   plt.suptitle(title,size=24)

   return 0

class plot_gf_kgrid_comp(plot_comp.plot_comp):
   def __init__(self, kgrid):
      for k in kgrid.__dict__.keys():
         setattr(self, k, getattr(kgrid, k))
      self.compgrid = zeros(kgrid.ncells, 'float')
      self.figure = plt.figure()
      self.ax = self.figure.add_subplot(111)
   
   def newfigure(self):
      self.figure = plt.figure()
      self.ax = self.figure.add_subplot(111)

   def show_grid(self, celldata, vmin=0., vmax=1.):
      cax = self.ax.imshow(celldata.swapaxes(0,1), vmin=vmin, vmax=vmax, 
                           extent=(0,self.ncells[0],0,self.ncells[1]))
      plt.colorbar(cax)
      self.ax.set_xticks(np.arange(self.ncells[1]+1))
      mag_ticks = np.arange(self.limits[0,0], self.limits[0,1]+self.cell_widths[0], self.cell_widths[0])
      mag_ticks = np.array(map(lambda x: "%.1f" % x, mag_ticks))
      mag_ticks[1::2] = ""
      self.ax.set_xticklabels(mag_ticks)
      self.ax.set_yticks(np.arange(self.ncells[0]+1))
      logre_ticks = np.arange(self.limits[1,0], self.limits[1,1]+self.cell_widths[1], self.cell_widths[1])
      self.ax.set_yticklabels(map(lambda x: "%.1f" % x, logre_ticks))
      self.ax.set_xlabel('magnitude', font_properties=Helvetica(16))
      self.ax.set_ylabel('log10(Re)', font_properties=Helvetica(16))
      plt.xlim(0, self.ncells[0])
      plt.ylim(0, self.ncells[1])

   def plot_gf_comp(self):
      """
      Plot the completeness due to poor GALFIT results.
      """
      for i in range(self.ncells[0]):
         for j in range(self.ncells[1]):
            kernel = self.kernels[(i,j)]
            self.compgrid[i,j] = float(kernel.nmodel_cell) / np.maximum(float(kernel.ndetect_cell), 1.)
            if kernel.nmodel_cell <= 4:
               self.compgrid[i,j] = 0.
      self.show_grid(self.compgrid)
      self.ax.set_title('GALFIT completeness\n%s' % (self.filename), 
                        font_properties=Helvetica(16))

   def plot_gf_ninput(self):
      self.ninput_grid = np.zeros(self.ncells)
      for i in range(self.ncells[0]):
         for j in range(self.ncells[1]):
            kernel = self.kernels[(i,j)]
            self.ninput_grid[i,j] = kernel.ninput_cell
      self.show_grid(self.ninput_grid, vmax=self.ninput_grid.max()*1.1)
      self.ax.set_title('Input Number of Galaxies\n%s' % self.filename, 
                        font_properties=Helvetica(16))

def tfkernels_bins2(kgrid1, kgrid2, x0grid=x0grid_def2, title="",
   colors=['red','blue','green'], pixscale=0.03):
   # show the 2x2 kernel images of both GALFIT & SExtractor kernel grids
   # display the two kernel grids side by side
   #mpl.rcParams['xtick.labelsize'] = 20
   #mpl.rcParams['ytick.labelsize'] = 20
   #mpl.rcParams['axes.labelsize']  = 22
   #mpl.rcParams['axes.titlesize'] = 24
   if len(x0grid[0]) != 4: raise ValueError, "please give x0 for a 2 by 2 subplot grid"
   fig = plt.figure(figsize=(14,6))
   axes = []
   def display_kernel(ax, k):
      prop = {'cmap':'hot','origin':'lower'}
      #ax.imshow(k.kernel.swapaxes(0,1), **prop)
      kmax = max(k.kernel.ravel())
      kseq = [0.9*kmax, 0.5*kmax, 0.1*kmax]
      ax.contour(k.kernel.swapaxes(0,1), kseq, colors=colors, linewidths=2.0)
      return ax
   gap = 0.005
   coords_GF = [[0.05,0.25+gap,0.05,0.25+gap], [0.5+gap,0.5+gap,0.1,0.1]]
   # start with GALFIT kernels (on the left)
   for i in range(4):
      ax0 = fig.add_axes([coords_GF[0][i], coords_GF[1][i], 0.2, 0.4])
      k0 = kernel_x0(kgrid1, x0grid[0][i], x0grid[1][i])
      ax0 = display_kernel(ax0, k0)
      axes += [ax0]
      ax0.text(0.1, 0.1, '[%.1f, %.1f]' % (k0.x0_cell[0], k0.x1_cell[0]),
         color='black', transform=ax0.transAxes,
         font_properties=Helvetica(18))
      r0bin_as = pixscale * (10.**k0.x0_cell[1])  # r0 of this bin in arcsec
      r1bin_as = pixscale * (10.**k0.x1_cell[1])
      ax0.text(0.1, 0.9, '[%.2f\", %.2f\"]' % (r0bin_as, r1bin_as),
         color='black', transform=ax0.transAxes,
         font_properties=Helvetica(18))
      ax0.plot([64.],[64.],'+', c='black', ms=20, mew=1.0)
      ax0.set_xlim(0, 128); ax0.set_ylim(0,128)
      ax0.set_yticks([14,64,114])
      ax0.tick_params(axis='y', color='black')
      if i in [0,2]:
         ax0.set_yticklabels([-1,0,1], color='black',
                             font_properties=Helvetica(16))     
         ax0.set_ylabel(r"$\delta \log R_e$", font_properties=Helvetica(24))
      else:
         ax0.set_yticklabels([])
         ax0.set_ylabel("")
      #if i in [1,3]:
      #   ax0.(r"$\delta \log R$")
      
      ax0.set_xticks([14,64,114])
      ax0.tick_params(axis='x', color='black')
      if i in [2,3]:
         ax0.set_xticklabels([-1,0,1],color='black',
                             font_properties=Helvetica(16))
         ax0.set_xlabel(r"$\delta m$", font_properties=Helvetica(24))
      else:
         ax0.set_xticklabels([])
         ax0.set_xlabel("")
   
   coords_SE = [[0.55,0.75+gap,0.55,0.75+gap], [0.5+gap,0.5+gap,0.1,0.1]]
   # now with SExtractor kernels (on the right)
   for i in range(4):
      ax0 = fig.add_axes([coords_SE[0][i], coords_SE[1][i], 0.2, 0.4])
      k0 = kernel_x0(kgrid2, x0grid[0][i], x0grid[1][i])
      ax0 = display_kernel(ax0, k0)
      axes += [ax0]
      ax0.text(0.1, 0.1, '[%.1f, %.1f]' % (k0.x0_cell[0], k0.x1_cell[0]),
         color='black', transform=ax0.transAxes,
         font_properties=Helvetica(18))
      r0bin_as = pixscale * (10.**k0.x0_cell[1])  # r0 of this bin in arcsec
      r1bin_as = pixscale * (10.**k0.x1_cell[1])
      ax0.text(0.1, 0.9, '[%.2f\", %.2f\"]' % (r0bin_as, r1bin_as),
         color='black', transform=ax0.transAxes,
         font_properties=Helvetica(18))
      ax0.plot([64.],[64.],'+', c='black', ms=20, mew=1.0)
      ax0.set_xlim(0, 128); ax0.set_ylim(0,128)
      ax0.set_yticks([14,64,114])
      ax0.tick_params(axis='y', color='black')
      if i in [0,2]:
         ax0.set_yticklabels([-1,0,1],color='black',
                             font_properties=Helvetica(16))     
         ax0.set_ylabel(r"$\delta \log R_e$", font_properties=Helvetica(24))
      else:
         ax0.set_yticklabels([])
         ax0.set_ylabel("")
      #if i in [1,3]:
      #   ax0.(r"$\delta \log R$")
      
      ax0.set_xticks([14,64,114])
      ax0.tick_params(axis='x', color='black')
      if i in [2,3]:
         ax0.set_xticklabels([-1,0,1],color='black',
                             font_properties=Helvetica(16))
         ax0.set_xlabel(r"$\delta m$", font_properties=Helvetica(24))
      else:
         ax0.set_xticklabels([])
         ax0.set_xlabel("")
   fig.text(0.22, 0.92, 'GALFIT', font_properties=Helvetica(24))
   fig.text(0.72, 0.92, 'SExtractor', font_properties=Helvetica(24))
   return axes


def nmodel_tot(kgrid,output="",wfits=False,vmin=5.,newfig=True):
   # shows how many points from simulation are used in each bin
   npix=(kgrid.limits[:,1]-kgrid.limits[:,0])/kgrid.pixdx
   ngrid=zeros(npix,"int")
   nx, ny = kgrid.ncells
   nmax = 0

   for i in range(nx):
      for j in range(ny):
         kern=kgrid.kernels[(i,j)]
         x0=kern.x0_cell
         x1=kern.x1_cell
         index0 = (x0-kgrid.limits[:,0])/kgrid.pixdx
         index0 = around(index0).astype('int')
         index1 = (x1-kgrid.limits[:,0])/kgrid.pixdx
         index1 = around(index1).astype('int')
         ngrid[index0[0]:index1[0],index0[1]:index1[1]]=kern.nmodel_cell
         if kern.nmodel_cell > nmax:
            nmax = kern.nmodel_cell
   if newfig: plt.figure()
   plt.imshow(ngrid.swapaxes(0,1),origin='lower',vmin=vmin,vmax=nmax)
   plt.xticks(arange(0,nx,2)*25,arange(0,nx,2)*0.5+kgrid.limits[0][0])
   plt.yticks(arange(0,ny,2)*10,arange(0,ny,2)*0.2+kgrid.limits[1][0])
   if wfits:
      mlutil.writefits2d(ngrid.swapaxes(0,1),output)
   return 0


def kimages(kgrid):
   nx,ny = kgrid.nbins
   for i in range(nx):
      for j in range(ny):
         kern = kgrid.kernels[(i,j)]
         sp_row = 10 - j
         sp_col = i + 1
         plt.subplot(ny,nx,10*(sp_row-1)+sp_col)
         plt.imshow(kern.kernel.swapaxes(0,1),origin='lower',vmax=1.e-3,vmin=0.)
         plt.xticks([])
         plt.yticks([])
   plt.subplots_adjust(hspace=0.,wspace=0.)
   plt.show()
   return 0


def kern_info(key,kgrid):
    kern = kgrid.kernels[key]
    nmodel = kern.nmodel-kern.nflagged
    maxval = max(kern.kernel.ravel())
    kernsum = sum(kern.kernel.ravel())
    completeness = kern.completeness
    #print "%10s  %5s  %7s  %7s  %7s" % ("key","nmod","max","sum","complt")
    print "%10s  %5d  %1.5g  %1.5g  %1.5g" % (str(key),nmodel,
           maxval,kernsum,completeness)


def compmap(kgrid,vmin=0.0,vmax=1.0,ax=None):
   # shows how many points from simulation are used in each bin
   cgrid=zeros(kgrid.ncells)

   for k in kgrid.kernels:
      i = k[0]
      j = k[1]
      kern=kgrid.kernels[(i,j)]
      cgrid[i,j] = kern.completeness
   if ax==None: 
      fig = plt.figure()
      ax = fig.add_subplot(111)
   cax = ax.imshow(cgrid.swapaxes(0,1),vmin=vmin,vmax=vmax,cmap=mpl.cm.hot,
                   extent=(0,kgrid.ncells[0]+1,0,kgrid.ncells[1]+1))
   plt.colorbar(cax, orientation='vertical')
   ax.set_xticks(arange(0, kgrid.ncells[0]+1)[::2])
   ax.set_xticklabels(arange(kgrid.limits[0][0],kgrid.limits[0][1]+1,
                      kgrid.cell_widths[0])[::2]) 
   ax.set_yticks(arange(0,kgrid.ncells[1]+1)[::2])
   ax.set_yticklabels(arange(kgrid.limits[1][0],kgrid.limits[1][1]+1,
                      kgrid.cell_widths[1])[::2])
   ax.set_title(kgrid.filename)
   return cgrid


def kernel_nummap(kgrid):
   nx = mnt(round((kgrid.limits[0,1]-kgrid.limits[0,0])/kgrid.pixdx[0]))
   ny = int(round((kgrid.limits[1,1]-kgrid.limits[1,0])/kgrid.pixdx[1]))
   nummap = zeros((nx,ny))
   for k in kgrid.kernels:
      kern = kgrid.kernels[k]
      ninput = kern.ninput_bin
      #if comp < 0: comp = 0.
      i0 = around((kern.x0_cell-kgrid.limits[:,0])/kgrid.pixdx).astype('int')
      i1 = around((kern.x1_cell-kgrid.limits[:,0])/kgrid.pixdx).astype('int')
      nummap[i0[0]:i1[0],i0[1]:i1[1]] = ninput
   plt.imshow(nummap.swapaxes(0,1),origin='lower')
   plt.xticks(arange(0,300,50),arange(22,28))
   plt.xlabel('magnitude')
   plt.yticks([0,50,100],[0.03,0.3,3.0])
   plt.ylabel('Re (arcsec)')
   plt.colorbar(orientation='horizontal')
   plt.show()
   return 0

