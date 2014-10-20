#!/usr/bin/env python

### For bivariate2

import numpy as np
import matplotlib.pyplot as plt
import plot_RLDist as pRL
from mpl_toolkits.axes_grid1 import make_axes_locatable
# import plot_gfmeasure as pgm
from bivariate2 import galaxy_samples as gs

__created__ = "2014/03/10"
__author__ = "Kuang-Han Huang"
# General philosophy: write a class that handles a field (a unique combination
# of depth in detection, color selection, and GALFIT measurement band). When
# showing all fields in one plot, divide the fields into tiers of depth in 
# the GALFIT measurement band.

class PlotRLDistFit(gs.LBGSample):
   def __init__(self, pars, paramfile, filtername, sample_name, **kwargs):
      self.sample_name = sample_name
      self.parameters = pars
      # Initialize all the correction kernels and everything for this LBG
      # Sample
      super(PlotRLDistFit, self).__init__(paramfile, filtername, **kwargs)
      # Now evaluate the distribution at the best-fit parameters
      self.loglikelihood_tot(pars, floor=1.e-50)
      self.limits = np.zeros((2, 2))
      # For plotting GALFIT-measured mag & Re
      self.titlesize = 20
      self.axlabelsize = 16
      self.ticklabelsize=14

   def set_labels(self, field, ax, dxlabel=1.0, dylabel=0.4):
      M = self.RLDist_factories[field].RLDist
      limits = np.array([M.xlimits, M.ylimits])
      # nxlabel = int(round((limits[0,1] - limits[0,0]) / dxlabel)) + 1
      # nylabel = int(round((limits[1,1] - limits[1,0]) / dylabel)) + 1
      # xticks = np.linspace(0, M.value.shape[0], num=nxlabel)
      # yticks = np.linspace(0, M.value.shape[1], num=nylabel)
      xticks = np.arange(0, M.value.shape[0], int(round(dxlabel/M.dx)))
      yticks = np.arange(0, M.value.shape[1], int(round(dylabel/M.dy)))
      xlabels = map(lambda x: '%.1f' % (x*self.dx+limits[0,0]), xticks)
      ylabels = map(lambda y: '%.1f' % (y*self.dy+limits[1,0]), yticks)
      ax.set_xticks(xticks)
      ax.set_yticks(yticks)
      ax.set_xticklabels(xlabels, size=self.ticklabelsize)
      ax.set_yticklabels(ylabels, size=self.ticklabelsize)
      ax.set_xlabel('%s magnitude' % self.filtername, size=self.axlabelsize)
      ax.set_ylabel('%s effective radius [arcsec]' % (self.filtername), 
                    size=self.axlabelsize)
      return ax

   def set_LF_labels(self, field, ax, dmag=1.0):
      M = self.RLDist_factories[field].RLDist
      nlabel = int(round((M.xlimits[1] - M.xlimits[0]) / dmag)) + 1
      # magticks = np.linspace(0, M.value.shape[0], num=nlabel)
      magticks = np.arange(0, M.value.shape[0], int(round(dmag/M.dx)))
      maglabels = map(lambda x: '%.1f' % (x*self.dx+M.xlimits[0]), magticks)
      ax.set_xticks(magticks)
      ax.set_xticklabels(maglabels, size=self.ticklabelsize)
      ax.set_xlabel('%s magnitude' % self.filtername, size=self.axlabelsize)
      return ax

   def set_SD_labels(self, field, ax, dlogr=0.4, direction='x'):
      M = self.RLDist_factories[field].RLDist
      nlabel = int(round((M.ylimits[1] - M.ylimits[0]) / dlogr)) + 1
      # rticks = np.linspace(0, M.value.shape[1], num=nlabel)
      rticks = np.arange(0, M.value.shape[1], int(round(dlogr/M.dy)))
      rlabels = map(lambda y: '%.1f' % (y*self.dy+M.ylimits[0]), rticks)
      if direction == 'x':
         ax.set_xticks(rticks)
         ax.set_xticklabels(rlabels, size=self.ticklabelsize)
         ax.set_xlabel('log10(Re) [arcsec]', size=self.axlabelsize)
      else:
         ax.set_yticks(rticks)
         ax.set_yticklabels(rlabels, size=self.ticklabelsize)
         ax.set_ylabel('log10(Re) [arcsec]', size=self.axlabelsize)
      return ax

   def show_RLDist(self, field_list=None, ax=None, cbar=False, dxlabel=1.0, dylabel=0.4, recalc=False):
      # Show the best-fit, corrected RL distribution
      if ax == None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
      if recalc:
         self.loglikelihood_tot(self.parameters, floor=1.e-50)
      if field_list == None:
         field_list = self.fields
      self.axCent = ax
      self.divider = make_axes_locatable(ax)
      # Combine the calculated RL distribution from several fields and show 
      # in one plot.
      # Make sure that the RL distributions for the fields to be combined have
      # the same array shape
      f0 = field_list[0]
      M0 = self.RLDist_factories[f0].RLDist
      dist_array = np.zeros(M0.value.shape)
      dx = 0.; dy = 0.
      if field_list == None:
         field_list = self.fields
      field_names = '+'.join(field_list).upper()
      for i in range(len(field_list)):
         M = self.RLDist_factories[field_list[i]].RLDist
         # note that here M is multiplied by (self.dx * self.dy)
         dist_array = dist_array + M.value
         self.limits[0] = M0.xlimits
         self.limits[1] = M0.ylimits
         self.dx = M0.dx
         self.dy = M0.dy
         self.xcenters = (M0.xcenters() - M0.xlimits[0]) / self.dx
         self.ycenters = (M0.ycenters() - M0.ylimits[0]) / self.dy
         assert (self.limits[0]==M.xlimits).all() & (self.limits[1]==M.ylimits).all(), "Fields have incompatible distribution limits."
      # Now show the combined RL likelihood distribution
      cax = ax.imshow(dist_array.T, 
                      extent=(0,dist_array.shape[0],0,dist_array.shape[1]))
      # Make axis labels
      ax = self.set_labels(field_list[0], ax, dxlabel=dxlabel, 
                           dylabel=dylabel)
      self.dist_array = dist_array
      if cbar:
         plt.colorbar(cax)
      return ax

   def show_mag_logRe(self, field_list=None, marker='o', size=4**2, color='black', ax=None, scatter_kw={'edgecolors':'none'}):
      # Show the GALFIT-measured magnitudes and logRe on top of the 
      # distribution
      _xlim = (0, 0)
      _ylim = (0, 0)
      if ax == None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
      else:
         _xlim = ax.get_xlim()
         _ylim = ax.get_ylim()
      if field_list == None:
         field_list = self.fields
      field_names = '+'.join(field_list).upper()
      for i in range(len(field_list)):
         f = field_list[i].lower()
         mag = self.galaxy_samples[f].mag
         mag = self.RLDist_factories[f].RLDist.get_xcoord(mag)
         logRe = self.galaxy_samples[f].logr
         logRe = self.RLDist_factories[f].RLDist.get_ycoord(logRe)
         # Re = 10.**(logRe)
         ax.scatter(mag, logRe, marker=marker, s=size, c=color, **scatter_kw)
         ax.set_xlabel('%s magnitude' % self.filtername.upper(), 
                       size=self.axlabelsize)
         ax.set_ylabel('log10(Re) [arcsec]', size=self.axlabelsize)
         ax.set_title(field_names, size=self.titlesize)
      ax.set_xlabel('%s magnitude' % self.filtername, size=self.axlabelsize)
      ax.set_ylabel('%s effective radius [arcsec]' % (self.filtername), 
                    size=self.axlabelsize)
      if _xlim != (0, 0):
         ax.set_xlim(*_xlim)
      if _ylim != (0, 0):
         ax.set_ylim(*_ylim)
      plt.draw()
      return ax

   def RLDist_boxcar(self, z0, field_list=None, dz=1.0):
      """
      Calculate the expected size-luminosity distribution if the selection 
      function is a boxcar of 100% within [z0-dz/2, z0+dz/2].
      """
      if field_list == None:
         field_list = self.fields
      M0 = self.RLDist_factories[field_list[0]].RLDist
      dist_array = np.zeros(M0.value.shape)
      # calls the RLDistributionFactory.RLdist_boxcar
      for f in field_list:
         factory = self.RLDist_factories[f]
         dist_array = dist_array + factory.RLdist_boxcar(self.parameters, 
                                   z0=z0, dz=dz)
      return dist_array

   def RLDist_dkgrid(self, field_list=None):
      """
      Calculate the combined size-luminosity distribution with only dropout-
      selection corrections.
      """
      if field_list == None:
         field_list = self.fields
      M0 = self.RLDist_factories[field_list[0]].RLDist
      dist_array = np.zeros(M0.value.shape)
      # calls the RLDistributionFactory.RLdist_dropout_selection
      for f in field_list:
         factory = self.RLDist_factories[f]
         DK = self.dropout_kgrids[f]
         M = factory.RLdist_dropout_selection(self.parameters, DK)
         dist_array = dist_array + M.value
      return dist_array

   def show_LF(self, z0, dz=1.0, field_list=None, axLF=None, lfbw=0.2, color='black', recalc=True, addtext=True, textsize=14):
      """
      Plot three different variants of LF:
      - Fully-corrected LF (should trace the observed number counts well)
      - Non-corrected LF with a boxcar P(z) within [z0-dz/2., z0+dz/2.]
      - Partially-corrected LF with realistic P(z) but no GALFIT smearing
      """
      if axLF == None:
         fig = plt.figure()
         axLF = fig.add_subplot(111)
      if field_list == None:
         field_list = self.fields
      if recalc:
         self.loglikelihood_tot(self.parameters, floor=1.e-50)
      self.axLF = axLF
      allmags = np.zeros(0)
      for f in field_list:
         allmags = np.concatenate([allmags, self.galaxy_samples[f].mag])
      M0 = self.RLDist_factories[field_list[0]].RLDist
      nbins = int(round((self.limits[0,1] - self.limits[0,0]) / lfbw))
      n , bins = np.histogram(allmags, 
                 np.linspace(self.limits[0,0], self.limits[0,1], nbins+1))
      nerr = [np.sqrt(n), np.sqrt(n)]
      # set a lower limit to the Poisson error bar if n == 1
      for i in range(len(n)):
         if n[i] == 1: nerr[0][i] = 1.-1.e-3
      xbincenters = (bins[:-1] - self.limits[0,0] + lfbw / 2.) / self.dx
      axLF.errorbar(xbincenters, n, yerr=nerr, fmt='.', ms=14., 
                    mfc=color, ls='None', mec=color, ecolor=color, capsize=6)
      LF = self.dist_array.sum(axis=1) # LF here contains the volume already
      # LFtot = np.sum(LF) * self.dx / lfbw
      # normalize the LF to predict the total number  of points
      normfactor = len(allmags) * lfbw / (np.sum(LF) * self.dx)
      LF = LF * normfactor
      axLF.semilogy(self.xcenters, LF, color=color, nonposy='mask', 
                    label='GALFIT TF', lw=3.)
      # Plot a RL distribution with boxcar P(z)
      dist0 = self.RLDist_boxcar(z0, field_list=field_list, dz=dz)
      LF0 = dist0.sum(axis=1)
      LF0 = LF0 * normfactor * self.phistar
      axLF.semilogy(self.xcenters, LF0, color=color, ls=':', lw=2, 
                    nonposy='mask',label='boxcar P(z)')
      # Plot a RL distribution with just P(z) but no GALFIT correction
      dist1 = self.RLDist_dkgrid(field_list=field_list)
      LF1 = dist1.sum(axis=1)
      LF1 = LF1 * normfactor * self.phistar
      axLF.semilogy(self.xcenters, LF1, color=color, ls='--', lw=2,
                    nonposy='mask',label='w/ dropout sel. kernel')
      axLF.set_xlim(self.axCent.get_xlim())
      axLF.set_xlabel('%s magnitude' % self.filtername.upper())
      axLF.set_ylabel(r'$dN/dm$', size=self.axlabelsize)
      if addtext:
         alpha, mstar = self.parameters[:2]
         phistar = self.phistar
         props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
         parstr = r'$\alpha=%.2f$' % alpha
         parstr = parstr + '\n' + r'$M^*=%.2f$' % mstar
         parstr = parstr + '\n' + r'$\phi^*=%.2f\times 10^{-3}$' % (phistar*1000.   )
         axLF.text(0.05, 0.95, parstr, bbox=props, va='top', ha='left', 
                   transform=axLF.transAxes, size=textsize)
      axLF.set_yticks([1., 1.e1, 1.e2])   
      axLF.set_ylim(1.e-1, max(n)*50.)
      plt.draw()
      return axLF

   def append_LF(self, z0, dz=1.0, field_list=None, axCent=None, divider=None, lfbw=0.2, color='black', recalc=True):
      """
      Append a panel on top of axCent to show the marginalized LF.
      """
      if axCent == None:
         axCent = self.axCent
      if divider == None:
         divider = self.divider
      axLF = divider.append_axes("top", size=1.5, pad=0.05)
      axLF = self.show_LF(z0, dz=dz, field_list=field_list, axLF=axLF, 
                          lfbw=lfbw, color=color, recalc=recalc, 
                          addtext=False)
      axLF.set_xticklabels([])
      axLF.set_xlabel("")
      return axLF

   def append_SD(self, z0, dz=1.0, field_list=None, axCent=None, divider=None, sdbw=0.2, color='black', recalc=True, textsize=14):
      """
      Append a panel on the right of axCent to show the marginalized size 
      distribution.
      Plot three different variants of size distribution:
      - Fully-corrected SD (should trace the observed number counts well)
      - Non-corrected SD with a boxcar P(z) within [z0-dz/2., z0+dz/2.]
      - Partially-corrected SD with realistic P(z) but no GALFIT smearing
      """
      if axCent == None:
         axCent = self.axCent
      if divider == None:
         divider = self.divider
      if field_list == None:
         field_list = self.fields
      if recalc:
         self.loglikelihood_tot(self.parameters, floor=1.e-50)
      axSD = divider.append_axes("right", size=1.35, pad=0.05)
      self.axSD = axSD
      all_logre = np.zeros(0)
      for f in field_list:
         all_logre = np.concatenate([all_logre, self.galaxy_samples[f].logr])
      M0 = self.RLDist_factories[field_list[0]].RLDist
      nbins = int(round((self.limits[1,1] - self.limits[1,0]) / sdbw))
      n , bins = np.histogram(all_logre, 
                 np.linspace(self.limits[1,0], self.limits[1,1], nbins+1))
      nerr = [np.sqrt(n), np.sqrt(n)]
      # set a lower limit to the Poisson error bar if n == 1
      for i in range(len(n)):
         if n[i] == 1: nerr[0][i] = 1.-1.e-3
      ybincenters = (bins[:-1] - self.limits[1,0] + sdbw / 2.) / self.dy
      axSD.errorbar(n, ybincenters, xerr=nerr, fmt='.', ms=14, mfc=color, 
                    ls='None', mec=color, ecolor=color, capsize=6)
      SD = self.dist_array.sum(axis=0)  # fully-corrected size distribution
      normfactor = len(all_logre) * sdbw / (np.sum(SD) * self.dy)
      SD = SD * normfactor
      axSD.semilogx(SD, self.ycenters, color=color, nonposx='mask', 
                    label='GALFIT TF', lw=3.)
      # Plot a RL distribution with boxcar P(z)
      dist0 = self.RLDist_boxcar(z0, field_list=field_list, dz=dz)
      SD0 = dist0.sum(axis=0)
      SD0 = SD0 * normfactor * self.phistar
      axSD.semilogx(SD0, self.ycenters, color=color, ls=':', lw=2, 
                    nonposx='mask',label='boxcar P(z)')
      # Plot a RL distribution with just P(z) but no GALFIT correction
      dist1 = self.RLDist_dkgrid(field_list=field_list)
      SD1 = dist1.sum(axis=0)
      SD1 = SD1 * normfactor * self.phistar
      axSD.semilogx(SD1, self.ycenters, color=color, ls='--', lw=2,
                    nonposx='mask',label='w/ dropout sel. kernel')
      # show size-distribution related parameters
      logr0, sigma, beta = self.parameters[2:]
      # parstr = r'$\log R_0=%.2f$ [arcsec]' % logr0
      # parstr = parstr + '\n' + r'$\sigma=%.2f$' % sigma
      # parstr = parstr + '\n' + r'$\beta=%.2f$' % beta
      # props = dict(boxstyle='round', color='wheat', alpha=0.5)
      # axSD.text(0.05, 0.95, parstr, transform=axSD.transAxes, va='top', 
      #           ha='left', size=textsize, bbox=props)
      axSD.set_ylim(self.axCent.get_ylim())
      axSD.set_yticklabels([])
      axSD.set_xticks([1., 1.e1, 1.e2])   
      axSD.set_xlim(1.e-1, max(n)*50.)
      axSD.set_xlabel(r'$dN/d\log R$', size=self.axlabelsize)
      plt.draw()
      return axSD

   def show_SD(self, z0, dz=1.0, field_list=None, axSD=None, sdbw=0.2, color='black', recalc=True, textsize=14, addtext=True):
      """
      Show the size distribution in a separate plot.
      Unfortunately, I can't reuse self.append_SD...
      """
      if axSD == None:
         fig = plt.figure()
         axSD = fig.add_subplot(111)
      self.axSD = axSD
      all_logre = np.zeros(0)
      if recalc:
         self.loglikelihood_tot(self.parameters, floor=1.e-50)
      for f in field_list:
         all_logre = np.concatenate([all_logre, self.galaxy_samples[f].logr])
      M0 = self.RLDist_factories[field_list[0]].RLDist
      nbins = int(round((self.limits[1,1] - self.limits[1,0]) / sdbw))
      n , bins = np.histogram(all_logre, 
                 np.linspace(self.limits[1,0], self.limits[1,1], nbins+1))
      nerr = [np.sqrt(n), np.sqrt(n)]
      # set a lower limit to the Poisson error bar if n == 1
      for i in range(len(n)):
         if n[i] == 1: nerr[0][i] = 1.-1.e-3
      logrbincenters = (bins[:-1] - self.limits[1,0] + sdbw / 2.) / self.dy
      axSD.errorbar(logrbincenters, n, yerr=nerr, fmt='.', ms=14, mfc=color, 
                    ls='None', mec=color, ecolor=color, capsize=6)
      SD = self.dist_array.sum(axis=0)  # fully-corrected size distribution
      normfactor = len(all_logre) * sdbw / (np.sum(SD) * self.dy)
      SD = SD * normfactor
      axSD.semilogy(self.ycenters, SD, color=color, nonposy='mask', 
                    label='GALFIT TF', lw=3.)
      # Plot a RL distribution with boxcar P(z)
      dist0 = self.RLDist_boxcar(z0, field_list=field_list, dz=dz)
      SD0 = dist0.sum(axis=0)
      SD0 = SD0 * normfactor * self.phistar
      axSD.semilogy(self.ycenters, SD0, color=color, ls=':', lw=2, 
                    nonposy='mask',label='boxcar P(z)')
      # Plot a RL distribution with just P(z) but no GALFIT correction
      dist1 = self.RLDist_dkgrid(field_list=field_list)
      SD1 = dist1.sum(axis=0)
      SD1 = SD1 * normfactor * self.phistar
      axSD.semilogy(self.ycenters, SD1, color=color, ls='--', lw=2,
                    nonposy='mask',label='w/ dropout sel. kernel')
      if addtext:
         # show size-distribution related parameters
         logr0, sigma, beta = self.parameters[2:]
         parstr = r'$\log R_0=%.2f$ [arcsec]' % logr0
         parstr = parstr + '\n' + r'$\sigma=%.2f$' % sigma
         parstr = parstr + '\n' + r'$\beta=%.2f$' % beta
         props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
         axSD.text(0.05, 0.95, parstr, transform=axSD.transAxes, va='top', 
                   ha='left', size=textsize, bbox=props)
      axSD.set_xlim(0, M0.value.shape[1])
      axSD.set_xticks(self.axCent.get_xticks())
      xticklabels = self.axCent.get_yticklabels()
      xticklabels = map(lambda x: x.get_text(), xticklabels)
      axSD.set_xticklabels(xticklabels)
      axSD.set_xlabel('log10(Re) [arcsec]', size=self.axlabelsize)
      axSD.set_yticks([1., 1.e1, 1.e2])   
      axSD.set_ylabel(r'$dN/d\log R$', size=self.axlabelsize)
      axSD.set_ylim(1.e-1, max(n)*50.)
      plt.draw()
      return axSD

class GOODS_PlotBdropsFit(PlotRLDistFit):
   def __init__(self, pars, paramfile, **kwargs):
      super(GOODS_PlotBdropsFit, self).__init__(pars, paramfile, 'f775w', 
                                                'B-dropouts', **kwargs)
   # to show fitting results for B-dropouts and V-dropouts in GOODS + HUDF
   def show_LF_all(self, z0, dz=1.0, lfbw=0.2, colors=['red','blue'], addtext=False):
      fig = plt.figure()
      axes = []
      fig2 = plt.figure()
      newax = fig2.add_subplot(111)
      for i in range(2):
         axes += [fig.add_subplot(2, 1, i+1)]
      for i in range(len(self.fields)):
         f = self.fields[i]
         print "Showing %s..." % f.upper()
         newax = self.show_RLDist(field_list=[f], ax=newax)
         ax = axes[i]
         self.show_LF(z0, dz=dz, field_list=[f], lfbw=lfbw, axLF=ax, 
                      color=colors[i], recalc=True, addtext=addtext)
         ax.text(0.05, 0.05, f.upper(), size=18, va='bottom', ha='left',
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
                 transform=ax.transAxes)
         self.set_LF_labels(f, ax)
      return axes

   def show_SD_all(self, z0, dz=1.0, sdbw=0.2, colors=['red','blue'], addtext=False):
      fig = plt.figure()
      axes = []
      fig2 = plt.figure()
      newax = fig2.add_subplot(111)
      for i in range(2):
         axes += [fig.add_subplot(2, 1, i+1)]
      for i in range(len(self.fields)):
         f = self.fields[i]
         print "Showing %s..." % f.upper()
         newax = self.show_RLDist(field_list=[f], ax=newax)
         ax = axes[i]
         self.show_SD(z0, dz=dz, field_list=[f], sdbw=sdbw, axSD=ax, 
                      color=colors[i], recalc=True, addtext=addtext)
         ax.text(0.95, 0.95, f.upper(), size=18, va='top', ha='right',
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
                 transform=ax.transAxes)
         self.set_SD_labels(f, ax, direction='x')
         M = self.RLDist_factories[f].RLDist
         ax.set_xlim(0, M.value.shape[1])
      return axes

class GDS_PlotLBGFit(PlotRLDistFit):
   def __init__(self, pars, paramfile, gfband, dropname, **kwargs):
      super(GDS_PlotLBGFit, self).__init__(pars, paramfile, gfband, dropname, 
                                           **kwargs)
   # to show fitting results in GOODS-S
   def show_LF_all(self, z0, dz=1.0, lfbw=0.2, colors=['purple','blue','green','red'], addtext=False):
      fig = plt.figure()
      axes = []
      fig2 = plt.figure()
      newax = fig2.add_subplot(111)
      for i in range(4):
         axes += [fig.add_subplot(2, 2, i+1)]
      for i in range(len(self.fields)):
         f = self.fields[i]
         print "Showing %s..." % f.upper()
         newax = self.show_RLDist(field_list=[f], ax=newax)
         ax = axes[i]
         self.show_LF(z0, dz=dz, field_list=[f], lfbw=lfbw, axLF=ax, 
                      color=colors[i], recalc=True, addtext=addtext)
         ax.text(0.05, 0.95, f.upper(), size=18, va='top', ha='left',
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
                 transform=ax.transAxes)
         self.set_LF_labels(f, ax)
      return axes

   def show_SD_all(self, z0, dz=1.0, sdbw=0.2, colors=['purple','blue','green','red'], addtext=False):
      fig = plt.figure()
      axes = []
      fig2 = plt.figure()
      newax = fig2.add_subplot(111)
      for i in range(4):
         axes += [fig.add_subplot(2, 2, i+1)]
      for i in range(len(self.fields)):
         f = self.fields[i]
         print "Showing %s..." % f.upper()
         newax = self.show_RLDist(field_list=[f], ax=newax)
         ax = axes[i]
         self.show_SD(z0, dz=dz, field_list=[f], sdbw=sdbw, axSD=ax, 
                      color=colors[i], recalc=True, addtext=addtext)
         ax.text(0.05, 0.95, f.upper(), size=18, va='top', ha='left',
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
                 transform=ax.transAxes)
         self.set_SD_labels(f, ax, direction='x')
         M = self.RLDist_factories[f].RLDist
         ax.set_xlim(0, M.value.shape[1])
      return axes

class GDS_PlotUdropstFit(GDS_PlotLBGFit):
   def __init__(self, pars, paramfile, **kwargs):
      super(GDS_PlotUdropstFit, self).__init__(pars, paramfile, 'f606w', 
                                               'U-dropout', **kwargs)

class GDS_PlotidropsFit(GDS_PlotLBGFit):
   def __init__(self, pars, paramfile, **kwargs):
      super(GDS_PlotidropsFit, self).__init__(pars, paramfile, 'f125w', 
                                               'i-dropout', **kwargs)

