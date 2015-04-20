#!/usr/bin/env python

import numpy as np
from sedtools import galtemp_colors as gc
import pysynphot as S
import starlib
import matplotlib.pyplot as plt
from dropout_selection import lbg_colorcrit as lcc
import os

biv2dir = '/Users/khuang/Dropbox/codes/mypytools/bivariate2'
lbgtemp = '%s/lbgtemp_z6_0.4Zsolar_constSFH_100Myr.sed' % biv2dir
magfiles_lowz = {'SB2': 'SB2colors_ext0.0_HST.dat',
                 'SB3': 'SB3colors_ext0.0_HST.dat',
                 'ES0': 'ES0colors_ext0.0_HST.dat',
                 'Sbc': 'Sbccolors_ext0.0_HST.dat',
                 'Scd': 'Scdcolors_ext0.0_HST.dat'}
starcolors_default = '/Users/khuang/Dropbox/codes/mypytools/pickles_stars_colors.cat'
sedtools = '/Users/khuang/Dropbox/codes/mypytools/sedtools'

# there are some observational evidence (see Castellanos et al. 2014) that 
# LBGs at z~3 have ~0.4 Z_solar and age of ~100-200 Myr

def LBG_colors(magfile, zcenter, sedfile=lbgtemp, zdr=1.5, dz=0.1, ebmv=0., extlaw='xgal', filters=gc.default_filters):
   z0 = zcenter - zdr
   z1 = zcenter + zdr
   gc.GalaxyTemplateColors(sedfile, magfile, z0=z0, z1=z1, dz=dz, ebmv=ebmv,
                           extlaw=extlaw, filters=filters)

class LBGColorPlotFactory(object):
   def __init__(self, magfiles_lbg, labels, zcenter, zdr=1.5, sedfile=lbgtemp, zmarks=None, starcolors=starcolors_default):
      # Plot color tracks
      # color1 (and same with color2) should be like ['f435w', 'f814w'] (two 
      # filter names in strings); this means the color is f435w - f814w
      self.plots = []
      self.magfiles = magfiles_lbg
      self.labels = labels
      self.zcenter = zcenter
      self.z0 = zcenter - zdr
      self.z1 = zcenter + zdr
      self.sedfile = sedfile
      self.starcolors = starcolors
      for i in range(len(self.magfiles)):
         P = gc.PlotGalaxyColors(self.magfiles[i])
         self.plots += [P]
      self.plots_lowz = {}
      
   def plot_lowz(self, color1, color2, ax=None, label=True, **plt_kwargs):
      if ax==None:
         fig = plt.figure(figsize=(10, 8))
         ax = fig.add_subplot(111)
      for g in magfiles_lowz.keys():
         P = gc.PlotLowZColors(sedtools+'/'+magfiles_lowz[g])
         self.plots_lowz[g] = P
         if label:
            gal_label = g
         else:
            if g == magfiles_lowz.keys()[0]:
               gal_label = "z<3 galaxies"
            else:
               gal_label = ""
         ax = self.plots_lowz[g].plot_colors(color1, color2, z1=self.z0-1., 
                                             label=gal_label, ax=ax, ls='--',
                                             **plt_kwargs)

      print self.plots_lowz.keys()
      # Now plot stars
      self.starlib = starlib.starlib(self.starcolors)
      c1 = self.starlib.colors(*color1)
      c2 = self.starlib.colors(*color2)
      ax.scatter(c2, c1, marker='*', s=8**2, label='stars', 
                 facecolor='yellow')
      return ax

   def plot_tracks(self, color1, color2, zmarks=None, zannotate=None, offset=(-0.5,0.5), legendloc=4, ax=None, title="", sizeannotate=14, zlabel=True, **plt_kwargs):
      if zmarks == None:
         self.zmarks = np.arange(self.z0, self.z1, 0.5)
      else:
         self.zmarks = zmarks
      if ax==None:
         fig = plt.figure(figsize=(10, 8))
         ax = fig.add_subplot(111)
      for i in range(len(self.plots)):
         ax = self.plots[i].plot_colors(color1, color2, z0=self.z0, 
                                        z1=self.z1, label=self.labels[i], 
                                        ax=ax, **plt_kwargs)
         if i == 0:
            ax = self.plots[i].mark_redshifts(self.zmarks, color1, color2, 
                               ax=ax, facecolor='none', edgecolor='red', 
                               s=10**2, marker='o',lw=2, zlabel=zlabel)
            if len(zannotate) > 0:
               c1 = self.plots[i].get_color(*color1)
               c2 = self.plots[i].get_color(*color2)
               for z in zannotate:
                  j = np.argsort(np.abs(self.plots[i].z - z))[0]
                  x = c2[j]
                  y = c1[j]
                  xtext = x + offset[0]
                  ytext = y + offset[1]
                  ax.annotate('z=%.1f' % z, xy=(x,y), 
                              xytext=(xtext, ytext), size=sizeannotate, 
                              arrowprops=dict(facecolor='black', arrowstyle='->')   )
      ax.legend(loc=legendloc, scatterpoints=1)
      ax.set_ylim(-1.,6.)
      ax.set_xlim(xmax=8.)
      ax.set_title(title)
      return ax

   def plot_colorcrit(self, colorcrit, ax, fc='blue', alpha=0.2, **kwargs):
      xmin, xmax = ax.get_xlim()
      ymin, ymax = ax.get_ylim()
      colorcrit.plotcrit_fill(xmin, ymax, ax=ax, xmax_plot=xmax, 
                              ymin_plot=ymin, fc=fc, alpha=alpha,
                              **kwargs)
      return ax

   def plot_galaxy(self, ax, color1, color1_err, color2, color2_err, **ebar_kwargs):
      # plot individual photometric points
      ylolim = False
      xlolim = False
      if color1_err < 0:
        ylolim = True
        color1_err = [[0.], [1.]]
      if color2_err < 0:
        xlolim = True
        color2_err = [[0.], [1.]]
      ax.errorbar(color2, color1, xerr=np.abs(color2_err), 
                  yerr=np.abs(color1_err), 
                  xlolims=xlolim, uplims=ylolim, **ebar_kwargs)
      ax.legend(scatterpoints=1, loc=4, numpoints=1)
      plt.draw()
      return ax

class UdropsPlotFactory(LBGColorPlotFactory):
   def __init__(self, drop='uvimos', droplabel='UBVY105', bands=['uvimos','f435w','f606w','f105w'], zcenter=3.3):
      udropsdir = '/Users/khuang/Dropbox/Research/bivariate/udrops_sample'
      magfiles = ['udrops_gds_colors_ext0.0.dat',
                  'udrops_gds_colors_ext0.1.dat',
                  'udrops_gds_colors_ext0.2.dat',
                  'udrops_gds_colors_ext0.3.dat']
      for i in range(len(magfiles)):
         magfiles[i] = os.path.join(udropsdir, magfiles[i])
      labels = ['E(B-V)=0','E(B-V)=0.1','E(B-V)=0.2','E(B-V)=0.3']
      self.zcenter = zcenter
      self.zmarks = [2.5, 3., 3.5, 4.]
      if len(bands) == 4:
        self.color1 = [bands[0], bands[1]]
        self.color2 = [bands[2], bands[3]]
      else:
        self.color1 = [bands[0], bands[1]]
        self.color2 = [bands[1], bands[2]]
      super(UdropsPlotFactory, self).__init__(magfiles, labels, self.zcenter, 
                                              zmarks=self.zmarks)
      # can also plot the dropout criteria
      self.lcc = lcc.colorcrit()
      self.lcc = self.lcc(drop, droplabel)      
      self.droplabel = droplabel

   def plot_tracks(self, **kwargs):
      ax = super(UdropsPlotFactory, self).plot_tracks(self.color1, 
                                     self.color2, 
                                     zmarks=self.zmarks, **kwargs)
      return ax

   def plot_all(self, xmin=-2, xmax=5, ymin=-0.5, ymax=6.5, **kwargs):
      ax = self.plot_tracks(**kwargs)
      ax.set_xlim(xmin, xmax)
      ax.set_ylim(ymin, ymax)
      ax = self.plot_colorcrit(self.lcc, ax=ax)
      ax.set_title('U-dropouts @ z~%.1f' % self.zcenter)
      ax.text(0.05, 0.95, self.droplabel, ha='left', va='top', 
              bbox=dict(boxstyle='round', facecolor='wheat'),
              transform=ax.transAxes, size='xx-large')
      return ax

class idropsPlotFactory(LBGColorPlotFactory):
   def __init__(self, drop='f775w', droplabel='Hua13'):
      idropsdir = '/Users/khuang/Dropbox/Research/bivariate/idrops_sample'
      magfiles = ['idrops_gds_colors_ext0.0.dat',
                  'idrops_gds_colors_ext0.1.dat',
                  'idrops_gds_colors_ext0.2.dat',
                  'idrops_gds_colors_ext0.3.dat']
      for i in range(len(magfiles)):
         magfiles[i] = os.path.join(idropsdir, magfiles[i])
      labels = ['E(B-V)=0','E(B-V)=0.1','E(B-V)=0.2','E(B-V)=0.3']
      self.zcenter = 6.0
      self.zmarks = [5.0, 5.5, 6.0, 6.5, 7.0]
      self.color1 = ['f775w', 'f850lp']
      self.color2 = ['f850lp', 'f125w']
      super(idropsPlotFactory, self).__init__(magfiles, labels, self.zcenter, 
                                              zmarks=self.zmarks)
      # can also plot the dropout criteria
      self.lcc = lcc.colorcrit()
      self.lcc = self.lcc(drop, droplabel)

   def plot_tracks(self, **kwargs):
      ax = super(idropsPlotFactory, self).plot_tracks(self.color1, 
                                     self.color2, 
                                     zmarks=self.zmarks, **kwargs)
      return ax

   def plot_all(self, xmin=-2, xmax=5, ymin=-0.5, ymax=6.5, **kwargs):
      ax = self.plot_tracks(**kwargs)
      ax.set_xlim(xmin, xmax)
      ax.set_ylim(ymin, ymax)
      ax = self.plot_colorcrit(self.lcc, ax=ax)
      ax.set_title('i-dropouts @ z~%.1f' % self.zcenter)
      return ax   

class BdropsPlotFactory(LBGColorPlotFactory):
   def __init__(self, drop='f435w', droplabel='Gia04'):
      idropsdir = '/Users/khuang/Dropbox/Research/bivariate/bdrops_sample'
      magfiles = ['bdrops_goods_colors_ext0.0.dat',
                  'bdrops_goods_colors_ext0.1.dat',
                  'bdrops_goods_colors_ext0.2.dat',
                  'bdrops_goods_colors_ext0.3.dat']
      for i in range(len(magfiles)):
         magfiles[i] = os.path.join(idropsdir, magfiles[i])
      labels = ['E(B-V)=0','E(B-V)=0.1','E(B-V)=0.2','E(B-V)=0.3']
      self.zcenter = 4.0
      self.zmarks = [3.0, 3.5, 4.0, 4.5, 5.0]
      self.color1 = ['f435w', 'f606w']
      self.color2 = ['f606w', 'f850lp']
      super(BdropsPlotFactory, self).__init__(magfiles, labels, self.zcenter, 
                                              zmarks=self.zmarks)
      # can also plot the dropout criteria
      self.lcc = lcc.colorcrit()
      self.lcc = self.lcc(drop, droplabel)

   def plot_tracks(self, **kwargs):
      ax = super(BdropsPlotFactory, self).plot_tracks(self.color1, 
                                     self.color2, 
                                     zmarks=self.zmarks, **kwargs)
      return ax

   def plot_all(self, xmin=-2, xmax=5, ymin=-0.5, ymax=6.5, **kwargs):
      ax = self.plot_tracks(**kwargs)
      ax.set_xlim(xmin, xmax)
      ax.set_ylim(ymin, ymax)
      ax = self.plot_colorcrit(self.lcc, ax=ax)
      ax.set_title('B-dropouts @ z~%.1f' % self.zcenter)
      return ax   

class VdropsPlotFactory(LBGColorPlotFactory):
   def __init__(self, drop='f606w', droplabel='Gia04'):
      idropsdir = '/Users/khuang/Dropbox/Research/bivariate/vdrops_sample'
      magfiles = ['vdrops_goods_colors_ext0.0.dat',
                  'vdrops_goods_colors_ext0.1.dat',
                  'vdrops_goods_colors_ext0.2.dat',
                  'vdrops_goods_colors_ext0.3.dat']
      for i in range(len(magfiles)):
         magfiles[i] = os.path.join(idropsdir, magfiles[i])
      labels = ['E(B-V)=0','E(B-V)=0.1','E(B-V)=0.2','E(B-V)=0.3']
      self.zcenter = 5.0
      self.zmarks = [4.0, 4.5, 5.0, 5.5, 6.0]
      self.color1 = ['f606w', 'f775w']
      self.color2 = ['f775w', 'f850lp']
      super(VdropsPlotFactory, self).__init__(magfiles, labels, self.zcenter, 
                                              zmarks=self.zmarks)
      # can also plot the dropout criteria
      self.lcc = lcc.colorcrit()
      self.lcc = self.lcc(drop, droplabel)

   def plot_tracks(self, **kwargs):
      ax = super(VdropsPlotFactory, self).plot_tracks(self.color1, 
                                     self.color2, 
                                     zmarks=self.zmarks, **kwargs)
      return ax

   def plot_all(self, xmin=-1.5, xmax=5.5, ymin=-1.0, ymax=5.0, **kwargs):
      ax = self.plot_tracks(**kwargs)
      ax.set_xlim(xmin, xmax)
      ax.set_ylim(ymin, ymax)
      ax = self.plot_colorcrit(self.lcc, ax=ax)
      ax.set_title('V-dropouts @ z~%.1f' % self.zcenter)
      return ax  