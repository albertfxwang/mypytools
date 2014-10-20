#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from sedtools import color_tracks as ct
import os, sys
import starlib
import lbg_colorcrit as lcc

class colortracks(ct.color_tracks, lcc.colorcrit):
   def __init__(self, drop, label, zarray=np.arange(0.1,10.,0.1)):
      lcc.colorcrit.__init__(self)
      lcc.colorcrit.__call__(self, drop, label)
      self.b1, self.b2 = self.bands[0], self.bands[1]
      if len(self.bands) == 3:
         self.b3, self.b4 = self.bands[1], self.bands[2]
      elif len(self.bands) == 4:
         self.b3, self.b4 = self.bands[2], self.bands[3]
      # ctname = "lbg_%s_%s_%s_%s_ext%.2f.ct"%(b1,b2,b3,b4,ebmv)
      # if os.path.exists(ctname):
      #    trackdata = ctname
      # else:
      #    trackdata = ""
      ct.color_tracks.__init__(self, zarray=zarray)
      self.fig = plt.figure()
      self.ax = self.fig.add_subplot(111)

   def reset_axes(self):
      self.fig.clf()
      self.ax = self.fig.add_subplot(111)

   def plotcrit_fill(self, xlo, yhi, **kwargs):
      lcc.colorcrit.plotcrit_fill(self, xlo, yhi, ax=self.ax, **kwargs)
      plt.draw()

   def LBG_colorbox(self, xlo, yhi, **kwargs):
      """
      Plot the color-selection region in shades. The criteria takes the 
      following form:
      i - J >= imj_lim  AND
      J - H <= jmh_lim  AND
      i - J >= coeff1 * (J - H) + coeff0
      """
      self.plotcrit_fill(xlo, yhi, ax=self.ax, **kwargs)
      self.ax.draw()
      # if self.ax==None:
      #    fig = plt.figure()
      #    self.ax = fig.add_subplot(111)
      # imj_array = np.arange(-1., 6.1, 0.1)
      # jmh_array = np.arange(-0.5, jmh_lim+0.05, 0.05)
      # def y(x):
      #    return np.maximum(np.ones(len(x))*imj_lim, coeff1 * x + coeff0)
      # self.ax.fill_between(jmh_array, np.ones(len(jmh_array))*imj_array[-1], 
      #                      y(jmh_array),
      #                      color='green', alpha=0.3)
      # self.ax.set_xlabel('F110W - F160W', size=18)
      # self.ax.set_ylabel('F814W - F110W', size=18)
      # self.ax.set_xlim(jmh_array[0],1.)
      # self.ax.set_ylim(-0.5,imj_array[-1])

   def calc_tracks(self, sp, ebmv=0.0, extlaw='xgal', lya_ew=0.0):
      """
      Given an unreddened model SED for LBG, calculate the color tracks.
      """
      ct.color_tracks.calc_tracks(self, sp, self.b1, self.b2, self.b3, self.b4,
                                  ebmv=ebmv, extlaw=extlaw, lya_ew=lya_ew,
                                  spkw='lbg')
      self.sedlabel='E(B-V)=%.2f' % ebmv

   def plot_color_tracks(self, annotate_z=[6.,7.,8.], **pltkw):
      ct.color_tracks.plot_tracks(self, ax=self.ax, setlabels=True, 
                                  zmarks=np.arange(2.,self.zarray.max(), 0.5), 
                                  drlim=[5.0, None],
                                  annotate_z=annotate_z, **pltkw)
   def plot_interlopers(self):
      galtypes = ['SB2', 'SB3', 'eso', 'imm', 'sbc', 'scd']
      for g in galtypes:
         # ctname = '%s_%s_%s_%s_%s_ext%.2f.ct' % (g, self.b1, self.b2, self.b3, self.b4, )
         # if not os.path.exists(ctname):
         ctg = ct.ctrack_local_galaxies(g, zarray=np.arange(0.,3.1,0.1))
         ctg.calc_tracks(self.b1, self.b2, self.b3, self.b4, ebmv=0.)
         # else:
         #    ctg = ct.ctrack_local_galaxies(g, zarray=np.arange(0.,3.1,0.1))
         ctg.sedlabel = g.upper()
         ctg.plot_tracks(ax=self.ax, drlim=[5.0, None])
      # Plot stars
      stars = starlib.starlib('/Users/khuang/Dropbox/codes/mypytools/pickles_stars_colors.cat')
      stars.plot_colorcolor(self.b1, self.b2, self.b3, self.b4, ax=self.ax, 
                            label="Stars")

def plot_zdrops_select(imj_lim, jmh_lim, coeff0, coeff1):
   za0 = zdrops_A370(ebmv=0.)
   za1 = zdrops_A370(ebmv=0.1)
   za0.LBG_colorbox(imj_lim, jmh_lim, coeff0, coeff1)
   za1.ax = za0.ax
   za0.plot_color_tracks()
   za1.plot_color_tracks(annotate_z=[])
   za1.plot_interlopers()


