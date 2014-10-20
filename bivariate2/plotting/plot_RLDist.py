#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from bivariate2 import RLDist
import pickle_util as pu

class PlotRLDist(RLDist.RLDistributionFactory):
   def __init__(self, maglimits, logrlimits, mag_convert, pixscale=0.03):
      super(PlotRLDist, self).__init__(maglimits, logrlimits, mag_convert, 
                                       pixscale)

   def show_dist(self, ax=None):
      if ax == None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
      shape = self.RLDist.value.shape
      ax.imshow(self.RLDist.value.T, extent=(0,shape[0],0,shape[1]))
      return ax

   def xticks(self, ax, dx=0.5, size=12, every=2):
      maglimits = self.RLDist.xlimits
      nticks = int(round((maglimits[1]-maglimits[0]) / dx))
      xlab = np.linspace(maglimits[0], maglimits[1], num=nticks+1)
      xticks = map(self.RLDist.get_xindex, xlab)
      ax.set_xticks(xticks[::every])
      ax.set_xticklabels(map(lambda x:'%.1f'%x, xlab[::every]), size=size)
      plt.draw()
      return ax

   def yticks(self, ax, dy=0.2, size=12, every=2):
      logrlimits = self.RLDist.ylimits
      nticks = int(round((logrlimits[1]-logrlimits[0]) / dy))
      ylab = np.linspace(logrlimits[0], logrlimits[1], num=nticks+1)
      yticks = map(self.RLDist.get_yindex, ylab)
      ax.set_yticks(yticks[::every])
      ax.set_yticklabels(map(lambda y:'%.2f'%y, 10.**(ylab[::every])), size=size)
      plt.draw()
      return ax

   def param_string(self, par):
      str_alpha = r"$\alpha=%.2f$" % par[0]
      str_mstar = r"$M^*=%.2f$" % par[1]
      str_r0 = r"$R_0=%.2f$" % (10.**par[2])
      str_sigma = r"$\sigma_{\ln R}=%.2f$" % par[3]
      str_beta = r"$\beta=%.2f$" % par[4]
      s = "\n".join([str_alpha,str_mstar,str_r0,str_sigma,str_beta])
      return s
