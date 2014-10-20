#!/usr/bin/env python

from numpy import *
import matplotlib
import matplotlib.pyplot as plt
import bivariate_lf as bl


sizemag_rel_bdrops = {'mbins':array([[24.,25.],[25.,26.],[26.,27.],[27.,28.5]]),
                      'logr0':array([0.8676, 0.8395, 0.7656, 0.8139]),
                      'sigma':array([0.7786, 0.8590, 0.5938, 0.7704])}

sizemag_rel_vdrops = {'mbins':array([[24.,25.],[25.,26.],[26.,27.],[27.,28.5]]),
                      'logr0':array([0.7300, 0.9153, 0.8620, 0.8365]),
                      'sigma':array([1.0600, 0.8888, 0.8498, 0.6227])}
# sigma is in lnR; convert to sigma in logR by dividing by ln(10)

param_bdrops = {'alpha':-1.6171, 'mstar':-20.5343, 'beta':0.2044}
param_vdrops = {'alpha':-1.6853, 'mstar':-20.5097, 'beta':0.2659}

def fp(fsize):
   f = matplotlib.font_manager.FontProperties(size=fsize)
   return f

def plot_bdrops_vdrops_sizemag(ax, dcolor='black', lcolors=['brown','green','red','orange']):
   linestyles = ['-','--',':','-.']
   mci = bl.mconvert('M1500_to_i.txt')
   mcz = bl.mconvert('M1500_to_z.txt')
   m0i = -21.0 + mci(4.0)    # the i-band apparent mag at z=4 corresponding to M0=-21 at 1500A
   m0z = -21.0 + mcz(5.0)    # the z-band apparent mag at z=5 corresponding to M0=-21 at 1500A
   logrpeak_b = zeros(4); logrpeak_v = zeros(4)
   mavg_b = zeros(4); mavg_v = zeros(4)
   yerrs_b = zeros((2,4)); yerrs_v = zeros((2,4))
   logsigma_b = sizemag_rel_bdrops['sigma'] / log(10.)
   logsigma_v = sizemag_rel_vdrops['sigma'] / log(10.)
   for i in range(4):
      # B-dropouts
      mavg_b[i] = average(sizemag_rel_bdrops['mbins'][i])
      logrpeak_b[i] = sizemag_rel_bdrops['logr0'][i] - 0.4 * param_bdrops['beta'] * \
         (mavg_b[i] - m0i)
      yerrs_b[0][i] = 10.**logrpeak_b[i] * (1. - 10.**((-1.)*logsigma_b[i]))
      yerrs_b[1][i] = 10.**logrpeak_b[i] * (10.**logsigma_b[i] - 1.)
      # V-dropouts
      mavg_v[i] = average(sizemag_rel_vdrops['mbins'][i])
      logrpeak_v[i] = sizemag_rel_vdrops['logr0'][i] - 0.4 * param_vdrops['beta'] * \
         (mavg_v[i] - m0z)
      yerrs_v[0][i] = 10.**logrpeak_v[i] * (1. - 10.**((-1.)*logsigma_v[i]))
      yerrs_v[1][i] = 10.**logrpeak_v[i] * (10.**logsigma_v[i] - 1.)
   #fig = plt.figure(figsize=(10,8))
   #ax = fig.add_subplot(111)
   ax.errorbar(mavg_b, 0.03 * 10.**logrpeak_b, yerr=0.03*yerrs_b, xerr=[0.5,0.5,0.5,0.75],\
      fmt='o', ms=10, ecolor=dcolor, mec=dcolor, mfc=dcolor, label='B-dropouts')
   ax.errorbar(mavg_v+0.1, 0.03 * 10.**logrpeak_v, yerr=0.03*yerrs_v, xerr=[0.5,0.5,0.5,0.75],\
      fmt='^', ms=10, ecolor=dcolor, mec=dcolor, mfc='none', label='V-dropouts', mew=1.8)
   ax.set_xlabel('apparent magnitude')
   ax.set_ylabel(r'$R_e$ [arcsec]')
   # plot various power-laws
   pindex = [0.1, 0.2, 0.3, 0.4]
   marr = arange(23., 29., 0.1)
   for i in range(len(pindex)):
      Larr = 10.**(marr / -2.5)
      rarr = Larr ** pindex[i]
      f = 0.03*10.**0.8396 / rarr[20]
      rarr = rarr * f
      ax.plot(marr, rarr, ls=linestyles[i], lw=2.0, label=r'$\bar{R} \propto L^{%.1f}$' % pindex[i],color=lcolors[i])
   ax.set_xlim(23.,28.8)
   ax.legend(loc=1, numpoints=1)
   return ax

