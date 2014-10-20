#!/usr/bin/env python

from numpy import *
import matplotlib as mpl
import matplotlib.pyplot as plt
from mypytools import cosmoclass
import scipy
from scipy import optimize
from my_mplfonts import Helvetica
import pmscolors

pms = pmscolors.pmscolors()

mpl.rcParams['lines.markersize'] = 8
cc0 = cosmoclass.cosmoclass(H0=70.2,omega_m=0.275,omega_l=0.725) # used in our work
cc1 = cosmoclass.cosmoclass(H0=70.,omega_m=0.3,omega_l=0.7) 
# used in Oesch+10, Ferguson+04, Bouwens+04
cc2 = cosmoclass.cosmoclass(H0=73.,omega_m=0.24,omega_l=0.76)
# used in Hathi+08 

cc_ours = cc1  # adopt cc1 as the cosmological parameters
zarr = [3.2, 4.0, 5.0, 6.0]
# r0arr_kpc = [2.38, 1.34, 1.19, 1.62]
r0arr_kpc = [2.39, 1.34, 1.19, 1.37]   # updated: 14/04/05
#logr0arr_pix = array([0.809, 0.799])
#r0arr_as = (10.**logr0arr_pix)*0.03
#r0arr_kpc = map(lambda z,r0: cc_ours.adsize(z,r0,unit='kpc'),zarr,r0arr_as)
# SIZE ERRORS AT Z~3.4 AND 6 NEED TO BE REVISED!
sigr0arr_kpc = [[0.2, 0.08, 0.15, 0.4], [0.2, 0.092, 0.28, 0.4]]  
# -row1 +row2

O10_size = {'z':array([5.0, 6.0, 6.8])+0.1,
            'r_mean':array([1.28, 0.9, 0.95]),
            'cc':cc1,
            'r_lo':array([1.05, 0.7, 0.7]),
            'r_hi':array([1.55, 1.1, 1.2]),
            'r_unit':'kpc'}

#O10_size_cc = {'z':O10_size['z'],
#               'r_mean':array([1.31, 0.92, 0.98]),
#               'r_lo':array([1.08, 0.72, 0.72]),
#               'r_hi':array([1.59, 1.13, 1.23])}
O10_size_cc = O10_size

F04_size = {'z':array([1.4, 2.3, 3.0, 4.0, 5.0])-0.1,
            'r_mean':[0.68, 0.375, 0.29, 0.25, 0.26],
            'cc':cc1,
            'r_lo':[0.6, 0.35, 0.275, 0.24, 0.25],
            'r_hi':[0.74, 0.4, 0.3, 0.255, 0.27],
            'r_unit':'arcsec'}

F04_size_cc = {'z':F04_size['z'],
               'r_mean':array(map(lambda z,r:cc_ours.adsize(z,r,unit='kpc'),F04_size['z'],F04_size['r_mean'])),
               'r_lo':array(map(lambda z,r:cc_ours.adsize(z,r,unit='kpc'),F04_size['z'],F04_size['r_lo'])),
               'r_hi':array(map(lambda z,r:cc_ours.adsize(z,r,unit='kpc'),F04_size['z'],F04_size['r_hi']))}
               
B04_size = {'z':array([2.5, 3.8, 4.9, 6.0])+0.2,
            'r_mean':array([1.99, 1.6, 1.18, 0.83]),
            'cc':cc1,
            'r_lo':array([1.85, 1.42, 1.1, 0.7]),
            'r_hi':array([2.13, 1.78, 1.26, 0.96]),
            'unit':'kpc'}

#B04_size_cc = {'z':B04_size['z'],
#               'r_mean':array([2.03, 1.64, 1.21, 0.85]),
#               'r_lo':array([1.89, 1.45, 1.13, 0.72]),
#               'r_hi':array([2.18, 1.82, 1.29, 0.99])}
B04_size_cc = B04_size
               
H08_size = {'z':array([3., 4., 5., 6.])-0.2,
            'r_mean':[0.32, 0.21, 0.14, 0.15],
            'cc':cc2,
            'r_lo':[0.26, 0.18, 0.13, 0.13],
            'r_hi':[0.38, 0.24, 0.15, 0.16],
            'unit':'arcsec'}

H08_size_cc = {'z':H08_size['z'],
               'r_mean':array(map(lambda z,r:cc_ours.adsize(z,r,unit='kpc'),H08_size['z'],H08_size['r_mean'])),
               'r_lo':array(map(lambda z,r:cc_ours.adsize(z,r,unit='kpc'),H08_size['z'],H08_size['r_lo'])),
               'r_hi':array(map(lambda z,r:cc_ours.adsize(z,r,unit='kpc'),H08_size['z'],H08_size['r_hi']))}
               #'r_mean':array([2.52, 1.50, 0.90, 0.88]),
               #'r_lo':array([2.05, 1.28, 0.84, 0.76]),
               #'r_hi':array([2.99, 1.71, 0.97, 0.94])}
               
S03_size = {'z':array([0.1]),
            'r_mean':array([4.]),
            'cc':None,
            'r_lo':array([2.96]),
            'r_hi':array([5.40]),
            'unit':'kpc'}
            

def convert_size_cc(their_cc, r_lit, z_lit, our_cc=cc0, unit='kpc'):
   # convert the size from other papers to the same cosmology as used in our paper
   if unit == 'kpc':
      # first solve for angular size, then calculate the size in kpc in our cosmology
      def dtheta(theta, r_lit, z, their_cc):
         r_try = their_cc.adsize(z, theta, unit='kpc')
         return abs(r_try - r_lit)
      theta_lit = optimize.fmin(dtheta, 0.1, args=[r_lit, z_lit, their_cc], disp=False)[0]
      r_cc = our_cc.adsize(z_lit, theta_lit, unit='kpc')
      return r_cc
   elif unit == 'arcsec':
      # convert from arcsec to size using our cosmology directly
      r_cc = our_cc.adsize(z_lit, r_lit, unit='kpc')
      return r_cc

def plot_sizeevol(zarr=zarr, r0_kpc=r0arr_kpc, sigr0_kpc=sigr0arr_kpc,
                  plot_S03=True, plot_UVgal=True):
   # plot the size evolution compiled from literature
   fig = plt.figure()
   ax = fig.add_subplot(111)
   ax.errorbar(zarr, r0_kpc, yerr=sigr0_kpc, fmt='*', mfc='blue', 
               ecolor='blue', mec='blue', label='This work', ms=20)
   # Plot Oesch+10
   ax.errorbar(O10_size_cc['z'], O10_size_cc['r_mean'],
      yerr=[O10_size_cc['r_mean']-O10_size_cc['r_lo'],
      O10_size_cc['r_hi']-O10_size_cc['r_mean']], fmt='^', color='green', 
      ecolor='green', mec='green', label='Oe10')
   # Plot Ferguson+04
   ax.errorbar(F04_size['z'], F04_size_cc['r_mean'],
      yerr=[F04_size_cc['r_mean']-F04_size_cc['r_lo'],
      F04_size_cc['r_hi']-F04_size_cc['r_mean']], fmt='s', color='red', 
      ecolor='red', mec='red', label='F04')
   # Plot Bouwens+04
   ax.errorbar(B04_size['z'], B04_size_cc['r_mean'],
      yerr=[B04_size_cc['r_mean']-B04_size_cc['r_lo'],
      B04_size_cc['r_hi']-B04_size_cc['r_mean']], fmt='D', color='purple',
      ecolor='purple', mec='purple', label='B04')
   # Plot Hathi+08
   ax.errorbar(H08_size['z'], H08_size_cc['r_mean'],
      yerr=[H08_size_cc['r_mean']-H08_size_cc['r_lo'],
      H08_size_cc['r_hi']-H08_size_cc['r_mean']], fmt='H', color='#104BA9', 
      ecolor='#104BA9', mec='#104BA9', label='H08')
   if plot_S03:
      # Plot Shen+03
      ax.errorbar(S03_size['z'], S03_size['r_mean'],
         yerr=[S03_size['r_mean']-S03_size['r_lo'],
         S03_size['r_hi']-S03_size['r_mean']], fmt='h', color='#FF7100', 
         ecolor='#FF7100', mec='#FF7100', label='S03 (optical)', ms=10)
   if plot_UVgal:
      ax.plot(0.1, 18.79, '^', color=pms.Blue_Purples, ms=12,
              label='local UV galaxies')
   # Plot size evolution scaling relation
   zarr_evol = arange(0.0, 7.1, 0.1)
   Hz = cc0.Hz(zarr_evol)
   size_ev1 = Hz ** (-1.)
   size_ev1 = size_ev1 * r0arr_kpc[2] / size_ev1[zarr_evol==5.0]
   ax.plot(zarr_evol, size_ev1, c='black', label=r'$H(z)^{-1}$', lw=2.0)
   size_ev2 = Hz ** (-2./3.)
   size_ev2 = size_ev2 * r0arr_kpc[2] / size_ev2[zarr_evol==5.0]
   ax.plot(zarr_evol, size_ev2, '--', c='black', label=r'$H(z)^{-2/3}$', lw=2.0)
   ax.legend(loc=3, numpoints=1)
   ax.set_yscale('log')
   xlims = [0, 7]
   ylims = [0.5, 8]
   ax.set_xlim(xlims)
   ax.set_ylim(ylims)
   xticks = arange(0, 8)
   yticks = arange(1, ylims[1]+1)
   ax.set_xticks(xticks)
   ax.set_xticklabels(xticks)
   ax.set_yticks(yticks)
   ax.set_yticklabels(yticks)
   ax.set_xlabel('Redshift')
   ax.set_ylabel('Size [kpc]')
   if plot_UVgal==True:
      ax.set_ylim(ymax=20.)
      ax.legend(loc=1, numpoints=1)
   
   return ax