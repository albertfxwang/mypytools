#!/usr/bin/env python

import numpy as np
from pygoods import Ftable
import matplotlib as mpl
import matplotlib.pyplot as plt
from my_mplfonts import Helvetica
import pmscolors

pms = pmscolors.pmscolors()


# Plot a simple size-luminosity distribution separated by Sersic indices

# The separation between late type (n < n_sep) and early-type (n > n_sep)
n_sep = 2.5

catdir = '/Users/khuang/Dropbox/Research/bivariate/bivariate_fit/catalogs'

def plot_RLdist(mag_array, re_array, n_array, dropname, colors=['blue','red'], 
                markers=['^','o'], ax=None, ms=10**2, gfband='F606W'):
   """
   Plot the raw RL distribution of the measured magnitudes & effective radii,
   but separated by Sersic indices. This is to see if there is any hint of
   separation in morphology that leads to the broad size distribution.
   """
   if ax==None:
      fig = plt.figure()
      ax = fig.add_subplot(111)
   logre_array = np.log10(re_array)
   disk = (n_array < n_sep)
   ellip = (n_array >= n_sep)
   # Now plot RL distribution for disks (n < n_sep)
   ax.set_yscale('log')
   ax.scatter(mag_array[disk], re_array[disk], s=ms, marker=markers[0],
              c=colors[0], edgecolor='none', label='$n < %.1f$'%n_sep)
   ax.scatter(mag_array[ellip], re_array[ellip], s=ms, marker=markers[1],
              c=colors[1], edgecolor='none', label='$n \geq %.1f$'%n_sep)
   ax.set_xlabel('%s magnitude' % gfband.upper(), 
                 font_properties=Helvetica(20))
   ax.set_ylabel('%s Re [pixels]' % gfband.upper(),
                 font_properties=Helvetica(20))
   plt.xticks(font_properties=Helvetica(14))
   plt.yticks(font_properties=Helvetica(14))
   #ax.set_title(dropname, font_properties=Helvetica(16))
   ax.text(0.05, 0.05, dropname, font_properties=Helvetica(24),
           ha='left', va='bottom', transform=ax.transAxes)
   ax.legend(loc=1, scatterpoints=1, prop=Helvetica(16))

def plot_RLdist_udrops(cat=catdir+'/udrops/udrops_goodss_ubvy_130517_vflags_galfit.fits',
                       markers=['^','o'], ax=None, ms=8**2,
                       colors=[pms.Bright_Blue,pms.Bright_Red]):
   c = Ftable(cat)
   mag_array1 = c.f606w_magout_gf[(c.udf_fit==True)&(c.f606w_magout_gf<=27.5)]
   re_array1 = c.f606w_reout_gf[(c.udf_fit==True)&(c.f606w_magout_gf<=27.5)]
   n_array1 = c.f606w_nout_gf[(c.udf_fit==True)&(c.f606w_magout_gf<=27.5)]

   mag_array2 = c.f606w_magout_gf[(c.deep_fit==True)&(c.f606w_magout_gf<=26.5)]
   re_array2 = c.f606w_reout_gf[(c.deep_fit==True)&(c.f606w_magout_gf<=26.5)]
   n_array2 = c.f606w_nout_gf[(c.deep_fit==True)&(c.f606w_magout_gf<=26.5)]

   mag_array3 = c.f606w_magout_gf[(c.ers_fit==True)&(c.f606w_magout_gf<=26.5)]
   re_array3 = c.f606w_reout_gf[(c.ers_fit==True)&(c.f606w_magout_gf<=26.5)]
   n_array3 = c.f606w_nout_gf[(c.ers_fit==True)&(c.f606w_magout_gf<=26.5)]

   mag_array4 = c.f606w_magout_gf[(c.wide_fit==True)&(c.f606w_magout_gf<=26.5)]
   re_array4 = c.f606w_reout_gf[(c.wide_fit==True)&(c.f606w_magout_gf<=26.5)]
   n_array4 = c.f606w_nout_gf[(c.wide_fit==True)&(c.f606w_magout_gf<=26.5)]

   mag_array = np.concatenate([mag_array1,mag_array2,mag_array3,mag_array4])
   re_array = np.concatenate([re_array1,re_array2,re_array3,re_array4])
   n_array = np.concatenate([n_array1,n_array2,n_array3,n_array4])

   plot_RLdist(mag_array, re_array, n_array, 'U-dropouts', colors=colors,
               markers=markers, ax=ax, ms=ms, gfband='F606W')

def plot_RLdist_bdrops(cat1=catdir+'/bdrops/bdrops_gf_v3.fits',
                       cat2=catdir+'/bdrops/bdrops_udf_gf_v3.fits',
                       markers=['^','o'], ax=None, ms=8**2,
                       colors=[pms.Bright_Blue,pms.Bright_Red]):
   c1 = Ftable(cat1)
   c2 = Ftable(cat2)
   gfcrit1 = (c1.f775w_gfflag==True)&(c1.magout<=26.5)
   gfcrit2 = (c2.f775w_gfflag==True)&(c2.magout<=28.5)
   mag_array1 = c1.magout[gfcrit1]
   re_array1 = c1.reout[gfcrit1]
   n_array1 = c1.nout[gfcrit1]
   mag_array2 = c2.magout[gfcrit2]
   re_array2 = c2.reout[gfcrit2]
   n_array2 = c2.nout[gfcrit2]

   mag_array = np.concatenate([mag_array1, mag_array2])
   re_array = np.concatenate([re_array1, re_array2])
   n_array = np.concatenate([n_array1, n_array2])
   plot_RLdist(mag_array, re_array, n_array, 'B-dropouts', colors=colors,
               markers=markers, ax=ax, ms=ms, gfband='F775W')

def plot_RLdist_vdrops(cat1=catdir+'/vdrops/vdrops_gf_v3.fits',
                       cat2=catdir+'/vdrops/vdrops_udf_gf_v3.fits',
                       markers=['^','o'], ax=None, ms=8**2,
                       colors=[pms.Bright_Blue,pms.Bright_Red]):
   c1 = Ftable(cat1)
   c2 = Ftable(cat2)
   gfcrit1 = (c1.f850lp_gfflag==True)&(c1.magout<=26.5)
   gfcrit2 = (c2.f850lp_gfflag==True)&(c2.magout<=28.5)
   mag_array1 = c1.magout[gfcrit1]
   re_array1 = c1.reout[gfcrit1]
   n_array1 = c1.nout[gfcrit1]
   mag_array2 = c2.magout[gfcrit2]
   re_array2 = c2.reout[gfcrit2]
   n_array2 = c2.nout[gfcrit2]

   mag_array = np.concatenate([mag_array1, mag_array2])
   re_array = np.concatenate([re_array1, re_array2])
   n_array = np.concatenate([n_array1, n_array2])
   plot_RLdist(mag_array, re_array, n_array, 'V-dropouts', colors=colors,
               markers=markers, ax=ax, ms=ms, gfband='F850LP')

def plot_RLdist_idrops(cat=catdir+'/idrops/idrops_goodss_130623_vflags_galfit.fits',
                       markers=['^','o'], ax=None, ms=8**2,
                       colors=[pms.Bright_Blue,pms.Bright_Red]):
   c = Ftable(cat)
   mag_array1 = c.f105w_magout_gf[(c.udf_fit==True)&(c.f105w_magout_gf<=27.5)]
   re_array1 = c.f105w_reout_gf[(c.udf_fit==True)&(c.f105w_magout_gf<=27.5)]
   n_array1 = c.f105w_nout_gf[(c.udf_fit==True)&(c.f105w_magout_gf<=27.5)]

   mag_array2 = c.f105w_magout_gf[(c.deep_fit==True)&(c.f105w_magout_gf<=26.5)]
   re_array2 = c.f105w_reout_gf[(c.deep_fit==True)&(c.f105w_magout_gf<=26.5)]
   n_array2 = c.f105w_nout_gf[(c.deep_fit==True)&(c.f105w_magout_gf<=26.5)]

   mag_array3 = c.f098m_magout_gf[(c.ers_fit==True)&(c.f098m_magout_gf<=26.5)]
   re_array3 = c.f098m_reout_gf[(c.ers_fit==True)&(c.f098m_magout_gf<=26.5)]
   n_array3 = c.f098m_nout_gf[(c.ers_fit==True)&(c.f098m_magout_gf<=26.5)]

   mag_array4 = c.f105w_magout_gf[(c.wide_fit==True)&(c.f105w_magout_gf<=26.5)]
   re_array4 = c.f105w_reout_gf[(c.wide_fit==True)&(c.f105w_magout_gf<=26.5)]
   n_array4 = c.f105w_nout_gf[(c.wide_fit==True)&(c.f105w_magout_gf<=26.5)]

   mag_array = np.concatenate([mag_array1,mag_array2,mag_array3,mag_array4])
   re_array = np.concatenate([re_array1,re_array2,re_array3,re_array4])
   n_array = np.concatenate([n_array1,n_array2,n_array3,n_array4])

   plot_RLdist(mag_array, re_array, n_array, 'i-dropouts', colors=colors,
               markers=markers, ax=ax, ms=ms, gfband='F105W/F098M')

def plot_RLdist_4drops():
   fig = plt.figure(figsize=(14,9))
   ax1 = plt.subplot(2,2,1)
   plot_RLdist_udrops(ax=ax1)
   ax2 = plt.subplot(2,2,2)
   plot_RLdist_bdrops(ax=ax2)
   ax3 = plt.subplot(2,2,3)
   plot_RLdist_vdrops(ax=ax3)
   ax4 = plt.subplot(2,2,4)
   plot_RLdist_idrops(ax=ax4)
   ax4.set_xlim(24., 28.)
   for ax in fig.axes:
      ax.set_ylim(0.1, 200.0)
   plt.subplots_adjust(top=0.95, bottom=0.07)
