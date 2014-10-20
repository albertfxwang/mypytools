#!/usr/bin/env python

import numpy as np
import bivRL as bl
from pygoods import Ftable
import matplotlib as mpl
import matplotlib.pyplot as plt
from my_mplfonts import Helvetica
import pmscolors

pms = pmscolors.pmscolors()


# Plot a simple size-luminosity distribution separated by Sersic indices

# The separation between late type (n < n_sep) and early-type (n > n_sep)
n_sep = 2.5
mag_lims = np.array([23.,29.])
logre_lims = np.array([-0.6,2.0])
magpix = np.arange(mag_lims[0], mag_lims[1], 0.02)

catdir = '/Users/khuang/Dropbox/Research/bivariate/bivariate_fit/catalogs'

def plot_RLdist(mag_array, re_array, dropname, logR0, beta, mean_kcorr,
                color='blue', marker='o', ax=None, ms=8**2, gfband='F606W'):
   """
   Plot the raw RL distribution of the measured magnitudes & effective radii,
   and fit a raw RL relation. 
   """
   if ax==None:
      fig = plt.figure()
      ax = fig.add_subplot(111)
   logre_array = np.log10(re_array)
   ax.set_yscale('log')
   ax.scatter(mag_array, re_array, s=ms, marker=marker,
              c=color, edgecolor='none')
   # Now fit a straight line to the points... just use np.polyfit (this assumes
   # there is only error in y, but what the heck)
   p = np.polyfit(mag_array, logre_array, 1)
   print "p[0], p[1]", p[0], p[1]
   # Plot the best-fit raw RL relation
   ax.plot(magpix, 10.**(p[0]*magpix+p[1]), lw=2., color='0.2', ls='--',
           label=r'$\beta_{\mathrm{raw}}=%.2f$'%(p[0]/-0.4))
   # Plot the best-fit beta after corrections
   beta_M = -0.4 * beta
   m0 = -21.0 + mean_kcorr
   logR = logR0 - 0.4 * beta * (magpix-m0)
   ax.plot(magpix, 10.**logR, lw=2., color='black', ls='-',
           label=r'$\beta=%.2f$'%beta) 
   ax.set_xlabel('%s magnitude' % gfband.upper(), 
                 font_properties=Helvetica(20))
   ax.set_ylabel('%s Re [pixels]' % gfband.upper(),
                 font_properties=Helvetica(20))
   xticks = np.arange(mag_lims[0], mag_lims[1])
   yticks = 10.**np.arange(-1., 3., 1.)
   ax.set_xticks(xticks)
   ax.set_xticklabels(map(lambda x:'%d'%x, xticks), 
                      font_properties=Helvetica(14))
   plt.yticks(font_properties=Helvetica(14))
   #ax.set_title(dropname, font_properties=Helvetica(16))
   ax.text(0.05, 0.05, dropname, font_properties=Helvetica(24),
           ha='left', va='bottom', transform=ax.transAxes)
   ax.legend(loc=1, scatterpoints=1, prop=Helvetica(16))
   ax.set_xlim(*mag_lims)
   ax.set_ylim(*(10.**logre_lims))

def bootstrap_error(mag_array, logre_array, nsamp=10000):
  """
  Estimate the uncertainty in the slope of the RL relation by bootstrap 
  sampling.
  """
  # p_array holds the slope for all resampled realizations
  p_array = np.zeros(nsamp)
  nobj = len(mag_array)
  for i in range(nsamp):
    index = np.random.choice(nobj, size=nobj, replace=True)
    mag_samp = mag_array.take(index)
    logre_samp = logre_array.take(index)
    p = np.polyfit(mag_samp, logre_samp, 1)
    p_array[i] = p[0]
  # return beta, not beta_M
  return p_array / (-0.4)


def plot_RLdist_udrops(logR0, beta, 
    cat=catdir+'/udrops/udrops_goodss_ubvy_130517_vflags_galfit.fits',
    marker='^', ax=None, ms=8**2,
    color=pms.Bright_Blue, bootstrap=False, nsamp=10000):
   c = Ftable(cat)
   mcv = bl.mconvert('kcorr/M1500_to_f606w_omega_m_0.3.txt')
   mean_kcorr = mcv(3.4)
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

   plot_RLdist(mag_array, re_array, 'U-dropouts', logR0, beta, mean_kcorr,
               color=color, marker=marker, ax=ax, ms=ms, gfband='F606W')
   if bootstrap == True:
      p_array = bootstrap_error(mag_array, np.log10(re_array), nsamp=nsamp)
      return p_array
   else:
      return 0

def plot_RLdist_bdrops(logR0, beta, 
                       cat1=catdir+'/bdrops/bdrops_gf_v3.fits',
                       cat2=catdir+'/bdrops/bdrops_udf_gf_v3.fits',
                       marker='o', ax=None, ms=8**2,
                       color=pms.Bright_Blue, bootstrap=False, nsamp=10000):
   c1 = Ftable(cat1)
   c2 = Ftable(cat2)
   mci = bl.mconvert('kcorr/M1500_to_f775w_omega_m_0.3.txt')
   mean_kcorr = mci(4.0)
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
   plot_RLdist(mag_array, re_array, 'B-dropouts', logR0, beta, mean_kcorr,
               color=color, marker=marker, ax=ax, ms=ms, gfband='F775W')
   if bootstrap == True:
      p_array = bootstrap_error(mag_array, np.log10(re_array), nsamp=nsamp)
      return p_array
   else:
      return 0

def plot_RLdist_vdrops(logR0, beta, 
                       cat1=catdir+'/vdrops/vdrops_gf_v3.fits',
                       cat2=catdir+'/vdrops/vdrops_udf_gf_v3.fits',
                       marker='o', ax=None, ms=8**2,
                       color=pms.Bright_Blue, bootstrap=False, nsamp=10000):
   c1 = Ftable(cat1)
   c2 = Ftable(cat2)
   mcz = bl.mconvert('kcorr/M1500_to_f850lp_omega_m_0.3.txt')
   mean_kcorr = mcz(5.0)
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
   plot_RLdist(mag_array, re_array, 'V-dropouts', logR0, beta, mean_kcorr,
               color=color, marker=marker, ax=ax, ms=ms, gfband='F850LP')
   if bootstrap == True:
      p_array = bootstrap_error(mag_array, np.log10(re_array), nsamp=nsamp)
      return p_array
   else:
      return 0

def plot_RLdist_idrops(logR0, beta,
    cat=catdir+'/idrops/idrops_goodss_130623_vflags_galfit.fits',
    marker='o', ax=None, ms=8**2,
    color=pms.Bright_Blue, bootstrap=False, nsamp=10000):
   c = Ftable(cat)
   mcy = bl.mconvert('kcorr/M1500_to_f105w_omega_m_0.3.txt')
   mean_kcorr = mcy(6.0)
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

   plot_RLdist(mag_array, re_array, 'i-dropouts', logR0, beta, mean_kcorr,
               color=color, marker=marker, ax=ax, ms=ms, gfband='F105W/F098M')
   if bootstrap == True:
      p_array = bootstrap_error(mag_array, np.log10(re_array), nsamp=nsamp)
      return p_array
   else:
      return 0

def plot_RLdist_4drops():
   fig = plt.figure(figsize=(14,9))
   ax1 = plt.subplot(2,2,1)
   plot_RLdist_udrops(0.729, 0.227, ax=ax1, marker='o', ms=5**2)
   ax2 = plt.subplot(2,2,2)
   plot_RLdist_bdrops(0.809, 0.218, ax=ax2, marker='o', ms=5**2)
   ax3 = plt.subplot(2,2,3)
   plot_RLdist_vdrops(0.799, 0.252, ax=ax3, marker='o', ms=5**2)
   ax4 = plt.subplot(2,2,4)
   plot_RLdist_idrops(0.675, 0.181, ax=ax4, marker='o', ms=5**2)
   #ax4.set_xlim(24., 28.)
   for ax in fig.axes:
      ax.set_ylim(0.1, 200.0)
   plt.subplots_adjust(top=0.95, bottom=0.07, left=0.07, right=0.97)
