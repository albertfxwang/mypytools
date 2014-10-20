#!/usr/bin/env python

from numpy import *
from pygoods import *
from bivariate import bivRL as bl
from bivariate import fit_lbg as fl
import matplotlib
import matplotlib.pyplot as plt
from bivariate import fit_lognormal as fitln
import scipy
from scipy import optimize
from my_mplfonts import Helvetica

# Best-fit parameters (031512)
parb = array([-1.68323375, -20.60007504, 0.80905733, 0.83207639, 0.21824641])
parv = array([-1.73995114, -20.53177958, 0.79928352, 0.90250446, 0.25213149])
# paru = array([-1.88230264, -20.67758783, 0.7270959,  0.72667011, 0.22459084]) # 130715_1
paru = array([-1.8391, -20.37539, 0.691768, 0.78518, 0.02624])
# pari = array([-1.84221235, -20.32331588, 0.6751357,  0.28752996, 0.1814661]) # 130715_1
pari = array([-1.9537, -21.27484, 0.59228, 0.76123, 0.84313])
limits1 = array([[21.0,26.5],[-1.0,2.0]])
limits2 = array([[23.0,28.5],[-1.0,2.0]])

catdir = '/Users/khuang/Dropbox/Research/bivariate/bivariate_fit/catalogs'
bivdir = '/Users/khuang/Dropbox/Research/bivariate/bivariate_fit'
logre_min = -0.5
logre_max = 2.0

def plot_sizedist(re_galfit, re_sex, dropname, sigma_best, colors=['blue', 'green'],
                  histls=['solid','dashed'], ls=['-','--'], ax=None,
                  binwd=0.1, pixdx=0.02, guess_gf=[5.0, 0.5], 
                  guess_se=[5.0, 0.5]):
   if ax==None:
      fig = plt.figure()
      ax = fig.add_subplot(111)
   logre_gf = log10(re_galfit)
   logre_se = log10(re_sex)
   logrebins = arange(logre_min, logre_max, binwd)
   logre_array = arange(logre_min, logre_max, pixdx)   
   hist_gf = ax.hist(logre_gf, logrebins, histtype='step', lw=1.5, 
             facecolor='none',
             edgecolor=colors[0], normed=True, label='GALFIT', ls=histls[0])
   hist_se = ax.hist(logre_se, logrebins, histtype='step', lw=1.5, 
             facecolor='none',
             edgecolor=colors[1], normed=True, label='SExtractor', ls=histls[1])
   p_galfit = fitln.fit_lognormal(re_galfit, *guess_gf)
   p_sex = fitln.fit_lognormal(re_sex, *guess_se)
   print "p_galfit", p_galfit
   print "p_sex", p_sex
   ln_galfit = fitln.lognormal(p_galfit[0], p_galfit[1], 10.**logre_array)
   ln_sex = fitln.lognormal(p_sex[0], p_sex[1], 10.**logre_array)
   ln_lambda = fitln.lognormal(p_galfit[0], 0.5, 10.**logre_array)
   ln_best = fitln.lognormal(p_galfit[0], sigma_best, 10.**logre_array)
   nf_galfit = hist_gf[0].max() / (ln_galfit*(10.**logre_array)).max() 
   nf_sex = hist_se[0].max() / (ln_sex*(10.**logre_array)).max()
   nf_lambda = hist_gf[0].max() / (ln_lambda*(10.**logre_array)).max()
   nf_best = hist_gf[0].max() / (ln_best*(10.**logre_array)).max()
   ax.plot(logre_array, (10.**logre_array)*ln_galfit*nf_galfit, 
           lw=3.0, ls=ls[0], 
           color=colors[0], label='$\sigma=%.2f$' % p_galfit[1])
   ax.plot(logre_array, (10.**logre_array)*ln_sex*nf_sex, 
           lw=3.0, ls=ls[1],
           color=colors[1], label='$\sigma=%.2f$' % p_sex[1])
   ax.plot(logre_array, (10.**logre_array)*ln_lambda*nf_lambda,
           lw=1.5, ls='-', color='black')
   ax.plot(logre_array, (10.**logre_array)*ln_best*nf_best,
           lw=1.5, ls='-', color='red')
   ax.set_xticklabels(ax.get_xticks(), font_properties=Helvetica(14))
   ax.set_yticklabels(ax.get_yticks(), font_properties=Helvetica(14))
   ax.set_xlabel('log10(Re) [pixel]', font_properties=Helvetica(18))
   ax.set_ylabel('P(Re)', font_properties=Helvetica(18))
   #ax.set_title(dropname, font_properties=Helvetica(16))
   ax.text(0.95, 0.5, dropname, transform=ax.transAxes, 
           horizontalalignment='right', verticalalignment='center',
           font_properties=Helvetica(18),
           bbox=dict(boxstyle='round', facecolor='#A5DB92', alpha=0.5))
   ax.legend(loc=9, ncol=2)
   ax.set_ylim(ymax=maximum(p_galfit.max(), p_sex.max())*1.2)


def plot_sizedist_bdrops(sigma_best, cat_goods=bivdir+'/bdrops_fitting/catalogs/bdrops_gf_v3.fits',
                         cat_udf=bivdir+'/bdrops_fitting/catalogs/bdrops_udf_gf_v3.fits',
                         colors=['blue','green'], histls=['solid','dashed'],
                         ls=['-','--'], ax=None, binwd=0.1, pixdx=0.02):
   c1 = Ftable(cat_goods)
   c2 = Ftable(cat_udf)
   mag_gf1 = c1.magout[c1.f775w_gfflag==True]
   re_gf1 = c1.reout[c1.f775w_gfflag==True]
   re_se1 = c1.i_flux_radius_1[c1.f775w_gfflag==True]
   mag_gf2 = c2.magout[c2.f775w_gfflag==True]
   re_gf2 = c2.reout[c2.f775w_gfflag==True]
   re_se2 = c2.i_flux_radius_1[c2.f775w_gfflag==True]
   re_gf_all = concatenate([re_gf1, re_gf2])
   re_se_all = concatenate([re_se1, re_se2])
   plot_sizedist(re_gf_all, re_se_all, 'z ~ 4', sigma_best, colors=colors,
                 histls=histls, ls=ls, ax=ax, binwd=binwd, pixdx=pixdx)

def plot_sizedist_vdrops(sigma_best, cat_goods=bivdir+'/vdrops_fitting/catalogs/vdrops_gf_v3.fits',
                         cat_udf=bivdir+'/vdrops_fitting/catalogs/vdrops_udf_gf_v3.fits',
                         colors=['blue','green'], histls=['solid','dashed'],
                         ls=['-','--'], ax=None, binwd=0.1, pixdx=0.02):
   c1 = Ftable(cat_goods)
   c2 = Ftable(cat_udf)
   mag_gf1 = c1.magout[c1.f850lp_gfflag==True]
   re_gf1 = c1.reout[c1.f850lp_gfflag==True]
   re_se1 = c1.i_flux_radius_1[c1.f850lp_gfflag==True]
   mag_gf2 = c2.magout[c2.f850lp_gfflag==True]
   re_gf2 = c2.reout[c2.f850lp_gfflag==True]
   re_se2 = c2.i_flux_radius_1[c2.f850lp_gfflag==True]
   re_gf_all = concatenate([re_gf1, re_gf2])
   re_se_all = concatenate([re_se1, re_se2])
   plot_sizedist(re_gf_all, re_se_all, 'z ~ 5', sigma_best, colors=colors,
                 histls=histls, ls=ls, ax=ax, binwd=binwd, pixdx=pixdx)

def plot_sizedist_udrops(sigma_best, cat=bivdir+'/udrops_fitting/catalogs/udrops_goodss_ubvy_130517_vflags_galfit.fits',
                         colors=['blue','green'], histls=['solid','dashed'],
                         ls=['-','--'], ax=None, binwd=0.1, pixdx=0.02):
   c = Ftable(cat)
   re_gf = c.f606w_reout_gf[c.f606w_gfflag==True]
   re_se = c.acs_f606w_flux_radius_1[c.f606w_gfflag==True]
   plot_sizedist(re_gf, re_se, 'z ~ 3.2', sigma_best, colors=colors,
                 histls=histls, ls=ls, ax=ax, binwd=binwd, pixdx=pixdx)

def plot_sizedist_idrops(sigma_best, cat=bivdir+'/idrops_fitting/catalogs/idrops_goodss_130623_vflags_galfit.fits',
                         colors=['blue','green'], histls=['solid','dashed'],
                         ls=['-','--'], ax=None, binwd=0.1, pixdx=0.02,
                         guess_gf=[5.0, 0.5], guess_se=[2.0, 0.4]):
   c = Ftable(cat)
   re_gf_f105w = c.f105w_reout_gf[(c.f105w_gfflag==True)&(c.ers==False)&(c.f105w_flux_radius_1>0)]
   re_se_f105w = c.f105w_flux_radius_1[(c.f105w_gfflag==True)&(c.ers==False)&(c.f105w_flux_radius_1>0)]
   print len(re_gf_f105w)
   plot_sizedist(re_gf_f105w, re_se_f105w, 'z ~ 6', sigma_best, colors=colors,
                 histls=histls, ls=ls, ax=ax, binwd=binwd, pixdx=pixdx,
                 guess_gf=guess_gf, guess_se=guess_se)
   #return re_se_f105w

def plot_sizedist_4drops():
   fig = plt.figure(figsize=(10,9))
   ax1 = fig.add_subplot(2,2,1)
   plot_sizedist_udrops(paru[3], ax=ax1)
   ax1.set_ylim(ymax=3.9)
   ax2 = fig.add_subplot(2,2,2)
   plot_sizedist_bdrops(parb[3], ax=ax2)
   ax2.set_ylim(ymax=3.3)
   ax3 = fig.add_subplot(2,2,3)
   plot_sizedist_vdrops(parv[3], ax=ax3)
   ax3.set_ylim(ymax=3.7)
   ax4 = fig.add_subplot(2,2,4)
   plot_sizedist_idrops(pari[3], ax=ax4)
   ax4.set_ylim(ymax=4.9)
   plt.subplots_adjust(bottom=0.07, top=0.95)

def plot_sizedist_bdrops_vdrops(c1, c2, c3, c4, colors=['blue','green','black']):
   # c1, c2 are the B-drops catalogs
   # c3, c4 are the V-drops catalogs
   # plot both the GALFIT and SExtractor sizes, as well as best-fit lognormal
   # functions to the raw distribution
   mod_lograrr = arange(-0.8, 1.8, 0.02)
   mod_lograrr2 = mod_lograrr+0.02
   mod_dr = 10.**mod_lograrr2 - 10.**mod_lograrr
   limits1 = bl.limits1.copy(); limits1[1] = array([-0.8, 1.8])
   limits2 = bl.limits2.copy(); limits2[1] = array([-0.8, 1.8])
   rlim = array([0., 30.])
   dh = 0.5
   fig = plt.figure(figsize=(8,10))
   
   # B-dropouts
   mci = bl.mconvert('M1500_to_i.txt')
   ax1 = fig.add_subplot(211)  # B-drops
   mag1, re1, crit1 = fl.cleandata('bdrops_gf_v2.cat', chisqnulim=0.4,
      limits=limits1, drop='b')
   mag2, re2, crit2 = fl.cleandata('bdrops_udf_gf_v2.cat', chisqnulim=5.0,
      limits=limits2, drop='b')
   r_se12 = concatenate((c1.i_flux_radius_1, c2.i_flux_radius_1))
   r_gf12 = concatenate((re1, re2))
   hse_1 = ax1.hist(r_se12, arange(rlim[0], rlim[1], dh), histtype='step', color=colors[0],
      label='SE size', normed=True,lw=1.5)
   hgf_1 = ax1.hist(r_gf12, arange(rlim[0], rlim[1], dh), histtype='step', color=colors[1],
      label='GF size', normed=True,lw=1.5)
   # plot the intrinsic models
   model01 = bl.bivariate_lf(parb, limits1, bl.pixdx, 'b', 'goods', kgrid=None, 
      zdgrid=None, mc=mci, meankcorr=mci(4.0))
   model02 = bl.bivariate_lf(parb, limits2, bl.pixdx, 'b', 'udf', kgrid=None,
      zdgrid=None, mc=mci, meankcorr=mci(4.0))
   sd01 = model01.model.sum(axis=0); sd01 = sd01 / (sum(sd01))
   sd02 = model02.model.sum(axis=0); sd02 = sd02 / (sum(sd02))
   sd012 = bl.A_GOODS * sd01 + bl.A_UDF * sd02
   sd012 = sd012 / (sum(sd012*0.02))
   sd012 = sd012 / (10.**mod_lograrr * log(10.))
   
   # fit lognormal distribution
   guess = array([5.0, 0.5])  # guess for peak and ln(sigma)
   xout1_se = optimize.fmin(fitln.logl_lognormal, guess, args=[r_se12])
   xout1_gf = optimize.fmin(fitln.logl_lognormal, guess, args=[r_gf12])
   rarr = arange(rlim[0], rlim[1], 0.2)
   ax1.plot(rarr, fitln.lognormal(xout1_se[0], xout1_se[1], rarr), '--', c=colors[0], lw=3.0,
      label=r'$\sigma=%.2f$'%xout1_se[1])
   ax1.plot(rarr, fitln.lognormal(xout1_gf[0], xout1_gf[1], rarr), '--', c=colors[1], lw=3.0,
      label=r'$\sigma=%.2f$'%xout1_gf[1])
   ax1.plot(10.**mod_lograrr, sd012, '-', c=colors[2], lw=2.0, 
      label=r'$\sigma=%.2f$'%(parb[3]))
   
   # V-dropouts
   ax2 = fig.add_subplot(212)
   mcz = bl.mconvert('M1500_to_z.txt')
   mag3, re3, crit3 = fl.cleandata('vdrops_gf_v2.cat', chisqnulim=0.5,
      limits=limits1, drop='v')
   mag4, re4, crit4 = fl.cleandata('vdrops_udf_gf_v2.cat', chisqnulim=5.0,
      limits=limits2, drop='v')
   r_se34 = concatenate((c3.z_flux_radius_1, c4.z_flux_radius_1))
   r_gf34 = concatenate((re3, re4))
   hse_2 = ax2.hist(r_se34, arange(rlim[0], rlim[1], dh), histtype='step', color=colors[0],
      label='SE size', normed=True, lw=1.5)
   hgf_2 = ax2.hist(r_gf34, arange(rlim[0], rlim[1], dh), histtype='step', color=colors[1],
      label='GF size', normed=True, lw=1.5)
   # plot the intrinsic model
   model03 = bl.bivariate_lf(parv, limits1, bl.pixdx, 'v', 'goods', kgrid=None,
      zdgrid=None, mc=mcz, meankcorr=mcz(5.0))
   model04 = bl.bivariate_lf(parv, limits2, bl.pixdx, 'v', 'udf', kgrid=None,
      zdgrid=None, mc=mcz, meankcorr=mcz(5.0))
   sd03 = model03.model.sum(axis=0); sd03 = sd03 / (sum(sd03))
   sd04 = model04.model.sum(axis=0); sd04 = sd04 / (sum(sd04))
   sd034 = bl.A_GOODS * sd03 + bl.A_UDF * sd04
   sd034 = sd034 / (sum(sd034*0.02))
   sd034 = sd034 / (10.**mod_lograrr * log(10.))

   xout2_se = optimize.fmin(fitln.logl_lognormal, guess, args=[r_se34])
   xout2_gf = optimize.fmin(fitln.logl_lognormal, guess, args=[r_gf34])
   rarr = arange(rlim[0], rlim[1], 0.1)
   ax2.plot(rarr, fitln.lognormal(xout2_se[0], xout2_se[1], rarr), '--', c=colors[0], lw=3.0,
      label=r'$\sigma=%.2f$'%xout2_se[1])
   ax2.plot(rarr, fitln.lognormal(xout2_gf[0], xout2_gf[1], rarr), '--', c=colors[1], lw=3.0,
      label=r'$\sigma=%.2f$'%xout2_gf[1])
   ax2.plot(10.**mod_lograrr, sd034, '-', c=colors[2], lw=2.0, 
      label=r'$\sigma=%.2f$'%(parv[3]))
   
   # Plot formatting
   ax1.set_xlabel('Size [pixels]')
   ax1.set_ylabel('P(R)')
   ax1.set_xlim(0, 30)
   ax1.legend(loc=1)
   #ax1.set_title('B-dropouts')
   ax1.text(0.45,0.9,'B-dropouts',size=18,transform=ax1.transAxes)
   ax2.set_xlabel('Size [pixels]')
   ax2.set_ylabel('P(R)')
   ax2.set_xlim(0, 30)
   ax2.legend(loc=1)
   #ax2.set_title('V-dropouts')
   ax2.text(0.45,0.9,'V-dropouts',size=18,transform=ax2.transAxes)
   
   return fig