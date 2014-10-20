#!/usr/bin/env python

from numpy import *
import matplotlib as mpl
import matplotlib.pyplot as plt
from mypytools import cosmoclass
import bivariate_lf as bl
import plot_sizemag as psm

cc = cosmoclass.cosmoclass(H0=70.5, omega_m=0.275, omega_l=0.725)
ad_1as_z4 = cc.adsize(4.0, 1.0, unit='kpc')
ad_1as_z5 = cc.adsize(5.0, 1.0, unit='kpc')

# Plot the size-luminosity relation of Huang+12, Shen+03, and Courteau+07
## Huang et al. 2012
sizemag_rel_bdrops_h12 = {'mbins':array([[24.,25.],[25.,26.],[26.,27.],[27.,28.5]]),
                         'logr0':array([0.8667, 0.8593, 0.7894, 0.8491]),
                         'sigma':array([0.7577, 0.8432, 0.6086, 0.7832])}

sizemag_rel_vdrops_h12 = {'mbins':array([[24.,25.],[25.,26.],[26.,27.],[27.,28.5]]),
                         'logr0':array([0.6422, 0.9382, 0.8569, 0.8486]),
                         'sigma':array([1.1617, 0.9055, 0.8264, 0.6017])}
mci = bl.mconvert('M1500_to_i.txt')
mcz = bl.mconvert('M1500_to_z.txt')
param_bdrops = {'alpha':-1.6171, 'mstar':-20.5343, 'logr0':0.815, 'sigma':0.766, 'beta':0.204}
param_vdrops = {'alpha':-1.6853, 'mstar':-20.5097, 'logr0':0.876, 'sigma':0.852, 'beta':0.266}
delta_logr0_bdrops = [0.037, -0.027]  # [+,-]
delta_logr0_vdrops = [0.14, -0.05]    # [+,-]
##

## Shen et al. 2003 (Fig. 6, Sersic Re fit)
## late-type: log10(Re[kpc]) = -0.4*alpha*M+(beta-alpha)*log10[1+10**(-0.4*(M-M0))]+gamma
## early-type: log10(Re[kpc]) = -0.4*a*M + b
shen03_late = {'alpha':0.25, 'beta':0.51, 'gamma':-1.71, 'M0':-20.91}
shen03_early = {'a':0.65, 'b':-5.06}

## Courteau et al. 2007 (Table 2 and Fig 4), 
courteau07_disk = {'beta':0.321, 'zp':-2.851}  # for all disk galaxies
courteau07_sd = {'beta':0.254, 'zp':-2.116}  # for Sd type only
# adopt the abs magnitude of Sun in I-band as 4.08 (Binney & Merrifield 1998)

## de Jong & Lacey 2000
dJL00_tot = {'Mstar':-22.17,'rstar':6.09,'sigma':0.28,'beta':0.253}
dJL00_disk = {'Mstar':-22.38, 'rstar':5.93, 'sigma':0.36, 'beta':0.214}

def fsize(size):
   fp = mpl.font_manager.FontProperties(size=size)
   return fp

def RL_LBG(logr0, beta, z, Marr, M0=-21.):
   # returns logR [kpc] as a function of abs mag
   adsize = cc.adsize(z, 1.0, unit='kpc')
   logr = log10(adsize*0.03) + logr0 - 0.4 * beta * (Marr - M0)
   return logr

def plot_LR():
   # first plot Huang et al. 2012
   fig = plt.figure(figsize=(10,6))
   ax = fig.add_subplot(111)
   m0i = -21.0; m0z = -21.0
   logrpeak_b = zeros(4); logrpeak_v = zeros(4)
   mavg_b = zeros(4); mavg_v = zeros(4)
   yerrs_b = zeros((2,4)); yerrs_v = zeros((2,4))
   logsigma_b = sizemag_rel_bdrops_h12['sigma'] / log(10.)
   logsigma_v = sizemag_rel_vdrops_h12['sigma'] / log(10.)
   #for i in range(4):
      # B-dropouts
   #   mavg_b[i] = average(sizemag_rel_bdrops_h12['mbins'][i]) - mci(4.0)
   #   logrpeak_b[i] = sizemag_rel_bdrops_h12['logr0'][i] - 0.4 * param_bdrops['beta'] * \
   #      (mavg_b[i] - m0i)
      # V-dropouts
   #   mavg_v[i] = average(sizemag_rel_vdrops_h12['mbins'][i]) - mcz(5.0)
   #   logrpeak_v[i] = sizemag_rel_vdrops_h12['logr0'][i] - 0.4 * param_vdrops['beta'] * \
   #      (mavg_v[i] - m0z)
   #ax.errorbar(mavg_b, log10(ad_1as_z4*0.03) + logrpeak_b, yerr=logrpeak_b, \
   #   xerr=[0.5,0.5,0.5,0.75],\
   #   fmt='o', ms=10, ecolor='black', mec='black', mfc='black', label='B-dropouts')
   #ax.plot(mavg_b, log10(ad_1as_z4*0.03)+logrpeak_b, 'o', ms=10, c='black')
   #ax.errorbar(mavg_v, log10(ad_1as_z5*0.03) + logrpeak_v, yerr=logrpeak_v, \
   #   xerr=[0.5,0.5,0.5,0.75],\
   #   fmt='^', ms=10, ecolor='black', mec='blue', mfc='blue', label='V-dropouts')
   #ax.plot(mavg_v, log10(ad_1as_z5*0.03)+logrpeak_v, '^', ms=10, c='blue')
   #ax.set_xlabel('absolute magnitude')
   #ax.set_ylabel(r'$R_e$ [kpc]')
   
   # plot the best-fit RL relations
   
   # define the absolute magnitude ranges probed for each curve
   Marr_bdrops = arange(-22.44, -17.3+0.01, 0.01)
   Marr_vdrops = arange(-22.65, -17.93+0.01, 0.01)
   Marr_shen03_etg = arange(-24., -19.+0.01, 0.01)
   Marr_shen03_ltg = arange(-24., -16.+0.01, 0.01)
   Marr_c07 = arange(-24.67, -17.17+0.01, 0.01)
   Marr_dJ00 = arange(-25.5, -16.5+0.01, 0.01)
   ax.semilogy(Marr_bdrops, 10.**RL_LBG(param_bdrops['logr0'], param_bdrops['beta'], 4.0, Marr_bdrops), c='black',\
      lw=3., label='z~4 (This work)')
   # plot what the z=4 RL relation will be given H(z)**-1 evolution
   f1 = sqrt(0.3*(1.+4.)**3+0.7)
   ax.semilogy(Marr_bdrops, 10.**RL_LBG(param_bdrops['logr0'], param_bdrops['beta'], 4.0, Marr_bdrops)*f1,
      c='0.5', lw=3.0, ls='--', label=r'z=0 from $H(z)^{-1}$')
   # plot what the z=4 RL relation will be given H(z)**(-2/3) evolution
   f2 = (0.3*(1.+4.)**3+0.7)**(1./3.)
   ax.semilogy(Marr_bdrops, 10.**RL_LBG(param_bdrops['logr0'], param_bdrops['beta'], 4.0, Marr_bdrops)*f2,
      c='0.5', ls='-', marker=None,  lw=3.0, label=r'z=0 from $H(z)^{-2/3}$')
   ax.semilogy(Marr_vdrops, 10.**RL_LBG(param_vdrops['logr0'], param_vdrops['beta'], 5.0, Marr_vdrops), c='blue', \
      lw=3., label='z~5 (This work)')
   xc1, yc1 = mpl.mlab.poly_between(Marr_bdrops,10.**RL_LBG(param_bdrops['logr0']+delta_logr0_bdrops[0],
      param_bdrops['beta'],4.0,Marr_bdrops),
      10.**RL_LBG(param_bdrops['logr0']+delta_logr0_bdrops[1],param_bdrops['beta'],4.0,Marr_bdrops))
   ax.fill(xc1, yc1, color='black', hatch='/', fill=False, lw=0.5)
   #ax.fill_between(Marr, y1=RL_LBG(param_bdrops['logr0']+delta_logr0_bdrops[0],\
   #   param_bdrops['beta'],4.0,Marr), 
   #   y2=RL_LBG(param_bdrops['logr0']+delta_logr0_bdrops[1],param_bdrops['beta'],4.0,Marr),
   #   color='blue', hatch='/')
   xc2, yc2 = mpl.mlab.poly_between(Marr_vdrops,10.**RL_LBG(param_vdrops['logr0']+delta_logr0_vdrops[0],
      param_vdrops['beta'],5.0,Marr_vdrops),
      10.**RL_LBG(param_vdrops['logr0']+delta_logr0_vdrops[1],param_vdrops['beta'],5.0,Marr_vdrops))
   ax.fill(xc2, yc2, color='blue', hatch='\\', fill=False, lw=0.5)
   #ax.fill_between(Marr, y1=RL_LBG(param_vdrops['logr0']+delta_logr0_vdrops[0],\
   #   param_vdrops['beta'],5.0,Marr),
   #   y2=RL_LBG(param_vdrops['logr0']+delta_logr0_vdrops[1],param_vdrops['beta'],5.0,Marr),
   #   color='red', alpha=0.2)
   # plot the width of the size distribution at the bottom left corner
   yerr1 = param_bdrops['sigma']/log(10.)
   yerr2 = param_vdrops['sigma']/log(10.)
   yerrpow = 1.2
   yerr_upper1 = 10.**(yerrpow+yerr1)-10.**(yerrpow)
   yerr_upper2 = 10.**(yerrpow+yerr2)-10.**(yerrpow)
   yerr_lower1 = 10.**(yerrpow)-10.**(yerrpow-yerr1)
   yerr_lower2 = 10.**(yerrpow)-10.**(yerrpow-yerr2)
   ax.set_yscale('log')
   ax.errorbar([-18.0], 10.**yerrpow, yerr=[[yerr_lower1],[yerr_upper1]], fmt=None, ecolor='black')
   ax.errorbar([-17.7], 10.**yerrpow, yerr=[[yerr_lower2],[yerr_upper2]], fmt=None, ecolor='blue')

   # plot the Shen et al. 2003 relations, when Sersic Re is used
   RL_s03_late = -0.4*shen03_late['alpha']*Marr_shen03_ltg+(shen03_late['beta']-shen03_late['alpha'])*log10(1.+10.**(-0.4*(Marr_shen03_ltg-shen03_late['M0'])))+shen03_late['gamma']
   ax.semilogy(Marr_shen03_ltg, 10.**RL_s03_late, '--', c='green', lw=2.0, label='S03 z~0 late')
   RL_s03_early = -0.4*shen03_early['a']*Marr_shen03_etg + shen03_early['b']
   ax.semilogy(Marr_shen03_etg, 10.**RL_s03_early, '--', c='red', lw=2.0, label='S03 z~0 early')

   # plot the Courteau et al. 2007 relation for disk galaxies
   RL_c07 = -0.4*courteau07_disk['beta']*(Marr_c07-4.08)+courteau07_disk['zp']
   ax.semilogy(Marr_c07, 10.**RL_c07, '-.', c='magenta', lw=2.0, label='C07 z~0')
   
   # plot de Jong & Lacey 2000 relation for disk galaxies (disk+bulge & disk-only)
   RL_dJL00_tot = log10(dJL00_tot['rstar'])-0.4*dJL00_tot['beta']*(Marr_dJ00-dJL00_tot['Mstar'])
   RL_dJL00_disk = log10(dJL00_disk['rstar'])-0.4*dJL00_disk['beta']*(Marr_dJ00-dJL00_disk['Mstar'])
   #print Marr, 10.**RL_dJL00_tot
   ax.semilogy(Marr_dJ00, 10.**RL_dJL00_tot, '-', c='orange', lw=2.0, label='dJL00 total')
   ax.semilogy(Marr_dJ00, 10.**RL_dJL00_disk, '--', c='orange', lw=2.0, label='dJL00 disk')
   
   #ytick = array([-1.,-0.5,0.,0.5,1.])
   ytick = array([0.1,1.,10.])
   #ax.set_yticks(minor=True)
   ax.set_yticklabels([0.1,1,10,100],[0.1,1,10,100],size=14)
   ax.set_xlim(-24., -17.); ax.set_ylim(0.1, 50.)
   xticks = ax.get_xticks()
   ax.set_xticks(xticks[1:])
   ax.set_xticklabels(xticks[1:])
   ax.legend(loc='lower center', ncol=3)
   ax.set_xlabel('Absolute Magnitude')
   ax.set_ylabel('Size [kpc]')
   return ax

