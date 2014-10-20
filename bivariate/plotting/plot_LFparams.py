#!/usr/bin/env python

from numpy import *
import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
from bivariate.plotting import plot_mcmc as pmc
from pygoods import *
from bivariate import fit_lbg as fl
from my_mplfonts import Helvetica
from bivariate import schechter
import pmscolors

pms = pmscolors.pmscolors()
## U-drops updated 14/01/04
paru_all = np.array([-2.05496823, -21.22223052, 9.776e-04, 0.66511482, 0.68345048, 0.20026905])
mcmc_udrops_cat = '/Users/khuang/Dropbox/Research/bivariate/bivariate_fit/udrops_fitting/mcmc_runs/udrops_mcmc_140103.fits'
## i-drops (130715_1)
pari_all = np.array([-1.84221235, -20.32331588, 3.406e-03, 0.67513577, 0.28752996, 0.1814661])
mcmc_idrops_cat = '/Users/khuang/Dropbox/Research/bivariate/bivariate_fit/idrops_fitting/mcmc_runs/idrops_mcmc_130730_100iter_2.fits'

"""
Plot the derived Schechter function parameters from previous studies at z=4 and 5
"""
z3par = {'S99':[-1.6, -21.2, 1.6],
         'vdB10':[-1.65, -20.93, 1.79],
         'AS00':[-1.57, -21.14, 4.4], # Adelberger & Steidel 2000
         'S06':[-1.43, -20.78, 1.7], # Sawicki & Thompson 2006
         'R09':[-1.73, -20.85, 1.71] # Reddy & Steidel 2009
         }
z3parerr = {'S99':[0.13, 0.15, 0],
            'vdB10':[[0.11,0.12], [0.13,0.14], [0.38,0.51]],
            'AS00':[0.11, 0.14, 1.0],
            'S06':[[0.09,0.17], [0.09,0.13], [0.32,0.59]],
            'R09':[0.13, 0.14, 0.53]
            }
z3_Maglim = {'H13':[-22.41, -18.12],
             'S99':[-22.75, -20.],
             'vdB10':[-23.20, -18.70],
             'AS00':[-22.85, -18.35],
             'S06':[-23., -18.],
             'R09':[-23.0, -18.5]}
z4par = {'B07':[-1.73,-20.98,1.3],
         'Y06':[-1.82,-21.14,1.5],
         'Be06':[-1.6,-20.7,1.3],
         'ST06':[-1.26,-21.0,0.9],
         'G05':[-1.64,-21.2,1.2],
         'Ou04':[-2.2,-21.0,1.2],
         'S99':[-1.6,-21.2,1.1],
         'vdB10':[-1.56,-20.84,1.36]}
z4parerr = {'B07':[0.05,0.1, 0.2],
            'Y06':[[0.09,0.09],[0.14,0.15],[0.4,0.4]],
            'Be06':[0.,0.,0.],
            'ST06':[[0.4,0.36],[0.4,0.5],[0.5,0.5]],
            'G05':[0.1,0.04,0.03],
            'Ou04':[0.2,0.1,0.2],
            'S99':[0.,0.,0.],
            'vdB10':[[0.08,0.08],[0.09,0.09],[0.23,0.2]]}
z4_Maglim = {'H13':[-22.5, -17.5],
             'B07':[-22.26, -16.26],
             'Y06':[-23.16, -19.16],
             'Be06':[-22.0, -16.39],
             'ST06':[-23.0, -18.5],
             'G05':[-23.0, -19.0],
             'Ou04':[-22.5, -20.0],
             'S99':[-23.0, -19.0],
             'vdB10':[-22.6, -18.7]}
z5par = {'B07':[-1.66,-20.64,1.0],
         'Y06':[-1.82,-20.72,1.2],
         'Be06':[-1.6,-20.55,0.9],
         'G05':[-1.51,-21.06,0.83],
         'Ou04':[-1.6,-20.3,2.4],
         'Oe07':[-1.54,-20.78,0.9],
         'I07':[-1.48,-21.28,0.41],
         'vdB10':[-1.65,-20.94,0.83]}
z5parerr = {'B07':[0.09,0.13,0.3],
            'Y06':[[0.,0.],[0.16,0.14],[0.4,0.3]],
            'Be06':[0.,0.,0.],
            'G05':[0.18,0.05,0.03],
            'Ou04':[0.,0.2,1.0],
            'Oe07':[0.1,0.16,0.3],
            'I07':[[0.38,0.32],[0.38,0.38],[0.29,0.3]],
            'vdB10':[[0.09,0.08],[0.1,0.11],[0.15,0.14]]}
z5_Maglim = {'H13':[-22.89, -17.89],
             'B07':[-21.66, -17.16],
             'Y06':[-22.84, -20.34],
             'Be06':[-22.39, -16.87],
             'G05':[-23.0, -19.5],
             'Ou04':[-22.5, -20.5],
             'Oe07':[-21.0, -17.0],
             'I07':[-23.08, -19.89],
             'vdB10':[-22.6, -19.0]}
z6par = {'B07':[-1.74,-20.23,1.4],
         'B04':[-1.6, -20.87, 0.23], # Bunker et al. 2004
         'Y04':[-1.8, -21.03, 0.46], # Yan & Windhorst 2004
         'M05':[-1.8, -20.83, 0.4], # Malhotra et al. 2005
         'S11':[-1.87, -20.25, 1.77], # Su et al. 2011
         'Be06':[-1.6, -21.3, 0.4],  # Beckwith et al. 2006
         'M09':[-1.71, -20.04, 1.8] # McLure et al. 2009
         }
z6parerr = {'B07':[0.19, 0.1, 0.2],
            'S11':[[0.14,0.14], [0.23,0.23], [0.49,0.62]],
            'M09':[0.11, 0.12, 0.5]}
z6_Maglim = {'H13':[-22.67, -19.17],
             'B07':[-22.13, -17.88],
             'B04':[-22.97, -17.97],
             'Y04':[-23.67, -17.0],
             'M05':[-21.61, -19.17],
             'S11':[-22., -17.],
             'Be06':[-21.97, -18.22],
             'M09':[-22.5, -20.6]}
plotkwarg = {'B07':{'label':'B07','marker':'^'},
             'Y06':{'label':'Y06','marker':'.'},
             'Be06':{'label':'Be06','marker':'p'},
             'ST06':{'label':'ST06','marker':'D'},
             'G05':{'label':'G05','marker':'H'},
             'Ou04':{'label':'Ou04','marker':'v'},
             'S99':{'label':'S99','marker':'d'},
             'Oe07':{'label':'Oe07','marker':'s'},
             'I07':{'label':'I07','marker':'o'},
             'vdB10':{'label':'vdB10','marker':'h'},
             'B04':{'label':'B04','marker':'<'},
             'Y04':{'label':'Y04','marker':'>'},
             'M05':{'label':'M05','marker':'1','mew':3},
             'S11':{'label':'S11','marker':'2','mew':3},
             'S06':{'label':'S06','marker':'3','mew':3},
             'R09':{'label':'R09','marker':'4','mew':3},
             'M09':{'label':'M09','marker':'+','mew':3},
             'AS00':{'label':'AS00','marker':'x','mew':3}
             }
colordic = {'B07':'blue',
            'Y06':'green',
            'Be06':'red',
            'ST06':'orange',
            'G05':'cyan',
            'Ou04':'#C44A2C',
            'S99':'#6C1B23',
            'Oe07':'#2976C8',
            'I07':'#637F3A',
            'vdB10':'#6B2138',
            'B04':'0.4',
            'Y04':'#B279C8',
            'M05':'#F98972',
            'S11':'#004731',
            'S06':'#77B6D0',
            'R09':'#679000',
            'M09':'#009878',
            'AS00':'#CBD34C'
            }
errbarkw = {'fmt':None, 'label':'_nolegend_'}
par_dict = {'z3':z3par, 'z4':z4par, 'z5':z5par, 'z6':z6par}
parerr_dict = {'z3':z3parerr, 'z4':z4parerr, 'z5':z5parerr, 'z6':z6parerr}

def fpsize(s):
   fp = matplotlib.font_manager.FontProperties(size=s)
   return fp


def mkerrbarkw(authorkey):
   kw = errbarkw.copy()
   kw['ecolor'] = colordic[authorkey]
   return kw


def mkplotkw(authorkey):
   kw = plotkwarg[authorkey].copy()
   kw['mfc'] = colordic[authorkey]
   kw['mec'] = colordic[authorkey]
   kw['linestyle'] = 'None'
   kw['ms'] = 10
   return kw

def plotparam_z(z, authorkey, i0, i1, ax=None):
   # i0, i1: index of array (alpha, Mstar, phistar)
   key = 'z%d' % z
   par = par_dict[key]
   #print "authorkey", authorkey
   parerr = parerr_dict[key]
   if ax == None:
      fig = plt.figure()
      ax = fig.add_subplot(111)
   if parerr.has_key(authorkey):
      if len(shape(parerr[authorkey])) == 1:
         ax.errorbar([par[authorkey][i0]], [par[authorkey][i1]],
                     xerr=[parerr[authorkey][i0]],
                     yerr=[parerr[authorkey][i1]],
                     **mkerrbarkw(authorkey))
      else:
         ax.errorbar([par[authorkey][i0]], [par[authorkey][i1]],
                     xerr=reshape(parerr[authorkey][i0],(2,1)),
                     yerr=reshape(parerr[authorkey][i1],(2,1)),
                     **mkerrbarkw(authorkey))
   ax.plot([par[authorkey][i0]], [par[authorkey][i1]],
           **(mkplotkw(authorkey)))
   return ax

def plotparam_z4(authorkey, ax, i0, i1):
   # i0, i1: index of array (alpha, Mstar, phistar)
   if len(shape(z4parerr[authorkey])) == 1:
      ax.errorbar([z4par[authorkey][i0]],[z4par[authorkey][i1]],xerr=[z4parerr[authorkey][i0]],
         yerr=[z4parerr[authorkey][i1]],**mkerrbarkw(authorkey))
   else:
      ax.errorbar([z4par[authorkey][i0]],[z4par[authorkey][i1]],
         xerr=reshape(z4parerr[authorkey][i0],(2,1)),
         yerr=reshape(z4parerr[authorkey][i1],(2,1)),
         **mkerrbarkw(authorkey))
   ax.plot([z4par[authorkey][i0]],[z4par[authorkey][i1]],**(mkplotkw(authorkey)))
   return ax


def plotparam_z5(authorkey, ax, i0, i1):
   # i0, i1: index of array (alpha, Mstar, phistar)
   if len(shape(z5parerr[authorkey])) == 1:
      ax.errorbar([z5par[authorkey][i0]],[z5par[authorkey][i1]],xerr=[z5parerr[authorkey][i0]],
         yerr=[z5parerr[authorkey][i1]],**mkerrbarkw(authorkey))
   else:
      ax.errorbar([z5par[authorkey][i0]],[z5par[authorkey][i1]],
         xerr=reshape(z5parerr[authorkey][i0],(2,1)),
         yerr=reshape(z5parerr[authorkey][i1],(2,1)),
         **mkerrbarkw(authorkey))
   ax.plot([z5par[authorkey][i0]],[z5par[authorkey][i1]],**mkplotkw(authorkey))
   return ax

def plot_LFparams(z, mcmc_cat, ax0, ax1, alpha, Mstar, phistar, ticksize=14,
                  plot_mcmc=False, 
                  paramlims=dict(alpha=array([-2.,-1.]), 
                                 mstar=array([-22.0,-20.0]), 
                                 phistar=array([0.1e-3,5.0e-3])),
                  parampix=dict(alpha=0.01,mstar=0.01,phistar=5.e-5),
                  paramscale=dict(alpha=0.5,mstar=0.2,phistar=200.),
                  htable=pmc.htable_b,
                  legendloc=4):
   """
   A general function to plot parameters at all redshifts.
   """
   # First read the proper parameter lists
   key = 'z%d' % z
   par = par_dict[key]
   parerr = parerr_dict[key]
   if plot_mcmc:
      c = Ftable(mcmc_cat)
   ax0.plot([alpha], [Mstar], '*', ms=18, c='black', label='This work')
   # plot confidence contours
   if plot_mcmc:
      xlim = paramlims['alpha']
      ylim = paramlims['mstar']
      dx   = parampix['alpha']
      dy   = parampix['mstar']
      sx   = paramscale['alpha']
      sy   = paramscale['mstar']
      pdf1 = pmc.scatter2pdf(xlim, ylim, dx, dy, c.alpha, c.mstar, sx, sy, 
                          htable[(0,1)], plot=False)
      # find 1- and 2-sigma levels
      levels1 = pmc.conflevel(pdf1)
      xgrid = arange(xlim[0], xlim[1], dx)
      ygrid = arange(ylim[0], ylim[1], dy)
      ax0.contour(xgrid, ygrid, pdf1, levels=levels1, colors='black')
   # Now plot the parameters alpha & Mstar from literature
   for k in par.keys():
      ax0 = plotparam_z(z, k, 0, 1, ax=ax0)
   ax0.set_xlabel(r'$\alpha$', font_properties=Helvetica(18)) 
   ax0.set_ylabel('$M^*$', font_properties=Helvetica(18))
   ax0.set_xticks(arange(-2.6,-0.9,0.3))
   ax0.legend(loc=legendloc, prop=Helvetica(14), numpoints=1)
   if z==3:
     ax0.text(0.1, 0.9, r'$z \sim 3.4$', transform=ax0.transAxes,
              color='black', font_properties=Helvetica(22))
   else:
    ax0.text(0.1, 0.9, r'$z \sim %d$' % z, transform=ax0.transAxes, 
            color='black', font_properties=Helvetica(22))
   # Now plot Mstar v.s. phistar
   ax1.plot([phistar*1.e3], [Mstar], '*', ms=18, c='black', label='This work')
   if plot_mcmc:
      xlim = paramlims['phistar']
      ylim = paramlims['mstar']
      dx   = parampix['phistar']
      dy   = parampix['mstar']
      sx   = paramscale['phistar']
      sy   = paramscale['mstar']
      pdf2 = pmc.scatter2pdf(xlim, ylim, dx, dy, c.phistar, c.mstar, sx, sy, 
                          htable[(1,2)], plot=False)
      levels2 = pmc.conflevel(pdf2); print levels2
      xgrid = arange(xlim[0], xlim[1], dx) * 1.e3
      ygrid = arange(ylim[0], ylim[1], dy)
      ax1.contour(xgrid, ygrid, pdf2, levels=levels2, colors='black')
   # Now plot values in the literature
   for k in par.keys():
      ax1 = plotparam_z(z, k, 2, 1, ax=ax1)
   ax1.set_ylabel('') 
   ax1.set_xlabel(r'$10^{-3}\ \phi^*\ \mathrm{Mpc}^{-3}$', 
                  font_properties=Helvetica(18))
   ax1.set_yticks(arange(-21.3,-20.5,0.2))
   ax1.set_yticklabels([])
   ax1.legend(loc=legendloc, prop=Helvetica(14), numpoints=1)
   if z==3:
     ax1.text(0.1, 0.9, r'$z \sim 3.4$', transform=ax1.transAxes,
              color='black', font_properties=Helvetica(22))
   else:
     ax1.text(0.1, 0.9, r'$z \sim %d$' % z, transform=ax1.transAxes, 
            color='black', font_properties=Helvetica(22))
   return ax0, ax1

def plot_LFparams_bdrops_vdrops(parb_all, parv_all):
   # parb_all, parv_all are in the format [alpha,mstar,phistar,logr0,sigma,beta]
   #phistar_b = fl.phistar(parb, 'b')
   #phistar_v = fl.phistar(parv, 'v')
   #matplotlib.rcParams['xtick.labelsize'] = 17
   #matplotlib.rcParams['ytick.labelsize'] = 17
   #matplotlib.rcParams['axes.labelsize'] = 19
   alpha_lim = array([-2.6,-0.4])
   mstar_lim = array([-21.8, -19.8])
   phistar_lim = array([0., 3.0])
   alpha_ticks = linspace(alpha_lim[0], alpha_lim[1], 7)[1:-1]
   mstar_ticks = linspace(mstar_lim[0], mstar_lim[1], 7)[1:-1]
   phistar_ticks = linspace(phistar_lim[0], phistar_lim[1], 7)[1:-1]
   phistar_b = parb_all[2]
   phistar_v = parv_all[2]
   fig = plt.figure(figsize=(10,9))
   ax0 = fig.add_axes([0.1, 0.55, 0.4, 0.4])
   ax1 = fig.add_axes([0.5, 0.55, 0.4, 0.4])
   ax2 = fig.add_axes([0.1, 0.08, 0.4, 0.4])
   ax3 = fig.add_axes([0.5, 0.08, 0.4, 0.4])
   bdrops_par = array([parb_all[0], parb_all[1], phistar_b])
   vdrops_par = array([parv_all[0], parv_all[1], phistar_v])
   print "bdrops_par", bdrops_par
   print "vdrops_par", vdrops_par
   ax0, ax1 = plot_LFparams(4, 'mcmc_runs/bdrops_mcmc_out.fits', ax0, ax1, 
                            *bdrops_par, plot_mcmc=True, 
                            paramlims=dict(alpha=alpha_lim,mstar=mstar_lim,
                                           phistar=phistar_lim/1.e3),
                            parampix=pmc.parampix_b, 
                            paramscale=pmc.paramscale_b,
                            htable=pmc.htable_b)
   ax2, ax3 = plot_LFparams(5, 'mcmc_runs/vdrops_mcmc_120906.fits', ax2, ax3,     
                            *vdrops_par, plot_mcmc=True,
                            paramlims=dict(alpha=alpha_lim,mstar=mstar_lim,
                                           phistar=phistar_lim/1.e3),
                            parampix=pmc.parampix_v,
                            paramscale=pmc.paramscale_v,
                            htable=pmc.htable_v)
   axes = [ax0,ax1,ax2,ax3]
   ax0.tick_params(top='off', right='off')
   ax1.tick_params(top='off', right='off')
   ax2.tick_params(top='off', right='off')
   ax3.tick_params(top='off', right='off')
   ax0.set_ylim(mstar_lim); ax1.set_ylim(mstar_lim)
   ax2.set_ylim(mstar_lim); ax3.set_ylim(mstar_lim)
   ax0.set_xlim(alpha_lim); ax1.set_xlim(phistar_lim)
   ax2.set_xlim(alpha_lim); ax3.set_xlim(phistar_lim)
   for ax in [ax1,ax3]:
      ax.set_xticks(phistar_ticks)
      ax.set_xticklabels(map(lambda x:'%.1f'%x, phistar_ticks), 
                         font_properties=Helvetica(14))
   for ax in [ax0,ax2]:
      ax.set_xticks(alpha_ticks)
      ax.set_xticklabels(map(lambda x:'%.1f'%x, alpha_ticks), 
                         font_properties=Helvetica(14))
   for ax in axes:
      ax.set_yticks(mstar_ticks)   
   for ax in [ax0, ax2]:
      ax.set_yticklabels(map(lambda y:'%.1f'%y, mstar_ticks), 
                         font_properties=Helvetica(14))
   return fig

def plot_LFparams_udrops_idrops(plot_mcmc=False):
   # parb_all, parv_all are in the format [alpha,mstar,phistar,logr0,sigma,beta]
   #phistar_b = fl.phistar(parb, 'b')
   #phistar_v = fl.phistar(parv, 'v')
   #matplotlib.rcParams['xtick.labelsize'] = 17
   #matplotlib.rcParams['ytick.labelsize'] = 17
   #matplotlib.rcParams['axes.labelsize'] = 19
   alpha_lim = [-2.2,-1.0]
   mstar_lim = [-21.5, -19.0]
   phistar_lim = [0., 8.0]
   alpha_ticks = linspace(alpha_lim[0], alpha_lim[1], 7)[1:-1]
   mstar_ticks = linspace(mstar_lim[0], mstar_lim[1], 7)[1:-1]
   phistar_ticks = linspace(phistar_lim[0], phistar_lim[1], 7)[1:-1]
   phistar_u = paru_all[2]
   phistar_i = pari_all[2]
   fig = plt.figure(figsize=(10,9))
   ax0 = fig.add_axes([0.13, 0.55, 0.4, 0.4])
   ax1 = fig.add_axes([0.53, 0.55, 0.4, 0.4])
   ax2 = fig.add_axes([0.13, 0.08, 0.4, 0.4])
   ax3 = fig.add_axes([0.53, 0.08, 0.4, 0.4])
   udrops_par = array([paru_all[0], paru_all[1], phistar_u])
   idrops_par = array([pari_all[0], pari_all[1], phistar_i])
   print "udrops_par", udrops_par
   print "idrops_par", idrops_par
   ax0, ax1 = plot_LFparams(3, mcmc_udrops_cat, 
                            ax0, ax1, 
                            *udrops_par, plot_mcmc=plot_mcmc,
                            paramlims=pmc.paramlims_u,
                            parampix=pmc.parampix_u,
                            paramscale=pmc.paramscale_u,
                            htable=pmc.htable_u,
                            legendloc=1)
   ax2, ax3 = plot_LFparams(6, mcmc_idrops_cat, 
                            ax2, ax3,     
                            *idrops_par, plot_mcmc=plot_mcmc,
                            paramlims=pmc.paramlims_i,
                            parampix=pmc.parampix_i,
                            paramscale=pmc.paramscale_i,
                            htable=pmc.htable_i)
   axes = [ax0,ax1,ax2,ax3]
   ax0.tick_params(top='off', right='off')
   ax1.tick_params(top='off', right='off')
   ax2.tick_params(top='off', right='off')
   ax3.tick_params(top='off', right='off')
   ax0.set_ylim(mstar_lim); ax1.set_ylim(mstar_lim)
   ax2.set_ylim(mstar_lim); ax3.set_ylim(mstar_lim)
   ax0.set_xlim(alpha_lim); ax1.set_xlim(phistar_lim)
   ax2.set_xlim(alpha_lim); ax3.set_xlim(phistar_lim)
   for ax in [ax1,ax3]:
      ax.set_xticks(phistar_ticks)
      ax.set_xticklabels(map(lambda x:'%.1f'%x, phistar_ticks), 
                         font_properties=Helvetica(16))
   for ax in [ax0,ax2]:
      ax.set_xticks(alpha_ticks)
      ax.set_xticklabels(map(lambda x:'%.1f'%x, alpha_ticks), 
                         font_properties=Helvetica(16))
   for ax in axes:
      ax.set_yticks(mstar_ticks)  
   for ax in [ax0, ax2]: 
      ax.set_yticklabels(map(lambda y:'%.1f'%y, mstar_ticks), 
                       font_properties=Helvetica(16))
   return fig


def plot_LFparams_z4(ax0, ax1, alpha, Mstar, phistar, ticksize=14):
   #fig = plt.figure(figsize=(12,6))
   #grid = ImageGrid(fig, 111,
   #       nrows_ncols=(1,2), axes_pad=0.0)
   #ax0 = grid[0]
   #c = sextractor('mcmc_runs/bdrops_mcmc_out.cat')
   c = Ftable('bdrops_mcmc_090612.fits')
   #ax0 = fig.add_subplot(1,2,1) # x: alpha; y:mstar
   ax0.plot([alpha], [Mstar], '*', ms=18, c='black', label='This work')
   # plot confidence contours
   xlim = pmc.paramlims_b['alpha']
   ylim = pmc.paramlims_b['mstar']
   dx   = pmc.parampix_b['alpha']
   dy   = pmc.parampix_b['mstar']
   sx   = pmc.paramscale_b['alpha']
   sy   = pmc.paramscale_b['mstar']
   pdf1 = pmc.scatter2pdf(xlim, ylim, dx, dy, c.alpha, c.mstar, sx, sy, pmc.htable_b[(0,1)],
      plot=False)
   # find 1- and 2-sigma levels
   levels1 = pmc.conflevel(pdf1)
   xgrid = arange(xlim[0], xlim[1], dx)
   ygrid = arange(ylim[0], ylim[1], dy)
   ax0.contour(xgrid, ygrid, pdf1, levels=levels1, colors='black')

   ax0 = plotparam_z4('B07', ax0, 0, 1)
   ax0 = plotparam_z4('Y06', ax0, 0, 1)
   ax0 = plotparam_z4('Be06', ax0, 0, 1)
   ax0 = plotparam_z4('ST06', ax0, 0, 1)
   ax0 = plotparam_z4('G05', ax0, 0, 1)
   ax0 = plotparam_z4('Ou04', ax0, 0, 1)
   ax0 = plotparam_z4('S99', ax0, 0, 1)
   ax0 = plotparam_z4('vdB10', ax0, 0, 1)

   ax0.set_xlabel(r'$\alpha$'); ax0.set_ylabel('$M^*$')
   ax0.set_xticks(arange(-2.6,-0.9,0.3))
   #ax0.set_xticklabels(ax0.get_xticklabels(), font_properties=fpsize(ticksize))
   #ax0.set_yticklabels(ax0.get_yticklabels(), font_properties=fpsize(ticksize))
   ax0.legend(loc=4, prop=fpsize(10.), numpoints=1)
   ax0.text(0.1, 0.9, r'$z \sim 4$', size=16, transform=ax0.transAxes, color='black')

   #ax1 = grid[1]
   #ax1 = fig.add_subplot(1,2,2)
   # x: mstar; y:phistar
   ax1.plot([phistar*1.e3], [Mstar], '*', ms=18, c='black', label='This work')
   # plot confidence intervals
   xlim = pmc.paramlims_b['phistar']
   ylim = pmc.paramlims_b['mstar']
   dx   = pmc.parampix_b['phistar']
   dy   = pmc.parampix_b['mstar']
   sx   = pmc.paramscale_b['phistar']
   sy   = pmc.paramscale_b['mstar']
   pdf2 = pmc.scatter2pdf(xlim, ylim, dx, dy, c.phistar, c.mstar, sx, sy, pmc.htable_b[(1,2)],
      plot=False)
   # find 1- and 2-sigma levels
   levels2 = pmc.conflevel(pdf2); print levels2
   xgrid = arange(xlim[0], xlim[1], dx) * 1.e3
   ygrid = arange(ylim[0], ylim[1], dy)
   ax1.contour(xgrid, ygrid, pdf2, levels=levels2, colors='black')

   ax1 = plotparam_z4('B07', ax1, 2, 1)
   ax1 = plotparam_z4('Y06', ax1, 2, 1)
   ax1 = plotparam_z4('Be06', ax1, 2, 1)
   ax1 = plotparam_z4('ST06', ax1, 2, 1)
   ax1 = plotparam_z4('G05', ax1, 2, 1)
   ax1 = plotparam_z4('Ou04', ax1, 2, 1)
   ax1 = plotparam_z4('S99', ax1, 2, 1)
   ax1 = plotparam_z4('vdB10', ax1, 2, 1)

   ax1.set_ylabel(''); ax1.set_xlabel(r'$10^{-3}\ \phi^*\ \mathrm{Mpc}^{-3}$')
   ax1.set_yticks(arange(-21.3,-20.5,0.2))
   ax1.set_yticklabels([])
   ax1.legend(loc=4, prop=fpsize(10.), numpoints=1)
   ax1.text(0.1, 0.9, r'$z \sim 4$', size=16, transform=ax1.transAxes, color='black')

   return ax0, ax1


def plot_LFparams_z5(ax0, ax1, alpha, Mstar, phistar, ticksize=14):
   #fig = plt.figure(figsize=(12,6))
   #grid = ImageGrid(fig, 111,
   #       nrows_ncols=(1,2), axes_pad=0.0)
   #ax0 = grid[0]
   #ax0 = fig.add_subplot(1,2,1)
   #c = sextractor('mcmc_runs/vdrops_mcmc_out.cat')
   c = Ftable('vdrops_mcmc_090612.fits')
   #ax0 = fig.add_subplot(1,2,1) # x: alpha; y:mstar
   ax0.plot([alpha], [Mstar], '*', ms=18, c='black', label='This work')
   # plot confidence contours
   xlim = pmc.paramlims_v['alpha']
   ylim = pmc.paramlims_v['mstar']
   dx   = pmc.parampix_v['alpha']
   dy   = pmc.parampix_v['mstar']
   sx   = pmc.paramscale_v['alpha']
   sy   = pmc.paramscale_v['mstar']
   pdf1 = pmc.scatter2pdf(xlim, ylim, dx, dy, c.alpha, c.mstar, sx, sy, pmc.htable_v[(0,1)],
      plot=False)
   # find 1- and 2-sigma levels
   levels1 = pmc.conflevel(pdf1)
   xgrid = arange(xlim[0], xlim[1], dx)
   ygrid = arange(ylim[0], ylim[1], dy)
   ax0.contour(xgrid, ygrid, pdf1, levels=levels1, colors='black')

   ax0 = plotparam_z5('B07', ax0, 0, 1)
   ax0 = plotparam_z5('Y06', ax0, 0, 1)
   ax0 = plotparam_z5('Be06', ax0, 0, 1)
   ax0 = plotparam_z5('G05', ax0, 0, 1)
   ax0 = plotparam_z5('Ou04', ax0, 0, 1)
   ax0 = plotparam_z5('Oe07', ax0, 0, 1)
   ax0 = plotparam_z5('I07', ax0, 0, 1)
   ax0 = plotparam_z5('vdB10', ax0, 0, 1)

   ax0.plot([alpha], [Mstar], '*', ms=18, c='black', label='This work')

   ax0.set_xlabel(r'$\alpha$'); ax0.set_ylabel('$M^*$')
   ax0.set_xticks(arange(-2.3,-0.8,0.2))
   ax0.set_yticks(arange(-22.,-19.,0.4))
   ax0.set_ylim(-21.8,-19.5)
   ax0.legend(loc=4, prop=fpsize(10.), numpoints=1)
   ax0.text(0.1, 0.9, r'$z \sim 5$', size=16, transform=ax0.transAxes, color='black')

   #ax1 = grid[1]
   #ax1 = fig.add_subplot(1,2,2)
   # plot confidence intervals
   xlim = pmc.paramlims_v['phistar']
   ylim = pmc.paramlims_v['mstar']
   dx   = pmc.parampix_v['phistar']
   dy   = pmc.parampix_v['mstar']
   sx   = pmc.paramscale_v['phistar']
   sy   = pmc.paramscale_v['mstar']
   pdf2 = pmc.scatter2pdf(xlim, ylim, dx, dy, c.phistar, c.mstar, sx, sy, pmc.htable_v[(1,2)],
      plot=False)
   # find 1- and 2-sigma levels
   levels2 = pmc.conflevel(pdf2); print levels2
   xgrid = arange(xlim[0], xlim[1], dx) * 1.e3
   ygrid = arange(ylim[0], ylim[1], dy)
   ax1.contour(xgrid, ygrid, pdf2, levels=levels2, colors='black')

   ax1 = plotparam_z5('B07', ax1, 2, 1)
   ax1 = plotparam_z5('Y06', ax1, 2, 1)
   ax1 = plotparam_z5('Be06', ax1, 2, 1)
   ax1 = plotparam_z5('G05', ax1, 2, 1)
   ax1 = plotparam_z5('Ou04', ax1, 2, 1)
   ax1 = plotparam_z5('Oe07', ax1, 2, 1)
   ax1 = plotparam_z5('I07', ax1, 2, 1)
   ax1 = plotparam_z5('vdB10', ax1, 2, 1)

   ax1.plot([phistar*1.e3], [Mstar], '*', ms=18, c='black', label='This work')

   ax1.set_ylabel(''); ax1.set_xlabel(r'$10^{-3}\ \phi^*\ \mathrm{Mpc}^{-3}$')
   ax1.set_yticks(arange(-21.8,-19.7,0.3))
   ax1.set_yticklabels([])
   ax1.set_ylim(-21.8, -19.7)
   ax1.set_xlim(0.0, 3.0)
   ax1.legend(loc=4, prop=fpsize(10.), numpoints=1)
   ax1.text(0.1, 0.9, r'$z \sim 5$', size=16, transform=ax1.transAxes, color='black')

   return ax0, ax1

def plot_LF_compare(drop, alpha, mstar, phistar, alpha_err=None, 
                    mstar_err=None, 
               phistar_err=None, ndraw=1000):
  """
  Plot the Schechter function from the literature.
  """
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.set_yscale('log')
  if drop == 'b':
    par_dic = z4par
    parerr_dic = z4parerr
    Maglim_dic = z4_Maglim
    dropname = 'B-dropouts'
  elif drop == 'v':
    par_dic = z5par
    parerr_dic = z5parerr
    Maglim_dic = z5_Maglim
    dropname = 'V-dropouts'
  elif drop == 'u':
    par_dic = z3par
    parerr_dic = z3parerr
    Maglim_dic = z3_Maglim
    dropname = 'U-dropouts'
  elif drop == 'i':
    par_dic = z6par
    parerr_dic = z6parerr
    Maglim_dic = z6_Maglim
    dropname = 'i-dropouts'
  M1500_lims = Maglim_dic['H13']
  M1500_arr = arange(M1500_lims[0], M1500_lims[1], 0.02)
  my_lf = phistar * schechter.mag_schechter(M1500_arr, mstar, alpha)
  ax.plot(M1500_arr, my_lf, label='This work', lw=3.0, c='black')
  for k in par_dic.keys():
    M1500_lims = Maglim_dic[k]
    M1500_arr = arange(M1500_lims[0], M1500_lims[1], 0.02)
    alpha_k, mstar_k, phistar_k = par_dic[k]
    phistar_k = phistar_k * 1.e-3
    lf_k = phistar_k * schechter.mag_schechter(M1500_arr, mstar_k, alpha_k)
    ax.plot(M1500_arr, lf_k, lw=1.0, label=k)
  plt.xticks(font_properties=Helvetica(14))
  plt.yticks(font_properties=Helvetica(14))
  plt.xlabel('$M_{1500}$', font_properties=Helvetica(18))
  plt.ylabel('$\phi(M)$ (Mpc$^{-3}$)', font_properties=Helvetica(18))
  plt.legend(loc=4, prop=Helvetica(13))
  plt.title(dropname, font_properties=Helvetica(18))
  if alpha_err != None:
    alpha_drawn = random.uniform(alpha_err[0], alpha_err[1], size=ndraw)
    mstar_drawn = random.uniform(mstar_err[0], mstar_err[1], size=ndraw)
    phistar_drawn = random.uniform(phistar_err[0], phistar_err[1], size=ndraw)
    M1500_lims = Maglim_dic['H13']
    M1500_arr = arange(M1500_lims[0], M1500_lims[1], 0.02)
    my_lf_drawn = []
    for i in range(ndraw):
      my_lf_drawn += [phistar_drawn[i]*schechter.mag_schechter(M1500_arr, mstar_drawn[i], alpha_drawn[i])]
    my_lf_drawn = array(my_lf_drawn)
    my_lf_max = amax(my_lf_drawn, axis=0)
    my_lf_min = amin(my_lf_drawn, axis=0)
    ax.fill_between(M1500_arr, my_lf_min, my_lf_max, color='black',
                    hatch='\\', facecolor='none', lw=2.0)
  z0_dic = dict(u=3.4,b=4.0,v=5.0,i=6.0)
  ax.text(0.05, 0.95, r'$z\sim %.1f$'%z0_dic[drop], 
          font_properties=Helvetica(22), ha='left', va='top',
          transform=ax.transAxes)

