#!/usr/bin/env python

from numpy import *
from pygoods import *
import matplotlib as mpl
import matplotlib.pyplot as plt
#import pymc
import KPDFadaptnumpy as KPDF
from bivariate import fit_lbg as fl
from bivariate import bivRL as bl
from bivariate import mlutil
from scipy import optimize
from multiprocessing import Queue, Process, Pool
import time
import os, glob
import cPickle
from my_mplfonts import Helvetica
import pmscolors

pms = pmscolors.pmscolors()
#mpl.rcParams['ytick.labelsize'] = 14
#mpl.rcParams['xtick.major.size'] = 6
#mpl.rcParams['ytick.major.size'] = 6
#mpl.rcParams['axes.titlesize'] = 18
parnames = ['alpha','mstar','phistar','logr0','sigma','beta']
labels = [r'$\alpha$',r'$M^*$',r'$\phi^*(\times 10^{-3})$',r'$\log R_0$',
          r'$\sigma_{\ln R}$',r'$\beta$']
par_dic = {'alpha':0, 'mstar':1, 'phistar':2, 'logr0':3, 'sigma':4, 'beta':5}
## Stuff for B-dropouts
#parb = array([-1.70, -20.683, 1.462e-3, 0.836, 0.698, 0.235])
ticks_b = {'alpha':[-1.5,-1.7,-1.9],
           'mstar':[-21.,-20.4],
           'phistar':[1.0e-3,2.0e-3],
           'logr0':[0.75,0.85],
           'sigma':[0.8,0.9],
           'beta':[0.1,0.2,0.3]}
paramlims_b = {'alpha':[-2.0,-1.4],
             'mstar':[-21.2,-20.0],
             'phistar':[0.5e-3,3.2e-3],
             'logr0':[0.7,0.95],
             'sigma':[0.70,1.0],
             'beta':[0.0,0.4]}
parampix_b = {'alpha':0.01,
              'mstar':0.01,
              'phistar':5.e-5,
              'logr0':0.005,
              'sigma':0.005,
              'beta':0.005}
paramscale_b = {'alpha':0.5,
                'mstar':0.2,
                'phistar':200.,
                'logr0':1.0,
                'sigma':1.0,
                'beta':1.0}
htable_b = {(0,1):0.015, (0,2):0.020, (0,3):0.018, (0,4):0.018, (0,5):0.015,
            (1,2):0.015, (1,3):0.018, (1,4):0.018, (1,5):0.020,
            (2,3):0.023, (2,4):0.025, (2,5):0.025,
            (3,4):0.015, (3,5):0.018,
            (4,5):0.020}                
## Stuff for V-dropouts
#parv = array([-1.67, -20.49, 1.117e-3, 0.868, 0.830, 0.276])
ticks_v = {'alpha':arange(-2.3,-1.,0.6),
           'mstar':[-21.3,-20.5,-19.7],
           'phistar':[1.5e-3,3.0e-3],
           'logr0':[0.6,0.8,1.0],
           'sigma':[0.8,1.0,1.2],
           'beta':arange(-0.4,0.6,0.4)}
paramlims_v = {'alpha':[-2.4,-0.8],
               'mstar':[-21.5,-19.5],
               'phistar':[0.,4.5e-3],
               'logr0':[0.5,1.1],
               'sigma':[0.6,1.4],
               'beta':[-0.5,0.7]}
parampix_v = {'alpha':0.01,
              'mstar':0.01,
              'phistar':2.e-5,
              'logr0':0.005,
              'sigma':0.005,
              'beta':0.005}
paramscale_v = {'alpha':0.6,
                'mstar':0.5,
                'phistar':250.,
                'logr0':1.4,
                'sigma':1.0,
                'beta':1.0}
htable_v = {(0,1):0.045, (0,2):0.043, (0,3):0.047, (0,4):0.047, (0,5):0.050,
            (1,2):0.032, (1,3):0.053, (1,4):0.053, (1,5):0.050,
            (2,3):0.055, (2,4):0.058, (2,5):0.060,
            (3,4):0.040, (3,5):0.050,
            (4,5):0.060}
## Stuff for U-dropouts
# 130715_1
par_u = array([-1.882, -20.678, 2.946e-03, 0.727, 0.727, 0.225])
ticks_u = {'alpha':arange(-2.3,-1.5,0.3),
           'mstar':arange(-21.0, -20.0, 0.5),
           'phistar':arange(1.5,6.,1.5)*1.e-3,
           'logr0':arange(0.6,0.9,0.1),
           'sigma':arange(0.6,0.85,0.1),
           'beta':arange(0.,0.5,0.1)}
paramlims_u = {'alpha':[-2.4,-1.3],
               'mstar':[-21.5,-19.9],
               'phistar':[0.,6.5e-3],
               'logr0':[0.55,0.9],
               'sigma':[0.55,0.9],
               'beta':[0.,0.45]}
parampix_u = {'alpha':0.01,
              'mstar':0.01,
              'phistar':2.e-5,
              'logr0':0.005,
              'sigma':0.005,
              'beta':0.005}
paramscale_u = {'alpha':0.6,
                'mstar':0.5,
                'phistar':250.,
                'logr0':1.4,
                'sigma':1.0,
                'beta':1.0}
htable_u = {(0,1):0.030, (0,2):0.040, (0,3):0.030, (0,4):0.035, (0,5):0.025,
            (1,2):0.045, (1,3):0.030, (1,4):0.025, (1,5):0.030,
            (2,3):0.035, (2,4):0.030, (2,5):0.040,
            (3,4):0.015, (3,5):0.020,
            (4,5):0.020}
## Stuff for i-dropouts
# 130715_1
par_i = array([-1.842, -20.323, 3.406e-03, 0.675, 0.288, 0.181])
ticks_i = {'alpha':[-2.4,-1.7,-1.],
           'mstar':[-21.5,-20.5,-19.5],
           'phistar':[3.0e-3, 6.0e-3],
           'logr0':[0.6,0.7,0.8],
           'sigma':[0.2,0.4,0.6],
           'beta':arange(-0.4,0.6,0.4)}
paramlims_i = {'alpha':[-3.5,-0.0],
               'mstar':[-22.0,-19.0],
               'phistar':[0.,1.e-2],
               'logr0':[0.5,0.9],
               'sigma':[0.0,0.7],
               'beta':[-0.5,0.7]}
parampix_i = {'alpha':0.02,
              'mstar':0.01,
              'phistar':5.e-5,
              'logr0':0.005,
              'sigma':0.005,
              'beta':0.005}
paramscale_i = {'alpha':0.6,
                'mstar':0.5,
                'phistar':250.,
                'logr0':1.4,
                'sigma':1.0,
                'beta':1.0}
htable_i = {(0,1):0.070, (0,2):0.200, (0,3):0.080, (0,4):0.080, (0,5):0.120,
            (1,2):0.150, (1,3):0.060, (1,4):0.060, (1,5):0.100,
            (2,3):0.070, (2,4):0.065, (2,5):0.100,
            (3,4):0.030, (3,5):0.050,
            (4,5):0.060}

def scatter2pdf(xlim, ylim, dx, dy, xarr, yarr, sx, sy, h, aspect=1.0, 
                plot=False):
   # plot confidence-interval contours from the distribution of points xarr and yarr
   # xgrid, ygrid are the pixel grids
   # sx, sy are the scales of x and y, so that their numerical values are
   # similar for numerical stability
   numpts = len(xarr)
   xgrid = arange(xlim[0], xlim[1], dx) * sx
   ygrid = arange(ylim[0], ylim[1], dy) * sy
   edata = array([xarr*sx, yarr*sy])
   edata = edata.swapaxes(0,1)
   print "shape(edata)", shape(edata)
   #print edata[:3]
   #edata = edata.reshape(numpts, 2)
   coord = meshgrid(xgrid, ygrid)
   X = coord[0]; X = X.reshape(len(X.ravel()))
   Y = coord[1]; Y = Y.reshape(len(Y.ravel()))
   gdata = array([X, Y]); gdata = gdata.swapaxes(0, 1)
   print "shape(gdata)", shape(gdata)
   #print gdata[0], gdata[-1]
   #gdata = KPDF.MPDF2DGrid2Array(xgrid, ygrid)
   #gdata = gdata.reshape((len(gdata)/2, 2))
   pdf = KPDF.MPDFGaussian(edata, gdata, h)
   pdf = pdf.reshape((len(ygrid), len(xgrid)))
   if plot:
      fig1 = plt.figure()
      ax1 = fig1.add_subplot(111)
      #ratio = (xgrid[-1]-xgrid[0]) / (ygrid[-1]-ygrid[0])
      ax1.imshow(pdf, origin='lower', aspect=aspect)
      ax1.set_title('h=%f' % h)
      fig1.show()
      return pdf, ax1
   else:
      return pdf

def getpdf(x, pdf_x, att_arr):
   if x<=att_arr[0]:
      return pdf_x[0]
   elif x>=att_arr[-1]:
      return pdf_x[-1]
   else:
      y = searchsorted(att_arr, x)
      p0 = pdf_x[y-1]; p1 = pdf_x[y]
      att0 = att_arr[y-1]; att1 = att_arr[y]
      p = p0 + (x-att0)*(p1-p0)/(att1-att0)
      return p

def dpdf(x, pdf_x, att_arr, level):
	return getpdf(x, pdf_x, att_arr) - level

def calc_confidence_interval(c, par0, target_levels=[0.682,0.954]):
   # given the table c with the traces from MCMC, calculate the confidence intervals for
   # each parameter
   pdfs = []
   att_arrs = []
   attributes = ['alpha','mstar','phistar','logr0','sigma','beta']
   low_vals = zeros(6)
   high_vals = zeros(6)
   for i in range(6):
      att = attributes[i]
      trace = c.__getitem__(att)
      h = KPDF.UPDFOptimumBandwidth(trace)
      att_arr = arange(min(trace)-2*h, max(trace)+2*h, h/10.)
      att_arrs += [att_arr]
      pdf_i = KPDF.UPDFEpanechnikov(trace, att_arr, h)
      pdfs += [pdf_i]
      levels = conflevel(pdf_i, target_levels=target_levels)
      # now find the 1-sigma
      low_val = optimize.brentq(dpdf, min(att_arr), par0[i], 
                args=(pdf_i,att_arr,levels[0]))  # 1-sigma lower bound
      high_val = optimize.brentq(dpdf, par0[i], max(att_arr), 
                 args=(pdf_i,att_arr,levels[0]))  # 1-sigma upper bound
      low_vals[i] = low_val
      high_vals[i] = high_val
   for i in range(6):
      print attributes[i], (low_vals[i]-par0[i]), (high_vals[i]-par0[i]), par0[i]
   return low_vals, high_vals

def conflevel(pdf, target_levels=[0.682, 0.954]):
   # given a pdf with a well-defined center, find the level that encloses
   # the percentages of probability given in levels. The default is 68.2% (1-sigma)
   # and 95.4% (2-sigma)
   #pdf_norm = pdf / sum(pdf.ravel())  # normalize the pdf if it's not yet normalized
   #ptot = sum(pdf_norm.ravel())
   def dp(pdflevel, pdf2, targetlevel):
      pdf_above = compress(pdf2.ravel()>=pdflevel, pdf2.ravel())
      d = sum(pdf_above) - targetlevel*sum(pdf2.ravel())
      # the difference between the target level
      # and the probability above pdflevel
      return abs(d)
   levels = zeros(len(target_levels))
   for i in range(len(levels)):
      tl = target_levels[i]
      guess = max(pdf.ravel()) / 5.
      f = optimize.fmin(dp, guess, args=[pdf, tl])
      levels[i] = f[0]
   return levels
   
def confcontours(drop, parnamex, parnamey, xarr, yarr, xbest, ybest, 
                 tlevels=[0.682,0.954], ax=None,
                 colors=['red','blue'], ms=10, mew=2.0):
   if ax==None:
      fig = plt.figure()
      ax = fig.add_subplot(111)
   # make a contour plot of confidence intervals
   parlist = (par_dic[parnamex], par_dic[parnamey])
   if parlist[0] > parlist[1]:
      parlist = (parlist[1], parlist[0])
   if drop == 'b':
      h = htable_b[parlist] 
      xlim = paramlims_b[parnamex]
      ylim = paramlims_b[parnamey]
      dx = parampix_b[parnamex]
      dy = parampix_b[parnamey]
      scalex = paramscale_b[parnamex]
      scaley = paramscale_b[parnamey]
      ticks = ticks_b
   elif drop == 'v':
      h = htable_v[parlist]
      xlim = paramlims_v[parnamex]
      ylim = paramlims_v[parnamey]
      dx = parampix_v[parnamex]
      dy = parampix_v[parnamey]
      scalex = paramscale_v[parnamex]
      scaley = paramscale_v[parnamey]
      ticks = ticks_v
   elif drop == 'u':
      h = htable_u[parlist]
      xlim = paramlims_u[parnamex]
      ylim = paramlims_u[parnamey]
      dx = parampix_u[parnamex]
      dy = parampix_u[parnamey]
      scalex = paramscale_u[parnamex]
      scaley = paramscale_u[parnamey]
      ticks = ticks_u
   elif drop == 'i':
      h = htable_i[parlist]
      xlim = paramlims_i[parnamex]
      ylim = paramlims_i[parnamey]
      dx = parampix_i[parnamex]
      dy = parampix_i[parnamey]
      scalex = paramscale_i[parnamex]
      scaley = paramscale_i[parnamey]
      ticks = ticks_i
   pdf = scatter2pdf(xlim, ylim, dx, dy, xarr, yarr, scalex, scaley, h, 
                     plot=False)
   # calculate probability density 
   levels = conflevel(pdf, target_levels=tlevels)
   xgrid = arange(xlim[0], xlim[1], dx)
   ygrid = arange(ylim[0], ylim[1], dy)
   ax.contour(xgrid, ygrid, pdf, levels=levels, colors=colors, 
              linewidths=2)
   ax.set_xticks(ticks[parnamex])
   ax.set_xticklabels(ticks[parnamex])
   ax.set_yticks(ticks[parnamey])
   ax.set_yticklabels(ticks[parnamey])
   ax.plot([xbest], [ybest], '*', ms=ms, mew=mew, c='black')
   return ax

def trace2pdf(trace):
   # convert from trace to PDF
   h = KPDF.UPDFOptimumBandwidth(trace)
   xarr = arange(min(trace)-2*h, max(trace)+2*h, h/10.)
   pdf = KPDF.UPDFEpanechnikov(trace, xarr, h)
   return pdf, xarr

def conf_interval_grids(c, drop, bestpars, ms=10, 
                        colors=[pms.Bright_Red,pms.Bright_Blue],
                        scatter=True, figsize=(10,9)):
   # plot a grid of 2D confidence contours for each pair of parameters
   #mpl.rcParams['xtick.labelsize'] = 18
   #mpl.rcParams['ytick.labelsize'] = 18
   #mpl.rcParams['axes.labelsize'] = 20
   #mpl.rcParams['axes.formatter.limits'] = [-2, 7]
   if drop == 'b':
   	ticks = ticks_b
   	paramlims = paramlims_b
   elif drop == 'v':
   	ticks = ticks_v
   	paramlims = paramlims_v
   elif drop == 'u':
    ticks = ticks_u
    paramlims = paramlims_u
   elif drop == 'i':
    ticks = ticks_i
    paramlims = paramlims_i
   pandx = 0.8/6.  # x width of each panel
   pandy = 0.8/6.  # y width of each panel
   fig = plt.figure(figsize=figsize)
   axes = {}
   for i in arange(0,6): # x-index of the panel
      for j in arange(1,6)[::-1]:
         if i >= j: continue
         ax = fig.add_axes([0.1+i*pandx, 0.1+(5-j)*pandy, pandx, pandy])
         parx = parnames[i]; pary = parnames[j]
         print parx, pary
         if scatter == True:
            ax.scatter(c.__getitem__(parx), c.__getitem__(pary), marker='o', 
                       s=2**2, c=colors[0])
            ax.plot(bestpars[i], bestpars[j], '*', ms=12, c=colors[0], lw=2)
            ax.set_xticklabels(ticks[parx])
            ax.set_yticklabels(ticks[pary])
            ax.set_xlim(paramlims[parx])
            ax.set_ylim(paramlims[pary])
         else:
            ax = confcontours(drop, parx, pary, c.__getitem__(parx), 
               c.__getitem__(pary), bestpars[i], bestpars[j], ms=ms, 
               colors=colors, ax=ax)
         ax.set_xticklabels('')
         ax.set_yticklabels('')
         if i == 0:
            ax.set_ylabel(labels[j], font_properties=Helvetica(16))
            ax.set_yticks(ticks[pary])
            if j == 2:  # for phistar
               ax.set_yticklabels(array(ticks[pary])*1.e3, 
                                  font_properties=Helvetica(12))
            else:
               ax.set_yticklabels(ticks[pary], font_properties=Helvetica(12))
         if j == 5:
            ax.set_xlabel(labels[i], font_properties=Helvetica(16))
            ax.set_xticks(ticks[parx])
            if i == 2:  # for phistar
               ax.set_xticklabels(array(ticks[parx])*1.e3,
                                  font_properties=Helvetica(12))
            else:
               ax.set_xticklabels(ticks[parx],
                                  font_properties=Helvetica(12))
         axes[(i,j)] = ax
      parx = parnames[i]
      pdf, xarr = trace2pdf(c.__getitem__(parx))
      ax = fig.add_axes([0.1+i*pandx, 0.1+(5-i)*pandy, pandx, pandy])
      ax.plot(xarr, pdf, c='black', lw=2)
      axes[(i,6)] = ax
      ax.set_yticks([])
      ax.set_yticklabels([])
      ax.set_xticks(ticks[parx])
      ax.set_xlim(paramlims[parx])
      if i==5:
         ax.set_xticklabels(ticks['beta'], font_properties=Helvetica(12))
         ax.set_xlabel(r'$\beta$', font_properties=Helvetica(16))
         ax.set_xticks(ticks['beta'])
         ax.set_xlim(paramlims['beta'])
      else:
         ax.set_xticklabels([])
   fig.show()
   return axes, fig
   
#def hist_par(traces, burn=0, thin=1, limits=dflim, bs=dfbins,\
#   title='B-dropouts'):
#   #r_alpha, tr_Mstar, tr_logr0, tr_sigma, tr_beta = read_db(dbfiles)
#   plt.figure(figsize=(14,10))
#   histarg = {'ec':'black', 'lw':0.5}
#   traces = list(traces)
#   for i in range(len(traces)):
#      traces[i] = traces[i][::thin]
#   # alpha
#   plt.subplot(2,3,1)
#   plt.hist(traces[0],arange(limits[0][0], limits[0][1], bs[0]),**histarg)
#   plt.title('alpha')
#   plt.xticks(arange(limits[0][0], limits[0][1], 0.2))
#   plt.xlim(*limits[0])
#   #plt.yticks(size=10); plt.xticks(size=10)
#   # Mstar
#   plt.subplot(2,3,2)
#   plt.hist(traces[1],arange(limits[1][0], limits[1][1], bs[1]),**histarg)
#   plt.title('Mstar')
#   plt.xticks(arange(limits[1][0], limits[1][1], 0.2))
#   plt.xlim(*limits[1])
#   # logr0
#   plt.subplot(2,3,3)
#   plt.hist(traces[2],arange(limits[2][0], limits[2][1], bs[2]), **histarg)
#   plt.title('logr0')
#   plt.xticks(arange(limits[2][0], limits[2][1], 0.05))
#   plt.xlim(*limits[2])
#   # sigma
#   plt.subplot(2,3,4)
#   plt.hist(traces[3],arange(limits[3][0], limits[3][1], bs[3]), **histarg)
#   plt.title('sigma')
#   plt.xticks(arange(limits[3][0], limits[3][1], 0.05))
#   plt.xlim(*limits[3])
#   # beta
#   plt.subplot(2,3,5)
#   plt.hist(traces[4],arange(limits[4][0], limits[4][1], bs[4]), **histarg)
#   plt.title('beta')
#   plt.xticks(arange(limits[4][0], limits[4][1], 0.05))
#   plt.xlim(*limits[4])
#   plt.suptitle(title)
#
#   return 0


def trace_par(traces,burn=0,thin=1):
   plt.figure(figsize=(14,10))
   traces = list(traces)
   for i in range(len(traces)):
      traces[i] = traces[i][::thin]
   N = len(traces[0])
   # alpha
   plt.subplot(2,3,1)
   plt.plot(range(N),traces[0])
   plt.title('alpha', size=10)
   # Mstar
   plt.subplot(2,3,2)
   plt.plot(range(N), traces[1])
   plt.title('Mstar', size=10)
   # logr0
   plt.subplot(2,3,3)
   plt.plot(range(N), traces[2])
   plt.title('logr0', size=10)
   # sigma
   plt.subplot(2,3,4)
   plt.plot(range(N), traces[3])
   plt.title('sigma', size=10)
   # beta
   plt.subplot(2,3,5)
   plt.plot(range(N), traces[4])
   plt.title('beta', size=10)
   plt.suptitle('B-dropouts; iter=%d, burn=%d, thin=%d' % (len(traces[0]),burn,thin))
   return 0


def scatter_par(traces, par, thin=1, title=''):
   # scatter plots for 2 pairs of parameters: to visualize the confidence regions
   traces = list(traces)
   for i in range(len(traces)):
      traces[i] = traces[i][::thin]
   plt.figure(figsize=(12,8))
   plt.subplot(1,2,1)
   plt.scatter(traces[0], traces[1], s=3, color='#009B95', marker='o')
   plt.plot(par[0], par[1], marker='x', ls='None', markersize=14, mew=2.0, color='black')
   plt.xlabel('alpha'); plt.ylabel('Mstar')
   plt.subplot(1,2,2)
   plt.scatter(traces[2], traces[3], s=3, color='#BF6F30', marker='o')
   plt.plot(par[2], par[3], marker='x', ls='None', markersize=14, mew=2.0, color='black')
   plt.xlabel('log_R0'); plt.ylabel('sigma')
   if len(title):
      plt.suptitle(title, size=14)
   return 0


def scatter_allpar(traces, par, thin=1, title=''):
   traces = list(traces)
   for i in range(len(traces)):
      traces[i] = traces[i][::thin]
   parname = ['alpha', 'Mstar', 'logr0', 'sigma', 'beta']
   fig = plt.figure(figsize=(19, 13))
   nsub = 1
   for i in range(5):
      for j in range(i+1, 5):
         plt.subplot(3, 4, nsub); nsub += 1
         plt.scatter(traces[i], traces[j], marker='o', s=3, color='#009B95')
         plt.plot(par[i], par[j], marker='x', ls='None', markersize=14, mew=2.0, color='black')
         plt.xlabel(parname[i]); plt.ylabel(parname[j])
   if len(title):
      plt.suptitle(title, size=14)
   return 0


def KPDFparams(xtrace, ytrace, xlim, ylim, dx, dy, thin=1, h=None, estimator="Gaussian"):
   # use Kernel PDF estimator to estimate the probability distribution function
   # in a 2-parameter plane in order to derive confidence region
   # if xtrace and ytrace were run with thin = 1, we can apply thinning here
   xtrace = xtrace[::thin]
   ytrace = ytrace[::thin]
   traces = array([xtrace, ytrace])
   traces = traces.swapaxes(0,1)
   x = arange(xlim[0], xlim[1], dx)
   y = arange(ylim[0], ylim[1], dy)
   npix = array([len(x), len(y)])
   coord = meshgrid(x, y)
   X = coord[0]; X = X.reshape(len(X.ravel()))
   Y = coord[1]; Y = Y.reshape(len(Y.ravel()))
   g = array([X, Y]); g = g.swapaxes(0, 1)
   if h == None:
      h = KPDF.MPDFOptimumBandwidth(traces)
   print "h = ", h
   if estimator == 'Gaussian':
      print "Gaussian estimator"
      PDF = KPDF.MPDFGaussian(traces, g, h)
   else:
      print "Epanechnikov estimator"
      PDF = KPDF.MPDFEpanechnikov(traces, g, h)
   
   PDF = PDF.reshape(npix[::-1])
   PDF = PDF / max(PDF.ravel())
   return PDF, coord[0], coord[1], npix


def autocor(x, k=1):
   xavg = average(x)
   xvar = (x - xavg)**2
   xvar = sum(xvar)
   xcov = (x[:-k] - xavg) * (x[k:] - xavg)
   xcov = sum(xcov)
   rho_k = xcov / xvar
   return rho_k

