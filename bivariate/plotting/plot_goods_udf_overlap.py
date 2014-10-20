#!/usr/bin/env python

from numpy import *
from pygoods import *
import matplotlib as mpl
import matplotlib.pyplot as plt


# Plot the GALFIT-measured mag & Re for sources that overlap between GOODS and UDF


def plot_bdrops_overlap(c1, c2, ax, tol=0.2, fsize=14):
   # c1: from GOODS; c2: from UDF; tol: matching distance in arcsec
   cullflags = [0, 1, 3, 14, 19]
   mag_goods = []
   re_goods = []
   mag_udf = []
   re_udf = []
   c1i = []; c2i = []
   gflags = []  # 1 means goods fit, 0 means poor fit
   n = 0
   for i in range(len(c1)):
      if c1.delta_j2000[i] < 0.:
         angdist = angsep.angsep(c1.alpha_j2000[i], c1.delta_j2000[i], c2.alpha_j2000,
            c2.delta_j2000)
         if min(angdist)*3600. <= tol:
            j = argsort(angdist)[0]
            if c2.cullflag[j] in cullflags:
               n += 1
               c1i += [i]; c2i += [j]
            if (c2.reout_err[j]/c2.reout[j] <= 0.6) & (c2.chisqnu[j]<5.0):
               if (c1.reout_err[i]/c1.reout[i] <= 0.6) & (c1.chisqnu[i]<0.4):
                  gflags += [1]  
               else: gflags += [0]
               
               mag_goods += [c1.magout[i]]
               re_goods += [c1.reout[i]]
               mag_udf += [c2.magout[j]]
               re_udf += [c2.reout[j]]
   mag_goods = array(mag_goods); re_goods = array(re_goods)
   mag_udf = array(mag_udf); re_udf = array(re_udf)
   gflags = array(gflags)
   #fig = plt.figure()
   #ax = fig.add_subplot(111)
   #ax.loglog(re_goods, re_udf, '.')
   ax.plot(compress(gflags==1,mag_goods-mag_udf), compress(gflags==1,log10(re_goods/re_udf)), 
      'o', ms=8, c='black', label='B-dropouts')
   #ax.plot(compress(gflags==0,mag_goods-mag_udf), compress(gflags==0,log10(re_goods/re_udf)),
   #   'x', ms=8, c='red', label='poorly-fit in GOODS')
   ax.set_xlim(-.5, .5); ax.set_ylim(-0.5, 0.5)
   yticks = [-0.5, -0.3, -0.1, 0.1, 0.3, 0.5]
   ax.set_xticks([-0.5, 0., 0.5]); ax.set_yticks(yticks)
   ax.set_xticklabels([-0.5, 0., 0.5],size=fsize); ax.set_yticklabels(yticks,size=fsize)
   #ax.set_xticklabels(array([2.,5.,10.])*0.03); ax.set_yticklabels(array([2.,5.,10.])*0.03)
   ax.set_xlabel('MAG_GOODS - MAG_HUDF', size=(fsize+3))
   ax.set_ylabel('log10(RE_GOODS/RE_HUDF)', size=(fsize+3))
   ax.plot([0.], [0.], marker='+', mew=2.0, c='green', ms=20)
   ax.plot([-0.5, 0.5], [0.1, 0.1], '--', c='black')
   ax.plot([-0.5, 0.5], [-0.1, -0.1], '--', c='black')
   print "A total of %d B-dropouts overlap between GOODS and UDF." % n
   return c1i, c2i, mag_goods, re_goods, mag_udf, re_udf


def plot_vdrops_overlap(c1, c2, ax, tol=0.2, fsize=14):
   # c1: from GOODS; c2: from UDF; tol: matching distance in arcsec
   cullflags = [0, 1, 3, 14, 19]
   mag_goods = []
   re_goods = []
   mag_udf = []
   re_udf = []
   c1i = []; c2i = []
   gflags = []  # 1 means goods fit, 0 means poor fit
   n = 0
   for i in range(len(c1)):
      if c1.delta_j2000[i] < 0.:
         angdist = angsep.angsep(c1.alpha_j2000[i], c1.delta_j2000[i], c2.alpha_j2000,
            c2.delta_j2000)
         if min(angdist)*3600. <= tol:
            print c1.cullflag[i]
            j = argsort(angdist)[0]
            if c2.cullflag[j] in cullflags:
               n += 1
               c1i += [i]; c2i += [j]
            if (c2.reout_err[j]/c2.reout[j] <= 0.6) & (c2.chisqnu[j]<5.0):
               if (c1.reout_err[i]/c1.reout[i] <= 0.6) & (c1.chisqnu[i]<0.5):
                  gflags += [1]  
               else: gflags += [0]
               
               mag_goods += [c1.magout[i]]
               re_goods += [c1.reout[i]]
               mag_udf += [c2.magout[j]]
               re_udf += [c2.reout[j]]
   mag_goods = array(mag_goods); re_goods = array(re_goods)
   mag_udf = array(mag_udf); re_udf = array(re_udf)
   gflags = array(gflags)
   #fig = plt.figure()
   #ax = fig.add_subplot(111)
   #ax.loglog(re_goods, re_udf, '.')
   ax.plot(compress(gflags==1,mag_goods-mag_udf), compress(gflags==1,log10(re_goods/re_udf)), 
      '^', ms=8, c='blue', label='V-dropouts')
   #ax.plot(compress(gflags==0,mag_goods-mag_udf), compress(gflags==0,log10(re_goods/re_udf)),
   #   'x', ms=8, c='red', label='poorly-fit in GOODS')
   ax.set_xlim(-.5, .5); ax.set_ylim(-0.5, 0.5)
   yticks = [-0.5, -0.3, -0.1, 0.1, 0.3, 0.5]
   ax.set_xticks([-0.5, 0., 0.5]); ax.set_yticks(yticks)
   ax.set_xticklabels([-0.5, 0., 0.5], size=fsize); ax.set_yticklabels(yticks, size=fsize)
   #ax.set_xticklabels(array([2.,5.,10.])*0.03); ax.set_yticklabels(array([2.,5.,10.])*0.03)
   ax.set_xlabel('MAG_GOODS - MAG_HUDF', size=(fsize+3))
   ax.set_ylabel('log10(RE_GOODS/RE_HUDF)', size=(fsize+3))
   ax.plot([0.], [0.], marker='+', mew=2.0, c='green', ms=20)
   ax.plot([-0.5, 0.5], [0.1, 0.1], '--', c='black')
   ax.plot([-0.5, 0.5], [-0.1, -0.1], '--', c='black')
   print "A total of %d V-dropouts overlap between GOODS and UDF." % n
   return c1i, c2i, mag_goods, re_goods, mag_udf, re_udf


def plot_all_overlap(cb1, cb2, cv1, cv2, tol=0.2, fsize=14):
   fig = plt.figure()
   ax = fig.add_subplot(111)
   x1 = plot_bdrops_overlap(cb1, cb2, ax, tol=tol, fsize=fsize)
   x2 = plot_vdrops_overlap(cv1, cv2, ax, tol=tol, fsize=fsize)
   # override the axis settings
   ax.set_xlim(-0.4, 0.4)
   ax.set_xticks([-0.4, 0., 0.4])
   ax.set_xticklabels([-0.4, 0., 0.4], size=fsize)
   ax.legend(loc=1, numpoints=1)
   return ax


def plot_diff_overlap(cb1, cb2, cv1, cv2, colorb='blue', colorv='red'):
   fig = plt.figure(figsize=(12,6))
   mpl.rcParams['lines.markersize'] = 15
   mpl.rcParams['axes.labelsize'] = 20
   mpl.rcParams['xtick.labelsize'] = 18
   mpl.rcParams['ytick.labelsize'] = 18
   n_bdrops_overlap = sum((cb1.cullflag==20)&(cb1.magout>0.))
   index_bdrops_overlap = compress((cb1.cullflag==20)&(cb1.magout>0.), arange(len(cb1.magout)))
   reout_bdrops_g = zeros(n_bdrops_overlap)
   reout_bdrops_u = zeros(n_bdrops_overlap)
   magout_bdrops_g = zeros(n_bdrops_overlap)
   magout_bdrops_u = zeros(n_bdrops_overlap)
   zmagauto_bdrops_u = zeros(n_bdrops_overlap)

   n_vdrops_overlap = sum((cv1.cullflag==20)&(cv1.magout>0.))
   index_vdrops_overlap = compress((cv1.cullflag==20)&(cv1.magout>0.), arange(len(cv1.magout)))
   reout_vdrops_g = zeros(n_vdrops_overlap)
   reout_vdrops_u = zeros(n_vdrops_overlap)
   magout_vdrops_g = zeros(n_vdrops_overlap)
   magout_vdrops_u = zeros(n_vdrops_overlap)
   zmagauto_vdrops_u = zeros(n_vdrops_overlap)
   
   for i in range(len(index_bdrops_overlap)):
      ii = index_bdrops_overlap[i]
      angdist = angsep.angsep(cb1.alpha_j2000[ii], cb1.delta_j2000[ii], cb2.alpha_j2000,
         cb2.delta_j2000)
      j = argsort(angdist)[0]
      reout_bdrops_g[i] = cb1.reout[ii]
      reout_bdrops_u[i] = cb2.reout[j]
      magout_bdrops_g[i] = cb1.magout[ii]
      magout_bdrops_u[i] = cb2.magout[j]
      zmagauto_bdrops_u[i] = cb2.z_mag_auto[j]

   for i in range(len(index_vdrops_overlap)):
      ii = index_vdrops_overlap[i]
      angdist = angsep.angsep(cv1.alpha_j2000[ii], cv1.delta_j2000[ii], cv2.alpha_j2000,
         cv2.delta_j2000)
      j = argsort(angdist)[0]
      reout_vdrops_g[i] = cv1.reout[ii]
      reout_vdrops_u[i] = cv2.reout[j]
      magout_vdrops_g[i] = cv1.magout[ii]
      magout_vdrops_u[i] = cv2.magout[j]
      zmagauto_vdrops_u[i] = cv2.z_mag_auto[j]

   #print magout_bdrops_g, magout_bdrops_u
   #print magout_vdrops_g, magout_vdrops_u
   #sp1 = fig.add_subplot(121)
   sp1 = fig.add_axes([0.1,0.1,0.35,0.8])
   sp1.plot(zmagauto_bdrops_u, log10(reout_bdrops_g/reout_bdrops_u), '.', 
      c=colorb,label='B-drops')
   sp1.plot(zmagauto_vdrops_u, log10(reout_vdrops_g/reout_vdrops_u), 'x',
      c=colorv, mew=2.0, label='V-drops')
   #sp2 = fig.add_subplot(122)
   sp2 = fig.add_axes([0.6,0.1,0.35,0.8])
   sp2.plot(zmagauto_bdrops_u, magout_bdrops_g - magout_bdrops_u, '.', 
      c=colorb, label='B-drops')
   sp2.plot(zmagauto_vdrops_u, magout_vdrops_g - magout_vdrops_u, 'x',
      c=colorv, mew=2.0, label='V-drops')
   sp1.set_xlabel('z-band MAG_AUTO in HUDF')
   sp1.set_ylabel('log10(Re_GOODS/Re_HUDF) [pixel]')
   sp1.set_xlim(23.0,27.5)
   sp2.set_xlim(23.0,27.5)
   xtick1 = sp1.get_xticks()
   sp1.set_xticks(arange(xtick1[1],xtick1[-1],1.0))
   sp1.set_xticklabels(arange(xtick1[1],xtick1[-1],1.0))
   ytick1 = sp1.get_yticks()
   sp1.set_yticks(arange(ytick1[1],ytick1[-1],0.1))
   sp1.set_yticklabels(arange(ytick1[1],ytick1[-1],0.1))
   sp2.set_xlabel('z-band MAG_AUTO in HUDF')
   sp2.set_ylabel('mag_GOODS - mag_HUDF')
   sp2.set_xticks(arange(xtick1[1],xtick1[-1],1.0))
   sp2.set_xticklabels(arange(xtick1[1],xtick1[-1],1.0))
   #sp2.set_yticklabels(sp2.get_yticks())
   sp2.set_ylim([-1.0, 1.0])
   sp1.legend(loc=2,numpoints=1)
   sp2.legend(loc=2,numpoints=1)
   
   return fig


