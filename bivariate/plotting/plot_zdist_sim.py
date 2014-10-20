#!/usr/bin/env python

from numpy import *
from pygoods import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import dropout_sel_sim as dss
import KPDFadaptnumpy as KPDF

fp11 = mpl.font_manager.FontProperties(size=11)
def fp(size):
   return mpl.font_manager.FontProperties(size=size)
   

def plot_bdrops_zdist(c):
   # plot the redshift distribution of selected B-drops in simulation
   # plot in 3 different M_1500 bins
   Marr = [-22., -21., -20.]; dM = 0.5
   hatch = ['/','\\','+']
   colors = ['blue','green','red']
   dropcrit = dss.bdrops_sel(c, zmagcut=26.5)[0]
   plt.figure()
   zarr = arange(3.0, 5.05, 0.05)
   for i in range(len(Marr)):
      withinbin = (c.m1500_input>=Marr[i]) & (c.m1500_input<(Marr[i]+dM))
      zsel = compress(withinbin & dropcrit, c.z_input)
      ztar = compress(withinbin & c.detect, c.z_input)
      print len(zsel), len(ztar)
      zselh = histogram(zsel, zarr)[0]
      ztarh = histogram(ztar, zarr)[0]
      zh = zselh.astype('float') / ztarh.astype('float')
      zh = zh / sum(zh)
      #h = KPDF.UPDFOptimumBandwidth(zh); print h
      #pz = KPDF.UPDFEpanechnikov(zh, zarr , h)
      #plt.plot(zarr, pz, label='%.1f<M1500<%.1f' % (Marr[i],(Marr[i]+dM)))
      plt.plot(zarr[:-1], zh, label='%.1f<M1500<%.1f' % (Marr[i],(Marr[i]+dM)))
   plt.xlabel('redshift')
   plt.ylabel('P(z)')
   plt.legend(loc='upper right', prop=fp)
   return 0

def show_zdist(zdgrid, zd_index, lw=1.5, ax=None):
   if ax==None:
      plt.clf()
      ax = plt.subplot(111)
   bbox_props = dict(boxstyle='round',fc='beige',ec='black',lw=1.2)
   zd = zdgrid.zdist[zd_index]
   plt.plot(zdgrid.zarr,zd.Pz,ls='steps-post',lw=lw, label='[%.1f:%.1f,%.1f:%.1f]'%(zd.M0,zd.M1,zd.logR0,zd.logR1))
   plt.xlabel('Input redshift',size=14)
   plt.ylabel('P(z)',size=14)
   plt.title("Index: %s" % str(zd_index))
   plt.text(0.05,0.5,"N_input=%d\nN_detect=%d\nN_select=%d"%(sum(zd.ninput),sum(zd.ndetect),sum(zd.nselect)),
      transform=ax.transAxes,bbox=bbox_props)
   plt.legend(loc=2)
   plt.ylim(0.,1.)

def plot_zdgrid(zdgrid, rindex=[6], Mindex=[4,6,8,10], ls=['steps','steps--','steps:','steps'],\
   lw=[2.0,1.5,1.5,1.0], color='black', title="", labelatt='M0', newfig=True, ax=None,
   legend_loc=1):
   # plot P(M, R, z) for different bands
   if newfig:
      fig = plt.figure()
      ax = fig.add_subplot(111)
   k = 0
   if len(Mindex) > 1:
      for i in range(len(Mindex)):
         j = 0
         Mi = Mindex[i]
         ri = rindex[j]
         zd = zdgrid.zdist[(Mi, ri)]
         print "M0=", zd.M0, "logR_{low}=", zd.logR0
         if labelatt == 'M0':
            label = r'$M_{low}$=%.1f' % (zd.M0)
         else:
            label = '%s=%.2f' % (labelatt, (zd.logR0))
         ax.plot(zdgrid.zarr, zd.Pz, ls=ls[k], lw=lw[k], color=color, label=label)
         k = k+1
   else:
      i = 0
      for j in range(len(rindex)):
         Mi = Mindex[i]
         ri = rindex[j]
         zd = zdgrid.zdist[(Mi, ri)]
         print "M0=", zd.M0, "logR_{low}=", zd.logR0
         if labelatt == 'M0':
            label = r'$\log_{10}R_{low}$=%.1f' % (zd.logR0)
         else:
            label = '%s=%.2f' % (labelatt, (zd.logR0))
         ax.plot(zdgrid.zarr, zd.Pz, ls=ls[k], lw=lw[k], color=color, label=label)
         k = k+1
   ax.legend(loc=legend_loc)
   ax.set_xlabel('Redshift')
   ax.set_ylabel(r'$P(M, R_e, z)$')
   ax.set_title(title)
   return ax


def plot_zdgrid_bdrops(zdgrid1, zdgrid1_nolya, fsize=12):
   # plot P(M, R, z) for B-dropouts with Lya [top row]
   # and w/o Lya [bottom row] in both magnitude bins [left column] and 
   # size bins [right column]
   fig = plt.figure(figsize=(12,8))
   axes = []
   
   # B-dropouts
   # B-dropouts in M bins w/ Lya
   ax1 = fig.add_axes([0.1, 0.5, 0.4, 0.4])
   ax1 = plot_zdgrid(zdgrid1, rindex=[6], Mindex=[4,6,8,10], labelatt='M0', newfig=False, 
      ax=ax1)
   axes += [ax1]
   # B-dropouts in logR bins w/ Lya
   ax2 = fig.add_axes([0.5, 0.5, 0.4, 0.4])
   ax2 = plot_zdgrid(zdgrid1, rindex=[4,6,8,10], Mindex=[8], labelatt=r'$\log_{10} R_{low}$', 
      newfig=False, ax=ax2)
   axes += [ax2]
   
   # SECOND ROW
   # B-dropouts in M bins w/o Lya
   ax3 = fig.add_axes([0.1, 0.1, 0.4, 0.4])
   ax3 = plot_zdgrid(zdgrid1_nolya, rindex=[6], Mindex=[4,6,8,10], labelatt='M0', newfig=False, 
      ax=ax3)
   axes += [ax3]
   # B-dropouts in logR bins w/o Lya
   ax4 = fig.add_axes([0.5, 0.1, 0.4, 0.4])
   ax4 = plot_zdgrid(zdgrid1_nolya, rindex=[4,6,8,10], Mindex=[8], labelatt=r'$\log_{10} R_{low}$', 
      newfig=False, ax=ax4)
   axes += [ax4]

   # formatting
   #axes[5].set_xticklabels(['',3.5,4.0,4.5,5.0,5.5,6.0])
   for i in range(len(axes)):
      #if i in [2, 3, 6, 7]:
      #   axes[i].legend(loc=2, prop=fp11)
      axes[i].set_xlim(3.0, 6.0)
      if i in [1, 3]:
         axes[i].set_yticklabels([])
         axes[i].set_ylabel("")
      if i in [0, 1]:
         axes[i].set_xticklabels([])
         axes[i].set_xlabel("")
      if i in [2, 3]:
         axes[i].set_xticklabels([3.,3.5,4.,4.5,5.,5.5,""], size=fsize)
      #if i in [2, 3]:
      #   axes[i].set_xticklabels(["",3.5,4.,4.5,5.,5.5,6.])
      if i in [0, 2]:
         axes[i].set_yticklabels(arange(0., 1., 0.2), size=fsize)
      #if i in [0, 1]:
      #   axes[i].set_title('B-drops')
      #if i in [2, 3]:
      #   axes[i].set_title('V-drops')
      
   

   # add dropout labels
   #fig.text(0.05, 0.85, r"B-dropouts with Ly$\alpha$", rotation='vertical', size=14)
   #fig.text(0.05, 0.55, r"B-dropouts without Ly$\alpha$", rotation='vertical', size=14)
   #fig.text(0.05, 0.25, r"V-dropouts with Ly$\alpha$", rotation='vertical', size=14)
   fig.text(0.03, 0.73, r"with Ly$\alpha$", rotation="vertical", size=14)
   fig.text(0.03, 0.35, r"without Ly$\alpha$", rotation="vertical", size=14)
   return axes

def plot_zdgrid_vdrops(zdgrid, zdgrid_nolya, fsize=12):
   # plot P(M, R, z) for V-dropouts with Lya [top row]
   # and w/o Lya [bottom row] in both magnitude bins [left column] and 
   # size bins [right column]
   fig = plt.figure(figsize=(12,8))
   axes = []
   ax1 = fig.add_axes([0.1, 0.5, 0.4, 0.4])
   ax1 = plot_zdgrid(zdgrid, rindex=[6], Mindex=[4,6,8,10], labelatt='M0', newfig=False, 
      ax=ax1)
   axes += [ax1]
   # V-dropouts in logR bins w/ Lya
   ax2 = fig.add_axes([0.5, 0.5, 0.4, 0.4])
   ax2 = plot_zdgrid(zdgrid, rindex=[4,6,8,10], Mindex=[8], labelatt=r'$\log_{10} R_{low}$', 
      newfig=False, ax=ax2)
   axes += [ax2]
   
   # V-dropouts in M bins w/o Lya
   ax3 = fig.add_axes([0.1, 0.1, 0.4, 0.4])
   ax3 = plot_zdgrid(zdgrid_nolya, rindex=[6], Mindex=[4,6,8,10], labelatt='M0', newfig=False, 
      ax=ax3)
   axes += [ax3]
   # V-dropouts in logR bins w/o Lya
   ax4 = fig.add_axes([0.5, 0.1, 0.4, 0.4])
   ax4 = plot_zdgrid(zdgrid_nolya, rindex=[4,6,8,10], Mindex=[8], labelatt=r'$\log_{10} R_{low}$', 
      newfig=False, ax=ax4)
   axes += [ax4]

   for i in range(len(axes)):
      axes[i].set_xlim(3.0, 6.0)
      axes[i].legend(loc=2, prop=fp(12))
      if i in [1, 3]:
         axes[i].set_yticklabels([])
         axes[i].set_ylabel("")
      if i in [0, 1]:
         axes[i].set_xticklabels([])
         axes[i].set_xlabel("")
      if i in [2, 3]:
         axes[i].set_xticklabels([3.,3.5,4.,4.5,5.,5.5,""], size=fsize)
      if i in [0, 2]:
         axes[i].set_yticklabels(arange(0., 1., 0.2), size=fsize)
   fig.text(0.03, 0.73, r"with Ly$\alpha$", rotation="vertical", size=(fsize+3))
   fig.text(0.03, 0.35, r"without Ly$\alpha$", rotation="vertical", size=(fsize+3))
   return axes

def plot_zdgrid_2(zdgrid):
   # plot two panels, one for P(M, R, z) in the same M bin, the other for P(M, R, z)
   # in the same R bin
   fig = plt.figure(figsize=(10, 5))
   ax1 = fig.add_axes([0.1, 0.1, 0.4, 0.8])
   ax2 = fig.add_axes([0.5, 0.1, 0.4, 0.8])
   plot_zdgrid(zdgrid, rindex=[6], Mindex=[4,6,8,10], labelatt='M0', newfig=False,\
      ax=ax1)
   plot_zdgrid(zdgrid, rindex=[4,6,8,10], Mindex=[8], labelatt=r'$\log_{10} R_{low}$',\
      newfig=False, ax=ax2)
   
   axes = [ax1, ax2]
   ax1.set_title(r'$\log_{10}(R_e) \in [%.2f, %.2f]$ [arcsec]' % \
      ((zdgrid.zdist[(6,6)].logR0+log10(0.03),\
      (zdgrid.zdist[(6,6)].logR1+log10(0.03)))))
   ax2.set_title(r'$M \in [%.1f, %.1f]$' % (zdgrid.zdist[(8,6)].M0, zdgrid.zdist[(9,6)].M0))
   ax2.set_yticklabels([])
   ax2.set_ylabel('')

   return axes

def plot_zdgrid_all(zdgrid_b, zdgrid_v):
   # plot P(M, R, z) for both B-dropouts [top row] and V-dropouts [bottom row]
   # in both magnitude bins [left column] and size bins [right column]
   # all without Lya
   fig = plt.figure(figsize=(12,8))
   axes = []
   mpl.rcParams['axes.labelsize'] = 20
   mpl.rcParams['xtick.labelsize'] = 18
   mpl.rcParams['ytick.labelsize'] = 18
   # B-dropouts
   # B-dropouts in M bins
   ax1 = fig.add_axes([0.1, 0.5, 0.4, 0.4])
   ax1 = plot_zdgrid(zdgrid_b, rindex=[6], Mindex=[4,6,8,10], labelatt='M0', newfig=False, 
      ax=ax1,legend_loc=1)
   axes += [ax1]
   # B-dropouts in logR bins
   ax2 = fig.add_axes([0.5, 0.5, 0.4, 0.4])
   ax2 = plot_zdgrid(zdgrid_b, rindex=[4,6,8,10], Mindex=[8], labelatt=r'$\log_{10} R_{low}$', 
      newfig=False, ax=ax2, legend_loc=1)
   axes += [ax2]
   
   # SECOND ROW
   # V-dropouts in M bins
   ax3 = fig.add_axes([0.1, 0.1, 0.4, 0.4])
   ax3 = plot_zdgrid(zdgrid_v, rindex=[6], Mindex=[4,6,8,10], labelatt='M0', newfig=False, 
      ax=ax3, legend_loc=2)
   axes += [ax3]
   # V-dropouts in logR bins
   ax4 = fig.add_axes([0.5, 0.1, 0.4, 0.4])
   ax4 = plot_zdgrid(zdgrid_v, rindex=[4,6,8,10], Mindex=[8], labelatt=r'$\log_{10} R_{low}$', 
      newfig=False, ax=ax4, legend_loc=2)
   axes += [ax4]

   # formatting
   #axes[5].set_xticklabels(['',3.5,4.0,4.5,5.0,5.5,6.0])
   for i in range(len(axes)):
      axes[i].set_xlim(3.0, 6.0)
      if i in [1, 3]:
         axes[i].set_yticklabels([])
         axes[i].set_ylabel("")
      if i in [0, 1]:
         axes[i].set_xticklabels([])
         axes[i].set_xlabel("")
      if i in [2, 3]:
         axes[i].set_xticklabels([3.,3.5,4.,4.5,5.,5.5,""])
      #if i in [2, 3]:
      #   axes[i].set_xticklabels(["",3.5,4.,4.5,5.,5.5,6.])
      if i == 0:
         axes[i].set_yticklabels(arange(0., 1., 0.2))
      if i == 2:
         axes[i].set_yticklabels(['']+list(arange(0.2,1.,0.2)))
      #if i in [0, 1]:
      #   axes[i].set_title('B-drops')
      #if i in [2, 3]:
      #   axes[i].set_title('V-drops')
      
   

   # add dropout labels
   #fig.text(0.05, 0.85, r"B-dropouts with Ly$\alpha$", rotation='vertical', size=14)
   #fig.text(0.05, 0.55, r"B-dropouts without Ly$\alpha$", rotation='vertical', size=14)
   #fig.text(0.05, 0.25, r"V-dropouts with Ly$\alpha$", rotation='vertical', size=14)
   fig.text(0.01, 0.73, r"B-dropouts", rotation="vertical", size=20)
   fig.text(0.01, 0.35, r"V-dropouts", rotation="vertical", size=20)
   return axes