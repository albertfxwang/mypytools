#!/usr/bin/env python

from numpy import *
from pygoods import *
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import bivariate_lf as bl
import fit_lbg as fl
#from fit_lbg import cleandata
import bivariate_fit as bf
import schechter
import mlutil
import cPickle
import zdist
from scipy import optimize
import fit_lognormal as fln

mci = bl.mconvert('M1500_to_i.txt')
mcz = bl.mconvert('M1500_to_z.txt')
limitsb1 = array([[22.5, 26.5], [-0.8, 1.8]])
limitsb2 = array([[24.5, 28.5], [-0.8, 1.8]])
limitsv1 = array([[23.0, 26.5], [-0.8, 1.8]])
limitsv2 = array([[24.5, 28.5], [-0.8, 1.8]])

#def bdrops_all():
#   mag1, re1, crit1 = fl.cleandata('bdrops_gf_v2.cat')
#   mag2, re2, crit2 = fl.cleandata('bdrops_udf_gf_v2.cat', chisqnulim=5.0, magautolim=28.5,\
#      limits=limits2)
#   plt.figure()
#   plt.scatter(mag1,log10(re1),marker='o',color='blue',s=6,label='GOODS ACS')
#   plt.scatter(mag2,log10(re2),marker='x',color='red',s=8,label='UDF')
#   plt.legend(loc='lower left')
#   plt.xlabel('i-band MAG_GALFIT')
#   plt.ylabel('i-band log10(Re) [pixels]')
#   plt.title('B-dropouts')
#   return 0


#def vdrops_all():
#   mag1, re1, crit1 = fl.cleandata('vdrops_gf.cat')
#   mag2, re2, crit2 = fl.cleandata('vdrops_udf_gf.cat', chisqnulim=5.0, magautolim=28.5,\
#      limits=limits2)
#   plt.figure()
#   plt.scatter(mag1,log10(re1),marker='o',color='blue',s=6,label='GOODS ACS')
#   plt.scatter(mag2,log10(re2),marker='x',color='red',s=8,label='UDF')
#   plt.legend(loc='lower left')
#   plt.xlabel('z-band MAG_GALFIT')
#   plt.ylabel('z-band log10(Re) [pixels]')
#   plt.title('V-dropouts')
#   return 0


def dropout_mag_logre_modgrid(mag,re,ax,limits,color='black'):
   #crit = (c.chisqnu<=0.4) & (c.reout_err/c.reout<=0.5)
   #crit = crit & ((c.cullflag==0)|(c.cullflag==1)|(c.cullflag==3)|(c.cullflag==14))
   #mag = compress(crit,c.magout)
   #logre = compress(crit,log10(c.reout))
   pixdx = array([0.02,0.02])
   logre = log10(re)
   x = (mag-limits[0][0])/pixdx[0]
   y = (logre-limits[1][0])/pixdx[1]
   ax.scatter(x,y,s=4,color=color)
   xmax = (limits[0][1]-limits[0][0])/pixdx[0]
   ymax = (limits[1][1]-limits[1][0])/pixdx[1]
   ax.set_xlim(0,xmax); ax.set_ylim(0,ymax)

   return ax


def findline(mag, logre):
   def linedist(a, b, x0, y0):
      d2 = (y0 - a*x0 - b)**2 / (a**2 + 1.)
      return d2
   def totaldist2(pl, m, r):
      d2 = linedist(pl[0], pl[1], m, r)
      d2 = sum(d2)
      return d2
   y = optimize.fmin(totaldist2, [-0.3, 0.5], args=[mag, logre])
   return y


#def observed_bdrops_beta(cat1, cat2):
#   # find the observed peak size-luminosity trend without correcting for incompleteness
#   mag1, re1, crit1 = fl.cleandata(cat1, reerr_ratiolim=0.6, chisqnulim=0.4, magautolim=26.5,\
#      zlo=3.0, drop='b', limits=limits1)
#   mag2, re2, crit2 = fl.cleandata(cat2, reerr_ratiolim=0.6, chisqnulim=5.0, magautolim=28.5,\
#      zlo=3.0, drop='b', limits=limits2)
#   mag = concatenate((mag1, mag2))
#   logre = log10(concatenate((re1, re2)))
#   y = findline(mag, logre)
#   return y


#def observed_vdrops_beta(cat1, cat2):
#   # find the observed peak size-luminosity trend without correcting for incompleteness
#   mag1, re1, crit1 = fl.cleandata(cat1, reerr_ratiolim=0.6, chisqnulim=0.5, magautolim=26.5,\
#      zlo=4.0, drop='v', limits=limits1)
#   mag2, re2, crit2 = fl.cleandata(cat2, reerr_ratiolim=0.6, chisqnulim=5.0, magautolim=28.5,\
#      zlo=4.0, drop='v', limits=limits2)
#   mag = concatenate((mag1, mag2))
#   logre = log10(concatenate((re1, re2)))
#   y = findline(mag, logre)
#   return y


def show_model(par, drop, field, newfig=True, axCent=None, fig1=None, lfbw=0.2, sdbw=0.2,
   colors=['blue','green','black','red']):
   if drop == 'b':
      zmean = 4.0
      mc = bl.mconvert('M1500_to_i.txt')
      if field == 'goods':
         cat = 'bdrops_gf_v3.cat'
         zdgrid = zdist.read_zdgrid('zdgrid/zdgrid_bdrops_nolya.p')
         limits = limitsb1
         magautolim=26.5; chisqnulim=0.4; reerrlim=0.6
         kgridfile = 'tfkernel/kernel_I.p'
         dataset = 'GOODS'
      elif field == 'udf':
         cat = 'bdrops_udf_gf_v3.cat'
         zdgrid = zdist.read_zdgrid('zdgrid/zdgrid_bdrops_udf_nolya.p')
         limits = limitsb2
         magautolim=28.5; chisqnulim=5.0; reerrlim=0.6
         kgridfile = 'tfkernel/kernel_I_udf.p'
         dataset = 'HUDF'
   elif drop == 'v':
      zmean = 5.0
      mc = bl.mconvert('M1500_to_z.txt')
      if field == 'goods':
         cat = 'vdrops_gf_v2.cat'
         zdgrid = zdist.read_zdgrid('zdgrid/zdgrid_vdrops_nolya_bston5.p')
         limits = limitsv1
         magautolim=26.5; chisqnulim=0.5; reerrlim=0.6
         kgridfile = 'tfkernel/kernel_Z.p'
         dataset = 'GOODS'
      elif field == 'udf':
         cat = 'vdrops_udf_gf_v2.cat'
         zdgrid = zdist.read_zdgrid('zdgrid/zdgrid_vdrops_udf_nolya_bston5.p')
         limits = limitsv2
         magautolim=28.5; chisqnulim=5.0; reerrlim=0.6
         kgridfile = 'tfkernel/kernel_Z_udf.p'
         dataset = 'HUDF'
         
   kgrid = mlutil.readkgrid(kgridfile)
   
   meankcorr = mc(zmean)
   fp = matplotlib.font_manager.FontProperties(size=9)
   if newfig:
      fig1 = plt.figure(figsize = (10, 15))
      axCent = plt.subplot(111)
   divider = make_axes_locatable(axCent)
   
   restlim = array([[-26.5, -15.5], [-0.5, 2.0]])
   pixdx = array([0.02, 0.02])

   if drop=='b': z0=3.0
   else: z0=4.0
   zd_flat = zdist.zdgrid(-25.0,-15.0,0.5,-0.6,1.8,0.2,z0,6.0,0.1,drop,zdgrid.area)
   zd_flat.flat_zdgrid(zlo=zmean-0.5, zhi=zmean+0.5)
   f = open('zdgrid_flat.p','w')
   cPickle.dump(zd_flat, f, 2)
   f.close()

   model0 = bl.bivariate_lf(par, limits, pixdx, drop, field, mc=mc, meankcorr=meankcorr,\
      zdgrid=None)
   V0 = zdgrid.dVdz[(zdgrid.zarr>=(zmean-0.5))&(zdgrid.zarr<=(zmean+0.5))]
   model0.model = model0.model * sum(V0)

   model_ds = bl.bivariate_lf(par, limits, pixdx, drop, field, kgrid=None, meankcorr=meankcorr,\
      zdgrid=zdgrid, mc=mc, norm=-1)
   model = bl.bivariate_lf(par, limits, pixdx, drop, field, kgrid=kgrid,
      meankcorr=meankcorr, M0=-21.0, zdgrid=zdgrid, mc=mc, norm=-1)
   # Show model + data points 
   #axCent.imshow(model.model.swapaxes(0,1), origin = 'lower', vmin = vmin, vmax = vmax,
   #   aspect = 'auto')
   mag,re,crit = fl.cleandata(cat, chisqnulim = chisqnulim,
      reerr_ratiolim = reerrlim, limits=limits, drop=drop)
   npts = len(mag); print npts
   axCent.scatter(mag,log10(re),s=4,color='black')
   axCent.contour(arange(limits[0][0],limits[0][1],pixdx[0]),\
      arange(limits[1][0],limits[1][1],pixdx[1]),model.model.swapaxes(0,1),6,colors=colors[0])
   
   yf = 0.5
   axCent.set_yticks(arange(limits[1][0], limits[1][1], 0.5))
   axCent.set_yticklabels(arange(limits[1][0], limits[1][1], 0.5))
   if drop == 'b':
      axCent.set_xlabel(r'GALFIT MAG in $i_{775}$')
      axCent.set_ylabel(r'GALFIT $\log_{10}(R_e)$ in $i_{775}$ [pixel]')
   elif drop == 'v': 
      axCent.set_xlabel(r'GALFIT MAG in $z_{850}$')
      axCent.set_ylabel(r'GALFIT $\log_{10}(R_e)$ in $z_{850}$ [pixel]')
   #plt.suptitle(r'$\vec{\mathbf{P}}=[%.2f, %.2f, %.2f$",$ %.2f, %.2f]$' % (par[0],par[1],
   #   10.**par[2]*0.03,par[3],par[4]),size=20)
   #axCent.text(0.1,0.1,"par = [%.2f,%.2f,%.2f,%.2f,%.2f]"%tuple(par),transform=axCent.transAxes,
   #   color='black')
   # plot a straight line through the observed points
   ydata = findline(mag, log10(re))
   xr = arange(limits[0][0], limits[0][1], 0.02)
   axCent.plot(xr, ydata[0] * xr + ydata[1], ':', c='black')
   # plot the straight line corresponding to the power-law relation with logR0 and beta
   m0 = -21.0 + mc(zmean)
   b = par[2] + 0.4 * par[4] * m0
   axCent.plot(xr, -0.4 * par[4] * xr + b, '--', lw=2.5, c=colors[3])
   
   axCent.text(0.1,0.9,dataset,transform=axCent.transAxes,color='black',size=14)

   # Show LF
   axLF = divider.append_axes("top", size=1.2, pad=0.0, sharex=axCent)
   n,bins = histogram(mag, arange(limits[0,0], limits[0,1]+lfbw, lfbw))
   nerr = [sqrt(n), sqrt(n)]
   for i in range(len(n)):
      if n[i]==1: nerr[0][i] = 1.-1.e-3
   axLF.errorbar(bins[:-1]+lfbw/2., n, yerr=nerr, fmt='.', ms=14.,\
      mfc='black', ls='None', mec='black', ecolor='black', capsize=6)
   LF = model.model.sum(axis=1) # LF here contains the volume already
   LFtot = sum(LF) * pixdx[0] / lfbw
   normfactor = npts / LFtot  # normalize the LF to predict the total number of points
   LF = LF * normfactor
   LF0 = model0.model.sum(axis=1)
   LF0 = LF0 * normfactor
   LF_ds = model_ds.model.sum(axis=1)
   LF_ds = LF_ds * normfactor
   axLF.semilogy(arange(limits[0,0],limits[0,1],pixdx[0]), LF, color = colors[2],\
      nonposy='mask', label='GALFIT TF')
   axLF.semilogy(arange(limits[0,0],limits[0,1],pixdx[0]), LF_ds, color = colors[1],\
      ls=':',lw=2,nonposy='mask',label='w/ dropout sel. kernel')
   axLF.semilogy(arange(limits[0,0],limits[0,1],pixdx[0]), LF0, color=colors[0],\
      ls='--',lw=2,nonposy='mask',label='Schechter')
   axLF.set_yticks([1.e-2,1.,1.e2])   
   axLF.set_ylim(1.e-2,max(n)*50.)
   xf = 1.0
   axCent.set_xticks(arange(limits[0][0], limits[0][1]+1., 1.))
   axCent.set_xticklabels(arange(limits[0][0], limits[0][1]+1., 1.))
   axCent.set_xlim(limits[0][0], limits[0][1])
   #axLF.legend(loc='upper left', prop=fp)
   

   # Show size distribution
   axSD = divider.append_axes("right", size="35%", pad=0.0, sharey=axCent)
   n, bins = histogram(log10(re), arange(limits[1,0], limits[1,1]+sdbw, sdbw))
   nerr = [sqrt(n), sqrt(n)]
   for i in range(len(n)):
      if n[i]==1: nerr[0][i] = 1.-1.e-3
   axSD.errorbar(n, bins[:-1]+sdbw/2., xerr=nerr, fmt='.', ms=14,\
      mfc='black', ls='None', mec='black', ecolor='black', capsize=6)
   SD = model.model.sum(axis=0)
   SDtot = sum(SD) * pixdx[1] / sdbw
   normfactor = npts / SDtot
   SD = SD * normfactor
   sizer = arange(limits[1,0], limits[1,1], pixdx[1])
   axSD.semilogx(SD, sizer, color=colors[2], label='GALFIT TF')
   SD0 = model0.model.sum(axis=0)
   SD0 = SD0 * normfactor
   SD_ds = model_ds.model.sum(axis=0)
   SD_ds = SD_ds * normfactor
   axSD.semilogx(SD_ds, sizer, color=colors[1], ls=':', lw=2, label='w/ dropout sel. kernel')
   axSD.semilogx(SD0, sizer, color=colors[0], ls='--', lw=2, label='lognormal')
   axSD.set_xticks([1.,10.,1.e2])
   axSD.set_xlim(1.e-2,max(SD)*5)
   #axSD.legend(loc='lower right', prop=fp)
   #axCent.set_ylim(limits[1][0], limits[1][1])
   axCent.set_ylim(-0.6, 1.8)
   plt.draw()
   fig1.show()

   for tl in axLF.get_xticklabels():
      tl.set_visible(False)
   for tl in axSD.get_yticklabels():
      tl.set_visible(False)
   return model,axCent,axSD,axLF,fig1


def plot_bdropsfit(par, colors=['0.0','0.2','0.4','0.6'], orientation='vertical'):
   if orientation == 'vertical':
      fig = plt.figure(figsize=(7,12))
      ax1 = fig.add_subplot(211)
   elif orientation == 'horizontal':
   	fig = plt.figure(figsize=(10,6))
   	ax1 = fig.add_subplot(121)
   p1 = show_model(par, 'b', 'goods', newfig=False, axCent=ax1, fig1=fig,
      colors=colors)
   if orientation == 'vertical':
      ax2 = fig.add_subplot(212)
   elif orientation == 'horizontal':
      ax2 = fig.add_subplot(122)
   p2 = show_model(par, 'b', 'udf', newfig=False, axCent=ax2, fig1=fig,
      colors=colors)
   return 0


def plot_vdropsfit(par, colors=['0.0','0.2','0.4','0.6']):
   if orientation == 'vertical':
      fig = plt.figure(figsize=(7,12))
      ax1 = fig.add_subplot(211)
   elif orientation == 'horizontal':
   	fig = plt.figure(figsize=(10,6))
   	ax1 = fig.add_subplot(121)
   p1 = show_model(par, 'v', 'goods', newfig=False, axCent=ax1, fig1=fig,
      colors=colors)
   if orientation == 'vertical':
      ax2 = fig.add_subplot(212)
   elif orientation == 'horizontal':
      ax2 = fig.add_subplot(122)
   p2 = show_model(par, 'v', 'udf', newfig=False, axCent=ax2, fig1=fig,
      colors=colors)
   return 0


def plot_sizedist(parb, parv):
   mod_lograrr = arange(bl.limits1[1][0], bl.limits1[1][1], 0.02)
   magb1, reb1, critb1 = fl.cleandata('bdrops_gf_v2.cat', drop='b', limits=bl.limits1, zlo=3.0)
   magb2, reb2, critb2 = fl.cleandata('bdrops_udf_gf_v2.cat', chisqnulim=5.0, magautolim=28.5,\
      limits=bl.limits2, drop='b', zlo=3.0)
   kgridb1 = mlutil.readkgrid('kernel_I.p')
   kgridb2 = mlutil.readkgrid('kernel_I_udf.p')
   mci = bl.mconvert('M1500_to_i.txt')
   reb = concatenate((reb1, reb2))
   modelb1 = bl.bivariate_lf(parb, bl.limits1, bl.pixdx, 'b', 'goods', kgrid=kgridb1, 
      zdgridfile='zdgrid_bdrops.p', mcfile='M1500_to_i.txt', meankcorr=mci(4.0))
   modelb2 = bl.bivariate_lf(parb, bl.limits2, bl.pixdx, 'b', 'udf', kgrid=kgridb2,
      zdgridfile='zdgrid_bdrops_udf.p', mcfile='M1500_to_i.txt', meankcorr=mci(4.0))
   sizedist_b = (modelb1.model.sum(axis=0) + modelb2.model.sum(axis=0))/ 10.**mod_lograrr
   magv1, rev1, critv1 = fl.cleandata('vdrops_gf_v2.cat', chisqnulim=0.5, drop='v',
      limits=bl.limits1, zlo=4.0)
   magv2, rev2, critv2 = fl.cleandata('vdrops_udf_gf_v2.cat', chisqnulim=5.0, magautolim=28.5,\
      limits=bl.limits2, drop='v', zlo=4.0)
   rev = concatenate((rev1, rev2))
   kgridv1 = mlutil.readkgrid('kernel_Z.p')
   kgridv2 = mlutil.readkgrid('kernel_Z_udf.p')
   mcz = bl.mconvert('M1500_to_z.txt')
   modelv1 = bl.bivariate_lf(parv, bl.limits1, bl.pixdx, 'v', 'goods', kgrid=kgridv1,
      zdgridfile='zdgrid_vdrops.p', mcfile='M1500_to_z.txt', meankcorr=mcz(5.0))
   modelv2 = bl.bivariate_lf(parv, bl.limits2, bl.pixdx, 'v', 'udf', kgrid=kgridv2,
      zdgridfile='zdgrid_vdrops_udf.p', mcfile='M1500_to_z.txt', meankcorr=mcz(5.0))
   sizedist_v = (modelv1.model.sum(axis=0) + modelb2.model.sum(axis=0))/ 10.**mod_lograrr
   # fit the uncorrected size distribution with lognormal function
   xout_b = fln.fit_lognormal(drop='b'); print xout_b
   xout_v = fln.fit_lognormal(drop='v'); print xout_v
   
   # plot
   fig = plt.figure(figsize=(8,10))
   ax1 = fig.add_subplot(211)
   rarr = arange(0.001, 41., 1.)
   pl_rarr = arange(0.001, 41., 0.01)
   
   h1 = ax1.hist(reb, rarr, color='gray', ec='none')
   fb = fln.lognormal(xout_b[0], xout_b[1], pl_rarr) 
   ax1.plot(pl_rarr, fb * max(h1[0]) / max(fb), color='black', label=r'$\sigma=%.2f$' %
      xout_b[1])
   ax1.plot(10.**mod_lograrr, sizedist_b * max(h1[0]) / max(sizedist_b), 
      color='red', label=r'$\sigma=%.2f$; corrected'%parb[3]) 

   ax2 = fig.add_subplot(212)
   h2 = ax2.hist(rev, rarr, color='gray', ec='none')
   fv = fln.lognormal(xout_v[0], xout_v[1], pl_rarr)
   ax2.plot(pl_rarr, fv * max(h2[0]) / max(fv), color='black', label=r'$\sigma=%.2f$' %
      xout_v[1])
   ax2.plot(10.**mod_lograrr, sizedist_v * max(h2[0]) / max(sizedist_v), 
      color='red', label=r'$\sigma=%.2f$; corrected'%parv[3])

   ax1.set_xlim(0, 40); ax2.set_xlim(0, 40)
   ax1.set_xlabel('Re [0.03" / pixel]')
   ax2.set_xlabel('Re [0.03" / pixel]')
   ax1.legend(loc=1)
   ax2.legend(loc=1)   
   ax1.set_title('B-dropouts (z~4)')
   ax2.set_title('V-dropouts (z~5)')

   return fig



