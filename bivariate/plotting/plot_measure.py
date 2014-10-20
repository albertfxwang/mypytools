#!/usr/bin/env python


from numpy import *
from pygoods import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import fit_bdrops as fb
#import fit_vdrops as fv
import bivRL as bl
import cosmoclass
import colorsys
from my_mplfonts import Helvetica


cc = cosmoclass.cosmoclass(H0=70.2,omega_m=0.27,omega_l=0.73)

def fprop(fsize):
	return mpl.font_manager.FontProperties(size=fsize)

class plot_measure_bdrops(fb.FitBdropsRL):
   def __init__(self, parfile):
      fb.FitBdropsRL.__init__(self, parfile, for_fit=False)
      self.c1 = Ftable(self.LBGCAT1)
      self.c2 = Ftable(self.LBGCAT2)
      self.z_mean = 4.0
      self.pixscale = 0.03
      self.crit1 = self.goodscrit
      self.crit2 = self.udfcrit

   def change_xtick_labelsize(self, ax, newsize=13, font=Helvetica):
      ax.set_xticklabels(ax.get_xticks(), ax.get_xticklabels(), 
                         font_properties=font(newsize))

   def change_ytick_labelsize(self, ax, newsize=13, font=Helvetica):
      ax.set_yticklabels(ax.get_yticks(), ax.get_yticklabels(),
                         font_properties=font(newsize))

   def plot_mag_Re(self, colorg='blue', coloru='red', absmag=False, kpc=False, 
                   outcat=None, ax=None):
      # plot the GALFIT measurement results
      mag1, re1 = self.datasets['goods']
      mag2, re2 = self.datasets['udf']
      if absmag:

         mag1 = mag1 - self.mc['goods'](self.z_mean)
         mag2 = mag2 - self.mc['udf'](self.z_mean)
      if kpc:
         re1 = re1 * cc.adsize(self.z_mean, self.pixscale, unit='kpc')
         re2 = re2 * cc.adsize(self.z_mean, self.pixscale, unit='kpc')
      if ax == None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
      ax.semilogy(mag1, re1, '.', color=colorg, label='GOODS')
      ax.semilogy(mag2, re2, 'x', color=coloru, label='HUDF')
      #yt = array([0.1,1.0,10.,100.])
      yt = array([-1.,0.,1.,2.])
      #if kpc==False:
      #	plt.yticks(yt, 0.03*yt)
      ax.legend(loc=1, numpoints=1)
      # plot the median error bars
      reerr1 = self.c1.reout_err[self.crit1] 
      reerr2 = cf2.reout_err[self.crit2]
      reerr1_plot = median(reerr1) 
      reerr2_plot = median(reerr2)
      magerr1 = self.c1.magout_err[self.crit1]
      magerr2 = self.c2.magout_err[self.crit2]
      magerr1_plot = median(magerr1) 
      magerr2_plot = median(magerr2)
      print "reerr1_plot, reerr2_plot = ", reerr1_plot, reerr2_plot
      print "magerr1_plot, magerr2_plot = ", magerr1_plot, magerr2_plot
      emag1 = 23.0; emag2 = 23.3
      ere1 = 0.5; ere2 = 0.5
      ax.set_yscale('log')
   
   def plot_nhist(self, nhi=10., dn=0.2, fraction=True, color='blue',
                  ax=None, label='', labelsize=18):
      """
      Plot GOODS + UDF together.
      """
      #mpl.rcParams['xtick.labelsize'] = 18
      #mpl.rcParams['xtick.labelsize'] = 18
      #mpl.rcParams['axes.labelsize'] = 20
      #mpl.rcParams['legend.fontsize'] = 16
      if ax == None:
         fig = plt.figure(figsize=(8,5))
         ax = fig.add_subplot(111)
      #c1b = sextractor(cat1b); c2b = sextractor(cat2b)
      #c1v = sextractor(cat1v); c2v = sextractor(cat2v)
      n_goods = self.c1.nout[self.crit1==True]
      n_udf = self.c2.nout[self.crit2==True]
      n_array = concatenate([n_goods, n_udf])
      if fraction == True:
         weights = ones(len(n_array))*(1./float(len(n_array)))
      x = ax.hist(n_array, arange(0., nhi, dn), histtype='step', 
                   edgecolor=color, lw=2.0, weights=weights, 
                   label=label)
      ax.set_xlabel('GALFIT Sersic index $n$', 
                    font_properties=Helvetica(labelsize))
      ax.set_xticks(arange(0.,nhi))
      if fraction:
         ax.set_ylabel('Fraction', font_properties=Helvetica(labelsize))
      else:
         ax.set_ylabel('Number', font_properties=Helvetica(labelsize))
      yticks = ax.get_yticks()
      ax.set_yticks(arange(0.,yticks[-1],0.04))
      
   def plot_ellipticity(self, d_ar=0.05, color='blue', ax=None, label=''):
      """
      Plot the ellipticity distribution.
      """
      axratio_array1 = self.c1.arout[self.crit1==True]
      axratio_array2 = self.c2.arout[self.crit2==True]
      axratio_array = concatenate([axratio_array1, axratio_array2])
      # define ellipticity array
      ell_array = 1. - axratio_array
      if ax == None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
      x = ax.hist(ell_array, bins=arange(0.,1.,d_ar), histtype='step', lw=2.0,
                  edgecolor=color, label=label)
      ax.set_xlabel('Ellipticity (1 - b/a)', font_properties=Helvetica(18))
      ax.set_ylabel('Number', font_properties=Helvetica(18))

   def plot_axis_ratio(self, d_ar=0.05, color='blue', ax=None, label=''):
      """
      Plot the axis ratio distribution.
      """
      axratio_array1 = self.c1.arout[self.crit1==True]
      axratio_array2 = self.c2.arout[self.crit2==True]
      axratio_array = concatenate([axratio_array1, axratio_array2])
      if ax == None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
      x = ax.hist(axratio_array, bins=arange(0.,1.,d_ar), histtype='step', 
                  lw=2.0, edgecolor=color, label=label)
      ax.set_xlabel('Axis Ratio (b/a)', font_properties=Helvetica(18))
      ax.set_ylabel('Number', font_properties=Helvetica(18))

class plot_measure_vdrops(plot_measure_bdrops):
   def __init__(self, parfile):
      plot_measure_bdrops.__init__(self,parfile)
      self.z_mean = 5.0


def plot_bdropouts_n(cat1, cat2, nlow=0., nhigh=8., s=0.5, v=0.5, alpha=0.8, deghi=240.):
	# plot the GALFIT measurements of B-dropouts, only to be color coded by their Sersic index
	# calculate the RGB color of each point from HSV system, where H is the hue calculated
	# by the Sersic index n. The range of hue is from 0 degree (red) to deghi degree,
	# with n=0 corresponding to purple and n>=8 corresponding to red
	mag1, re1, crit1 = fl.cleandata(cat1, reerr_ratiolim=0.6, chisqnulim=0.4,
		zlo=-1.0, drop='b', limits=bl.limits1)
	mag2, re2, crit2 = fl.cleandata(cat2, reerr_ratiolim=0.6, chisqnulim=5.0,
		zlo=-1.0, drop='b', limits=bl.limits2)
	c1 = sextractor(cat1); c2 = sextractor(cat2)
	n1 = c1.nout[crit1]
	n2 = c2.nout[crit2]
	n1_cap = minimum(n1, nhigh); n2_cap = minimum(n2, nhigh)  # enforce maximum value of n
	print len(n1_cap), len(n2_cap)
	delta_n = nhigh - nlow; print delta_n
	hue1 = deghi - (n1_cap / delta_n * deghi)  # convert n=0 to purple (315 degree)
	hue2 = deghi - (n2_cap / delta_n * deghi)
	rgbcolor1 = map(colorsys.hsv_to_rgb, hue1, array([s]).repeat(len(mag1)), array([v]).repeat(len(mag1)))
	rgbcolor2 = map(colorsys.hsv_to_rgb, hue2, array([s]).repeat(len(mag2)), array([v]).repeat(len(mag2)))
	fig = plt.figure()
	plt.scatter(mag1, log10(re1), s=8**2, marker='o', c=rgbcolor1, edgecolor='none', alpha=alpha)
	plt.scatter(mag2, log10(re2), s=8**2, marker=7, c=rgbcolor2, edgecolor='none', alpha=alpha)
	plt.xlabel('GALFIT F775W magnitude')
	plt.ylabel('GALFIT log10(Re) [pixels]')
	plt.title('B-dropouts')
	

def plot_bdropouts_color(cat1, cat2, nlow=0., nhigh=8., s=0.5, v=0.5, alpha=0.8, deghi=240.):
	# plot the GALFIT measurements of B-dropouts, only to be color coded by their Sersic index
	# calculate the RGB color of each point from HSV system, where H is the hue calculated
	# by the Sersic index n. The range of hue is from 0 degree (red) to deghi degree,
	# with n=0 corresponding to purple and n>=8 corresponding to red
	mag1, re1, crit1 = fl.cleandata(cat1, reerr_ratiolim=0.6, chisqnulim=0.4,
		zlo=-1.0, drop='b', limits=bl.limits1)
	mag2, re2, crit2 = fl.cleandata(cat2, reerr_ratiolim=0.6, chisqnulim=5.0,
		zlo=-1.0, drop='b', limits=bl.limits2)
	c1 = sextractor(cat1); c2 = sextractor(cat2)
	n1 = c1.nout[crit1]
	n2 = c2.nout[crit2]
	n1_cap = minimum(n1, nhigh); n2_cap = minimum(n2, nhigh)  # enforce maximum value of n
	print len(n1_cap), len(n2_cap)
	delta_n = nhigh - nlow; print delta_n
	hue1 = deghi - (n1_cap / delta_n * deghi)  # convert n=0 to purple (315 degree)
	hue2 = deghi - (n2_cap / delta_n * deghi)
	rgbcolor1 = map(colorsys.hsv_to_rgb, hue1, array([s]).repeat(len(mag1)), array([v]).repeat(len(mag1)))
	rgbcolor2 = map(colorsys.hsv_to_rgb, hue2, array([s]).repeat(len(mag2)), array([v]).repeat(len(mag2)))
	fig = plt.figure()
	plt.scatter(mag1, log10(re1), s=8**2, marker='o', c=rgbcolor1, edgecolor='none', alpha=alpha)
	plt.scatter(mag2, log10(re2), s=8**2, marker=7, c=rgbcolor2, edgecolor='none', alpha=alpha)
	plt.xlabel('GALFIT F775W magnitude')
	plt.ylabel('GALFIT log10(Re) [pixels]')
	plt.title('B-dropouts')
	
	
def plot_vdropouts(c1, c2, colorg='blue', coloru='red', absmag=False, kpc=False, outcat=None):
   # plot the GALFIT measurement results
   mag1, re1, crit1 = fl.cleandata(c1, reerr_ratiolim=0.6, chisqnulim=0.5,
      zlo=-1.0, drop='v', limits=bl.limits1)
   mag2, re2, crit2 = fl.cleandata(c2, reerr_ratiolim=0.6, chisqnulim=5.0,
      zlo=-1.0, drop='v', limits=bl.limits2)
   print "sum(crit1), sum(crit2)", sum(crit1), sum(crit2)
   if absmag:
   	mag1 = mag1 - bl.meankcorr_z5_z   # convert to absolute magnitude
   	mag2 = mag2 - bl.meankcorr_z5_z
   if kpc:
   	re1 = re1 * cc.adsize(5.0, 0.03, unit='kpc')
   	re2 = re2 * cc.adsize(5.0, 0.03, unit='kpc')
   cf1 = sextractor(c1); cf2 = sextractor(c2)
   plt.semilogy(mag1, re1, '.', color=colorg, label='GOODS')
   plt.semilogy(mag2, re2, 'x', color=coloru, label='HUDF')
   # plot the median error bars
   reerr1 = cf1.reout_err[crit1]; reerr2 = cf2.reout_err[crit2]
   reerr1_plot = median(reerr1); reerr2_plot = median(reerr2)
   magerr1 = cf1.magout_err[crit1]; magerr2 = cf2.magout_err[crit2]
   magerr1_plot = median(magerr1); magerr2_plot = median(magerr2)
   print "reerr1_plot, reerr2_plot =", reerr1_plot, reerr2_plot
   print "magerr1_plot, magerr2_plot =", magerr1_plot, magerr2_plot
   emag1 = 23.0; emag2 = 23.3
   ere1 = 0.5; ere2 = 0.5
   plt.yscale('log')
   #plt.errorbar([emag1],[ere1],xerr=magerr1_plot,yerr=reerr1_plot,fmt='.',elinewidth=2.0,
   #	ecolor=colorg,color=colorg)
   #plt.errorbar([emag2],[ere2],xerr=magerr2_plot,yerr=reerr2_plot,fmt='x',elinewidth=2.0,
  # 	ecolor=coloru,color=coloru)
   if kpc == False:
		#yt = array([0.1,1.0,10.,100.])
		yt = array([-1.,0.,1.,2.])
		#plt.yticks(yt, 0.03*yt)
   plt.legend(loc='upper right', numpoints=1)
   if outcat:
   	f = open(outcat,'w')
   	f.write('# 1 mag\n')
   	f.write('# 2 Re\n')
   	for i in range(len(mag1)):
   		f.write('%f  %f  ' % (mag1[i], re1[i]))
   		f.write('\n')
   	for i in range(len(mag2)):  # write B-dropouts in HUDF
   		f.write('%f  %f  ' % (mag2[i], re2[i]))
   		f.write('\n')
   	f.close()
   return 0


def plot_all(c11, c12, c21, c22, colorg='blue', coloru='red'):
   # plot both B-dropouts and V-dropouts
   # c11 & c12: B-dropouts catalogs in GOODS & UDF
   # c21 & c22: V-dropouts catalogs in GOODS & UDF
   mpl.rcParams['xtick.labelsize'] = 20
   mpl.rcParams['ytick.labelsize'] = 20
   mpl.rcParams['lines.markersize'] = 10
   mpl.rcParams['axes.labelsize'] = 22
   if type(c11) != type('abc'):
      print "Please provide catalog file names."
      raise
   fig = plt.figure(figsize=(13,6))
   #plt.subplot(1,2,1)
   ax1 = fig.add_axes([0.1,0.1,0.35,0.8])
   plot_bdropouts(c11, c12, colorg=colorg, coloru=coloru)
   plt.title('$B$-dropouts', size=22)
   plt.xlabel('GALFIT i-band mag')
   #plt.ylabel('Re [arcsec]')
   plt.ylabel(r'$R_e$ [pixel]')
   plt.xlim(22.5,28.5)
   plt.ylim(0.1, 100.)
   #plt.subplot(1,2,2)
   ax2 = fig.add_axes([0.6,0.1,0.35,0.8])
   plot_vdropouts(c21, c22, colorg=colorg, coloru=coloru)
   plt.title('$V$-dropouts', size=22)
   plt.xlabel('GALFIT z-band mag')
   #plt.ylabel('Re [arcsec]')
   plt.ylabel(r'$R_e$ [pixel]')
   plt.xlim(22.5,28.5)
   plt.ylim(0.1, 100.)
   return 0


def plot_bdrops_nhist(cat1, cat2, axes, nhi=10., dn=0.2, color='blue'):
   # plot the GALFIT-measured Sersic index distribution 
   mag1, re1, crit1 = fl.cleandata(cat1, reerr_ratiolim=0.6, chisqnulim=0.4, magautolim=26.5,\
      zlo=3.0, drop='b', limits=bl.limits1)
   mag2, re2, crit2 = fl.cleandata(cat2, reerr_ratiolim=0.6, chisqnulim=5.0, magautolim=28.5,\
      zlo=3.0, drop='b', limits=bl.limits2)
   c1 = sextractor(cat1); c2 = sextractor(cat2)
   narr1 = compress(crit1, c1.nout)
   narr2 = compress(crit2, c2.nout)
   narr12 = concatenate((narr1,narr2))
   x = axes.hist(narr12, arange(0., nhi, dn), histtype='step', lw=1.5, color=color)
   return axes, narr12





def plot_bdrops_nsep(cat1, cat2, nsep=2.5, title='B-dropouts'):
   # plot GALFIT-measure mag & Re of B-dropouts separately for n<=2.5 and n>2.5 objects
   fig = plt.figure()
   ax1 = fig.add_subplot(111)
   mag1, re1, crit1 = fl.cleandata(cat1, reerr_ratiolim=0.6, chisqnulim=0.4, magautolim=26.5,\
      zlo=3.0, drop='b', limits=bl.limits1)
   mag2, re2, crit2 = fl.cleandata(cat2, reerr_ratiolim=0.6, chisqnulim=5.0, magautolim=28.5,\
      zlo=3.0, drop='b', limits=bl.limits2)
   c1 = sextractor(cat1); c2 = sextractor(cat2)
   mag1_disk = compress(crit1 & (c1.nout<=2.5), c1.magout)
   mag1_devauc = compress(crit1 & (c1.nout>2.5), c1.magout)
   re1_disk = compress(crit1 & (c1.nout<=2.5), c1.reout)
   re1_devauc = compress(crit1 & (c1.nout>2.5), c1.reout)
   mag2_disk = compress(crit2 & (c2.nout<=2.5), c2.magout)
   mag2_devauc = compress(crit2 & (c2.nout>2.5), c2.magout)
   re2_disk = compress(crit2 & (c2.nout<=2.5), c2.reout)
   re2_devauc = compress(crit2 & (c2.nout>2.5), c2.reout)
   ax1.plot(mag1_disk, log10(re1_disk), marker='s', ls='none', mec='none', mfc='blue', ms=4,
      label=r'GOODS; $n \leq 2.5$')
   ax1.plot(mag2_disk, log10(re2_disk), marker='s', ls='none', mec='none', mfc='red', ms=4,
      label=r'UDF; $n \leq 2.5$')
   ax1.plot(mag1_devauc, log10(re1_devauc), marker='o', ls='none', mec='blue', mfc='none',
      ms=6, mew=1.1, label=r'GOODS; $n > 2.5$')
   ax1.plot(mag2_devauc, log10(re2_devauc), marker='o', ls='none', mec='red', mfc='none',
      ms=6, mew=1.1, label=r'UDF; $n > 2.5$')
   ax1.legend(loc=3, prop=fp10)
   ax1.set_xlabel('i-band GALFIT mag')
   ax1.set_yticks([-1.,0.,1.,2.])
   ax1.set_yticklabels([0.003,0.03,0.3,3.0])
   ax1.set_ylabel('i-band GALFIT Re [arcsec]')
   ax1.set_title(title)
   return ax1


def plot_vdrops_nsep(cat1, cat2, nsep=2.5, title='V-dropouts'):
   # plot GALFIT-measure mag & Re of B-dropouts separately for n<=2.5 and n>2.5 objects
   fig = plt.figure()
   ax1 = fig.add_subplot(111)
   mag1, re1, crit1 = fl.cleandata(cat1, reerr_ratiolim=0.6, chisqnulim=0.5, magautolim=26.5,\
      zlo=4.0, drop='v', limits=bl.limits1)
   mag2, re2, crit2 = fl.cleandata(cat2, reerr_ratiolim=0.6, chisqnulim=5.0, magautolim=28.5,\
      zlo=4.0, drop='v', limits=bl.limits2)
   c1 = sextractor(cat1); c2 = sextractor(cat2)
   mag1_disk = compress(crit1 & (c1.nout<=2.5), c1.magout)
   mag1_devauc = compress(crit1 & (c1.nout>2.5), c1.magout)
   re1_disk = compress(crit1 & (c1.nout<=2.5), c1.reout)
   re1_devauc = compress(crit1 & (c1.nout>2.5), c1.reout)
   mag2_disk = compress(crit2 & (c2.nout<=2.5), c2.magout)
   mag2_devauc = compress(crit2 & (c2.nout>2.5), c2.magout)
   re2_disk = compress(crit2 & (c2.nout<=2.5), c2.reout)
   re2_devauc = compress(crit2 & (c2.nout>2.5), c2.reout)
   ax1.plot(mag1_disk, log10(re1_disk), marker='s', ls='none', mec='none', mfc='blue', ms=4,
      label=r'GOODS; $n \leq 2.5$')
   ax1.plot(mag2_disk, log10(re2_disk), marker='s', ls='none', mec='none', mfc='red', ms=4,
      label=r'UDF; $n \leq 2.5$')
   ax1.plot(mag1_devauc, log10(re1_devauc), marker='o', ls='none', mec='blue', mfc='none',
      ms=6, mew=1.1, label=r'GOODS; $n > 2.5$')
   ax1.plot(mag2_devauc, log10(re2_devauc), marker='o', ls='none', mec='red', mfc='none',
      ms=6, mew=1.1, label=r'UDF; $n > 2.5$')
   ax1.legend(loc=3, prop=fp10)
   ax1.set_xlabel('z-band GALFIT mag')
   ax1.set_yticks([-1.,0.,1.,2.])
   ax1.set_yticklabels([0.003,0.03,0.3,3.0])
   ax1.set_ylabel('z-band GALFIT Re [arcsec]')
   ax1.set_title(title)
   return ax1
