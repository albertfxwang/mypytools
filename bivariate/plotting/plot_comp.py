#!/usr/bin/env python

from numpy import *
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cPickle
from pygoods import *
#import fit_lbg as fl
from my_mplfonts import Helvetica

cmap_hot = mpl.cm.hot
cmap_jet = mpl.cm.jet
cmap_hsv = mpl.cm.hsv

class plot_comp(object):
   def showmap(self, marray, bw=array([0.5,0.2]), 
      dtype='float', vrange=None, xdata=[], ydata=[],
      newfig=True, fig=None, ax=None, cmap=cmap_hot, tickw=array([1.0, 0.4]), 
      cbar_extend='neither', aspect='auto'):
      """
      Show the values of marray in each pixel on a map
      can decide to mark the axes at every nth bin
      Assume the object has the following attributes:
      pixdx
      limits
      """
      # nbins = (self.limits[:,1] - self.limits[:,0]) / bw
      # nbins = nbins.astype('int')
      # if not (shape(marray) == nbins).all():
      #    print "shape of input array is [%d, %d], incompatible with number of bins [%d, %d]" % (
      #       shape(marray)[0], shape(marray)[1], nbins[0], nbins[1])
      #    raise ValueError
      # xbins = arange(limits[0][0], limits[0][1], bw[0])
      # ybins = arange(limits[1][0], limits[1][1], bw[1])
      npix = (self.limits[:,1] - self.limits[:,0]) / self.pixdx
      xticks = arange(0, npix[0], tickw[0]/self.pixdx[0])
      xticks = xticks.astype('int')
      yticks = arange(0, npix[1], tickw[1]/self.pixdx[1])
      yticks = yticks.astype('int')
      xtickvals = self.limits[0][0] + xticks * self.pixdx[0]
      ytickvals = self.limits[1][0] + yticks * self.pixdx[1]
      #print len(xticks), len(xtickvals), len(yticks), len(ytickvals)
      #print yticks, ytickvals
      xticklabels = []; yticklabels = []
      for k in range(len(xticks)):
         xticklabels += ['%.1f' % xtickvals[k]]
      for k in range(len(yticks)):
         yticklabels += ['%.1f' % ytickvals[k]]
                  
      if newfig == True:
         fig = plt.figure()
         ax = fig.add_subplot(111)
      if vrange==None:
         vrange=[min(marray.ravel()), max(marray.ravel())]
      cax = ax.imshow(marray.swapaxes(0,1), origin='lower', vmin=vrange[0],
         vmax=vrange[1], cmap=cmap, aspect=aspect)
      cbar = fig.colorbar(cax, ax=ax, extend=cbar_extend)
      ax.set_xticks(xticks); ax.set_xticklabels(xticklabels)
      ax.set_yticks(yticks); ax.set_yticklabels(yticklabels)
      # plot the points overlaid on top
      if len(xdata):
         xpixdata = (xdata - self.limits[0][0]) / self.pixdx[0]
         ypixdata = (ydata - self.limits[1][0]) / self.pixdx[1]
         ax.plot(xpixdata, ypixdata, 'o', mec='white', mfc='green', ms=3)
         ax.set_xlim(0, npix[0]); ax.set_ylim(0, npix[1])
      self.fig = fig
      self.ax = ax
      self.cbar = cbar

class plot_zdgrid_comp(plot_comp):
   def __init__(self, zdgrid, pixdx=array([0.02, 0.02])):
      for k in zdgrid.__dict__.keys():
         setattr(self, k, getattr(zdgrid, k))
      self.pixdx = pixdx
      self.limits = array([[self.M0, self.M1], [self.logR0, self.logR1]])
      nx = (self.M1 - self.M0) / self.pixdx[0]
      self.nx = int(round(nx))
      ny = (self.logR1 - self.logR0) / self.pixdx[1]
      self.ny = int(round(ny))

   def zdgrid_detcomp(self, show=True, **kwargs):
      
      compgrid = zeros((self.nx, self.ny), 'float')
   
      for zd in self.zdist:
         kern = self.zdist[zd]
         ix0 = (kern.M0 - self.M0) / self.pixdx[0]
         ix0 = int(round(ix0))
         ix1 = (kern.M1 - self.M0) / self.pixdx[0]
         ix1 = int(round(ix1))
         iy0 = (kern.logR0 - self.logR0) / self.pixdx[1]
         iy0 = int(round(iy0))
         iy1 = (kern.logR1 - self.logR0) / self.pixdx[1]
         iy1 = int(round(iy1))
         detcomp = float(sum(kern.ndetect)) / float(sum(kern.ninput))
         compgrid[ix0:ix1,iy0:iy1] = detcomp
      # self.detcompgrid = compgrid
      if show:
         self.showmap(compgrid, **kwargs)
         self.ax.set_title("Detection completeness for \n%s" % os.path.split(self.filename)[-1], 
                           size=18)

   def zdgrid_ninput(self, show=True, **kwargs):
      # returns the ninput map
      numgrid = zeros((self.nx, self.ny), 'float')
   
      for zd in self.zdist:
         kern = self.zdist[zd]
         ix0 = (kern.M0 - self.M0) / self.pixdx[0]
         ix0 = int(round(ix0))
         ix1 = (kern.M1 - self.M0) / self.pixdx[0]
         ix1 = int(round(ix1))
         iy0 = (kern.logR0 - self.logR0) / self.pixdx[1]
         iy0 = int(round(iy0))
         iy1 = (kern.logR1 - self.logR0) / self.pixdx[1]
         iy1 = int(round(iy1))
         numgrid[ix0:ix1,iy0:iy1] = sum(kern.ninput)
      if show:
         self.showmap(numgrid, **kwargs)
         self.ax.set_title("Number of Input Galaxies for \n%s" % os.path.split(self.filename)[-1], 
                           size=18)


   def zdgrid_selcomp(self, z0, z1, show=True, **kwargs):
      compgrid = zeros((self.nx, self.ny), 'float')
   
      for zd in self.zdist:
         kern = self.zdist[zd]
         ix0 = (kern.M0 - self.M0) / self.pixdx[0]
         ix0 = int(round(ix0))
         ix1 = (kern.M1 - self.M0) / self.pixdx[0]
         ix1 = int(round(ix1))
         iy0 = (kern.logR0 - self.logR0) / self.pixdx[1]
         iy0 = int(round(iy0))
         iy1 = (kern.logR1 - self.logR0) / self.pixdx[1]
         iy1 = int(round(iy1))
         iz0 = searchsorted(self.zarrold, z0)
         iz1 = searchsorted(self.zarrold, z1)
         # print iz0, iz1
         if sum(kern.ndetect)>0:
            print zd, sum(kern.nselect), sum(kern.ninput)
            detcomp = float(sum(kern.nselect[iz0:iz1])) / float(maximum(1.0, sum(kern.ninput[iz0:iz1])))
            if sum(kern.ninput[iz0:iz1]) <= 4:
               detcomp = 0.
            compgrid[ix0:ix1,iy0:iy1] = detcomp
      if show:
         self.showmap(compgrid, **kwargs)
         self.ax.set_title("Selection Completeness for \n%s" % os.path.split(self.filename)[-1], 
                           size=18)


def sexcompmap(c, maglim, logrelim, magbw, logrebw):
	# the SExtractor detection completeness in zmag-logRe bins
	# calculated from the dropout-simulations catalog (NOT the GALFIT simulations)
	amag = arange(maglim[0], maglim[1]+magbw, magbw)
	alogre = arange(logrelim[0], logrelim[1]+logrebw, logrebw)
   # use histogram2d
	mag_input = c.z_mag_in
	logre_input = log10(c.re_in)
	mag_detect = compress(c.detect, mag_input)
	logre_detect = compress(c.detect, logre_input)
	N_input = histogram2d(mag_input, logre_input, bins=[amag,alogre])[0]
	N_detect = histogram2d(mag_detect, logre_detect, bins=[amag,alogre])[0]
	compgrid = where(N_input>0, N_detect.astype('float')/N_input.astype('float'), 0.)
	#compgrid = zeros((len(amag), len(alogre)), 'float')
	#for i in range(len(amag)):
	#   for j in range(len(alogre)):
	#      crit = (c.magin >= amag[i]) & (c.magin < (amag[i] + magbw))
	#      crit = crit & (log10(c.rein) >= alogre[j]) & (log10(c.rein) < (alogre[j] + logrebw))
	#      adetect = compress(crit, c.recovered_se)
	#      n_input = sum(crit)
	#      n_detect = sum(adetect)
	#      if n_input > 0:
	#         compgrid[i, j] = float(n_detect) / float(n_input)
	return compgrid


def ninputmap(c, maglim, logrelim, magbw, logrebw):
   amag = arange(maglim[0], maglim[1], magbw)
   alogre = arange(logrelim[0], logrelim[1], logrebw)
   numgrid = zeros((len(amag), len(alogre)), 'float')
   for i in range(len(amag)):
      for j in range(len(alogre)):
         crit = (c.magin >= amag[i]) & (c.magin < (amag[i] + magbw))
         crit = crit & (log10(c.rein) >= alogre[j]) & (log10(c.rein) < (alogre[j] + logrebw))
         n_input = sum(crit)
         numgrid[i, j] = n_input
   return numgrid


def gfcompmap(c, maglim, logrelim, magbw, logrebw):
   amag = arange(maglim[0], maglim[1], magbw)
   alogre = arange(logrelim[0], logrelim[1], logrebw)
   compgrid = zeros((len(amag), len(alogre)), 'float')
   for i in range(len(amag)):
      for j in range(len(alogre)):
         crit = (c.magin >= amag[i]) & (c.magin < (amag[i] + magbw))
         crit = crit & (log10(c.rein) >= alogre[j]) & (log10(c.rein) < (alogre[j] + logrebw))
         adetect = compress(crit, c.recovered_galfit)
         n_input = sum(crit)
         n_detect = sum(adetect)
         if n_input > 0:
            compgrid[i, j] = float(n_detect) / float(n_input)
   return compgrid


def gfcompmap2(kgrid):
   # plot the completeness map of the GALFIT transfer function kernel grid
   npix = (kgrid.limits[:,1] - kgrid.limits[:,0]) / kgrid.pixdx 
   npix = npix.astype('int')
   compgrid = zeros(npix)
   for k in kgrid.kernels:
      kern = kgrid.kernels[k]
      i0 = (kern.x0_bin - kgrid.limits[:,0]) / kgrid.pixdx
      i0 = around(i0).astype('int')
      i1 = (kern.x1_bin - kgrid.limits[:,0]) / kgrid.pixdx
      i1 = around(i1).astype('int')
      compgrid[i0[0]:i1[0], i0[1]:i1[1]] = kern.completeness
   compgrid = maximum(compgrid, 0.)

   return compgrid


def gfnummap(kgrid):
   # returns the input numbers of each bin
   npix = (kgrid.limits[:,1] - kgrid.limits[:,0]) / kgrid.pixdx 
   npix = npix.astype('int')
   numgrid = zeros(npix)
   for k in kgrid.kernels:
      kern = kgrid.kernels[k]
      i0 = (kern.x0_bin - kgrid.limits[:,0]) / kgrid.pixdx
      i0 = around(i0).astype('int')
      i1 = (kern.x1_bin - kgrid.limits[:,0]) / kgrid.pixdx
      i1 = around(i1).astype('int')
      numgrid[i0[0]:i1[0], i0[1]:i1[1]] = kern.ninput_bin

   return numgrid


def dkgrid_comp(dkgrid):
   # returns the completeness map
   nx = (dkgrid.M1 - dkgrid.M0) / 0.02
   nx = int(round(nx))
   ny = (dkgrid.logR1 - dkgrid.logR0) / 0.02
   ny = int(round(ny))
   compgrid = zeros((nx, ny), 'float')
   
   for k in dkgrid.dkcorr:
      kern = dkgrid.dkcorr[k]
      ix0 = (kern.M0 - dkgrid.M0) / 0.02
      ix0 = int(round(ix0))
      ix1 = (kern.M1 - dkgrid.M0) / 0.02
      ix1 = int(round(ix1))
      iy0 = (kern.logR0 - dkgrid.logR0) / 0.02
      iy0 = int(round(iy0))
      iy1 = (kern.logR1 - dkgrid.logR0) / 0.02
      iy1 = int(round(iy1))
      compgrid[ix0:ix1,iy0:iy1] = kern.completeness

   return compgrid
      

def dkgrid_ntarget(dkgrid):
   # returns the ntarget map
   nx = (dkgrid.M1 - dkgrid.M0) / 0.02
   nx = int(round(nx))
   ny = (dkgrid.logR1 - dkgrid.logR0) / 0.02
   ny = int(round(ny))
   numgrid = zeros((nx, ny), 'float')

   for k in dkgrid.dkcorr:
      kern = dkgrid.dkcorr[k]
      ix0 = (kern.M0 - dkgrid.M0) / 0.02
      ix0 = int(round(ix0))
      ix1 = (kern.M1 - dkgrid.M1) / 0.02
      ix1 = int(round(ix1))
      iy0 = (kern.logR0 - dkgrid.logR0) / 0.02
      iy0 = int(round(iy0))
      iy1 = (kern.logR1 - dkgrid.logR1) / 0.02
      iy1 = int(round(iy1))
      numgrid[ix0:ix1,iy0:iy1] = kern.n_target

   return numgrid



def read_dkgrid(dkgridfile):
   f = open(dkgridfile)
   dkgrid = cPickle.load(f)
   f.close()
   return dkgrid

def plot_galfit_qf_comp_magbin(c_sim, c_drop, limits=[22.0, 26.5], chisqnulim=0.4):
   # plot the GALFIT completeness for simulations and real dropout sample
   # binned in GALFIT-measured magnitudes
   magbins = arange(limits[0], limits[1]+0.5, 0.5)
   qfcrit_sim = (c_sim.magout>0.)&((c_sim.reouterr/c_sim.reout)<=0.6) & (c_sim.chisqnu<=chisqnulim)
   qfcrit_drop = (c_drop.magout>0.)&((c_drop.reout_err/c_drop.reout)<=0.6) & (c_drop.chisqnu<=chisqnulim)
   magout_qf_sim = compress(qfcrit_sim, c_sim.magout)
   magout_qf_drop = compress(qfcrit_drop, c_drop.magout)
   n_sim = histogram(compress(c_sim.magout>0.,c_sim.magout), magbins)[0]
   n_qf_sim = histogram(compress(qfcrit_sim,c_sim.magout), magbins)[0]
   n_drop = histogram(compress(c_drop.magout>0.,c_drop.magout), magbins)[0]
   n_qf_drop = histogram(compress(qfcrit_drop, c_drop.magout), magbins)[0]
   print n_qf_sim
   n_qf_sim = maximum(n_qf_sim, 1)
   n_qf_drop = maximum(n_qf_drop, 1)   # enforce at least 1 in the denominator
   frac_qf_sim = array(n_qf_sim).astype('float') / array(n_sim).astype('float')
   frac_qf_drop = array(n_qf_drop).astype('float') / array(n_drop).astype('float')
   fig = plt.figure()
   magarr = arange(limits[0], limits[1], 0.5)
   plt.plot(magarr, frac_qf_sim, ls='steps--', label='simulation')
   plt.plot(magarr, frac_qf_drop, ls='steps', label='dropout')
   return frac_qf_sim, frac_qf_drop

	
def show_galfit_qf_comp3(kgrid1, kgrid2, chisqnulim=[0.4,5.0], drop='b', 
	colors=['0.5','0.0'], fields=['GOODS','HUDF']):
	limits1 = around(kgrid1.limits, 1)
	limits2 = around(kgrid2.limits, 1)
	print kgrid1.dxgrid, kgrid2.dxgrid
	nbins = []
	nbins += [around((limits1[:,1]-limits1[:,0])/kgrid1.dxgrid).astype('int')]
	nbins += [around((limits2[:,1]-limits2[:,0])/kgrid2.dxgrid).astype('int')]
	if drop == 'b':
		c1 = Ftable('bdrops_gf_v3.fits')
		c2 = Ftable('bdrops_udf_gf_v3.fits')
	elif drop == 'v':
		c1 = Ftable('vdrops_gf_v2.fits')
		c2 = Ftable('vdrops_udf_gf_v2.fits')
	kgarr = [kgrid1, kgrid2]
	carr = [c1, c2]
	if drop == 'b':
		chisqnulim = [0.4, 5.0]
	elif drop == 'v':
		chisqnulim = [0.5, 5.0]
	plt.rc('lines', linewidth=2.0)
	plt.rc('axes', labelsize=18, titlesize=22)
	plt.rc('xtick', labelsize=13)
	plt.rc('ytick', labelsize=13)
	fp = mpl.font_manager.FontProperties(size=15)
	n_gf = []
	n_gf_qf_sim = []
	n_gf_qf_real = []
	plt.figure(figsize=(12,10))
	for i in range(2):
		# i is iterating over (GOODS, HUDF) fields
		kg = kgarr[i]
		c = carr[i]
		qfcomp = zeros(nbins[i])
		print "shape(qfcomp)", shape(qfcomp)
		for k in kg.kernels:
			#print k
			qfcomp[k[0],k[1]] = maximum(kg.kernels[k].completeness, 0.)
		mbins = arange(kg.limits[0,0], kg.limits[0,1]+kg.dxgrid[0], kg.dxgrid[0])
		rbins = arange(kg.limits[1,0], kg.limits[1,1]+kg.dxgrid[1], kg.dxgrid[1])
		dm = kg.dxgrid[0] / 2.
		dlogr = kg.dxgrid[1] / 2.
		#nbins = around((kg.limits[:,1]-kg.limits[:,0])/kg.dxgrid).astype('int')
		ndrops_gf = np.histogram2d(c.magout[c.reout>0.], log10(c.reout[c.reout>0.]),
			bins=[mbins, rbins])  # number of dropouts in this GALFIT (mag, logRe) bin
		ndrops_gf_qf_sim = ndrops_gf[0] * qfcomp
		#print "shape(ndrops_gf[0])", shape(ndrops_gf[0])
		#print "shape(qfcomp)", shape(qfcomp)
		#print "shape(ndrops_gf_qf_sim", shape(ndrops_gf_qf_sim)
		# number of dropouts expected from simulation that will satisfy our GALFIT quality
		# criteria in this GALFIT (mag, logRe) bin
		qfcrit = (c.reout>0.) & (c.reout_err/c.reout <= 0.6) & (c.chisqnu <= chisqnulim[i])
		ndrops_gf_qf_real = np.histogram2d(c.magout[qfcrit], log10(c.reout[qfcrit]),
			bins=[mbins,rbins])
		# number of real dropouts that satisfy our GALFIT quality criteria in this GALFIT
		# (mag, logRe) bin
		# plot over magnitude bins
		# subplot layout:
		# mag_bin_GOODS   mag_bin_HUDF
		# Re_bin_GOODS    Re_bin_HUDF
		plt.subplot(2, 2, i+1)  # mag bins in GOODS/HUDF
		ndrops_gf_sim_mbins = ndrops_gf_qf_sim.sum(axis=1)
		ndrops_gf_real_mbins = ndrops_gf_qf_real[0].sum(axis=1)
		plt.plot(mbins[:-1], ndrops_gf_sim_mbins, drawstyle='steps-post', 
			label='simulation', linestyle='--', color=colors[0])
		plt.plot(mbins[:-1], ndrops_gf_real_mbins, drawstyle='steps-post',
			label='sample', linestyle='-', color=colors[1])
		# plot Poissonian error bars for the number of real dropouts
		plt.errorbar(mbins[:-1]+dm, ndrops_gf_real_mbins, yerr=sqrt(ndrops_gf_real_mbins), 
			fmt=None, ecolor=colors[1], capsize=8)
		plt.xticks(mbins[:-1][::2])
		plt.xlabel('GALFIT magnitude')
		yticks = plt.yticks()[0]
		plt.yticks(yticks[::2])
		plt.title(fields[i])
		plt.xlim(mbins[0], mbins[-2])
		#if i == 0:
		plt.ylabel('Number with good fits')
		plt.legend(loc=2, prop=fp)
		# ******************************************** #
		plt.subplot(2, 2, i+3)  # Re bins in GOODS/HUDF
		ndrops_gf_sim_rbins = ndrops_gf_qf_sim.sum(axis=0)
		ndrops_gf_real_rbins = ndrops_gf_qf_real[0].sum(axis=0)
		plt.plot(rbins[:-1], ndrops_gf_sim_rbins, drawstyle='steps-post',
			label='simulation', linestyle='--', color=colors[0])
		plt.plot(rbins[:-1], ndrops_gf_real_rbins, drawstyle='steps-post',
			label='sample', linestyle='-', color=colors[1])
		# plot Poissonian error bars for the number of real dropouts
		plt.errorbar(rbins[:-1]+dlogr, ndrops_gf_real_rbins, 
			yerr=sqrt(ndrops_gf_real_rbins), fmt=None, ecolor=colors[1], capsize=8)
		plt.xticks(rbins[:-1][::2])
		plt.xlabel('GALFIT log10(Re) [pixels]')
		yticks = plt.yticks()[0]
		plt.yticks(yticks[::2])
		plt.xlim(rbins[0], rbins[-2])
		#if i == 0:
		plt.ylabel('Number with good fits')
		plt.legend(loc=2, prop=fp)
		n_gf += [ndrops_gf[0]]
		n_gf_qf_sim += [ndrops_gf_qf_sim]
		n_gf_qf_real += [ndrops_gf_qf_real[0]]
	#plt.rcdefaults()
	return n_gf, n_gf_qf_sim, n_gf_qf_real
	

def show_galfit_qf_comp(c_sim, c_drop, limits=array([[22.0, 26.5], [-1.0, 2.0]]), \
   maglim=26.5, chisqnulim=0.4, cmap=mpl.cm.gray):
   # calculate the fraction of sources thrown away due to poor GALFIT quality
   # first, only pick the bins that cover the real dropouts
   fig = plt.figure(figsize=(7,10))
   palette = plt.cm.gray
   palette.set_over('cyan', 1.0)
   palette.set_under('red', 1.0)
   palette.set_bad('black', 1.0)
   sp1 = fig.add_axes([0.1,0.65,0.8,0.25])
   sp2 = fig.add_axes([0.1,0.375,0.8,0.25])
   sp3 = fig.add_axes([0.1,0.10,0.8,0.25])
   bw = array([0.5, 0.2])
   nbins = (limits[:,1] - limits[:,0]) / bw
   nbins = nbins.astype('int')
   mlims = arange(limits[0,0], limits[0,1]+bw[0], bw[0])
   rlims = arange(limits[1,0], limits[1,1]+bw[1], bw[1])
   binedge = [mlims, rlims]
   qf = (c_drop.reout_err/c_drop.reout<=0.6)&(c_drop.chisqnu<=chisqnulim)&(c_drop.magout<=maglim)
   #ndrops = zeros(nbins, 'int')
   ndrops_tot = histogram2d(c_drop.magout, log10(c_drop.reout), bins=binedge)[0]
   ndrops_qf  = histogram2d(compress(qf,c_drop.magout), compress(qf,log10(c_drop.reout)),
      bins=binedge)[0]
   #ndrops_tot = maximum(ndrops_tot, 1)
   comp_qf = ndrops_qf.astype('float') / maximum(ndrops_tot.astype('float'), 1.)
   nsim_tot = histogram2d(compress(c_sim.recovered_galfit,c_sim.magout),\
      compress(c_sim.recovered_galfit,log10(c_sim.reout)), bins=binedge)[0]
   nsim_tot = maximum(nsim_tot, 1)
   qf_sim = (c_sim.reouterr/c_sim.reout<=0.6)&(c_sim.chisqnu<=chisqnulim)&(c_sim.recovered_galfit==1)
   nsim_qf  = histogram2d(compress(qf_sim,c_sim.magout),compress(qf_sim,log10(c_sim.reout)),
      bins=binedge)[0]
   #place(nsim_qf, ndrops_tot<1, 0)
   comp_sim_qf = nsim_qf.astype('float') / nsim_tot.astype('float')

   x1 = showmap(comp_qf, limits=limits, cmap=cmap, newfig=False, fig=fig, ax=sp1,
      vrange=[0.,1.])
   #x1[2].set_xlabel('GALFIT mag_out')
   x1[2].set_ylabel('log10(Re_out) [pixel]')
   x1[2].set_xticklabels([])  # unset the tick labels
   x1[3].set_ticks([0.0,0.25,0.5,0.75,1.0])
   x2 = showmap(comp_sim_qf, limits=limits, cmap=cmap, newfig=False, fig=fig, ax=sp2,
      vrange=[0.,1.])
   #x2[2].set_xlabel('GALFIT mag_out')
   x2[2].set_ylabel('log10(Re_out) [pixel]')
   x2[2].set_xticklabels([])
   x2[3].set_ticks([0.0,0.25,0.5,0.75,1.0])
   comp_sim_frac = comp_qf / maximum(comp_sim_qf, 1.e-5)
   putmask(comp_sim_frac, comp_qf<1.e-5, 0.)
   comp_sim_frac_masked = ma.masked_where(comp_sim_frac<0.05, comp_sim_frac)
   print max(comp_sim_frac.ravel()), min(comp_sim_frac.ravel())
   x3 = showmap(comp_sim_frac_masked, limits=limits, cmap=palette, newfig=False, fig=fig,
      ax=sp3, vrange=[0.8,1.2], mask=comp_sim_frac_masked.mask, cbar_extend='both')
   x3[2].set_xlabel('GALFIT mag_out')
   x3[2].set_ylabel('log10(Re_out) [pixel]')
   x3[3].set_ticks([0.8,1.0,1.2])
   
   return x1, x2, x3




def plot_det_comp(c, limits, bw=[0.5,0.2], pixdx=array([0.02,0.02])):
	# calculate the SExtractor detection completeness in z-band
	nbins = (limits[:,1]-limits[:,0]) / bw
	nbins = nbins.astype('int')
	m_within = (c.z_magin>=limits[0,0]) & (c.z_magin<limits[0,1])
	r_within = (log10(c.re_input)>=limits[1,0]) & (log10(c.re_input)<limits[1,1])
	within = m_within & r_within
	ninput_bins = histogram2d(compress(within,c.z_magin), 
		compress(within,log10(c.re_input)), bins=nbins)[0]
	#print ninput_bins
	ndet_bins = histogram2d(compress((c.detect==1)&within,c.z_magin),
		compress((c.detect==1)&within,log10(c.re_input)), 
		bins=nbins)[0]
	print ndet_bins
	det_comp = where(ninput_bins>0, 
		array(ndet_bins).astype('float') / array(ninput_bins).astype('float'),
		0.)
	showmap(det_comp, limits, bw=bw, pixdx=pixdx)
	return det_comp
	
def plot_gf_det_comp(c, limits, bw=array([0.5,0.2]), pixdx=array([0.02,0.02])):
	# calculate the SExtractor detection completeness in z-band
	nbins = (limits[:,1]-limits[:,0]) / bw
	nbins = nbins.astype('int')
	m_within = (c.magin>=limits[0,0]) & (c.magin<limits[0,1])
	r_within = (log10(c.rein)>=limits[1,0]) & (log10(c.rein)<limits[1,1])
	within = m_within & r_within
	ninput_bins = histogram2d(compress(within,c.magin), 
		compress(within,log10(c.rein)), bins=nbins)[0]
	#print ninput_bins
	ndet_bins = histogram2d(compress((c.recovered_se==1)&within,c.magin),
		compress((c.recovered_se==1)&within,log10(c.rein)), 
		bins=nbins)[0]
	print ndet_bins
	det_comp = where(ninput_bins>0, 
		array(ndet_bins).astype('float') / array(ninput_bins).astype('float'),
		0.)
	showmap(det_comp, limits, bw=bw, pixdx=pixdx)
	return det_comp
	
def plot_galfit_qf_26(kgrid1, kgrid2):
	# compare the GALFIT completeness near mag=26.0 (between 25.5 and 26.5 mag)
	# in GOODS-depth and UDF-depth images (in i-band or z-band)
	# kgrid1 is the kernels in the GOODS-depth image; kgrid2 is the kernels in the UDF-depth image
	ninput_tot1 = 0; ninput_tot2 = 0
	ngood_tot1 = 0; ngood_tot2 = 0
	# first the GOODS-depth kernels
	for k in kgrid1.kernels:
		k1 = kgrid1.kernels[k]
		if (round(k1.x0_bin[0],1)==25.5) or (round(k1.x1_bin[0],1)==26.5):
			ninput_tot1 += k1.ninput_bin
			ngood_tot1 += k1.ngood_bin
	# then the UDF-depth kernels
	for k in kgrid2.kernels:
		k2 = kgrid2.kernels[k]
		if (round(k2.x0_bin[0],1)==25.5) or (round(k2.x1_bin[0],1)==26.5):
			ninput_tot2 += k2.ninput_bin
			ngood_tot2 += k2.ngood_bin
	# output
	print "GOODS-depth kernel: n_input=%d, n_good=%d (%.1f%%)" % (ninput_tot1, ngood_tot1, float(ngood_tot1)/float(ninput_tot1)*100.)
	print "UDF-depth kernel: n_input=%d, n_good=%d (%.1f%%)" % (ninput_tot2, ngood_tot2, float(ngood_tot2)/float(ninput_tot2)*100.)
	return ninput_tot1, ngood_tot1, ninput_tot2, ngood_tot2