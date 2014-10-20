#!/usr/bin/env python

import numpy as N
import matplotlib.pyplot as plt
from pygoods import Ftable, sextractor
import pyfits
from stats import running

class plot_gfsim(running.running):
	def __init__(self,catalog):
		self.c = Ftable(catalog)
		self.gfmeasure = (self.c.recovered_galfit==1) & (self.c.re_out_err/self.c.re_out<=0.6) & (self.c.mag_auto<99.0) & (self.c.magerr_auto>0.)
		self.dmag = (self.c.mag_out-self.c.mag_in)[self.gfmeasure]
		self.dlogRe = (N.log10(self.c.re_out)-N.log10(self.c.re_in))[self.gfmeasure]
		self.mag_auto = self.c.mag_auto[self.gfmeasure]
		self.magerr_auto = self.c.magerr_auto[self.gfmeasure]
		self.SN_auto = (1.0857/self.c.magerr_auto)[self.gfmeasure]
		self.mag_in = (self.c.mag_in)[self.gfmeasure]
		self.re_in = (self.c.re_in)[self.gfmeasure]
		self.mag_out = (self.c.mag_out)[self.gfmeasure]
		self.re_out = (self.c.re_out)[self.gfmeasure]
		self.robust = False
		self.fig = None
		self.ax = None

	def newfig(self):
		self.fig = plt.figure()
		self.ax = self.fig.add_subplot(111)

	def newax(self):
		self.ax = self.fig.add_subplot(111)

	def SN_dmag(self,nsamp=50,interval=50,sn_floor=3.0,robust=True,scatter=True):
		running.running.__init__(self,self.SN_auto[self.SN_auto>=sn_floor],self.dmag[self.SN_auto>=sn_floor],
			nsamp=nsamp,interval=interval)
		if self.ax == None:
			self.newfig()
		if scatter:
			self.ax.scatter(self.xval,self.yval,marker='x',s=2**2,c='0.5')
		if robust:
			if self.robust==False:
				self.calc_robust()
				self.robust = True
			self.ax.errorbar(self.x_bins,self.loc_robust_bins,yerr=self.scale_robust_bins,marker='^')
		else:
			self.mean()
			self.std()
			self.ax.errorbar(self.x_bins,self.mean_bins,yerr=self.std_bins)
		self.ax.set_xlabel('S/N auto',size=14)
		self.ax.set_xscale('log')
		self.ax.set_ylabel(r'$\Delta$mag',size=14)
		self.ax.plot([sn_floor,max(self.SN_auto)],[0.,0.],'--',c='black',lw=1.2)
		self.ax.set_xlim(sn_floor,max(self.x_bins+30))
		self.ax.set_title(self.c.filename+' robust='+str(robust),size=18)

	def SN_dlogre(self,nsamp=50,interval=50,sn_floor=3.0,robust=True,scatter=True):
		running.running.__init__(self,self.SN_auto[self.SN_auto>=sn_floor],self.dlogRe[self.SN_auto>=sn_floor],
			nsamp=nsamp,interval=interval)
		if self.ax == None:
			self.newfig()
		if scatter:
			self.ax.scatter(self.xval,self.yval,marker='x',s=2**2,c='0.5')
		if robust:
			if self.robust==False:
				self.calc_robust()
				self.robust = True
			self.ax.errorbar(self.x_bins,self.loc_robust_bins,yerr=self.scale_robust_bins,marker='^')
		else:
			self.mean()
			self.std()
			self.ax.errorbar(self.x_bins,self.mean_bins,yerr=self.std_bins)
		self.ax.set_xlabel('S/N auto',size=14)
		self.ax.set_xscale('log')
		self.ax.set_ylabel(r'$\Delta\log(R_e)$',size=14)
		self.ax.plot([sn_floor,max(self.SN_auto)],[0.,0.],'--',c='black',lw=1.2)
		self.ax.set_xlim(sn_floor,max(self.x_bins+30))
		self.ax.set_title(self.c.filename+' robust='+str(robust),size=18)

	def SN_dmag_scale(self,nsamp=50,interval=50,sn_floor=3.0,robust=True,label=""):
		running.running.__init__(self,self.SN_auto[self.SN_auto>=sn_floor],self.dmag[self.SN_auto>=sn_floor],
			nsamp=nsamp,interval=interval)
		if self.ax == None:
			self.newfig()
		if robust==True:
			self.calc_robust()
			self.robust = True
			self.ax.plot(self.x_bins,self.scale_robust_bins,marker='^',label=label	)
		else:
			self.mean()
			self.std()
			self.ax.plot(self.x_bins,self.std_bins,marker='^')
		self.ax.set_xlabel('S/N auto',size=14)
		self.ax.set_xscale('log')
		self.ax.set_ylabel(r'$\sigma(\Delta\mathrm{mag})$',size=14)
		self.ax.plot([sn_floor,max(self.SN_auto)],[0.,0.],'--',c='black',lw=1.2)
		self.ax.set_xlim(sn_floor,max(self.x_bins+30))
		self.ax.set_title(self.c.filename+' robust='+str(robust),size=18)

	def SN_dlogre_scale(self,nsamp=50,interval=50,sn_floor=3.0,robust=True,label=""):
		running.running.__init__(self,self.SN_auto[self.SN_auto>=sn_floor],self.dlogRe[self.SN_auto>=sn_floor],
			nsamp=nsamp,interval=interval)
		if self.ax == None:
			self.newfig()
		if robust:
			self.calc_robust()
			self.robust = True
			self.ax.plot(self.x_bins,self.scale_robust_bins,marker='^',label=label	)
		else:
			self.mean()
			self.std()
			self.ax.plot(self.x_bins,self.std_bins,marker='^')
		self.ax.set_xlabel('S/N auto',size=14)
		self.ax.set_xscale('log')
		self.ax.set_ylabel(r'$\sigma[\Delta\log(R_e)]$',size=14)
		self.ax.plot([sn_floor,max(self.SN_auto)],[0.,0.],'--',c='black',lw=1.2)
		self.ax.set_xlim(sn_floor,max(self.x_bins+30))
		self.ax.set_title(self.c.filename+' robust='+str(robust),size=18)