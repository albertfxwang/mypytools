#!/usr/bin/env python

from numpy import *
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from bivariate import bivRL as bl
from mpl_toolkits.axes_grid1.axes_grid import ImageGrid
from mpl_toolkits.axes_grid1 import make_axes_locatable
from bivariate import zdist, mlutil
from my_mplfonts import Helvetica
from pygoods import Ftable, sextractor, parseconfig
import cPickle
import copy
from bivariate import fit_bdrops as fb
from bivariate import fit_udrops as fu
from bivariate import fit_idrops as fi
import time
import yaml

mci = bl.mconvert('/Users/khuang/Dropbox/Research/bivariate/bivariate_fit/kcorr/M1500_to_f775w_omega_m_0.3.txt')
rfmt = FormatStrFormatter('%.3f')
ndic_candels = {'udf':1, 'deep':2, 'ers':3, 'wide':4}
ndic_goods = {'goods':1, 'udf':2}

def gfcrit_dropout(c, gfband, reerr_lim=0.6, chi2nu_lim=2.0):
   gfcrit = ((getattr(c,'%s_magflag'%gfband) == 0) & \
             (getattr(c,'%s_reflag'%gfband) == 0))
   gfcrit = (gfcrit & (getattr(c,'%s_nflag'%gfband) == 0))
   gfcrit = (gfcrit & (getattr(c,'%s_reout_err_gf'%gfband)/getattr(c,'%s_reout_gf'%gfband)<=reerr_lim))
   gfcrit = (gfcrit & (getattr(c,'%s_chi2nu_gf'%gfband) <= chi2nu_lim))
   vfcrit = (c.vflag != 10) & (c.vflag != 4)
   gfcrit = gfcrit & vfcrit
   return gfcrit

class plot_bivariate_RL(object):
   """
   A class for plotting bivariate size-luminosity distribution fitting 
   results. This should be a parent class for more specialized classes.
   """
   # def __init__(self, parfile, field, drop, gfband, nproc_model=2,
   #              z_mean=3.2, expand=False):
   #    """
   #    Consider getting rid of self.__init__ function...
   #    """
   #    bl.bivariate_RL_class.__init__(self)
   #    #def __init__(self, parfile, field, nproc_model=2):
   #    if drop in ['f435w','f606w']:
   #       nfield = ndic_goods[field.lower()]
   #    else:
   #       nfield = ndic_candels[field.lower()]
   #    # c = parseconfig(parfile)
   #    c = yaml.load(open(parfile))
   #    self.limits = array(c['LIMITS%d'%nfield]).reshape((2,2))
   #    self.pixdx = array(c['PIXDX'])
   #    self.npix = (self.limits[:,1] - self.limits[:,0]) / self.pixdx
   #    self.npix = around(self.npix).astype('int')
   #    if c.has_key('MCFILE1'):
   #       self.mcfile = c['MCFILE%d'%nfield]
   #    else:
   #       self.mcfile = c['MCFILE']
   #    self.zdgridfile = c['ZDGRIDFILE%d'%nfield]
   #    self.kgridfile = c['KGRIDFILE%d'%nfield]
   #    self.sdfile = c['SDFILE%d'%nfield]
   #    self.intfracfile = c['INTFRACFILE%d'%nfield]
   #    self.drop = drop
   #    self.field = field
   #    self.mag_lolim = c['MAG_LOLIMS'][nfield-1]
   #    self.nproc_model = nproc_model
   #    self.gfband = gfband
   #    self.npix = around((self.limits[:,1]-self.limits[:,0])/self.pixdx)
   #    self.npix = self.npix.astype('int')
   #    self.mc = bl.mconvert(self.mcfile)
   #    self.zdgrid = cPickle.load(open(self.zdgridfile))
   #    self.kgrid = cPickle.load(open(self.kgridfile))
   #    self.kgrid.interpolate_completeness(self.limits, self.pixdx)
   #    self.z_mean = z_mean
   #    self.M0 = c['M0']
   #    self.expand = expand
   #    if 'REERR_LIM' in c.keys():
   #       self.reerr_lim = c['REERR_LIM']
   #    else:
   #       self.reerr_lim = 0.6
   #    if 'CHI2NU_LIM' in c.keys():
   #       self.chi2nu_lim = c['CHI2NU_LIM']
   #    else:
   #       self.chi2nu_lim = 1000.
   #    if len(self.sdfile) > 0:
   #       self.sd = cPickle.load(open(self.sdfile))
   #    else:
   #       self.sd = None
   #    if len(self.intfracfile) > 0:
   #       self.intfrac = Ftable(self.intfracfile)
   #    else:
   #       self.intfrac = None
   #    if expand:
   #      limits2 = self.limits.copy()
   #      limits2[0][1] += 0.5
   #      bl.Pijk_dVdz(self.zdgrid, self.mc, limits2, self.pixdx)
   #    else:
   #      bl.Pijk_dVdz(self.zdgrid, self.mc, self.limits, self.pixdx)

   def plot_rawmodel(self, par, ax=None, colorbar=True, show_LF=False,
                     show_sizedist=False):
      self.bivariate_RL(par, self.limits, self.pixdx, self.drop, 
                        self.field, mc=self.mc, add_interloper=False,
                        nproc_model=self.nproc_model)
      self.raw_model = self.model.copy()
      self.vmin = 0.
      self.vmax = max(self.raw_model.ravel())
      self.npix = shape(self.raw_model)
      if ax==None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
      cax = ax.imshow(self.raw_model.swapaxes(0,1), origin='lower',
                extent=(0,self.npix[0],0,self.npix[1]), vmin=self.vmin,
                vmax=self.vmax)
      if colorbar:
         plt.colorbar(cax)
      #plt.colorbar(cax)
      ax.set_xlim(0, self.npix[0])
      ax.set_xticks(arange(0, self.npix[0], 50))
      ax.set_xticklabels(map(lambda x:'%.1f'%x, 
                         arange(self.limits[0][0], self.limits[0][1], 1.)),
                         font_properties=Helvetica(14))
      ax.set_yticks(arange(0.,self.npix[1], 20)+10.)
      ax.set_yticklabels(map(lambda y:'%.1f'%y, 
                     arange(self.limits[1][0], self.limits[1][1], 0.4)+0.2),
                     font_properties=Helvetica(14))
      ax.set_xlabel('$M_{1500}$', font_properties=Helvetica(18))
      ax.set_ylabel('log(Re) [pixels]', font_properties=Helvetica(18))
      #fig.show()
      #return raw_model
      divider = make_axes_locatable(ax)
      if show_LF:
            axLF = divider.append_axes("top", size=1.2, pad=0.0, sharex=ax)
            axLF.semilogy(arange(self.limits[0,0],self.limits[0,1],self.pixdx[0]), 
                          self.raw_model.sum(axis=1), color='black', 
                          nonposy='mask', lw=2.0)
            axLF.set_ylim(1.e-2, 1.e2)
            axLF.set_yticks([1.e-2, 1., 1.e2])


   def plot_zdgridmodels(self, par, ax=None):
      """ 
      Model with P(z) corrections.
      """
      self.bivariate_RL(par, self.limits, self.pixdx, self.drop, self.field,
                        zdgrid=self.zdgrid, mc=self.mc, add_interloper=False,
                        nproc_model=self.nproc_model)
      
      if ax==None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
      self.zdgridmodel = self.model.copy()
      self.vmin = 0.
      self.vmax = max(self.zdgridmodel.ravel())
      cax = ax.imshow(self.model.swapaxes(0,1), origin='lower',
                extent=(0,self.npix[0],0,self.npix[1]), vmin=self.vmin,
                vmax=self.vmax)
      plt.colorbar(cax)
      ax.set_xlim(0, self.npix[0])
      ax.set_xticks(arange(0, self.npix[0], 50))
      ax.set_xticklabels(map(lambda x:'%.1f'%x, 
                         arange(self.limits[0][0], self.limits[0][1], 1.)),
                         font_properties=Helvetica(14))
      ax.set_yticks(arange(0.,self.npix[1], 20)+10.)
      ax.set_yticklabels(map(lambda y:'%.1f'%y,
                     arange(self.limits[1][0], self.limits[1][1], 0.4)+0.2),
                     font_properties=Helvetica(14))
      ax.set_xlabel('$M_{1500}$', font_properties=Helvetica(18))
      ax.set_ylabel('log(Re) [arcsec]', font_properties=Helvetica(18))

   def plot_interlopermodels(self, par, interlopers_only=True, ax=None):
      if ax==None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
      self.bivariate_RL(par, self.limits, self.pixdx, self.drop, self.field,
                        zdgrid=self.zdgrid, mc=self.mc, add_interloper=True,
                        sd=self.sd, intfrac=self.intfrac, 
                        mag_lolim=self.mag_lolim, nproc_model=self.nproc_model)
      if interlopers_only:
         intmodel = self.add_interlopers(self.intfrac, self.sd, 
                                         self.mag_lolim)
         cax = ax.imshow(intmodel.model.swapaxes(0,1), origin='lower',
                   extent=(0,self.npix[0],0,self.npix[1]))
         plt.colorbar(cax)
      else:         
         cax = ax.imshow(self.model.swapaxes(0,1), origin='lower',
                extent=(0,self.npix[0],0,self.npix[1]))
         plt.colorbar(cax)
      ax.set_xlim(0, self.npix[0])
      ax.set_xticks(arange(0, self.npix[0], 50))
      ax.set_xticklabels(map(lambda x:'%.1f'%x, 
                         arange(self.limits[0][0], self.limits[0][1], 1.)),
                         font_properties=Helvetica(14))
      ax.set_yticks(arange(0.,self.npix[1], 20)+10.)
      ax.set_yticklabels(map(lambda y:'%.1f'%y,
                     arange(self.limits[1][0], self.limits[1][1], 0.4)+0.2),
                     font_properties=Helvetica(14))
      ax.set_xlabel('$M_{1500}$', font_properties=Helvetica(18))
      ax.set_ylabel('log(Re) [arcsec]', font_properties=Helvetica(18))

   def plot_fullmodel(self, par, field, ax=None, add_interloper=False, 
                      mag0=22.0, mag1=26.5, cbar=True, **imshow_kw):
      if ax==None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
      model = self.models[field]
      npix0 = copy.deepcopy(self.npix[field])
      new_npix = [int(round((mag1-mag0)/self.PIXDX[0])), self.npix[field][1]]
      new_model = zeros(new_npix)
      i0 = (self.limits[field][0][0] - mag0) / self.PIXDX[0]
      i0 = int(i0)
      i1 = (self.limits[field][0][1] - mag0) / self.PIXDX[0]
      i1 = int(i1)
      #if sum(self.zdgrid.Pk.ravel()) == 0.:
      #   bl.Pijk_dVdz(self.zdgrid, self.mc, self.limits, self.pixdx)
      model.bivariate_RL(par, self.limits[field], self.pixdx, self.DROP, field,
                        zdgrid=self.zdgrids[field], mc=self.mc[field], 
                        add_interloper=self.ADD_INTERLOPER,
                        sd=self.sd[field], intfrac=self.intfrac[field], 
                        mag_lolim=self.mag_lolims[field], 
                        kgrid=self.tfkgrids[field], nproc_model=self.nproc_model,
                        expand=self.expand)
      #fullmodel = self.model
      new_model[i0:i1,:] = model.model.copy()
      if 'extent' not in imshow_kw.keys():
        imshow_kw['extent'] = (0,new_npix[0],0,new_npix[1])
      cax = ax.imshow(new_model.swapaxes(0,1), **imshow_kw)
      if cbar==True:
         plt.colorbar(cax)
      #ax.set_xlim(0, npix0[0])
      #ax.set_xticks(arange(0, npix0[0]+50, 50))
      ax.set_xticks(arange(0, new_npix[0]+50, 50))
      ax.set_xticklabels(map(lambda x:'%.1f'%x, 
                         #arange(self.limits[0][0], self.limits[0][1]+1., 1.)),
                         arange(mag0, mag1+1., 1.)),
                         font_properties=Helvetica(14))
      ax.set_yticks(arange(0.,npix0[1], 20)+10.)
      ax.set_yticklabels(map(lambda y:'%.1f'%y,
                     arange(self.limits[field][1][0], self.limits[field][1][1], 0.4)+0.2),
                     font_properties=Helvetica(14))
      ax.set_xlabel('$M_{1500}$', font_properties=Helvetica(18))
      ax.set_ylabel('log(Re) [arcsec]', font_properties=Helvetica(18))
      ax.set_xlim(0, new_npix[0])
      ax.set_ylim(0, new_npix[1])
      return ax

   def plot_step1(self, par, figsize=(8,10)):
      """
      Plot raw model + model with P(z) corrections.
      """
      fig = plt.figure(figsize=figsize)
      ax1 = fig.add_subplot(2,1,1)
      self.plot_rawmodel(par, ax=ax1)
      ax2 = fig.add_subplot(2,1,2)
      self.plot_zdgridmodels(par, ax=ax2)
      return 0

   def plot_step2(self, par, figsize=(8,10)):
      """
      Plot interloper models and model with P(z) corrections & interlopers.
      """
      fig = plt.figure(figsize=figsize)
      ax1 = fig.add_subplot(2,1,1)
      self.plot_interlopermodels(par, interlopers_only=True, ax=ax1)
      ax2 = fig.add_subplot(2,1,2)
      self.plot_interlopermodels(par, interlopers_only=False, ax=ax2)
      return 0
      
   def plot_step3(self, par, figsize=(8,10), mag0=22.0, mag1=26.5):
      """
      Plot the model with P(z) & interlopers and the model also with GALFIT
      transfer function kernels.
      """
      fig = plt.figure(figsize=figsize)
      ax1 = fig.add_subplot(2,1,1)
      self.plot_interlopermodels(par, ax=ax1, interlopers_only=False)
      ax2 = fig.add_subplot(2,1,2)
      self.plot_fullmodel(par, ax=ax2, mag0=mag0, mag1=mag1)

   def plot_data(self, field, mag0, ax=None, color='white', 
                 markersize=5**2):
      """
      ci is the FITS table containing the GALFIT measurements of all dropouts.
      It is assumed that attribute self.datasets exist for each field.
      """
      if ax==None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
      # crit = (getattr(ci,'%s_fit'%field)==True)
      # crit = crit & (getattr(ci,'%s_magout_gf'%self.gfband)<=self.limits[0][1])
      # mag_out = getattr(ci,'%s_magout_gf'%self.gfband)[crit]
      # logre_out = log10(getattr(ci,'%s_reout_gf'%self.gfband))[crit]
      mag_out = self.datasets[field][0]
      logre_out = self.datasets[field][1]
      x = (mag_out - mag0) / self.pixdx[0]
      y = (logre_out - self.limits[field][1][0]) / self.pixdx[1]
      ax.scatter(x, y, s=markersize, color=color, label=field)
      plt.draw()
      return ax

   def plot_schechterLF(self, ci, params, phistar, z0, z1, field, ax=None, 
                        add_interloper=False,
                        binwidth=0.5, marker='o', mcolor='black'):
      """
      Plot Schechter function of the convolved model.
      """
      if ax == None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
      recalc_model = 1
      #if (hasattr(self,'model')==False): 
      #   recalc_model = 1
      #elif (sum(self.model.ravel()) <= 0.):
      #   recalc_model = 1
      # Now plot the number of points
      ax.set_yscale('log')
      #crit = (getattr(ci,'%s_fit'%self.field)==True)
      # crit = (getattr(ci,'%s_magout_gf'%self.gfband)<=self.limits[0][1])
      # crit = crit & (getattr(ci,self.field)==True)
      # gfcrit = gfcrit_dropout(ci, self.gfband, self.reerr_lim, self.chi2nu_lim)
      # mag_out = getattr(ci,'%s_magout_gf'%self.gfband)[crit&gfcrit]
      mag_out = self.datasets[field][0]
      mag_bins = arange(self.limits[field][0][0], 
                        self.limits[field][0][1]+binwidth, 
                        binwidth)
      Nobs = histogram(mag_out, bins=mag_bins)[0] / binwidth
      Nobs_err = sqrt(Nobs)
      Nobs_high = Nobs + Nobs_err
      Nobs_low = maximum(Nobs - Nobs_err, 1.e-10)
      Nobs_err_low = Nobs - Nobs_low
      mask = (Nobs==0)
      ax.errorbar((mag_bins[:-1]+binwidth/2.)[mask==False], 
                  Nobs[mask==False], 
                  yerr=[Nobs_err_low[mask==False],Nobs_err[mask==False]], 
                  fmt=marker, ms=10, mfc='none', mec=mcolor,
                  ecolor=mcolor, capsize=8, capthick=1.3, elinewidth=1.3,
                  label=field.upper())
      if recalc_model:
         # calculate uncorrected model
         #model0 = self.model_z(params, self.z_mean, None, self.mc)
         model = self.models[field]
         model.bivariate_RL(params, self.limits[field], self.pixdx, self.drop, 
                        field, zdgrid=self.zdgrids[field], mc=self.mc[field], 
                        add_interloper=self.ADD_INTERLOPER,
                        sd=self.sd[field], intfrac=self.intfrac[field], 
                        mag_lolim=self.mag_lolims[field], 
                        kgrid=self.tfkgrids[field], nproc_model=self.nproc_model, 
                        delta_z=self.delta_z, expand=self.expand)
         model_full = model.model.copy()
         # The normalization should be 
         # LF = phistar * sum(model_full) * pixdx[1]
         print "z0, z1", z0, z1
         i0 = searchsorted(self.zdgrids[field].zarr, z0)
         i1 = searchsorted(self.zdgrids[field].zarr, z1)
         #V0 = sum(self.zdgrid.dVdz[i0:i1])
         LF0 = zeros(shape(model_full)[0])
         for i in range(i0, i1):
            #if self.zdgrid.Pz_tot[i] >= 0.01:
            zi = self.zdgrids[field].zarr[i] + self.zdgrids[field].dz / 2.
            mz = self.models[field].model_z(params, zi, None, self.mc[field])
            LF0 = LF0 + mz.sum(axis=1) * self.zdgrids[field].dVdz[i]
         LF0 = LF0 * phistar * self.pixdx[1]
         #LF0 = model0.sum(axis=1) * self.pixdx[1] * V0 * phistar
         LF = model_full.sum(axis=1) * self.pixdx[1] #* self.zdgrid.Veff 
         LF = LF * phistar
         #normfactor = sum(model_full.ravel()) * self.pixdx[1]
         #normfactor = len(mag_out) / normfactor
         #LF = model_full.sum(axis=1) * normfactor
         #LF = model_full.sum(axis=1) / sum(self.zdgrid.dVdz)
         # normalize by the total effective volume
         mag_arr = arange(self.limits[field][0][0], self.limits[field][0][1], 
                          self.pixdx[0])
         ax.plot(mag_arr, LF, label='corrected LF', lw=2.0, color=mcolor)
         ax.plot(mag_arr, LF0, label='uncorrected LF', color='0.2', lw=1.0)
         ax.set_xlabel('%s magnitude' % self.gfband.upper(), 
                       font_properties=Helvetica(16))
         ax.set_ylabel('Number per magnitude', font_properties=Helvetica(16))
      ax.set_ylim(ymin=min(model_full.sum(axis=1) / sum(self.zdgrids[field].dVdz)))
      ax.legend(loc=4, prop=Helvetica(13))
      return ax   

   def plot_sizedist(self, ci, params, phistar, z0, z1, field, ax=None, 
                     add_interloper=True,
                     binwidth=0.2, marker='o', mcolor='black'):
      """
      Plot the log-normal size distribution.
      """
      if ax == None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
      ax.set_yscale('log')
      # crit = (getattr(ci,'%s_fit'%field)==True)
      # logre_out = log10(getattr(ci,'%s_reout_gf'%self.gfband.lower())[crit])
      logre_bins = arange(self.limits[field][1][0], 
                          self.limits[field][1][1]+binwidth,
                          binwidth)
      Nobs = histogram(self.datasets[field][1], bins=logre_bins)[0] / binwidth
      Nobs_err = sqrt(Nobs)
      Nobs_high = Nobs + Nobs_err
      Nobs_low = maximum(Nobs - Nobs_err, 1.e-10)
      Nobs_err_low = Nobs - Nobs_low
      mask = (Nobs==0)
      ax.errorbar((logre_bins[:-1]+binwidth/2.)[mask==False],
                  Nobs[mask==False],
                  yerr=[Nobs_err_low[mask==False],Nobs_err[mask==False]],
                  fmt=marker, ms=10, mfc='none', mec=mcolor,
                  ecolor=mcolor, capsize=8, capthick=1.3, elinewidth=1.3,
                  label=field.upper())
      print "z0, z1", z0, z1
      logr_arr = arange(self.limits[field][1][0], self.limits[field][1][1], 
                        self.pixdx[1])
      #self.bivariate_RL(params, self.limits, self.pixdx, self.drop, 
      #                  self.field,
      #                  zdgrid=None, mc=self.mc, 
      #                  add_interloper=False,
      #                  sd=None, intfrac=None, 
      #                  mag_lolim=self.mag_lolim, kgrid=None,
      #                  nproc_model=self.nproc_model, delta_z=self.delta_z,
      #                  expand=self.expand)
      #SD0 = self.model.sum(axis=0) * Vtot * phistar * self.pixdx[0]
      # Now calculate the full corrected size distribution
      self.models[field].bivariate_RL(params, self.limits[field], self.pixdx, 
                        self.drop, field, zdgrid=self.zdgrids[field], 
                        mc=self.mc[field], add_interloper=self.add_interloper,
                        sd=self.sd[field], intfrac=self.intfrac[field], 
                        mag_lolim=self.mag_lolims[field], 
                        kgrid=self.tfkgrids[field], 
                        nproc_model=self.nproc_model, delta_z=self.delta_z,
                        expand=self.expand)
      SD = self.models[field].model.sum(axis=0) * self.pixdx[0] * phistar
      i0 = searchsorted(self.zdgrids[field].zarr, z0)
      i1 = searchsorted(self.zdgrids[field].zarr, z1)
      # Calculated un-corrected model at z=self.z_mean
      SD0 = zeros(shape(self.models[field].model)[1])
      for i in range(i0, i1):
         zi = self.zdgrids[field].zarr[i] + self.zdgrids[field].dz / 2.
         mz = self.models[field].model_z(params, zi, None, self.mc[field])
         SD0 = SD0 + mz.sum(axis=0) * self.zdgrids[field].dVdz[i]
      SD0 = SD0 * phistar * self.pixdx[1]
      # The normalization should be 
      # SD = phistar * sum(model_full) * pixdx[0]
      # Now start plotting 
      ax.plot(logr_arr, SD0, label='uncorrected', color='0.2', lw=1.0)
      ax.plot(logr_arr, SD, label='corrected', lw=2.0, color=mcolor)
      ax.set_xlabel('log10(Re) [pixels]', font_properties=Helvetica(16))
      ax.set_ylabel('Number per dex in log10(Re)', font_properties=Helvetica(16))
      ax.set_ylim(ymin=min(self.models[field].model.sum(axis=0) / sum(self.zdgrids[field].dVdz)))
      ax.legend(loc=4, prop=Helvetica(13))
      return ax

# def plot_udrops_4fields_goodss_RL(par, c_drop, figsize=(12,8),
#                                   add_interloper=False,
#                                   datacolor='white', plot_data=True,
#                                   parfile='udrops_goodss_fit_RL.par.yml',
#                                   markersize=5**2,
#                                   mag0=22.0, mag1=27.5, 
#                                   expand=True, **imshow_kw):
class plot_goodss_4fields(plot_bivariate_RL, fu.FitUdropsRL, fi.FitIdropsRL):
   """
   Plot 2D models for the parameters & the measured points in each depth
   in CANDELS GOODS-S: UDF, Deep, ERS, and Wide.
   This is the production function.
   """
   def __init__(self, c_drop, parfile, drop, figsize=(12,8),
                datacolor='white', plotdata=True,
                markersize=5**2,
                mag0=22.0, mag1=27.5,
                **imshow_kw):
      if drop == 'uvimos':
         fu.FitUdropsRL.__init__(self, parfile)
      elif drop == 'f775w':
         fi.FitIdropsRL.__init__(self, parfile)
      self.c = c_drop
      self.figsize = figsize
      self.datacolor = datacolor
      self.plotdata = plotdata
      self.markersize = markersize
      self.mag0 = mag0
      self.mag1 = mag1
      self.imshow_kw = imshow_kw
      self.npix = {}
      for i in range(4):
         f = self.FIELDS[i]
         limits = self.limits[f]
         self.npix[f] = around((limits[:,1]-limits[:,0])/self.PIXDX).astype('int')
   
      # fields = ['udf', 'deep', 'ers', 'wide']
      self.figures = []
      self.axes = {}
      if drop == 'uvimos':
         self.gfband = 'F606W'
      else:
         self.gfband = 'F125W'
   
   def plot_RL(self, par, cbar=False):
      """
      Plot the RL distribution along with data points used in the fit, 
      separately for each field.
      """
      fig = plt.figure(figsize=self.figsize)
      self.figures += [fig]
      self.axes['RL'] = []

      for i in range(4):
         # p = plot_bivariate_RL(parfile, fields[i], 'uvimos', 'f606w', 
         #                       expand=expand)
         # p = fu.FitUdropsRL(parfile)
         ax_i = fig.add_subplot(2, 2, i+1)
         self.axes['RL'] += [ax_i]
         #x0 = round((mag0 - p.limits[0][0]) / p.pixdx[0])
         #x0 = int(x0)
         #x1 = round((mag1 - p.limits[0][0]) / p.pixdx[0])
         #x1 = int(x1)
         #imshow_kw = dict(extent=(x0,x1,0,p.npix[1]))
         self.plot_fullmodel(par, self.FIELDS[i], ax=ax_i,
                             add_interloper=self.ADD_INTERLOPER,
                             mag0=self.mag0, mag1=self.mag1, 
                             cbar=cbar, **self.imshow_kw)
         if self.plotdata:
           self.plot_data(self.FIELDS[i], self.mag0, 
                          ax=ax_i, color=self.datacolor,
                          markersize=self.markersize)
         #ax_i.set_xlim(0, shape(p.model)[0])
         #ax_i.set_ylim(0, shape(p.model)[1])
         ax_i.text(0.15, 0.8, self.FIELDS[i].upper(), color='white', 
                   font_properties=Helvetica(16), transform=ax_i.transAxes,
                   horizontalalignment='center')
         ax_i.set_xlabel('%s magnitude' % self.gfband.upper(), 
                         font_properties=Helvetica(18))
      if not cbar:
         plt.subplots_adjust(wspace=0.15, hspace=0.15, left=0.05, right=0.975)
      parstr = "[ %.2f %.2f %.2f %.2f %.2f ]" % \
         (par[0],par[1],par[2],par[3],par[4])
      plt.suptitle("U-drops; "+parstr, font_properties=Helvetica(20))

# def plot_udrops_4fields_goodss_schechter_LF(par, cu, phistar, figsize=(12,8),
#                                      add_interloper=False, 
#                                      mag0=22.0, mag1=27.5, z0=3.0, z1=4.0,
#                                      parfile='udrops_goodss_fit_RL.par',
#                                      mcolors=['blue','green','red','black'],
#                                      delta_z=1.5):
   def plot_LF(self, par, phistar, z0=3.1, z1=4.1,
               mcolors=['blue','green','red','black']):
      """
      Plot the Schechter LF corrected for selection effects and measurement 
      errors.
      """      
      markers = ['o', '^', 'p', 's']
      fig = plt.figure(figsize=self.figsize)
      self.figures += [fig]
      self.axes['LF'] = []
      pms = []
      for i in range(4):
         field = self.FIELDS[i]
         ax = fig.add_subplot(2, 2, i+1)
         self.axes['LF'] += [ax]
         # p = plot_bivariate_RL(parfile, fields[i], 'uvimos', 'f606w')
         # pms += [p]
         # p.delta_z = delta_z
         self.plot_schechterLF(self.c, par, phistar, z0, z1, field, ax=ax, 
                            add_interloper=False,
                            binwidth=0.5, marker=markers[i], mcolor=mcolors[i])
         ax.text(0.1, 0.9, field.upper(), transform=ax.transAxes, ha='center',
                 va='center', font_properties=Helvetica(16))
         ax.set_ylim(ymin=1.e-1, ymax=10.**3.0)
         ax.set_xlim(self.mag0, self.mag1)
      #model0 = pms[0].bivariate_RL_M(par, pms[0].limits, pms[0].pixdx, 
      #                          M0=-21.0+pms[0].mc(6.0))
      #LF0 = model0.sum(axis=1)
      # mag_arr = arange(pms[0].limits[0][0], pms[0].limits[0][1], pms[0].pixdx[0])
      par_str = map(lambda x:'%.2f'%x, par)
      plt.suptitle('U-dropouts\n%s' % par_str, font_properties=Helvetica(20))
      #ax.plot(mag_arr, LF0, '--', lw=2.0, color='black')
      # plt.savefig('temp.pdf')

   def plot_SD(self, par, phistar, 
                     figsize=(12,9),
                     z0=3.1, z1=4.1,
                     mcolors=['blue','green','red','black']):
      """
      Plot the log-normal size distribution.
      """
      fields = ['udf', 'deep', 'ers', 'wide']
      markers = ['o', '^', 'p', 's']
      fig = plt.figure(figsize=self.figsize)
      self.figures += [fig]
      self.axes['sizedist'] = []
      for i in range(4):
         ax = fig.add_subplot(2, 2, i+1)
         self.axes['sizedist'] += [ax]
         # p = plot_bivariate_RL(parfile, fields[i], 'uvimos', 'f606w')
         # p.delta_z = delta_z
         # pms += [p]
         self.plot_sizedist(self.c, par, phistar, z0, z1, self.fields[i], 
                            ax=ax, add_interloper=self.add_interloper,
                            mcolor=mcolors[i])
         ax.text(0.1, 0.9, self.fields[i].upper(), transform=ax.transAxes, 
                 ha='center', va='center', font_properties=Helvetica(16))
         ax.set_ylim(ymin=1.e-1, ymax=10.**3.0)
         #ax.set_xlim(p.limits[1][0], p.limits[1][1])
         ax.set_xlim(-0.6, 1.8)
      par_str = map(lambda x:'%.2f'%x, par)
      plt.suptitle('U-dropouts\n%s' % par_str, font_properties=Helvetica(20))


def plot_idrops_4fields_goodss_RL(par, c_drop, figsize=(12,9),
                                  add_interloper=False,
                                  datacolor='white', plot_data=True,
                                  parfile='idrops_goodss_fit_RL.par',
                                  markersize=5**2,
                                  mag0=22.0, mag1=27.5, **imshow_kw):
   """
   Plot 2D models for the parameters & the measured points in each depth
   in CANDELS GOODS-S: UDF, Deep, ERS, and Wide.
   This is the production function.
   """
   fields = ['udf', 'deep', 'ers', 'wide']
   axes = []
   fig = plt.figure(figsize=figsize)
   for i in range(4):
      if fields[i]=='ers':
         gfband = 'f098m'
      else:
         gfband = 'f105w'
      p = plot_bivariate_RL(parfile, fields[i], 'f775w', gfband)
      ax_i = fig.add_subplot(2, 2, i+1)
      axes += [ax_i]
      p.plot_fullmodel(par, ax=ax_i, add_interloper=add_interloper,
                       mag0=mag0, mag1=mag1, **imshow_kw)
      if plot_data:
         p.plot_data(c_drop, fields[i], mag0, ax=ax_i, color=datacolor, 
                     markersize=markersize)
      #ax_i.set_xlim(0, shape(p.model)[0])
      #ax_i.set_ylim(0, shape(p.model)[1])
      ax_i.text(0.15, 0.8, fields[i].upper(), color='white', 
                font_properties=Helvetica(16), transform=ax_i.transAxes,
                horizontalalignment='center')
      if fields[i] == 'ers':
         ax_i.set_xlabel('F098M magnitude', font_properties=Helvetica(18))
      else:
         ax_i.set_xlabel('F105W magnitude', font_properties=Helvetica(18))
   plt.subplots_adjust(wspace=0.15, hspace=0.15, left=0.1, right=0.925)
   parstr = "[ %.2f %.2f %.2f %.2f %.2f ]" % \
      (par[0],par[1],par[2],par[3],par[4])
   plt.suptitle("i-drops; "+parstr, font_properties=Helvetica(20))
   return axes



def plot_idrops_4fields_goodss_schechter_LF(par, ci, phistar, figsize=(12,9),
                                     parfile='idrops_goodss_fit_RL.par',
                                     mag0=22.0, mag1=27.5,
                                     add_interloper=False,
                                     delta_z=1.5):
   """
   Plot the Schechter LF corrected for selection effects and measurement 
   errors.
   """
   fields = ['udf', 'deep', 'ers', 'wide']
   mcolors = ['blue', 'green', 'red', 'black']
   markers = ['o', '^', 'p', 's']
   fig = plt.figure(figsize=figsize)
   axes = []
   pms = []
   for i in range(4):
      if fields[i] == 'ers':
         gfband = 'f098m'
      else:
         gfband = 'f105w'
      ax = fig.add_subplot(2,2,i+1)
      axes += [ax]
      p = plot_bivariate_RL(parfile, fields[i], 'f775w', gfband)
      pms += [p]
      p.delta_z = delta_z
      p.plot_schechterLF(ci, par, phistar, 5.3, 6.8, ax=ax, 
                         add_interloper=add_interloper,
                        binwidth=0.5, marker=markers[i], mcolor=mcolors[i])
      ax.text(0.1, 0.9, fields[i].upper(), transform=ax.transAxes, ha='center',
              va='center', font_properties=Helvetica(16))
      ax.set_ylim(ymin=1.e-4)
   #model0 = pms[0].bivariate_RL_M(par, pms[0].limits, pms[0].pixdx, 
   #                          M0=-21.0+pms[0].mc(6.0))
   #LF0 = model0.sum(axis=1)
   mag_arr = arange(pms[0].limits[0][0], pms[0].limits[0][1], pms[0].pixdx[0])
   #ax.plot(mag_arr, LF0, '--', lw=2.0, color='black')
   par_str = map(lambda x:'%.2f'%x, par)
   plt.suptitle('i-dropouts\n%s' % par_str, font_properties=Helvetica(20))

def show_modelz(model, z0, z1, dz, t=1.0, vmax=1.e4):
   """
   Display the model contribution at each redshift. Flash the model for a given
   amount of time.
   """
   fig = plt.figure()
   ax = fig.add_subplot(111)
   cbar = 0
   for z in arange(z0, z1, dz):
      mz = model.modelz_all['%.1f'%z]
      cax = ax.imshow(mz.swapaxes(0,1), extent=(0,shape(mz)[0],0,shape(mz)[1]),
                      vmin=0., vmax=vmax)
      if cbar==0:
         plt.colorbar(cax)
         cbar = 1
      ax.set_title('z = %.1f' % z, size=18)
      plt.draw()
      time.sleep(t)

def plot_idrops_4fields_goodss_sizedist(par, ci, phistar, parfile, 
                                      figsize=(12,9),
                                      add_interloper=False, z0=5.3, z1=6.8,
                                      mcolors=['blue','green','red','black'],
                                      delta_z=1.5):
   """
   Plot the log-normal size distribution.
   """
   fields = ['udf', 'deep', 'ers', 'wide']
   markers = ['o', '^', 'p', 's']
   fig = plt.figure(figsize=figsize)
   axes = []
   pms = []
   for i in range(4):
      ax = fig.add_subplot(2, 2, i+1)
      axes += [ax]
      if fields[i] == 'ers':
         gfband = 'f098m'
      else:
         gfband = 'f105w'
      p = plot_bivariate_RL(parfile, fields[i], 'f775w', gfband)
      p.delta_z = delta_z
      pms += [p]
      p.plot_sizedist(ci, par, phistar, z0, z1, ax=ax, 
                    add_interloper=add_interloper,
                    mcolor=mcolors[i])
      ax.text(0.1, 0.9, fields[i].upper(), transform=ax.transAxes, ha='center',
              va='center', font_properties=Helvetica(16))
      ax.set_ylim(ymin=1.e-1, ymax=10.**3.0)
      #ax.set_xlim(p.limits[1][0], p.limits[1][1])
      ax.set_xlim(-0.6, 1.8)
   par_str = map(lambda x:'%.2f'%x, par)
   plt.suptitle('i-dropouts\n%s' % par_str, font_properties=Helvetica(20))

def plot_restmodel(par, mc, rlimits, binwidth=array([0.5,0.2])):
   # plot the model in (M_1500, logR) plane
   fig = plt.figure()
   ax = plt.subplot(111)
   pixdx = array([0.02,0.02])
   rmodel = bl.bivariate_RL_class()
   model0 = rmodel.bivariate_RL_M(par, rlimits, pixdx)
   #plt.hot()
   ax.imshow(model0.swapaxes(0,1), origin='lower')
   nbins = (rlimits[:,1]-rlimits[:,0]) / binwidth
   nbins = around(nbins).astype('int')
   npix = (rlimits[:,1]-rlimits[:,0]) / pixdx
   npix = around(npix).astype('int')
   xticks = arange(0, npix[0]+1, int(binwidth[0]/pixdx[0]))[::2]
   yticks = arange(0, npix[1]+1, int(binwidth[1]/pixdx[1]))[::5]
   ax.set_xticks(xticks)
   ax.set_yticks(yticks)
   #ax.grid(color='white', linestyle='-', linewidth=1.2)
   ax.set_xticklabels(rlimits[0,0]+xticks*pixdx[0], font_properties=Helvetica(14))
   ax.set_yticklabels(rlimits[1,0]+yticks*pixdx[1], font_properties=Helvetica(14))
   ax.set_xlabel(r'$M_{1500}$', font_properties=Helvetica(20))
   ax.set_ylabel(r'$\log R_e$ [pixel]', font_properties=Helvetica(20))
   #plt.colorbar(mappable=rmodel.model,ax=ax)
   divider =    divider = make_axes_locatable(ax)
   axLF = divider.append_axes("top", size="60%", pad=0.0)
   Marr = arange(0, (rlimits[0,1]-rlimits[0,0])/pixdx[0])
   axLF.semilogy(Marr, model0.sum(axis=1), 
                 color='black', nonposy='mask', lw=2.0)
   axLF.set_xticks([])
   axLF.set_yticks([])
   axLF.text(0.5, 0.2, 'Schechter Function', font_properties=Helvetica(18),
             horizontalalignment='center', verticalalignment='bottom',
             transform=axLF.transAxes)
   axSD = divider.append_axes("right", size="35%", pad=0.0)
   logRarr = arange(0, (rlimits[1,1]-rlimits[1,0])/pixdx[1])
   axSD.semilogx(model0.sum(axis=0), logRarr, color='black', lw=2, 
                 nonposx='mask')
   axSD.set_ylim(0, len(logRarr))
   axSD.set_yticklabels([])
   axSD.set_xticklabels([])
   axSD.text(0.05, 0.5, 'log-normal', font_properties=Helvetica(18), 
             horizontalalignment='left', verticalalignment='center',
             transform=axSD.transAxes)
   return ax, fig
   

def plot_no_tf_model(par, limits, drop, field, zdgrid):
   fig = plt.figure()
   ax = plt.subplot(111)
   pixdx = array([0.02,0.02])
   binwidth = array([0.5,0.2])
   if drop=='b':
      mc = bl.mconvert('M1500_to_i.txt')
      meankcorr=mc(4.0)
   elif drop == 'v':
      mc = bl.mconvert('M1500_to_z.txt')
      meankcorr=mc(5.0)
   model = bl.bivariate_lf(par, limits, pixdx, drop, field, kgrid=None, zdgrid=zdgrid,
      add_interloper=False, meankcorr=meankcorr, mc=mc)
   ax.imshow(model.model.swapaxes(0,1), origin='lower')
   nbins = (limits[:,1]-limits[:,0]) / binwidth
   nbins = around(nbins).astype('int')
   npix = (limits[:,1]-limits[:,0]) / pixdx
   npix = around(npix).astype('int')
   xticks = arange(0, npix[0]+1, int(binwidth[0]/pixdx[0]))
   yticks = arange(0, npix[1]+1, int(binwidth[1]/pixdx[1]))
   ax.set_xticks(xticks)
   ax.set_yticks(yticks)
   ax.grid(color='white', linestyle='-', linewidth=1.2)
   ax.set_xticklabels(limits[0,0]+xticks*pixdx[0])
   ax.set_yticklabels(limits[1,0]+yticks*pixdx[1])
   ax.set_xlabel(r'$i_{775}$-band magnitude')
   ax.set_ylabel(r'$\log R_e$ [pixel]')
   return ax, fig

def plot_full_model(par, limits):
   # quick & dirty full model: for B-dropouts in GOODS
   fig = plt.figure()
   ax = plt.subplot(111)
   pixdx = array([0.02,0.02])
   binwidth = array([0.5,0.2])
   zdgrid = zdist.read_zdgrid('zdgrid/zdgrid_bdrops_nolya.p')
   kgrid = mlutil.readkgrid('tfkernel/kernel_I.p')
   mc = bl.mconvert('M1500_to_i.txt')
   model = bl.bivariate_lf(par, limits, pixdx, 'b', 'goods', kgrid=kgrid, 
                           zdgrid=zdgrid, mc=mc, meankcorr=mc(4.0), 
                           add_interloper=False)
   ax.imshow(model.model.swapaxes(0,1), origin='lower', cmap=mpl.cm.hot)
   nbins = (limits[:,1]-limits[:,0]) / binwidth
   nbins = around(nbins).astype('int')
   npix = (limits[:,1]-limits[:,0]) / pixdx
   npix = around(npix).astype('int')
   xticks = arange(0, npix[0]+1, int(binwidth[0]/pixdx[0]))
   yticks = arange(0, npix[1]+1, int(binwidth[1]/pixdx[1]))
   ax.set_xticks(xticks)
   ax.set_yticks(yticks)
   ax.grid(color='white', linestyle='-', linewidth=1.2)
   ax.set_xticklabels(limits[0,0]+xticks*pixdx[0])
   ax.set_yticklabels(limits[1,0]+yticks*pixdx[1])
   ax.set_xlabel(r'$i_{775}$-band magnitude')
   ax.set_ylabel(r'$\log R_e$ [pixel]')
   ax.set_xlim(0,npix[0])
   ax.set_ylim(0,npix[1])
   return ax, fig
      
def plot_bdrops_wmodel(c, limits, ax):
   # plot the B-dropouts GALFIT mag & Re on the same model grid
   pixdx = array([0.02,0.02])
   npix = (limits[:,1]-limits[:,0]) / pixdx
   npix = around(npix).astype('int')
   cullcrit = in1d(c.cullflag, [0,1,3,12])
   galfitqf = ((c.reout_err/c.reout)<=0.6) & (c.chisqnu<=0.5)
   mags = c.magout[cullcrit&galfitqf]
   logre = log10(c.reout[cullcrit&galfitqf])
   x = (mags - limits[0,0]) / pixdx[0]
   y = (logre - limits[1,0]) / pixdx[1]
   #fig = plt.figure()
   #ax = fig.add_axes([0.1,0.1,0.78,0.65])
   ax.plot(x, y, '.', ms=10, c='green')
   #ax.set_xticklabels([])
   #ax.set_yticklabels([])
   ax.set_xlim(0,npix[0])
   ax.set_ylim(0,npix[1])
   return ax
   

def plot_3models(par, zdgridfile, kgrid):
   fig = plt.figure(figsize=(6,8))
   # set up axes grid
   grid = ImageGrid(fig, 111, nrows_ncols=(3,1), share_all=True, axes_pad=0.05, 
      aspect=1.5)
   raw_model = bl.bivariate_lf(par, bl.limits1, bl.pixdx, meankcorr=mc(4.0))
   zdmodel = bl.bivariate_lf(par, bl.limits1, bl.pixdx, zdgridfile=zdgridfile,
      mcfile='M1500_to_i.txt', meankcorr=mc(4.0))
   finalmodel = bl.bivariate_lf(par, bl.limits1, bl.pixdx, kgrid=kgrid,
      zdgridfile=zdgridfile, mcfile='M1500_to_i.txt', meankcorr=mc(4.0))
   vmax = max(zdmodel.model.ravel())
   # now show raw model
   plt.gray()
   grid[0].imshow(raw_model.model.swapaxes(0,1), origin='lower')
   grid[1].imshow(zdmodel.model.swapaxes(0,1), origin='lower', vmax=vmax)
   grid[2].imshow(finalmodel.model.swapaxes(0,1), origin='lower', vmax=vmax)
   for i in range(3):
      grid[i].set_xticks(arange(50, 300, 50))
      grid[i].set_xticklabels(arange(22., 27., 1.),size=10)
      grid[i].tick_params(axis='x', color='white')
      grid[i].set_xlim(50, 275)
      grid[i].set_yticks(arange(0, 200, 50))
      grid[i].set_yticklabels(0.03 * 10.**arange(-1.,3.,1.),size=10)
      grid[i].tick_params(axis='y', color='white')
      grid[i].axis["left"].label.set_text('Re [arcsec]')
      grid[i].set_ylim(0, 130)
   grid[2].axis["bottom"].label.set_text(r"$m$")
   return fig, grid

def plot_2models_thesis(par, parfile, drop='f435w', gfband='f775w',
                        field='goods', cmap=mpl.cm.hot):
   #fig = plt.figure(figsize=(8,9))
   #ax1 = fig.add_subplot(2,1,1)
   #ax2 = fig.add_subplot(2,1,2)

   pmb = plot_bivariate_RL(parfile, field, drop, gfband, z_mean=4.0)
   #pmb.plot_rawmodel(par, ax=ax1)
   #pmb.plot_zdgridmodels(par, ax=ax2)
   pmb.plot_step1(par, figsize=(8,9))
   plt.subplots_adjust(left=0.13,right=0.92)
   pmb.plot_step2(par, figsize=(8,9))
   pmb.plot_step3(par, figsize=(8,9), mag0=21.0, mag1=26.5)

def plot_2models(par, zdgrid, kgrid, limits, drop='f435w', 
                 field='goods', cmap=mpl.cm.hot, normalize=True):
   """
   The production script that makes Figure 7 in Huang et al. 2013.
   """
   fig = plt.figure(figsize=(8,10))
   # set up axes grid
   #grid = ImageGrid(fig, 111, nrows_ncols=(2,1), share_all=True, axes_pad=0.05, 
   #   aspect=1.5)
   raw_limits = limits.copy()
   if drop == 'f435w':
      mcfile = 'kcorr/M1500_to_f775w.txt'
      z0 = 4.0
   elif drop == 'f606w':
      mcfile = 'kcorr/M1500_to_f850lp.txt'
      z0 = 5.0
   mc = bl.mconvert(mcfile)
   raw_limits[0] = limits[0] - mc(z0)
   raw_limits[0] = around(raw_limits[0], 1)
   raw_model = bl.bivariate_lf_M(par, raw_limits, bl.pixdx)
   finalmodel = bl.bivariate_RL(par, limits, bl.pixdx, drop, field, 
                                kgrid=kgrid, zdgrid=zdgrid, mc=mc,
                                add_interloper=False)
   # normalize the models
   if normalize:
   	raw_model.model = raw_model.model / sum(raw_model.model.ravel())
   	finalmodel.model = finalmodel.model / sum(finalmodel.model.ravel())
   shape1 = shape(raw_model.model)
   shape2 = shape(finalmodel.model)
   nxr, nyr = shape(raw_model.model)
   nx, ny = shape(finalmodel.model)
   dm = limits[0][1] - limits[0][0]
   dlogr = limits[1][1] - limits[1][0]
   nxticks = int(dm / 0.5)
   nyticks = int(dlogr / 0.2)
   yticklabels = arange(limits[1][0]+0.2, limits[1][1], 0.4)
   yticklabels = around(yticklabels, 1)
   yticks = (yticklabels - limits[1][0]) / bl.pixdx[0]
   yticks = yticks.astype('int')
   #for i in range(nyticks):
   #   if i % 2 == 0:
   #      #yticks += ['%.3f' % (0.03 * 10.**(limits[1][0]+i*0.2))]
   #      yticks += ['%.1f' % (limits[1][0]+i*0.2)]
   #   else:
   #      yticks += ['']
   #vmax = max(raw_model.model.ravel())
   # now show raw model
   plt.hot()
   ax1 = fig.add_subplot(211)
   ax2 = fig.add_subplot(212)
   axes = [ax1, ax2]
   img1 = ax1.imshow((raw_model.model.swapaxes(0,1)), origin='lower', 
                     cmap=cmap, aspect=0.8)
   cax1 = ax1.figure.colorbar(img1,ax=ax1)
   cax_ticks1 = linspace(0., max(raw_model.model.ravel()), 5)
   cax1.set_ticks(cax_ticks1)
   cax1.set_ticklabels(map(lambda x:'%.1e'%x, cax_ticks1))
   img2 = ax2.imshow((finalmodel.model.swapaxes(0,1)), origin='lower', 
                     cmap=cmap, aspect=0.8)
   cax2 = ax2.figure.colorbar(img2,ax=ax2)
   cax_ticks2 = linspace(0., max(finalmodel.model.ravel()), 5)
   cax2.set_ticks(cax_ticks2)
   cax2.set_ticklabels(map(lambda x:'%.1e'%x, cax_ticks2))
   #grid[0].imshow(raw_model.model.swapaxes(0,1), origin='lower')
   #grid[1].imshow(finalmodel.model.swapaxes(0,1), origin='lower')
   for i in range(2):
      axes[i].tick_params(axis='x', color='white')
      axes[i].set_xlim(0, nx)
      #grid[i].yaxis.set_major_formatter(rfmt)
      axes[i].set_yticks(yticks)
      axes[i].set_yticklabels(yticklabels, font_properties=Helvetica(14))
      axes[i].tick_params(axis='y', color='white')
      axes[i].set_ylim(0, ny)
      axes[i].set_ylabel(r"$\log_{10}(R_e)$ [pixel]",
                         font_properties=Helvetica(20))
      #axes[i].grid(b=True, which='major', linestyle='-', color='white')
   axes[0].set_xticks(arange(0, nxr+nxr/nxticks, nxr/nxticks))
   axes[0].set_xticklabels(arange(raw_limits[0][0], raw_limits[0][1]+0.5, 
                           0.5), font_properties=Helvetica(14))
   axes[0].set_xlabel(r"$M_{1500}$", font_properties=Helvetica(20))
   axes[0].text(0.1,0.1,"Original",color="white",transform=axes[0].transAxes,
                font_properties=Helvetica(18))
   axes[1].set_xticks(arange(0, nx+nx/nxticks, nx/nxticks))
   axes[1].set_xticklabels(arange(limits[0][0], limits[0][1]+0.5, 0.5),
                           font_properties=Helvetica(14))
   axes[1].set_xlabel(r"$m_{\mathrm{F775W}}$",
                      font_properties=Helvetica(20))
   axes[1].text(0.1,0.1,"Transformed",color="white",
                transform=axes[1].transAxes, font_properties=Helvetica(18))
   #plt.colorbar()
   return fig, axes

