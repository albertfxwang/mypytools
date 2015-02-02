## Tools for processing/displaying EAZY results
import numpy as np
import os, subprocess
from pygoods import Ftable, sextractor
from sedtools.bands import filters
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz
from PhotomTools import photomUtils as pu
import pysynphot as S
import cosmoclass
from stats import distributions

eazy_path = '/Users/khuang/bin/eazy'
# a symbol for magnitude upper limits
downarrow = [(-2,0),(2,0),(0,0),(0,-4),(-2,-2),(0,-4),(2,-2),(0,-4),(0,0)]
pc_cm = 3.0857e18   # parsec in cm
area_10pc = 4 * np.pi * (10*pc_cm)**2  # in cm^2
L_solar = 3.826e33  # in erg/s
bc03_phys = '/Users/khuang/eazy-1.00/templates/BC03/bc03_m42_tau.phys'
bc03_phys_dustfree = '/Users/khuang/eazy-1.00/templates/BC03/bc03_m42_tau_dustfree.phys'
cosmo = cosmoclass.cosmoclass(70., 0.3, 0.7)  # default cosmology parameters

def scinote2exp(scinote, nprec=3):
   # convert scientific notation in the format of 1.0e10 into (1.0, 10)
   n = scinote.split('e')
   fltstr = "%%.%df" % nprec
   if np.abs(int(n[1])) <= 2:
      return fltstr % float(scinote)
   else:
      # return float(n[0]), int(n[1])
      return "%.2f \\times 10^{%d}" % (float(n[0]), int(n[1]))

class EAZY(object):
   def __init__(self, tempfluxunit='fnu', physfile=bc03_phys, paramfile='zphot.param', with_specz=False):
      # catalog is the file containing all the fluxes; first column in catalog
      # is the object ID.
      with open(paramfile, 'rb') as f1:
         lines = f1.readlines()
         for l in lines:
            if l.startswith('CATALOG_FILE'):
               catalog = l.split()[1]
      self.input_catalog = catalog
      self.c = sextractor(catalog)
      if with_specz:
         assert len(self.c._d) % 2 == 0, "Number of columns doesn't seem right... missing some fluxes or flux errors?\n(or do you mean with_specz=True?)"
      else:
         assert len(self.c._d) % 2 == 1, "Number of columns doesn't seem right... missing some fluxes or flux errors?\n(or do you mean with_specz=True?)"
      self.nbands = (len(self.c._d) - 1) / 2
      self.objid = self.c._1
      # also assumes that the first line of catalog is the header line
      f = open(catalog)
      hdr = f.readline()
      f.close()
      self.filter_names = []
      # read filter names; they should be within the keys in 
      # sedtools.bands.filters
      if os.path.exists('OUTPUT/photz.zout'):
         self.read_output()
      for b in hdr.split():
         if b.lower().startswith('f_'):
            self.filter_names += [b.lower()[2:]]
      if tempfluxunit == 'Lsolar':
         self.normfactor = area_10pc / L_solar
      else:
         self.normfactor = 1.0
      self.physfile = physfile
      self.phys = sextractor(physfile)

   def read_output(self):
      # Read in the output photo-z file
      out = sextractor("OUTPUT/photz.zout")
      hdr = open("OUTPUT/photz.zout", "rb").readline()
      hdr_list = hdr.split()[1:]
      # print hdr_list
      for i in range(1, len(hdr_list)+1):
         setattr(self, hdr_list[i-1].lower(), getattr(out, "_%d" % i))

   def read_photz_errors(self, mode='a', nsig=1):
      """
      Read phot-z errors (relative to the peak of P(z)) for plotting with 
      plt.errorbar.
      """
      self.read_output()
      if mode == 'a':
         zbest = self.z_a
      elif mode == '1':
         zbest = self.z_1
      else:
         zbest = self.z_2
      if nsig == 1:
         zlo = self.l68
         zhi = self.u68
      else:
         zlo = self.l95
         zhi = self.u95
      error_lo = zbest - zlo
      error_hi = zhi - zbest
      return [zbest, error_lo, error_hi]

   def read_SED(self, objid, mode='1'):
      """
      Read the best-fit SED and determine the best-fit physical properties.
      """
      temp = sextractor("OUTPUT/%s.temp_sed" % str(objid))
      # temp flux is in uJy, the same as catalog flux
      # Also read the best-fit template number
      header = temp._header.split('\n')
      h2 = header[1].split('templ:norm')[-1].split(':')
      MOD_BEST = int(h2[0])
      ## MOD_BEST is the model number in TEMPLATES_FILE; the file order has
      ## to match between TEMPLATES_FILE and physfile for the best-fit physical
      ## Calculate the normalization factor, from template to observed flux
      ## The goal is to match their L_nu in rest-frame
      sp_obs = S.ArraySpectrum(temp._1, temp._2, fluxunits='fnu')  # in uJy
      sp_obs = sp_obs * 1.e-29  # convert from uJy to erg/s/cm^2/Hz
      sp_mod = S.FileSpectrum(os.path.split(self.physfile)[0] + '/' + self.phys.filename[MOD_BEST-1])  # in flam
      sp_mod.convert('fnu')  # convert flux unit to fnu
      if mode == '1':
         zpeak = self.z_1[self.id==objid][0]
      elif mode == 'a':
         zpeak = self.z_a[self.id==objid][0]
      else:
         zpeak = self.z_2[self.id==objid][0]
      LUMDIST = cosmo.lumdist(zpeak, unit='Mpc') * 1e6  # in parsec
      FACTOR = (1./(1.+zpeak)) * (LUMDIST / 10.)**2  
      # a factor in front of template flux ratios
      # Use the average at rest-frame 2000 A and 6000 A for the normalization factor
      A1 = FACTOR * (sp_obs.sample(2000*(1+zpeak)) / sp_mod.sample(2000))
      A2 = FACTOR * (sp_obs.sample(10000*(1+zpeak)) / sp_mod.sample(10000))
      # print "A1, A2: %.2e, %.2e" % (A1, A2)
      NORM = (A1 + A2) / 2.
      ## parameters to be correct.
      mass_best = self.phys.smass[MOD_BEST-1] * NORM
      log_age_best = self.phys.log_age[MOD_BEST-1] 
      sfr_best = self.phys.sfr[MOD_BEST-1] * NORM
      tau = self.phys.tau[MOD_BEST-1]
      ebmv = self.phys.ebmv[MOD_BEST-1]
      return temp, zpeak, mass_best, log_age_best, sfr_best, tau, ebmv, MOD_BEST

class PlotEAZY(EAZY):
   def __init__(self, tempfluxunit='fnu', physfile=bc03_phys, wavelim=[3000.,8e4], with_specz=False):
      EAZY.__init__(self, tempfluxunit=tempfluxunit, physfile=physfile, with_specz=with_specz)
      self.lambda_max = 30000.  # default value for maximum lambda in plots
      self.lambda_min = 3000.
      self.flux_ujy_max = 1.0
      self.physfile = physfile
      self.phys = sextractor(physfile)
      self.wavelim = wavelim

   def plot_Pz(self, objid, ax=None, savefig=True, mode='a', specz=-1, **plot_kwargs):
      assert objid in self.objid, "Object ID %s does not exist in catalog." % str(objid)
      pz = sextractor("OUTPUT/%s.pz" % str(objid))
      if ax == None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
      pz._2 = pz._2 - pz._2.min()
      print pz._2.min(), pz._2.max()
      prob = np.exp(-0.5 * pz._2)
      print prob.min(), prob.max()
      probIntegrated = cumtrapz(prob, x=pz._1)[-1]
      prob = prob / probIntegrated
      ax.plot(pz._1, prob, **plot_kwargs)
      if specz > 0:
         ax.plot([specz, specz], [0., ax.get_ylim()[1]], lw=2, ls='--', 
                 c='navy', label=r'$z_{\mathrm{spec}}=%.3f$' % specz)
      ax.legend(loc=0)
      ax.set_xlabel('Redshift')
      ax.set_ylabel(r'$P(z)$')
      if mode == 'a':
         z_peak = self.z_a[self.id==objid][0]
      elif mode=='1':
         z_peak = self.z_1[self.id==objid][0]
      else:
         z_peak = self.z_2[self.id==objid][0]
      ax.set_title('Object %s [z_peak = %.3f]' % (str(objid), z_peak))
      if savefig:
         fig.savefig("OUTPUT/%s_Pz.png" % str(objid))
      return ax, z_peak

   def plot_multiple_Pz(self, objids, ax=None, mode='a', **plot_kwargs):
      """
      Show the P(z) for multiple objects in the same plot.
      """
      if ax == None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
      for x in objids:
         ax, z_peak = self.plot_Pz(x, ax=ax, savefig=False, mode=mode, **plot_kwargs)
         ax.lines[-1].set_label('%s [z_peak = %.2f]' % (x, z_peak))
      ax.set_title("")
      ax.legend()
      return ax

   def plot_photom(self, objid, ax=None, savefig=False, ebar_kwargs={'ms':10, 'mec':'black', 'mfc':'none', 'ecolor':'black', 'mew':1.5, 'capsize':8, 'capthick':1.5}):
      # Plot the observed fluxes and the best-fit SED
      if ax == None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
         ax.set_yscale('log')
         ax.set_xscale('log')
      obs = sextractor("OUTPUT/%s.obs_sed" % str(objid))
      # First, plot the detections
      detect = (obs._2 > 0)
      ax.errorbar(obs._1[detect], obs._2[detect], yerr=obs._3[detect], 
                  fmt='s', **ebar_kwargs)
      ymin = ax.get_ylim()[0]
      # Then plot upper limits
      # yerr_bottom = obs._2[detect==False] - 1.e-10
      # ax.errorbar(obs._1[detect==False], obs._2[detect==False], 
      #             yerr=[yerr_bottom,obs._3[detect==False]], uplims=True, 
      #             **ebar_kwargs)
      ax.scatter(obs._1[detect==False], obs._3[detect==False],
                 marker=downarrow, edgecolor=ebar_kwargs['ecolor'],
                 s=(ebar_kwargs['ms']*2)**2, linewidths=1.2)
      # ax.set_ylim(ymin=ymin)
      ylims = ax.get_ylim()
      mag_lo = pu.uJy2ABmag(ylims[0])
      mag_hi = pu.uJy2ABmag(ylims[1])
      yticks_mag = np.arange(np.ceil(mag_lo), mag_hi-2.0, -1.0)
      yticks_flux = pu.ABmag2uJy(yticks_mag)
      ax.set_yticks(yticks_flux)
      ax.set_yticklabels(["%d" % int(x) for x in yticks_mag])
      ax.set_ylabel('Observed magnitude')
      self.lambda_max = obs._1.max() * 1.2
      self.lambda_min = obs._1.min() * 0.8
      self.flux_ujy_max = obs._2.max() * 10.
      ax.set_xlim(self.lambda_min, self.lambda_max)
      # xlim will be overridden if also plots SED
      ax.set_ylim(ymax=self.flux_ujy_max)
      return ax

   def plot_SED(self, objid, xmax=None, ax=None, savefig=False, plot_kwargs={'color':'black'}, mode='1'):
      if ax == None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
         ax.set_yscale('log')
         ax.set_xscale('log')
      (temp,z_best,mass_best,log_age_best,sfr_best,tau,ebmv,mod_best) = self.read_SED(objid, mode=mode)
      print "Best-fit parameters:"
      print "Stellar mass = %.2e [M_solar]" % mass_best
      print "Age = %.2e [yrs]" % 10.**log_age_best
      print "SFR = %.2f [M_solar/yr]" % sfr_best
      print "tau = %.2f [Gyr]" % tau
      print "E(B-V) = %.2f" % ebmv
      # Also try to figure out E(B-V) and tau...
      l = open("OUTPUT/%s.temp_sed" % str(objid)).readlines()[1]
      l = l.split()
      for i in range(len(l)):
         if l[i] == 'z=':
            z_temp = float(l[i+1])
            break
      ax.plot(temp._1, temp._2, **plot_kwargs)
      # ax.set_xlim(self.lambda_min, self.lambda_max)
      ax.set_xlim(*self.wavelim)
      # lambda_ticks = [5000.,10000.,20000.,30000.,40000.,50000.]
      lambda_ticks = [5.e3, 1.e4, 2.e4, 5.e4]
      ax.set_xticks(lambda_ticks)
      ax.set_xticklabels(map(lambda x: "%.1f" % (x/1.e4), lambda_ticks))
      ax.set_xlabel(r"$\lambda$ [$\mu$m]")
      if not len(ax.get_ylabel()):
         ax.set_ylabel(r"$F_{\nu}$ [$\mu$Jy]")
      if mode == 'a':
         z_peak = self.z_a[self.id==objid][0]
      elif mode == '1':
         z_peak = self.z_1[self.id==objid][0]
      else:
         z_peak = self.z_2[self.id==objid][0]
      ax.set_title('Object %s [z_peak = %.3f]' % (str(objid), z_peak), weight='bold')
      if savefig:
         plt.savefig('OUTPUT/%s_SED.png' % str(objid))
      return ax, z_peak, [mass_best, 10.**log_age_best, sfr_best, tau, ebmv]


   def plot_all(self, objid, objname="", savefig=True, legend_loc=1, mode='1', specz=-1):
      fig = plt.figure(figsize=(10,12))
      ax1 = fig.add_subplot(211)
      ax1.set_yscale('log')
      ax1.set_xscale('log')
      ax1 = self.plot_photom(objid, ax=ax1)
      ax1, zp, props = self.plot_SED(objid, ax=ax1, mode=mode,
                                     plot_kwargs=dict(color='blue',lw=2))
      ax2 = fig.add_subplot(212)
      ax2, zp = self.plot_Pz(objid, ax=ax2, savefig=False, mode=mode, lw=2,
                             specz=specz)
      ax2.set_title("")
      if len(objname):
         ax1.set_title('Object %s [z_peak = %.3f]' % (objname, zp))
         ax2.set_title('Object %s [z_peak = %.3f]' % (objname, zp))
      plt.draw()
      if savefig:
         plt.savefig("%s_SED_Pz.png" % str(objid))
      return fig

   def plot_all_1panel(self, objid, objname="", savefig=False, legend_loc=1, mode='1', specz=-1, xlow2=0.65, xsize2=0.3, ylow2=0.12, ysize2=0.3, legendfontsize='medium', pztickfontsize='large', pztextfontsize='large', zmax=6):
      fig = plt.figure(figsize=(10, 8))
      ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
      print "\nPlotting %s..." % objid
      ax1.set_yscale('log')
      ax1.set_xscale('log')
      ax1 = self.plot_photom(objid, ax=ax1)
      ax1, zp, props = self.plot_SED(objid, ax=ax1, mode=mode,
                                     plot_kwargs=dict(color='blue',lw=2))
      # Now plot P(z)
      bbox1 = ax1.get_position()
      xlow2, ylow2 = bbox1.p0 + bbox1.size * np.array([xlow2, ylow2])
      xsize2, ysize2 = bbox1.size * np.array([xsize2, ysize2])
      # print "Boundaries for ax1:"
      # print np.concatenate([bbox1.p0, bbox1.size])
      # print "New values for the boundaries of ax2:"
      print xlow2, ylow2, xsize2, ysize2
      ax2 = fig.add_axes([xlow2, ylow2, xsize2, ysize2])
      ax2, zp = self.plot_Pz(objid, ax=ax2, savefig=False, mode=mode, lw=2,
                             specz=specz)
      ax2.set_title("")
      ax2.set_xlim(-0.1, zmax)
      if ax2.legend():
         ax2.legend().set_visible(False)  # turn off P(z) legend
      ax2.set_yticklabels([])
      ax2.set_ylabel(ax2.get_ylabel(), size=pztickfontsize)
      ax2.legend(fontsize=pztextfontsize)
      if len(objname):
         ax1.set_title('Object %s [z_peak = %.3f]' % (objname, zp))
      plt.draw()
      if savefig:
         plt.savefig("%s_SED_Pz.png" % str(objid))
      return fig


   def plot_multiple_all(self, objids, objnames=[], ncols=2, savefig=False, legend_loc=1, mode='a', **plt_kwargs):
      """
      Plot SED fits and P(z) in subplots. One subplot for each SED fit, and 
      one subplot showing all the P(z) curves.
      """
      plt_kwargs_def = dict(color='blue',lw=2)
      for k in plt_kwargs.keys():
         plt_kwargs_def[k] = plt_kwargs[k]
      fig = plt.figure(figsize=(14,12))
      n_plots = len(objids) + 1
      if n_plots % ncols == 0:
         nrows = n_plots / ncols
      else:
         nrows = n_plots / ncols + 1
      print "nrows, ncols:", nrows, ncols
      axes = []
      # Plot individual SED fits
      for i in range(len(objids)):
         ax = fig.add_subplot(nrows, ncols, i+1)
         ax.set_yscale('log')
         ax.set_xscale('log')
         ax = self.plot_photom(objids[i], ax=ax)
         ax, zp, props = self.plot_SED(objids[i], ax=ax, mode=mode,
                                       plot_kwargs=plt_kwargs_def)
         ax.text(0.1, 0.9, objids[i], transform=ax.transAxes, ha='left', 
                 va='top', bbox=dict(boxstyle='round',facecolor='none'),
                 size='x-large')
         ax.set_title("")
         axes.append(ax)
      # Now plot all P(z) in the same panel
      ax = fig.add_subplot(nrows, ncols, i+2)
      ax = self.plot_multiple_Pz(objids, ax=ax, mode=mode, lw=2)
      ax.legend(loc=legend_loc, fontsize='large')
      axes.append(ax)
      return axes

   def plot_photz_vs_specz(self, specz, zmax, ax=None, mode='a', nsig=1, **ebar_kwargs):
      """
      Plot phot-z v.s. spec-z comparison.
      """
      ebar_kwargs_def = dict(fmt='s', elinewidth=1.5, capsize=8)
      for k in ebar_kwargs.keys():
         ebar_kwargs_def[k] = ebar_kwargs[k]
      zbest, error_lo, error_hi = self.read_photz_errors(mode=mode, nsig=nsig)
      assert len(specz) == len(zbest), "Please provide the same number of spec-z's as the number of objects in the catalog."
      if ax == None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
      ax.errorbar(specz, zbest, yerr=[error_lo, error_hi], **ebar_kwargs_def)
      ax.plot([0, zmax], [0, zmax], lw=2, ls='--', c='black')
      ax.set_xlim(0, zmax)
      ax.set_ylim(0, zmax)
      ax.set_xlabel('spec. redshift')
      ax.set_ylabel('phot. redshift')
      return ax


class PlotEAZY_wspecz(PlotEAZY):
   def __init__(self, tempfluxunit='fnu', physfile=bc03_phys, wavelim=[3000.,8e4]):
      PlotEAZY.__init__(self, tempfluxunit=tempfluxunit, physfile=physfile, 
                        wavelim=wavelim, with_specz=True)
   def plot_all(self, objid, xmax=None, ax=None, savefig=False, textbox=True, plot_kwargs={'color':'black'}, mode='1'):
      fig = plt.figure(figsize=(10,8))
      ax1 = fig.add_subplot(111)
      ax1.set_yscale('log')
      ax1.set_xscale('log')
      ax1 = self.plot_photom(objid, ax=ax1)
      ax1, zp, props = self.plot_SED(objid, ax=ax1, mode=mode, plot_kwargs=plot_kwargs)
      if mode == '1':
         chi2 = self.chi_1
      elif mode == 'a':
         chi2 = self.chi_a
      else:
         chi2 = self.chi_2
      ax1.set_title('Object %s [z_spec = %.2f]' % (objid, zp))
      if textbox:
         mass, age, sfr, tau, ebmv = props
         masstxt = scinote2exp('%e' % mass)
         sedtext = "With IRAC:\n"
         sedtext = sedtext + "$M_{\mathrm{star}} = %s/\\mu\ \mathrm{[M_{\odot}]}$\n" % masstxt
         sedtext = sedtext + "$E(B-V) = %.2f$\n" % ebmv
         sedtext = sedtext + "$\mathrm{SFR} = %.2f/\\mu\ \mathrm{[M_{\odot}/yr]}$\n" % sfr
         ssfr = sfr / mass * 1.e9  # in units of Gyr^-1
         sedtext = sedtext + "$\mathrm{sSFR} = %.2f\ [\mathrm{Gyr}^{-1}]$\n" % (ssfr)
         age = scinote2exp('%e' % age)
         sedtext = sedtext + "$\mathrm{Age} = %s\ \mathrm{[yrs]}$\n" % age
         sedtext = sedtext + "$\\tau = %.1f\ \mathrm{[Gyrs]}$" % tau
         chi2_nu = chi2[self.id==objid][0] / self.nfilt[self.id==objid][0]
         sedtext = sedtext + "\n$\chi^2 = %.2f$" % (chi2_nu)
         ax1.text(0.95, 0.05, sedtext, size='large', ha='right', va='bottom',
                 transform=ax1.transAxes, multialignment='left',
                 bbox=dict(boxstyle='round,pad=0.3',facecolor='LightPink'))
      plt.draw()
      if savefig:
         plt.savefig("%s_SED_Pz.png" % str(objid))
      return ax1

   def plot_all_multifits(self, objids, objname="", savefig=False, plot_kwargs={'lw':2}, mode='1'):
      """
      Plot multiple fits of the same object. Each fit is represented by a 
      different objid in the EAZY catalog, but they have the same photometry.
      This method doesn't check if the photometry is the same---it just uses
      the photometry from the first objid. The user should make sure that
      all objids are from the same photometry.
      """
      fig = plt.figure(figsize=(10,8))
      ax = fig.add_subplot(111)
      ax.set_yscale('log')
      ax.set_xscale('log')
      if mode == '1':
         chi2 = self.chi_1
      elif mode == 'a':
         chi2 = self.chi_a
      else:
         chi2 = self.chi_2
      ax = self.plot_photom(objids[0], ax=ax)
      for x in objids:
         ax, zp, props = self.plot_SED(x, ax=ax, mode=mode,
                                       plot_kwargs=plot_kwargs)
         chi2_nu = chi2[self.id==x][0] / self.nfilt[self.id==x][0]
         ax.lines[-1].set_label(r"$z=%.2f$; $\chi^2_{\nu}=%.2f$" % (zp, chi2_nu))
         print ax.lines[-1].get_label()
         ax.set_title("")
      ax.legend(loc=2, fontsize='x-large')
      if objname:
         ax.set_title(objname, size='xx-large')
      plt.draw()
      return ax


def plot_HST_IRAC_SED(objid, ax=None, colors=['blue', 'red'], savefig=True, legend_loc=2, mode='a'):
   curdir = os.getcwd()
   print "Plotting object %s..." % str(objid)
   os.chdir('hst_only')
   # Read catalog name
   with open('zphot.param', 'rb') as f1:
      lines = f1.readlines()
      for l in lines:
         if l.startswith('CATALOG_FILE'):
            catalog_hst = l.split()[1]
   p1 = PlotEAZY(catalog_hst)
   if mode == 'a':
      zp_1 = p1.z_a[p1.id==objid][0]
   elif mode == '1':
      zp_1 = p1.z_1[p1.id==objid][0]
   else:
      zp_1 = p1.z_2[p1.id==objid][0]
   os.chdir(curdir+'/with_irac')
   with open('zphot.param', 'rb') as f2:
      lines = f2.readlines()
      for l in lines:
         if l.startswith('CATALOG_FILE'):
            catalog_irac = l.split()[1]
   p2 = PlotEAZY(catalog_irac)
   if mode == 'a':
      zp_2 = p2.z_a[p2.id==objid][0]
   elif mode == '1':
      zp_2 = p2.z_1[p2.id==objid][0]
   else:
      zp_2 = p2.z_2[p2.id==objid][0]
   if ax == None:
      fig = plt.figure()
      ax = fig.add_axes([0.12,0.12,0.78,0.78])
   ax.set_yscale('log')
   ax.set_xscale('log')
   ax = p2.plot_photom(objid, ax=ax, savefig=False)
   os.chdir('../hst_only')
   print "Plotting HST_ONLY:"
   ax, z_peak, props_hst = p1.plot_SED(objid, ax=ax, savefig=False, mode=mode,
                    plot_kwargs={'color':colors[0],'label':'HST only' % zp_1,'lw':1.0, 'ls':'--'}) 
   os.chdir(curdir+'/with_irac')
   print "Plotting WITH_IRAC:"
   ax, z_peak, props_irac = p2.plot_SED(objid, ax=ax, savefig=False, mode=mode,
                    plot_kwargs={'color':colors[1],'label':'With IRAC' % zp_2,'lw':2.0})
   ax.set_xlim(xmin=p2.lambda_min * 0.8, xmax=p2.lambda_max * 1.5)
   ax.legend(loc=legend_loc) 
   os.chdir(curdir)
   if savefig:
      plt.savefig('%s_SED.png' % str(objid))
   return ax, z_peak, props_hst, props_irac

def plot_HST_IRAC_Pz(objid, ax=None, colors=['blue', 'maroon'], savefig=True, legend_loc=1, mode='a'):
   curdir = os.getcwd()
   print "Plotting P(z) for object %s..." % str(objid)
   os.chdir('hst_only')
   # Read catalog name
   with open('zphot.param', 'rb') as f1:
      lines = f1.readlines()
      for l in lines:
         if l.startswith('CATALOG_FILE'):
            catalog_hst = l.split()[1]
   p1 = PlotEAZY(catalog_hst)
   if mode == 'a':
      zp_1 = p1.z_a[p1.id==objid][0]
   elif mode == '1':
      zp_1 = p1.z_1[p1.id==objid][0]
   os.chdir(curdir+'/with_irac')
   with open('zphot.param', 'rb') as f2:
      lines = f2.readlines()
      for l in lines:
         if l.startswith('CATALOG_FILE'):
            catalog_irac = l.split()[1]
   p2 = PlotEAZY(catalog_irac)
   if mode == 'a':
      zp_2 = p2.z_a[p2.id==objid][0]
   elif mode == '1':
      zp_2 = p2.z_1[p2.id==objid][0]
   if ax == None:
      fig = plt.figure()
      ax = fig.add_axes([0.12,0.12,0.78,0.78])
   os.chdir('../hst_only')
   ax, z_peak = p1.plot_Pz(objid, ax=ax, savefig=False, color=colors[0],
                   label='HST only (z_peak=%.3f)' % zp_1, ls='--', lw=1.5,
                   mode=mode) 
   os.chdir(curdir+'/with_irac')
   ax, z_peak = p2.plot_Pz(objid, ax=ax, savefig=False, color=colors[1],
                   label='With IRAC (z_peak=%.3f)' % zp_2, lw=3.0,
                   mode=mode)
   ax.legend(loc=legend_loc) 
   os.chdir(curdir)
   if savefig:
      plt.savefig('%s_Pz.png' % str(objid))
   return ax, z_peak

def plot_HST_IRAC_all(objid, colors=['blue','maroon'], outputfile="", legend_loc=2, mode='1', textbox=True):
   # Plot both P(z) and SED for the object with ID objid
   fig = plt.figure(figsize=(10,12))
   ax1 = fig.add_subplot(2, 1, 1)
   ax1, z_peak, props_hst, props_irac = plot_HST_IRAC_SED(objid, ax=ax1, 
         colors=colors, savefig=False, legend_loc=legend_loc, mode=mode)
   title = 'Object %s [z_peak = %.3f]' % (str(objid), z_peak)
   ax1.set_title(title)
   ax2 = fig.add_subplot(2, 1, 2)
   ax2, z_peak = plot_HST_IRAC_Pz(objid, ax=ax2, colors=colors, savefig=False,
                                 legend_loc=legend_loc, mode=mode)
   if textbox:
      mass, age, sfr, tau, ebmv = props_irac
      masstxt = scinote2exp('%e' % mass)
      sedtext = "With IRAC:\n"
      sedtext = sedtext + "$M_{\mathrm{star}} = %s/\\mu\ \mathrm{[M_{\odot}]}$\n" % masstxt
      sedtext = sedtext + "$E(B-V) = %.2f$\n" % ebmv
      sedtext = sedtext + "$\mathrm{SFR} = %.2f/\\mu\ \mathrm{[M_{\odot}/yr]}$\n" % sfr
      ssfr = sfr / mass * 1.e9  # in units of Gyr^-1
      sedtext = sedtext + "$\mathrm{sSFR} = %.2f\ [\mathrm{Gyr}^{-1}]$\n" % (ssfr)
      age = scinote2exp('%e' % age)
      sedtext = sedtext + "$\mathrm{Age} = %s\ \mathrm{[yrs]}$\n" % age
      sedtext = sedtext + "$\\tau = %.1f\ \mathrm{[Gyrs]}$" % tau
      ax1.text(0.95, 0.05, sedtext, size='large', ha='right', va='bottom',
              transform=ax1.transAxes, multialignment='left',
              bbox=dict(boxstyle='round,pad=0.3',facecolor='LightPink'))
   # if objname:
      # ax2.set_title(title)
   ax2.set_title("")
   if len(outputfile):
      plt.savefig(outputfile)
   plt.show()
   return [ax1, ax2]

def plot_dustfree(objid, outputfile="", mode="1", textbox=True):
   """
   Plot EAZY fitting results for dust-free templates only. Only consider the 
   WITH_IRAC case.
   Should be called inside the dust-free EAZY directory so that the OUTPUT 
   directory exists.
   """
   P = PlotEAZY(physfile=bc03_phys_dustfree)
   fig = P.plot_all(objid, objname=objid, savefig=False, mode=mode)
   # Now add text box showing the best-fit properties
   ax1, ax2 = fig.axes
   props_df = P.read_SED(objid)
   if textbox:
      mass = props_df[2]
      age = 10.**props_df[3]
      sfr = props_df[4]
      tau = props_df[5]
      ebmv = props_df[6]
      masstxt = scinote2exp('%e' % mass)
      sedtext = "With IRAC:\n"
      sedtext = sedtext + "$M_{\mathrm{star}} = %s/\\mu\ \mathrm{[M_{\odot}]}$\n" % masstxt
      sedtext = sedtext + "$E(B-V) = %.2f$\n" % ebmv
      sedtext = sedtext + "$\mathrm{SFR} = %.2f/\\mu\ \mathrm{[M_{\odot}/yr]}$\n" % sfr
      ssfr = sfr / mass * 1.e9  # in units of Gyr^-1
      sedtext = sedtext + "$\mathrm{sSFR} = %.2f\ [\mathrm{Gyr}^{-1}]$\n" % (ssfr)
      age = scinote2exp('%e' % age)
      sedtext = sedtext + "$\mathrm{Age} = %s\ \mathrm{[yrs]}$\n" % age
      sedtext = sedtext + "$\\tau = %.1f\ \mathrm{[Gyrs]}$" % tau
      ax1.text(0.95, 0.05, sedtext, size='large', ha='right', va='bottom',
              transform=ax1.transAxes, multialignment='left',
              bbox=dict(boxstyle='round,pad=0.3',facecolor='LightPink'))
   ax2.set_title("")
   plt.draw()
   if len(outputfile):
      fig.savefig(outputfile)
   return fig

def plot_HST_IRAC_all_1panel(objid, ax=None, colors=['blue','maroon'], outputfile="", legend_loc=2, SEDtype='default', mode='1', textbox=True, title=False, xlow2=0.65, xsize2=0.3, ylow2=0.12, ysize2=0.3, magmin=30, magmax=23, legendfontsize='x-large', pztickfontsize='large', pztextfontsize='large',qsodir='/Users/khuang/Dropbox/Research/surfsup_dropbox/HIGHZ_ALL/eazy/QSO',stardir='/Users/khuang/Dropbox/Research/surfsup_dropbox/HIGHZ_ALL/eazy/STAR'):
   curdir = os.getcwd()
   # Plot both P(z) and SED for the object with ID objid
   if ax == None:
      fig = plt.figure(figsize=(8,6))
      ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
   else:
      fig = ax.get_figure()
   # First plot the best-fit SED as usual
   print "\nPlotting %s..." % objid
   ax1 = ax
   axes = []
   ax1, z_peak, props_hst, props_irac = plot_HST_IRAC_SED(objid, ax=ax1, 
         colors=colors, savefig=False, legend_loc=legend_loc, mode=mode)
   print "z_peak = %.2f" % z_peak
   # Plot best-fit QSO spectrum
   print "\nPlotting QSO spectrum...\n"
   os.chdir(qsodir)
   plt_qso = PlotEAZY()
   out_qso = plt_qso.plot_SED(objid, xmax=ax1.get_xlim()[1], ax=ax1, 
      plot_kwargs=dict(color='SeaGreen', ls='--', lw=1, label='QSO'))
   # Plot best-fit STAR spectrum
   print "\nPlotting STAR spectrum...\n"
   os.chdir(stardir)
   plt_star = PlotEAZY(with_specz=True)
   out_star = plt_star.plot_SED(objid, xmax=ax1.get_xlim()[1], ax=ax1,
      plot_kwargs=dict(color='0.1', ls=':', lw=1, label='STAR'))
   # re-set HST only line width
   for j in range(len(ax1.lines)):
      if ax1.lines[j].get_label() == 'HST only':
         # set line width for HST_only to 2
         ax1.lines[j].set_linewidth(2)  
         print "Line label: %s" % ax1.lines[j].get_label()
         break
   os.chdir(curdir)
   if title:
      title = 'Object %s [z_peak = %.3f]' % (str(objid), z_peak)
      ax1.set_title(title)
   else:
      ax1.set_title("")
   ax1.set_xlim(3000, 8.e4)
   ymin = pu.ABmag2uJy(magmin)
   ymax = pu.ABmag2uJy(magmax)
   ax1.set_ylim(ymin, ymax)
   # Set the y-axis tick labels
   mag_array = np.arange(magmin-1, magmax, -1)
   yticks = [pu.ABmag2uJy(x) for x in mag_array]
   ax1.set_yticks(yticks)
   ax1.set_yticklabels(['%d'%y for y in mag_array])
   ax1.set_xlabel('Observed $\lambda$ [$\mu\mathrm{m}$]')
   legend = ax1.legend(title=objid, loc=2, fontsize=legendfontsize)
   plt.setp(legend.get_title(), fontsize=legendfontsize)
   axes.append(ax1)
   plt.draw()  ## Important for setting the right axes boundaries for ax1!!
   print "Plotted SED."
   # Now plot P(z)
   # Read bbox properties AFTER plotting SEDs, because if ax is one of the 
   # AxesGrid axes, its boundaries won't be set properly before plotting the 
   # SED!
   bbox1 = ax1.get_position()
   xlow2, ylow2 = bbox1.p0 + bbox1.size * np.array([xlow2, ylow2])
   xsize2, ysize2 = bbox1.size * np.array([xsize2, ysize2])
   # print "Boundaries for ax1:"
   # print np.concatenate([bbox1.p0, bbox1.size])
   # print "New values for the boundaries of ax2:"
   print xlow2, ylow2, xsize2, ysize2
   ax2 = fig.add_axes([xlow2, ylow2, xsize2, ysize2])
   ax2, z_peak = plot_HST_IRAC_Pz(objid, ax=ax2, colors=colors, savefig=False,
                                 legend_loc=legend_loc, mode=mode)
   ## Plot P(z) for QSO templates as well...
   # os.chdir(qsodir)
   # ax2, z_peak_qso = plt_qso.plot_Pz(objid, ax=ax2, savefig=False, 
   #                   mode=mode, ls='--', lw=1, color='SeaGreen')
   # print "\nz_peak = %.3f (%.3f) for QSO (galaxy) templates.\n" % (z_peak_qso,z_peak)
   # os.chdir(curdir)
   ax2.legend().set_visible(False)  # turn off P(z) legend
   ax2.set_yticklabels([])
   ax2.set_ylabel(ax2.get_ylabel(), size=pztickfontsize)
   xlims2 = ax2.get_xlim()
   if z_peak < 9:
      xlims2 = (xlims2[0], 10.)
      xticks2 = [0, 2, 4, 6, 8, 10]
   else:
      xticks2 = [0, 3, 6, 9, 12]
   ax2.set_xlim(xlims2)
   # xticks2 = np.linspace(*ax2.get_xlim(), num=5)
   ax2.set_xticks(xticks2)
   ax2.set_xticklabels(['%d'%x for x in xticks2], size=pztickfontsize)
   ax2.set_xlabel('z', labelpad=0.01, size=pztickfontsize)
   ax2.set_title('')
   xtext = 0.1
   ytext = 0.9
   ha = 'left'
   if z_peak < 5:
      xtext = 0.9
      ha = 'right'
   ax2.text(xtext, ytext, '$z=%.1f$'% z_peak, size=pztextfontsize, ha=ha, 
            va='top', transform=ax2.transAxes)
   # NO best-fit SED texts
   if len(outputfile):
      plt.savefig(outputfile)
   axes.append(ax2)
   plt.draw()
   return axes

class MCEAZY(EAZY):
   def __init__(self, tempfluxunit='fnu', physfile=bc03_phys, with_specz=False):
      EAZY.__init__(self, tempfluxunit=tempfluxunit, physfile=physfile, 
                    paramfile='zphot.param.orig', with_specz=with_specz)
      # os.system('cp zphot.param zphot.param.orig') 
      # make a copy of the original zphot.param
      # Also read zphot.param
      self.param_lines = {}
      self.with_specz = with_specz
      with open('zphot.param.orig', 'rb') as f:
         lines = f.readlines()
         for i in range(len(lines)):
            if (lines[i][0] != "#") & (len(lines[i]) > 1):
               l = lines[i].split()
               self.param_lines[l[0]] = l[1:]

   def perturb_flux(self, niter):
      """
      Perturb the fluxes with errors and write a new catalog.
      niter: number of the desired iterations; because EAZY is too slow
             interpolating templates each time it runs, I'll just create a 
             giant catalog that has niter * nobj entries in it and fit just 
             once.
      """
      root = os.path.splitext(self.input_catalog)[0]
      new_catalog = root + "_%diters.cat" % niter
      # new_catalog = root + "_iter%d.cat" % niter
      # Now process the IDs to make them unique
      objids = []
      for i in range(niter):
         objids += [x + '-%d' % i for x in self.c._1]
      grand_output = [objids]
      if self.with_specz:
         grand_output += [np.tile(self.c._2, niter)]
      # grand_output = [self.c._1]  # the ID column
      # new_fluxerr = np.zeros((len(self.c), self.nbands))
      for i in range(self.nbands):
         if self.with_specz:
            flux_col = i * 2 + 3
            fluxerr_col = flux_col + 1
         else:
            flux_col = i * 2 + 2
            fluxerr_col = flux_col + 1
         old_flux = getattr(self.c, "_%d" % flux_col)
         old_fluxerr = getattr(self.c, "_%d" % fluxerr_col)
         old_fluxerr_floor = np.maximum(1.e-6, old_fluxerr)
         # Now duplicate the arrays niter times
         old_flux = np.tile(old_flux, niter)
         old_fluxerr = np.tile(old_fluxerr, niter)
         old_fluxerr_floor = np.tile(old_fluxerr_floor, niter)
         # this is to make np.random.normal work; objects who get this fluxerr
         # floor will not have their perturbed flux propogated to output
         perturb = (np.array(old_flux > 0))
         new_flux = np.where(perturb, 
                             np.random.normal(old_flux, old_fluxerr_floor), 
                             old_flux)
         grand_output += [["%e" % flux for flux in new_flux]]
         grand_output += [["%e" % fluxerr for fluxerr in old_fluxerr]]
         # new_flux[:,i] = np.where(perturb, 
         #                          np.random.normal(old_flux, old_fluxerr), 
         #                          old_flux)
         # new_fluxerr[:,i] = old_fluxerr
      # Now write to output
      grand_output = np.array(grand_output)
      np.savetxt(new_catalog, np.transpose(grand_output), "%s", 
                 header=self.c._header.strip(), delimiter="  ",
                 comments="")
      return new_catalog

   def rewrite_param(self, new_catalog):
      """
      Replace the old catalog with the new catalog in zphot.param.
      """
      print "Updating zphot.param..."
      with open('zphot.param', 'wb') as f:
         for k in self.param_lines:
            if k == 'CATALOG_FILE':
               value = self.param_lines[k]
               value[0] = new_catalog
            elif k in ["OBS_SED_FILE","POFZ_FILE"]:
               value = self.param_lines[k]
               value[0] = "n"
            else:
               value = self.param_lines[k]
            f.write(k + '  ' + ' '.join(value) + '\n')
      print "Done."

   def MCSampling(self, niter):
      """
      Runs Monte Carlo sampling of the input photometry catalog, runs EAZY
      for each realization, and collect the set of values for stellar
      population properties.
      We don't need to read the best-fit parameters for the input photometry.
      At the end of each iteration, read the best-fit parameters.
      """
      new_catalog = self.perturb_flux(niter)
      self.rewrite_param(new_catalog)
      # Run EAZY
      subprocess.call(["eazy"])
      # Now read output
      self.read_output()
      # for objid in self.id:
      print "Reading EAZY output..."
      output = {}
      for objid in self.c._1:
         print "Reading output for %s..." % objid
         if objid not in output.keys():
            output[objid] = dict(mass=[], log_age=[], sfr=[], tau=[], ebmv=[],
                                 mod_best=[], zbest=[])
         for i in range(niter):
            try:
               newid = objid + "-%d" % i
               (temp,zbest,mass,log_age,sfr,tau,ebmv,mod_best) = self.read_SED(newid)
               zbest = self.z_1[self.id==newid][0]
               output[objid]['mass'].append(mass)
               output[objid]['log_age'].append(log_age)
               output[objid]['sfr'].append(sfr)
               output[objid]['tau'].append(tau)
               output[objid]['ebmv'].append(ebmv)
               output[objid]['mod_best'].append(mod_best)
               output[objid]['zbest'].append(zbest)
            except:
               print "Could not read output for object %s!!" % newid
               continue
         print "Writing output for %s..." % objid
         output_file = objid + ".mc"
         if not os.path.exists(output_file):
            with open(output_file, 'wb') as f:
               f.write("# 1 NITER\n")
               f.write('# 2 ZBEST\n')
               f.write('# 3 SMASS\n')
               f.write('# 4 LOG_AGE\n')
               f.write('# 5 SFR\n')
               f.write('# 6 TAU\n')
               f.write('# 7 EBMV\n')
               f.write('# 8 MOD_BEST\n')
         if len(output[objid]['mass']):
            x = output[objid]
            with open(output_file, 'ab') as f:
               for j in range(len(output[objid]['mass'])):
                  line = '%d %f %e ' % (niter, x['zbest'][j], x['mass'][j])
                  line += '%e %f %f ' % (x['log_age'][j], x['sfr'][j], x['tau'][j])
                  line += '%f %d ' % (x['ebmv'][j], x['mod_best'][j])
                  f.write(line + '\n')
      # print "Finish iteration %d." % i

def fix_header(fname):
   """
   Fix the MC output catalog header.
   """
   header = ['# 1 NITER\n','# 2 ZBEST\n','# 3 SMASS\n','# 4 LOG_AGE\n',
             '# 5 SFR\n','# 6 TAU\n','# 7 EBMV\n','# 8 MOD_BEST\n']
   with open(fname) as f:
      lines = f.readlines()
   with open(fname, 'wb') as f2:
      for i in range(len(header)):
         f2.write(header[i])
      for j in range(len(lines)):
         if lines[j][0] != "#":
            f2.write(lines[j])
   print "Done."

def read_MCresults(objname, p=0.68, ebmv=-1, xgrid=200, with_specz=False, mu=1.0):
   """
   Read the Monte Carlo simulation results and calculate confidence intervals.
   """
   mc = sextractor('MonteCarlo/' + objname + '.mc')
   bf = EAZY(with_specz=with_specz)
   bf.read_output()
   R = bf.read_SED(objname)
   #bestprops = temp, z_best, mass_best, log_age_best, sfr_best, tau, ebmv, MOD_BEST
   bestprops = dict(zbest=R[1], smass=R[2], log_age=R[3], sfr=R[4])
   confint = {}
   # Also filter by E(B-V) if desired
   if ebmv == 'best':
      ebmv = R[6]
   if ebmv < 0:  # all E(B-V) values
      filt = np.ones(len(mc), 'bool')
   elif ebmv < 0.05:  # E(B-V) = 0
      filt = (mc.ebmv < 0.05)
   elif ebmv < 0.15:  # E(B-V) = 0.1
      filt = ((mc.ebmv > 0.05) & (mc.ebmv < 0.15))
   elif ebmv < 0.25:  # E(B-V) = 0.2
      filt = ((mc.ebmv > 0.15) & (mc.ebmv < 0.25))
   else:  # E(B-V) = 0.3
      filt = (mc.ebmv > 0.25)
   print "======================================="
   print "%d out of %d simulated points are used." % (np.sum(filt), len(filt))
   print "======================================="
   print "Best-fit E(B-V) = %.2f" % R[6]
   print ""
   if with_specz:
      properties = ['smass', 'log_age', 'sfr']
   else:
      properties = ['zbest','smass','log_age','sfr']
   print "Using magnification factor mu = %.2f:" % mu
   for x in properties:      
      if x == 'smass':
         d = distributions.Distribution1D(mc.smass[filt] / mu)
         print "%s: %f * 10^9 M_solar" % (x, bestprops[x] * 1.e-9 / mu)
         scale = 1.e-9
         x0 = bestprops[x] / mu
      elif x == 'log_age':
         d = distributions.Distribution1D(10.**(mc.log_age[filt]))
         print "%s: %f Myr" % (x, 10.**bestprops[x] / 1.e6)
         scale = 1.e-6
         x0 = (10.**bestprops[x])
      elif x == 'sfr':
         # d = distributions.Distribution1D(getattr(mc, x)[filt])
         d = distributions.Distribution1D(mc.sfr[filt] /mu)
         print "%s: %f" % (x, bestprops[x] / mu)
         scale = 1.0
         x0 = bestprops[x] / mu
      limits = d.conf_interval(xgrid=100, p=p, x0=x0, print_it=True,
                               scale=scale)
      confint[x] = limits
      print ""
   # Also print specific SFR
   d = distributions.Distribution1D((mc.sfr / mc.smass * 1.e9)[filt])
   x0 = bestprops['sfr']/bestprops['smass']*1.e9
   print "sSFR: %f M_solar/Gyr" % (x0)
   scale = 1.0
   limits = d.conf_interval(xgrid=xgrid, p=p, x0=x0, print_it=True,
                            scale=scale)
   print ""
   return confint

class Plot_MCEAZY(object):
   def __init__(self, objname):
      assert os.path.exists('hst_only')
      assert os.path.exists('with_irac')
      self.c1 = sextractor('hst_only/MonteCarlo/%s.mc' % objname)
      self.c2 = sextractor('with_irac/MonteCarlo/%s.mc' % objname)
      self.objname = objname

   def hist_MCEAZY(self, prop, bins1, bins2, ax=None, logprop=False, xlabel="", **hist_kwargs):
      """
      Plots the histograms of MC sampling results. Should be called in the 
      directory above both hst_only and with_irac.
      """
      if ax == None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
      assert hasattr(self.c1, prop), "Column '%s' not in the output catalog." % prop
      x1 = getattr(self.c1, prop)
      x2 = getattr(self.c2, prop)
      if logprop:
         x1 = np.log10(np.maximum(x1, 1e-3))
         x2 = np.log10(np.maximum(x2, 1e-3))
      h1 = ax.hist(x1, bins1, lw=2, color='0.7', label='hst_only',
                   histtype='stepfilled', normed=True, **hist_kwargs)
      h2 = ax.hist(x2, bins2, lw=2, color='red', label='with_irac',
                   histtype='step', normed=True, **hist_kwargs)
      ax.set_xlabel(xlabel)
      ax.set_title(self.objname)
      ax.legend(loc=2)
      return ax

   def hist_all(self, **hist_kwargs):
      """
      Plots the histogram for photo-z (top-left), stellar mass (top-right),
      SFR (bottom-left), and age (bottom-right).
      """
      fig = plt.figure(figsize=(12,9))
      ax1 = fig.add_subplot(221)
      self.hist_MCEAZY('zbest', np.arange(0.,11.,0.2), ax=ax1, 
                       xlabel='Photo-z')
      ax2 = fig.add_subplot(222)
      self.hist_MCEAZY('smass', np.arange(7.,11.,0.1), logprop=True, ax=ax2,
                       xlabel='log10(Stellar Mass)', **hist_kwargs)
      ax3 = fig.add_subplot(223)
      self.hist_MCEAZY('sfr', np.arange(-1.,4.,0.1), logprop=True, ax=ax3,
                       xlabel='log10(SFR)', **hist_kwargs)
      ax3.legend(loc=1)
      ax4 = fig.add_subplot(224)
      self.hist_MCEAZY('log_age', np.arange(7.,10.,0.1), ax=ax4,
                       xlabel='log10(Age)', **hist_kwargs)
      ax4.legend(loc=1)
      for ax in [ax1, ax2, ax3, ax4]:
         ax.set_title("")
      plt.suptitle(self.objname, size=28)
      return [ax1,ax2,ax3,ax4]

# class Plot_MCEAZY_wspaec(object):
#    def __init__(self, objname):
#       self.c1 = sextractor('MonteCarlo/%s.mc' % objname)
#       self.objname = objname

