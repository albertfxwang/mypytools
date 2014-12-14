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

eazy_path = '/Users/khuang/bin/eazy'
# a symbol for magnitude upper limits
downarrow = [(-2,0),(2,0),(0,0),(0,-4),(-2,-2),(0,-4),(2,-2),(0,-4),(0,0)]
pc_cm = 3.0857e18   # parsec in cm
area_10pc = 4 * np.pi * (10*pc_cm)**2  # in cm^2
L_solar = 3.826e33  # in erg/s
bc03_phys = '/Users/khuang/eazy-1.00/templates/BC03/bc03_m42_tau.phys'
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
         assert len(self.c._d) % 2 == 0, "Number of columns doesn't seem right... missing some fluxes or flux errors?"
      else:
         assert len(self.c._d) % 2 == 1, "Number of columns doesn't seem right... missing some fluxes or flux errors?"
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

   def read_SED(self, objid):
      """
      Read the best-fit SED and determine the best-fit physical properties.
      """
      temp = sextractor("OUTPUT/%s.temp_sed" % str(objid))
      # temp flux is in uJy, the same as catalog flux
      # Also read the best-fit template number
      header = temp._header.split('\n')
      h2 = header[1].split('templ:norm')[-1].split(':')
      MOD_BEST = int(h2[0])
      # print "MOD_BEST: ", (MOD_BEST)
      ## MOD_BEST is the model number in TEMPLATES_FILE; the file order has
      ## to match between TEMPLATES_FILE and physfile for the best-fit physical
      ## Calculate the normalization factor, from template to observed flux
      ## The goal is to match their L_nu in rest-frame
      sp_obs = S.ArraySpectrum(temp._1, temp._2, fluxunits='fnu')  # in uJy
      sp_obs = sp_obs * 1.e-29  # convert from uJy to erg/s/cm^2/Hz
      sp_mod = S.FileSpectrum(os.path.split(self.physfile)[0] + '/' + self.phys.filename[MOD_BEST-1])  # in flam
      sp_mod.convert('fnu')  # convert flux unit to fnu
      zpeak = self.z_1[self.id==objid][0]
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
      return temp, mass_best, log_age_best, sfr_best, tau, ebmv, MOD_BEST

class PlotEAZY(EAZY):
   def __init__(self, tempfluxunit='fnu', physfile=bc03_phys, wavelim=[3000.,8e4], with_specz=False):
      EAZY.__init__(self, tempfluxunit=tempfluxunit, physfile=physfile, with_specz=with_specz)
      self.lambda_max = 30000.  # default value for maximum lambda in plots
      self.lambda_min = 3000.
      self.flux_ujy_max = 1.0
      self.physfile = physfile
      self.phys = sextractor(physfile)
      self.wavelim = wavelim

   def plot_Pz(self, objid, ax=None, savefig=True, mode='a', **plot_kwargs):
      assert objid in self.objid, "Object ID %s does not exist in catalog." % str(objid)
      pz = sextractor("OUTPUT/%s.pz" % str(objid))
      if ax == None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
      prob = np.exp(-0.5 * pz._2)
      probIntegrated = cumtrapz(prob, x=pz._1)[-1]
      prob = prob / probIntegrated
      ax.plot(pz._1, prob, **plot_kwargs)
      ax.set_xlabel('Redshift')
      ax.set_ylabel(r'$P(z)$')
      if mode == 'a':
         z_peak = self.z_a[self.id==objid][0]
      elif mode=='1':
         z_peak = self.z_1[self.id==objid][0]
      ax.set_title('Object %s [z_peak = %.3f]' % (str(objid), z_peak))
      if savefig:
         fig.savefig("OUTPUT/%s_Pz.png" % str(objid))
      return ax, z_peak

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
      (temp,mass_best,log_age_best,sfr_best,tau,ebmv,mod_best) = self.read_SED(objid)
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
      lambda_ticks = [5000.,10000.,20000.,30000.,40000.,50000.]
      ax.set_xticks(lambda_ticks)
      ax.set_xticklabels(map(lambda x: "%.1f" % (x/1.e4), lambda_ticks))
      ax.set_xlabel(r"$\lambda$ [$\mu$m]")
      if not len(ax.get_ylabel()):
         ax.set_ylabel(r"$F_{\nu}$ [$\mu$Jy]")
      if mode == 'a':
         z_peak = self.z_a[self.id==objid][0]
      elif mode == '1':
         z_peak = self.z_1[self.id==objid][0]
      ax.set_title('Object %s [z_peak = %.3f]' % (str(objid), z_peak), weight='bold')
      if savefig:
         plt.savefig('OUTPUT/%s_SED.png' % str(objid))
      return ax, z_peak, [mass_best, 10.**log_age_best, sfr_best, tau, ebmv]

   def plot_all(self, objid, objname="", savefig=True, legend_loc=1, mode='a'):
      fig = plt.figure(figsize=(10,12))
      ax1 = fig.add_subplot(211)
      ax1.set_yscale('log')
      ax1.set_xscale('log')
      ax1 = self.plot_photom(objid, ax=ax1)
      ax1, zp, props = self.plot_SED(objid, ax=ax1, mode=mode)
      ax2 = fig.add_subplot(212)
      ax2, zp = self.plot_Pz(objid, ax=ax2, savefig=False, mode=mode)
      if len(objname):
         ax1.set_title('Object %s [z_peak = %.3f]' % (objname, zp))
         ax2.set_title('Object %s [z_peak = %.3f]' % (objname, zp))
      plt.draw()
      plt.savefig("%s_SED_Pz.png" % str(objid))

class PlotEAZY_wspecz(PlotEAZY):
   def __init__(self, tempfluxunit='fnu', physfile=bc03_phys, wavelim=[3000.,8e4]):
      PlotEAZY.__init__(self, tempfluxunit=tempfluxunit, physfile=physfile, 
                        wavelim=wavelim, with_specz=True)
   def plot_all(self, objid, z_spec, xmax=None, ax=None, savefig=True, textbox=True, plot_kwargs={'color':'black'}):
      fig = plt.figure(figsize=(10,8))
      ax1 = fig.add_subplot(111)
      ax1.set_yscale('log')
      ax1.set_xscale('log')
      ax1 = self.plot_photom(objid, ax=ax1)
      ax1, zp, props = self.plot_SED(objid, ax=ax1, mode='1')
      ax1.set_title('Object %s [z_spec = %.2f]' % (objid, z_spec))
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
         ax1.text(0.95, 0.05, sedtext, size='large', ha='right', va='bottom',
                 transform=ax1.transAxes, multialignment='left',
                 bbox=dict(boxstyle='round,pad=0.3',facecolor='LightPink'))
      plt.draw()
      if savefig:
         plt.savefig("%s_SED_Pz.png" % str(objid))
      return ax1

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

def plot_HST_IRAC_all(objid, colors=['blue','maroon'], outputfile="", legend_loc=2, SEDtype='default', mode='1', textbox=True):
   # Plot both P(z) and SED for the object with ID objid
   fig = plt.figure(figsize=(10,12))
   ax1 = fig.add_subplot(2, 1, 1)
   ax1, z_peak, props_hst, props_irac = plot_HST_IRAC_SED(objid, ax=ax1, 
         colors=colors, savefig=False, legend_loc=legend_loc, mode=mode)
   title = 'Object %s [z_peak = %.3f]' % (str(objid), z_peak)
   ax1.set_title(title)
   # if objname:
   #    ax1.set_title(title + '\n(%s SED)' % SEDtype)
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

class MCEAZY(EAZY):
   def __init__(self, tempfluxunit='fnu', physfile=bc03_phys):
      EAZY.__init__(self, tempfluxunit=tempfluxunit, physfile=physfile, 
                    paramfile='zphot.param.orig')
      # os.system('cp zphot.param zphot.param.orig') 
      # make a copy of the original zphot.param
      # Also read zphot.param
      self.param_lines = {}
      with open('zphot.param.orig', 'rb') as f:
         lines = f.readlines()
         for i in range(len(lines)):
            if (lines[i][0] != "#") & (len(lines[i]) > 1):
               l = lines[i].split()
               self.param_lines[l[0]] = l[1:]

   def perturb_flux(self, niter):
      """
      Perturb the fluxes with errors and write a new catalog.
      """
      root = os.path.splitext(self.input_catalog)[0]
      new_catalog = root + "_iter%d.cat" % niter
      # new_flux = np.zeros((len(self.c), self.nbands))
      # new_fluxerr = np.zeros((len(self.c), self.nbands))
      grand_output = [self.c._1]  # the ID column
      # new_fluxerr = np.zeros((len(self.c), self.nbands))
      for i in range(self.nbands):
         flux_col = i * 2 + 2
         fluxerr_col = flux_col + 1
         old_flux = getattr(self.c, "_%d" % flux_col)
         old_fluxerr = getattr(self.c, "_%d" % fluxerr_col)
         old_fluxerr_floor = np.maximum(1.e-6, old_fluxerr)
         # this is to make np.random.normal work; objects who get this fluxerr
         # floor will not have their perturb flux propogated to output
         perturb = (old_flux > 0)
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
            else:
               value = self.param_lines[k]
            f.write(k + '  ' + ' '.join(value) + '\n')
      print "Done."

   def MCSampling(self, nstart=1, nfinish=5):
      """
      Runs Monte Carlo sampling of the input photometry catalog, runs EAZY
      for each realization, and collect the set of values for stellar
      population properties.
      We don't need to read the best-fit parameters for the input photometry.
      At the end of each iteration, read the best-fit parameters.
      """
      for i in range(nstart, nfinish+1):
         print "Start iteration %d..." % i
         new_catalog = self.perturb_flux(i)
         self.rewrite_param(new_catalog)
         # Run EAZY
         subprocess.call(["eazy"])
         # Now read output
         self.read_output()
         for objid in self.id:
            output_file = objid + ".mc"
            try:
               (temp,mass,log_age,sfr,tau,ebmv,mod_best) = self.read_SED(objid)
               zbest = self.z_1[self.id==objid][0]
            except:
               print "Could not read output for object %s!!" % objid
               continue
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
            with open(output_file, 'ab') as f:
               f.write('%s  %f  %e  %e  %f  %f  %f  %d\n' % (i, zbest, mass, log_age, sfr, tau, ebmv, mod_best)) 
         print "Finish iteration %d." % i

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

class Plot_MCEAZY(object):
   def __init__(self, objname):
      assert os.path.exists('hst_only')
      assert os.path.exists('with_irac')
      self.c1 = sextractor('hst_only/MonteCarlo/%s.mc' % objname)
      self.c2 = sextractor('with_irac/MonteCarlo/%s.mc' % objname)
      self.objname = objname

   def hist_MCEAZY(self, prop, bins, ax=None, logprop=False, xlabel="", **hist_kwargs):
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
      h1 = ax.hist(x1, bins, lw=2, color='blue', label='hst_only',
                   histtype='step', normed=True, **hist_kwargs)
      h2 = ax.hist(x2, bins, lw=2, color='red', label='with_irac',
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
