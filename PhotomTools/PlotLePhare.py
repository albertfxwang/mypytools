#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from PhotomTools import LePhare as LP

class PlotLePhare(LP.LePhare):
   ## Le Phare output is object by object
   def __init__(self, objid, paramfile='', specdir='.'):
      curdir = os.getcwd()
      self.paramfile = paramfile
      self.readObjSpec(objid, specdir=specdir)  # read the best-fit SED; I don't trust the best-fit parameters here, though!
      # os.chdir(specdir)
      # self.readCatOut() # read the best-fit physical parameters
      self.objIndex = np.arange(len(self.data['IDENT']))[self.data['IDENT']==objid][0]
      # os.chdir(curdir)

   def plot_Pz(self, sedtype='GAL-1', outputdir='.', ax=None, savefig=True, xbox=0.05, ybox=0.95, txtPreFix="", txtProp={}, hatch='\\', **plot_kwargs):
      plot_kwargs_copy = default_plot_kwargs.copy()
      if len(plot_kwargs.keys()):
         for k in plot_kwargs.keys():
            plot_kwargs_copy[k] = plot_kwargs[k]
      if ax == None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
      txtProp_copy = default_txtProp.copy()
      if len(txtProp.keys()):
         for k in txtProp.keys():
            txtProp_copy[k] = txtProp[k]
      line = ax.plot(self.zarray, self.Pz, **plot_kwargs_copy)[0]
      ax.set_xlabel('LePhare Redshift')
      ax.set_ylabel(r'$P(z)$')
      # z_peak = self.photz  # from SPEC output; NEED VERIFICATION!!
      z_peak = self.data['Z_BEST'][self.objIndex]
      ax.set_title('Object %s [phot-z = %.3f]' % (self.objid, z_peak))
      # show the probability within the "1-sigma" range
      # bestprop = self.bestfitProps[sedtype] # read from SPEC output
      # zlow = bestprop['Zinf']
      # zhigh = bestprop['Zsup']
      zlow = self.data['Z_BEST68_LOW'][self.objIndex]
      zhigh = self.data['Z_BEST68_HIGH'][self.objIndex]
      pztext = 'Phot-z=%.3f (%.2f-%.2f)' % (self.photz, zlow, zhigh)
      if len(txtPreFix):
         pztext = txtPreFix + '\n' + pztext
      ax.text(xbox, ybox, pztext, transform=ax.transAxes,
              **txtProp_copy)
      # Also draw a shaded area within the "1-sigma" range
      ymax = ax.get_ylim()[1]
      if len(hatch):
         ax.fill_between([zlow, zhigh], [ymax, ymax], 0, hatch=hatch, 
                         facecolor='none', edgecolor=line.get_color())
      # ax.plot([zlow, zlow], [0, ymax], ls='--', lw=1.2, color='black')
      # ax.plot([zhigh, zhigh], [0, ymax], ls='--', lw=1.2, color='black')
      # fill_between somehow screws up the y-limits... so reset the limits
      ax.set_ylim(0, ymax)
      if savefig:
         fig.savefig("%s/%s_Pz.png" % (output_dir, self.objid))
      return ax

   def plot_photom(self, outputdir='.', ax=None, savefig=False, ebar_kwargs={}):
      if ax == None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
         ax.set_yscale('log')
         ax.set_xscale('log')
      if len(ebar_kwargs):
         for k in ebar_kwargs:
            default_ebar_kwargs[k] = ebar_kwargs[k]
      # First, plot the detections
      include = (self.photom['emag'] < 1.0) & (self.photom['Mag'] > 0)
      detect = np.logical_and((self.photom['emag'] > 0), include)
      nondetect = np.logical_and((self.photom['emag'] <= 0), include)
      ax.errorbar(self.photom['Lbd_mean'][detect], 
                  self.photom['fnu'][detect], 
                  yerr=self.photom['fnu_err'][detect], 
                  xerr=self.photom['Lbd_width'][detect]/2.,
                  fmt='s', 
                  **default_ebar_kwargs)
      # Then plot the upper limits
      ax.errorbar(self.photom['Lbd_mean'][nondetect],
                  self.photom['fnu'][nondetect],
                  xerr=self.photom['Lbd_width'][nondetect]/2.,
                  fmt=None, **default_ebar_kwargs)
      ax.scatter(self.photom['Lbd_mean'][nondetect], 
                 self.photom['fnu'][nondetect],
                 marker=downarrow, edgecolor=default_ebar_kwargs['ecolor'],
                 facecolor=default_ebar_kwargs['ecolor'],
                 s=(default_ebar_kwargs['ms']*2)**2, linewidths=1.2)
      ax.set_xlim(self.lambda_min, self.lambda_max)
      ax.set_ylim(ymax=self.photom['fnu'][include].max() * 10.)
      ax.set_xlabel('Wavelength (A)')
      ax.set_ylabel(r'$F_{\nu}$ [$\mu$Jy]')
      if savefig:
         fig.savefig("%s/%s_photom.png" % (output_dir, self.objid))
      return ax

   def plot_SED(self, sedtype='GAL-1', outputdir='.', ax=None, savefig=True, xbox=0.05, ybox=0.95, txtPreFix="", txtProp={}, plotGAL2=False, **plot_kwargs):
      plot_kwargs_copy = default_plot_kwargs.copy()
      txtProp_copy = default_txtProp.copy()
      if ax == None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
         ax.set_yscale('log')
         ax.set_xscale('log')
      if len(plot_kwargs.keys()):
         for k in plot_kwargs.keys():
            plot_kwargs_copy[k] = plot_kwargs[k]
      # replace box properties if new keywords are provided
      if len(txtProp.keys()):
         for k in txtProp.keys():
            txtProp_copy[k] = txtProp[k]
      l1 = ax.plot(self.wave, self.flux, 
                   label='z=%.3f' % self.bestfitProps['GAL-1']['Zphot'], 
                   **plot_kwargs_copy)
      if len(self.wave2) and plotGAL2:
         plot_kwargs_copy['ls'] = '--'
         l2 = ax.plot(self.wave2, self.flux2, 
                      label='z=%.3f' % self.bestfitProps['GAL-2']['Zphot'],
                      **plot_kwargs_copy)
         ax.legend(loc=4)
      nxticks = np.minimum(7, self.nfilters)
      log_lambda_ticks = np.linspace(np.log10(self.lambda_min / 0.8), 
                                     np.log10(self.lambda_max / 1.2), 
                                     nxticks+1)
      lambda_ticks = 10. ** log_lambda_ticks
      ax.set_xticks(lambda_ticks)
      # convert x-ticks into microns
      ax.set_xticklabels(map(lambda x: "%.1f" % (x / 1.e4), lambda_ticks))
      ax.set_xlabel(r"Wavelength [$\mu$m]")
      ax.set_ylabel(r'$F_{\nu}$ [$\mu$Jy]')
      ax.set_title('Object %s [phot-z = %.3f]' % (self.objid, self.photz), 
                   weight='bold')
      # Also display the best-fit SED properties... 
      for modtype in self.bestfitProps.keys():
         if (self.bestfitProps[modtype]['Chi2'] > 0) and (self.bestfitProps[modtype]['Chi2'] < self.bestfitProps['GAL-1']['Chi2']):
            print "****************************************"
            print "Warning: best-fit SED type for ID=%s is %s!!" % (self.objid, modtype)
            print "****************************************"
      # bestprop = self.bestfitProps[sedtype]
      logmass = self.data['MASS_BEST'][self.objIndex]
      # mass = scinote2exp('%e' % 10.**bestprop['Mass'])
      mass = scinote2exp('%e' % 10.**logmass)
      sedtext = "$M_{\mathrm{star}} = %s/\\mu\ \mathrm{[M_{\odot}]}$\n" % mass
      # sedtext = sedtext + "$E(B-V) = %.2f$\n" % bestprop['EB-V']
      sedtext = sedtext + "$E(B-V) = %.2f$\n" % self.data['EBV_BEST'][self.objIndex]
      # age = scinote2exp('%e' % 10.**bestprop['Age'])
      age = self.data['AGE_BEST'][self.objIndex]
      age = scinote2exp('%e' % age, nprec=2)
      sedtext = sedtext + "$\mathrm{Age} = %s\ \mathrm{[yrs]}$\n" % age
      # num_model = bestprop['Model']
      num_model = self.data['MOD_BEST'][self.objIndex]
      # tau = tau_array[num_model % len(tau_array)]
      # sedtext = sedtext + "$\\tau = %.1f\ \mathrm{[Gyrs]}$\n" % tau
      # sfr = scinote2exp('%e' % 10.**(bestprop['SFR']), nprec=2)
      logsfr = self.data['SFR_BEST'][self.objIndex]
      sfr = scinote2exp('%e' % 10.**logsfr)
      sedtext = sedtext + "$\mathrm{SFR} = %s/\\mu\ \mathrm{[M_{\odot}/yr]}$" % sfr
      # sedtext = sedtext + "$\chi_{\\nu}^2$ = %.2f" % bestprop['Chi2']
      # ------------------------------------------------------------
      # Find the best-fit metallicity... only works for BC03 models!
      # if (num_model <= 9):
      #    Z = 0.2  # 0.2 Z_solar (m42 models)
      # elif (num_model <= 18):
      #    Z = 0.4  # 0.4 Z_solar (m52 models)
      # else:
      #    Z = 1.0  # Z_solar (m62 models)
      # sedtext = sedtext + "$Z = %.2f\ Z_{\odot}$" % Z
      # ------------------------------------------------------------
      if len(txtPreFix):
         sedtext = txtPreFix +'\n' + sedtext
      ax.text(xbox, ybox, sedtext, transform=ax.transAxes, **txtProp_copy)
      if savefig:
         fig.savefig("%s/%s_SED.png" % (output_dir, self.objid))
      return ax

   def plot_all(self, axes=None, outputdir='.', objname="", savefig=True, xbox=0.05, ybox=0.95, txtPreFix="", txtProp={}, ebar_kwargs={}, SED_plot_kwargs={'lw':1.4, 'color':'black'}, Pz_plot_kwargs={'lw':2}):
      if axes == None:
         fig = plt.figure(figsize=(10,12))
         # Top panel: photometry + SED
         ax1 = fig.add_subplot(211)
         ax2 = fig.add_subplot(212)
         ax1.set_yscale('log')
         ax1.set_xscale('log')
      else:
         assert len(axes)==2, "Please provide a list of two matplotlib axes for the axes argument."
         ax1, ax2 = axes
      ax1 = self.plot_photom(outputdir=outputdir, ax=ax1, savefig=False,
                             ebar_kwargs=ebar_kwargs)
      ax1 = self.plot_SED(outputdir=outputdir, ax=ax1, savefig=False, 
                          txtProp=txtProp, txtPreFix=txtPreFix, 
                          xbox=xbox, ybox=ybox, **SED_plot_kwargs)
      ax2 = self.plot_Pz(outputdir=outputdir, ax=ax2, savefig=False,
                         xbox=xbox, ybox=ybox, txtProp=txtProp, 
                         txtPreFix=txtPreFix, **Pz_plot_kwargs)
      if len(objname):
         ax1.set_title('Object %s [photo-z = %.3f]' % (objname, self.photz))
         ax2.set_title('Object %s [photo-z = %.3f]' % (objname, self.photz))
      plt.draw()
      plt.savefig("%s/%s_SED_Pz_lephare.png" % (outputdir, self.objid))
      return ax1, ax2

def plot_HST_IRAC_all(objid, cluster_name, objname="", colors=['blue','red'], savefig=True, legend_loc=2, outputdir='.', outputname=""):
   # Use objid to find the LePhare output spec file (Idxxxxxxxxxx.spec)
   # objname will appear as the name in the figure title
   # colors[0] for HST_only, and colors[1] for with_IRAC
   # Must call this in the directory above hst_only/ and with_irac/
   specfile = "Id%09d.spec" % objid
   objid_sq = int('10' + str(objid))
   specfile_sq = "Id%09d.spec" % objid_sq
   fig = plt.figure(figsize=(10,12))
   ax1 = fig.add_subplot(211)
   ax1.set_xscale('log')
   ax1.set_yscale('log')
   ax2 = fig.add_subplot(212)
   print "Reading hst_only..."
   hstPlot = PlotLePhare(objid, 'lephare_%s_hst_only_bc03.param'%cluster_name, 
                         specdir='hst_only')
   # hstPlot_sq = PlotLePhare(objid_sq, 'hst_only')
   # hstPlot = PlotLePhare('hst_only/%s' % specfile)
   print "Reading with_irac..."
   iracPlot = PlotLePhare(objid,'lephare_%s_with_irac_bc03.param'%cluster_name, 
                          specdir='with_irac')
   print "Reading star/qso..."
   iracPlot_sq = PlotLePhare(objid_sq, 
                             'lephare_%s_with_irac_bc03_starqso.param'%cluster_name, 
                             specdir='with_irac')
   # iracPlot = PlotLePhare('with_irac/%s' % specfile)
   # First plot the photometry points with iracPlot only 
   ax1 = iracPlot.plot_photom(ax=ax1, savefig=False, 
                              ebar_kwargs={'mec':'black','ecolor':'black'})
   # Now plot the best-fit SED from HST_only
   ax1 = hstPlot.plot_SED(ax=ax1, savefig=False, xbox=0.05, ybox=0.95, 
                          txtPreFix='HST ONLY:', 
                          txtProp={'size':'large','bbox':dict(boxstyle='round,pad=0.5',facecolor='LightCyan')}, 
                          plotGAL2=False, color=colors[0], ls='--')
   # Plot the best-fit SED from with_IRAC
   ax1 = iracPlot.plot_SED(ax=ax1, savefig=False, xbox=0.95, ybox=0.05, 
                           txtPreFix='WITH IRAC:',
                           txtProp={'va':'bottom','ha':'right','size':'large','bbox':dict(boxstyle='round,pad=0.5',facecolor='LightPink')},
                           color=colors[1])
   # Plot the P(z) curves from HST only
   ax2 = hstPlot.plot_Pz(ax=ax2, savefig=False, xbox=0.05, ybox=0.95, 
                         txtPreFix='HST ONLY:', 
                         txtProp={'size':'x-large',
                         'bbox':dict(boxstyle='round,pad=0.3',facecolor='LightCyan')},
                         hatch='', color=colors[0], lw=2)
   Pzmax_hst = hstPlot.Pz.max()
   # Plot the P(z) curves from with_IRAC
   ax2 = iracPlot.plot_Pz(ax=ax2, savefig=False, xbox=0.05, ybox=0.8,
                          txtPreFix='WITH IRAC:',
                          txtProp={'size':'x-large', 'bbox':dict(boxstyle='round,pad=0.3',facecolor='LightPink')},
                          hatch='', color=colors[1], lw=2)
   # Check if STAR or QSO have better chi2 for this object than GAL-1
   for sedtype in iracPlot_sq.bestfitProps.keys():
      if iracPlot_sq.bestfitProps[sedtype]['Nline'] > 0:
         if iracPlot_sq.bestfitProps[sedtype]['Chi2'] < iracPlot.bestfitProps['GAL-1']['Chi2']:
            print "**********************************************************"
            print "Warning: %s has lower chi2 than galaxy for object %d!!" % (   sedtype, objid)
            print "**********************************************************"
   Pzmax_irac = iracPlot.Pz.max()
   ymax = np.maximum(Pzmax_hst, Pzmax_irac) * 1.2
   ax2.set_ylim(0, ymax)
   # Put hatch around the confidence interval
   bp_hst = hstPlot.bestfitProps['GAL-1']
   bp_irac = iracPlot.bestfitProps['GAL-1']
   ax2.fill_between([bp_hst['Zinf'], bp_hst['Zsup']], ymax*5, hatch='\\',
                  edgecolor=colors[0], facecolor='none')
   ax2.fill_between([bp_irac['Zinf'], bp_irac['Zsup']], ymax*5, hatch='/',
                    edgecolor=colors[1], facecolor='none')
   ax2.set_ylim(0, ymax)
   if len(objname):
      ax1.set_title('Object %s' % (objname), size=28)
      # ax2.set_title('Object %s' % (objname))
   else:
      ax1.set_title('Object ID=%d' % objid, size=28)
      # ax2.set_title('Object ID=%d' % objid)
   # ax1, ax2 = hstPlot.plot_all(axes=[ax1,ax2], objname=objname, savefig=False,
   #                         ebar_kwargs={'mec':colors[0],'ecolor':colors[0]},
   #                         SED_plot_kwargs={'color':colors[0]})
   
   # ax1, ax2 = iracPlot.plot_all(axes=[ax1,ax2],objname=objname,savefig=False,
   #                         ebar_kwargs={'mec':colors[1],'ecolor':colors[1]},
   #                         SED_plot_kwargs={'color':colors[1]})
   plt.draw()
   if savefig:
      if len(outputname):
         outloc = '%s/%s' % (outputdir, outputname)
      else:
         outloc = '%s/obj%d_lephare_hst_irac.png' % (outputdir, objid)
      plt.savefig(outloc)
   return ax1, ax2

