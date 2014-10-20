#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import photomUtils as pu
from scipy.integrate import cumtrapz

# Some default matplotlib stuff
default_ebar_kwargs={'ms':12, 'marker':'.', 'mec':'black', 'mfc':'black', 
                     'ecolor':'black', 
                     'mew':1.2, 'capsize':9, 'capthick':1.5}
default_plot_kwargs={'lw':1.3}
default_txtProp = {'va':'top', 'ha':'left', 'multialignment':'left',
                  'bbox':dict(boxstyle='round,pad=0.3',facecolor='LightCyan'),
                  'size':'x-large'}
# a symbol for magnitude upper limits
downarrow = [(-2,0),(2,0),(0,0),(0,-4),(-2,-2),(2,-2),(0,-4),(0,0)]

def scinote2exp(scinote):
   # convert scientific notation in the format of 1.0e10 into (1.0, 10)
   n = scinote.split('e')
   if np.abs(int(n[1])) <= 2:
      return "%.3f" % float(scinote)
   else:
      # return float(n[0]), int(n[1])
      return "%.2f \\times 10^{%d}" % (float(n[0]), int(n[1]))

class PlotLePhare(object):
   ## Le Phare output is object by object
   def __init__(self, specfile):
      # read the Le Phare output *.spec file
      f = open(specfile, 'rb')
      lines = f.readlines()
      f.close()
      self.bestfitProps = {}
      readPz = False
      readPhotom = False
      SEDstart = len(lines)
      nphotom = 0
      self.photom = {}
      photomKeys = []
      for i in range(len(lines)):
         if lines[i].startswith('# Ident'):
            # First two lines are the Object ID, spec-z (if available), and photo-z
            print "Read object ID & redshifts..."
            l_identity = lines[i+1].split()
            self.objid = l_identity[0]
            self.specz = float(l_identity[1])
            self.photz = float(l_identity[2])
            continue
         elif lines[i].startswith('# Mag'):
            nphotom = len(lines[i].split()[1:])
            self.photomKeys = lines[i].split()[1:]
            for k in self.photomKeys:
               self.photom[k] = []
            # Also add Fnu in micro-Jansky
            self.photom['fnu'] = []
            continue
         elif lines[i].startswith('FILTERS'):
            print "Read the number of filters..."
            # Then read the number of filters
            self.nfilters = int(lines[i].split()[1])
            continue
         elif lines[i].startswith('PDF'):
            print "Read the number of P(z) steps..."
            # Read the number of P(z) steps
            self.zsteps = int(lines[i].split()[1])
            continue
         elif lines[i].startswith('# Type'):
            print "Reading best-fit properties..."
            print lines[i]
            attributes = lines[i].split()[1:]
            nattr = len(attributes)
            # read the values for each SED type
            j = 1
            while 1:
               l_attr = lines[i+j].split()
               if len(l_attr) != nattr:
                  break
               self.bestfitProps[l_attr[0]] = dict(zip(attributes[1:], l_attr[1:]))
               # convert into either integer or float
               for k in self.bestfitProps[l_attr[0]]:
                  if k.lower() in ['nline', 'model', 'library', 'nband', 'extlaw']:
                     self.bestfitProps[l_attr[0]][k] = int(self.bestfitProps[l_attr[0]][k])
                  else:
                     self.bestfitProps[l_attr[0]][k] = float(self.bestfitProps[l_attr[0]][k])
               j += 1
            print "A total of %d SED types read." % j
            continue
         elif (len(lines[i].split()) == nphotom) and (readPhotom == False):
            # Read the object photometry
            if not readPhotom:
               print "Read object photometry in %d filters..." % self.nfilters
               readPhotom = True
            for j in range(i, i+self.nfilters):
               photomList = [float(x) for x in lines[j].split()]
               for k in range(nphotom):
                  self.photom[self.photomKeys[k]] += [photomList[k]]
            # convert lists into numpy arrays
            for k in self.photom:
               self.photom[k] = np.array(self.photom[k])
            # calculate fnu in micro-Jansky from AB mag
            self.photom['fnu'] = pu.ABmag2uJy(self.photom['Mag'])
            # also calculate flux errors --- first calculate S/N
            # If mag_err < 0, then mag is upper limit
            SN = pu.magerr2sn(self.photom['emag'])
            self.photom['fnu_err'] = np.where(SN > 0, self.photom['fnu'] / SN,
                                              self.photom['fnu'])
            # Estimate the range in wavelength to show in the plots
            self.lambda_max = (self.photom['Lbd_mean'][-1] + self.photom['Lbd_width'][-1])  * 1.2
            self.lambda_min = (self.photom['Lbd_mean'][0] - self.photom['Lbd_width'][0]) * 0.8
            continue
         elif not readPz and len(lines[i].split()) == 2:
            print "Reading P(z)..."
            # Read P(z) for the next self.zsteps lines
            PzLines = [lines[j].split() for j in range(i, i+self.zsteps)]
            PzBlock = [[float(l[0]), float(l[1])] for l in PzLines]
            PzBlock = np.array(PzBlock)
            self.zarray = PzBlock[:,0]
            self.Pz = PzBlock[:,1]
            SEDstart = i + self.zsteps
            probIntegrated = cumtrapz(self.Pz, x=self.zarray)[-1]
            # prob = self.Pz / probIntegrated
            self.Pz = self.Pz / probIntegrated
            readPz = True
         elif (readPz == True) and (i == SEDstart):
            print "Reading SED for GAL-1..."
            # Read the best-fit SED flux (in AB mag)
            self.wave = []
            self.flux = []
            self.wave2 = []
            self.flux2 = []
            self.fluxunit = 'uJy'
            j = 0
            # First read GAL-1
            # while i + j < len(lines):
            while j < self.bestfitProps['GAL-1']['Nline']:
               if len(lines[i+j].split()) == 2:
                  l = lines[i+j].split()
                  self.wave += [float(l[0])]
                  self.flux += [pu.ABmag2uJy(float(l[1]))]
                  j = j + 1
            if self.bestfitProps['GAL-2']['Nline'] > 0:
               # Also read a second galaxy SED
               print "Reading SED for GAL-2..."
               j2 = 0
               while j2 < self.bestfitProps['GAL-2']['Nline']:
                  if len(lines[i+j+j2].split()) == 2:
                     # print i+j+j2
                     l = lines[i+j+j2].split()
                     self.wave2 += [float(l[0])]
                     self.flux2 += [pu.ABmag2uJy(float(l[1]))]
                     j2 = j2 + 1
            self.wave = np.array(self.wave)
            self.flux = np.array(self.flux)
            if len(self.wave2):
               self.wave2 = np.array(self.wave2)
               self.flux2 = np.array(self.flux2)
            continue

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
      z_peak = self.photz  # NEED VERIFICATION!!
      ax.set_title('Object %s [phot-z = %.3f]' % (self.objid, z_peak))
      # show the probability within the "1-sigma" range
      bestprop = self.bestfitProps[sedtype]
      zlow = bestprop['Zinf']
      zhigh = bestprop['Zsup']
      prob1sig = bestprop['PDF']
      pztext = 'Phot-z=%.3f (%.2f-%.2f)' % (self.photz, zlow, zhigh)
      # pztext = 'Prob(%.2f-%.2f) = %.2f' % (zlow, zhigh, (prob1sig/100.))
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
      bestprop = self.bestfitProps[sedtype]
      mass = scinote2exp('%e' % 10.**bestprop['Mass'])
      sedtext = "$M_{\mathrm{star}} = %s/\\mu$ [$M_{\odot}$]\n" % mass
      sedtext = sedtext + "$E(B-V) = %.2f$\n" % bestprop['EB-V']
      age = scinote2exp('%e' % 10.**bestprop['Age'])
      sedtext = sedtext + "$\mathrm{Age} = %s$ [yrs]\n" % age
      sfr = scinote2exp('%e' % bestprop['SFR'])
      sedtext = sedtext + "$\mathrm{SFR} = %s/\\mu$ [$M_{\odot}$/yr]\n" % sfr
      sedtext = sedtext + "$\chi_{\\nu}^2$ = %.2f" % bestprop['Chi2']
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

def plot_HST_IRAC_all(objid, objname="", colors=['blue','red'], savefig=True, legend_loc=2, outputdir='.', outputname=""):
   # Use objid to find the LePhare output spec file (Idxxxxxxxxxx.spec)
   # objname will appear as the name in the figure title
   # colors[0] for HST_only, and colors[1] for with_IRAC
   # Must call this in the directory above hst_only/ and with_irac/
   specfile = "Id%09d.spec" % objid
   fig = plt.figure(figsize=(10,12))
   ax1 = fig.add_subplot(211)
   ax1.set_xscale('log')
   ax1.set_yscale('log')
   ax2 = fig.add_subplot(212)
   hstPlot = PlotLePhare('hst_only/%s' % specfile)
   iracPlot = PlotLePhare('with_irac/%s' % specfile)
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
   Pzmax_irac = iracPlot.Pz.max()
   ymax = np.maximum(Pzmax_hst, Pzmax_irac) * 1.2
   print "Pz max:", ymax
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

