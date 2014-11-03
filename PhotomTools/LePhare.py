#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import photomUtils as pu
from scipy.integrate import cumtrapz
from pygoods import sextractor
import copy, os, sys, subprocess

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
# a list of model tau's (for exponentially declining SFH)
tau_array = [0.1, 0.3, 1.0, 2.0, 3.0, 5.0, 10.0, 15.0, 30.0]
numTau = len(tau_array)

def scinote2exp(scinote, nprec=3):
   # convert scientific notation in the format of 1.0e10 into (1.0, 10)
   n = scinote.split('e')
   fltstr = "%%.%df" % nprec
   if np.abs(int(n[1])) <= 2:
      return fltstr % float(scinote)
   else:
      # return float(n[0]), int(n[1])
      return "%.2f \\times 10^{%d}" % (float(n[0]), int(n[1]))

# Write a parent class for LePhare-related tasks, then subclasses for plotting 
# and running MC sampling of photometry
class LePhare(object):
   def __init__(self, paramfile):
      self.paramfile = paramfile
      self.getCatOut()

   def getCatOut(self):
      assert len(self.paramfile) > 0, "No parameter file provided..."
      with open(self.paramfile, 'rb') as f:
         lines = f.readlines()
         for line in lines:
            if line.startswith('CAT_OUT'):
               l = line.split()
               self.catout = l[1]

   def readCatOut(self):
      """
      Read the output catalog (CAT_OUT).
      """
      assert len(self.catout) > 0
      columns = []
      data = []
      self.data = {}
      with open(self.catout, 'rb') as f:
         lines = f.readlines()
      # n_hdr = 0
      # n_line = 0
      for i in range(len(lines)):
         if lines[i].startswith('# AB'):
            # read the filter names
            bands = lines[i].strip().split()[2:]
         elif lines[i].startswith('# Output format'):
            # identified the start of header lines
            n_hdr = i + 1
            break
      # print "bands:", bands
      for i in range(n_hdr, len(lines)):
         # read the column names
         if lines[i].startswith('# '):
            l_hdr = lines[i][1:].strip()
            l_hdr = l_hdr.split(',')[:-1]
            for l in l_hdr:
               colname, colnum = l.split()
               colnum = int(colnum)
               # test if the column name jumps in number
               if colname == 'STRING_INPUT':
                  pass
               elif colname != "MAG_ABS()":
                  columns += [colname]
               else:
                  # these are the absolute magnitudes in all the filters
                  for b in bands:
                     columns += ['ABSMAG_' + b]
         elif lines[i].startswith('###'):
            # this is the line separating header and data
            n_data = i + 1
            break
      # Now start reading data
      for i in range(n_data, len(lines)):
         l = lines[i].split()  # split the fields for each object
         if len(l) == 0:
            break
         outlist = []
         for entry in l:
            try:
               x = int(entry)
            except:
               x = float(entry)
            outlist += [x]
         data += [outlist]
      # Now collect the results
      for j in range(len(columns)):
         self.data[columns[j]] = np.array(map(lambda x: x[j], data))

   def getSpecFile(self, objid):
      specfile = "Id%09d.spec" % objid
      return specfile

   def readObjSpec(self, objid, specdir='.'):
      # construct the output spec file name from object ID, and delegate the 
      # work to self.readSpec
      self.objid = objid
      curdir = os.getcwd()
      os.chdir(specdir)
      specfile = self.getSpecFile(objid)
      self.readSpec(specfile)
      self.getCatOut()
      self.readCatOut()
      os.chdir(curdir)

   def readSpec(self, specfile):
      """
      Read the Le Phare output *.spec file; will set the following attributes:
      self.bestfitProps  --- best-fit stellar population properties
      self.photom        --- the photometry from the input catalog
      self.objid
      self.specz         --- spec-z if there was any in the input catalog
      self.photz         --- best-fit photo-z (peak of P(z)?)
      self.photomKeys    --- filter names(?)
      self.nfilters      --- number of filters
      self.zsteps        --- number of redshift steps
      self.zarray        --- the redshift grid over which P(z) is defined
      self.Pz            --- the P(z) curve
      self.lambda_max
      self.lambda_min    --- the range of wavelength for plotting
      self.wave 
      self.flux 
      self.wave2
      self.flux2         --- the wavelengths & fluxes of the 1st and 2nd 
                             best SED models
      """
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

class MCLePhare(LePhare):
   """
   Runs Monte Carlo sampling to determine error bars for SED fitting.
   """
   def __init__(self, paramfile, inputcat):
      self.paramfile = paramfile
      self.inputcat = inputcat
      self.catHeader = ""
      self.catColumns = []
      # self.getCatOut()

   def readCatalog(self):
      # Read the input photometry catalog to LePhare, because we will 
      # perturb the magnitudes later
      # read the filters in the catalog
      with open(self.inputcat, 'rb') as f:
         lines = f.readlines()
         for l in lines:
            if l[0] == '#':
               l2 = l.split()
               if l2[1] == 'ID':
                  # this is the header line
                  self.catHeader = copy.copy(l)
                  for l3 in l2[1:]:
                     self.catColumns += [l3.lower()]
                  break
      c = sextractor(self.inputcat)
      for i in range(len(self.catColumns)):
         setattr(self, self.catColumns[i], getattr(c, '_%d'%(i+1)))

   def MCSampling(self, N):
      """
      Runs Monte Carlo sampling of the input photometry catalog, runs LePhare
      for each realization, and collect the set of values for stellar
      population properties.
      We don't need to read the best-fit parameters for the input photometry.
      """
      # read the photometry of input catalog
      self.readCatalog()
      MainDIR = os.getcwd()
      root = os.path.splitext(self.inputcat)[0]
      assert len(self.catHeader) > 0
      if not os.path.exists('MonteCarlo'):
         os.mkdir('MonteCarlo')
      for objid in self.id:
         objout = 'MonteCarlo/MCOutput_OBJ%d.txt' % objid
         if not os.path.exists(objout):
            with open(objout, 'wb') as log:
               print >> log,"# ITER  ZPHOT  CHI2  logAGE  EBMV  logSMASS  SFR"
      for niter in range(N):
         # create a new input catalog with perturbed photometry
         # Now perturb detected magnitudes; do NOTHING to upper limits
         newcat = '%s_MC.cat' % (root)
         with open('MonteCarlo/' + newcat, 'wb') as f:
            f.write(self.catHeader)
            for i in range(len(self.id)):
               newline = "%d  " % self.id[i]
               for j in range(1, len(self.catColumns))[::2]:
                  mag = getattr(self, self.catColumns[j])[i]
                  mag_err = getattr(self, self.catColumns[j+1])[i]
                  if (mag < 0) or (mag > 90):
                     # no photometry in this filter
                     pass
                  elif (mag_err < 0):
                     # upper limit in this filter
                     pass
                  else:
                     mag = np.random.normal(mag, mag_err)
                  newline += "%.3f  %.3f  " % (mag, mag_err)
               newline += "\n"
               f.write(newline)
         # update the parameter file with the perturbed catalog
         newparamfile = os.path.splitext(self.paramfile)[0]+"_MC.param" 
         with open(self.paramfile, 'rb') as f2:
            lines2 = f2.readlines()
         # Now replace input values with values for each iteration of MC
         for i2 in range(len(lines2)):
            if lines2[i2].startswith('CAT_IN'):
               lines2[i2] = lines2[i2].replace(self.inputcat, newcat)
            elif lines2[i2].startswith('CAT_OUT'):
               l2 = lines2[i2].split()
               l2[1] = "%s.out" % (os.path.splitext(l2[1])[0]+'_MC')
               lines2[i2] = ' '.join(l2) + '\n'
            elif lines2[i2].startswith('PARA_OUT'):
               l2 = lines2[i2].split()
               print l2
               if not os.path.exists('MonteCarlo/%s' % l2[1]):
                  os.system('cp %s MonteCarlo/%s' % (l2[1], l2[1]))
               # l2[1] = "%s" % l2[1]
               # lines2[i2] = ' '.join(l2)
            # elif lines2[i2].startswith('PDZ_OUT'):
            #    l2 = lines2[i2].split()
            #    l2[1] = "%s" % (l2[1] + '_run%d' % niter)
            #    lines2[i2] = ' '.join(l2) + '\n'
         # write new parameter file
         print "Write new parameter file %s..." % newparamfile
         os.chdir('MonteCarlo')
         with open(newparamfile, 'wb') as f2n:
            for i2 in range(len(lines2)):
               f2n.write(lines2[i2])
         # run LePhare fitting
         print "Now run LePhare for iteration %d..." % niter
         rcode = subprocess.call(["zphota", "-c", newparamfile])
         # read the best-fit parameters, store the results
         # I need the following properties:
         # - photo z
         # - stellar mass
         # - age
         # - E(B-V)
         # - SFR
         print "Collect the best-fit stellar pop values..."
         newrun = LePhare(newparamfile)
         newrun.readCatOut()
         for objid in self.id:
            # self.readObjSpec(objid, specdir='.')
            # bestProps = self.bestfitProps['GAL-1']  
            # Read the stellar pop. parameters from the output catalog
            ## restrict to the first galaxy component
            with open('MCOutput_OBJ%d.txt' % objid, 'ab') as log:
               j = np.arange(len(newrun.data['IDENT']))[newrun.data['IDENT']==objid][0]
               # the columns are ITER  PHOTZ  CHI2  AGE  EBMV  SMASS  SFR
               photz = newrun.data['Z_BEST'][j]
               chi2 = newrun.data['CHI_BEST'][j]
               age = newrun.data['AGE_BEST'][j]
               ebmv = newrun.data['EBV_BEST'][j]
               smass = newrun.data['MASS_BEST'][j]
               sfr = newrun.data['SFR_BEST'][j]
               outline = "%d  %.6f  %.6f  %.6e  %.6f  %.6e  %.6e" % (niter, photz, chi2, age, ebmv, smass, sfr)
               print >> log, outline
         # go back to the previous directory...
         print "Finished iteration %d." % niter
         os.chdir(MainDIR)
      print "Monte Carlo simulations all done!"

class PlotLePhare(LePhare):
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

