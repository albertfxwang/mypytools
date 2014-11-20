#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import photomUtils as pu
from scipy.integrate import cumtrapz
from pygoods import sextractor
import copy, os, sys, subprocess
from stats import distributions as dist


# a list of model tau's (for exponentially declining SFH)
tau_array = [0.1, 0.3, 1.0, 2.0, 3.0, 5.0, 10.0, 15.0, 30.0]
numTau = len(tau_array)

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

   def readObjOutput(self, objid):
      try:
         j = np.arange(len(self.data['IDENT']))[self.data['IDENT']==objid][0]
      except:
         self.readCatOut()
         j = np.arange(len(self.data['IDENT']))[self.data['IDENT']==objid][0]
      output_dict = {}
      for k in self.data.keys():
         output_dict[k] = self.data[k][j]
      return output_dict

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
            self.wave = np.array(self.wave)
            self.flux = np.array(self.flux)
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
               self.wave2 = np.array(self.wave2)
               self.flux2 = np.array(self.flux2)
               j = j + j2
            if self.bestfitProps['QSO']['Nline'] > 0:
               # read QSO component
               self.wave3 = []
               self.flux3 = []
               print "Reading SED for QSO model..."
               j3 = 0
               while j3 < self.bestfitProps['QSO']['Nline']:
                  if len(lines[i+j+j3].split()) == 2:
                     # print i+j+j2
                     l = lines[i+j+j3].split()
                     self.wave3 += [float(l[0])]
                     self.flux3 += [pu.ABmag2uJy(float(l[1]))]
                     j3 = j3 + 1
               self.wave3 = np.array(self.wave3)
               self.flux3 = np.array(self.flux3)
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
                     if l3.lower() not in self.catColumns:
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
               print >> log,"# ITER  ZPHOT  CHI2  AGE  EBMV  logSMASS  SFR"
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
               ebmv = newrun.data['EBV_BEST'][j]
               smass = newrun.data['MASS_BEST'][j]
               sfr = newrun.data['SFR_BEST'][j]
               age = newrun.data['AGE_BEST'][j]
               # age = smass / sfr  ## A KLUDGE THAT ONLY WORKS FOR CSF!!
               outline = "%d  %.6f  %.6f  %.6e  %.6f  %.6e  %.6e" % (niter, photz, chi2, age, ebmv, smass, sfr)
               print >> log, outline
         # go back to the previous directory...
         print "Finished iteration %d." % niter
         os.chdir(MainDIR)
      print "Monte Carlo simulations all done!"

   def readMCResults(self, objid):
      assert os.path.exists('MonteCarlo'), "The directory MonteCarlo does not exist."
      c = sextractor('MonteCarlo/MCOutput_OBJ%d.txt' % objid)
      self.MCresults = {}
      self.MCresults['photz'] = c._2
      self.MCresults['chi2'] = c._3
      self.MCresults['log_age'] = np.log10(c._4)
      self.MCresults['ebmv'] = c._5
      self.MCresults['log_mass'] = c._6
      self.MCresults['log_sfr'] = c._7
      self.MCresults['log_ssfr'] = c._7 - c._6
      # This only works with constant SFH!!
      # self.MCresults['log_age'] = self.MCresults['log_mass'] - self.MCresults['log_sfr']

   def calc_conf_intervals(self, objid, xgrid=200, p=0.68, print_it=False):
      """
      Calculate the confidence intervals for the following properties:
      - photo-z
      - log10(AGE)
      - log10(MASS)
      - log10(SFR)
      - log10(sSFR)
      """
      try:
         objout_dict = self.readObjOutput(objid)
      except:
         self.getCatOut()
         self.readCatOut()
         objout_dict = self.readObjOutput(objid)  
         # THE LePhare output for this object
      self.readMCResults(objid)
      confInt = {}
      # confInt stores the best-fit values and distances to the upper and lower
      # bounds of the confidence interval
      # For example, if the confidence interval is between 5.5 and 6.5, with
      # the best-fit value being 5.9, then 
      # confInt['photz'] = [5.9, 0.6, 0.4]

      # phot-z
      photz_dist = dist.Distribution1D(self.MCresults['photz'])
      photz_best = objout_dict['Z_BEST']
      if print_it: print "Confidence interval for phot-z:"
      photz_lo, photz_hi = photz_dist.conf_interval(xgrid=xgrid, p=p, 
                                                    x0=photz_best)
      confInt['photz'] = [photz_best, photz_hi-photz_best, photz_best-photz_lo]
      if print_it: print ""

      # log10(MASS)  [M_solar]
      log_mass_best = objout_dict['MASS_BEST']
      log_mass_dist = dist.Distribution1D(self.MCresults['log_mass'])
      if print_it: print "Confidence interval for log10(MASS):"
      log_mass_lo, log_mass_hi = log_mass_dist.conf_interval(xgrid=xgrid,
                                 p=p, x0=log_mass_best, print_it=print_it)
      confInt['log_mass'] = [log_mass_best, log_mass_hi-log_mass_best, log_mass_best-log_mass_lo]
      if print_it: print ""

      # log10(SFR)   [M_solar/yr]
      log_sfr_best = objout_dict['SFR_BEST']
      log_sfr_dist = dist.Distribution1D(self.MCresults['log_sfr'])
      if print_it: print "Confidence interval for log10(SFR):"
      log_sfr_lo, log_sfr_hi = log_sfr_dist.conf_interval(xgrid=xgrid, p=p,
                                 x0=log_sfr_best, print_it=print_it)
      confInt['log_sfr'] = [log_sfr_best, log_sfr_hi-log_sfr_best, log_sfr_best-log_sfr_lo]
      if print_it: print ""

      # log10(sSFR)   [yr^-1]
      log_ssfr_best = log_sfr_best - log_mass_best
      log_ssfr_dist = dist.Distribution1D(self.MCresults['log_sfr'] - self.MCresults['log_mass'])
      if print_it: print "Confidence interval for log10(sSFR):"
      log_ssfr_lo, log_ssfr_hi = log_ssfr_dist.conf_interval(xgrid=xgrid,
                                 p=p, x0=log_ssfr_best, print_it=print_it)
      confInt['log_ssfr'] = [log_ssfr_best, log_ssfr_hi-log_ssfr_best, log_ssfr_best-log_ssfr_lo]
      if print_it: print ""

      # log10(AGE)  [yr]
      ## AGAIN: THIS ONLY WORKS FOR CONSTANT SFH!!!
      log_age_best = np.log10(objout_dict['AGE_BEST'])
      # log_age_best = objout_dict['MASS_BEST'] - objout_dict['SFR_BEST']
      log_age_dist = dist.Distribution1D(self.MCresults['log_age'])
      if print_it: print "Confidence interval for log10(AGE):"
      log_age_lo, log_age_hi = log_age_dist.conf_interval(xgrid=xgrid, p=p, 
                                 x0=log_age_best, print_it=print_it)
      confInt['log_age'] = [log_age_best, log_age_hi-log_age_best, log_age_best-log_age_lo]
      if print_it: print ""

      return confInt

