#!/usr/bin/env python

"""
Simulations that insert artificial objects into images, and perform photometry
using SExtractor.
Last Updated: 14/06/19
Update Log:
- 140619: improved how the method collect_results handles the input columns;
          now read input columns (except the essential ones) from the parameter
          file (OTHERCOLNAMES).
"""

import numpy as np
import os, sys, glob, time, string
from sextutils import sextractor
import yaml
from pyraf import iraf
iraf.stsdas()
iraf.artdata()
iraf.images()
import pyfits
import fake_galaxies
from fake_galaxies import gtype_str
_hconvolve = True
try:
   from hconvolve import hconvolve, imhconvolve
except:
   print "No hconvolve... use scipy.signal.fftconvolve instead."
   _hconvolve = False
from scipy.signal import fftconvolve  # just in case

version= "0.1"
match_rad = 3.0   # positional match radius (in pixels) for artificial galaxies

class SExtractorSim(object):
   def __init__(self, parfile):
      assert parfile.endswith('.sim'), "Simulation input file must have extension .sim"
      c = yaml.load(open(parfile))
      self.c = c
      self.start = c['NSTART']
      self.nstop = c['NSTOP']
      self.n = c['NSTART']
      assert (c['SAVE'] and (c['NSTOP'] > 1))==False, "Error: You must use SAVE = no for NSTOP > 1"
      print "Simulation parameters from %s:" % (parfile)
      for k in c.keys():
         print "%-20s %s" % (k, c[k])
      self.detect_band = c['DETECTBAND']
      self.bands = c['BANDS']
      self.realimages = dict(zip(self.bands, c['REALIMAGE']))
      self.flagimages = dict(zip(self.bands, c['FLAGIMAGE']))
      self.rmsimages = dict(zip(self.bands, c['RMSIMAGE']))
      try:
         self.wht_type = c['WEIGHT_TYPE']
      except:
         self.wht_type = 'MAP_RMS'
      self.gains = dict(zip(self.bands, c['GAIN']))
      self.zeropoints = dict(zip(self.bands, c['MAGZPT']))
      self.psffiles = dict(zip(self.bands, c['PSFFILE']))
      self.detect_image = self.realimages[self.detect_band]
      self.fakeimages = {}
      self.segimages = {}
      self.catalogs = {}
      self.newcatalogs = {}
      self.noiselessimages = {}
      self.sextparfile = c['SEXTPARFILE']
      self.sexfile = c['SEXFILE']
      # Can take PSF-matching kernels... but require that they are specified
      # as dictionaries in the parameter file, to avoid confusion
      if 'PSFMATCH' not in c.keys():
         self.psfmatch = False
      else:
         self.psfmatch = c['PSFMATCH']
      if self.psfmatch:
         self.psfmatch_kernels = c['PSFMATCH_KERNELS']
      print "Artdata parameters:"
      iraf.lpar('artdata')
      self.root = os.path.splitext(parfile)[0]   # e.g. root = "run1m"
      hdr = pyfits.getheader(self.detect_image)
      self.xmax = hdr['naxis1']
      self.ymax = hdr['naxis2']
      # put in this try/except clause so that TPHOT sims will work...
      try:
         self.logrmin = np.log10(c['RMIN'] / c['SCALE'])  # in pixels
         self.logrmax = np.log10(c['RMAX'] / c['SCALE'])
         self.pixscale = c['SCALE']
      except:
         self.logrmin = np.log10(c['RMIN'] / c['SCALE'][0])
         self.logrmax = np.log10(c['RMAX'] / c['SCALE'][0])
         self.pixscale = c['SCALE'][0]
      self.rdistfunc = c['RADIUS_DISTRIBUTION']
      self.lognormal_mag0 = c['LOGNORMAL_MAG0']
      self.lognormal_peak = c['LOGNORMAL_PEAK']
      self.lognormal_sigma = c['LOGNORMAL_SIGMA']
      self.lognormal_beta = c['LOGNORMAL_BETA']
      self.diskfrac = c['DISKFRAC']
      self.ngal = c['NGALAXIES']
      self.maglow = c['MAGLOW']
      self.maghigh = c['MAGHIGH']
      self.flagmax = c['FLAGMAX']
      self.save = c['SAVE']
      if c.has_key('OTHERCOLNAMES'):
         self.othercolnames = c['OTHERCOLNAMES']
      else:
         self.othercolnames = []
      if c.has_key('MAGFILE'):
         self.magfile = c['MAGFILE']
      else:
         self.magfile = None
      for b in self.bands:
         broot = self.root+'_'+b
         if not os.path.exists(broot):  
            # create output directory if needed
            os.mkdir(broot)
      # Create list of input galaxies in the detect band directory
      igalfile = "%s_%s/%s.allgal" % (self.root, self.detect_band, self.root)
      igfile = open(igalfile, 'a')
      igfile.write("# %s \n" % (time.ctime()))
      igfile.write("# %s \n" % (version))
      igfile.close()
      self._finished_sextractor = False
      # Fake galaxies...
      fg_args = [self.realimages, self.flagimages, self.bands]
      fg_kwargs = dict(ngal=self.ngal, diskfrac=self.diskfrac, 
                  magfile=self.magfile, rdist=self.rdistfunc, mag0=self.maglow,
                  mag1=self.maghigh, logr0=self.logrmin, logr1=self.logrmax, 
                  lognormal_beta=self.lognormal_beta, 
                  lognormal_mag0=self.lognormal_mag0, 
                  lognormal_peak=self.lognormal_peak,
                  lognormal_sigma=self.lognormal_sigma, 
                  flagmax=self.flagmax,
                  othercols_list=self.othercolnames)
      self.fg = fake_galaxies.fake_galaxies(*fg_args, **fg_kwargs)
   
   def makenoiselessimage(self, band, galfile, magz, psffile, save=0, gain=1.0):
      """
      Creates a noiseless convolved image
      """
      assert os.path.exists(psffile), "PSF image %s does not exist." % psffile
      broot = self.root + '_%s' % band
      outfile = broot + '_sim.fits'
      outfile_nonoise = broot + '_sim_noiseless.fits'
      # print outfile,xmax,ymax,galfile,magz,gain
      iraf.unlearn('artdata')
      iraf.unlearn('mkobjects')
      iraf.artdata.dynrange=1.e5
      print "Running iraf.mkobject for %s..." % band
      iraf.mkobjects(outfile_nonoise, output="", title="", ncols=self.xmax, 
         nlines=self.ymax, header="", background=0.0, objects=galfile,
         xoffset=0., yoffset=0., star="gaussian", radius=0.1,
         beta=2.5, ar=1., pa=0., distance=1., exptime=1., 
         magzero=magz, gain=gain, rdnoise=0., poisson=0,
         seed=2, comments=1)
      print "Convolving with PSF..."
      if _hconvolve:
         imhconvolve(outfile_nonoise, psffile, outfile, overwrite=True)
      else:
         print "No hconvolve..."
         outimage = pyfits.getdata(outfile_nonoise)
         psfimage = pyfits.getdata(psffile)
         outimage2 = fftconvolve(outimage, psfimage, mode='same')
         h = pyfits.open(outfile, mode='update')
         h[0].data = outimage2
         h.flush()
         h.close()
      self.noiselessimages[band] = outfile
      os.remove(outfile_nonoise)

   def addsimulated(self, band, save=0):
      """
      Add the noiseless images of artificial galaxies to the real images.
      """
      # simulation = root+'_sim.fits'
      assert os.path.exists(self.noiselessimages[band]), \
         "Noiseless image with artificial galaxies not calculated."
      broot = self.root + '_%s' % band
      outimage = broot + '.fits'
      if os.path.exists(outimage):
         os.remove(outimage)
      noiseless_img = pyfits.getdata(self.noiselessimages[band])
      realimage_img = pyfits.getdata(self.realimages[band])
      hdr = pyfits.getheader(self.realimages[band])
      simulated_img = realimage_img + noiseless_img
      if self.psfmatch:
         if band == self.detect_band:
            pass
         else:
            assert band in self.psfmatch_kernels.keys(), "PSF-match kernel for %s does not exist." % band
            kernel = pyfits.getdata(self.psfmatch_kernels[band])
            print "Convolving with PSF-match kernel in %s..." % band
            # if _hconvolve:
            #    simulated_img = hconvolve(simulated_img, kernel)
            # else:
            # Have to use scipy for this one... otherwise there is some weird
            # artifacts from convolution
            simulated_img = fftconvolve(simulated_img, kernel, mode='same')
      pyfits.append(outimage, simulated_img, hdr)
      self.fakeimages[band] = outimage
      # iraf.imcalc(realimage+","+simulation,outimage,"im1+im2")
      if not save:
         os.remove(self.noiselessimages[band])

   def insert_fake_sources(self):
      """
      Generate attributes for artificial galaxies, and then insert artificial 
      galaxies into real images. Calls self.makenoiselessimage and 
      self.addsimulated.
      """
      ### Need to be updated for multi-band, input-SED case
      # if no specified list of fake galaxies
      # make artdata file of input fake galaxies
      if glob.glob('*.list'):
         os.system('rm *.list')
      # args = [self.realimages, self.flagimages, self.bands]
      # kwargs = dict(ngal=self.ngal, diskfrac=self.diskfrac, magfile=self.magfile,
      #               rdist=self.rdistfunc, mag0=self.maglow, mag1=self.maghigh,
      #               logr0=self.logrmin, logr1=self.logrmax, lognormal_beta=self.lognormal_beta,
      #               lognormal_mag0=self.lognormal_mag0, lognormal_peak=self.lognormal_peak,
      #               lognormal_sigma=self.lognormal_sigma, flagmax=self.flagmax,
      #               othercols=self.othercolnames)
      # self.fg = fake_galaxies.fake_galaxies(*args, **kwargs)
      self.fg.spawn_galaxies()
      # Will update later to use fake_galaxies_sexsim!
      for b in self.bands:
         self.fg.makegals_multiband(self.flagimages[b], bands=[b])  
         # write input file for mkobjects
         self.makenoiselessimage(b, self.fg.artfiles[b], self.zeropoints[b],
                                 self.psffiles[b])
         self.addsimulated(b)
      # designate a detection image
      self.fake_detect_image = self.fakeimages[self.detect_band]

   def run_sextractor(self, sex_exe='cex'):
      """
      Run SExtractor through all bands.
      """
      assert hasattr(self, 'fg'), "Please generate fake sources first."
      for b in self.bands:
         broot = self.root + '_' + b
         self.catalogs[b] = "%s_%d.cat" % (broot, self.n)
         self.segimages[b] = '%s_segnew.fits' % broot         
         self.newcatalogs[b] = "%s_%d.newcat" % (broot, self.n)
         n_args = "%s,%s -c %s" % (self.fake_detect_image, self.fakeimages[b], self.sexfile)
         n_args = n_args + " -CATALOG_NAME %s" % (self.newcatalogs[b])
         n_args = n_args + " -MAG_ZEROPOINT %9.4f" % (self.zeropoints[b])
         n_args = n_args + " -GAIN %12.4f" % (self.gains[b])
         n_args = n_args + " -FLAG_IMAGE %s" % (self.flagimages[b])
         n_args = n_args + " -PARAMETERS_NAME %s" % (self.sextparfile)
         n_args = n_args + " -CHECKIMAGE_TYPE SEGMENTATION"
         n_args = n_args + " -CHECKIMAGE_NAME %s" % (self.segimages[b])
         n_args = n_args + " -WEIGHT_TYPE %s,%s" % (self.wht_type, self.wht_type)
         if hasattr(self.c, 'WEIGHT_THRESH'):
            wt = self.c['WEIGHT_THRESH']
            n_args = n_args + " -WEIGHT_THRESH %f,%f" % (wt, wt)
         n_args = n_args + " -WEIGHT_IMAGE %s,%s" % (\
                           self.rmsimages[self.detect_band], self.rmsimages[b])
         print "%s %s" % (sex_exe, n_args)
         sys.stdout.flush()
         fpipe = os.popen("%s %s" % (sex_exe, n_args))
         fpipe.close()
         # Identify fake galaxies
         cnew = sextractor(self.newcatalogs[b])
         f = open(self.catalogs[b], 'wb')
         f.write(cnew._header)
         ncolumns = len(cnew._colnames)
         default_columns = ['X_IN', 'Y_IN', 'MAG_IN', 'RE_IN', 
                           'GTYPE_IN [devauc=%d, expdisk=%d' % (fake_galaxies.devauc, fake_galaxies.disk), 
                           'AXIS_RATIO_IN', 'PA_IN', 'ID_THISRUN']
         for j in range(len(default_columns)):
            ncolumns += 1
            f.write('# %d %s\n' % (ncolumns, default_columns[j]))
         # f.write('# %d X_IN\n' % (ncolumns+1))
         # f.write('# %d Y_IN\n' % (ncolumns+2))
         # f.write('# %d MAG_IN\n' % (ncolumns+3))
         # f.write('# %d RE_IN\n' % (ncolumns+4))
         # # f.write('# %d RE_IN_ARCSEC\n' % (ncolumns+5))
         # f.write('# %d GTYPE_IN [devauc=%d, expdisk=%d]\n' % ((ncolumns+6), fake_galaxies.devauc, fake_galaxies.disk))
         # f.write('# %d AXIS_RATIO_IN\n' % (ncolumns+7))
         # f.write('# %d PA_IN\n' % (ncolumns+8))
         # f.write('# %d ID_THISRUN\n' % (ncolumns+9))
         for j in range(len(self.fg.othercolnames)):
            ncolumns += 1
            f.write('# %d %s\n' % ((ncolumns), self.fg.othercolnames[j].upper()))
         if self.magfile:
            f.write('# %d IGAL\n' % (ncolumns+1))
         n_fake_gals = 0
         for i in range(self.ngal):
            dist = np.sqrt((self.fg.x[i]-cnew.x_image)**2 + (self.fg.y[i]-cnew.y_image)**2)
            if dist.min() > match_rad:
               continue
            j = np.argsort(dist)[0]  # index in cnew that matches this fake galaxy
            f.write(' '.join(cnew._colentries[j]))  # write the entry from cnew
            f.write(' %.2f %.2f ' % (self.fg.x[i], self.fg.y[i]))
            f.write(' %.2f %.2f ' % (self.fg.mag[b][i], self.fg.re[i]))
            f.write(' %4d  %.2f ' % (self.fg.gtype[i], self.fg.axis_ratio[i]))
            f.write(' %.2f  %4d ' % (self.fg.position_angle[i], i))
            for k in range(len(self.fg.othercolnames)):
               f.write(' %s ' % str(self.fg.othercols[self.fg.othercolnames[k]][i]))
            if self.magfile:
               f.write(' %s ' % str(self.fg.igals[i]))
            f.write('\n')
            n_fake_gals += 1
            self.fg.detected[i] = 1
         # Append non-detected galaxies in the end; substitute all values by -1
         for i in range(self.ngal):
            if self.fg.detected[i] == 0:
               f.write(' '.join([repr(-1)]*len(cnew._colnames)))
               f.write(' %.2f %.2f ' % (self.fg.x[i], self.fg.y[i]))
               f.write(' %.2f %.2f ' % (self.fg.mag[b][i], self.fg.re[i]))
               f.write(' %4d  %.2f ' % (self.fg.gtype[i], self.fg.axis_ratio[i]))
               f.write(' %.2f  %4d ' % (self.fg.position_angle[i], i))
               for k in range(len(self.fg.othercolnames)):
                  f.write(' %s ' % str(self.fg.othercols[self.fg.othercolnames[k]][i]))
               if self.magfile:
                  f.write(' %s ' % str(self.fg.igals[i]))
               f.write('\n')
         f.close()
         print "%d fake galaxies identified." % n_fake_gals
         os.system("mv %s %s/run%d.cat" % (self.catalogs[b], broot, self.n))
         os.system("mv %s %s/run%d.newcat" % (self.newcatalogs[b], broot, self.n)) 
         os.system("mv %s %s/glart%d.list" % (self.fg.artfiles[b], broot, self.n))                
         os.system("cp %s %s" % (self.psffiles[b], broot))

      self._finished_sextractor = True
   

   def cleanup(self):
      # if not self.save:
      # os.system('rm obj*.fits')
      print "removing %s*.fits" % self.root
      os.system('rm %s*.fits' % self.root)

   def collect_results(self, output_name, colname_file=None, careful_match=False):
      # allruns = {}
      bandcat = {}
      # first collect all SExtractor runs in the detection band
      os.chdir('%s_%s' % (self.root, self.detect_band))
      allruns = glob.glob('run*[0-9].cat')
      os.chdir('../')
      # for b in self.bands:
      #    os.chdir('%s_%s' % (self.root, b))
      #    allruns[b] = glob.glob('run*[0-9].cat')
      #    os.chdir('../')
      # then check if certain runs are missing in some filters
      # if yes, remove the catalog from the list for detection band
      print "checking if a run exists in all filters..."
      # for r in allruns[self.detect_band]:
      for r in allruns:
         for b in self.bands:
            if b == self.detect_band:
               pass
            else:
               # if r not in allruns[b]:
               if not os.path.exists('%s_%s/%s' % (self.root, b, r)):
                  print r
                  allruns.remove(r)
                  # allruns[b].remove(r)
      # Now merge catalogs for each filter, including all columns from the 
      # SExtractor runs. Each filter has a catalog named [root]_[band].cat.
      for b in self.bands:
         print "Merging catalog for %s..." % b
         broot = '%s_%s' % (self.root, b)
         os.chdir(broot)
         bandcat[b] = '%s.cat' % (broot)
         f = open(bandcat[b], 'wb')
         _wrote_header = 0 
         for r in allruns:
         # for r in allruns[self.detect_band]:
            sys.stdout.write("Processing %s_%s/%s... \r" % (self.root, b, r)) 
            sys.stdout.flush()
            # c = sextractor(r)
            # if _wrote_header == 0:
            #    colnames = c._colnames
            # ncols = len(colnames)
            # objstr = ""
            # if _wrote_header == 0:
            #    f.write(c._header)
            #    _wrote_header = 1
            # for i in range(len(c)):
            #    objstr += ' '.join(map(lambda x: str(c.__getattribute__(x)), colnames))
            #    objstr += '\n'
            # f.write(objstr)
            fc = open(r)
            lines = fc.readlines()
            if _wrote_header == 0:
               # f.write(c._header)
               for i in range(len(lines)):
                  if lines[i].startswith('# '):
                     f.write(lines[i])
               _wrote_header = 1
            for i in range(len(lines)):
               if not lines[i].startswith('#'):
                  f.write(lines[i])
         f.close()
         os.chdir('../')
      # Finally, merge catalogs from each filter, only including columns
      # specified by colname_file.
      if colname_file == None:
         colname_file = '/Users/khuang/Dropbox/codes/mypytools/sexsim/output_colnames.yml'
      cname = yaml.load(open(colname_file))
      # Build headers
      columns = {}
      colnames = []
      for b in self.bands:
         broot = '%s_%s' % (self.root, b)
         cb = sextractor('%s/%s' % (broot, bandcat[b]))
         if b == self.detect_band:
            cb0 = cb
         for cn in cname['SEXTRACTOR_COLUMNS']:
            name = '%s_%s' % (b.upper(), cn.upper())
            columns[name] = getattr(cb, cn.lower())
            colnames += [name]
         for cn in cname['INPUT_BAND_COLUMNS']:
            name = '%s_%s' % (b.upper(), cn.upper())
            columns[name] = getattr(cb, cn.lower())
            colnames += [name]
         # record a DETECT column
         if b == self.detect_band:
            colnames += ['DETECT']
            columns['DETECT'] = np.where(getattr(cb, 'x_image')>0, 1, 0)
      # following are input columns that do not depend on filters
      # revert back to detection band

      for cn in cname['INPUT_ALL_COLUMNS'] + self.othercolnames:
         name = cn.upper()
         columns[name] = getattr(cb0, cn.lower())
         colnames += [name]
      # write headers
      print colnames
      f = open(output_name, 'wb')
      for i in range(len(colnames)):
         f.write('# %d %s\n' % (i+1, colnames[i]))
      # write results
      for i in range(len(cb)):
         row = map(lambda x: str(columns[x][i]), colnames)
         f.write(' '.join(row))
         f.write('\n')
      self.columns = columns
      f.close()

def fix_redshift_col(c, colname='redshift', orderCorrect=103, copy=False):
   # Fix the column order 
   colOrder = c._d[colname]
   if copy:
      os.system('cp %s %s' % (c._fname, c._fname+'.copy'))
   if colOrder != orderCorrect:
      print "Fix %s..." % c._fname
      for name in c._colnames:
         # identify which column occupies the desired column order
         if c._d[name] == orderCorrect:
            otherColname = name
            break
      # swap their order
      c._d[colname] = orderCorrect
      c._d[otherColname] = colOrder
      # now re-write the catalog
      f = open(c._fname, 'wb')
      # first, write headers
      newcolnames = []
      for i in range(len(c._colnames)):
         for cn in c._colnames:
            if c._d[cn] == i+1:
               f.write('# %d %s\n' % ((i+1), cn.upper()))
               newcolnames.append(cn)
      # Now write values
      objstr = ""
      for j in range(len(c)):
         objstr += ' '.join(map(lambda x: str(c.__getattribute__(x)[j]), newcolnames))
         objstr += '\n'
      f.write(objstr)
      f.close()

def check_all(allruns, checked=None):
   # To check the column orders in all SExtractor sim output
   if checked==None:
      checked = [False] * len(allruns)
   try:
      for i in range(len(allruns)):
         n = os.path.split(allruns[i])[-1]
         n = os.path.splitext(n)[0][3:]
         n = int(n)
         if not checked[n]:
            c = sextractor(allruns[i])
            fix_redshift_col(c)
            checked[n] = True
   except:
      print "Something weird happened... Last checked was %s" % allruns[i]
   return checked

class sextractor_sim(SExtractorSim):
   pass

def calc_maglim(catalog, band, root, fluxerr_factor, apercol='iso', sigma=3, dsigma=0.2, func='median'):
   """
   Calculate the N-sigma magnitude limit from the merged simulation catalog.
   """
   c = sextractor(catalog)
   fluxcol = '%s_flux_%s' % (band.lower(), apercol.lower())
   fluxerrcol = '%s_fluxerr_%s' % (band.lower(), apercol.lower())
   flux = getattr(c, fluxcol)
   fluxerr = fluxerr_factor * getattr(c, fluxerrcol)
   sn = flux / fluxerr
   # choose only within a narrow range of S/N
   sn_calc = (sn >= (sigma-dsigma)) & (sn <= (sigma+dsigma))
   # ensure that the fake source is detected
   x_image = getattr(c, '%s_x_image' % band.lower())
   y_image = getattr(c, '%s_y_image' % band.lower())
   detect = (x_image > 0) & (y_image > 0)
   calc = np.logical_and(sn_calc, detect)
   assert np.sum(calc) >= 20, "Fewer than 20 (only %d) objects in the calculation..." % np.sum(calc)
   flux_sn = flux[calc]
   # Now read the zeropoint from the simulation parameter file
   zp = -1.0
   par = yaml.load(open('%s.sim' % root))
   for i in range(len(par['BANDS'])):
      if par['BANDS'][i] == band.lower():
         zp = par['MAGZPT'][i]
   assert zp > 0, "Unable to read zeropoint from the parameter file..."
   if func == 'median':
      maglim_sn = np.median(zp - 2.5 * np.log10(flux_sn))
   else:
      maglim_sn = np.average(zp - 2.5 * np.log10(flux_sn))
   return maglim_sn

      