#!/usr/bin/env python

"""
Simulations that insert artificial objects into images, and perform photometry
using SExtractor.
"""

import numpy as np
import os, sys, glob, time, string
from pygoods import *
import yaml
from pyraf import iraf
iraf.stsdas()
iraf.artdata()
iraf.images()
import fake_galaxies
from fake_galaxies import gtype_str
from hconvolve import imhconvolve

version= "0.1"
match_rad = 3.0   # positional match radius (in pixels) for artificial galaxies

class sextractor_sim(object):
   def __init__(self, parfile):
      assert parfile.endswith('.sim'), "Simulation input file must have extension .sim"
      try:
         c = yaml.load(open(parfile))
      except:
         c = parseconfig(parfile)
      self.c = c
      self.niter = c['NITER']
      self.n = c['NSTART']
      print c['SAVE']
      assert (c['SAVE'] and (c['NITER'] > 1))==False, "Error: You must use SAVE = no for NITER > 1"
      print "Simulation parameters from %s:" % (parfile)
      for k in c.keys():
         print "%-20s %s" % (k, c[k])
      self.detect_band = c['DETECTBAND']
      self.bands = c['BANDS']
      self.realimages = dict(zip(self.bands, c['REALIMAGE']))
      self.flagimages = dict(zip(self.bands, c['FLAGIMAGE']))
      self.rmsimages = dict(zip(self.bands, c['RMSIMAGE']))
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
      print "Artdata parameters:"
      iraf.lpar('artdata')
      self.root = os.path.splitext(parfile)[0]   # e.g. root = "run1m"
      hdr = pyfits.getheader(self.detect_image)
      self.xmax = hdr['naxis1']
      self.ymax = hdr['naxis2']
      self.logrmin = np.log10(c['RMIN'] / c['SCALE'])  # in pixels
      self.logrmax = np.log10(c['RMAX'] / c['SCALE'])
      self.pixscale = c['SCALE']
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
         
      # otherpars=[]
      # # sparfile = self.root+'.sextpar'
      # if (len(magfile) > 0): 
      #    print "Retrieving magnitudes from %s" % (magfile)
      #    (inputmags,input_otherpars,notherpars)=getdata(magfile,bands)
      #    # Update VECTOR_ASSOC in the SExtractor parameter file
      #    sfi = open(sextparfile,'r')
      #    sfo = open(sparfile,'w')
      #    lines = sfi.readlines()
      #    for l in lines:
      #      if string.find(l,'VECTOR_ASSOC') == -1:
      #        sfo.write(l)
      #    sfo.write('VECTOR_ASSOC(%d)\n' % (notherpars+7))
      #    sfi.close()
      #    sfo.close()
      # else:
      # shutil.copyfile(sextparfile,sparfile)
   
   def makenoiselessimage(self, band, galfile, magz, psffile, save=0, gain=1.0):
      """
      Creates a noiseless convolved image
      """
      assert os.path.exists(psffile), "PSF image %s does not exist." % psffile
      broot = self.root + '_%s' % band
      outfile = broot + '_sim.fits'
      # print outfile,xmax,ymax,galfile,magz,gain
      iraf.unlearn('artdata')
      iraf.unlearn('mkobjects')
      iraf.artdata.dynrange=1.e5
      print "Running iraf.mkobject for %s..." % band
      iraf.mkobjects(outfile, output="", title="", ncols=self.xmax, 
         nlines=self.ymax, header="", background=0.0, objects=galfile,
         xoffset=0., yoffset=0., star="gaussian", radius=0.1,
         beta=2.5, ar=1., pa=0., distance=1., exptime=1., 
         magzero=magz, gain=gain, rdnoise=0., poisson=0,
         seed=2, comments=1)
      imhconvolve(outfile, psffile, outfile, overwrite=True)
      self.noiselessimages[band] = outfile

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
      pyfits.append(outimage, simulated_img, hdr)
      self.fakeimages[band] = outimage
      # iraf.imcalc(realimage+","+simulation,outimage,"im1+im2")
      if not save:
         os.remove(self.noiselessimages[band])

   def insert_fake_sources(self):
      """
      Generate attributes for artificial galaxies, and then insert artificial 
      galaxies into real images.
      """
      # if no specified list of fake galaxies
      # make artdata file of input fake galaxies
      if glob.glob('*.list'):
         os.system('rm *.list')
      args = [self.realimages, self.flagimages, self.bands]
      kwargs = dict(ngal=self.ngal, diskfrac=self.diskfrac, 
                    rdist=self.rdistfunc, mag0=self.maglow, mag1=self.maghigh,
                    logr0=self.logrmin, logr1=self.logrmax, lognormal_beta=self.lognormal_beta,
                    lognormal_mag0=self.lognormal_mag0, lognormal_peak=self.lognormal_peak,
                    lognormal_sigma=self.lognormal_sigma, flagmax=self.flagmax)
      self.fg = fake_galaxies.fake_galaxies_gfsim(*args, **kwargs)
      # Will update later to use fake_galaxies_sexsim!
      for b in self.bands:
         self.fg.makegals_multiband(self.flagimages[b], bands=[b])  
         # write input file for mkobjects
         self.makenoiselessimage(b, self.fg.artfiles[b], self.zeropoints[b],
                                 self.psffiles[b])
         self.addsimulated(b)

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
         n_args = "%s,%s -c %s" % (self.fakeimages[self.detect_band], self.fakeimages[b], self.sexfile)
         n_args = n_args + " -CATALOG_NAME %s" % (self.newcatalogs[b])
         n_args = n_args + " -MAG_ZEROPOINT %9.4f" % (self.zeropoints[b])
         n_args = n_args + " -GAIN %12.4f" % (self.gains[b])
         n_args = n_args + " -FLAG_IMAGE %s" % (self.flagimages[b])
         n_args = n_args + " -PARAMETERS_NAME detectnew.param"
         n_args = n_args + " -CHECKIMAGE_TYPE SEGMENTATION"
         n_args = n_args + " -CHECKIMAGE_NAME %s" % (self.segimages[b])
         n_args = n_args + " -WEIGHT_TYPE MAP_RMS,MAP_RMS"
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
         f.write('# %d X_IN\n' % (ncolumns+1))
         f.write('# %d Y_IN\n' % (ncolumns+2))
         f.write('# %d MAG_IN\n' % (ncolumns+3))
         f.write('# %d RE_IN\n' % (ncolumns+4))
         f.write('# %d GTYPE_IN [devauc=%d, expdisk=%d]\n' % ((ncolumns+5), fake_galaxies.devauc, fake_galaxies.disk))
         f.write('# %d AXIS_RATIO_IN\n' % (ncolumns+6))
         f.write('# %d PA_IN\n' % (ncolumns+7))
         f.write('# %d ID_THISRUN\n' % (ncolumns+8))
         n_fake_gals = 0
         jlist = []
         for i in range(self.ngal):
            dist = np.sqrt((self.fg.x[i]-cnew.x_image)**2 + (self.fg.y[i]-cnew.y_image)**2)
            if dist.min() > match_rad:
               continue
            j = np.argsort(dist)[0]  # index in cnew that matches this fake galaxy
            if j not in jlist:
               # to avoid double-counting objects in cnew as matches to fake galaxies
               jlist += [j]
               f.write(' '.join(cnew._colentries[j]))  # write the entry from cnew
               f.write(' %.2f %.2f ' % (self.fg.x[i], self.fg.y[i]))
               f.write(' %.2f %.2f ' % (self.fg.mag[i], self.fg.re[i]))
               f.write(' %4d  %.2f ' % (self.fg.gtype[i], self.fg.axis_ratio[i]))
               f.write(' %.2f  %4d ' % (self.fg.position_angle[i], i))
               f.write('\n')
               n_fake_gals += 1
         f.close()
         print "%d fake galaxies identified." % n_fake_gals
         os.system("mv %s %s/run%d.cat" % (self.catalogs[b], broot, self.n))
         os.system("mv %s %s/run%d.newcat" % (self.newcatalogs[b], broot, self.n)) 
         os.system("mv %s %s/glart%d.list" % (self.fg.artfiles[b], broot, self.n))                
         os.system("cp %s %s" % (self.psffiles[b], broot))

      self._finished_sextractor = True
   
   def cleanup(self):
      if not self.save:
         # os.system('rm obj*.fits')
         print "removing %s*.fits" % self.root
         os.system('rm %s*.fits' % self.root)

   def collect_results(self):
      raise NotImplementedError
      # for b in bands:
      #    broot = root+'_'+b
      #    fakeimage = broot+".fits"
      #    catalog = "%s_%d.cat" % (broot,n)
      #    print catalog
      #    try:
      #       cat = sextractor(catalog)
      #       ncat = len(cat)
      #    except:
      #       ncat = 0  # no SExtractor detection of fake objects
      #       continue

   
