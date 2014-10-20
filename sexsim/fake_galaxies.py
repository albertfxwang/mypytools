#!/usr/bin/env python

import numpy as np
from pyraf import iraf
import pyfits
from pygoods import sextractor

"""
Generate attributes of fake galaxies randomly, following specified distributions
if necessary.
"""
# Axis ratio distribution parameters... from Ferguson et al. 2004
qpar = {"devauc": (0,0.3, 0.9), "expdisk": (1,0.05, 0.01)}
devauc = 0
disk = 1
devauc_correction = 0.
def gtype_str(gtype):
   if gtype == devauc:
     return "devauc"
   elif gtype == disk:
     return "expdisk"
   else:
     print "Galaxy type %d not recognized." % gtype
     return None

def axratio_devauc(ngal):
   qp = qpar['devauc']
   q = np.random.uniform(qp[1], qp[2], size=ngal) + qp[0]
   return q

def axratio_disk(ngal):
   qp = qpar['expdisk']
   q = np.random.normal(qp[1],qp[2], size=ngal)
   q = np.maximum(q, 0.01) 
   q = np.minimum(q, 1.)
   return q

# Tests to see if galaxy is on the image
def test_flag(x, y, flagimg):
   """Tests to see if galaxy is (mostly) on the image"""
   flg_max = flagimg[y-6:y+4, x-6:x+4].max()
   return flg_max

class FakeGalaxies(object):
   def __init__(self, realimages, flagimages, bands, magfile=None, 
                othercols_list=[], ngal=40, 
                diskfrac=0.5, rdist='loguniform', mag0=22.0, mag1=27.0, 
                logr0=-1.0, logr1=2.0, lognormal_beta=0.3, lognormal_mag0=25.0, 
                lognormal_peak=0.7, lognormal_sigma=0.5, flagmax=1):
      """
      A class that will generate random values for attributes of fake galaxies.
      The attributes are the same in all measurement images; in GALFIT sims we 
      do not care about the SEDs of individual sources.
      Arguments realimages, flagimages should be dictionaries.
      """
      self.ngal = ngal
      self.bands = bands
      self.realimages = realimages
      self.flagimages = flagimages
      self.diskfrac = diskfrac
      self.rdist = rdist
      self.magfile = magfile
      if self.magfile:
         self.c = sextractor(self.magfile)
         print "The first N columns of magfile have the input magnitudes of the sources,"
         print "where N = len(bands). Starting from column N+1, additional columns"
         print "will be stored."
         print "The first N columns are input magnitudes in the following order:"
         print " ".join(self.bands)
      self.mag0 = mag0
      self.mag1 = mag1
      self.logr0 = logr0
      self.logr1 = logr1
      self.lognormal_beta = lognormal_beta
      self.lognormal_peak = lognormal_peak
      self.lognormal_sigma = lognormal_sigma
      self.lognormal_mag0 = lognormal_mag0
      self.flagmax = flagmax
      hdr = pyfits.getheader(realimages[bands[0]])
      self.xmax = hdr['naxis1']
      self.ymax = hdr['naxis2']
      self.othercols_list = othercols_list
      # Initialize attributes...
      # If there is an input magnitude list, use the magnitudes for input.
      # The input magnitude file might also carry additional columns, which 
      # will be specified by the othercols parameters.
      self.spawn_galaxies()

   def spawn_galaxies(self):
      """
      Generate fake galaxies.
      """
      self.mag = {}
      self.othercols = {}
      self.othercolnames = []
      if self.magfile:
         # IMPORTANT
         # randomly draws self.ngal galaxies from magfile
         # So self.mag0, self.mag1 has NO effects
         igals = np.random.random_integers(0, len(self.c)-1, size=self.ngal)
         self.igals = igals
         # randomly generate detection-band magnitude
         # mag_b0 = np.random.uniform(self.mag0, self.mag1, size=self.ngal)
         for i in range(len(self.bands)):
            b = self.bands[i]
            self.mag[b] = getattr(self.c, '_%d' % (i+1)).take(igals)
         ## DO NOT SCALE THE MAGNITUDES IF ONE CARES ABOUT THE INPUT REDSHIFT!!
         # Now scale the magnitudes to match the detection-band magnitudes
         # dmag_b0 = mag_b0 - self.mag[self.bands[0]]
         # for b in self.bands:
         #    self.mag[b] = self.mag[b] + dmag_b0
         # Now collect other columns
         for j in range(len(self.othercols_list)):
            name = self.othercols_list[j]
            self.othercols[name] = getattr(self.c, '_%d' % (len(self.bands)+1+j)).take(igals)
      else:   
         for b in self.bands:
            self.mag[b] = np.random.uniform(self.mag0, self.mag1, size=self.ngal)
      self.othercolnames = self.othercols.keys()
      print "Other column names: ", self.othercolnames
      if self.rdist == 'lognormal':
         self.logre = self.rlognormal()
      else:
         # uniform distribution of input log(Re)
         self.logre = np.random.uniform(self.logr0, self.logr1, size=self.ngal)
      self.re = 10.**(self.logre)
      if np.__version__ >= '1.7.0':
         self.gtype = np.random.choice([devauc, disk], size=self.ngal, 
                                       p=[1.-self.diskfrac, self.diskfrac])
      else:
         r = np.random.uniform(0., 1., size=self.ngal)
         self.gtype = np.where(r>=self.diskfrac, devauc, disk)
      self.axis_ratio = self.get_axialratio()
      self.position_angle = np.random.uniform(0., 360., size=self.ngal)
      self.x, self.y = self.get_xy(flagmax=self.flagmax)
      self.detected = np.zeros(self.ngal)  # to store which galaxies are detected
      self.artfiles = {}

   def get_xy(self, edgebuffer=60, flagmax=1):
      """
      Randomly determines the input (x, y) positions of each artificial source.
      This method could be overwritten by child classes.
      """
      xarr = np.zeros(self.ngal)
      yarr = np.zeros(self.ngal)
      flag_img = pyfits.getdata(self.flagimages[self.bands[0]])
      for i in range(self.ngal):
         offimage = 1
         while offimage:
            x = np.random.uniform(edgebuffer/2., self.xmax-edgebuffer/2.)
            y = np.random.uniform(edgebuffer/2., self.ymax-edgebuffer/2.)
            if test_flag(x, y, flag_img) < flagmax:
               offimage = 0
               xarr[i] = x
               yarr[i] = y
      return xarr, yarr

   def rlognormal(self, b=None):
      """
      Randomly generates a radius drawn from lognormal distribution
      """
      # default to using the magnitude in the first band
      if b == None:
         b = self.bands[0]
      lumratio = 10**((self.mag[b]-self.lognormal_mag0)/-2.5)
      mu = self.lognormal_peak + self.lognormal_beta * np.log(lumratio)
      val = np.random.normal(mu, self.lognormal_sigma)
      return val

   def get_axialratio(self):
      """
      Gets a random axial ratio from a distribution of inclinations & true
      axial ratios
      """
      sini = np.random.rand();
      # Calculate intrinsic axis ratio q
      condlist = [self.gtype==devauc, self.gtype==disk]
      choicelist = [axratio_devauc(self.ngal), axratio_disk(self.ngal)]
      q = np.where(self.gtype==devauc, axratio_devauc(self.ngal), 
                   axratio_disk(self.ngal))
      cosi = np.sqrt(1-sini**2)
      ba = np.sqrt(sini*sini + (q*cosi)**2)
      return ba

   def update_mag_re(self):
      """
      Update the way one draws magnitudes and Re; to be implemented later.
      """
      raise NotImplementedError

   def makegals_multiband(self, flagimage, igalfile="", bands=None):
      """Makes the galaxy list files """
      
      print "in makegals"      
      # Write the galaxy parameters out to files
      if bands==None:
         bands = self.bands
      for b in bands:
         self.artfiles[b] = "glart_%s.list" % (b)  # input file for iraf.mkobject
         f = open(self.artfiles[b], "wb")
         for i in range(self.ngal):
            # Write lines for the mkobjects galaxy list
            if self.gtype[i] == devauc:
               f.write("%10.2f %10.2f %8.3f %12s %8.3f %6.2f %6.2f no " \
                  % (self.x[i], self.y[i], self.mag[b][i]+devauc_correction, \
                    gtype_str(self.gtype[i]), self.re[i], self.axis_ratio[i], \
                    self.position_angle[i]))
            else:
               f.write("%10.2f %10.2f %8.3f %12s %8.3f %6.2f %6.2f no " 
                  % (self.x[i], self.y[i], self.mag[b][i], gtype_str(self.gtype[i]), \
                     self.re[i], self.axis_ratio[i], self.position_angle[i]))
            f.write("\n")
       
         f.close()
      # Write out galaxies to igalfile if desired
      if len(igalfile) > 0:
         igfile = open(igalfile, 'a')  # the *.allgal file in the detection-band directory
         for i in range(self.ngal):
            outstring = "%10.2f %10.2f " % (self.x[i], self.y[i])
            outstring = outstring + "%d %8.3f %6.2f %6.2f" % (self.gtype[i], \
                        self.re[i], self.axis_ratio[i], self.position_angle[i])
            for b in self.bands:
               outstring = outstring + "%8.3f " % (self.mag[b][i])
            if len(self.othercols):
               keys = self.othercols.keys()
               for k in keys:
                  outstring = '%s ' % self.othercols[k][i]
            igfile.write("%s\n" % outstring)
         igfile.close()
      print "finish makegals"

class fake_galaxies(FakeGalaxies):
   pass

class FakeGalaxiesCustomXY(FakeGalaxies):
   def __init__(self, realimages, flagimages, bands, posfile, **kwargs):
      super(FakeGalaxiesCustomXY, self).__init__(realimages, flagimages, bands,
                                                 **kwargs)
      # now overwrite the input x, y
      self.x, self.y = self.get_custom_xy(posfile)

   def get_custom_xy(self, posfile):
      # posfile is a text file with allowed input positions (in pixel coord.)
      # this method will randomly drawn ngal of them as input positions
      # posfile needs to be in SExtractor format with columns x and y
      f = sextractor(posfile)
      assert hasattr(f, 'x')
      assert hasattr(f, 'y')
      nf = len(f)
      indices = np.random.choice(np.arange(self.ngal), size=self.ngal, 
                                 replace=False)
      xarr = self.x.take(indices)
      yarr = self.y.take(indices)
      return xarr, yarr

