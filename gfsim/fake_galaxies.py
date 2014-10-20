#!/usr/bin/env python

import numpy as np
from pyraf import iraf
import pyfits

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

class fake_galaxies_gfsim(object):
   def __init__(self, realimages, flagimages, bands, ngal=40, diskfrac=0.5, 
                rdist='loguniform', mag0=22.0, mag1=27.0, logr0=-1.0, 
                logr1=2.0, lognormal_beta=0.3, lognormal_mag0=25.0, 
                lognormal_peak=0.7, lognormal_sigma=0.5,
                edgebuffer=60, flagmax=1):
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
      flag_img = pyfits.getdata(flagimages[bands[0]])
      self.diskfrac = diskfrac
      self.rdist = rdist
      self.mag0 = mag0
      self.mag1 = mag1
      self.logr0 = logr0
      self.logr1 = logr1
      self.lognormal_beta = lognormal_beta
      self.lognormal_peak = lognormal_peak
      self.lognormal_sigma = lognormal_sigma
      self.lognormal_mag0 = lognormal_mag0
      hdr = pyfits.getheader(realimages[bands[0]])
      xmax = hdr['naxis1']
      ymax = hdr['naxis2']
      # Initialize attributes...
      self.mag = np.random.uniform(mag0, mag1, size=ngal)
      if rdist == 'lognormal':
         self.logre = self.rlognormal()
      else:
         self.logre = np.random.uniform(logr0, logr1, size=ngal)
      self.re = 10.**(self.logre)
      self.gtype = np.random.choice([devauc, disk], size=ngal, 
                                    p=[1.-diskfrac, diskfrac])
      self.axis_ratio = self.get_axialratio()
      self.position_angle = np.random.uniform(0., 360., size=ngal)
      self.x = np.zeros(ngal)
      self.y = np.zeros(ngal)
      for i in range(ngal):
         offimage = 1
         while offimage:
            x = np.random.uniform(edgebuffer/2., xmax-edgebuffer/2.)
            y = np.random.uniform(edgebuffer/2., ymax-edgebuffer/2.)
            if test_flag(x, y, flag_img) < flagmax:
               offimage = 0
               self.x[i] = x
               self.y[i] = y
      self.artfiles = {}

   def rlognormal(self):
      """
      Randomly generates a radius drawn from lognormal distribution
      """
      lumratio = 10**((self.mag-self.lognormal_mag0)/-2.5)
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
      q = np.where(self.gtype==devauc, axratio_devauc(self.ngal), axratio_disk(self.ngal))
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
                  % (self.x[i], self.y[i], self.mag[i]+devauc_correction, \
                    gtype_str(self.gtype[i]), self.re[i], self.axis_ratio[i], \
                    self.position_angle[i]))
            else:
               f.write("%10.2f %10.2f %8.3f %12s %8.3f %6.2f %6.2f no " 
                  % (self.x[i], self.y[i], self.mag[i], gtype_str(self.gtype[i]), \
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
               outstring = outstring + "%8.3f " % (self.mag[i])
            igfile.write("%s\n" % outstring)
         igfile.close()
      print "finish makegals"

# class fake_galaxies_sexsim(fake_galaxies_gfsim):
#    """
#    Mostly the same as fake_galaxies_gfsim, only now need to take care of each
#    galaxy's SED.
#    """
#    # Not implemented yet... will update later
#    pass