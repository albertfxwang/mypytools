#!/usr/bin/env python

import numpy as np
import pyfits
import pywcs
import os
from hconvolve import hconvolve, imhconvolve
from pyraf import iraf  
iraf.artdata()

"""
Classes that represent FITS images. Also contains special classes for images
from different instruments.
"""
### Functions that cannot be written as methods.
### ......

def compare_crvals(fitsimage1, fitsimage2, tol=0.01):
   """
   Compare the CRVALs of two images within a tolerance (given in arcsec).
   """
   tol_degree = tol / 3600.
   same_crval1 = (abs(fitsimage1.hdr['crval1'] - fitsimage2.hdr['crval1']) <= tol_degree)
   same_crval2 = (abs(fitsimage1.hdr['crval2'] - fitsimage2.hdr['crval2']) <= tol_degree)
   return (same_crval1 & same_crval2)

def match_images_multires(image1, image2):
   """
   Find the WCS footprint of image1, and make a cutout from image2 that matches
   the footprint of image1.
   Make sure both images have the same CRVALs.
   """
   # First shift their CRPIXs
   fitsimage1 = fitsimage(image1)
   fitsimage2 = fitsimage(image2)
   ready1 = fitsimage1.check_crpix_half_integer()
   ready2 = fitsimage2.check_crpix_half_integer()
   if (ready1==False):
      # shift CRPIXs of image1 to half-integers, and record the new CRVALs
      # and then shift CRVALs of image2 as well.
      crpix1_1 = fitsimage1.hdr['crpix1']
      crpix1_2 = fitsimage1.hdr['crpix2']
      new_crpix1_1 = int(crpix1_1) + 0.5
      new_crpix1_2 = int(crpix1_2) + 0.5
      new_crval1, new_crval2 = fitsimage1.shift_crpix(new_crpix1_1, new_crpix1_2)
      fitsimage1 = fitsimage(image1)
      # Now shift CRPIXs of image2
      fitsimage2.shift_crval(new_crval1, new_crval2)
      fitsimage2 = fitsimage(image2)
      new_crpix2_1 = fitsimage2.hdr['crpix1']
      new_crpix2_2 = fitsimage2.hdr['crpix2']
      assert fitsimages.check_crpix_half_integer(), \
         "CRPIX for %s: [%.2f, %.2f] still not half-integers..." % (image2, new_crpix2_1, new_crpix2_2)
      print "Both images now have half-integer CRPIXs:"
      print "%s has CRPIXs = [%.2f, %.2f]" % (image1, fitsimage1.hdr['crpix1'], fitsimage1.hdr['crpix2'])
      print "%s has CRPIXs = [%.2f, %.2f]" % (image2, fitsimage2.hdr['crpix1'], fitsimage2.hdr['crpix2'])
   elif (ready2==False) | (compare_crvals(fitsimage1, fitsimage2)==False):
      # Now update the WCS of image2 only
      crval1_1 = fitsimage1.hdr['crval1']
      crval1_2 = fitsimage1.hdr['crval2']
      new_crpix2_1, new_crpix2_2 = fitsimage2.shift_crval(crval1_1, crval1_2)
      fitsimage2 = fitsimage(image2)
      assert fitsimage2.check_crpix_half_integer(), \
         "CRPIX for %s: [%.2f, %.2f] still not half-integers..." % (image2, new_crpix2_1, new_crpix2_2)
   # Make a cutout from image2 that matches the WCS footprint of image1
   ncols1 = fitsimage1.hdr['naxis1']
   nlines1 = fitsimage1.hdr['naxis2']
   ra_min, dec_min = fitsimage1.wcs.wcs_pix2sky([[1,1]], 1)[0]
   ra_max, dec_max = fitsimage1.wcs.wcs_pix2sky([[ncols1, nlines1]], 1)[0]
   x2_min, y2_min = np.around(fitsimage2.wcs.wcs_sky2pix([[ra_min, dec_min]], 1)[0])
   x2_max, y2_max = np.around(fitsimage2.wcs.wcs_sky2pix([[ra_max, dec_max]], 1)[0])
   print "pixel range for image2 to match image1: [%d:%d,%d:%d]" % (x2_min, x2_max, y2_min, y2_max)
   return x2_min, x2_max, y2_min, y2_max

class fitsimage(object):
   """
   A base class for images. Will define info that applies to all images here.
   Instrument-specific attributes will be defined in subclasses.
   For now, assume that the images are 2D (have I ever used a 3D image?)
   """
   def __init__(self, filename, make_wcs=True):
      self.filename = filename
      hdr = pyfits.getheader(filename)
      self.hdr = hdr
      self.data = pyfits.getdata(filename)
      # self.psffile = psffile  # store the file name for its PSF
      self._has_wcs = make_wcs
      if self._has_wcs:
         self.wcs = pywcs.WCS(hdr)
      else:
         self.wcs = None

   def export(self, name=None):
      """
      Export current header and data to a file.
      """
      assert name, "Invalid file name to export to."
      if os.path.exists(name):
         h = pyfits.open(name, mode='update')
         h[0].header = self.hdr
         h[0].data = self.data
         h.flush()
         h.close()
      else:
         pyfits.append(name, self.data, self.hdr)

   def check_crpix_half_integer(self):
      """
      Check if CRPIXs are half integers.
      """
      x1 = (abs((self.hdr['CRPIX1'] - int(self.hdr['CRPIX1']))) == 0.5)
      if not x1:
         print "CRPIX1 %.2f is not a half integer." % (self.hdr['CRPIX1'])
      x2 = (abs((self.hdr['CRPIX1'] - int(self.hdr['CRPIX1']))) == 0.5)
      if not x2:
         print "CRPIX2 %.2f is not a half integer." % (self.hdr['CRPIX2'])
      return (x1 and x2)

   def shift_crpix(self, new_crpix1, new_crpix2, export=True):
      """
      Shift WCS to use new values for CRPIX1 and CRPIX2, and write the updated
      values to the header.
      """
      assert self._has_wcs, "self.wcs is not yet defined!"
      new_crvals = self.wcs.wcs_pix2sky([[new_crpix1, new_crpix2]], 1)[0]
      self.hdr['crpix1'] = new_crpix1
      self.hdr['crpix2'] = new_crpix2
      self.hdr['crval1'] = new_crvals[0]
      self.hdr['crval2'] = new_crvals[1]
      if export:
         self.export(self.filename)
      return new_crvals

   def shift_crval(self, new_crval1, new_crval2, export=True):
      """
      Shift WCS to use new values for CRVAL1 and CRVAL2, and write the updated
      values to the header.
      """
      assert self._has_wcs, "self.wcs is not yet defined!"
      new_crpixs = self.wcs.wcs_sky2pix([[new_crval1, new_crval2]], 1)[0]
      self.hdr['crval1'] = new_crval1
      self.hdr['crval2'] = new_crval2
      self.hdr['crpix1'] = new_crpixs[0]
      self.hdr['crpix2'] = new_crpixs[1]
      if export:
         self.export(self.filename) 
      return new_crpixs     

   def display(self):
      os.system('ds9 %s & ' % self.filename)

   def convolve(self, kernel, output):
      """
      Convolves self with a kernel, and save the output.
      """
      kernelimg = pyfits.getdata(kernel)
      convimg = hconvolve(self.data, kernelimg)
      if os.path.exists(output):
         os.remove(output)
      pyfits.append(output, convimg, self.hdr)


class HST_image_4sim(fitsimage):
   """
   Special methods that mostly applies to HST images.
   """
   def __init__(self, filename, make_wcs=True):
      fitsimage.__init__(self, filename, make_wcs)
      self.noiseless = None

   def calc_zeropoint(self):
      zpt = -2.5*np.log10(self.hdr['photflam']) - 5.*np.log10(self.hdr['photplam']) - 2.408
      return zpt

   def makenoiselessimage(self, galfile, simroot, magz=99.0, save=0, gain=1.0, 
                          psffile=""):
      """
      Creates a noiseless convolved image
      """
      assert os.path.exists(psffile), "PSF image %s does not exist." % psffile
      if magz >= 99.0:
         magz = self.calc_zeropoint()
      outfile=simroot+'_sim.fits'
      print outfile,xmax,ymax,galfile,magz,gain
      iraf.unlearn('artdata')
      iraf.unlearn('mkobjects')
      iraf.artdata.dynrange=1.e5
      iraf.mkobjects(outfile, output="", title="", ncols=self.hdr['naxis1'], 
         nlines=self.hdr['naxis2'],header="", background=0.0, objects=galfile,
         xoffset=0., yoffset=0., star="gaussian", radius=0.1,
         beta=2.5, ar=1., pa=0., distance=1., exptime=1., 
         magzero=magz, gain=gain, rdnoise=0., poisson=0,
         seed=2, comments=1)
      imhconvolve(outfile, psffile, outfile, overwrite=True)
      self.noiseless = outfile

   def addsimulated(self, galfile,root,realimage, gain=1.0, psffile="", save=0):
      # simulation = root+'_sim.fits'
      assert self.noiseless != None, "Noiseless image with artificial galaxies not calculated."
      outimage = root+'.fits'
      if os.path.exists(outimage):
         os.remove(outimage)
      # makenoiselessimage(galfile,root+'_sim',magz,xmax,ymax,
      #    save=save,gain=gain,psffile=psffile)
      noiseless_img = pyfits.getdata(self.noiseless)
      simulated_img = self.data + noiseless_img
      pyfits.append(outimage, simulated_img, self.hdr)
      # iraf.imcalc(realimage+","+simulation,outimage,"im1+im2")
      if not save:
        os.remove(self.noiseless)

