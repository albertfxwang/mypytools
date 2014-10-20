#!/usr/bin/env python

import numpy as np
import pyfits 
import pywcs
import os
import pyraf_utils
from PhotomTools import hst_zeropoints

class FITSImage(object):
   # a base class for FITS images
   def __init__(self, filename):
      self.filename = filename
      h = pyfits.open(filename, mode='update')
      self.data = h[0].data
      self.header = h[0].header
      update = 0
      if 'CPDIS1' in self.header.keys():
         h[0].header.remove('cpdis1')
         update = 1
      if 'CPDIS2' in self.header.keys():
         h[0].header.remove('cpdis2')
         update = 1
      if update:
         h.flush()
         h.close()
         h = pyfits.open(filename)
         self.header = h[0].header
      self.wcs = pywcs.WCS(header=self.header, fobj=h)
      self.scale = self.pixscale()

   def pixscale(self):
      """
      Returns pixel scale in arcsec/pixel.
      """
      cd1_1 = self.header['cd1_1']
      try:
         cd1_2 = self.header['cd1_2']
      except:
         cd1_2 = 0.
      s = np.sqrt((cd1_1*3600.)**2 + (cd1_2*3600.)**2)
      return s

   def cut(self, xc, yc, xw, yw, newfile, unit='pixel', overwrite=False):
      """      
      make a cutout around center (xc, yc)
      If unit == 'pixel', xc and yc are pixel coordinates; otherwise, 
      xc, yc should be RA and DEC.
      xw, yw are the *widths* of the cutout; if unit=='pixel', xw, yw are also
      in pixels; if unit=='radec', xw, yw are in arcsec.
      """      
      assert unit in ['pixel', 'radec']
      if os.path.exists(newfile) & (overwrite==False):
         print "File %s already exists! Use overwrite=True to overwrite it." % newfile
         return 0
      # convert to pixel coordinates
      if unit == 'radec':
         xc, yc = self.wcs.wcs_sky2pix([[xc, yc]], 1)[0]
         xw = wx / self.scale
         yw = yw / self.scale
      xlo = int(round(xc - xw / 2.))
      xhi = int(round(xc + xw / 2.))
      ylo = int(round(yc - yw / 2.))
      yhi = int(round(yc + yw / 2.))
      if (overwrite == True) & (os.path.exists(newfile)):
         os.remove(newfile)
      pyraf_utils.imcopy(self.filename, [xlo,xhi,ylo,yhi], newfile)

   def blocksum(self, blksize):
      """
      Calculate the block-summed image. Block size is specified by blksize.
      """
      assert (self.data.shape[0] % blksize == 0), "number of image columns not an integer multiple of blksize..."
      assert (self.data.shape[1] % blksize == 0), "number of image rows not an integer multiple of blksize..."
      step1 = np.add.reduceat(self.data, np.arange(self.data.shape[0])[::blksize], axis=0)
      step2 = np.add.reduceat(step1, np.arange(self.data.shape[1])[::blksize], axis=1)
      # print "Old image shape: ", self.data.shape
      # print "New image shape: ", step2.shape
      self.data = step2
      self.scale = self.scale * blksize
      # Also need to update WCS keywords
      self.header['crpix1'] = self.header['crpix1'] / float(blksize)
      self.header['crpix2'] = self.header['crpix2'] / float(blksize)
      self.header['cd1_1'] = self.header['cd1_1'] * blksize
      self.header['cd2_2'] = self.header['cd2_2'] * blksize
      if self.header.has_key('cd1_2'):
         self.header['cd1_2'] = self.header['cd1_2'] * blksize
      if self.header.has_key('cd2_1'):
         self.header['cd2_1'] = self.header['cd2_1'] * blksize

   def writeto(self, name=None, overwrite=False):
      # update the FITS image with the current values of data and header
      if name == None:
         # overwrite the FITS file
         h = pyfits.open(self.filename, mode='update')
         h[0].data = self.data
         h[0].header = self.header
         h.flush()
         h.close()
      else:
         if os.path.exists(name):
            if (overwrite==False):
               raise ValueError, "%s already exists. Set overwrite to True?" % name
            else:
               os.remove(name)
         pyfits.append(name, self.data, self.header)

class HSTImage(FITSImage):
   def ABzeropoint(self):
      """
      Calculate zeropoint in AB mag using header keywords PHOTFLAM and PHOTPLAM.
      """
      photflam = self.header['photflam']
      photplam = self.header['photplam']
      zpt = hst_zeropoints.wfc3_abmag_zeropoint(photflam, photplam)
      return zpt

   def effgain(self):
      """
      Calculate the effective gain (ccdgain * exptime) for images in counts/sec.
      """
      ccdgain = self.header['ccdgain']
      exptime = self.header['exptime']
      g = ccdgain * exptime
      return g


class WeightImage(FITSImage):
   def __init__(self, filename):
      super(WeightImage, self).__init__(filename)

   def make_rms(self, rmsimage, cval=1.e10, norm=False):
      """
      Make an RMS image out of a weight image.
      """
      if norm:
         wht = self.data / self.header['exptime']  
         # remember to divide by exposure time, otherwise the unit of rms will 
         # be wrong
      else:
         wht = self.data
      rms = np.where(wht>0, 1./np.sqrt(wht), cval)
      if os.path.exists(rmsimage):
         os.remove(rmsimage)
      pyfits.append(rmsimage, rms, self.header)

   def make_flg(self, flgimage, flgvalue=1):
      """
      Make a flag image out of a weight image.
      """
      flg = np.where(self.data > 0, 0, 1).astype('int16')
      if os.path.exists(flgimage):
         os.remove(flgimage)
      pyfits.append(flgimage, flg, self.header)
      
