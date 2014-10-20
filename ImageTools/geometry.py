#!/usr/bin/env python

import numpy as np
from pyraf import iraf
iraf.images()
import pyfits
import os, sys, subprocess
from pygoods import sextractor
from drizzlepac import updatehdr
import pywcs

### Performs astrometric alignments of input images to reference images using
### SExtractor catalogs.

class AlignImage(object):
   def __init__(self, input_image, ref_image, input_rms, ref_rms, input_flg, ref_flg, coofile):
      ### Initialization
      self.input_image = input_image
      self.ref_image = ref_image
      self.input_rms = input_rms
      self.ref_rms = ref_rms
      self.input_seg = input_rms.replace('rms', 'geo_seg')
      self.ref_seg = ref_rms.replace('rms', 'geo_seg')
      self.input_flg = input_flg
      self.ref_flg = ref_flg
      self.input_cat = os.path.splitext(self.input_image)[0] + '.cat'
      self.ref_cat = os.path.splitext(self.ref_image)[0] + '.cat'
      # initialize WCS info
      self.hdr_ref = pyfits.getheader(self.ref_image)
      self.wcs_ref = pywcs.WCS(self.hdr_ref)
      if self.hdr_ref.has_key('cd1_1'):
         cd1_1 = self.hdr_ref['cd1_1']
         if self.hdr_ref.has_key('cd1_2'):
            cd1_2 = self.hdr_ref['cd1_2']
         else:
            cd1_2 = 0.
         self.ref_pixscale = 3600. * np.sqrt(cd1_1**2 + cd1_2**2)
      else:
         self.ref_pixscale = np.abs(self.hdr_ref['cdelt1']) * 3600.
      # Also find x0ref_old, y0ref_old
      self.hdr_in = pyfits.getheader(self.input_image)
      cd1_1 = self.hdr_in['cd1_1']
      if self.hdr_in.has_key('cd1_2'):
         cd1_2 = self.hdr_in['cd1_2']
      else:
         cd1_2 = 0.
      self.input_pixscale = np.sqrt(cd1_1**2 + cd1_2**2) * 3600.
      self.wcs_in = pywcs.WCS(self.hdr_in)
      self.coofile = coofile

   def run_sextractor(self, sexconf, magzero=25.0, detect_thresh=10., sexcmd='cex'):
      ### Run SExtractor on both input and reference images to produce catalogs
      ### that iraf.geomap will use to calculate the transformations. SExtractor
      ### should be configured to only detect bright sources.
      sexcmd_input = [sexcmd, self.input_image, '-c', sexconf]
      sexcmd_input += ['-CATALOG_NAME', self.input_cat]
      sexcmd_input += ['-DETECT_THRESH', str(detect_thresh)]
      sexcmd_input += ['-MAG_ZEROPOINT', str(magzero)]
      sexcmd_input += ['-FLAG_IMAGE', self.input_flg]
      sexcmd_input += ['-FLAG_TYPE', 'OR']
      sexcmd_input += ['-WEIGHT_TYPE', 'MAP_RMS']
      sexcmd_input += ['-WEIGHT_THRESH', str(1.e10)]
      sexcmd_input += ['-WEIGHT_IMAGE', self.input_rms]
      sexcmd_input += ['-CHECKIMAGE_TYPE', 'SEGMENTATION']
      sexcmd_input += ['-CHECKIMAGE_NAME', self.input_seg]
      print "Run SExtractor on input images..."
      print ' '.join(sexcmd_input)
      subprocess.call(sexcmd_input)
      ## Now run SExtractor on the reference image as well
      sexcmd_ref = [sexcmd, self.ref_image, '-c', sexconf]
      sexcmd_ref += ['-CATALOG_NAME', self.ref_cat]
      sexcmd_ref += ['-DETECT_THRESH', str(detect_thresh)]
      sexcmd_ref += ['-MAG_ZEROPOINT', str(magzero)]
      sexcmd_ref += ['-FLAG_IMAGE', self.ref_flg]
      sexcmd_ref += ['-FLAG_TYPE', 'OR']
      sexcmd_ref += ['-WEIGHT_TYPE', 'MAP_RMS']
      sexcmd_ref += ['-WEIGHT_THRESH', str(1.e10)]
      sexcmd_ref += ['-WEIGHT_IMAGE', self.ref_rms]
      sexcmd_ref += ['-CHECKIMAGE_TYPE', 'SEGMENTATION']
      sexcmd_ref += ['-CHECKIMAGE_NAME', self.ref_seg]
      print "Run SExtractor on refernce images..."
      print ' '.join(sexcmd_ref)
      subprocess.call(sexcmd_ref)

   def get_xy_input(self, numbers):
      # get the (x, y) coordinates for a list of objects that will be used later
      # for geomap
      N = len(numbers)
      # self.x_input, self.y_input should be pais of (x, y) coordinates
      # they should match self.x_ref, self.y_ref defined later.
      c = sextractor(self.input_cat)
      self.x_input = np.zeros(N)
      self.y_input = np.zeros(N)
      for i in range(len(numbers)):
         self.x_input[i] = c.x_image[c.number==numbers[i]][0]
         self.y_input[i] = c.y_image[c.number==numbers[i]][0]

   def get_xy_ref(self, numbers):
      # get the (x, y) coordinates for a list of objects that will be used later
      # for geomap
      N = len(numbers)
      # self.x_ref, self.y_ref should be pais of (x, y) coordinates
      # they should match self.x_input, self.y_input defined later.
      c = sextractor(self.ref_cat)
      self.x_ref = np.zeros(N)
      self.y_ref = np.zeros(N)
      for i in range(len(numbers)):
         self.x_ref[i] = c.x_image[c.number==numbers[i]][0]
         self.y_ref[i] = c.y_image[c.number==numbers[i]][0]

   def write_xy_map(self, coofile=""):
      # write a coordinate file that matches the input image to the reference image
      # Note that our definition of input v.s. reference is reversed from the 
      # definition in IRAF documentation
      if not coofile: coofile = self.coofile
      f = open(coofile, 'wb')
      if not hasattr(self, 'x_ref'):
         raise ValueError, "run self.get_xy_ref and self.get_xy_input first!"
      for i in range(len(self.x_ref)):
         f.write('%.2f  %.2f  ' % (self.x_input[i], self.y_input[i]))
         f.write('%.2f  %.2f  ' % (self.x_ref[i], self.y_ref[i]))
         f.write('\n')
      f.close()

   def geomap(self, database, xmax=-1, ymax=-1):
      if not self.coofile:
         raise ValueError, "Specify self.coofile first!"
      if os.path.exists(database):
         os.remove(database)
      if xmax < 0: xmax = self.hdr_in['naxis1']
      if ymax < 0: ymax = self.hdr_in['naxis2']
      iraf.geomap(self.coofile, database, 1, xmax, 1, ymax, interactive='no')
      self.database = database

   def geotran(self, database, coofile):
      # Run geotran on self.input_image... is this what I want to do?
      output_image = os.path.splitext(self.input_image)[0] + '_geotran.fits'
      if os.path.exists(output_image):
         os.remove(output_image)
      h = pyfits.getheader(self.input_image)
      iraf.geotran(self.input_image, output_image, database, coofile, 
                   geometry='geometric', 
                   xscale=1.0, yscale=1.0, interpolant='spline3')

   def read_shifts(self, database):
      # Read xshift and yshift from database. These shifts are in units of 
      # REFERENCE IMAGE pixels. Here we will convert the shift into degrees,
      # which we will apply to shift the CRVALs of the input image. Also read
      # the amount of rotation.
      # Essentially, xshift and yshift are the pixel coordinate values in the 
      # REFERENCE image for the ORIGIN (1, 1) of the INPUT image, or 
      # (x0ref, y0ref), after alignment. 
      # So the true "shift" should be the difference between (x0ref, y0ref)
      # and (x0ref_old, y0ref_old).
      # Find the values for xshift, yshift, and xrot
      f = open(database)
      lines = f.readlines()
      for line in lines:
         l = line.split()
         if l[0] == 'xshift':
            x0ref = float(l[1])
         elif l[0] == 'yshift':
            y0ref = float(l[1])
         elif l[0] == 'xrotation':
            xrot = float(l[1])
      # enforce xrot to between -180. and 180.
      if xrot > 180.:
         xrot = xrot - 360.
      # convert xshift, yshift to arcsec
      # Need to find the pixel scale (in arcsec) of the reference image
      ra0, dec0 = self.wcs_in.wcs_pix2sky([[1,1]], 1)[0]
      x0ref_old, y0ref_old = self.wcs_ref.wcs_sky2pix([[ra0, dec0]], 1)[0]
      # Now calculate xshift and yshift in arcsec
      xshift = (x0ref - x0ref_old) * self.ref_pixscale
      yshift = (y0ref - y0ref_old) * self.ref_pixscale
      return xshift, yshift, xrot
      
   def apply_shift(self, database, input_image="", suffix='geo', force_north_up=True):
      # Apply the shifts to input image WCS header keywords
      # Once the database has been calculated, one can use this to shift many
      # other images.
      # First read the shifts (in arcsec) and rotation (in degrees)
      xshift, yshift, xrot = self.read_shifts(database)
      # Then calculate the shifts in INPUT IMAGE PIXELS
      xshift_pix = xshift / self.input_pixscale
      yshift_pix = yshift / self.input_pixscale
      print "xshift_pix, yshift_pix", xshift_pix, yshift_pix
      # Now call updatehdr.updatewcs_with_shift
      if not input_image:
         input_image = self.input_image
      output_image = os.path.splitext(input_image)[0] + '_%s.fits' % suffix
      if os.path.exists(output_image):
         print "Removing %s..." % output_image
         os.remove(output_image)
      os.system('cp %s %s' % (input_image, output_image))
      updatehdr.updatewcs_with_shift(output_image, input_image, rot=xrot,
                                     scale=1.0, xsh=-xshift_pix, ysh=-yshift_pix, 
                                     force=True, verbose=True)
      if force_north_up:
         h = pyfits.open(output_image, mode='update')
         pixscale = np.sqrt(h[0].header['cd1_1']**2 + h[0].header['cd1_2']**2) # in degrees
         h[0].header['cd1_1'] = np.sign(h[0].header['cd1_1']) * pixscale
         h[0].header['cd1_2'] = 0.
         h[0].header['cd2_1'] = 0.
         h[0].header['cd2_2'] = np.sign(h[0].header['cd2_2']) * pixscale
         h[0].header['orientat'] = 0.
         h.flush()
         h.close()

   def run_alignment(self, reflist, database, suffix='geo'):
      # Run all alignment steps AFTER running SExtractor and identified a list
      # of bright objects to calculate shifts.
      print "If you have not run SExtractor on both the input and the reference"
      print "image... do it now and identify a list of objects to calculate the"
      print "shifts. Then write the segmentation numbers of these objects to "
      print "the file as the reflist argument."
      r = sextractor(reflist)
      self.get_xy_input(r.segnum_in)
      self.get_xy_ref(r.segnum_ref)
      self.write_xy_map(self.coofile)
      self.geomap(database)
      self.apply_shift(database, suffix=suffix)

def paste_image(source, dest, output):
   # Paste the image data from one file onto the pixel grid defined by another
   # file. The user should make sure that the source and destination images
   # have the same pixel scales and orientations.
   import images
   s = images.FITSImage(source)  # source of imaging data
   d = images.FITSImage(dest)    # defines the destination pixel grid
   # Check if the two images have the same pixel scale
   assert round(s.scale, 3) == round(d.scale, 3), "Pixel scales of source and desination are different!"
   # calculate the pixel coordinates (defined by dest) of the edges of the
   # source image
   ra0, dec0 = s.wcs.wcs_pix2sky([[1, 1]], 1)[0]
   x0d, y0d = d.wcs.wcs_sky2pix([[ra0, dec0]], 1)[0]
   # Make sure to correct CRPIX1 and CRPIX2 as well
   dx = x0d - np.round(x0d)
   dy = y0d - np.round(y0d)
   print "x0d, y0d before rounding up:", x0d, y0d
   x0d = int(np.round(x0d))
   y0d = int(np.round(y0d))
   print "x0d, y0d", x0d, y0d
   assert ((x0d>0) & (y0d>0)), "Source image extends outside of destination grid."
   nx, ny = s.header['naxis1'], s.header['naxis2']
   # ra1, dec1 = s.wcs.wcs_pix2sky([[nx, ny]], 1)[0]
   # x1d, y1d = d.wcs.wcs_sky2pix([[ra1, dec1]], 1)[0]
   # print "x1d, y1d before rounding up:", x1d, y1d
   # x1d = np.floor(x1d); y1d = np.floor(y1d)
   x1d = x0d + nx
   y1d = y0d + ny
   print "x1d, y1d", x1d, y1d
   assert ((x1d<=d.header['naxis1']) & (y1d<=d.header['naxis2'])), "Source image extends outside of destination grid."
   # Now paste the image data
   newdata = np.zeros(d.data.shape)
   newdata[y0d-1:y1d-1,x0d-1:x1d-1] = s.data.copy()
   newheader = s.header.copy()
   # Change the WCS keywords to match the destination
   for k in ['cd1_1', 'cd1_2', 'cd2_1', 'cd2_2']:
      if d.header.has_key(k):
         newheader[k] = d.header[k]
   newwcs = pywcs.WCS(d.header)
   newcrpix1 = d.header['crpix1'] - dx
   newcrpix2 = d.header['crpix2'] - dy
   newcrval1, newcrval2 = newwcs.wcs_pix2sky([[newcrpix1, newcrpix2]], 1)[0]
   newheader['crpix1'] = newcrpix1
   newheader['crpix2'] = newcrpix2
   # newheader['crval1'] = newcrval1
   # newheader['crval2'] = newcrval2
   newheader['crval1'] = d.header['crval1']
   newheader['crval2'] = d.header['crval2']
   pyfits.writeto(output, newdata, header=newheader, clobber=True)
