#!/usr/bin/env python

import numpy as np
from pyraf import iraf
iraf.stsdas()
import os, sys
import pyfits
from scipy import ndimage
from scipy.ndimage.measurements import center_of_mass
from scipy.ndimage.interpolation import shift, rotate, zoom
from ImageTools import image_moments
from pygoods import sextractor

def make_iracpsf(psfname, oversample=10, interp='sinc', iracpix=0.6,
                 radius=3.0):
   """
   Make an oversample IRAC PSF given the PSF image in the IRAC pixel scale, 
   for use in TFIT.
   psfname: the input IRAC PSF
   oversample: the oversampling factor (default is 10)
   interp: the algorithm of interpolation (default is linear)
           other options are nearest, poly3, poly5, spline3, sinc, lsinc, 
           and drizzle.
   iracpix: the pixel scale (in arcsec) of the input IRAC PSF (default is 0.6)
   radius: the radius of the circularized PSF in arcsec (beyond which is set to 0)
   """
   newscale = int(round(iracpix / oversample * 1000))
   # the new pixel scale in milli-arcsec
   interp = interp.replace('[','_')
   interp = interp.replace(']','_')
   psfroot = os.path.splitext(psfname)[0]
   psfroot = psfroot + '_%dmas' % newscale
   psfroot = psfroot + '_%s' % interp
   size0 = pyfits.getdata(psfname).shape[0]
   size1 = size0 * oversample
   # First, run iraf.imlintran
   factor = 1. / oversample; print "factor", factor
   if os.path.exists('temp1.fits'):
      os.remove('temp1.fits')
   iraf.imlintran(psfname, 'temp1.fits', 0., 0., factor, factor, 
                  ncols=size1, nlines=size1, interpolant=interp,
                  fluxconserve='yes')
   # Second, circularize the PSF, assume that the input PSF has equal numbers 
   # of rows and columns
   imgsize = pyfits.getdata('temp1.fits').shape[0]
   imgcenter = (imgsize + 1) / 2
   print "imgcenter", imgcenter
   pixrad = radius / iracpix * oversample
   pixrad = np.round(pixrad)
   mathstr = 'if ((x-%d)**2 + (y-%d)**2) < %f**2 then im1 else 0' % (imgcenter, imgcenter, pixrad)
   print "mathstr", mathstr
   if os.path.exists('temp2.fits'):
      os.remove('temp2.fits')
   iraf.imcalc('temp1.fits', 'temp2.fits', mathstr)
   # Third, trim the borders; shed border with units of 100 pixels
   nborder = ((imgsize - 2 * pixrad) / 2) / 100
   print "nborder:", nborder
   if nborder > 0:
      begin = nborder * 100 + 1
      end = imgsize - nborder * 100 - 1
      iraf.imcopy('temp2.fits[%d:%d,%d:%d]' % (begin,end,begin,end),
                  psfroot+'.fits')
   else:
      os.system('mv temp2.fits %s' % (psfroot+'.fits'))
   os.system('rm temp*.fits')

def sharpen_iracpsf(psfname, mag, interpolant='sinc'):
   """
   Make PSF sharper by shrinking pixels!
   mag: pixel ratio (mag > 1 means shrinking pixels and make PSF sharper)
   """
   newpsf = os.path.splitext(psfname)[0] + '_mag%.2f.fits' % mag
   if os.path.exists(newpsf):
      os.remove(newpsf)
   oldpsf_array = pyfits.getdata(psfname)
   ncols, nrows = oldpsf_array.shape
   iraf.imlintran(input=psfname, output=newpsf, xrotation=0., yrotation=0., 
                  xmag=mag, ymag=mag, ncols=ncols, nlines=nrows, 
                  interpolant=interpolant)

def sharpen_iracpsf_zoom(psfname, zoom_factor):
   """
   Sharpen PSF using scipy.ndimage.interpolation.zoom.
   zoom < 1.0 makes PSF SHARPER.
   """
   newpsf = os.path.splitext(psfname)[0] + '_zoom%.2f.fits' % zoom_factor
   psf = pyfits.getdata(psfname)
   hdr = pyfits.getheader(psfname)
   zoomed = zoom(psf, zoom_factor)
   if os.path.exists(newpsf):
      os.remove(newpsf)
   pyfits.append(newpsf, zoomed, hdr)

### Below are newer tasks using pyfits and scipy, therefore side-stepping IRAF

def circular_mask(xc, yc, ncol, nlines, r):
   """
   Create an array that masks out the region beyond the circle with center
   (xc, yc) and radius r.
   All of xc, yc, ncol, nlines need to be integers. Radius r does not need to 
   be an integer.
   """
   rmax = np.minimum(ncol-xc, nlines-yc)
   r = np.minimum(r, rmax)
   y, x = np.ogrid[-xc:ncol-xc,-yc:nlines-yc]
   mask = (x**2 + y**2 <= r**2)
   print "shape(mask)", mask.shape
   return mask

def zoom_iracpsf(psfname, oversample=10., spline_order=3, radius=5.0,
                 iracpix=0.6, outname=None, filter_width=5.0, inner_rad=2.0):
   """
   Use scipy.ndimage.interpolate.zoom to oversample the PSF. Use this to 
   make a PSF on the high-res pixel grid from the low-res pixel grid.
   iracpix: input IRAC pixel scale in arcsec.
   radius: radius of the desired PSF image in arcsec
   """
   psf0 = pyfits.getdata(psfname)
   hdr0 = pyfits.getheader(psfname)
   print psf0.shape
   pixrad = radius / iracpix * oversample
   pixrad = np.round(pixrad)
   psf1 = zoom(psf0, oversample, order=spline_order)   
   print "shape(psf1)", psf1.shape
   xc, yc = np.array(psf1.shape) / 2.
   xc = int(np.floor(xc)); yc = int(np.floor(yc))
   print "xc, yc ", xc, yc
   # Center the PSF again
   # In case the low-level noise skews the image moment, we filter the PSF
   # image first
   # Make a circular mask around the center, and only calculate the image
   # moments for the central part; again this is to guard against PSF wings
   # skewing the center of mass
   inner_pixrad = inner_rad / iracpix * oversample
   inner_pixrad = np.round(inner_pixrad)
   cmask1 = circular_mask(xc, yc, psf1.shape[0], psf1.shape[1], inner_pixrad)
   psf1m = np.where(cmask1==True, psf1, 0.)
   psf1_filtered = ndimage.filters.gaussian_filter(psf1m, filter_width)
   # Now calculate the center of mass of the filtered PSF
   cm1 = center_of_mass(psf1_filtered)
   print "CM of the filtered PSF: (%.2f, %.2f)" % tuple(cm1)
   xshift = xc - cm1[0]
   yshift = yc - cm1[1]
   print "xshift, yshift:", xshift, yshift
   psf1 = shift(psf1, [xshift, yshift], order=1,
                                      mode='wrap')
   cmask = circular_mask(xc, yc, psf1.shape[0], psf1.shape[1], pixrad)
   psf1 = np.where(cmask==True, psf1, 0.)
   print "Shifted PSF center: ", center_of_mass(psf1)
   # assume that CDELT1 is in arcsec/pix
   hdr0['cdelt1'] = hdr0['cdelt1'] / oversample
   hdr0['cdelt2'] = hdr0['cdelt2'] / oversample
   # mas_str = '%2d' % abs(int(round(hdr0['cdelt1']*1000.)))
   mas_str = '%2d' % abs(int(round(iracpix / oversample * 1000.)))
   if outname==None:
      outname = os.path.splitext(psfname)[0] + '_%2smas.fits' % (mas_str)
   if os.path.exists(outname):
      os.remove(outname)
   # Trim the borders if there is any
   yc2, xc2 = np.array(psf1.shape) / 2.
   xc2 = int(np.floor(xc2)); yc2 = int(np.floor(yc2))
   xmin = np.maximum(0, xc2 - pixrad*1.2)
   xmax = np.minimum(psf1.shape[1], xc2 + pixrad*1.2)
   ymin = np.maximum(0, yc2 - pixrad*1.2)
   ymax = np.minimum(psf1.shape[0], yc2 + pixrad*1.2)
   print "xmin, xmax, ymin, ymax", xmin, xmax, ymin, ymax
   psf2 = psf1[ymin:ymax,xmin:xmax]
   print "shape(psf2)", np.shape(psf2)
   psf2 = psf2 / psf2.sum()
   pyfits.append(outname, psf2, hdr0)
   return psf1

def combine_rotated_modelpsf(modelpsf, channel, angle_file, norm_const=96.8):
   """
   Performs a weighted sum of the model PSF rotated at different angles based
   on the AOR. Weights are the exposure times of each frame, and both exptime 
   and position angle are read from the CROTA2 and EXPTIME columns from 
   angle_file.
   Note that scipy.ndimage.interpolation.rotate rotates clockwise.
   After this step, one can re-scale the PSF image into the desired pixel 
   scales.
   Definition of CROTA2: "CROTA2 The simplest way is to give CROTA2, the angle 
   between the North and the second axis of the image counted positive to the 
   East, and the coordinate increments CDELT1, CDELT2. The last two are the 
   per pixel increment along RA and DEC."
   So the rotation angle should be -1 * CROTA2.
   """
   c = sextractor(angle_file)
   exptime = c.exptime / norm_const  # most frames have exptime=norm_const
   angles = -1. * c.crota2
   # angles = 360. - c.crota2
   psf0 = pyfits.getdata(modelpsf)
   hdr0 = pyfits.getheader(modelpsf)
   # Now center the image at the measured center of mass
   cm0 = center_of_mass(psf0)
   xc = psf0.shape[0] / 2.
   yc = psf0.shape[0] / 2.
   xshift = xc - cm0[0]
   yshift = yc - cm0[1]
   psf0_centered = shift(psf0, [xshift, yshift])
   psf0_summed = np.zeros(psf0.shape)
   for i in range(len(angles)):
      psf0_rot = rotate(psf0_centered, angles[i], reshape=False)
      psf0_summed = psf0_rot * exptime[i]
   # normalize the PSF
   psf0_summed = psf0_summed / psf0_summed.sum()
   outname = '%s_modelpsf_rotsum.fits' % channel
   if os.path.exists(outname):
      os.remove(outname)
   pyfits.append(outname, psf0_summed, hdr0)


def zoom_model_iracpsf(modelpsf, pixscale=0.06, radius=3.0):
   """
   From the extended model PSF/PRF, calculate the over-sampled IRAC PSF in the 
   desired channel. The resulting pixel scale of the PSF file is given by the
   pixscale parameter in arcsec. The circularization radius (in arcsec) is 
   given by the parameter radius.
   """
   pixrad = np.round(radius / pixscale)
   h = pyfits.open(modelpsf)
   # assume that CDELT1 and CDELT2 are in deg/pix
   pixscale0 = abs(h[0].header['cdelt1'] * 3600.) # convert from deg/pix to arcsec/pi
   print "pixscale0", pixscale0
   # h[0].header['cdelt1'] = -1.*pixscale0
   # h[0].header['cdelt2'] = pixscale0
   # h.flush()
   h.close()
   zoomfactor = pixscale0 / pixscale
   print "zoomfactor", zoomfactor
   
   # Now zoom the PSF
   zoom_iracpsf(modelpsf, oversample=zoomfactor, radius=radius, iracpix=pixscale0)
   # psf0 = pyfits.getdata(modelpsf)
   # hdr = pyfits.getheader(modelpsf)
   # hdr['cdelt1'] = -1. * pixscale / 3600.
   # hdr['cdelt2'] = pixscale / 3600.
   # psf1 = ndimage.interpolation.zoom(psf0, zoomfactor, order=3)
   # newname = os.path.splitext(modelpsf)[0] + '_%2dmas.fits' % int(round(pixscale*1000.))
   # pyfits.append(newname, psf1, hdr)

def restore_modelpsf_cdelt(modelpsf, defvalue=6.7963e-5):
   h = pyfits.open(modelpsf, mode='update')
   h[0].header['cdelt1'] = -1.*defvalue
   h[0].header['cdelt2'] = defvalue
   h.flush()
   h.close()

