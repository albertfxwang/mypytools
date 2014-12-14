#!/usr/bin/env python

import numpy as np
import pyfits, pywcs
import os
import hconvolve 
import pyfits
from scipy.ndimage import filters

def flatten_pixcoord(image):
   """
   Construnct an array representing the pixel coordinates of every pixel in 
   the image, in the shape of (Npix, 2).
   """
   hdr = pyfits.getheader(image)
   naxis1 = hdr['naxis1']
   naxis2 = hdr['naxis2']
   pixcoord = np.mgrid[1:naxis1+1,1:naxis2+1]
   pixcoord = pixcoord.swapaxes(0,1).swapaxes(1,2)
   pixcoord = pixcoord.ravel().reshape(naxis1*naxis2, 2)
   return pixcoord

def transform_pixcoord(source, dest):
   """
   Transform the pixel coordinates from the source image to the dest(ination)
   image.
   """
   wcs_src = pywcs.WCS(pyfits.getheader(source))
   wcs_dest = pywcs.WCS(pyfits.getheader(dest))
   coord_src = flatten_pixcoord(source)
   sky_src = wcs_src.wcs_pix2sky(coord_src, 1)
   coord_dest = wcs_dest.wcs_sky2pix(sky_src, 1)
   return np.around(coord_dest).astype('int')

def drizzle_mask(source, dest, output):
   """
   ** The final task in making the mask for a light map. **
   Given a source segmentation map, make a new mask on the same pixel grid as 
   the destination image. The source should be a mask image with value 1 for 
   the pixels to keep, and 0 for the pixels to be masked out.
   """
   if os.path.exists(output):
      os.system('rm %s' % output)
   coord_src = flatten_pixcoord(source)
   print "Transforming coordinates..."
   coord_dest = transform_pixcoord(source, dest)
   hdr_dest = pyfits.getheader(dest)
   mask_in = pyfits.getdata(source)
   mask_out = np.zeros((hdr_dest['naxis1'], hdr_dest['naxis2']))
   # Now this might take a lot of time... is there a more clever way?
   print "Start looping all the pixels, might take a while..."
   for i in range(coord_src.shape[0]):
      x0, y0 = coord_src[i]
      if mask_in[y0-1,x0-1] > 0:
         x1, y1 = coord_dest[i]
         mask_out[y1-1,x1-1] = 1
   # write output
   pyfits.append(output, data=mask_out, header=hdr_dest)

def blur_mask(mask, blur_kernel, threshold=0.1):
   """
   Given an input mask image (numpy array) of values 0 or 1 (1 being accepted 
   pixels), convolve the mask by blur_kernel (normalized to 1) and set pixels 
   above the threshold to value 1 and then the rest 0.
   """
   k = pyfits.getdata(blur_kernel)
   k = k / k.sum()
   mask = hconvolve.hconvolve(mask, k)
   mask = np.where(mask >= threshold, 1, 0).astype('int')
   return mask


def filter_segmap(segimage, id_keep, output, blur_kernel="", threshold=0.1):
   """
   Specify a list of ID numbers to keep, and zero-out the rest of the 
   segmentation map.
   """
   seg = pyfits.getdata(segimage)
   mask = np.zeros(seg.shape, 'int')
   # Loop through all IDs... is there a better way??
   for x in id_keep:
      mask = np.where(seg==x, 1, mask)
   seg_masked = np.where(mask==1, 1, 0)
   if os.path.exists(output):
      os.system('rm %s' % output)
   # Now convolve with a blurring kernel if desired
   if len(blur_kernel):
      mask = blur_mask(mask, blur_kernel, threshold=threshold)
      # k = pyfits.getdata(blur_kernel)
      # mask = hconvolve.hconvolve(mask, )
   pyfits.append(output, data=seg_masked, header=pyfits.getheader(segimage))
   return mask

def make_irac_lightmap(id_keep, hr_segmap, hr_mask, irac_psf, irac_drz, irac_output, blur_threshold=0.1, sigma=1.0):
   """
   Make an IRAC map for cluster members (or an arbitrary set of ID numbers).
   id_keep: ID numbers in the high-res segmentation map to keep
   hr_segmap: high-res segmentation map
   hr_mask: high-res mask image (an intermediate product)
   irac_psf: PSF in IRAC, used to blur the high-res mask image
   irac_drz: IRAC science image, onto which we drizzle the high-res mask image
   irac_output: file name of the output IRAC light map
   blur_threshold: the threshold (between 0 and 1) in the step of blurring the 
                   high-res mask image.
   """
   # First step, zero-out the non cluster members
   mask = filter_segmap(hr_segmap, id_keep, hr_mask, blur_kernel=irac_psf, 
                        threshold=blur_threshold)
   # Now we have a mask image in high-res, drizzle the pixels onto the low-res
   # pixel grid
   if os.path.exists("irac_mask.fits"):
      os.system('rm irac_mask.fits')
   drizzle_mask(hr_mask, irac_drz, "irac_mask.fits")
   irac_input = pyfits.getdata(irac_drz)
   irac_mask = pyfits.getdata("irac_mask.fits")
   irac_map = np.where(irac_mask>0, irac_input, 0.)
   # Also smooth the output light map with a Gaussian kernel
   if sigma > 0:
      print "Smoothing the IRAC mask..."
      irac_map = filters.gaussian_filter(irac_map, sigma)
   irac_hdr = pyfits.getheader(irac_drz)
   os.system('rm %s' % irac_output)
   pyfits.append(irac_output, data=irac_map, header=irac_hdr)
   print "Done."



