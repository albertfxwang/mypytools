#!/usr/bin/env python

import numpy as np
import pyfits
import os, sys
import pywcs
import yaml

def nearest_half_integer(a):
   """
   Return the nearest half-integer.
   """
   a_dec = a - int(a)
   diff_a_dec = np.sign(a) * 0.5 - a_dec
   return round(a + diff_a_dec, 1)

def wcs_swarp(image, crvals, pixscale):
   """pixscale is the desired pixel scale in arcsec."""
   hdr = pyfits.getheader(image)
   wcs = pywcs.WCS(hdr)
   cdmatrix = np.array([[-1.*pixscale/3600., 0.], [0., pixscale/3600.]])
   wcs.wcs.cd = cdmatrix
   crpix = wcs.wcs_sky2pix([crvals], 1)
   return crpix[0] 

def update_crvals_crpixs(image, crvals, crpixs):
   hdu = pyfits.open(image, mode='update')
   hdu[0].header['crval1'] = crvals[0]
   hdu[0].header['crval2'] = crvals[1]
   hdu[0].header['crpix1'] = crpixs[0]
   hdu[0].header['crpix2'] = crpixs[1]
   hdu.flush()
   hdu.close()

def swarp_images(swarp_params):
   c = yaml.load(open(swarp_params))
   
   # Step 1: re-calculate the CRVALs that correspond to the nearest half-integer
   # pixels to the current CRPIXs. Therefore we do not SWarp the low-res images.
   # To ensure that CRPIXs are positive for the high-res image, we should start
   # with the CRVALs of the high-res image, find the image pixel coordinates
   # of these CRVALs in the low-res image, then find the nearest half-integer
   # coordinates in the low-res image, and finally go back to the high-res 
   # image.
   lores_hdr = pyfits.getheader(os.path.join(c['lores_dir'], c['lores_drz']))
   hires_hdr0 = pyfits.getheader(os.path.join(c['hires_dir'], c['hires_input_drz']))
   # print "Original CRVALs: %f, %f" % (lores_hdr['crval1'], lores_hdr['crval2'])
   # First, get the CRVALs in the high-res image
   print "Original CRVALs: %f, %f" % (hires_hdr0['crval1'], hires_hdr0['crval2'])
   crvals0 = [hires_hdr0['crval1'], hires_hdr0['crval2']]
   lores_wcs = pywcs.WCS(lores_hdr)
   # Calculate the pixel coordinates in the low-res frame for the CRVALs of the
   # high-res frame
   crpixs0_lr = lores_wcs.wcs_sky2pix([crvals0], 1)[0]
   lr_crpix1 = crpixs0_lr[0]
   lr_crpix2 = crpixs0_lr[1]
   # lr_crpix1 = lores_hdr['crpix1']
   # lr_crpix2 = lores_hdr['crpix2']
   # Now calculate the nearest half-integer CRPIXs for the low-res frame
   lr_crpix1_new = nearest_half_integer(lr_crpix1)
   lr_crpix2_new = nearest_half_integer(lr_crpix2)
   crvals_new = lores_wcs.wcs_pix2sky([[lr_crpix1_new, lr_crpix2_new]], 1)[0]
   print "New CRVALs: %f, %f" % (crvals_new[0], crvals_new[1])
   lr_crpixs_new = [lr_crpix1_new, lr_crpix2_new]
   # Now update the headers of all lo-res images
   lr_output_drz = os.path.join(c['lores_dir'], c['lores_drz'])
   # lr_output_drz_shift = lr_output_drz[:-5] + '_shift.fits'
   lr_output_wht = os.path.join(c['lores_dir'], c['lores_unc'])
   # lr_output_wht_shift = lr_output_wht[:-5] + '_shift.fits'
   # for f in [lr_output_drz_shift, lr_output_wht_shift]:
      # if os.path.exists(f):
      #    os.remove(f)
   # os.system('cp %s %s' % (lr_output_drz, lr_output_drz_shift))
   # os.system('cp %s %s' % (lr_output_wht, lr_output_wht_shift))
   # update_crvals_crpixs(lr_output_drz_shift, crvals_new, lr_crpixs_new)
   # update_crvals_crpixs(lr_output_wht_shift, crvals_new, lr_crpixs_new)
   update_crvals_crpixs(lr_output_drz, crvals_new, lr_crpixs_new)
   update_crvals_crpixs(lr_output_wht, crvals_new, lr_crpixs_new)

   # Step 2: Write the .head file for the SWarp output of the high-res image
   hr_input_drz = os.path.join(c['hires_dir'], c['hires_input_drz'])
   hr_input_wht = os.path.join(c['hires_dir'], c['hires_input_wht'])
   hr_output_drz = os.path.join(c['hires_dir'], c['hires_output_drz'])
   hr_output_wht = os.path.join(c['hires_dir'], c['hires_output_wht'])
   for f in [hr_output_drz, hr_output_wht]:
      if os.path.exists(f):
         os.remove(f)
   hires_input_hdr = pyfits.getheader(hr_input_drz)
   # hr_crpix1 = hires_input_hdr['crpix1']
   # hr_crpix2 = hires_input_hdr['crpix2']
   hr_wcs = pywcs.WCS(hires_input_hdr)
   hr_crpix = hr_wcs.wcs_sky2pix([crvals_new], 1)[0]
   print hr_crpix[0], hr_crpix[1]
   hr_crpix1_new = nearest_half_integer(hr_crpix[0])
   hr_crpix2_new = nearest_half_integer(hr_crpix[1])
   print "New CRPIX's for high-res image: %.1f, %.1f" % (hr_crpix1_new, hr_crpix2_new)
   headfile = hr_output_drz[:-4] + 'head'
   f = open(headfile, 'wb')
   # Remember to reserve 8 character spaces for FITS header keywords. So DO NOT
   # put the "=" sign before the 9th character of each line!
   # f.write("""RADECSYS= 'FK5     '
   f.write("""EQUINOX =  2000.0
CTYPE1  = 'RA---TAN'
CUNIT1  = 'deg     '
CTYPE2  = 'DEC--TAN'
CUNIT2  = 'deg     '
""")
   f.write("CRVAL1  =  %.10f\n" % crvals_new[0])
   f.write("CRVAL2  =  %.10f\n" % crvals_new[1])
   f.write("CRPIX1  =  %.1f\n" % hr_crpix1_new)
   f.write("CRPIX2  =  %.1f\n" % hr_crpix2_new)
   f.write("CD1_1   =  %.10f\n" % (-1.*c['hr_scale']/3600.))
   f.write("CD1_2   =  0\n")
   f.write("CD2_1   =  0\n")
   f.write("CD2_2   =  %.10f\n" % (c['hr_scale']/3600.))
   f.write("END\n")
   f.close()

   # Step 3: Run SWarp
   swarp_cmd = 'swarp %s ' % hr_input_drz
   swarp_cmd += '-c %s ' % c['swarp_file']
   swarp_cmd += '-IMAGEOUT_NAME %s ' % hr_output_drz
   swarp_cmd += '-WEIGHTOUT_NAME %s ' % hr_output_wht
   swarp_cmd += '-WEIGHT_IMAGE %s ' % hr_input_wht
   print swarp_cmd
   os.system(swarp_cmd)

   # Now check the swarped images
   # os.system('ds9 %s %s &' % (hr_output_drz, lr_output_drz_shift))
   os.system('ds9 %s %s &' % (hr_output_drz, lr_output_drz))
   # lr_hdr_new = pyfits.getheader(lr_output_drz_shift)
   lr_hdr_new = pyfits.getheader(lr_output_drz)
   hr_hdr_new = pyfits.getheader(hr_output_drz)
   print "Low-res image CRVALs: %f, %f" % (lr_hdr_new['crval1'], lr_hdr_new['crval2'])
   print "Low-res image CRPIXs: %.1f, %.1f" % (lr_hdr_new['crpix1'], lr_hdr_new['crpix2'])
   print "High-res image CRVALs: %f, %f" % (hr_hdr_new['crval1'], hr_hdr_new['crval2'])
   print "High-res image CRPIXs: %.1f, %.1f" % (hr_hdr_new['crpix1'], hr_hdr_new['crpix2'])
   print "Done!"
   # return lr_output_drz_shift, lr_output_wht_shift

def shift_crpix(hires_drz, hires_wht, lores_drz, lores_unc, hr_dir='.', lr_dir='.'):
   """
   Shift the CRPIX's to half-integers, re-calculate the CRVALS, and then 
   match the CRVALS of hi- and lo-res images.
   If hires_drz and lores_drz have the same CRVALs, and their pixel sizes
   are integer multiples of each other, then making CRPIXs half-integer in the
   low-res image should automatically make the CRPIXs of high-res image half-
   integers.
   """
   hi = pyfits.open(os.path.join(hr_dir, hires_drz), mode='update')
   hi_wht = pyfits.open(os.path.join(hr_dir, hires_wht), mode='update')
   lo = pyfits.open(os.path.join(lr_dir, lores_drz), mode='update')
   lo_unc = pyfits.open(os.path.join(lr_dir, lores_unc), mode='update')
   hi_hdr = hi[0].header
   lo_hdr = lo[0].header
   crpix1_lo = lo_hdr['crpix1']
   crpix2_lo = lo_hdr['crpix2']
   print "Current CRPIXs for the low-res image: %.2f, %.2f" % (crpix1_lo, crpix2_lo)
   # calculate the new half-integer CRPIX's
   crpix1_lo_new = nearest_half_integer(crpix1_lo)
   crpix2_lo_new = nearest_half_integer(crpix2_lo)
   print "New CRPIXs for the low-res image will be: %.2f, %.2f" % (crpix1_lo_new, crpix2_lo_new)
   # construct WCS instances 
   wcs_hires = pywcs.WCS(hi_hdr)
   wcs_lores = pywcs.WCS(lo_hdr)
   # Re-calculate CRVALs
   crvals_lo_new = wcs_lores.wcs_pix2sky([[crpix1_lo_new, crpix2_lo_new]], 1)[0]
   print "New CRVALs: %.10f, %.10f" % (crvals_lo_new[0], crvals_lo_new[1])
   # Re-calculate CRPIXs for high-res image
   crpix_hi_new = wcs_hires.wcs_sky2pix([crvals_lo_new], 1)[0]
   print "New CRPIXs for hi-res image: %.2f, %.2f" % (crpix_hi_new[0], crpix_hi_new[1])
   if np.abs(crpix_hi_new[0] - nearest_half_integer(crpix_hi_new[0])) > 0.01:
      crpix_hi_new[0] = nearest_half_integer(crpix_hi_new[0])
      print "Enforcing new CRPIX1 to %.2f..." % (crpix_hi_new[0])
   if np.abs(crpix_hi_new[1] - nearest_half_integer(crpix_hi_new[1])) > 0.01:
      crpix_hi_new[1] = nearest_half_integer(crpix_hi_new[1])
      print "Enforcing new CRPIX2 to %.2f..." % (crpix_hi_new[1])
   # Now update the headers
   hi[0].header['crval1'] = crvals_lo_new[0]
   hi[0].header['crval2'] = crvals_lo_new[1]
   hi[0].header['crpix1'] = crpix_hi_new[0]
   hi[0].header['crpix2'] = crpix_hi_new[1]
   hi_wht[0].header['crval1'] = crvals_lo_new[0]
   hi_wht[0].header['crval2'] = crvals_lo_new[1]
   hi_wht[0].header['crpix1'] = crpix_hi_new[0]
   hi_wht[0].header['crpix2'] = crpix_hi_new[1]
   lo[0].header['crval1'] = crvals_lo_new[0]
   lo[0].header['crval2'] = crvals_lo_new[1]
   lo[0].header['crpix1'] = crpix1_lo_new
   lo[0].header['crpix2'] = crpix2_lo_new
   lo_unc[0].header['crval1'] = crvals_lo_new[0]
   lo_unc[0].header['crval2'] = crvals_lo_new[1]
   lo_unc[0].header['crpix1'] = crpix1_lo_new
   lo_unc[0].header['crpix2'] = crpix2_lo_new
   hi.flush()
   hi.close()
   lo.flush()
   lo.close()
