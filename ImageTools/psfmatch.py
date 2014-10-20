#!/usr/bin/env python

import numpy as np
import pyfits
# from scipy import signal, ndimage
from hconvolve import fftdeconvolve
import os
from scipy import ndimage

def psfmatch(psf_ref, psf_2match, kernelname):
   """
   Derive the kernel that matches psf_2match to psf_ref, so that psf_2match,
   when convolved with the kernel, gives psf_ref.
   Make sure that both PSF images have the same pixel scales, are centered, and
   have the same image size.
   """
   psf1 = pyfits.getdata(psf_ref)
   psf2 = pyfits.getdata(psf_2match)
   kernel = fftdeconvolve(psf1, psf2)   
   # normalize the kernel
   kernel = kernel / kernel.sum()
   hdr2 = pyfits.getheader(psf_2match)
   if os.path.exists(kernelname):
      os.remove(kernelname)
   pyfits.append(kernelname, kernel.real, hdr2)

def zoom_psf(psfname, zoom=None, cutpad=0, newpixscale=0.06, psfsize=3.9, norm=True):
   """
   Rescale ACS PSF using scipy.ndimage.interpolation.zoom to match the 
   WFC3/IR pixel scale.
   zoom < 1.0 makes PSF SHARPER.
   Assume that psf image is square.
   """
   psf = pyfits.getdata(psfname)
   hdr = pyfits.getheader(psfname)
   if zoom == None:
      oldsize = psf.shape[0]
      newsize = int(round(psfsize / newpixscale))
      zoom = float(newsize) / oldsize
   zoomed = ndimage.interpolation.zoom(psf, zoom)
   if norm:
      zoomed = zoomed / zoomed.sum()
   print zoomed.shape
   newpsf = os.path.splitext(psfname)[0] + '_zoomed.fits'
   if cutpad > 0:
      zoomed = zoomed[cutpad:-cutpad,cutpad:-cutpad]
   if os.path.exists(newpsf):
      os.remove(newpsf)
   pyfits.append(newpsf, zoomed, hdr)