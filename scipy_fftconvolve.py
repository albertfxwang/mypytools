#!/usr/bin/env python

import numpy as np
from scipy.signal import fftconvolve
import pyfits
import os

# Use scipy.signal.fftconvolve to do convolution, and then return an array
# with the same shape as array1.

def scipy_fftconvolve(array1, array2):
   output = fftconvolve(array1, array2, mode='same')
   m, n = array2.shape
   output = output[m/2:-m/2,n/2:-n/2]
   return output

def scipy_imfftconvolve(image1, image2, output_image):
   # provide FITS image file names
   array1 = pyfits.getdata(image1)
   array2 = pyfits.getdata(image2)
   output = scipy_fftconvolve(array1, array2)
   if os.path.exists(output_image):
      os.remove(output_image)
   hdr1 = pyfits.getheader(image1)
   pyfits.append(output_image, output, hdr1)
   