#!/usr/bin/env python

import numpy as np
try:
   import pyfftw
   from pyfftw.interfaces.numpy_fft import fftshift, ifftshift
   from pyfftw.interfaces.numpy_fft import fftn, ifftn
   _pyfftw = 1
except:
   _pyfftw = 0
import anfft
import pyfits, glob, os
from scipy import fftpack, signal
import multiprocessing

"""
Imports anfft that utilizes FFTW3 to do Fast Fourier Transform really fast.
Hopefully will be faster than any implementation in numpy and scipy!
Convolution function based on the script by Tamas Haraszti:
http://www.rzuser.uni-heidelberg.de/~ge6/Programing/convolution.html
Upate 2014/3/6: calls anfft.rfftn for DFT of *real* input arrays. This gains a
another factor of 3 in speed over the more general function anfft.fftn.
"""

__author__ = "Kuang-Han Huang"
__date__ = "Mar 6, 2014"
__version__ = "1.1"


def padzero2d_k(a, nr, nc):
   """
   a: array to be padded with zero
   nr: number of rows of the output
   nc: number of columns of the output
   do padding so that kernel is centered at the same place as image
   AND also wrap around the kernel
   """
   
   s = list(a.shape)
   z = np.zeros((nr, nc), a.dtype.char)
   if s[0] % 2 == 0:
      rpos = s[0] / 2 + 1 
      rneg = s[0] / 2 - 1
   else:
      rpos = (s[0] + 1) / 2
      rneg = (s[0] - 1) / 2
   if s[1] % 2 == 0:
      cpos = s[1] / 2 + 1
      cneg = s[1] / 2 - 1
   else:
      cpos = (s[1] + 1) / 2
      cneg = (s[1] - 1) / 2
      
   # Wrap around: it's tricky when kernel size is odd
   # Do this by copying quarters of kernel
   z[0:rpos,0:cpos] = a[-rpos:,-cpos:]
   z[-rneg:,0:cpos] = a[0:rneg,-cpos:]
   z[0:rpos,-cneg:] = a[-rpos:,0:cneg]
   z[-rneg:,-cneg:] = a[0:rneg,0:cneg]

   return z
   

def padzero2d_i(a, nr, nc):
   """
   a: array to be padded with zero
   nr: number of rows of the output
   nc: number of columns of the output
   pad zeros to the end of image
   """
   s = list(a.shape)

   index = [slice(None)] * len(s)
   index[0] = slice(0, s[0])
   index[1] = slice(0, s[1])
   s[0] = nr
   s[1] = nc
   z = np.zeros(s, a.dtype.char)
   z[index] = a
   return z	


def hconvolve(image, kernel, pad=True, threads=multiprocessing.cpu_count()):
   """ Not so simple convolution """
   # The size of the image and kernel
   r1, c1 = image.shape
   r2, c2 = kernel.shape

   # Pad zeros of half the size of the kernel
   if pad:
      if _pyfftw:
         # for some reason, pyfftw requires a slightly different padding width
         # if calling rfftn and irfftn; if calling fftn and ifftn, no such change
         # is necessary
         if r2 % 2 == 0:
            r = r1 + r2/2 
         else: 
            r = r1 + (r2 + 1) / 2
         if c2 % 2 == 0:
            c = c1 + c2/2
         else:
            c = c1 + (c2) / 2
      else:
         if r2 % 2 == 0:
            r = r1 + r2/2 
         else: 
            r = r1 + (r2 + 1) / 2
         if c2 % 2 == 0:
            # c = c1 + c2/2 + 1
            c = c1 + c2 / 2
         else:
            # c = c1 + (c2 + 1) / 2
            c = c1 + (c2 / 2)
      
   # end of if pad

   # Does padding:
   # pad zeros on the END of image
   image_p = padzero2d_i(image, r, c)
   #image_p = image.copy()
   # pad zeros on the SIDES of kernel SYMMETRICALLY and then WRAP AROUND
   kernel_p = padzero2d_k(kernel, r, c)

   if _pyfftw:
      f1 = pyfftw.interfaces.numpy_fft.rfftn(image_p, threads=threads)
      f2 = pyfftw.interfaces.numpy_fft.rfftn(kernel_p, threads=threads)
      fftimage = f1 * f2
      if pad:
         conved = pyfftw.interfaces.numpy_fft.irfftn(fftimage, threads=threads)[:r1,:c1].real
      else:
         conved = pyfftw.interfaces.numpy_fft.irfftn(fftimage, threads=threads).real
   else:
      fftimage = anfft.rfftn(image_p) * anfft.rfftn(kernel_p)
      if pad:
         conved = anfft.irfftn(fftimage)[:r1,:c1].real
      else:
         conved = anfft.irfftn(fftimage).real

   return conved

def imhconvolve(imagefile, kernelfile, outputfile, pad=True, overwrite=False):
   """
   Arguments are FITS images instead of just arrays. It is a wrapper
   around the function hconvolve.
   """
   assert os.path.exists(imagefile), "File %s does not exist." % imagefile
   if os.path.exists(outputfile):
      if not overwrite:
         raise NameError, "Image %s already exists; set overwrite=True to overwrite it." % outputfile
      else:
         os.remove(outputfile)
   image = pyfits.getdata(imagefile)
   kernel = pyfits.getdata(kernelfile)
   header = pyfits.getheader(imagefile)
   conved = hconvolve(image,kernel,pad=pad)
   pyfits.append(outputfile, conved, header)

def imfftconvolve(imagefile, kernelfile, outputfile):
   """
   Use scipy.signal.fftconvolve instead of anfft or pyfftw. It is slower, but 
   it gives more robust results...
   """
   assert os.path.exists(imagefile)
   if os.path.exists(outputfile):
      if not overwrite:
         raise NameError, "Image %s already exists; set overwrite=True to overwrite it." % outputfile
      else:
         os.remove(outputfile)
   image = pyfits.getdata(imagefile)
   kernel = pyfits.getdata(kernelfile)
   header = pyfits.getheader(imagefile)
   conved = signal.fftconvolve(image, kernel, mode='same')
   pyfits.append(outputfile, conved, header)

def gausskern(n, sigma=1.0):
   # n needs to be an odd integer
   if n%2 == 0:
      raise ValueError, "n needs to be an odd positive integer"
   k = np.zeros((n,n), dtype='float')
   center = (n+1)/2 - 1
   for i in range(n):
      for j in range(n):
         dist2 = (i - center)**2 + (j - center)**2
         g = np.exp(-dist2 / (2. * sigma**2))
         k[i,j] = g
   return k


def Convolve(image1, image2, MinPad=False, pad=True):
	""" Not so simple convolution """

	#Just for comfort:
	FFt = np.fft.fft2
	iFFt = np.fft.ifft2

	#The size of the images:
	r1,c1 = image1.shape
	r2,c2 = image2.shape

	#MinPad results simpler padding,smaller images:
	if MinPad:
		r = r1+r2
		c = c1+c2
	else:
		#if the Numerical Recipies says so:
		r = 2*max(r1,r2)
		c = 2*max(c1,c2)
	 
	#For nice FFT, we need the power of 2:
	if pad:
		pr2 = int(np.log(r)/np.log(2.0) + 1.0 )
		pc2 = int(np.log(c)/np.log(2.0) + 1.0 )
		rOrig = r
		cOrig = c
		r = 2**pr2
		c = 2**pc2
	#end of if pad
	
	#numpy fft has the padding built in, which can save us some steps
	#here. The thing is the s(hape) parameter:
	fftimage = FFt(image1,s=(r,c)) * FFt(image2[::-1,::-1],s=(r,c))

	#return fftimage.real
	if pad:
		return (iFFt(fftimage))[:rOrig,:cOrig].real
		#return (iFFt(fftimage)).real
	else:
		return (iFFt(fftimage)).real

def fftdeconvolve(image, psf):
   """
   De-convolution by directly dividing the DFT... may not be numerically 
   desirable for many applications. Noise could be an issue.
   Use scipy.fftpack for now; will re-write for anfft later...
   Taken from this post on stackoverflow.com: 
   http://stackoverflow.com/questions/17473917/is-there-a-equivalent-of-scipy-signal-deconvolve-for-2d-arrays
   """
   if not _pyfftw:
      raise NotImplementedError
   image = image.astype('float')
   psf = psf.astype('float')

   # image_fft = fftpack.fftshift(fftpack.fftn(image))
   # psf_fft = fftpack.fftshift(fftpack.fftn(psf))
   image_fft = fftshift(fftn(image))
   psf_fft = fftshift(fftn(psf))
   kernel = fftshift(ifftn(ifftshift(image_fft / psf_fft)))
   return kernel

def RLdeconvolve(imagefile, psffile, deconvfile, maxiter=20, tol=1.e-3):
   """
   To solve the kernel such that deconv convolves with psf = image. Or in 
   other words, image1 is the "reference" image.
   Using the iterative Richardson-Lucy algorithm, following Binney and 
   Merrifield, Galactic Astronomy, Appendix C.
   *** Input images are strongly preferred to have odd number of rows and 
   columns.***
   """
   image = pyfits.getdata(imagefile)
   assert image.min() > 0, "Input image has to be positive!"
   psf = pyfits.getdata(psffile)
   ncols, nlines = image.shape
   ncols_psf, nlines_psf = psf.shape
   if (ncols_psf<ncols) & (nlines_psf<nlines):
      width = (ncols - ncols_psf) / 2
      psf_padded = np.pad(psf, width, mode='constant')
   else:
      psf_padded = psf
   psf_flip = psf_padded[::-1,::-1]
   # if image1.shape != image2.shape:
   #    raise ValueError, "image1 and image2 should have the same dimensions."
   hdr = pyfits.getheader(imagefile)
   assert np.abs(psf.sum() - 1.0) <= 1.e-5, "PSF file is not normalized."
   # enforces the normalization of image1 and image2
   # should I record the normalization constant?
   image = image / image.sum()
   psf = psf / psf.sum()
   # initial guess of kernel
   last_deconv = image.mean() * np.ones(image.shape)
   last_deconv = last_deconv / last_deconv.sum()
   last_image = signal.fftconvolve(last_deconv, psf, mode='same')  # f_i
   niter = 0
   while niter < maxiter:
      niter += 1
      relative_blur = image / last_image
      error_est = signal.fftconvolve(relative_blur, psf_flip, mode='same')
      last_deconv = last_deconv * error_est
      last_image = signal.fftconvolve(last_deconv, psf, mode='same')
      # last_L = last_deconv / last_image * psf_padded
      # new_deconv = signal.fftconvolve(image, last_L, mode='same')
      # last_image = signal.fftconvolve(last_deconv, psf, mode='same')
      if np.max(np.abs((last_image - image) / image)) <= tol:
         print "Converged in %d iterations." % niter
         # last_deconv = new_deconv.copy()
         break
         # last_deconv = last_deconv * signal.fftconvolve(image/last_image, psf,
         #                                                mode='same')
         # last_deconv = last_deconv / last_deconv.sum()
         # last_image = signal.fftconvolve(last_deconv, psf, mode='same')
      # new_deconv = new_deconv / new_deconv.sum()
      # last_image = signal.fftconvolve(new_deconv, psf, mode='same')
      # print last_image.max()
      # if np.abs((new_deconv - last_deconv) / last_deconv).max() <= tol:
      
      # last_deconv = new_deconv.copy()

   if niter == maxiter:
      print "Max iterations (%d) reached." % (maxiter)
      print "Last iteration has mean deviation of %f" % (np.max(np.abs((last_image - image) / image)))
   if os.path.exists(deconvfile):
      os.remove(deconvfile)
   pyfits.append(deconvfile, last_deconv, hdr)


