"""
nirspec_extract.py

Program to extract 1D spectra from reduced 2D NIRSPEC spectra.  These reduced
 spectra should be rectified and sky-subtracted, with the dispersion running
 along the x-axis.  The niredux code in the nirspec directory will produce
 such files.
"""

import numpy
from scipy import interpolate,ndimage
from mostools import spectools as st
import spec_simple as ss
import ccdredux as c
import pyfits as p
import pylab as plot
from math import sqrt,log

#-----------------------------------------------------------------------

def clear_all():
   """
   Clears all three figure windows.
   """

   plot.figure(1)
   plot.clf()
   plot.figure(2)
   plot.clf()
   plot.figure(3)
   plot.clf()

#-----------------------------------------------------------------------

def read_nirspec_spec(filename, x1=0, x2=0, y1=0, y2=0, informat='new',
                      verbose=True):
   """
   Reads in a 2D spectrum produced by Matt Auger's niredux code.
   New format (default):
      hdu0 = 2D reduced spectrum
      hdu1 = 2D variance spectrum
   Old format:
      hdu0 = 2D reduced spectrum
      hdu1 = 1D wavelength vector
      hdu2 = 1D sky vector

   Inputs:
      filename - name of file containing reduced spectrum
      x1, etc. - numbers defining the section to trim out of the 2D spectrum,
                 if non-zero.  Defaults = 0
      informat - format of input file (see above). Default = 'new'
      verbose  - set to False to eliminate output
   """

   if(verbose):
      print "Reading file %s." % filename

   try:
      hdulist = p.open(filename)
   except:
      hdulist = p.open(filename,ignore_missing_end=True)

   hdulist.info()

   """ Trim the data if requested """
   xt1,xt2,yt1,yt2 = c.define_trimsec(hdulist[0],x1,x2,y1,y2)
   d = hdulist[0].data[yt1:yt2,xt1:xt2].copy()

   if informat=='old':
      w = hdulist[1].data.copy()
      v = hdulist[2].data.copy()
   else:
      hdr = hdulist[0].header
      v = numpy.median(hdulist[1].data[yt1:yt2,xt1:xt2].copy(),axis=0)
      w = 1.0*numpy.arange(v.size) + 1.0*xt1
      w *= hdr['cd1_1']
      w += hdr['crval1']

   hdulist.close()

   return d,w,v

#-----------------------------------------------------------------------

def nir_extract(filename, outname, x1=0, x2=0, y1=0, y2=0,
                informat='new',outformat='text', apmin=-4., apmax=4.,
                muorder=3, sigorder=3, fitrange=None, weight='gauss', 
                owave=None, do_plot=True, do_subplot=True, 
                stop_if_nan=True):
   """
   Given the calibrated and sky-subtracted 2d spectrum, extracts the 1d 
    spectrum, taking into account things specific to the form of the 2d 
    NIRSPEC spectra.
   """ 

   """ Set up NIRSPEC detector characteristics """
   gain = 4.     # Value in e-/ADU
   rdnoise = 25. # Value in e-

   """ Read the 2d spectrum """
   print ""
   d,w,v = read_nirspec_spec(filename,x1,x2,y1,y2,informat=informat)
   if informat=='old':
      s = v
   else:
      s = numpy.sqrt(v)

   """ Show the 2d spectrum """
   mp,sp = c.sigma_clip(d)
   if do_plot:
      plot.figure(1)
      plot.clf()
      plot.imshow(d,origin='lower',cmap=plot.cm.gray,vmin=mp-sp,vmax=mp+4*sp)

   """ Find the trace """
   if(do_plot):
      # plot.figure(2)
      plot.figure(2, figsize=(15,12))
      plot.clf()
   mu0,sig0 = ss.find_trace(d,apmin=apmin,apmax=apmax,do_plot=do_plot,
                            do_subplot=do_subplot)

   """ Trace the object down the 2d spectrum """
   mupoly,sigpoly = ss.trace_spectrum(d,mu0,sig0,fitrange=fitrange,
                                      muorder=muorder,sigorder=sigorder,
                                      do_plot=do_plot,do_subplot=do_subplot)

   """ Extract the 1d spectrum """
   spec,varspec = ss.extract_spectrum(d,mupoly,sigpoly,apmin=apmin,apmax=apmax,
                                      weight=weight,sky=s,gain=gain,
                                      rdnoise=rdnoise,do_plot=do_plot,
                                      do_subplot=do_subplot)
   if informat=='new':
      varspec = v

   """ Check for NaN values in the returned spectra """
   spec_nan = False
   
   if ((numpy.isnan(spec)).sum() > 0 or (numpy.isnan(varspec)).sum() > 0):
      spec_nan = True

   if (stop_if_nan and spec_nan):
      print ""
      print "WARNING: There is at least one NaN in the extracted spectrum "
      print " or variance spectrum.  Stopping here. "
      print "If you understand the cause of this/these NaN value(s) and want"
      print " to continue, re-run nir_extract with stop_if_nan=False"
      print ""
      return
   elif (spec_nan):
      spec[numpy.isnan(spec)] = 0.0
      varspec[numpy.isnan(varspec)] = 2.0 * numpy.nanmax(varspec)

   """ 
   Linearize the wavelength scale if needed (taken from code by Matt Auger) 
     Old format: explicitly linearize here
     New format: wavelength scale is alread linearized by niredux
   """

   if owave is None:
      if informat=='old':
         w0 = w[0]
         w1 = w[-1]
         owave = numpy.linspace(w0,w1,w.size)
      else:
         owave = w
   
   tmpwave,outspec = ss.resample_spec(w,spec,owave)
   tmpwave,outvar  = ss.resample_spec(w,varspec,owave)
   tmpwave,sky     = ss.resample_spec(w,s,owave)

   """ Plot the linearized spectrum, with rms """
   rms = numpy.sqrt(outvar)
   if do_plot:
      if (do_subplot):
         plot.figure(3)
      else:
         # Re-use fig 4, which showed an earlier version of the extracted 
         # spectrum
         plot.figure(4)
      plot.clf()
      plot.axhline(color='k')
      plot.plot(owave,outspec,linestyle='steps',label='spectrum')
      plot.plot(owave,rms,'r',linestyle='steps',label='RMS')
      plot.xlabel('Wavelength')
      plot.ylabel('Relative Flux Density')
      plot.title('Output Spectrum')
      plot.legend(loc='best')
      plot.xlim([owave[0],owave[-1]])

   """ Write the output spectrum in the requested format """
   if(outformat == 'mwa'):
      st.make_spec(outspec,outvar,owave,outname,clobber=True)
   else:
      ss.save_spectrum(outname,owave,outspec,outvar)

