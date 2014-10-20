#!/usr/bin/python

######## Generate the bandpass objects and calculate zeropoints in pysynphot#########
### K-H. Huang
### TO BE RUN ON SCIENCE LINUX CLUSTER
### -Version 1.1bv: includes zeropoint fluxes calculated by pysynphot by two ways:
### 1. Renormalize a blackbody spectrum to the zeropoint magnitude and calculate the flux through a given
###    bandpass using pysynphot.integrate directly (method 1)
### 2. Use pysynphot to implement the equivalent of "effstim" function in synphot (method 2) -- 09/23/2008
### -Version 1.2bv: rewrite the calculation as function calc_zeroflux -- 09/23/2008
### -Version 1.3: calculate zeropoint fluxes directly by converting ABmag to micro-Jansky -- 09/28/08
### -Version 1.3.1: changed bandpass into a dictionary, and will be referenced by bandpass name


import os
import numpy as np
import pysynphot as S
import matplotlib.pyplot as plt
#from pyraf import iraf
#from iraf import stsdas,hst_calib,ttools,synphot

rootdir = os.getenv('mypytools')
filters = {}

filters['ukpno'] = S.FileBandpass(rootdir+'/sedtools/bandpass/CANDELS_FILTERS/MOSAIC/U_kpno_mosaic.bp')   # "Mosaic U band in TFIT catalog
filters['uctio'] = S.FileBandpass(rootdir+'/sedtools/bandpass/CANDELS_FILTERS/CTIO/U_ctio.bp')
filters['uvimos'] = S.FileBandpass(rootdir+'/sedtools/bandpass/U_vimos.bp')
filters['f435w'] = S.ObsBandpass('acs,wfc2,f435w')
filters['f475w'] = S.ObsBandpass('acs,wfc2,f475w')
filters['f555w'] = S.ObsBandpass('acs,wfc2,f555w')
filters['f606w'] = S.ObsBandpass('acs,wfc2,f606w')
filters['f625w'] = S.ObsBandpass('acs,wfc2,f625w')
filters['f775w'] = S.ObsBandpass('acs,wfc2,f775w')
filters['f814w'] = S.ObsBandpass('acs,wfc2,f814w')
filters['f850lp'] = S.ObsBandpass('acs,wfc2,f850lp')
filters['f098m'] = S.ObsBandpass('wfc3,ir,f098m')
filters['f105w'] = S.ObsBandpass('wfc3,ir,f105w')
filters['f110w'] = S.ObsBandpass('wfc3,ir,f110w')
filters['f125w'] = S.ObsBandpass('wfc3,ir,f125w')
filters['f140w'] = S.ObsBandpass('wfc3,ir,f140w')
filters['f160w'] = S.ObsBandpass('wfc3,ir,f160w')
filters['irac1'] = S.FileBandpass(rootdir+'/sedtools/bandpass/CANDELS_FILTERS/IRAC/irac_ch1.dat')     # IRAC ch1 (vice versa)
filters['irac2'] = S.FileBandpass(rootdir+'/sedtools/bandpass/CANDELS_FILTERS/IRAC/irac_ch2.dat')
filters['irac3'] = S.FileBandpass(rootdir+'/sedtools/bandpass/CANDELS_FILTERS/IRAC/irac_ch3.dat')
filters['irac4'] = S.FileBandpass(rootdir+'/sedtools/bandpass/CANDELS_FILTERS/IRAC/irac_ch4.dat')

for k in filters.keys():
   filters[k].name = k.lower()
   
#obsband = ['U_mosaic','f435w','f606w','f775w','f850lp','J_isaac','H_isaac','Ks_isaac',\
#           'irac1.fits','irac2.fits','irac3.fits','irac4.fits']
#zeropoint = [31.251, 25.65288, 26.49341, 25.64053, 24.84315, 26.0, 26.0, 26.0,\
#             22.416, 22.195, 20.603, 21.781]
# zeropoint magnitudes in ABmag
#bb = S.BlackBody(10000)

#def zeroflux(zparray):
    # convert an array of zeropoint magnitudes into micro-Jansky
#    zfarray = M.ABmag_to_uJy(zparray)
#    return zfarray

def write_eazy_filter_file(output='FILTERS.EAZY.res'):
   # write my own filter file for EAZY
   f = open(output, 'wb')
   n = 0
   for b in filters.keys():
      n += 1
      print b
      band = filters[b]
      f.write("# %d %s F%d \n" % (len(band.wave), band.name, n))
      for i in range(len(band.wave)):
         w = band.wave[i]
         t = band.sample(w)
         f.write('%d  %18.6g  %18.6g \n' % ((i+1), w, t))
   f.close()

def resampleFilter(bp, output, w0, w1, dw=20.):
   """
   Re-sample a filter curve between w0 and w1 (angstroms) with step dw.
   Usually we use this to re-sample the filter curve to lower resolutions 
   (e.g., to be used by Le Phare).
   First argument, bp, should be a pysynphot filter instance.
   """
   nwaves = (w1 - w0) // dw + 1
   nwaves = int(nwaves)
   print "Number of resampled wavelength steps for %s: " % bp.name, nwaves
   if nwaves > 350:
      print "Warning: more than 350 wavelength steps (%d)!!" % nwaves
   w1 = w0 + dw * (nwaves - 1)  
   # adjust w1 to make the wavelength range an integer multiple of dw
   wave = np.linspace(w0, w1, nwaves)
   throughput = bp.sample(wave)
   f = open(output, 'wb')
   f.write('#  wavelength (A)   throughput  \n')
   # Pad the first and last wavelength with an entry with 0 throughput
   f.write('%12.3f  0.\n' % (wave[0] - dw))
   for i in range(nwaves):
      f.write('%12.3f  %12.4e\n' % (wave[i], throughput[i]))
   f.write('%12.3f  0.\n' % (wave[-1] + dw))
   f.close()

def resampleDefaultFilters(nw=3, nwaves=300):
   """
   Re-sample all default filters to low resolution.
   """
   for b in filters:
      bp = filters[b]
      w0 = bp.pivot() - nw * bp.equivwidth()
      w1 = bp.pivot() + nw * bp.equivwidth()
      dw = (w1 - w0) / nwaves
      dw = np.floor(dw)
      dw = np.maximum(dw, 1.)
      output = '%s_lr.bp' % bp.name
      resampleFilter(bp, output, w0, w1, dw=dw)


def plotFilter(bp, **plot_kwargs):
   plt.figure()
   plt.plot(bp.wave, bp.throughput, **plot_kwargs)
   w0 = bp.pivot() - bp.equivwidth() * 2
   w1 = bp.pivot() + bp.equivwidth() * 2
   plt.xlim(w0, w1)
   plt.xlabel('Wavelength (A)')
   plt.ylabel('Throughput')
   if hasattr(bp, 'name'):
      plt.title(bp.name.upper())


