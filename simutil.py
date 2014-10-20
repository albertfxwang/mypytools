#!/usr/bin/env python

# Utilities for manipulating simulation files done by K-H Huang Oct 2009
# Changed default z range -- 10/14/09 K-H Huang

from numpy import *
from pygoods import *
import os, glob
import KPDFadaptnumpy as KPDF
import pysynphot as S

bands = ['b','v','i','z']

def bdrop(bmag,vmag,zmag):
    """Usage: bdrop(bmag,vmag,zmag)
       returns True if input mag combination satisfies 
       B-dropout criteria from K-S. Lee et al. 2006
    """
    c1 = (bmag-vmag)>=1.2+1.4*(vmag-zmag)
    c2 = (bmag-vmag)>=1.2
    c3 = (vmag-zmag)<=1.2
    return (c1 * c2 * c3) # equiv to (c1 and c2 and c3)

def vdrop(vmag,imag,zmag):
    """Usage: vdrop(vmag,imag,zmag)
       returns True if mag combination satisfies V-dropout
       criteria in K-S. Lee et al. 2006
    """
    c1 = (vmag-imag)>1.5+0.9*(imag-zmag)
    c2 = (vmag-imag)>2.0
    c3 = c1 + c2 # c3>0 if any of c1 or c2 is true; equiv to c1 or c2
    c4 = (vmag-imag)>=1.2
    c5 = (imag-zmag)<=1.3
    return (c3 * c4 * c5) # equiv to (c3 and c4 and c5)

    
class dropouts:
    def __init__(self,Mlim,rband,obands,mags,zarray,ebmv,completeness):
        self.Mmin = Mlim[0]
        self.Mmax = Mlim[1]
        self.rband = rband
        for i in range(len(obands)):
            setattr(self,obands[i],mags[i])
        self.zarray = zarray
        self.zrange = arange(2.0,6.02,0.02) # default zrange
        if len(zarray):
            self.zdist = calczdist(zarray,self.zrange)
        else:
            self.zdist = zeros(len(self.zrange))
        self.ebmv = ebmv
        self.completeness = completeness
    def newzrange(self,zrange):
        self.zrange = zrange
        self.zdist = calczdist(self.zarray,self.zrange)

def plot_bdropcrit():
    vmzr = arange(-2.,5.)
    bmvr = arange(-2.,10.)
    bmv1 = 1.2 + 1.4 * vmzr
    line1 = plot(vmzr,bmv1,color='black')
    bmv2 = ones(len(vmzr))*1.2
    line2 = plot(vmzr,bmv2,color='black')
    vmz3 = ones(len(bmvr))*1.2
    line3 = plot(vmz3,bmvr,color='black')
        

def add_lyman_alpha(sp, obslam, eqwidth):
   # adds Lyman alpha line flux to the SED at a single wavelength
   # use equivalent width information
   # since equivalent width is observed, it is affected by dust extinction
   # and also is in observed wavelength. Therefore one needs to add Lyman alpha
   # only after applying redshift, IGM, and dust extinction
   # both obslam and eqwidth are in angstrom
   wave = sp.wave.copy()
   # first resample the SED at the observed wavelength of Lya
   if obslam not in wave:
      wave = concatenate([wave, [obslam]])
      wave = sort(wave)
   if eqwidth >= 0.:
      # for Lya in emission: add the line flux to the flam @ a single wavelength
      sp = sp.resample(wave)
      f0 = compress(sp.wave == obslam, sp.flux)[0]  # the continuum level in units of flam
      # make a second SED containing zero flux except at obslam
      i = searchsorted(wave, obslam)
      wave2 = wave[i-5:i+5]
      flux2 = zeros(len(wave2))
      # calculate line flux of Lya
      df = eqwidth * f0 * 2 / (wave2[6]-wave2[4])
      flux2[5] = df
      sp2 = S.ArraySpectrum(wave=wave2, flux=flux2, waveunits='angstrom',\
         fluxunits='flam')
      sp2 = sp2.resample(wave)
      sp3 = sp + sp2
   else:
      # Lya in absorption: zero out f within the appropriate range
      # need to be careful about the range of wavelength that we set flux to zero
      # in order for trapezoid integration to get the correct absorbed flux
      eqwidth = abs(eqwidth)
      sp = sp.resample(wave)
      wave2 = wave.copy()
      lam1 = obslam - eqwidth/2.  # real edges of the absorption trough
      lam2 = obslam + eqwidth/2.
      i1 = searchsorted(sp.wave, lam1)
      i2 = searchsorted(sp.wave, lam2)
      if lam1 not in wave2:   # add edges in the wavelength array
         wave2 = concatenate([wave2, [lam1]])
      if lam2 not in wave2:
         wave2 = concatenate([wave2, [lam2]])
      lam_0 = lam1 - 1.e-6
      lam_1 = lam2 + 1.e-6
      if lam_0 not in wave2:  # add edges in the wavelength array
         wave2 = concatenate([wave2, [lam_0]])
      if lam_1 not in wave2:
         wave2 = concatenate([wave2, [lam_1]])
      wave2 = sort(wave2)
      sp2 = sp.resample(wave2)  # add in wavelength entries corresponding to the absorption edge
      flux2 = sp2.flux.copy()
      i11 = searchsorted(sp2.wave, lam1)
      i22 = searchsorted(sp2.wave, lam2)
      flux2[i11:i22+1] = 0.
      sp3 = S.ArraySpectrum(wave=wave2, flux=flux2, waveunits='angstrom',\
         fluxunits='flam')
   return sp3


def test_lya(sp, obsband, obslam, eqwidth):
   sp3 = add_lyman_alpha(sp, obslam, eqwidth)
   f1 = (sp * obsband).integrate(fluxunits='flam')
   f2 = (sp3 * obsband).integrate(fluxunits='flam')
   ew = (f2 - f1) / sp.sample(obslam)
   print "input equiv width; calculated equiv width:"
   print eqwidth, ew


      

def make_lyaew_dist(xr, fdist):
   c = sextractor('Lyaew_list.dat')
   h = KPDF.UPDFOptimumBandwidth(c.lya_ew)
   yr = KPDF.UPDFEpanechnikov(c.lya_ew, xr, h)
   yr = yr / max(yr)
   f = open(fdist, 'w')
   f.write('# 1 Lya_EW\n')
   f.write('# 2 PDF\n')
   for i in range(len(xr)):
      f.write('%f  %f\n' % (xr[i], yr[i]))
   f.close()
   return xr, yr
