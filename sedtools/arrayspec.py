# Subclasses for taking in spectra and spectral elements from arrays
### Last updated by Kuang-Han Huang on 2/22/08 ###

import math
import pysynphot as S
from pysynphot.spectrum import *
from pysynphot import units
from pysynphot import spectrum
#from pysynphot import refs
#from pysynphot.spectrum import FlatSpectrum
#from numprint import *
from numpy import *
import os

pytool_path = os.getenv('mypytools')

try:
   defwaveset = refs._default_waveset
except:
   defwaveset = spectrum._computeDefaultWaveset()


class ArraySourceSpectrum(TabularSourceSpectrum):
    '''Class for a source spectrum that is read in from a wavelength and flux array
    '''
    def __init__(self, wave, flux, waveunits="angstrom", fluxunits="flam"):
        self._wavetable = wave
        self._fluxtable = flux
        self.waveunits = units.Units(waveunits)
        self.fluxunits = units.Units(fluxunits)
        self.filename = None
        self.ToInternal()


class ArraySpectralElement(TabularSpectralElement):
    def __init__(self, wave, throughput, waveunits="angstrom"):
        self.name = None
        self.wavetable  = wave
        self.throughputtable = throughput
        self.waveunints = units.Units(waveunits)
        self.throughputunits = 'none'
    

def ABmag(spec,bandpass):
    if bandpass == None:
        bandpass = S.ArrayBandpass(defwaveset,ones(len(defwaveset)))
    stflux = 10.**(48.6/-2.5)
    abu = S.ArraySpectrum(defwaveset,ones(len(defwaveset))*10.**(48.6/-2.5),
        fluxunits='fnu')
    #abu = spectrum.FlatSpectrum(10.**(48.6/-2.5),fluxunits='fnu')
    if bandpass.wave[-1] > spec.wave[-1]: 
        # if bandpass wavelenth range is longer than spectrum    
        n = searchsorted(bandpass.wave,spec.wave[-1])
        merge = concatenate((spec.wave,bandpass.wave[n:]))
        merge = sort(unique(merge))
        flux_standard = S.ArraySpectrum(merge,ones(len(merge))*stflux,
            fluxunits='fnu')
    else:
        flux_standard = S.ArraySpectrum(spec.wave,ones(len(spec.wave))*stflux,
            fluxunits='fnu') 
        # Sets the wavelength array of abu
#   print spec.flux.min(),spec.flux.max()
#   print spec.flux
#   l = format("%10.1f %10.2e",spec.wave,spec.flux)
#   print l
#   print (spec*bandpass).integrate()
#   print (flux_standard*bandpass).integrate()	
    numerator = (spec*bandpass).integrate()
    denominator = (flux_standard*bandpass).integrate()
    ratio = numerator/denominator
    if ratio <= 0.:
        abmag = 999.
    else:
        abmag = -2.5*math.log10(ratio)
    return abmag
    
def Vegamag(spec,band):       # added by Kuang-Han Huang
    vega = S.FileSpectrum(pytool_path+'/sedtools/alpha_lyr_stis_003.fits')
    # make sure that this vega spectrum is in the same directory
    vegaflux = (vega*band).integrate()
    # vegaflux corresponds to vegamag 0 in any band
    flux = (spec*band).integrate()
    vegamag = -2.5*log10(flux/vegaflux)
    return vegamag

def monoband(lambda_central,width=100):
    """construct a uniform bandpass @ lambda_central & width"""
    # can just use pysynphot.Box...
    wlo = lambda_central - width + 1
    w0 = lambda_central - width/2. + 1
    w1 = lambda_central + width/2. + 1
    whi = lambda_central + width
    waveset = arange(wlo,whi)
    fluxset1 = zeros(len(arange(wlo,w0)))
    fluxset2 = ones(len(arange(w0,w1)))
    fluxset3 = zeros(len(arange(w1,whi)))
    fluxset = concatenate([fluxset1,fluxset2,fluxset3])
    uni = S.ArrayBandpass(waveset,fluxset)
    #uni = uni.compute(waveset)
    return uni
   
