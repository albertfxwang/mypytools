#!/usr/bin/env python

from numpy import *
import pysynphot as S
from pysynphot import extinction
from pysynphot.extinction import _ExtinctionLaw,_buildDefaultWaveset
from pysynphot import spectrum
#from pysynphot.spectrum import _computeDefaultWaveset
from pysynphot import units

_waveset = _buildDefaultWaveset()
_waveset = 1.0/_waveset

def computeTrans_C00():
   """
   Calculate the attenuation curve k(lambda) for unit E(B-V).
   """
   wave_um = _waveset
   wave0 = compress((wave_um<0.12),wave_um)
   wave1 = compress((wave_um>=0.12)*(wave_um<0.63),wave_um)
   wave2 = compress((wave_um>=0.63)*(wave_um<2.2),wave_um)
   wave3 = compress(wave_um>=2.2,wave_um)
   
   k1 = 2.659 * (-2.156 + 1.509/wave1 - 0.198/wave1**2 + 0.011/wave1**3) + 4.05
   k2 = 2.659 * (-1.857 + 1.040/wave2) + 4.05
   k3 = ones(len(wave3),'float')*k2[-1]
   k0 = ones(len(wave0),'float')*k1[0]
   k = concatenate([k0,k1,k2,k3])
   return k

class Calzetti00(_ExtinctionLaw):
   citation = 'Calzetti et al. 2000 (ApJ, 533, 682)'
   name = 'CALZETTI00'
   def __init__(self, ebmv):
      """
      Calculate the dust attenuation following the Calzetti+2000 relation,
      parameterized by the E(B-V) of the stellar continuum. The relation in 
      E(B-V) between the stellar and nebular attenuation is roughly
      E(B-V)_stellar = (0.44 +/- 0.03) E(B-V)_neb
      as given in Calzetti+2000.
      """
      self._wavetable = _waveset.copy() * 10000.
      # self.A_lam = 0.4 * ebmv * computeTrans_C00()
      self.transparencytable = 10.**(-1 * self.A_lam)

# class Extinction(spectrum.ArraySpectralElement):
#    """extinction = Extinction(extinction in magnitudes,
#    'gal1|smc|lmc reddening laws)"""
#    def __init__(self, extval):
#       ''' Extinction mimics as a spectral element.
#       '''
#       law = Calzetti00(extval)
#       #self._wavetable = law._wavetable
#       #self._throughputtable = law.transparencytable
#       self.extinction = S.ArrayBandpass(law._wavetable,law.transparencytable)
#       self.name=law.name
#       self.citation=law.citation
#       self.waveunits=units.Units('angstrom')
#       self.isAnalytic=False
#       self.warnings={}
#       self.extval = extval
#    #def __call__(self):
      
