#!/usr/bin/env python

import numpy as np
from pygoods import sextractor, Ftable
import pyfits

bands_goodsn = ['kpno_u', 'acs_f435w', 'acs_f606w', 'acs_f775w', 'acs_f850lp',
                'wfc3_f105w', 'wfc3_f125w', 'wfc3_f140w', 'wfc3_f160w',
                'cfht_ks', 'irac_ch1', 'irac_ch2', 'irac_ch3', 'irac_ch4']

bands_goodss = ['ctio_u', 'vimos_u', 'acs_f435w', 'acs_f606w', 'acs_f775w', 
                'acs_f814w', 'acs_f850lp', 'wfc3_f098m', 'wfc3_f105w',
                'wfc3_f125w', 'wfc3_f140w', 'wfc3_f160w', 'isaac_ks',
                'hawki_ks', 'irac_ch1', 'irac_ch2', 'irac_ch3', 'irac_ch4']

def ujy_to_abmag(flux):
   return 23.9 - 2.5 * np.log10(flux)

def ujy_to_abmag_1sig(flux, fluxerr):
   mag = np.where(flux/fluxerr > 1., ujy_to_abmag(flux), ujy_to_abmag(fluxerr))
   mag = np.where(fluxerr <= 0., 99.0, mag)
   return mag

class CANDELSMWcat(object):
   def __init__(self, filename, mode='readonly'):
      self.h = pyfits.open(filename, mode=mode)
      self.mode = mode
      self.Columns = []
      for i in range(len(self.h[1].data.columns)):
         self.Columns += [self.h[1].data.columns[i]]

   def add_mag_1sig(self, bands):
      # Add the 1-sigma magnitude columns to the table
      for b in bands:
         print b
         flux = self.h[1].data['%s_flux' % b.lower()]
         fluxerr = self.h[1].data['%s_fluxerr' % b.lower()]
         mag_1sig = ujy_to_abmag_1sig(flux, fluxerr)
         column = pyfits.Column(name='%s_mag_1sig' % b.lower(),
                                array=mag_1sig, 
                                format='D')
         self.Columns += [column]
      # Add the new columns
      # self.h[1].data.columns = self.h[1].data.columns + newcolumns
   
   def writeto(self, newfilename, clobber=True):
      tbhdu = pyfits.new_table(self.Columns)
      tbhdu.writeto(newfilename, clobber=clobber)

