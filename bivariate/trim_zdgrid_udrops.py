#!/usr/bin/env python

"""
To trim fat from zdgrid_udrops*.p
"""

import os, cPickle

# trimlist = ['acs_f435w', 'acs_f606w', 'wfc3_f105w', 'wfc3_f098m', 'wfc3_f160w',
#             'vimos_u']
trimlist = ['acs_f435w', 'acs_f606w', 'wfc3_f105w', 'wfc3_f098m', 'wfc3_f160w']

class zdgrid_sup(object):
   def __init__(self, field):
      self.field = field


def trim_zdgrid_udrops(filename, field, write=True):
   f = open(filename)
   root = os.path.splitext(filename)[0]
   supname = root + '_sup.p'
   zdgrid = cPickle.load(f)
   if write:
      zdgrid2 = zdgrid_sup(field)
   f.close()
   for t in trimlist:
      for k in zdgrid.__dict__.keys():
         if k.startswith(t):
            if write:
               setattr(zdgrid2, k, zdgrid.__dict__[k])
            delattr(zdgrid, k)
   if hasattr(zdgrid, 'sn_lolim_crit'):
      delattr(zdgrid, 'sn_lolim_crit')
   if hasattr(zdgrid, 'sn_hilim_crit'):
      delattr(zdgrid, 'sn_hilim_crit')
   if hasattr(zdgrid, 'lcc'):
      delattr(zdgrid, 'lcc')
   if write:
      f = open(supname, 'wb')
      cPickle.dump(zdgrid2, f, 2)
      f.close()
   f = open(filename, 'wb')
   cPickle.dump(zdgrid, f, 2)
   f.close()
