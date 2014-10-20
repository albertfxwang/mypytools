#!/usr/bin/env python

import numpy as N
import fitsutil
import pyfits
import os, shutil
from pygoods import Ftable

def write_logre(gfsimcat):
   """
   Take a GALFIT simulation catalog, add two additional columns 
   logre_in and logre_out from existing columns re_in and re_out.
   """
   c = Ftable(gfsimcat)
   root, ext = os.path.splitext(gfsimcat)
   newroot = root + '_4kgrid'
   newname = newroot + ext  
   # New file name
   logre_in = N.log10(c.re_in)
   logre_out = N.where(c.re_out>0., N.log10(c.re_out), -99.0)
   # Calculate output log10(re_out)
   colnames = ['logre_in', 'logre_out']
   newarrays = [logre_in, logre_out]
   formats = ['D', 'D']
   shutil.copy(gfsimcat, newname)
   fitsutil.add_columns(newname, colnames, newarrays, formats)
   os.remove('%s.OLD' % newname)
