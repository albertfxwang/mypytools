#!/usr/bin/env python

import numpy as np
import os, sys
from pygoods import sextractor
import pyfits

def make_tfitcat(sexcat, flagimage=None, segimage=None, snfloor=0.0,
                 idlist=[], outname=""):
   """
   From SExtractor catalog, make a catalog that can be used for TFIT, with the
   following columns:
   NUMBER
   X_IMAGE
   Y_IMAGE
   XMIN_IMAGE
   YMIN_IMAGE
   XMAX_IMAGE
   YMAX_IMAGE
   BACKGROUND
   FLUX_ISO
   """
   # Optional input flag image and then make an output catalog for TFIT that
   # only contains the objects not flagged out
   c = sextractor(sexcat)
   if flagimage == None:
      all_objects = c.number
   else:
      # Now determine which objects are not completely flagged out
      flagimg = pyfits.getdata(flagimage)
      segimg = pyfits.getdata(segimage)
      segimg_flagged = np.where(flagimg == 0, segimg, 0)
      # zero-out the pixels flagged out (with value >0 in the flag image)
      all_objects = np.unique(segimg_flagged.ravel())
      all_objects = all_objects[all_objects > 0]
      # remove 0 from the array
      # Now apply a S/N floor to get rid of very low S/N detections (in ISO aperture)
      if snfloor > 0:
         objid_sncrit = c.number[(c.flux_iso/c.fluxerr_iso)>=snfloor]
         all_objects = np.intersect1d(all_objects, objid_sncrit)
   if len(idlist):
      all_objects = np.intersect1d(all_objects, idlist)
   
   print "number of all objects", len(all_objects)

   root = os.path.splitext(sexcat)[0]
   if len(outname) == 0:
      outname = root + '_tfit.cat'
   f = open(outname, 'wb')
   f.write('# 1 NUMBER\n')
   f.write('# 2 X_IMAGE\n')
   f.write('# 3 Y_IMAGE\n')
   f.write('# 4 XMIN_IMAGE\n')
   f.write('# 5 YMIN_IMAGE\n')
   f.write('# 6 XMAX_IMAGE\n')
   f.write('# 7 YMAX_IMAGE\n')
   f.write('# 8 BACKGROUND\n')
   f.write('# 9 FLUX_ISO\n')
   f.write('# 10 FLUXERR_ISO\n')
   all_index = np.arange(len(c))
   for i in range(len(all_objects)):
      j = all_index[c.number == all_objects[i]][0]
      f.write('%d %f %f ' % (c.number[j], c.x_image[j], c.y_image[j]))
      f.write('%d %d ' % (c.xmin_image[j], c.ymin_image[j]))
      f.write('%d %d ' % (c.xmax_image[j], c.ymax_image[j]))
      f.write('%f %f ' % (c.background[j], c.flux_iso[j]))
      f.write('%f ' % c.fluxerr_iso[j])
      f.write('\n')

   f.close()

if __name__ == "__main__":
   sexcat = sys.argv[1]
   make_tfitcat(sexcat)
   