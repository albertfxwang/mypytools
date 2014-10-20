#!/usr/bin/env python

import readgfheader as rh
import numpy as N
from pygoods import sextractor, Ftable
import os
import pyfits
import glob, os

def merge_galfit(outputfile, outputdir):
   # If have not merged everything yet... merge now
   rh.write_gf_output(outputdir,'obj*out.fits',outputfile)
   # this will be a SExtractor-format catalog, with FILENAME as the first column...

def add_crashes(outputfile, inputdir):
   """
   Add entries for objects that GALFIT crashed on. Look for GALFIT input files
   with file names obj[xxx], where [xxx] is the object ID. Use default 
   non-measurement values for these objects that GALFIT failed.
   """
   inputfiles = glob.glob(os.path.join(inputdir, 'obj*[0-9]'))
   # Assume that the GALFIT input files have file names obj[ID]
   c0 = sextractor(outputfile)
   newlines = []
   for f in inputfiles:
      f = os.path.split(f)[-1]
      objid = int(f.split('_')[0][3:])
      fout = 'obj%d_out.fits' % objid
      if fout not in c0.filename:
         l = ['obj%d_crashed' % objid, -1., -1., -1., -1., 99.0, 99.0,
            -1., -1., -1., -1., -1., -1., 999.0, 999.0, -1., 1, 1, 1, 1, 1, 1, 
            1, 0]
         newlines += [' '.join(str(x) for x in l) + '\n']
   f = open(outputfile)
   oldlines = f.readlines()
   f.close()
   f = open(outputfile, 'wb')
   lines = oldlines + newlines
   for i in range(len(lines)):
      f.write(lines[i])
   f.close()

def add_objid(outputfile):
   """
   Add a column OBJID to the catalog.
   """
   c = sextractor(outputfile)
   ncols = len(c._colnames)
   if 'objid' not in c._colnames:
      c._header = c._header + '# %d OBJID\n' % (ncols + 1)
      objid = []
      for i in range(len(c)):
         num = c.filename[i].split('_')[0][3:]
         num = int(num)
         objid += [num]
      objid = N.array(objid)
      f = open(outputfile, 'wb')
      f.write(c._header)
      for i in range(len(c)):
         for j in range(len(c._colnames)):
            f.write('%s ' % c._colentries[i][j])
         f.write('%d ' % objid[i])
         f.write('\n')
      f.close()

def add_header_prefix(catalog, prefix):
   # Add a prefix (e.g., filter name) to all columns in the catalog
   c = sextractor(catalog)
   nobj = len(c)
   colnames = c._colnames
   for i in range(len(colnames)):
      colnames[i] = '_'.join([prefix, colnames[i]])
   os.system('cp %s %s.copy' % (catalog, catalog))
   f = open(catalog, 'wb')
   for i in range(len(colnames)):
      f.write('# %d %s\n' % (i+1, colnames[i].upper()))
   for i in range(nobj):
      for j in range(len(c._d)):
         f.write('%s ' % c._colentries[i][j])
      f.write('\n')
   f.close()

def step2():
   # now re-format this and make object IDs as the first column
   c = sextractor('udrops_goodss_h13_130305_galfit.cat')
   objectid = N.zeros(len(c),'int')
   for i in range(len(c)):
      root = os.path.splitext(c.filename[i])[0]  # in the format of objXXXX, where XXXX is the object ID
      objectid[i] = int(root[3:-4])  # the object ID
   header = c._header.replace('FILENAME','OBJECTID')
   # Now write to the catalog
   f = open('udrops_goodss_h13_130305_galfit2.cat','wb')
   f.write(header)
   for i in range(len(c)):
      f.write('%8d ' % objectid[i])
      for j in range(1,len(c._colnames)):
         f.write('%s ' % c._colentries[i][j])
      f.write('\n')
   f.close()

def step3(udropscat,galfitcat,outputcat):
   # Now merge the GALFIT results with original U-dropout catalog
   # Both catalogs should be FITS tables
   # Output catalog should 
   if os.path.exists(outputcat):
      os.remove(outputcat)
   c1 = Ftable(udropscat)  # the U-dropout catalog with the same column as the GOODS-S MW catalog
   c2 = Ftable(galfitcat)  # the merged GALFIT catalog
   colnames = c1.Columns + c2.Columns[1:]  # bypass c2.objectid (redundant)
   formats = c1.d.formats + c2.d.formats[1:]
   # Find matching indices
   #i1 = N.arange(len(c1.d),'int')  # indices in c1
   i2 = N.ones(len(c1.d),'int') * (-1) # to store indices in c2 corresponding to i1
   # e.g., c2.objectid[i2[0]] = c1.id_1[0]
   for i1 in range(len(c1.d)):
      j = N.arange(len(c2.d))[c2.objectid==c1.id_1[i1]]
      if len(j) > 0:
         i2[i1] = j[0]
   # Now initialize the arrays
   arrays = {}
   # First collect all columns from c1
   for col in c1.Columns:
      arrays[col] = getattr(c1,col)
   for i in range(len(c2.Columns)):
      if c2.Columns[i] == 'objectid': 
         pass
      elif c2.d.formats[i] == 'D':
         arrays[c2.Columns[i]] = N.ones(len(c1.d),'float')*(-99.0)
      else:
         arrays[c2.Columns[i]] = N.ones(len(c1.d),'int')*(-1)
   print arrays.keys()
   # Now copy the arrays in c2
   for i1 in range(len(c1.d)):
      if i2[i1] >= 0:
         for col in c2.Columns[1:]:
            arrays[col][i1] = getattr(c2,col)[i2[i1]]
   # Now convert arrays into a list
   arrays2 = []
   allcolumns = []
   for col in colnames:
      arrays2 += [arrays[col]]
   # Now construct pyfits.ColDefs
   allcolumns = []
   for i in range(len(colnames)):
      allcolumns += [pyfits.Column(name=colnames[i],format=formats[i],array=arrays2[i])]
   allcols = pyfits.ColDefs(allcolumns)
   tbhdu = pyfits.new_table(allcols)
   tbhdu.writeto(outputcat)
