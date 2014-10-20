#!/usr/bin/env python

from numpy import *
import pyfits
import os, glob
import string

def readgfheader(img):
   """assumes that extension 2 is the model image"""
   delim = '+/-'
   h = pyfits.open(img)
   h2 = h[2]
   hxcent1 = h2.header['1_XC']
   xcent1,xcent1err = hxcent1.split(delim)
   if '*' in xcent1: xflag=1
   else: xflag=0
   hycent1 = h2.header['1_YC']
   ycent1,ycent1err = hycent1.split(delim)
   if '*' in ycent1: yflag=1
   else: yflag=0
   hmag1 = h2.header['1_MAG']
   mag1,mag1err = hmag1.split(delim)
   if '*' in mag1: magflag=1
   else: magflag=0
   hre1 = h2.header['1_RE']
   re1,re1err = hre1.split(delim)
   if '*' in re1: reflag=1
   else: reflag=0
   hn1 = h2.header['1_N']
   n1,n1err = hn1.split(delim)
   if '*' in n1: nflag=1
   else: nflag=0
   har1 = h2.header['1_AR']
   ar1,ar1err = har1.split(delim)
   if '*' in ar1: arflag=1
   else: arflag=0
   hpa1 = h2.header['1_PA']
   pa1,pa1err = hpa1.split(delim)
   if '*' in pa1: paflag=1
   else: paflag=0
   hchisqnu = h2.header['CHI2NU']
   chisqnu = hchisqnu  # this is float
   try:
      ncomp = int(h2.header['NCOMP'])
   except:
      ncomp = -1
   results = [xcent1,xcent1err,ycent1,ycent1err,mag1,mag1err,re1,re1err,n1,n1err,ar1,ar1err,pa1,pa1err]
   #print results
   results = map(string.replace,results,['*']*14,['']*14)  # replace '*' by ''
   results = results + [chisqnu]
   results = map(float,results) # convert to float
   results = results + [xflag,yflag,magflag,reflag,nflag,arflag,paflag,ncomp]
   return results


def write_gf_output(dir, fpattern, output_filename):
   """
   Read the GALFIT output from all the images in the directory dir.
   The model image is in the 2nd extension.
   Then write the catalog to output_filename.
   fpattern --- the file pattern used in glob to search for GALFIT output images; should include .fits extension
   """
   print "fpattern:", fpattern
   images = glob.glob(dir+'/'+fpattern)   # assume that all the FITS images in dir are
                                       	# GALFIT output images
   print "%d images found; now write merged GALFIT catalog." % len(images)
   nobj = len(images)
   xobj = zeros(nobj)
   xobj_err = zeros(nobj)
   yobj = zeros(nobj)
   yobj_err = zeros(nobj)
   mag = zeros(nobj)
   mag_err = zeros(nobj)
   Reff = zeros(nobj)
   Reff_err = zeros(nobj)
   sersicn = zeros(nobj)
   sersicn_err = zeros(nobj)
   ar = zeros(nobj)
   ar_err = zeros(nobj)
   pa = zeros(nobj)
   pa_err = zeros(nobj)
   chi2nu = zeros(nobj)
   xflag = zeros(nobj,'int')
   yflag = zeros(nobj,'int')
   magflag = zeros(nobj,'int')
   reflag = zeros(nobj,'int')
   nflag = zeros(nobj,'int')
   arflag = zeros(nobj,'int')
   paflag = zeros(nobj,'int')
   ncomps = zeros(nobj, 'int')

   for i in range(nobj):
      result = readgfheader(images[i])
      xobj[i] = result[0]
      xobj_err[i] = result[1]
      yobj[i] = result[2]
      yobj_err[i] = result[3]
      mag[i] = result[4]
      mag_err[i] = result[5]
      Reff[i] = result[6]
      Reff_err[i] = result[7]
      sersicn[i] = result[8]
      sersicn_err[i] = result[9]
      ar[i] = result[10]
      ar_err[i] = result[11]
      pa[i] = result[12]
      pa_err[i] = result[13]
      chi2nu[i] = result[14]
      xflag[i] = result[15]
      yflag[i] = result[16]
      magflag[i] = result[17]
      reflag[i] = result[18]
      nflag[i] = result[19]
      arflag[i] = result[20]
      paflag[i] = result[21]
      ncomps[i] = result[22]

   # now write the output
   f = open(output_filename, 'w')
   f.write('# 1 FILENAME\n')
   f.write('# 2 X_OUT\n')
   f.write('# 3 X_OUT_ERR\n')
   f.write('# 4 Y_OUT\n')
   f.write('# 5 Y_OUT_ERR\n')
   f.write('# 6 MAG_OUT\n')
   f.write('# 7 MAG_OUT_ERR\n')
   f.write('# 8 RE_OUT\n')
   f.write('# 9 RE_OUT_ERR\n')
   f.write('# 10 N_OUT\n')
   f.write('# 11 N_OUT_ERR\n')
   f.write('# 12 AR_OUT\n')
   f.write('# 13 AR_OUT_ERR\n')
   f.write('# 14 PA_OUT\n')
   f.write('# 15 PA_OUT_ERR\n')
   f.write('# 16 CHI2NU\n')
   f.write('# 17 XFLAG\n')
   f.write('# 18 YFLAG\n')
   f.write('# 19 MAGFLAG\n')
   f.write('# 20 REFLAG\n')
   f.write('# 21 NFLAG\n')
   f.write('# 22 ARFLAG\n')
   f.write('# 23 PAFLAG\n')
   f.write('# 24 NCOMPS\n')
   for i in range(nobj):
      f.write('%s ' % os.path.split(images[i])[-1])
      f.write('%f %f ' % (xobj[i], xobj_err[i]))
      f.write('%f %f ' % (yobj[i], yobj_err[i]))
      f.write('%f %f ' % (mag[i], mag_err[i]))
      f.write('%f %f ' % (Reff[i], Reff_err[i]))
      f.write('%f %f ' % (sersicn[i], sersicn_err[i]))
      f.write('%f %f ' % (ar[i], ar_err[i]))
      f.write('%f %f ' % (pa[i], pa_err[i]))
      f.write('%f ' % chi2nu[i])
      f.write('%d %d ' % (xflag[i],yflag[i]))
      f.write('%d %d ' % (magflag[i],reflag[i]))
      f.write('%d %d ' % (nflag[i],arflag[i]))
      f.write('%d ' % paflag[i])
      f.write('%d ' % ncomps[i])
      f.write('\n')

   f.close()
	