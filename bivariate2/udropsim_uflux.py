#!/usr/bin/env python

from numpy import *
from pygoods import sextractor, Ftable
import cPickle
import KPDFadaptnumpy as KPDF
import os
import fitsutil
from stats import draw_from_pdf as dfp

#uvimos_magzero = 25.128
uvimos_magzero = 26.158
# tfitsim_catalog = '/Users/khuang/Dropbox/Research/bivariate/bivariate_fit/udrops_fitting/simcatalogs/run1_udf_130625.fits'
# c = Ftable(tfitsim_catalog)
# uvimos_fluxerr = uvimos_magzero - 2.5 * log10(c.uvimos_fitquerr)


def calc_umag_pdf(tfitsimcatalog, umag_array=arange(21.,35.0,0.5),
                  pdfname='uvimos_mag_pdf.p', b_snlim=3.0):
   """
   Calculate the PDF of output magnitude around 
   input value of U-band magnitudes when S/N>=1 in U-band.
   """
   c = Ftable(tfitsimcatalog)
   include = (c.uvimos_mag_out>0.)&(c.detect==True)&(c.uvimos_mag_out<99.)&\
             (c.f435w_magerr_iso < (1.0857/b_snlim)) & \
             (c.uvimos_fitqty > 0.)
   print "sum(include)=", sum(include)
   umag_in = c.uvimos_mag_in[include==True]
   umag_out = c.uvimos_mag_out[include==True]
   dmag_array = arange(-10.,6.05,0.05)  # x-coordinates of the PDF
   dumag = umag_array[1]-umag_array[0]  # the bin width of U-band input magnitude
   kpdf_umag = {}
   kpdf_umag['h'] = {}
   kpdf_umag['dmag_array'] = dmag_array
   for umag in umag_array[:-1]:
      mag0 = umag
      mag1 = mag0+dumag
      print mag0,mag1
      bincrit = (umag_in>=mag0)&(umag_in<mag1)
      mag_in_bin = umag_in[bincrit]
      mag_out_bin = umag_out[bincrit]
      if sum(bincrit)>=5:
         h = KPDF.UPDFOptimumBandwidth(mag_out_bin)
         kpdf_umag['h']['%.1f'%mag0] = h
         pdf = KPDF.UPDFEpanechnikov((mag_out_bin-mag_in_bin),dmag_array,h)
      else:
         pdf = ones(len(dmag_array))
         pdf = pdf / sum(pdf)
      kpdf_umag['%.1f'%mag0] = pdf
   f = open(pdfname,'wb')
   cPickle.dump(kpdf_umag,f,2)
   f.close()

def calc_ulimmag_pdf(tfitsimcatalog, pdfname, 
                     limmag_array=arange(26.,31.,0.02)):
   c = Ftable(tfitsimcatalog)
   uvimos_1sig_mag = uvimos_magzero - 2.5 * log10(c.uvimos_fitquerr)
   umag_crit = uvimos_1sig_mag>0
   uvimos_1sig_mag = uvimos_1sig_mag[umag_crit]
   kpdf_ulimmag = {}
   kpdf_ulimmag['dmag_array'] = limmag_array
   h = KPDF.UPDFOptimumBandwidth(uvimos_1sig_mag)
   kpdf_ulimmag['h'] = h
   pdf = KPDF.UPDFEpanechnikov(uvimos_1sig_mag, limmag_array, h)
   kpdf_ulimmag['pdf'] = pdf
   # Also calculate the fraction that have S/N < 1
   dlimmag = limmag_array[1] - limmag_array[0]
   limmag_bins = concatenate([limmag_array, [limmag_array[-1]+dlimmag]])
   uvimos_ston = (c.uvimos_fitqty/c.uvimos_fitquerr)[umag_crit]
   # n1 = histogram(uvimos_1sig_mag[uvimos_ston<1.], limmag_bins)[0]
   # n2 = histogram(uvimos_1sig_mag, limmag_bins)[0]
   # limfrac = n1.astype('float') / maximum(1.,n2.astype('float'))
   # kpdf_ulimmag['limfrac'] = limfrac
   f = open(pdfname, 'wb')
   cPickle.dump(kpdf_ulimmag, f, 2)
   f.close()

def kpdf_from_array(x, xarr, pdf):
   """
   Given an array of coordinates xarr and the PDF defined at these coordinates, return the value
   of PDF(x) at the location x. Linearly interpolate between the discrete sampling of PDF.
   xarr must be a sorted, ascending array of coordinates. 
   pdf must be a non-negative array.
   """
   i = xarr.searchsorted(x)
   if i==0 or i==len(xarr):
      raise ValueError, "x (%.2f) is not within the range of [%.2f,%.2f]." % (x,xarr[0],xarr[-1])
   xarr0 = xarr[i-1]; pdf0 = pdf[i-1]
   xarr1 = xarr[i]; pdf1 = pdf[i]
   p = pdf0 + (x-xarr0)*(pdf1-pdf0)/(xarr1-xarr0)  # linear interpolation
   return p

def draw_umag(u_mag_in, pdffile, onesigma=False, ndraw=3, mag0=21.0, 
              mag1=34.0):
   """
   Draw VIMOS U-band magnitudes given U-band input magnitudes.
   umag_out = umag_in + dmag_drawn
   (Not implemented yet) For each U-band input magnitude, draw ndraw magnitudes from the PDF.
   """
   dmagbin = 0.5  # the bin width of the input magnitudes in the PDF file
   if onesigma==True:
      dmag0 = 25.0
      dmag1 = 29.95
   else:
      dmag0 = -10.
      dmag1 = 6.
   # First read the PDF file
   f = open(pdffile,'rb')
   kpdf = cPickle.load(f)
   f.close()
   # Now read the U-dropout simulation catalog
   #c = Ftable(udropsimcatalog)
   #u_mag_in = c.u_mag_in
   nobj = len(u_mag_in)
   #dmag_drawn = zeros(len(u_mag_in))
   u_mag_drawn = zeros(len(u_mag_in))
   # draw dmag efficiently by binning the input magnitudes
   index_sort = argsort(u_mag_in)  # the index of each element in the sorted array from the unsorted array
   u_mag_in_sort = sort(u_mag_in)  # e.g., u_mag_in_sort[0] = u_mag_in[index_sort[0]]
   # A quick fix right now for input U-band magnitudes brighter than mag0 or fainter than mag1:
   # if U_mag_in < mag0, use the bin in mag0 to draw dmag
   # if U_mag_in > mag1, use the bin in mag1 to draw dmag
   # First, draw magnitudes for U_mag_in < mag0
   #i00 = u_mag_in_sort.searchsorted(mag0)
   #if i00 > 0:
   #   pdf = kpdf[mag0]
   #   mag_in_bin = u_mag_in_sort[:i0]
   #   index_sort_bin = index_sort[:i0]
   #   def pdf_func(x):
   #      return kpdf_from_array(x,kpdf['dmag_array'],pdf)
   #   dmag_drawn_bin = dfp.draw_from_pdf_func(i0,pdf_func,[dmag0,dmag1])

   marr = arange(mag0,mag1,dmagbin)
   if min(u_mag_in) < mag0:
      marr = concatenate([[min(u_mag_in)],marr])
   if max(u_mag_in) > mag1:
      marr = concatenate([marr,[max(u_mag_in)]])
   for m in marr:
      if m < mag0:
         i0 = 0
         i1 = u_mag_in_sort.searchsorted(mag0)
         k = "%.1f" % mag0
      elif m > mag1:
         i0 = u_mag_in_sort.searchsorted(mag1)
         i1 = len(u_mag_in_sort)
         #k = "%.1f" % mag1
         # For u_mag_in > mag1: directly assign 99 to drawn magnitude (but not the 1-sigma)
         if onesigma==False:
            index_sort_bin = index_sort[i0:i1]
            for j in range(i1-i0):
               jj = index_sort_bin[j]
               u_mag_drawn[jj] = 99.0
            continue
      else:
         i0 = u_mag_in_sort.searchsorted(m)
         i1 = u_mag_in_sort.searchsorted(m+dmagbin)
         k = "%.1f" % m
      #if i1 > i0:
      pdf = kpdf[k]
      mag_in_bin = u_mag_in_sort[i0:i1]
      index_sort_bin = index_sort[i0:i1]
      #def pdf_func(x):
      #   return kpdf_from_array(x,kpdf['dmag_array'],pdf)
      dmag_drawn_bin = dfp.draw_from_pdf((i1-i0), pdf, kpdf['dmag_array'])
      #dmag_drawn_bin = dfp.draw_from_pdf_func((i1-i0),pdf_func,[dmag0,dmag1])
      # Now map the drawn magnitudes back to the unsorted array u_mag_in
      for j in range(i1-i0):
         jj = index_sort_bin[j]
         if onesigma==True:
            u_mag_drawn[jj] = dmag_drawn_bin[j]
         else:
            u_mag_drawn[jj] = dmag_drawn_bin[j]+u_mag_in[jj]
      print "drawn %s" % k
   
   return u_mag_drawn

def draw_ulimmag(nobj, pdffile):
   """
   Randomly draw a 1-sigma limiting magnitudes for VIMOS U-band, and this is
   independent of the input U-band magnitude (more dependent on the position?)
   """
   f = open(pdffile, 'rb')
   x = cPickle.load(f)
   f.close()
   ulimmag_drawn = dfp.draw_from_pdf(nobj, x['pdf'], x['dmag_array'])
   return ulimmag_drawn

