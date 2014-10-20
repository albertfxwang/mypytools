#!/usr/bin/env python
### NEED TO CONSOLIDATE WITH mkinput_photometry.py !!!

from numpy import *
import mypytools
import pysynphot as S
import os, sys
from mypytools.sedtools.arrayspec import ABmag,monoband
from mypytools.sedtools import magredshift as mr
from mypytools.sedtools import meiksin
from pygoods import *
from mypytools.sedtools import simkcorr
from mypytools.sedtools.bands import filters
from multiprocessing import Pool, Process, Queue, Manager
import pyfits
import yaml

# Definitions
#sp = S.FileSpectrum('mylbg_sfr10.sed')
# bband = S.ObsBandpass('acs,wfc2,f435w')
# vband = S.ObsBandpass('acs,wfc2,f606w')
# iband = S.ObsBandpass('acs,wfc2,f775w')
# i814band = S.ObsBandpass('acs,wfc2,f814w')
# zband = S.ObsBandpass('acs,wfc2,f850lp')
#uband = S.FileBandpass('U_vimos.bp') # VIMOS U-band
#uni1500 = S.FileBandpass('uni1500.bp')
#uni1500 = S.Box(1500., 100.)
# yband = S.ObsBandpass('wfc3,ir,f105w')
# jband = S.ObsBandpass('wfc3,ir,f125w')
# hband = S.ObsBandpass('wfc3,ir,f160w')

# band_defs = {'acsf435w':bband, 'acsf606w':vband, 'acsf775w':iband, 'acsf850lp':zband,
# 	'wfc3f105w':yband, 'wfc3f125w':jband, 'wfc3f160w':hband, 'acsf814w':i814band}
highz_bands = ['f435w', 'f475w', 'f555w', 'f606w', 'f625w', 'f775w', 'f814w', 
               'f850lp', 'f098m', 'f105w', 'f110w', 'f125w', 'f140w', 'f160w']
galtemplate_dir = '/Users/khuang/Dropbox/codes/mypytools/sedtools/galtemplates'
SED_LBG1 = galtemplate_dir+'/lbgtemp_z3_Zsolar_constSFH_200Myr.sed'
SED_LBG2 = galtemplate_dir+'/lbgtemp_z6_0.4Zsolar_constSFH_100Myr.sed'


def draw_ebmv_uniform(extlo,exthi,size):
   """draw from a uniform distribution of E(B-V)"""
   ebmv = random.uniform(extlo,exthi,size=size)
   return ebmv

def draw_ebmv_gaussian(extpeak,extsigma,size):
   """draw from a gaussian distribution of E(B-V), but enforce E(B-V)>0
   for negative values of E(B-V), redraw from a flat distribution between 0 and 0.5"""
   ebmv = random.normal(extpeak,extsigma,size=size)
   flat_ebmv = random.uniform(0., max(ebmv), size=size)
   ebmv = where(ebmv>0., ebmv, flat_ebmv)
   #ebmv = maximum(ebmv, 0.)
   return ebmv

def draw_ebmv_lognormal(extpeak,extsigma,size):
	"""draw from a lognormal distribution of E(B-V)"""
	logextpeak = log(extpeak)
	logebmvr = random.normal(logextpeak, extsigma, size=size)
	ebmvr = exp(logebmvr)
	return ebmvr
	
def draw_z_uniform(zlo,zhi,size):
   """draw from a uniform distribution of redshift"""
   z = random.uniform(zlo,zhi,size=size)
   return z

def draw_M_uniform(Mlo,Mhi,size):
   """draw from a uniform distribution of M_1500"""
   M = random.uniform(Mlo,Mhi,size=size)
   return M

def draw_age_uniform(size):
   """draw age from [50, 100, 150] Myr old pool"""
   ages = zeros(size,'int')
   for i in range(len(ages)):
      a = random.uniform(0.,3.)
      if (a<1.0): ages[i] = 50
      elif (a<2.0): ages[i] = 100
      else: ages[i] = 150
   return ages

def draw_sfr_uniform(size):
   """draw SFR from [10,50,100] M_solar/yr pool"""
   sfr = zeros(size,'int')
   for i in range(len(sfr)):
      s = random.uniform(0.,3.)
      if (s<1.0): sfr[i] = 10
      elif (s<2.0): sfr[i] = 50
      else: sfr[i] = 100
   return sfr

def draw_lya(ngal, xpdf, ypdf):
   """
   Draw Lya equivalent width from the given distribution. The distribution is given by the
   xpdf and ypdf parameters.
   xpdf: the value of the random variables
   ypdf: the probability density at the given value of xpdf
   """
   lya_ew = zeros(ngal, 'float')
   c = sextractor('Lya_EW_dist.txt')
   i = 0
   while i < ngal:
      ew_r = random.uniform(-100., 100.)
      j = searchsorted(xpdf, ew_r)
      pdf_r = ypdf[j-1] + (ypdf[j] - ypdf[j-1])/(xpdf[j] - xpdf[j-1]) * (ew_r - xpdf[j-1])
      r = random.uniform(0., 1.)
      if r <= pdf_r:
         lya_ew[i] = ew_r
         i += 1
   return lya_ew

def mkinput(ngal, Mlo, Mhi, Mdist, extpar, zlo, zhi, bands,
            igmroutine=meiksin.meiksin, spec=SED_LBG2, extdist='uniform', 
            lyadist='Lya_EW_dist.txt', w0=0., w1=20., restwave=1500.):
   """
   Generate input catalog for LBG simulation
   Use randomly drawn values of E(B-V), z, absolute magnitude at the specified
   rest-frame wavelength, and Ly-alpha equivalent width.
   """
   sp = S.FileSpectrum(spec)
   random.seed()
   # randomly draw E(B-V), z, and M
   # use either a flat distribution ('uniform') or a Gaussian distribution ('gaussian')
   if extdist == 'uniform':
      ebmvr = draw_ebmv_uniform(extpar[0],extpar[1],size=ngal)
   elif extdist == 'gaussian':
      ebmvr = draw_ebmv_gaussian(extpar[0],extpar[1],size=ngal)
   elif extdist == 'lognormal':
      ebmvr = draw_ebmv_lognormal(extpar[0],extpar[1],size=ngal)
   zr = draw_z_uniform(zlo,zhi,size=ngal)   # draw a redshift from a flat distribution between zlo and zhi

   if Mdist == 'uniform':
      Mr = draw_M_uniform(Mlo,Mhi,size=ngal)  
      # draw an M_1500 from a flat distribution between Mlo and Mhi
   elif Mdist == 'single':
      # Use a single value for input absolute magnitude M=-21.0, therefore 
      # will ignore MLO and MHI.
      # This is useful if one is interested in colors only; apparent magnitudes
      # can be scaled later with respect to the reference absolute magnitude
      Mr = ones(ngal) * -21.0
   assert len(Mr) == ngal, "Input absolute magnitudes are not generated."
   #?? Should I draw from a reasonbly assumed LF?? 
 
   # draw equivalent width from Lya EW distribution at z~1-3 provided by Naveen Reddy 
   if lyadist == 'uniform':
      lya_ew = random.uniform(w0, w1, size=ngal)
   elif len(lyadist) > 0:
      ld = sextractor(lyadist)
      lya_ew = draw_lya(ngal, ld.lya_ew, ld.pdf)
   else:
      lya_ew = zeros(ngal)
   if bands == 'highz_bands':
      bands = {}
      for b in highz_bands:
         bands[b] = filters[b]
   else:
      bandnames = copy.copy(bands)
      bands = {}
      for b in bandnames:
         bands[b] = filters[b]

   mags = zeros((len(bands), ngal))
   kc = zeros((len(bands), ngal))
   restband = S.Box(restwave, 100.)

   for i in range(ngal):
      if i % 500 == 0:
         print "%d done." % i
         sys.stdout.flush()
      ext = S.Extinction(ebmvr[i], 'xgal')  
      spec_ext = sp * ext  # apply dust extinction to rest-frame SED
      # calculate observed magnitudes (normalized to M_1500 = Mr[i]) of each band, using the drawn redshift & 
      # Lya EW
      for j in range(len(bands)):
         # b = bands[j]
         b = bands.keys()[j]
         mags[j][i], kc[j][i] = simkcorr.simkcorr(spec_ext, Mr[i], restband, 
                                                  bands[b], zr[i], lya_ew[i])
   results = []
   for i in range(len(bands)):
      results += [mags[i]]
   results += [zr, ebmvr, Mr, lya_ew]
   return array(results)


def write_inputcat_sex(catname, bands, zr, ebmvr, Mr, lya_ew):
   f = open(catname,'w')
   for i in range(len(bands)):
   	f.write('# %d %s\n' % ((i+1), bands[i]))
   f.write('# %d z\n' % (len(bands)+1))
   f.write('# %d EBMV\n' % (len(bands)+2))
   f.write('# %d M [Restframe %.1f A mag]\n' % ((len(bands)+3),restwave))
   f.write('# %d Lya_EW\n' % (len(bands)+4))
   for j in range(len(ebmvr)):
      for i in range(len(bands)):
         f.write('%.3f ' % mags[i][j])
      f.write("%10.7f  %10.7f  %10.7f  " % (zr[j], ebmvr[j], Mr[j]))
      f.write('%.2f ' % (lya_ew[j]))
      f.write("\n")
   f.close()

def write_inputcat_fits(fitsname, colnames, colvalues, colformats):
   # write to a FITS table for quick read
   columns = []
   for i in range(len(colnames)):
      col = pyfits.Column(name=colnames[i], format=colformats[i], array=colvalues[i])
      columns += [col]
   coldefs = pyfits.ColDefs(columns)
   tbhdu = pyfits.new_table(coldefs)
   tbhdu.writeto(fitsname)
   print "Written to %s" % fitsname

def worker(ngal, Mlo, Mhi, Mdist, extpar, zlo, zhi, lyadist, w0, w1, bands, 
           restwave, q_out):
   output = mkinput(ngal, Mlo, Mhi, Mdist, extpar, zlo, zhi, bands, 
      igmroutine=meiksin.meiksin,spec=spec, extdist=extdist, lyadist=lyadist,
      w0=w0, w1=w1, restwave=restwave)
   q_out.put(output)

if __name__ == "__main__":
   parfile = sys.argv[1]
   # c = parseconfig(parfile)
   c = yaml.load(open(parfile))
   N = c['NUMBER']; print 'NUMBER:', N
   Mlo = c['MLO']; print 'MLO:', Mlo
   Mhi = c['MHI']; print 'MHI:', Mhi
   Mdist = c['MDIST']; print 'MDIST:', Mdist
   extdist = c['EXTDIST']; print 'EXTDIST:', extdist
   if extdist == 'uniform':
      extlo = c['EXTLO']; print 'EXTLO:', extlo
      exthi = c['EXTHI']; print 'EXTHI:', exthi
      extpar = [extlo,exthi]
   elif extdist == 'gaussian':
      extpeak = c['EXTPEAK']; print 'EXTPEAK:', extpeak
      extsigma = c['EXTSIGMA']; print 'EXTSIGMA:', extsigma
      extpar = [extpeak,extsigma]
   elif extdist == 'lognormal':
   	extpeak = c['EXTPEAK']; print 'EXTPEAK:', extpeak
   	extsigma = c['EXTSIGMA']; print 'EXTSIGMA:', extsigma
   	extpar = [extpeak,extsigma]
   zlo = c['ZLO']; zhi = c['ZHI']
   spec = c['SPEC']; print spec
   catname = c['CATNAME']
   bands = c['BANDS']
   try:
      lyadist = c['LYADIST']
   except:
      lyadist = ""
   if lyadist == 'uniform':
      w0 = c['LYA_EW0']
      w1 = c['LYA_EW1']
   else: 
      w0 = 0.
      w1 = 0.
   if 'NPROC' in c.keys():
      nproc = c['NPROC']
   else:
   	nproc = 1
   try:
      restwave = c['RESTWAVE']
   except:
      restwave = 1500.
   if nproc > 1:
      # use multiprocessing
      processes = []
      q = Queue()
      results = zeros(0)
      nblock = N / nproc
      for i in range(nproc):
         n0 = i * nblock
         n1 = min((i+1)*nblock, N)
         blocksize = n1 - n0
         print "block %d has size %d" % (i, blocksize)
         p = Process(target=worker, args=(blocksize, Mlo, Mhi, Mdist, extpar, 
            zlo, zhi, lyadist, w0, w1, bands, restwave, q))
         processes += [p]
         p.start()
      # now get results
      for i in range(nproc):
         x = q.get()
         x = array(x)
         if len(results) == 0:
            results = x
         else:
            results = concatenate((array(results), x), axis=1)
         #print shape(results)
      # join the processes
      for i in range(nproc):
         processes[i].join()
      print shape(results)
   else:
   	results = mkinput(N,Mlo,Mhi,Mdist,extpar,zlo,zhi,bands,spec=spec,
   		extdist=extdist,lyadist=lyadist)
   if bands == 'highz_bands':
      bands = highz_bands
   colnames = bands +['z_input','ebmv','M_input','lya_ew']
   print colnames
   colformats = ['D'] * (len(bands)+4)
   colvalues = results
   # write to a FITS table
   
   write_inputcat_fits(catname, colnames, colvalues, colformats)
   h = pyfits.open(catname, mode='update')
   h[0].header.update('INPUTSED', os.path.split(c['SPEC'])[-1])
   h.flush()
   h.close()
   print "DONE"
