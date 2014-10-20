#!/usr/bin/env python
### NEED TO CONSOLIDATE WITH dropsim_mkinput.py !!!

import numpy as np
import pysynphot as S
import os, sys, time
from sedtools import calzetti00 as c00
from sedtools.bands import filters
from sedtools.arrayspec import ABmag,monoband
from sedtools import magredshift as mr
from sedtools import meiksin
from mypytools.sedtools.bands import filters
from pygoods import *
from sedtools import simkcorr
import pyfits, fitsutil, time
from multiprocessing import Queue, Process, Pool 
# New scripts for calculating colors (6/15/2014)
# following sedtools/galtemp_colors.py
from sedtools import FileSED
import cosmoclass

highz_filters = ['uvimos', 'f435w', 'f475w', 'f555w', 'f606w', 'f625w', 
                  'f775w', 'f814w', 'f850lp', 'f098m', 'f105w', 'f110w', 
                  'f125w', 'f140w', 'f160w']
cosmo_def = cosmoclass.cosmoclass(H0=70., omega_m=0.3, omega_l=0.7)
uni1500 = S.Box(1500.,100.)  # a uniform bandpass of width 100 A around 1500 A

class SimInputPhotomFactory(FileSED.GalaxySEDFactory):
   def __init__(self, sed_filename, cosmo=cosmo_def, sedproperties={}, normmag=-21., normband=1500., ebmv_range=[0.,0.3], ebmv_distlaw='gaussian', ebmv_peak=0.15, ebmv_sigma=0.05, z_range=[3.,6.], z_distlaw='flat', m1500_range=[-25.,-15.], m1500_distlaw='flat'):
      super(SimInputPhotomFactory, self).__init__(sed_filename, 
            cosmo=cosmo, normmag=normmag, normband=normband)
      self.ebmv_range = ebmv_range
      self.ebmv_distlaw = ebmv_distlaw
      self.ebmv_peak = ebmv_peak
      self.ebmv_sigma = ebmv_sigma
      self.z_range = z_range
      self.z_distlaw = z_distlaw
      self.m1500_range = m1500_range
      self.m1500_distlaw = m1500_distlaw
      self.normband = normband

   def calcMagnitudes(self, normmag, z, ebmv, filters):
      """
      Calculate the expected magnitudes for each simulated source in the 
      given list of filters, given their model SED and a distribution in 
      M1500, z, E(B-V), and Lya EW. The filters need to be pysynphot 
      bandpass objects.
      Steps:
      1. Apply dust extinction
      2. Normalize SED to some M_1500
      3. Add Ly-alpha emission (if included)
      4. Redshift
      5. Calculate magnitudes in all filters
      """
      ## NOTE: apply dust attenuation BEFORE normalizing M1500!!
      ## We are interested in the UV luminosity that is NOT obscured by dust
      ## so we need to normalize the SEDs to match the emergent UV luminosity
      ## Therefore, dust only modifies the UV color
      self.add_dust(ebmv, law='xgal')
      self.normalize_abmag(normmag, self.normband)  # e.g., M1500
      # -------------------------------------
      # optional: add Ly-alpha emission here (to be implemented)
      # -------------------------------------
      sp = self.redshift(z)  # Include IGM attenuation
      mags = [sp.ABmag(f) for f in filters]
      self.reset()
      return np.array(mags)

   def run(self, N, filters, output_filename, append=False):
      
      # First, randomly generates M1500, z, and E(B-V) (and maybe Lya EW)
      if self.m1500_distlaw == 'flat':
         m1500 = self.draw_m1500_uniform(N)
      if self.z_distlaw == 'flat':
         z = self.draw_z_uniform(N)
      if self.ebmv_distlaw == 'gaussian':
         ebmv = self.draw_ebmv_gaussian(N)
      else:
         ebmv = self.draw_ebmv_uniform(N)
      # Maybe draw Ly-alpha EW here...
      # initialize the array that stores the output magnitudes
      # mags = np.zeros((len(filters), N))  
      mags = np.zeros((N, len(filters)))
      for i in range(N):
         if i % 10 == 0:
            sys.stdout.write( "Fake source number %d:  \r" % i )
            sys.stdout.flush()
         mags[i] = self.calcMagnitudes(m1500[i], z[i], ebmv[i], filters)
         # self.add_dust(ebmv[i], law='xgal')
         # self.normalize_abmag(m1500[i], uni1500)
         # # print "After normalizing to M1500 = %.2f:" % m1500[i]
         # # print "%s = %.2f" % (filters[0].name, self.copy.ABmag(filters[0]))
         
         # sp = self.redshift(z[i])
         # for j in range(len(filters)):
         #    f = filters[j]
         #    mags[j][i] = sp.ABmag(f)
         # self.reset()
      # write to output
      mags = mags.T
      print "Finished calculation. Now flush output..."
      if append:
         print "Warning: make sure the catalog to append to has the same colums (in the same order) as the new run!"
         f = open(output_filename, 'ab')
      else:
         f = open(output_filename, 'wb')
         # write header
         f.write('## Template SED: %s \n' % self.sp.name)
         f.write('## Calculated on %s \n' % time.ctime())
         for j in range(len(filters)):
            f.write('# %d %s \n' % ((j+1), filters[j].name))
         f.write('# %d M1500 \n' % (len(filters) + 1))
         f.write('# %d REDSHIFT \n' % (len(filters) + 2))
         f.write('# %d EBMV \n' % (len(filters) + 3))
      for i in range(N):
         objstr1 = "  ".join(map(lambda s:"%.4f" % s, mags[:,i]))
         objstr2 = "  ".join(map(lambda s:"%.2f" % s,(m1500[i],z[i],ebmv[i])))
         f.write(objstr1 + "  " + objstr2 + '\n')
      f.close()

   def draw_ebmv_uniform(self, N):
      """
      draw from a uniform distribution of E(B-V)
      """
      # ebmv = random.uniform(extlo,exthi,size=size)
      ebmv = np.random.uniform(*self.ebmv_range, size=N)
      return ebmv

   def draw_ebmv_gaussian(self, N):
      """
      draw from a gaussian distribution of E(B-V)
      """
      #ebmv = random.normal(extpeak,extsigma,size=size)
      # Generate one value at a time... must be a more efficient way of doing 
      # this!
      ebmv = np.zeros(N)
      for i in range(N):
         drawn = 0
         while drawn == 0:
            x = np.random.normal(self.ebmv_peak, self.ebmv_sigma)
            if x > 0:
               ebmv[i] = x
               drawn = 1
               break
      return ebmv

   def draw_ebmv_lognormal(self, N):
      """
      draw from a gaussian distribution of E(B-V)
      self.ebmv_peak becomes log10(peak_ebmv)
      """
      #ebmv = random.normal(extpeak,extsigma,size=size)
      log_ebmv = random.normal(np.log10(self.ebmv_peak), self.ebmv_sigma,
                               size=N)
      ebmv = 10. ** log_ebmv
      return ebmv

   def draw_z_uniform(self, N):
      """
      draw from a uniform distribution of redshift
      """
      z = np.random.uniform(*self.z_range, size=N)
      return z

   def draw_m1500_uniform(self, N):
      """
      draw from a uniform distribution of M_1500
      """
      m1500 = np.random.uniform(*self.m1500_range, size=N)
      return m1500

   def draw_lya(self, ngal, xpdf, ypdf):
      raise NotImplementedError
      # """
      # Draw Lya equivalent width from the given distribution. The distribution is    given by the
      # xpdf and ypdf parameters.
      # xpdf: the value of the random variables
      # ypdf: the probability density at the given value of xpdf
      # """
      # lya_ew = zeros(ngal, 'float')
      # c = sextractor('Lya_EW_dist.txt')
      # i = 0
      # while i < ngal:
         # ew_r = random.uniform(-100., 100.)
         # j = searchsorted(xpdf, ew_r)
         # pdf_r = ypdf[j-1] + (ypdf[j] - ypdf[j-1])/(xpdf[j] - xpdf[j-1]) * (ew_r    - xpdf[j-1])
         # r = random.uniform(0., 1.)
         # if r <= pdf_r:
            # lya_ew[i] = ew_r
            # i += 1
      # return lya_ew

# def draw_age_uniform(size):
#    """draw age from [50, 100, 150] Myr old pool"""
#    ages = zeros(size,'int')
#    for i in range(len(ages)):
#       a = random.uniform(0.,3.)
#       if (a<1.0): ages[i] = 50
#       elif (a<2.0): ages[i] = 100
#       else: ages[i] = 150
#    return ages

# def draw_sfr_uniform(size):
#    """draw SFR from [10,50,100] M_solar/yr pool"""
#    sfr = zeros(size,'int')
#    for i in range(len(sfr)):
#       s = random.uniform(0.,3.)
#       if (s<1.0): sfr[i] = 10
#       elif (s<2.0): sfr[i] = 50
#       else: sfr[i] = 100
#    return sfr



# def mkinput_sed(ngal, M_limits, M_dist, extpar, z_limits, bands, igmroutine=meiksin.meiksin,
#    spec = 'mylbg_sfr10.sed', extdist='uniform', lyadist='Lya_EW_dist.txt',nproc=1):
#    """Generate input catalog for GOODS simulation using template SED.
#       Use randomly drawn values of E(B-V), z, M_1500, etc.
#    """
   
#    sp = S.FileSpectrum(spec)
#    zlo, zhi = z_limits
#    Mlo, Mhi = M_limits

#    # randomly draw E(B-V), z, and M1500
#    # use either a flat distribution ('uniform') or a Gaussian distribution ('gaussian')
#    print "Start drawing E(B-V)"
#    if extdist == 'uniform':
#       ebmvr = draw_ebmv_uniform(extpar[0],extpar[1],size=ngal)
#    elif extdist == 'gaussian':
#       ebmvr = draw_ebmv_gaussian(extpar[0],extpar[1],size=ngal)
#    elif extdist == 'lognormal':
#       extpeak = log10(extpar[0]); extsigma = extpar[1]
#       logebmvr = draw_ebmv_lognormal(extpeak,extsigma,size=ngal)
#       ebmvr = 10.**logebmvr
#    print "finished drawing E(B-V)"
#    zr = draw_z_uniform(zlo,zhi,size=ngal)   # draw a redshift from a flat distribution between zlo and zhi

#    if M_dist == 'uniform':
#       Mr = draw_mag_uniform(Mlo,Mhi,size=ngal)  # draw an M_1500 from a flat distribution between Mlo and Mhi
#    #?? Should I draw from a reasonbly assumed LF?? 
 
#    # draw equivalent width from Lya EW distribution at z~1-3 provided by Naveen Reddy 
#    if lyadist:
#       ld = sextractor(lyadist)
#       lya_ew = draw_lya(ngal, ld.lya_ew, ld.pdf)
#    else:
#       lya_ew = zeros(ngal)
   
#    magnitudes = {}
#    kcorr = {}
#    for b in bands:
#       magnitudes[b] = zeros(ngal,'float')
#       kcorr[b] = zeros(ngal,'float')
#    def worker(n0,n1,qout_dic,qindex_dic):
#       mag_proc = {}
#       kcorr_proc = {}
#       Mr_proc = Mr[n0:n1]
#       zr_proc = zr[n0:n1]
#       lya_ew_proc = lya_ew[n0:n1]
#       ebmvr_proc = ebmvr[n0:n1]
#       for b in bands:
#          mag_proc[b] = zeros(n1-n0)
#          kcorr_proc[b] = zeros(n1-n0)
#       for i in range(n1-n0):
#          if i % 500 == 0:
#             print "%d done." % i
#             sys.stdout.flush()
#          ext = S.Extinction(ebmvr_proc[i],'xgal')
#          spec_ext = sp * ext  # apply dust extinction to rest-frame SED
#          # calculate observed magnitudes (normalized to M_1500 = Mr[i]) of each band, using the drawn redshift & 
#          # Lya EW
#          for b in bands:
#             mag_proc[b][i], kcorr_proc[b][i] = simkcorr.simkcorr(spec_ext, Mr_proc[i], uni1500, filters[b], zr_proc[i], lya_ew_proc[i])
#       # Now paste the results into Queue
#       for b in bands:
#          qout_dic[b].put(mag_proc[b])
#          qindex_dic[b].put([n0,n1])

#    chunksize = float(ngal) / nproc
#    qout_dic = {}
#    qindex_dic = {}
#    for b in bands:
#       qout_dic[b] = Queue()
#       qindex_dic[b] = Queue()
#    procs = []
#    for j in range(nproc):
#       j0 = int(round(chunksize*j))
#       j1 = int(round(chunksize*(j+1)))
#       print j0, j1
#       p = Process(target=worker,args=(j0,j1,qout_dic,qindex_dic))
#       procs += [p]
#       p.start()

#    # now retrieve results from the queues
#    for i in range(nproc):
#       for b in bands:
#          i0,i1 = qindex_dic[b].get()
#          magnitudes[b][i0:i1] = qout_dic[b].get()
   
#    return magnitudes, ebmvr, zr, Mr, lya_ew

   
# # def write_inputcat(catname, bands, magnitudes, zr, ebmvr, Mr, lya_ew):
#    """
#    Write a FITS table instead of a SExtractor-format catalog.
#    """
#    columns = []
#    # first record the magnitudes
#    for b in bands:
#       colname = '%s_mag' % b
#       columns += [pyfits.Column(name=colname,array=magnitudes[b],format='D')]
#    columns += [pyfits.Column(name='z_input',array=zr,format='D')]
#    columns += [pyfits.Column(name='ebmv_input',array=ebmvr,format='D')]
#    columns += [pyfits.Column(name='M1500_input',array=Mr,format='D')]
#    columns += [pyfits.Column(name='Lya_EW_input',array=lya_ew,format='D')]
#    coldefs = pyfits.ColDefs(columns)
#    tbhdu = pyfits.new_table(coldefs)
#    tbhdu.writeto(catname)


# if __name__ == "__main__":
#    parfile = sys.argv[1]
#    c = parseconfig(parfile)
#    N = c['NUMBER']; print 'NUMBER:', N
#    M_limits = [c['MLO'],c['MHI']]; print "M_limits:", M_limits
#    Mdist = c['MDIST']; print 'MDIST:', Mdist
#    extdist = c['EXTDIST']; print 'EXTDIST:', extdist
#    if extdist == 'uniform':
#       extlo = c['EXTLO']; print 'EXTLO:', extlo
#       exthi = c['EXTHI']; print 'EXTHI:', exthi
#       extpar = [extlo,exthi]
#    elif (extdist == 'gaussian') | (extdist=='lognormal'):
#       extpeak = c['EXTPEAK']; print 'EXTPEAK:', extpeak
#       extsigma = c['EXTSIGMA']; print 'EXTSIGMA:', extsigma
#       extpar = [extpeak,extsigma]; print extpar

#    z_limits = [c['ZLO'],c['ZHI']]
#    spec = c['SPEC']
#    bands = c['BANDS']; print "BANDS:", bands
#    catname = c['CATNAME']
#    nproc = c['NPROC']
#    if os.path.exists(catname):
#       os.remove(catname)
#    try:
#       lyadist = c['LYADIST']
#    except:
#       lyadist = ""
#    t1 = time.time()
#    magnitudes, ebmvr, zr, Mr, lya_ew = mkinput_sed(N, M_limits, Mdist, extpar,\
#       z_limits, bands, extdist=extdist, spec=spec, lyadist=lyadist, nproc=nproc)
#    t2 = time.time()
#    print "%.2f seconds passed." % (t2-t1)
#    write_inputcat(catname, bands, magnitudes, zr, ebmvr, Mr, lya_ew)
#    print "DONE"
