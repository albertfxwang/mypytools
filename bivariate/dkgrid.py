#!/usr/bin/env python

import numpy as N
from pygoods import sextractor
import cPickle
import dropout_sel_sim as dss
import KPDFadaptnumpy as KPDF
from bivRL import mconvert, A_GOODS, A_UDF
import os
import cosmoclass

"""
Calculates the redshift distribution of a given dropout sample from simulation.
Uses Kernel PDF estimator to calculate density distribution. Then construct a
zdist class object that stores all the redshift distribution information for 
use in bivariate_lf.py in order to construct bivariate model.
"""

def calcdist(zarray, z0, z1, dz):
   """
   Calculates redshift distribution, or equivalently, the completeness function
   in each M_1500 bin C(z) (selection completeness in terms of redshift)
   """
   grid = N.arange(z0, z1, dz)
   if len(zarray) > 4:
      h = KPDF.UPDFOptimumBandwidth(zarray)  # get optimal bandwidth
      p = KPDF.UPDFEpanechnikov(zarray, grid, h)
      p = p / sum(p) # normalize p so that sum(p) = 1.0
   else:
      p = N.zeros(len(grid))
   return p


def contaminant_fraction(c, dropcrit, z0, z1, zmagcut=26.5):
   """
   Calculates the contaminant franction (objects selected as dropouts which are outside
   the target redshift range) for a given set of redshift boundaries.
   """
   zcrit = (c.z_input >= z0) & (c.z_input < z1)
   #n_all = sum(dropcrit & (c.z_mag_auto <= zmagcut))
   n_all = sum(dropcrit)
   n_sel = sum(zcrit & dropcrit)
   print n_all, n_sel
   return (float(n_all-n_sel)/float(n_all))



def makedkcorr(c, Mbincrit, logrbincrit, dropcrit, dmz0, dmz1, ddmz, drop, z0=3.5, z1=4.5):
   """
   Calculates differential K-correction kernel and completeness for a given 
   M1500 bin.
   The completeness is given as the function C(z), which is the completeness of dropout
   selection as a function of redshift in a given magnitude bin.
   
   """
   bincrit = logrbincrit & Mbincrit  # to select objects in the mag & size bin
   zcrit = (c.z_input >= z0) & (c.z_input < z1)
   selcrit = (zcrit & bincrit & dropcrit)  
   # to select all dropouts in this M1500 bin
   zarr = N.compress(selcrit, c.z_input)  # the array of the redshifts of all dropouts in this M1500 bin
   karr = N.zeros(len(zarr))   
   # zcent: reference redshift (zcent = 4.0 for B-dropouts; zcent = 5.0 for V-dropouts)
   if drop == 'b':
      mc = mconvert('M1500_to_i.txt')
      zcent = 4.0
   elif drop == 'v':
      mc = mconvert('M1500_to_z.txt')
      zcent = 5.0
   #n_ztarget = sum(zcrit & bincrit & c.detect)  # NOT including SE detection incompleteness
   n_ztarget = sum(zcrit & bincrit)  # including SE detection incompleteness
   n_zselect = sum(zcrit & bincrit & dropcrit)
      
   # calculate completeness as follows:
   # completeness = (num. of objects that are detected and selected as dropouts) / 
   #                (num. of INPUT objects in this bin)
   if n_ztarget > 5:   
      completeness = float(n_zselect) / float(n_ztarget)
   else: completeness = 0.
   # calculate the K-correction of all redshifts
   for i in range(len(zarr)):
      karr[i] = mc(zarr[i])
   dkarr = karr - mc(zcent)  # the differential K-correction
   grid = N.arange(dmz0, dmz1+ddmz, ddmz)  # the grid over which differential K-correction will be calculated
   if n_zselect > 4:
      h = KPDF.UPDFOptimumBandwidth(dkarr)
      p_dk = KPDF.UPDFEpanechnikov(dkarr, grid, h)
      if sum(p_dk) > 0:
         p_dk = p_dk / sum(p_dk)
      else: p_dk = N.zeros(len(grid))
   else:
      p_dk = N.zeros(len(grid))
   # include zarr in the kernel for inspection purposes
   return completeness, p_dk, n_ztarget, n_zselect, zarr
   

def makedkcorr2(c, Mbincrit, dropcrit, dmzarray, zarray, z0, z1):
   """
   Use the completeness function C(z) in each magnitude bin to construct the 
   'differential K-correction' kernel.
   z0, z1: the redshift boundary defining the dropout sample in question
   """
   Fz = N.zeros(len(zarray)-1)  # fraction of total input objects within each redshift bin
   Cz = N.zeros(len(zarray)-1)  # completeness of selected object in each redshift bin
   nin_array = N.zeros(len(zarray)-1, 'int')
   zcrit = ((c.z_input >= z0) & (c.z_input < z1))
   n_all = sum(zcrit & Mbincrit)
   for i in range(len(Fz)-1):
      if zarray[i] > z1: continue
      elif zarray[i+1] < z0: continue
      else:
         zlo = zarray[i]; zhi = zarray[i+1]
         if zarray[i] < z0: zlo = z0
         if zarray[i+1] > z1: zhi = z1
         zbincrit = (c.z_input>=zlo) & (c.z_input<zhi)
         n_zbin = sum(zbincrit & Mbincrit)
         nsel_zbin = sum(zbincrit & dropcrit & Mbincrit)
         if n_all > 0:
            Fz[i] = float(n_zbin) / float(n_all)
         else: Fz[i] = 0.
         if n_zbin > 0:
            Cz[i] = float(nsel_zbin) / float(n_zbin)
         else: Cz[i] = 0.
         nin_array[i] = n_zbin
   #FzCz = Fz * Cz
   return Fz, Cz, zarray, nin_array


class dkcorr:
   """
   The differential K-correction kernel object.
   Terminologies:
   K-correction here is defined as the difference between apparent magnitude m at the first
   band and the absolute magnitude at the second band M. For example, the K-correction between
   the observed z-band and rest-frame 1500 A is K = m_z - M_1500. This is a different definition
   from the official definition of K-correction (which takes away distance modulus), but it's
   convenient for the purpose here (or I'm just lazy to get all the terminologies right). 
   Apparently K depends on the redshift of the source, the observed and rest-frame bands, 
   and the cosmology model used.
   Differential K-correction means the K-correction at a given redshift z minus the K-correction
   value at the reference redshift z0, using the same bands and the same cosmology models.
   """
   def __init__(self, M0, M1, logR0, logR1, completeness, p, n_target, n_select, zarr):
      self.M0 = M0
      self.M1 = M1
      self.logR0 = logR0
      self.logR1 = logR1
      self.completeness = completeness
      self.dkcorr = p
      self.n_target = n_target
      self.n_select = n_select
      self.zarr = zarr
      # the probability that the differential K-correction at dmz_k
      # is p_k, where p_k is the kth element of p
   #def __init__(self, M0, M1, dmzarray, Fz, Cz, nin_array):
   #   self.M0 = M0
   #   self.M1 = M1
   #   self.dmzarray = dmzarray # redshifts where completeness are defined
   #   self.Fz = Fz  # the F(z) function, the fraction of total in each redshift bin
   #   self.Cz = Cz  # the C(z), or redshift completeness function
   #   self.nin_array = nin_array


class dkgrid:
   """ 
   The grid of differential K-correction kernels at different M1500 bins.
   """
   def __init__(self, M0, M1, dM, logR0, logR1, dlogR, dmz0, dmz1, ddmz, drop, z0, z1):
      # dmz: the bin width of the differential K-correction
      self.M0 = M0
      self.M1 = M1
      self.dM = dM
      self.logR0 = logR0
      self.logR1 = logR1
      self.dlogR = dlogR
      self.dmz0 = dmz0   # lower bound of differential K-correction
      self.dmz1 = dmz1  # higher bound of differential K-correction
      self.ddmz = ddmz  # step size in differential K-correction, should be the same as model
                      # pixel width in the magnitude direction
      self.dmzarray = N.arange(dmz0, dmz1+ddmz, ddmz)  # the bin boundaries 
      self.z0 = z0
      self.z1 = z1
      if drop == 'b': 
         zcent = 4.0
         mc = mconvert('M1500_to_i.txt')
      elif drop == 'v': 
         zcent = 5.0  # the fiducial central redshift of the given dropout sample
         mc = mconvert('M1500_to_z.txt')
      zarray = N.zeros(len(self.dmzarray))  
      # the redshifts corresponding to the magnitude of the differential K-corr dmz
      for i in range(len(self.dmzarray)):
         zarray[i] = mc.revert2z(self.dmzarray[i] + mc(zcent))
      self.zarray = zarray  # zarray could be larger than the [z0, z1] interval... be careful
      
      self.Mlolims = N.arange(M0, M1, dM)
      self.Mhilims = self.Mlolims + dM
      self.logRlolims = N.arange(logR0, logR1, dlogR)
      self.logRhilims = self.logRlolims + dlogR
      if self.Mhilims[-1] > M1: self.Mhilims[-1] = M1
      self.dkcorr = {}
      if drop == 'b':
         self.drop = 'b'
         self.dropout = "B-dropout"
      if drop == 'v':
         self.drop = 'v'
         self.dropout = "V-dropout"

   def __call__(self, simcat, field, imaflags_thresh=4, bstonlim=5.0):
      # use simcat to make differential K-correction kernels in each M1500 bin
      c = sextractor(simcat)
      # do dropout selection
      if field == 'goods':
         zmagcut=26.5
      elif field == 'udf':
      	zmagcut=28.5
      if self.drop == 'b':
         dropcrit, bmags, vmags, zmags = dss.bdrops_sel_sim(c, field=field,\
            imaflags_thresh=imaflags_thresh)
      if self.drop == 'v':
         dropcrit, vmags, imags, zmags = dss.vdrops_sel_sim(c, field=field,\
            imaflags_thresh=imaflags_thresh, bstonlim=bstonlim)
      
      # iterates through each M1500 bin
      for i in range(len(self.Mlolims)):
         for j in range(len(self.logRlolims)):
            Mbincrit = ((c.m1500_input >= self.Mlolims[i]) & (c.m1500_input < self.Mhilims[i]))
            logrbincrit = ((N.log10(c.re_input) >= self.logRlolims[j]) & (N.log10(c.re_input) <\
               self.logRhilims[j]))
            # calls makedkcorr() to make kernels
            completeness, p_dk, n_ztarget, n_zselect, zarr = makedkcorr(c, Mbincrit,\
               logrbincrit,dropcrit, self.dmz0, self.dmz1, self.ddmz, self.drop,  z0=self.z0,\
               z1=self.z1)
            
            #Fz, Cz, zarray, nin_array = makedkcorr2(c, Mbincrit, dropcrit, self.dmzarray,\
            #   self.zarray, self.z0, self.z1)
            #print i, self.Mlolims[i], self.Mhilims[i], max(Cz)
            
            #k = dkcorr(self.Mlolims[i], self.Mhilims[i], self.dmzarray, Fz, Cz, nin_array)
            k = dkcorr(self.Mlolims[i], self.Mhilims[i], self.logRlolims[j],\
               self.logRhilims[j], completeness, p_dk, n_ztarget, n_zselect, zarr)
            self.dkcorr[(i,j)] = k

   def delta_grid(self):
      # make a "perfect" kernel grid, for testing
      for i in range(len(self.Mlolims)):
         for j in range(len(self.logRlolims)):
            dmzarray = N.arange(self.dmz0, self.dmz1+self.ddmz, self.ddmz)
            p_dk = N.zeros(len(dmzarray))
            nc = len(p_dk) / 2
            p_dk[nc] = 1.0
            k = dkcorr(self.Mlolims[i], self.Mhilims[i], self.logRlolims[j],\
               self.logRhilims[j], 1.0, p_dk, 1, 1, [0.0])
            self.dkcorr[(i,j)] = k

   def uniform_grid(self):
      # make a uniform z distribution kernel grid, for testing
      for i in range(len(self.Mlolims)):
         for j in range(len(self.logRlolims)):
            dmzarray = N.arange(self.dmz0, self.dmz1+self.ddmz, self.ddmz)
            p_dk = N.ones(len(dmzarray)) / float(len(dmzarray))
            k = dkcorr(self.Mlolims[i], self.Mhilims[i], self.logRlolims[j],\
              self.logRhilims[j], 1.0, p_dk, 1, 1, [0.0])
            self.dkcorr[(i,j)] = k


def make_dkgrid(simcat, drop, filename, ddmz=0.02, zmagcut=26.5, z0=3.6, z1=4.4,\
   imaflags_thresh=4, bstonlim=2.0):
   print "Remember to use the proper zmagcut for the HUDF sample!!!"
   # B-drops redshift range: 3.6 < z < 4.4
   # V-drops redshift range: 4.6 < z < 5.6
   M0 = -25.0
   M1 = -15.0
   dM = 0.5
   logR0 = -0.6
   logR1 = 1.8
   dlogR = 0.2
   dmz0 = -1.0
   dmz1 = 1.0
   print "zmagcut =", zmagcut
   print "dmz0 =", dmz0
   print "dmz1 =", dmz1
   print "ddmz =", ddmz
   print "z0 =", z0
   print "z1 =", z1
   
   dk = dkgrid(M0, M1, dM, logR0, logR1, dlogR, dmz0, dmz1, ddmz, drop, z0, z1)
   dk(simcat, zmagcut=zmagcut, imaflags_thresh=imaflags_thresh, bstonlim=bstonlim)
   write_zdgrid(dk, filename)
   print "done."


def make_all_dkgrid_goods(simcat):
   os.system('cp dkgrid_bdrops.p dkgrid_bdrops.p.OLD')
   make_dkgrid(simcat, 'b', 'dkgrid_bdrops.p', z0=3.6, z1=4.4) # B-drops in GOODS
   os.system('cp dkgrid_vdrops.p dkgrid_vdrops.p.OLD')
   make_dkgrid(simcat, 'v', 'dkgrid_vdrops.p', z0=4.6, z1=5.6, bstonlim=5.0)
   # V-drops in GOODS
   return 0
