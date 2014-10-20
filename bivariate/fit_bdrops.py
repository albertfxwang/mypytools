#!/usr/bin/env python

from numpy import *
import fit_lbg as fl
import cPickle, os, time, sys
from pygoods import *
import bivRL as bl
import bivRM as bm

__version__ = "1.0 05 June 2013"
__author__ = "Kuang-Han Huang; JHU/STScI"

"""
Defines class FitBdropsRL (a super class of fit_lbg.FitBivariateRL)
to fit the size-luminosity distribution of B-dropouts in GOODS and UDF
"""

def read_pickled(pfile):
   f = open(pfile, 'rb')
   p = cPickle.load(f)
   f.close()

def read_kgrid(kgridfile):
   f = open(kgridfile, 'rb')
   kgrid = cPickle.load(f)
   f.close()
   return kgrid

def read_zdgrid(zdgridfile):
   # the same mechanism as reading kgrid
   return read_kgrid(zdgridfile) 

class FitBdropsRL(fl.FitBivariateRL):
   """
   A class to fit B-dropouts in GOODS & UDF.
   """
   def __init__(self, parfile, for_fit=True):
      """
      Defines the data sets, model limits, parameter boundaries.
      catalog1: B-dropout catalog in GOODS
      catalog2: B-dropout catalog in UDF
      GOODS = 1; UDF = 2
      vflags: acceptable visual cull flags in the sample
      """
      fl.FitBivariateRL.__init__(self, parfile)
      c1 = Ftable(self.LBGCAT1)
      c2 = Ftable(self.LBGCAT2)
      self.fields = ['goods', 'udf']

      # Create bl.bivariate_RL_class instances
      self.models = {}
      for f in self.fields:
         self.models[f] = bl.bivariate_RL_class()

      # Define limits
      self.limits = {}
      self.limits['goods'] = array(self.LIMITS1).reshape((2,2))
      self.limits['udf'] = array(self.LIMITS2).reshape((2,2))
      self.delta_z = 1.5
      
      #print self.limits['goods']
      # create self.dataset
      self.datasets = {}
      # GOODS criteria
      # To be included in fitting, objects need to satisfy 1. the vflag 
      # criteria and 2. the GALFIT quality criteria, and 3. the magnitude limit
      goodscrit = in1d(c1.cullflag, self.CULLFLAGS)
      goodscrit = goodscrit & ((c1.reout_err/c1.reout)<=self.REERR_RATIOLIM1)
      goodscrit = goodscrit & (c1.chisqnu<=self.CHISQNULIM1)
      goodscrit = goodscrit & (c1.magout<=self.MAGAUTOLIM1)
      self.goodscrit = goodscrit
      mag_goods = c1.magout[goodscrit==True]
      logre_goods = log10(c1.reout[goodscrit==True])
      self.datasets['goods'] = array([mag_goods, logre_goods])

      # UDF criteria
      udfcrit = in1d(c2.cullflag, self.CULLFLAGS)
      udfcrit = udfcrit & ((c2.reout_err/c2.reout)<=self.REERR_RATIOLIM2)
      udfcrit = udfcrit & (c2.chisqnu<=self.CHISQNULIM2)
      udfcrit = udfcrit & (c2.magout<=self.MAGAUTOLIM2)
      self.udfcrit = udfcrit
      mag_udf = c2.magout[udfcrit==True]
      logre_udf = log10(c2.reout[udfcrit==True])
      self.datasets['udf'] = array([mag_udf, logre_udf])

      if for_fit:
         # Define & read GALFIT transfer function kernel grids
         self.tfkgrids = {}
         self.tfkgrids['goods'] = read_kgrid(self.KERNEL1)
         self.tfkgrids['udf'] = read_kgrid(self.KERNEL2)
         #self.tfkgrids['goods'].interpolate_completeness(self.limits['goods'],
         #                                                self.pixdx)
         #self.tfkgrids['udf'].interpolate_completeness(self.limits['udf'],
         #                                              self.pixdx)

         # Define & read P(z) kernels
         # Also calculate zdgrid.Pk
         self.zdgrids = {}
         self.zdgrids['goods'] = read_zdgrid(self.ZDGRIDFILE1)
         self.zdgrids['udf'] = read_zdgrid(self.ZDGRIDFILE2)

         # Calculate Pk
         for f in self.fields:
            bl.Pijk_dVdz(self.zdgrids[f], self.mc[f], self.limits[f], self.pixdx)

         # Define parameter boundaries
         self.boundary = [self.BOUND_ALPHA, self.BOUND_MSTAR, self.BOUND_R0, \
            self.BOUND_SIGMA, self.BOUND_BETA]
         for i in range(5):
            self.boundary[i] = fl.replaceNone(self.boundary[i])

         # Figure out interloper fraction stuff
         if self.add_interloper==True:
            # First read the size distributions
            self.sd['goods'] = cPickle.load(open(self.SDFILES[0]))
            self.sd['udf'] = cPickle.load(open(self.SDFILES[1]))
            # Then read the interloper fractions
            self.intfrac['goods'] = Ftable(self.INTFRAC_FILES[0])
            self.intfrac['udf'] = Ftable(self.INTFRAC_FILES[1])
            self.mag_lolims['goods'] = self.MAG_LOLIMS[0]
            self.mag_lolims['udf'] = self.MAG_LOLIMS[1]

         # Define mlfunc_args
         self.mlfunc_args = [self.datasets, self.limits, self.pixdx,
            self.tfkgrids, self.zdgrids, self.M0,
            self.mc, self.z_mean, self.drop,
            self.add_interloper, self.sd, self.intfrac, self.mag_lolims,
            self.nproc_model, None, None, None, None]

class FitBdropsRM(FitBdropsRL,fl.FitBivariateRM):
   """
   Performs size-stellar mass distribution fitting. Most of the initialization
   is the same as the RL distribution case, just need an extra argument
   M2M_kernel in the fitting.
   """
   def __init__(self, parfile):
      FitBdropsRL.__init__(self, parfile)
      self.m2m_kernel = read_kgrid(self.M2M_KERNELFILE)
      #self.limits_m = {}
      #self.limits_m['goods'] = array(self.LIMITS_M1).reshape(2,2)
      #self.limits_m['udf'] = array(self.LIMITS_M2).reshape(2,2)
      self.phistar_m = []
      self.logM0 = self.LOGM0
      for f in self.fields:
         self.models[f] = bm.bivariate_RM_class()

      # define mlfunc_args for self.mlfunc_RM
      self.mlfunc_args = [self.datasets, self.limits, 
                          self.pixdx, self.tfkgrids, self.zdgrids, self.logM0,
                          self.mc, self.z_mean, self.drop, self.add_interloper,
                          self.sd, self.intfrac, self.mag_lolims,
                          self.m2m_kernel, self.nproc_model]

def do_fit_bdrops_RL(parfile, picklefile=None):
   FB = FitBdropsRL(parfile)
   FB.mlfit()
   if picklefile != None:
      FB.pickle(picklefile)

def do_fit_bdrops_RM(parfile, picklefile=None):
   FB = FitBdropsRM(parfile)
   FB.mlfit()
   if picklefile != None:
      FB.pickle(picklefile)

if __name__ == "__main__":
   parfile = sys.argv[1]
   mode = sys.argv[2]
   if len(sys.argv) > 3:
      picklefile = sys.argv[3]
   else:
      picklefile = None
   if mode == 'RL':
      do_fit_bdrops_RL(parfile, picklefile=picklefile)
   else:
      do_fit_bdrops_RM(parfile, picklefile=picklefile)
