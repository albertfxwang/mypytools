#!/usr/bin/env python

from numpy import *
from bivariate import fit_lbg as fl
import cPickle, os, time
from pygoods import *
import cPickle
from bivariate import bivRL as bl

__version__ = "1.0 22 May 2013"
__author__ = "Kuang-Han Huang; JHU/STScI"

"""
Defines class FitUdropsRL (a super class of fit_lbg.FitBivariateRL)
to fit the size-luminosity distribution of U-dropouts in CANDELS
GOODS-S.
"""
# Define the weight limits for each depth
logwhtlims = {}
logwhtlims['udf'] = [5.0, 10.0]
logwhtlims['deep'] = [4.2, 5.0]
logwhtlims['wide'] = [0.0, 4.2]
# ERS does not really need to use log(weight)

def read_kgrid(kgridfile):
   f = open(kgridfile, 'rb')
   kgrid = cPickle.load(f)
   f.close()
   return kgrid

def read_zdgrid(zdgridfile):
   # the same mechanism as reading kgrid
   return read_kgrid(zdgridfile) 

class FitUdropsRL(fl.FitBivariateRL):
   """
   A class to fit U-dropouts in CANDELS GOODS-S.
   """
   def __init__(self, parfile, expand=False):
      """
      Defines the data sets, model limits, parameter boundaries.
      """
      fl.FitBivariateRL.__init__(self, parfile)
      self.catalog = self.DATAFILE
      c = Ftable(self.catalog)
      # U-dropout catalog with GALFIT results
      self.fields = ['udf','deep','ers','wide']

      # Create bl.bivariate_RL_class instances
      self.models = {}
      for f in self.fields:
         self.models[f] = bl.bivariate_RL_class()

      # Define limits
      self.limits = {}
      self.limits['udf'] = array(self.LIMITS1).reshape((2,2))
      self.limits['deep'] = array(self.LIMITS2).reshape((2,2))
      self.limits['ers'] = array(self.LIMITS3).reshape((2,2))
      self.limits['wide'] = array(self.LIMITS4).reshape((2,2))
      self.expand = expand 
      print "self.expand", self.expand

      #print self.limits['deep']
      # create self.dataset
      self.datasets = {}
      self.ra = {}
      self.dec = {}
      self.objid = {}
      # UDF criteria
      udfcrit = (c.wfc3_f160w_logwht>=5.0)
      udfcrit = udfcrit & fl.between(c.f606w_magout_gf,
                                     self.limits['udf'][0])
      udfcrit = udfcrit & fl.between(log10(c.f606w_reout_gf),
                                     self.limits['udf'][1])
      # ERS criteria
      erscrit = (c.wfc3_f098m_weight > 0)
      erscrit = erscrit & fl.between(c.f606w_magout_gf,
                                     self.limits['ers'][0])
      erscrit = erscrit & fl.between(log10(c.f606w_reout_gf),
                                     self.limits['ers'][1])
      # Deep criteria
      deepcrit = (c.wfc3_f098m_weight==0)
      deepcrit = deepcrit & (c.wfc3_f160w_logwht>=4.2) & \
         (c.wfc3_f160w_logwht<5.0)
      deepcrit = deepcrit & fl.between(c.f606w_magout_gf,
                                       self.limits['deep'][0])
      deepcrit = deepcrit & fl.between(log10(c.f606w_reout_gf),
                                       self.limits['deep'][1])
      # Wide criteria
      widecrit = (c.wfc3_f160w_logwht<4.2) & (c.wfc3_f098m_weight==0)
      widecrit = widecrit & fl.between(c.f606w_magout_gf,
                                       self.limits['wide'][0])
      widecrit = widecrit & fl.between(log10(c.f606w_reout_gf),
                                       self.limits['wide'][1])
      vfcrit = (c.vflag != 10) & (c.vflag != 4)
      # Visual classification criteria
      gfcrit = ((c.f606w_magflag == 0) & (c.f606w_reflag == 0))
      gfcrit = (gfcrit & (c.f606w_nflag == 0))
      gfcrit = (gfcrit & (c.f606w_reout_err_gf/c.f606w_reout_gf<=self.REERR_LIM))
      gfcrit = (gfcrit & (c.f606w_reout_err_gf>0.))
      gfcrit = (gfcrit & (c.f606w_chi2nu_gf <= self.CHI2NU_LIM))
      # Data set for UDF
      self.udf_fit = (udfcrit&vfcrit&gfcrit)
      mag_udf = c.f606w_magout_gf[udfcrit&vfcrit&gfcrit]
      logre_udf = log10(c.f606w_reout_gf[udfcrit&vfcrit&gfcrit])
      self.datasets['udf'] = array([mag_udf, logre_udf])
      self.ra['udf'] = c.ra[udfcrit&vfcrit&gfcrit]
      self.dec['udf'] = c.dec[udfcrit&vfcrit&gfcrit]
      self.objid['udf'] = c.number[udfcrit&vfcrit&gfcrit]
      # Data set for GOODS-S Deep
      self.deep_fit = (deepcrit&vfcrit&gfcrit)
      mag_deep = c.f606w_magout_gf[deepcrit&vfcrit&gfcrit]
      logre_deep = log10(c.f606w_reout_gf[deepcrit&vfcrit&gfcrit])
      self.datasets['deep'] = array([mag_deep, logre_deep])
      self.ra['deep'] = c.ra[deepcrit&vfcrit&gfcrit]
      self.dec['deep'] = c.dec[deepcrit&vfcrit&gfcrit]
      self.objid['deep'] = c.number[deepcrit&vfcrit&gfcrit]
      # Data set for ERS
      self.ers_fit = (erscrit&vfcrit&gfcrit)
      mag_ers = c.f606w_magout_gf[erscrit&vfcrit&gfcrit]
      logre_ers = log10(c.f606w_reout_gf[erscrit&vfcrit&gfcrit])
      self.datasets['ers'] = array([mag_ers, logre_ers])
      self.ra['ers'] = c.ra[erscrit&vfcrit&gfcrit]
      self.dec['ers'] = c.dec[erscrit&vfcrit&gfcrit]
      self.objid['ers'] = c.number[erscrit&vfcrit&gfcrit]
      # Data set for wide
      self.wide_fit = (widecrit&vfcrit&gfcrit)
      mag_wide = c.f606w_magout_gf[widecrit&vfcrit&gfcrit]
      logre_wide = log10(c.f606w_reout_gf[widecrit&vfcrit&gfcrit])
      self.datasets['wide'] = array([mag_wide, logre_wide])
      self.ra['wide'] = c.ra[widecrit&vfcrit&gfcrit]
      self.dec['wide'] = c.dec[widecrit&vfcrit&gfcrit]
      self.objid['wide'] = c.number[widecrit&vfcrit&gfcrit]

      print "Total numbers used in fitting:"
      for f in self.fields:
         print "%12s  %6d" % (f.upper(), len(self.datasets[f][0]))

      # Define & read GALFIT transfer function kernel grids
      # Also interpolate the GALFIT completeness map.
      self.tfkgrids = {}
      self.tfkgrids['udf'] = read_kgrid(self.KGRIDFILE1)
      self.tfkgrids['udf'].interpolate_completeness(self.limits['udf'],
                                                    self.pixdx)
      self.tfkgrids['deep'] = read_kgrid(self.KGRIDFILE2)
      self.tfkgrids['deep'].interpolate_completeness(self.limits['deep'],
                                                     self.pixdx)
      self.tfkgrids['ers'] = read_kgrid(self.KGRIDFILE3)
      self.tfkgrids['ers'].interpolate_completeness(self.limits['ers'],
                                                    self.pixdx)
      self.tfkgrids['wide'] = read_kgrid(self.KGRIDFILE4)
      self.tfkgrids['wide'].interpolate_completeness(self.limits['wide'],
                                                     self.pixdx)

      # Define & read P(z) kernels
      # Also calculate zdgrid.Pk
      self.zdgrids = {}
      self.zdgrids['udf'] = read_zdgrid(self.ZDGRIDFILE1)
      self.zdgrids['deep'] = read_zdgrid(self.ZDGRIDFILE2)
      self.zdgrids['ers'] = read_zdgrid(self.ZDGRIDFILE3)
      self.zdgrids['wide'] = read_zdgrid(self.ZDGRIDFILE4)
      # Calculate Pk
      for f in self.fields:
         if self.expand:
            limits_f = self.limits[f].copy()
            limits_f[0][1] += 0.5
            bl.Pijk_dVdz(self.zdgrids[f], self.mc[f], limits_f, self.pixdx)
         else:
            bl.Pijk_dVdz(self.zdgrids[f], self.mc[f], self.limits[f], self.pixdx)

      # Define parameter boundaries
      self.boundary = [self.BOUND_ALPHA, self.BOUND_MSTAR, self.BOUND_R0, \
         self.BOUND_SIGMA, self.BOUND_BETA]
      for i in range(5):
         self.boundary[i] = fl.replaceNone(self.boundary[i])

      # Controls the redshift range around self.z_mean over which to 
      # include P(z)
      self.delta_z = 1.5
      print "self.delta_z = ", self.delta_z

      if self.add_interloper==True:
         self.sd['udf'] = cPickle.load(open(self.SDFILE1))
         self.sd['deep'] = cPickle.load(open(self.SDFILE2))
         self.sd['ers'] = cPickle.load(open(self.SDFILE3))
         self.sd['wide'] = cPickle.load(open(self.SDFILE4))
         self.intfrac['udf'] = Ftable(self.INTFRACFILE1)
         self.intfrac['deep'] = Ftable(self.INTFRACFILE2)
         self.intfrac['ers'] = Ftable(self.INTFRACFILE3)
         self.intfrac['wide'] = Ftable(self.INTFRACFILE4)
         self.mag_lolims['udf'] = self.MAG_LOLIMS[0]
         self.mag_lolims['deep'] = self.MAG_LOLIMS[1]
         self.mag_lolims['ers'] = self.MAG_LOLIMS[2]
         self.mag_lolims['wide'] = self.MAG_LOLIMS[3]


      # Define mlfunc_args
      self.mlfunc_args = [self.datasets, self.limits, self.pixdx,
         self.tfkgrids, self.zdgrids, self.M0,
         self.mc, self.z_mean, self.drop,
         self.add_interloper, self.sd, self.intfrac, self.mag_lolims,
         self.nproc_model, None, None, None, None]
   def set_mlfunc(self):
      # Define/refresh mlfunc_args
      self.mlfunc_args = [self.datasets, self.limits, self.pixdx,
         self.tfkgrids, self.zdgrids, self.M0,
         self.mc, self.z_mean, self.drop,
         self.add_interloper, self.sd, self.intfrac, self.mag_lolims,
         self.nproc_model, None, None, None, None]


def do_fit_udrops_RL(parfile, picklefile=None):
   FU = FitUdropsRL(parfile)
   FU.mlfit()
   if picklefile != None:
      FU.pickle(picklefile)

def do_fit_udrops_RM(parfile, picklefile=None):
   FU = FitUdropsRM(parfile)
   FU.mlfit()
   if picklefile != None:
      FU.pickle(picklefile)

def run_mcmc_udrops_RL(parfile, parameters, chain_catalog):
   FU = FitUdropsRL(parfile)
   FU.build_sampler(FU.NWALKERS)
   FU.mcmc_sample(parameters, iterations=FU.NITER_MCMC, 
                  progress_file=chain_catalog)
   # Now save the results
   # FU.pickle(picklefile)
   # f = open(chain_catalog, 'wb')
   # f.write('# 1 ALPHA\n')
   # f.write('# 2 MSTAR\n')
   # f.write('# 3 LOGR0\n')
   # f.write('# 4 SIGMA\n')
   # f.write('# 5 BETA\n')
   # f.write('# 6 PHISTAR\n')
   # f.write('# ')


if __name__ == "__main__":
   parfile = sys.argv[1]
   mode = sys.argv[2]
   # if len(sys.argv) > 3:
   #    picklefile = sys.argv[3]
   # else:
   #    picklefile = None
   if mode == 'RL':
      do_fit_udrops_RL(parfile, picklefile=picklefile)
   elif mode == 'mcmc_RL':
      chain_catalog = sys.argv[3]
      parameters = sys.argv[4:]
      print "parameters:", parameters
      parameters = array(parameters).astype('float')
      run_mcmc_udrops_RL(parfile, parameters, chain_catalog)
   else:
      do_fit_udrops_RM(parfile, picklefile=picklefile)