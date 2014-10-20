#!/usr/bin/env python

from numpy import *
import fit_lbg as fl
import cPickle, os, time
from pygoods import *
import cPickle
import bivRL as bl
import yaml

__version__ = "1.0 22 May 2013"
__author__ = "Kuang-Han Huang; JHU/STScI"

"""
Defines class FitIdropsRL (a super class of fit_lbg.FitBivariateRL)
to fit the size-luminosity distribution of i-dropouts in CANDELS
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

class FitIdropsRL(fl.FitBivariateRL):
   """
   A class to fit i-dropouts in CANDELS GOODS-S.
   """
   def __init__(self, parfile, reerr_lim=0.6, udf_maglim=30.0, expand=False):
      """
      Defines the data sets, model limits, parameter boundaries.
      udf_maglim: the magnitude limit in UDF to which we extend the fitting
      to; beyond self.limits['udf'][0][1] I only fit 1D magnitude distribution
      using SExtractor corrections.
      """
      fl.FitBivariateRL.__init__(self, parfile)
      self.catalog = self.DATAFILE
      c = Ftable(self.catalog)
      # i-dropout catalog with GALFIT results
      self.fields = ['udf','deep','ers','wide']
      #self.fields = ['deep', 'ers', 'wide']
      self.delta_z = 1.5
      self.expand = expand
      
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

      # Define mconvert instances
      self.mc = {}
      self.mc['udf'] = bl.mconvert(self.MCFILE1)
      self.mc['deep'] = bl.mconvert(self.MCFILE2)
      self.mc['ers'] = bl.mconvert(self.MCFILE3)
      self.mc['wide'] = bl.mconvert(self.MCFILE4)

      #print self.limits['deep']
      # create self.dataset
      self.datasets = {}
      self.ra = {}
      self.dec = {}
      self.objid = {}
      # UDF criteria
      # udfcrit = (c.udf==True)
      udfcrit = (c.wfc3_f160w_logwht>=5.0)
      # udfcrit = udfcrit & fl.between(c.f105w_magout_gf,
      #                                self.limits['udf'][0])
      # udfcrit = udfcrit & fl.between(log10(c.f105w_reout_gf),
      #                                self.limits['udf'][1])
      udfcrit = udfcrit & fl.between(c.f125w_magout_gf,
                                     self.limits['udf'][0])
      udfcrit = udfcrit & fl.between(log10(c.f125w_reout_gf),
                                     self.limits['udf'][1])
      udfcrit = udfcrit & (c.interloper_flag_udf==False)
      # ERS criteria
      # erscrit = (c.ers==True)
      erscrit = (c.wfc3_f098m_weight > 0)
      # erscrit = erscrit & fl.between(c.f098m_magout_gf,
      #                                self.limits['ers'][0])
      # erscrit = erscrit & fl.between(log10(c.f098m_reout_gf),
      #                                self.limits['ers'][1])
      erscrit = erscrit & fl.between(c.f125w_magout_gf,
                                     self.limits['ers'][0])
      erscrit = erscrit & fl.between(log10(c.f125w_reout_gf),
                                     self.limits['ers'][1])
      erscrit = erscrit & (c.interloper_flag_ers==False)
      # Deep criteria
      # deepcrit = (c.deep==True) 
      deepcrit = (c.wfc3_f098m_weight==0) & ((c.wfc3_f160w_logwht<5.0))
      deepcrit = deepcrit & (c.wfc3_f160w_logwht>=4.2) & \
         (c.wfc3_f160w_logwht<5.0)
      # deepcrit = deepcrit & fl.between(c.f105w_magout_gf,
      #                                  self.limits['deep'][0])
      # deepcrit = deepcrit & fl.between(log10(c.f105w_reout_gf),
      #                                  self.limits['deep'][1])
      deepcrit = deepcrit & fl.between(c.f125w_magout_gf,
                                       self.limits['deep'][0])
      deepcrit = deepcrit & fl.between(log10(c.f125w_reout_gf),
                                       self.limits['deep'][1])
      deepcrit = deepcrit & (c.interloper_flag_deep==False)
      # Wide criteria
      # widecrit = (c.wide==True)
      widecrit = (c.wfc3_f160w_logwht<4.2) & (c.wfc3_f098m_weight==0)
      # widecrit = widecrit & fl.between(c.f105w_magout_gf,
      #                                  self.limits['wide'][0])
      # widecrit = widecrit & fl.between(log10(c.f105w_reout_gf),
      #                                  self.limits['wide'][1])
      widecrit = widecrit & fl.between(c.f125w_magout_gf,
                                       self.limits['wide'][0])
      widecrit = widecrit & fl.between(log10(c.f125w_reout_gf),
                                       self.limits['wide'][1])
      widecrit = widecrit & (c.interloper_flag_wide==False)
      # Below is probably redundant, but just to be safe...
      vfcrit = (c.vflag != 10) & (c.vflag != 4)
      # Visual classification criteria
      gfcrit = ((c.f125w_magflag == 0) & (c.f125w_reflag == 0))
      gfcrit = (gfcrit & (c.f125w_nflag == 0))
      gfcrit = (gfcrit & (c.f125w_reout_err_gf/c.f125w_reout_gf<=self.REERR_LIM))
      gfcrit = (gfcrit & (c.f125w_chi2nu_gf<=self.CHI2NU_LIM))
      # Data set for UDF
      # mag_udf = c.f105w_magout_gf[udfcrit]
      # logre_udf = log10(c.f105w_reout_gf[udfcrit])
      mag_udf = c.f125w_magout_gf[udfcrit&gfcrit&vfcrit]
      logre_udf = log10(c.f125w_reout_gf[udfcrit&gfcrit&vfcrit])
      self.datasets['udf'] = array([mag_udf, logre_udf])
      self.ra['udf'] = c.ra[udfcrit&vfcrit&gfcrit]
      self.dec['udf'] = c.dec[udfcrit&vfcrit&gfcrit]
      self.objid['udf'] = c.number[udfcrit&vfcrit&gfcrit]
      # Data set for GOODS-S Deep
      # mag_deep = c.f105w_magout_gf[deepcrit]
      # logre_deep = log10(c.f105w_reout_gf[deepcrit])
      mag_deep = c.f125w_magout_gf[deepcrit&gfcrit&vfcrit]
      logre_deep = log10(c.f125w_reout_gf[deepcrit&gfcrit&vfcrit])
      self.datasets['deep'] = array([mag_deep, logre_deep])
      self.ra['deep'] = c.ra[deepcrit&vfcrit&gfcrit]
      self.dec['deep'] = c.dec[deepcrit&vfcrit&gfcrit]
      self.objid['deep'] = c.number[deepcrit&vfcrit&gfcrit]
      # Data set for ERS
      # mag_ers = c.f098m_magout_gf[erscrit]
      # logre_ers = log10(c.f098m_reout_gf[erscrit])
      mag_ers = c.f125w_magout_gf[erscrit&gfcrit&vfcrit]
      logre_ers = log10(c.f125w_reout_gf[erscrit&gfcrit&vfcrit])
      self.datasets['ers'] = array([mag_ers, logre_ers])
      self.ra['ers'] = c.ra[erscrit&vfcrit&gfcrit]
      self.dec['ers'] = c.dec[erscrit&vfcrit&gfcrit]
      self.objid['ers'] = c.number[erscrit&vfcrit&gfcrit]
      # Data set for wide
      # mag_wide = c.f105w_magout_gf[widecrit]
      # logre_wide = log10(c.f105w_reout_gf[widecrit])
      mag_wide = c.f125w_magout_gf[widecrit&gfcrit&vfcrit]
      logre_wide = log10(c.f125w_reout_gf[widecrit&gfcrit&vfcrit])
      self.datasets['wide'] = array([mag_wide, logre_wide])
      self.ra['wide'] = c.ra[widecrit&vfcrit&gfcrit]
      self.dec['wide'] = c.dec[widecrit&vfcrit&gfcrit]
      self.objid['wide'] = c.number[widecrit&vfcrit&gfcrit]

      print "Total numbers used in fitting:"
      for f in self.fields:
         print "%12s  %6d" % (f.upper(), len(self.datasets[f][0]))

      # Define & read GALFIT transfer function kernel grids
      self.tfkgrids = {}
      self.tfkgrids['udf'] = read_kgrid(self.KGRIDFILE1)
      self.tfkgrids['deep'] = read_kgrid(self.KGRIDFILE2)
      self.tfkgrids['ers'] = read_kgrid(self.KGRIDFILE3)
      self.tfkgrids['wide'] = read_kgrid(self.KGRIDFILE4)
      for f in self.fields:
         self.tfkgrids[f].interpolate_completeness(self.limits[f], self.pixdx)

      # Define & read P(z) kernels
      # Also calculate zdgrid.Pk
      self.zdgrids = {}
      self.zdgrids['udf'] = read_zdgrid(self.ZDGRIDFILE1)
      self.zdgrids['deep'] = read_zdgrid(self.ZDGRIDFILE2)
      self.zdgrids['ers'] = read_zdgrid(self.ZDGRIDFILE3)
      self.zdgrids['wide'] = read_zdgrid(self.ZDGRIDFILE4)
      # Calculate Pk
      for f in self.fields:
         bl.Pijk_dVdz(self.zdgrids[f], self.mc[f], self.limits[f], self.pixdx)

      # Define parameter boundaries
      self.boundary = [self.BOUND_ALPHA, self.BOUND_MSTAR, self.BOUND_R0, \
         self.BOUND_SIGMA, self.BOUND_BETA]
      for i in range(5):
         self.boundary[i] = fl.replaceNone(self.boundary[i])

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

      # Specifically to extend UDF magnitude limit to constrain faint-end
      # slope for i-dropouts
      self.zdist_mag_udf = cPickle.load(open(self.ZDIST_MAG_UDF))
      limits_ext = [self.limits['udf'][0][1]-0.5, self.MAGLIM_UDF]
      self.zdist_mag_udf.calc_Pk(limits_ext, self.pixdx[0], 
                                 self.mc['udf'])
      self.kgrid1d_udf = cPickle.load(open(self.KGRID1D_UDF))
      mag_array_udf = c.wfc3_f125w_mag[(c.udf==True) & \
                        (c.wfc3_f125w_mag>=self.limits['udf'][0][1]) & \
                        (c.wfc3_f125w_mag<self.MAGLIM_UDF) & \
                        ((c.wfc3_f125w_flux/c.wfc3_f125w_fluxerr)>=3.0)]
      self.mag_array_udf = mag_array_udf
      print "In addition, extend UDF down to %.1f magnitude using SExtractor transfer function only; %d more objects." % \
         (self.MAGLIM_UDF, len(mag_array_udf))
      # Define mlfunc_args
      self.mlfunc_args = [self.datasets, self.limits, self.pixdx,
         self.tfkgrids, self.zdgrids, self.M0,
         self.mc, self.z_mean, self.drop,
         self.add_interloper, self.sd, self.intfrac, self.mag_lolims,
         self.nproc_model, self.zdist_mag_udf, self.kgrid1d_udf,
         self.mag_array_udf, self.MAGLIM_UDF]

def do_fit_idrops_RL(parfile, picklefile=None):
   FU = FitIdropsRL(parfile)
   FU.mlfit()
   if picklefile != None:
      FU.pickle(picklefile)

def do_fit_idrops_RM(parfile, picklefile=None):
   FU = FitIdropsRM(parfile)
   FU.mlfit()
   if picklefile != None:
      FU.pickle(picklefile)

if __name__ == "__main__":
   parfile = sys.argv[1]
   mode = sys.argv[2]
   if len(sys.argv) > 3:
      picklefile = sys.argv[3]
   else:
      picklefile = None
   if mode == 'RL':
      do_fit_idrops_RL(parfile, picklefile=picklefile)
   else:
      do_fit_idrops_RM(parfile, picklefile=picklefile)