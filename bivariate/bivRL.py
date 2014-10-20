#!/usr/bin/env python
##########  IMPORTANT!!!  ##############
##! bivRL.py (formerly bivariate_lf.py)
########################################
# v1.0.0 H. Ferguson 9/22/03
# v2.1.0 H. Ferguson 1/6/05	- corrected application of radius
#       distribution so that does not renormalize the total
#       to 1 within the observed size range. Instead it is
#       normalized to 1 integrated over the whole size distribution  
#
# Last Updated: K-H Huang 02/05/2013
"""
   Fit a bivariate luminosity-size relation. The LF is assumed to be
   a Schechter function. The size distribution is assumed to be
   log-normal. The size-luminosity relation is assumed to be a powerlaw.

   Usage: bivariate_lf configfile
   configfile - parameters governing the fit
       DATAFILE     # input data file (sextractor catalog format)
       TECHNIQUE    # simplex | grid | anneal (minimization technique) 
       TRANSFER_FCT # Observational transfer function (pickle file)
       MAG_LIMITS   # Magnitude limits
       LOGR_LIMITS  # Limits in np.log(r)
       ALPHA        # Initial guess for alpha
       MSTAR        # Initial guess for Mstar
       RPEAK        # Initial guess for peak of size distribution
       SIGMA        # Initial guess for width of size distribution (e folding)
       M0           # Reference absolute mag for size - luminosity relation
       BETA         # Size-luminosity powerlaw 

   The implicit assumption here is that the bins are finely grided. So 
   the integral of the distribution functions is just the sum over bins * 
   bin size.
   A word about the transfer function. The mlfit routines are completely
   oblivious to this. The transfer-function file gets read in by the
   bivariate_lf subroutine and the transfer function gets applied there
   to each model. The transfer-function file name is a global variable.
"""
__version__ = "3.0"
__author__ = "H. C. Ferguson, STScI / K.-H. HUANG, JHU"

import time
import numpy as np
from transfunc.transfunc_apply import apply_transfer_2d
#import apply_transfer_nd
#from transfunc import transfunc_apply, transfunc_calc
from transfunc.bivmodel import bivmodel
from pygoods import sextractor
#from pygoods.numprint import *
import interlopers
import multiprocessing as mp
import cPickle
import hconvolve
from scipy import signal

# Globals... for GOODS + UDF data set (update later for more general cases)
_first_call = 1
meankcorr_z4_i = 46.013  
# mean mag-correction from rest-frame 1500 A to i-band @ z=4
meankcorr_z5_z = 46.401  
# mean mag-correction from rest-frame 1500 A to z-band @ z=5
# distance modulus + K-correction
limits1 = np.array([[21.0, 26.5], [-2.0, 3.0]])
limits2 = np.array([[23.0, 28.5], [-2.0, 3.0]])
pixdx = np.array([0.02, 0.02])
p1 = np.array([-1.6, -21.0, 0.7, 0.7, 0.3])  # just a sample set of parameters
# total volume surveyed (comoving volume)
A_GOODS = 2.6146e-5  
# total area of GOODS (320 arcmin^2) in steradian, MINUS THE AREA OF UDF
A_UDF = 9.308e-7       
# total area of UDF (11 arcmin^2) in steradian
V_B_GOODS = 1.2519e7  
# in Mpc^3, from z=3.6 to z=4.4 (B-drops) over 320 arcmin^2
V_B_UDF = 4.30355e5    
# in Mpc^3, from z=3.6 to z=4.4 (B-drops) over 11 arcmin^2 of HUDF
V_V_GOODS = 1.1292e7  
# in Mpc^3, from z=4.5 to z=5.5 (V-drops) over 320 arcmin^2
V_V_UDF = 3.88148e5     
# in Mpc^3, from z=4.5 to z=5.5 (B-drops) over 11 arcmin^2 of HUDF

zc_dic = {'f435w':4.0,'f606w':5.0,'f775w':6.0,'uvimos':3.4}  
# keys are the dropout bands

#M0 = -21.		# Reference point of size-magnitude relation


class bivariate_RL_class(bivmodel):
   def __init__(self, MDR=7.0):
      #self.limits = limits
      # limits in magnitude and logR dimensions
      #self.pixdx = pixdx
      #self.M0 = M0
      self.MDR = MDR

   def bivariate_RL_M(self, (alpha,Mstar,logr0,sigma,beta), limits, pixdx,
                      M0=-21.0, MDR=7.):
      """
      Compute the bivariate size-luminosity distribution model on the M-R grid, 
      where M is absolute magnitude (at e.g., rest-frame 1500A).
  
      Arguments:
      alpha,Mstar    -- Schechter-function parameters
      logr0,sigma    -- Lognormal distribution peak, width
                        sigma in ln(r) (to be consistent with what people usually 
                        cite for lambda distribution), 
                        logr0 is in log10(R)
      beta           -- Size-luminosity power-law index
      limits         -- grid limits in ABSOLUTE mag,size
      M0             -- the nominal luminosity where logr0 is evaluated (peak of 
                        the size distribution)
      MDR            -- absolute

      The reference magnitude for the size-luminosity relation is taken
      from the global variable M0
      """
      dmag = pixdx[0]
      dlogr = pixdx[1]  # base 10
      npix = np.around((limits[:,1]-limits[:,0])/pixdx).astype('int')
      M = np.arange(limits[0,0],limits[0,1],dmag)  # array in absolute magnitude
      if abs(M[-1]-limits[0,1])<=dmag/1.e3: 
         # M[-1] shouldn't be upper limit in M
         M = M[:-1]
      logr = np.arange(limits[1][0],limits[1][1],dlogr) 
      # logr = log10(r)
      if abs(logr[-1]-limits[1][1])<=dlogr/1.e3:
         logr = logr[:-1]
      # Compute the distribution functions at the grid points and insert 
      # into array
      result = np.zeros(npix,"float")
      # directly computing Schechter function without importing schechter.py
      delta_M = M - Mstar
      delta_M = np.maximum(delta_M, -1.*MDR)
      delta_M = np.minimum(delta_M, 1.*MDR)
      lf = 0.4*np.log(10.)*10.**(-0.4*(alpha+1.)*delta_M)*\
         np.exp(-10.**(-0.4*delta_M))
      # this is the part of LF without the normalization constant
      for i in range(len(M)):
         expo = -1.*(np.log(10.))**2
         expo = expo*((logr-logr0+0.4*beta*(M[i]-M0)))**2
         expo = expo/(2.*sigma**2)
         sd = np.log(10.)/(np.sqrt(2.*np.pi)*sigma)*np.exp(expo)
         result[i] = lf[i]*sd

      # enforce dynamical range in M so that LF does not extend to infinity
      M_uplim = Mstar + MDR; M_lolim = Mstar - MDR
      index_lolim = (M_lolim - limits[0,0]) / pixdx[0]
      index_lolim = int(round(index_lolim))
      if index_lolim < 0: pass
      elif index_lolim >= npix[0]: 
         result = np.zeros(npix)
      else:
         result[:index_lolim,:] = 0.0
      index_uplim = (M_uplim - limits[0,0]) / pixdx[0] 
      index_uplim = int(round(index_uplim))
      if index_uplim >= npix[0]: pass
      elif index_uplim < 0:
         result = np.zeros(npix)
      else:
         result[index_uplim:,:] = 0.0
      return result

   def model_z(self, params, z, zdgrid, mc):
      """
      A method that returns the RL model contribution at redshift z,
      corrected for P(z) if provided.
      This is the method that should be replaced by sub-classes.
      """
      if zdgrid != None:
         k = np.searchsorted(zdgrid.zarr, z) - 1
         k = np.maximum(k, 0)
         zk = zdgrid.zarr[k] + zdgrid.dz/2.
         m0 = self.M0 + mc(zk)
         park = params.copy()
         park[1] = params[1] + mc(zk)
         MODEL_k = self.bivariate_RL_M(park, self.limits, self.pixdx, m0)
         MODEL_k = MODEL_k * zdgrid.Pk[k]
      else:
         m0 = self.M0 + mc(z)
         p2 = params.copy()
         p2[1] = params[1] + mc(z)
         MODEL_k = self.bivariate_RL_M(p2, self.limits, self.pixdx, m0)
      return MODEL_k

   def convolve_zdist(self, params, zdgrid, mc, nproc_model=1):
      model_zdgrid = np.zeros(self.npix)
      manager = mp.Manager()
      q_out = manager.Queue()
      # create a queue to collect all results
      if zdgrid != None:
         self.z0 = np.maximum(self.z0, zdgrid.zarr[0])
         self.z1 = np.minimum(self.z1, zdgrid.zarr[-1])
         zarr = np.arange(self.z0, self.z1, zdgrid.dz)
         zarr_chunks = np.array_split(zarr, nproc_model)
         processes = []
         for i in range(nproc_model):
            p = mp.Process(target=convolve_dist_worker,
                          args=(self, zdgrid, zarr_chunks[i][0], 
                          zarr_chunks[i][-1], params, 
                          mc, q_out))
            processes += [p]
            p.start()
         # First wait for processes to finish
         for i in range(nproc_model):
            processes[i].join()
         # Now get output
         while not q_out.empty():
            model_zdgrid = model_zdgrid + q_out.get()
      else:
         m0 = self.M0 + self.meankcorr
         p2 = params.copy()
         p2[1] = params[1] + self.meankcorr
         model_zdgrid = self.model_z(params, self.z_mean, zdgrid, self.limits, mc)
      return model_zdgrid

   def add_interlopers(self, intfrac, sd, mag_lolim):
      intmodel = interlopers.interloper_model(self.limits, self, 
                  intfrac, sd, pixdx=self.pixdx, mag_lolim=mag_lolim)
      return intmodel

   def bivariate_RL(self, (alpha,Mstar,logr0,sigma,beta),
                  modellimits,pixdx,drop,field,kgrid=None,verbose=0,
                  zdgrid=None, mc=None, 
                  add_interloper=True, M0=-21.0,
                  sd=None, intfrac=None, mag_lolim=23.0,
                  norm=-1, expand=False, nproc_model=1,
                  delta_z=1.5,
                  floor=1.e-50):
      """
      Usage:
      bivariate_RL((alpha,Mstar,logr0,sigma,beta),
                  modellimits,pixdx,drop,field,kgrid=None,verbose=0,
                  zdgrid=None,mc=None,z0=3.0,z1=6.0,
                  MDR=7.0,M0=-21.0,add_interloper=True,
                  norm=-1)
      modellimits, pixdx: analogous, but in observed frame (m,logr)
      drop: 'b' or 'v' (B-drops or V-drops); add more options later...
      field: 'goods' or 'udf' (also add more options later...)
      kgrid: transfer function (for size-magnitude measurement errors) -- 
      if decide to apply transfer function kernels here
      zdgrid: kernel grids containing P(z) in magnitude and size bins
      mc: K-correction kernel file (from M to m)
      M0: the absolute magnitude where logr0 is evaluated
      MDR: absolute magnitude dynamic range of the model around Mstar
      add_interloper: whether to include interloper contributions or not 
      (need to be more flexible later)
      """
      # enforce positive sigma & beta!
      sigma = abs(sigma); beta = abs(beta)
      self.params = np.array([alpha,Mstar,logr0,sigma,beta])
      limits0 = modellimits
      self.limits = modellimits.copy()
      if expand:
         limits[0][1] += 0.5 
      # add 0.5 mag to the magnitude limit to take into account the number
      # density that is scattered into the magnitude range from fainter 
      # magnitudes?
      self.pixdx = pixdx
      self.npix = np.around((self.limits[:,1]-self.limits[:,0])/self.pixdx).astype('int')
      self.npix0 = np.around((limits0[:,1]-limits0[:,0])/self.pixdx).astype('int')
      # STEP 1: create the model distribution in (m, logr) space. If no dropout-selection
      # kernels are provided (zdgridfile==None), then use the mean redshift to transform
      # the model in the (M, logr) plane to the (m, logr) plane (using the mean redshift).
      # If dropout-selection kernels are provided, apply them to do the conversion from
      # (M, logr) plane to (m, logr) plane.      
      #mc = mconvert(mcfile)
      self.z_mean = zc_dic[drop]
      self.meankcorr = mc(self.z_mean)  # the central redshift of the dropout sample
      self.z0 = self.z_mean - delta_z
      self.z1 = self.z_mean + delta_z
      self.M0 = M0
      self.field = field
      # determine how to convert from absolute mag to apparent mag   
      model_zdgrid = self.convolve_zdist(self.params,zdgrid, mc,
                                         nproc_model=nproc_model) 
      bivmodel.__init__(self, modellimits, pixdx, model_zdgrid) 
      # STEP 2: add the interloper contributions if add_interloper==True.
      # Need revision for more general applications.
      if add_interloper == True:
         intmodel = self.add_interlopers(intfrac, sd, mag_lolim)
         if expand == True:
            self.model[:self.npix0[0],:] = self.model[:self.npix0[0],:] + intmodel.model

      # STEP 3: convolve with GALFIT transfer function kernels if kgrid != None
      if kgrid != None:
         if verbose:
            print "Smearing with transfer function: ",time.ctime()
         smeared_model = apply_transfer_2d(self, kgrid, 
                                           nproc_model=nproc_model)
         self.model = smeared_model.model
         if verbose:
            print "Normalizing: ",time.ctime()
         if verbose:
            print "Returning normalized model: ",time.ctime()
      biv_modelarray = np.maximum(self.model, floor)
      #biv_modelarray = maximum(model.model,1.e-50)
      if norm > 0.:
         modsum = sum(biv_modelarray.ravel()) * self.pixdx[0] * self.pixdx[1]
         biv_modelarray = biv_modelarray * (norm / modsum)
      if expand:
         biv_modelarray = biv_modelarray[:self.npix0[0],:]
         # self.limits[field] = self.limits0.copy()
      self.model = biv_modelarray

def convolve_dist_worker(RLmodel, zdgrid, z0, z1, params, mc, q_out):   
   for z in np.arange(z0, z1, zdgrid.dz):
         MODEL_z = RLmodel.model_z(params, z+zdgrid.dz/2., zdgrid, mc)
         q_out.put(MODEL_z)
   

class univariate_LF_class(bivariate_RL_class):
   def __init__(self, limits0, maglim, zarray, pixdx=np.array([0.02,0.02]), 
                drop='f775w', field='udf', M0=-21.0):
      bivariate_RL_class.__init__(self)
      self.maglim = maglim
      self.limits0 = limits0
      self.pixdx = pixdx
      self.drop = drop
      self.field = field
      self.zarray = zarray
      self.M0 = M0
      self.limits = self.limits0.copy()
      # use 0.5 mag brighter limit to calculate the contribution to the faint-
      # end. Should I experiment with this option?
      self.limits[0][0] = self.limits0[0][1]-0.5
      self.limits[0][1] = self.maglim
      self.npix = np.around((self.limits[:,1]-self.limits[:,0])/self.pixdx)
      self.npix = self.npix.astype('int')
      self.magarr = np.arange(self.limits[0][0],self.limits[0][1],self.pixdx[0])
      self.univariate_LF = np.zeros((len(zarray),self.npix[0]))
      self.univariate_LF_sum = None
      self.LF0 = np.zeros((len(zarray),self.npix[0]))

   def univariate_LF_all(self, params, mc, zdist_mag):
      """
      Calculate the extended LF beyond the bivariate RL magnitude limit for 
      all redshifts.
      """
      self.params = params
      for i in range(len(self.zarray)):
         z = self.zarray[i]
         # The Schechter function at redshift z
         M = self.model_z(params, z, None, mc).sum(axis=1)*self.pixdx[1]
         self.LF0[i] = M.copy()*zdist_mag.dVdz[i]
         # Now multiply by the completeness in each M1500 bin
         # assume that zdist_mag.Pk() has been called
         j = np.searchsorted(zdist_mag.zarr, z)
         M = M * zdist_mag.Pk[j]
         self.univariate_LF[i] = M
      self.univariate_LF_sum = self.univariate_LF.sum(axis=0)
      self.LF0_sum = self.LF0.sum(axis=0)

   def apply_transfer_1d(self, kgrid1d, floor=1.e-50):
      dm = kgrid1d.binwidth
      model_LF_sum = np.zeros(self.npix[0])
      model_all = []
      for k in kgrid1d.kernels.keys():
         m0 = float(k)
         m1 = m0 + dm
         if m1 < self.limits[0][0]:
            continue
         elif m0 > self.limits[0][1]:
            continue
         i0 = np.searchsorted(self.magarr, m0)
         i1 = np.searchsorted(self.magarr, m1)
         if i1 > i0:
            ki = kgrid1d.kernels[k]
            hw = len(ki.kernel) / 2
            model_k = np.zeros(self.npix[0])
            model_k[i0:i1] = self.univariate_LF_sum[i0:i1].copy()
            #print model_k
            model_k = signal.fftconvolve(model_k, ki.kernel, 
                                         mode='full')[hw:-hw]
            model_all += [model_k]
            model_LF_sum = model_LF_sum + model_k
      self.univariate_LF_sum = np.maximum(model_LF_sum, floor)
      #return model_all

   def loglikelihood(self, mag_array):
      """
      Calculate the -loglikelihood!
      """
      mag_pix = np.concatenate([self.magarr, [self.magarr[-1]+self.pixdx[0]]])
      mag_counts = np.histogram(mag_array, mag_pix)[0]
      mag_counts_positive = np.compress(mag_counts>0, mag_counts)
      model_LF = np.compress(mag_counts>0, self.univariate_LF_sum)
      logl = (-np.log10(model_LF)*mag_counts_positive).sum()
      return logl


def whichMbin(M,dropdic):
   # search in which bin M is
   # zdist must be a dictionary of dropouts class objects
   Mbin = -1
   for k in dropdic.keys():
      Mmax = dropdic[k].Mmax
      Mmin = dropdic[k].Mmin
      if (M<Mmax) & (M>=Mmin):
         Mbin = k
         break
   return Mbin 

class mconvert:
   # defines magnitude conversion class (from M_1500 to m)
   def __init__(self, mconvfile):
      c = sextractor(mconvfile)
      self.zarr = c.z
      self.dmarr = c.dm
      self.__doc__ = """
   The first column of mconvfile needs to be redshift, and second column needs to be 
   dm (where M_1500 + dm = m_app).
   """
   def __call__(self, z):
      # does interpolation of magnitude conversion factor dm
      i = np.searchsorted(self.zarr, z)
      if i == 0:
         dz = z - self.zarr[i]  # dz <= 0
         ddm = (self.dmarr[1] - self.dmarr[0]) / (self.zarr[1] - 
                                                  self.zarr[0])*dz
         dm = self.dmarr[0] - ddm
      elif i == len(self.zarr):
         dz = z - self.zarr[i-1]
         ddm = (self.dmarr[-1] - self.dmarr[-2]) / (self.zarr[-1] - 
                                                    self.zarr[-2])*dz
         dm = self.dmarr[-1] + ddm
      else:    
         dz = z - self.zarr[i-1]
         ddm = (self.dmarr[i] - self.dmarr[i-1]) / (self.zarr[i] - 
                                                    self.zarr[i-1])*dz
         dm = self.dmarr[i-1] + ddm
      return dm

   def revert2z(self, dm):
      """
      Given dm, find the corresponding redshift z. Only works when dm is also 
      a monotonically increasing function of z (fortunately it is).
      """
      i = np.searchsorted(self.dmarr, dm)
      if i == 0:
         ddm = dm - self.dmarr[0]
         dz = (self.zarr[1] - self.zarr[0]) / (self.dmarr[1] - 
                                               self.dmarr[0]) * ddm
         z = self.zarr[0] - dz
      elif i == len(self.dmarr):
         ddm = dm - self.dmarr[-1]
         dz = (self.zarr[-1] - self.zarr[-2]) / (self.dmarr[-1] - 
                                                 self.dmarr[-2]) * ddm
         z = self.zarr[-1] + dz
      else:
         ddm = dm - self.dmarr[i]
         dz = (self.zarr[i] - self.zarr[i-1]) / (self.dmarr[i] - 
                                                 self.dmarr[i-1]) * ddm
         z = self.zarr[i] + dz
      return z
      

def Pijk_dVdz(zdgrid, mc, modellimits, pixdx):
   """
   Calculate P(M, R, z)*dV/dz for all redshift bins. Uses the redshift 
   selection function from zdgrid.
   Returns a 3D array Pk, which contains a 2D model in the (m, logR) space for 
   every redshift bin zk. At each redshift zk, the 2D model has the value 
   P_ijk(z)*dV/dz; since P_ijk is divided into many (m, logR) cells, 
   Pk[zk] is a 2D distribution that has constant values in each (m,logR) cell 
   (with the values being P_ijk(z)). Pk[zk] is multiplied (pixel by pixel) by 
   the bivariate distribution phi(M,logR) later on and integrate over M to 
   get the corrected bivariate distribution phi'(m,logR). The conversion 
   between absolute mag M and apparent mag m is done when "pasting" P_ijk(z)
   onto the appropriate (m,logR) cell.
   """
   zarr = zdgrid.zarr + zdgrid.dz / 2. # the central redshift of each bin
   #print zdgrid.zarr,zdgrid.dz
   npix = np.around((modellimits[:,1] - modellimits[:,0]) / pixdx)
   npix = npix.astype('int')
   Pk = np.zeros((len(zarr), npix[0], npix[1]))  # a 3D array, with 2D arrays the same shape as the (m, logR) model in each redshift slice
   #print "shape(Pk)", np.shape(Pk)
   # calculate Pk in the (m, logR) plane
   #def model_zk(k)
   def paste_zd(zdkey, k):
      """
      For each (M1500, logR) cell, transform from M1500 to observed magnitudes.
      """
      zd = zdgrid.zdist[zdkey]
      mbin0 = zd.M0 + mck
      mbin1 = zd.M1 + mck
      # calculate the corresponding limits of each bin limit (M,logR)
      # need to take care of "edge effects"
      ix0 = (mbin0 - modellimits[0,0]) / pixdx[0]; ixc0 = int(np.ceil(ix0))
      if ixc0 >= npix[0]: 
         return 0  # this bin is too faint
      ixc0 = np.maximum(0, ixc0)
      #if ixc0 <= 0: ixc0 = 0
      ix1 = (mbin1 - modellimits[0,0]) / pixdx[0]; ixf1 = int(np.floor(ix1))
      if ix1 < 0: 
         return 0  # this bin is too bright
      ixf1 = np.minimum(npix[0],ixf1)
      #if ixf1 >= npix[0]: ixf1 = npix[0]
      iy0 = round((zd.logR0 - modellimits[1,0]) / pixdx[1]); iy0 = int(iy0)
      if iy0 >= npix[1]: 
         return 0  # this bin is at the size that's too large
      iy1 = round((zd.logR1 - modellimits[1,0]) / pixdx[1]); iy1 = int(iy1)
      if iy1 < 0: 
         return 0  # this bin is at the size that's too small
      # within this cell in (m, logR) space, at redshift zk, this is P(m, logR, zk)
      #Pk[k,ixc0:ixf1,iy0:iy1] = zd.Pz[k]  # paste P_ijk onto Pk
      #try:
      Pk[k,ixc0:ixf1,iy0:iy1] += zd.Pz[k]  
      # paste P_ijk onto Pk; should be += or =?
      #except:
      #   print "[%d,%d:%d,%d:%d]" % (k,ixc0,ixf1,iy0,iy1)
      #   raise ValueError

      # now consider edges of each bin
      if ixc0 != 0:
         Pk[k,ixc0-1,iy0:iy1] += (float(ixc0) - ix0) * zd.Pz[k]
      if ixf1 != npix[0]: 
         Pk[k,ixf1,iy0:iy1] += (ix1 - float(ixf1)) * zd.Pz[k]
      return 1


   #print "start mapping through redshift bins..."
   for k in range(len(zarr)):
      zk = zarr[k]
      mck = mc(zk)
      xout = map(paste_zd, zdgrid.zdist.keys(), 
                 np.ones(len(zdgrid.zdist.keys()),'int')*k)
      Pk[k] = Pk[k] * zdgrid.dVdz[k]
      # zdgrid.dVdz[k] is the volumne enclose by the surveyed area around
      # redshift zk
   #print "finish mapping through redshift bins"
         
   #return Pk
   zdgrid.Pk = Pk  
   # Set attribute for zdgrid
