#!/usr/bin/env python

from numpy import *
import bivRL as bl
from transfunc.bivmodel import bivmodel
import gauss
from hconvolve import hconvolve
import time
from transfunc.transfunc_apply import apply_transfer_2d
import interlopers
import cPickle



"""
Construct a bivariate size-stellar mass (RM) distribution, then transform into
size-luminosity distribution using reasonable M/L ratios.
The RM distribution is defined on the log(M_star)--log(R) plane.
The transformation from stellar mass to luminosity (M1500) is via the relation
M1500 = log10(M_star) + X
where X is the conversion factor. X could follow a distribution and does not
have to be a single value. p(X) could also be different for different 
log(M_star) bins.
"""

def schechter_M(logm_arr, alpha_m, logMstar_m):
   """
   Return the Schechter function in logM.
   The Schechter function in logM is

   phi(logM)dlogM = phistar_m * log(10) * 10.**[(alpha_m+1)*(logM-logMstar_m)]
                    * exp[-10.**(logM-logMstar_m)]

   without the normalization constant phistar_m.
   """
   p = log(10) * 10.**((logm_arr-logMstar_m)*(1.+alpha_m))
   p = p * exp(-10.**((logm_arr-logMstar_m)))
   return p

def lognormal(logr_arr, logr_peak, sigma):
   """
   Return a log-normal size distribution.
   """
   p = log(10.)/(sigma*sqrt(2.*pi))
   p = p * exp(-1.*(log(10.)*(logr_arr-logr_peak))**2/(2.*sigma**2))
   return p

def searchsorted_descending(x, element):
   """
   Find the index in the array x that element>x[i] but element<x[i+1].
   x is sorted in descending order, thus different from the built-in 
   searchsorted function.
   """
   x_asc = x[::-1]
   i = searchsorted(x_asc, element)
   j = minimum(len(x)-i, len(x))
   return j

class bivariate_RM_class(bl.bivariate_RL_class):
   def __init__(self, logMDR=7.0):
      self.logMDR = logMDR
      # stellar mass dynamic range

   def bivariate_RM1500_0(self, (alpha_m,logMstar_m,logr0,sigma,beta_m), 
                          limits, pixdx=array([0.02,0.02]), 
                          logM0=10.0, 
                          a=2.5, b=-11.0, kcorr=0.):
      """
      Construct a size-stellar mass distribution, but tranformed onto the 
      (M1500-R) plane. The SMF is assumed to be a Schechter function:
      phi(logM)dlogM = phistar_m * log(10) * 10.**((alph_m+1)(logM-logMstar_m))
                       * exp[-10**(logM-logMstar_m)] * dlogM
      Then, use the variable transformation from logM to M1500 using the 
      equation

      M1500 = (-2.5/a) * logM + b

      to transform phi(logM)dlogM into phi(M1500)dM1500. If kcorr!=0, then 
      further transform into apparent magnitude:

      m = M1500 + kcorr

      Then use log-normal size distribution and a power-law mass-size relation: 
      
      logR_peak = logR0 + beta_m * (logM-logM0)

      To construct the bivariate (M1500, logR) distribution given the 
      parameters of the Schechter SMF (alpha_m, logMstar_m), log-normal size
      distribution parameters (logR0, sigma), and the mass-size relation 
      beta_m.

      The parameters are expressed in stellar mass, as well as the pivot 
      stellar mass logM0 (where logR0 is defined). However, limits and pixdx
      are defined in terms of magnitude (either absolute magnitude at 1500 A 
      or apparent magnitude). This could be confusing, but I can't think of a
      better way of doing this.

      alpha_m:    low-mass end of the SMF
      logMstar_m: knee of the SMF in stellar mass
      logr0:      peak of the size distribution at logM0(=10.0 by default)
      sigma:      width of the size distribution at a given stellar mass
      beta_m:     slope of the size-stellar mass relation
      limits:     model limits in (M1500, logR) or (m, logR)
      pixdx:      pixel size (default to [0.02, 0.02])
      logM0:      the pivot stellar mass where logR0 is defined
      MDR:        dynamic range of the Schechter function in log(M_star)
      a:          scaling factor between M1500 and logM; use a=2.5 for logM,
                  and other factors for different scalings between M1500 and 
                  logM.
      kcorr:      if logMstar is in apparent magnitude, it is the shift in 
                  magnitude between M1500 and the apparent magnitude

      Another use of this method is to convert all M1500 into apparent 
      magnitudes; this way one can work directly in (m, logR) plane, which is
      what one observes.
      """
      dmag = pixdx[0]
      dlogr = pixdx[1]
      npix = around((limits[:,1]-limits[:,0])/pixdx).astype('int')
      npix = abs(npix)
      mag_arr = arange(limits[0][0], limits[0,1], dmag)
      # The magnitude array
      # M1500_arr is an ASCENDING array
      if len(mag_arr) > npix[0]:
         mag_arr = mag_arr[:npix[0]]
      logr_arr = arange(limits[1][0], limits[1][1], dlogr)
      if len(logr_arr) > npix[1]:
         logr_arr = logr_arr[:npix[1]]
      # Now perform variable transformation from M1500 back to logM 
      # (stellar mass)
      # if kcorr != 0: M1500_arr - kcorr is the absolute magnitude
      logM_arr = -0.4 * a * (mag_arr - kcorr - b)
      # Now work in (stellar mass, logR) space!!
      # logM_arr is a DESCENDING array, but all elements should be positive
      # Compute the distribution functions at the grid points and insert into 
      # array
      SMF = schechter_M(logM_arr, alpha_m, logMstar_m)
      # SMF without phistar_m
      logr_peak = logr0 + beta_m * (logM_arr-logM0)
      sd = map(lognormal, tile(logr_arr,(npix[0],1)), logr_peak, 
            ones(npix[0])*sigma)
      # shape(sd)[0] should equal len(logm_arr)
      sd = array(sd).swapaxes(0,1)  
      # flip array so that the 2nd dimension is the same length as logm_arr
      result = sd * SMF
      result = result.swapaxes(0,1)
      # swap the axes back
      # enforce dynamical range in logM so that SMF does not extend to infinity
      logM_uplim = logMstar_m + self.logMDR 
      logM_lolim = logMstar_m - self.logMDR
      index_lolim = searchsorted_descending(logM_arr, logM_lolim)
      if index_lolim >= len(logM_arr): 
         pass
      elif index_lolim == 0: 
         result = zeros(npix)
      else:
         result[index_lolim:,:] = 0.0
      index_uplim = searchsorted_descending(logM_arr, logM_uplim)
      if index_uplim == 0: 
         pass
      elif index_uplim >= len(logM_arr):
         result = zeros(npix)
      else:
         result[:index_uplim,:] = 0.0
      return result

   def model_z(self, params, z, zdgrid, mc):
      """
      Returns the model contribution from redshift z as input to 
      self.convolve_zdist.
      It is a wrapper around self.M2M, which does all the heavy lifting. This
      method only takes care of P(z) grid manipulations.
      """
      if zdgrid != None:
         k = searchsorted(zdgrid.zarr, z)
         model_zdgrid = zeros(self.npix)
         zk = zdgrid.zarr[k] + zdgrid.dz/2.
         kcorr = mc(zk)
         MODEL_k = self.M2M(params, kcorr=kcorr)
         MODEL_k = MODEL_k * zdgrid.Pk[k]
      else:
         kcorr = mc(z)
         MODEL_k = self.M2M(params, kcorr=kcorr)
      return MODEL_k

      pass

   def M2M(self, params, kcorr=0.):
      """
      Calculate RM distribution, then apply the transformation from 
      log(M_star) directly to m_obs.
      self.M2M_kernel should be a dictionary, whose keys are (logMlo,logMhi)
      denoting the boundary of each stellar mass bin, and the values are
      a tuple (X0, p(X) or sigma_X)---sigma_X being the scatter of a 
      Gaussian around X0 *in magnitudes*. The conversion between M1500 and 
      logM is

      M1500 = (-2.5/a) * logM + X
      
      params are in (alpha_m, logMstar_m, logR0, sigma, beta_m)
      """
      # Strategy: convert R-M distribution into R-m (apparent magnitude) 
      # distribution using the "X factor" and the magnitude correction at 
      # each redshift z.
      bivmodel_RL = zeros(self.npix)
      # the kernel width in the X dimension
      kw = self.M2M_kernel['kwidth']       
      # kw is in pixels
      X_ref = self.M2M_kernel['xref']
      # a reference value for X. I will shift RM distribution by X_ref first,
      # and then shift the kernel in each bin to center around the respective
      # X values.
      # if 'type' == 'distribution': not yet implemented...
      if self.M2M_kernel['type'] == 'distribution':
         # the correction is M1500 = -log10(M_star) + X, with X a random
         # variable following probability distribution p(X)
         # First, convert logM into m (apparent magnitude) and calculate the
         # uncorrected R-m distribution converted from RM distribution.
         # Pixel size is unchanged.
         raise ValueError, "distribution-type M2M kernels not yet implemented..."
         mstar_ref = -1.*params[1] + X_ref + kcorr
         m0_ref = -1.*self.logM0 + X_ref + kcorr
         params_ref = params.copy()
         params_ref[1] = mstar_ref
         #print "params_ref", params_ref
         #print "m0_ref", m0_ref
         model0 = self.bivariate_RM_0(params_ref, self.limits, 
                                      pixdx=self.pixdx, logM0=m0_ref,
                                      xdirection=1, a=1)
         #print "max(model0)", max(model0.ravel())
         # Now apply relative shifts for each kernel (relative to X_ref) as
         # well as smearing around X0 in each bin.
         for k in self.M2M_kernel['kernels'].keys():
            if k[1] < self.limits_m[0][1]:
               continue
            elif k[0] > self.limits_m[0][0]:
               continue
            #print "k", k
            X0 = self.M2M_kernel['kernels'][k][0]
            # the reference X of this bin
            logMlo = k[0]
            logMhi = k[1]
            # The limits in M1500 of this stellar mass bin before smearing
            mobs_hi = -1.*logMlo + X0 + kcorr
            mobs_lo = -1.*logMhi + X0 + kcorr
            j0 = round((mobs_lo-self.limits[0][0])/self.pixdx[0])
            j0 = maximum(int(j0), 0)
            j1 = round((mobs_hi-self.limits[0][0])/self.pixdx[0])
            j1 = minimum(int(j1), self.npix[0]-1)
            #print "k, j0, j1", k, j0, j1
            if j1 <= j0:
               continue
            model_RL_k = model0.copy()
            # Copy the section in this log(M_star) bin onto the M1500 grid
            #model_RL_k[j0:j1] = model_RM_k[i0:i1].copy()
            model_RL_k[:j0,:] = 0.
            model_RL_k[j1:,:] = 0.
            # Now smear the model
            X_shift = X0 - X_ref
            # the relative shift between X0 and X_ref---this will be the mean
            # of the Gaussian distribution of this kernel
            X_shift_pix = X_shift / self.pixdx[0]
            #if type(self.M2M_kernel[k][1]) == type(1.0):
            # A Gaussian kernel, with sigma given in magnitude
            sigma_X = self.M2M_kernel['kernels'][k][1]
            # sigma_X is in magnitudes
            sigma_X_pix = sigma_X / self.pixdx[0]
            # sigma_X_pix is in pixels
            if kw % 2 == 1:
               xarr = arange(-(kw-1)/2., (kw+1)/2.)
            else:
               xarr = arange(kw) - kw/2.
            if sigma_X > 0:
               kernel = gauss.gauss(xarr, X_shift_pix, sigma_X_pix)
            else:
               # a delta-function kernel
               kernel = zeros(kw)
               kernel[int(kw/2 + X_shift_pix)] = 1.0
            #else:
            #   # the 1-D smearing kernel
            #   kernel = self.M2M_kernel[k][1].copy()
            # Make kernel into a 2D array (pad with 0 on both sides) for 
            # convolution
            kernel2d = zeros((len(kernel), 3))
            kernel2d[:,1] = kernel
            # smear with p(X)
            model_RL_k = hconvolve(model_RL_k, kernel2d)
            # paste back onto bivmodel_M1500
            bivmodel_RL = bivmodel_RL + model_RL_k
      elif self.M2M_kernel['type'] == 'equation':
         # Uses a linear equation M1500 = -a * log(M_star) + b, and a Gaussian
         # spread around this equation (the uncertainty is the same throughout)
         # all stellar mass bins.
         # If dm is given by self.pixdx[0], then dlogM needs to be rescaled 
         # to self.pixdx[0]/a
         a, b, sigma_X = self.M2M_kernel['kernels']
         sigma_X_pix = sigma_X / self.pixdx[0]
         # if M2M_kernel['type']=='equation', then M2M_kernel['kernels']
         # should contain a list of three numbers: a, b, and sigma_X. The
         # transformation now is M1500 = -a * log(M_star) + b, and sigma_X
         # smears the model in the magnitude direction.
         mstar_ref = (-2.5 / a) * params[1] + b + kcorr
         m0_ref = (-2.5 / a) * self.logM0 + b + kcorr
         params_ref = params.copy()
         params_ref[1] = mstar_ref
         #print "params_ref", params_ref
         #print "m0_ref", m0_ref
         model0 = self.bivariate_RM1500_0(params, self.limits, 
                                      pixdx=self.pixdx, logM0=self.logM0,
                                      a=abs(a), b=b, kcorr=kcorr)
         #return model0
         # Now apply the smearing due to uncertain M/L
         if kw % 2 == 1:
            xarr = arange(-(kw-1)/2., (kw+1)/2.)
         else:
            xarr = arange(kw) - kw/2.
         if sigma_X > 0:
            kernel = gauss.gauss(xarr, 0., sigma_X_pix)
         else:
            # a delta-function kernel
            kernel = zeros(kw)
            kernel[int(kw/2)] = 1.0
         kernel2d = zeros((len(kernel), 3))
         kernel2d[:,1] = kernel
         # smear with p(X)
         bivmodel_RL = hconvolve(model0, kernel2d)
         
      return bivmodel_RL


   def bivariate_RM(self, (alpha_m,logMstar_m,logr0,sigma,beta_m), 
                    limits, pixdx, drop, field, 
                    kgrid=None,verbose=0,
                    zdgrid=None, mc=None, M2M_kernel=None,
                    add_interloper=True,
                    sd=None, intfrac=None, mag_lolim=23.0,
                    norm=-1, logM0=10.0, nproc_model=1):
      ###################### NOTES ON M2M_kernel #############################
      # M2M_kernel should be the key to convert stellar mass into M_1500. 
      # It should be a dictionary with the following entries:
      # 'kwidth': width (length) of the transformation kernel array, in pixels
      # 'type': 'distribution' or 'equation'; if 'type'=='distribution', the
      #         kernel in each mass bin is specified by a central value X0 
      #         (in the equation M1500 = -log(M_star) + X, X in mag) and
      #         a number (in
      #         magnitudes) sppecifying the RMS of the Gaussian distribution
      #         centered around X0; if 'type'=='equation', then kernel is 
      #         specified by three numbers a, b, and sigma_X as in 
      #         M1500 = -2.5/a * log(M_star) + b (note the minus sign in front
      #         of a), with sigma_X being the RMS around this linear equation.
      # 'kernels': it is another dictionary with keys (logM0, logM1) being the 
      #            lower and upper limits of each stellar mass bin,
      #            and each bin could have a different kernel (or one can use
      #            a large range between logM0, logM1 to encompass the entire
      #            parameter space).
      ########################################################################
      self.params = array([alpha_m, logMstar_m, logr0, sigma, beta_m])
      # limits in log(M_star) v.s. log(Re)
      #self.limits_m = limits_m
      # limits in m_obs (apparent magnitude) v.s. log(Re)
      self.limits = limits
      self.pixdx = pixdx
      self.npix = around((self.limits[:,1]-limits[:,0]) / \
                         self.pixdx).astype('int')
      # the pivot log(M_star)
      self.logM0 = logM0
      self.M0 = -21.0
      self.M2M_kernel = M2M_kernel
      self.z_mean = bl.zc_dic[drop]
      self.meankcorr = mc(self.z_mean)  
      self.z0 = self.z_mean-1.5
      self.z1 = self.z_mean+1.5
      # the central redshift of the dropout sample
      # step 1: calculate the contribution to the RL distribution at redshift
      # z from the RM distribution
      model_zdgrid = self.convolve_zdist(self.params,zdgrid, mc,
                                         nproc_model=nproc_model) 
      bivmodel.__init__(self, limits, pixdx, model_zdgrid)
      # STEP 2: add the interloper contributions if add_interloper==True.
      # Need revision for more general applications.
      if add_interloper == True:
         intmodel = self.add_interlopers(intfrac, sd, mag_lolim)
         self.model = self.model + intmodel.model

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
      biv_modelarray = maximum(self.model, 0.)
      #biv_modelarray = maximum(model.model,1.e-50)
      if norm > 0.:
         modsum = sum(biv_modelarray.ravel()) * self.pixdx[0] * self.pixdx[1]
         biv_modelarray = biv_modelarray * (norm / modsum)
      self.model = biv_modelarray

def SMF_M1500(M1500_arr, alpha_m, Mstar_1500, a=2.5):
   """
   Returns the SMF function in Schechter function form, but converted into
   M1500.
   """
   p = log(10) * (0.4 * a) * 10.**(-0.4*(alpha_m+1.)*(M1500_arr - Mstar_1500))
   p = p * exp(-10.**(-0.4 * a * (M1500_arr - Mstar_1500)))
   return p



