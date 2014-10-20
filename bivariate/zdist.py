#!/usr/bin/env python

import numpy as N
from pygoods import sextractor, Ftable
import cPickle
import gzip
import KPDFadaptnumpy as KPDF
from bivRL import mconvert, A_GOODS, A_UDF
import os
import cosmoclass
import udropsim_uflux as uu
import matplotlib.pyplot as plt
from dropout_selection import lbg_colorcrit as lcc
import scipy
import numexpr as ne
try:
   from scipy.interpolate import interp1d, RectBivariateSpline, UnivariateSpline
   interp = 1
except ImportError:
   print "Cannot import scipy.interpolate."
   interp = 0

"""
Calculates the redshift distribution of a given dropout sample from simulation.
Uses Kernel PDF estimator to calculate density distribution. Then construct a
zdist class object that stores all the redshift distribution information for 
use in bivariate_lf.py in order to construct bivariate model.
"""

cc = cosmoclass.cosmoclass(70.0, 0.3, 0.7)

def make_Pz(c, bincrit, dropcrit, zarr, n_repeat=1, interpolate=False, 
            dznew=0.1, nfloor=4):
   """
   Calculate completeness & redshift distribution.
   Then construct the zdist instance element.
   Note: zarr is the array of the *lower bound* of each redshift bin.
   If interpolate==True, interpolate in the *REDSHIFT* direction only.
   bincrit and dropcrit should already been "tiled".
   """
   z_in = N.tile(c.z_in,n_repeat)
   detect = N.tile(c.detect,n_repeat)
   dz = zarr[1] - zarr[0]  # native redshift bin width
   # dznew is the bin width of the redshift array at which to calculate P(z); 
   # if interpolate==False, dznew=0.1 is the native bin width (the bin width 
   # used to tally selection fraction from simulations).
   if interpolate==False:
      dznew = zarr[1]-zarr[0]  
      # redshift interval between successive P(z) values
   zarr2 = N.concatenate([zarr,[zarr[-1]+dz]])
   #print "zarr", zarr  
   #print "zarr2", zarr2
   # add the upper boundary to zarr
   selcrit = ne.evaluate("(bincrit==True) & (dropcrit==True)")
   detcrit = ne.evaluate("(bincrit==True) & (detect==True)")
   Pz = N.zeros(len(zarr2))
   zarr_bin = z_in[bincrit==True]
   zarr_det = z_in[detcrit]
   zarr_sel = z_in[selcrit]
   ninput = N.histogram(zarr_bin,zarr2)[0]  
   # number counts in each redshift interval for all objects in 
   # bins (M,logR)-bin
   ndetect = N.histogram(zarr_det,zarr2)[0]
   nselect = N.histogram(zarr_sel,zarr2)[0]  
   # number counts in each redshift interval for all objects in 
   # this (M,logR)-bin 
   # selected as dropout
   Pz = nselect.astype('float') / N.maximum(ninput.astype('float'),1.0)
   # Enforce Pz=0 if ninput < nfloor
   Pz = N.where(ninput < (nfloor*n_repeat), 0., Pz)
   if interpolate==True:
      if interp == 0:
         raise ImportError, "Could not import scipy.interpolate."
      zarr2 = zarr2[:-1]+dz/2.  
      # zarr2 now is the centers of each redshift bin
      #print "zarr2[0], zarr2[-1]", zarr2[0], zarr2[-1]
      fnew = interp1d(zarr2,Pz,kind='cubic')  # cubic interpolation
      #fnew = UnivariateSpline(zarr2, Pz, k=3)
      #print "zarr2[:-1]", zarr2[:-1]
      #print "zarr2[0], zarr2[-1]", zarr2[0], zarr2[-1]
      zarr2 = N.around(zarr2, 2)
      dznew = N.round(dznew, 2)
      zarr2new = N.arange(zarr2[1],zarr2[-1],dznew)
      #print "zarr2new[-1]", zarr2new[-1]
      # Now zarr2new are the centers of the new redshift bins;
      # zarr2new does not go beyond zarr2 
      #print "zarr2new", zarr2new 
      # the redshift coordinate at which to calculate the interpolated P(z)
      Pz = fnew(zarr2new)
      Pz = N.maximum(Pz, 0.)  # enforce P(z) to be positive or zero
      Pz = N.minimum(Pz, 1.)  # enforce P(z) to be <= 1
      # Now, need to return lower-edges of the new redshift bins
      return Pz, ndetect, nselect, ninput, dznew, zarr2new-dznew/2.
   else:
      return Pz, ndetect, nselect, ninput, dz, zarr


class zdist:
   def __init__(self,M0, M1, logR0, logR1, zarr, Pz, ndetect, nselect, ninput):
      # create a zdist instance, but does not calculate P(z)
      self.M0 = M0
      self.M1 = M1
      self.logR0 = logR0
      self.logR1 = logR1
      z0 = zarr[0]
      dz = zarr[1] - zarr[0]
      self.z0 = z0
      self.z1 = zarr[-1] + dz
      self.dz = dz
      self.Pz = Pz
      self.zarr = zarr
      self.ndetect = ndetect
      self.nselect = nselect
      self.ninput = ninput
     

class zdgrid(object):
   def __init__(self, M0, M1, dM, logR0, logR1, dlogR, z0, z1, dz, drop, area):
      """
      Makes P(z) in a grid of M_1500 and logR bins.
      """
      #print "dz =", dz
      # area: surveyed area in steradians;
      self.M0 = M0
      self.M1 = M1
      self.dM = dM
      self.logR0 = logR0
      self.logR1 = logR1
      self.dlogR = dlogR
      self.z0 = z0  # lower z bound of P(z)
      self.z1 = z1  # higher z bound of P(z)
      self.dz = dz
      self.zarr = N.arange(z0, z1, dz)
      self.Mlolims = N.arange(M0, M1, dM)
      self.Mhilims = self.Mlolims + dM
      self.logRlolims = N.arange(logR0, logR1, dlogR)
      self.logRhilims = self.logRlolims + dlogR
      if self.Mhilims[-1] > M1: 
         self.Mhilims[-1] = M1
      self.Mbins = N.concatenate([self.Mlolims,[self.Mhilims[-1]]])
      self.logrbins = N.concatenate([self.logRlolims,[self.logRhilims[-1]]])
      self.zdist = {}
      self.drop = drop
      self.area = area
      #if drop == 'b':
      #   self.drop = 'b'
      #if drop == 'v':
      #   self.drop = 'v'
      # Store the cosmology parameters used
      self.cc = cc
      # calculate dVdz
      self.calc_dVdz()
      self.Pk = None

   def calc_dVdz(self):
      dVdz = N.zeros(len(self.zarr))
      zarr_aug = N.concatenate([self.zarr, [self.zarr[-1]+self.dz]])
      def calc_comoving_V(z):
         V = self.cc.comoving_volume(z, unit='Mpc3')
         return V
      comoving_V = N.array(map(calc_comoving_V, zarr_aug))
      self.dVdz = (comoving_V[1:] - comoving_V[:-1]) * self.area / (4.*N.pi)
      #for i in range(len(self.zarr)):
      #   V1 = cc.comoving_volume(self.zarr[i], unit='Mpc3') * \
      #      self.area / (4.*N.pi)
      #   V2 = cc.comoving_volume(self.zarr[i]+self.dz, unit='Mpc3') * \
      #      self.area / (4.*N.pi)
      #   dVdz[i] = V2 - V1
      #self.dVdz = dVdz
      #self.area = area

   def __call__(self, c, interpolate=False, dznew=0.1, plot_diag=False, 
                ktest=None, n_repeat=1, expand=[0.,0.], mag_in_col='m1500_in',
                re_in_col='re_in'):
      """
      Catalog instance c needs the following two attributes:
      z_in: input redshift array
      detect: array of whether an object is detected or not
      Relegates the determination of self.dropcrit to super classes.
      """
      mag_in = c.__getitem__(mag_in_col)
      re_in = c.__getitem__(re_in_col)
      self.mag_in_col = mag_in_col
      self.re_in_col = re_in_col
      if hasattr(self,'dropcrit')==False:
         raise ValueError, 'Color selection not performed yet.'
      print "expand:", expand
      self.expand = expand
      # Now calculate P(z) for the entire simulation---does not 
      # interpolate here
      self.n_repeat = n_repeat
      self.dznew = dznew
      print "n_repeat", n_repeat
      zarr_old = self.zarr.copy()
      dz_old = zarr_old[1]-zarr_old[0]
      P = make_Pz(c,N.ones(len(c.d)*n_repeat),self.dropcrit,zarr_old,
                  n_repeat=n_repeat,interpolate=False)
      self.Pz_tot,self.ndetect_tot,self.nselect_tot,self.ninput_tot,\
         self.dz_old,self.zarrnew = P
      if plot_diag==True:
         plt.hist(c.z_in[self.dropcrit==True],N.arange(z0,z1+dz,dz),
                  histtype='step',lw=2.0)
      self.dzold = dz_old
      self.zarrold = zarr_old
      if interpolate==True:
         P = make_Pz(c,N.ones(len(c.d)*n_repeat),self.dropcrit,self.zarr,
                     interpolate=interpolate,dznew=dznew,
                     n_repeat=n_repeat)
         self.Pz_tot,self.ndetect_tot,self.nselect_tot,self.ninput_tot,\
            self.dz,self.zarr = P
         # self.zarr are the lower-edges of each redshift bin
         #zarrnew = N.arange(zarr_old[0],zarr_old[-1],dznew)
         #self.zarr = zarrnew   
         # need to change this in order for 
         # bivRL.Pijk_dVdz to work
         self.dz = dznew
         self.dVdz = N.zeros(len(self.zarr))
         # Calculate dVdz
         print "Re-calculate dVdz..."
         self.calc_dVdz()
         #print self.zarr
         #print "self.dznew", self.dznew
         #for i in range(len(self.zarr)):
         #   V1 = cc.comoving_volume(self.zarr[i], unit='Mpc3') * \
         #      self.area / (4.*N.pi)
         #   V2 = cc.comoving_volume(self.zarr[i]+self.dznew, unit='Mpc3') * \
         #      self.area / (4.*N.pi)
         #   self.dVdz[i] = V2 - V1
      # is self.Veff used anywhere else?
      self.Veff = N.sum(self.dVdz*self.Pz_tot)
      # Now calculate P(z) for each (M1500, logR)-bin
      logRe_in = N.log10(getattr(c, self.re_in_col))
      # bin edges need to contain the upper edges
      mindex = N.searchsorted(self.Mbins, mag_in) - 1
      rindex = N.searchsorted(self.logrbins, logRe_in) - 1
      for i in range(len(self.Mlolims)):
         for j in range(len(self.logRlolims)):
            if (ktest!=None) & ((i,j)!=ktest): 
               continue  
               # skip other bins if ktest is specified
            print "(i, j) = (%d, %d)" % (i, j)
            print "M0, logR0 = ", self.Mlolims[i], self.logRlolims[j]
            print "M1, logR1 = ", self.Mhilims[i], self.logRhilims[j]
            bincrit = ((mag_in >= (self.Mlolims[i]-expand[0])) & \
                       (mag_in < (self.Mhilims[i]+expand[0])))
            bincrit = bincrit & (logRe_in >= (self.logRlolims[j]-expand[1])) & \
                     (logRe_in < (self.logRhilims[j]+expand[1]))
            #bincrit = ne.evaluate("(mindex==i) & (rindex==j)")
            print "Total number of points:", N.sum(bincrit)
            if n_repeat>1:
               bincrit = N.tile(bincrit,n_repeat)
            # self.zarr and self.dz already set to zarrnew and dznew if 
            # interpolate==True
            Pz, ndetect, nselect, ninput, dznew, zarrnew = make_Pz(c, bincrit, 
               self.dropcrit, self.zarrold,
               n_repeat=n_repeat,interpolate=interpolate,dznew=dznew)
            k = zdist(self.Mlolims[i],self.Mhilims[i],self.logRlolims[j],
                      self.logRhilims[j],self.zarr,Pz,
                      ndetect.astype('float')/n_repeat, 
                      nselect.astype('float')/n_repeat, 
                      ninput.astype('float')/n_repeat)
            # Record S/N upper & lower limits
            #k.n_snlolim = N.sum(self.sn_lolim_crit[bincrit])
            #k.n_snhilim = N.sum(self.sn_hilim_crit[bincrit])
            #k.n_colorcrit = N.sum(self.colorcrit[bincrit])
            #print "k.n_snlolim", k.n_snlolim
            #print "k.n_snhilim", k.n_snhilim
            #print "k.n_colorcrit", k.n_colorcrit
            if interpolate==True:
               #k.zarrold = k.zarr.copy()
               #k.dzold = k.dz
               k.zarr = zarrnew  # the new lower bounds of each redshift bin
               k.dz = dznew
            self.zdist[(i,j)] = k
      # self.M1500_in = c.m1500_in.copy()
      setattr(self, mag_in_col, mag_in.copy())
      # self.re_in = c.re_in.copy()
      setattr(self, re_in_col, re_in.copy())
      self.z_in = c.z_in.copy()
      if hasattr(self,'c'):
         del self.c  # need to delete the Ftable instance for cPickle to work!


   def delta_zdgrid(self, zc=4.0):
      # create a grid with P(zc) = 1.0 and P(z)=0 elsewhere. For testing purposes.
      for i in range(len(self.Mlolims)):
         for j in range(len(self.logRlolims)):
            k = N.searchsorted(self.zarr, zc)
            Pz = N.zeros(len(self.zarr)); Pz[k] = 1.
            ninput = N.zeros(len(self.zarr),'int'); ninput[k]=1
            ndetect = N.zeros(len(self.zarr),'int'); ndetect[k]=1
            nselect = N.zeros(len(self.zarr),'int'); nselect[k]=1
            kern = zdist(self.Mlolims[i], self.Mhilims[i], self.logRlolims[j],\
               self.logRhilims[j], self.zarr, Pz, ndetect, nselect, ninput)
            self.zdist[(i,j)] = kern

   def flat_zdgrid(self, zlo, zhi):
      # create a grid of P(z) = 1.0 within zlo < z < zhi and P(z)=0 elsewhere
      for i in range(len(self.Mlolims)):
         for j in range(len(self.logRlolims)):
            klo = N.searchsorted(self.zarr, zlo)
            khi = N.searchsorted(self.zarr, zhi)
            ninput = N.zeros(len(self.zarr),'int'); ninput[klo:khi]=1
            ndetect = N.zeros(len(self.zarr),'int'); ndetect[klo:khi]=1
            nselect = N.zeros(len(self.zarr),'int'); nselect[klo:khi]=1
            Pz = N.zeros(len(self.zarr)); Pz[klo:khi]=1.0
            kern = zdist(self.Mlolims[i], self.Mhilims[i], self.logRlolims[j],\
               self.logRhilims[j], self.zarr, Pz, ndetect, nselect, ninput)
            self.zdist[(i,j)] = kern

   # def calc_sel_comp(self,z0=2.5,z1=4.0,show=True):
   #    """
   #    Calculate selection completeness in (M1500, logR) plane.
   #    """
   #    zcrit = N.tile((self.z_in>=z0)&(self.z_in<z1), self.n_repeat)
   #    M1500_in_zr = N.tile(self.M1500_in, self.n_repeat)[zcrit]
   #    #logR_in_zr = N.log10(self.re_in)[(self.z_in>=z0)&(self.z_in<z1)]
   #    logR_in_zr = N.tile(N.log10(self.re_in), self.n_repeat)[zcrit]
   #    #M1500_in_sel_zr = M1500_in_zr[(self.z_in>=z0)&(self.z_in<z1)&\
   #    #   (self.dropcrit[:len(self.M1500_in)]==True)]
   #    selcrit = (zcrit & self.dropcrit)
   #    M1500_in_sel_zr = N.tile(self.M1500_in, self.n_repeat)[selcrit]
   #    logR_in_sel_zr = N.tile(N.log10(self.re_in), self.n_repeat)[selcrit]
   #    #logR_in_sel_zr = N.log10(self.re_in)[(self.z_in>=z0)&(self.z_in<z1)&\
   #    #   (self.dropcrit[:len(self.M1500_in)]==True)]
   #    #logR_in_sel_zr = N.tile(logR_in_sel_zr, self.n_repeat)
   #    n_input_zr = N.histogram2d(M1500_in_zr,logR_in_zr,
   #                               bins=[self.Mbins,self.logrbins])[0]
   #    n_select_zr = N.histogram2d(M1500_in_sel_zr,logR_in_sel_zr,
   #                                bins=[self.Mbins,self.logrbins])[0]
   #    self.sel_comp_map = n_select_zr.astype('float')/N.maximum(n_input_zr,1.0)
   #    if show:
   #       fig = plt.figure()
   #       ax = fig.add_subplot(111)
   #       cax = ax.imshow(self.sel_comp_map.swapaxes(0,1),origin='lower',vmin=0.,vmax=1.0, 
   #          extent=(0,N.shape(n_input_zr)[0],0,N.shape(n_input_zr)[1]))
   #       ax.set_xticks(N.arange(len(self.Mbins))[::2])
   #       ax.set_xticklabels(self.Mbins[::2])
   #       ax.set_yticks(N.arange(len(self.logrbins))[::2])
   #       ax.set_yticklabels(self.logrbins[::2])
   #       ax.set_xlabel('input magnitude',size=14)
   #       ax.set_ylabel('input logRe', size=14)
   #       #ax.set_title('U-dropout selection completeness map in rest-frame 1500 A')
   #       ax.text(0.75,0.9,'z=[%.1f,%.1f]'%(z0,z1),transform=ax.transAxes,color='white',size=12)
   #       cbar = fig.colorbar(cax)

   # def calc_det_comp_M1500(self,z0=2.5,z1=4.0,show=True):
   #    """
   #    Calculate the detection completeness in (M1500,logR) plane within a 
   #    given redshift range
   #    """
   #    self.Mbins = N.concatenate([self.Mlolims,[self.Mhilims[-1]]])
   #    self.logrbins = N.concatenate([self.logRlolims,[self.logRhilims[-1]]])
   #    M1500_in_zr = self.M1500_in[(self.z_in>=z0)&(self.z_in<z1)]
   #    logR_in_zr = N.log10(self.re_in)[(self.z_in>=z0)&(self.z_in<z1)]
   #    M1500_in_det_zr = self.M1500_in[(self.z_in>=z0)&(self.z_in<z1)&(self.detect==True)]
   #    logR_in_det_zr = N.log10(self.re_in)[(self.z_in>=z0)&(self.z_in<z1)&(self.detect==True)]
   #    n_input_zr = N.histogram2d(M1500_in_zr,logR_in_zr,bins=[self.Mbins,self.logrbins])[0]
   #    n_detect_zr = N.histogram2d(M1500_in_det_zr,logR_in_det_zr,bins=[self.Mbins,self.logrbins])[0]
   #    self.det_comp_map = n_detect_zr.astype('float')/N.maximum(n_input_zr,1.0)
   #    if show:
   #       fig = plt.figure()
   #       ax = fig.add_subplot(111)
   #       cax = ax.imshow(self.det_comp_map.swapaxes(0,1),origin='lower',
   #                       vmin=0.,vmax=1.0, 
   #                       extent=(0,N.shape(n_input_zr)[0],
   #                               0,N.shape(n_input_zr)[1]))
   #       ax.set_xticks(N.arange(len(self.Mbins))[::2])
   #       ax.set_xticklabels(self.Mbins[::2])
   #       ax.set_yticks(N.arange(len(self.logrbins))[::2])
   #       ax.set_yticklabels(self.logrbins[::2])
   #       ax.set_xlabel('input magnitude',size=14)
   #       ax.set_ylabel('input logRe', size=14)
   #       #ax.set_title('U-dropout detection completeness map in rest-frame 1500 A')
   #       ax.text(0.75,0.9,'z=[%.1f,%.1f]'%(z0,z1),transform=ax.transAxes,color='white',size=12)
   #       cbar = fig.colorbar(cax)

   def calc_det_comp_m(self, detband='acs_f850lp', show=True):
      """
      Calculate the detection completeness in (m, logR) plane, where m is the
      apparent magnitude in the detection band.
      """
      self.mbins = N.arange(20.,30.5,0.5)
      self.logrbins = N.concatenate([self.logRlolims,[self.logRhilims[-1]]])
      m = getattr(self,'%s_mag'%detband)[:len(self.re_in)]
      logR = N.log10(self.re_in)
      m_detected = m[self.detect==True]
      logR_detected = N.log10(self.re_in)[self.detect==True]
      n_input = N.histogram2d(m, logR, bins=[self.mbins, self.logrbins])[0]
      n_detect = N.histogram2d(m_detected, logR_detected, 
                               bins=[self.mbins, self.logrbins])[0]
      self.det_comp_map_m = n_detect.astype('float') / N.maximum(n_input, 1.0)
      if show:
         fig = plt.figure()
         ax = fig.add_subplot(111)
         cax = ax.imshow(self.det_comp_map_m.swapaxes(0,1), origin='lower',
                         vmin=0., vmax=1.0, 
                         extent=(0,N.shape(n_input)[0],
                                 0,N.shape(n_input)[1]))
         ax.set_xticks(N.arange(len(self.mbins))[::2])
         ax.set_xticklabels(self.mbins[::2])
         ax.set_yticks(N.arange(len(self.logrbins))[::2])
         ax.set_yticklabels(self.logrbins[::2])
         ax.set_xlabel('Input magnitude in %s' % detband, size=14)
         ax.set_ylabel('Input logRe', size=14)
         cbar = fig.colorbar(cax)


   def index(self, M, logR):
   	i = N.searchsorted(self.Mlolims, M)
   	j = N.searchsorted(self.logRlolims, logR)
   	return (i-1,j-1)
   	
   def Pz_bin(self, M, logR):
   	i, j = self.index(M, logR)
   	return self.zdist[(i,j)].Pz

def calc_snlolim_frac(zdgrid):
   s = [len(zdgrid.Mlolims), len(zdgrid.logRlolims)]
   snlolim_frac = N.zeros(s)
   for k in zdgrid.zdist.keys():
      k0 = k[0]
      k1 = k[1]
      snlolim_frac[k0,k1] = float(zdgrid.zdist[k].n_snlolim)/float(N.sum(zdgrid.zdist[k].ninput))
   return snlolim_frac

def calc_snhilim_frac(zdgrid):
   s = [len(zdgrid.Mlolims), len(zdgrid.logRlolims)]
   snhilim_frac = N.zeros(s)
   for k in zdgrid.zdist.keys():
      k0 = k[0]
      k1 = k[1]
      snhilim_frac[k0,k1] = float(zdgrid.zdist[k].n_snhilim)/float(N.sum(zdgrid.zdist[k].ninput))
   return snhilim_frac

def calc_colorcrit_frac(zdgrid):
   s = [len(zdgrid.Mlolims), len(zdgrid.logRlolims)]
   snhilim_frac = N.zeros(s)
   for k in zdgrid.zdist.keys():
      k0 = k[0]
      k1 = k[1]
      snhilim_frac[k0,k1] = float(zdgrid.zdist[k].n_colorcrit)/N.maximum(float(N.sum(zdgrid.zdist[k].ndetect)),1.0)
   return snhilim_frac

def write_zdgrid(zdgrid,outname):
   f = open(outname,'w')
   cPickle.dump(zdgrid,f,2)
   f.close()

def read_zdgrid(filename):
   if os.path.splitext(filename)[1] == '.gz':
      # It is a gzipped pickled file.
      f = gzip.open(filename, 'rb')
   else:
      # It is a regular pickeld file.
      f = open(filename, 'rb')
   zdgrid=cPickle.load(f)
   f.close()
   return zdgrid

def make_zdgrid(c, drop, filename, z0=3.0, z1=6.0, dz=0.1, field='goods',\
   bstonlim=2.0, dlogR=0.2):
   # an example workflow of how to produce a z distribution grid
   # parameters
   # c: the SExtractor catalog instance of the simulation output catalog
   M0 = -25.0
   M1 = -15.0
   dM = 0.5
   logR0 = -0.6
   logR1 = 1.8
   if field == 'goods':
   	area = A_GOODS
   elif field == 'udf':
   	area = A_UDF
   
   zd = zdgrid(M0, M1, dM, logR0, logR1, dlogR, z0, z1, dz, drop, area)
   print zd.dVdz
   zd(c, field=field, bstonlim=bstonlim)
   write_zdgrid(zd,filename)
   
   print "done."

class delta_zdgrid(zdgrid):
   def __call__(self,zc=4.0):
      for i in range(len(self.Mlolims)):
         print i, self.Mlolims[i], self.Mhilims[i] 
         #bincrit = ((c.m1500_input >= self.Mlolims[i]) & (c.m1500_input < self.Mhilims[i]))
         #completeness, p = makezdist(c, bincrit, dropcrit, self.z0, self.z1, self.dz)
         #zarr = N.arange(self.z0, self.z1 + self.dz, self.dz)
         p = N.zeros(len(self.zarr))
         zi = (zc - self.z0) / self.dz
         zi = int(round(zi))
         p[zi] = 1.0
         k = zdist(self.Mlolims[i], self.Mhilims[i], 1.0, p)
         self.zdist[i] = k

class interpolate_zdgrid(zdgrid):
   def __init__(self, zdgridfile):
      """
      A sub-class of zdgrid for interpolating P(z) at finer M1500 and logRe
      bins. The interpolation is performed AFTER the P(z) grids are calculated
      at native bin widths of dM1500 = 0.5 mag, dlogRe = 0.2.
      """
      x = read_zdgrid(zdgridfile)
      self.x = x
      self.filename = zdgridfile
      self.M1500_in = x.M1500_in.copy()
      self.re_in = x.re_in.copy()
      self.z_in = x.z_in.copy()
      self.detect = x.detect.copy()
      self.dropcrit = x.dropcrit.copy()
      self.dVdz = x.dVdz.copy()
      self.Veff = x.Veff
      M0_arr = []
      logR0_arr = []
      for k in x.zdist.keys():
         zdk = x.zdist[k]
         if zdk.M0 not in M0_arr:
            M0_arr += [zdk.M0]
         if zdk.logR0 not in logR0_arr:
            logR0_arr += [zdk.logR0]
      # self.M0_arr and self.logR0_arr are the lower edges of each native bin
      self.dM = x.dM
      self.dlogR = x.dlogR
      self.M0_arr = x.Mlolims.copy()
      self.logR0_arr = x.logRlolims.copy()
      self.Mc_arr = self.M0_arr + self.dM / 2. 
      self.logRc_arr = self.logR0_arr + self.dlogR / 2.
      self.Mlolims = x.Mlolims.copy()
      self.Mhilims = x.Mhilims.copy()
      self.logRlolims = x.logRlolims.copy()
      self.logRhilims = x.logRhilims.copy()
      # self.Mc_arr, self.logRc_arr are the CENTERS of each bin
      self.Mlims = [self.Mlolims[0], self.Mhilims[-1]]
      # The lower & upper limits of M1500
      self.logRlims = [self.logRlolims[0], self.logRhilims[-1]]
      # The lower & upper limits of logR
      self.zarr = x.zarr.copy()
      self.dz = x.dz

      # Now construct the 3-D array of P(z) in [M1500, logR, z] bins
      self.nbins = [len(self.Mlolims), len(self.logRlolims)]
      self.Pijk = N.zeros((self.nbins[0], self.nbins[1], len(self.zarr)))
      # self.Pijk is the 3-D P(z) array
      for k in x.zdist.keys():
         for l in range(len(self.zarr)):
            self.Pijk[k[0], k[1], l] = x.zdist[k].Pz[l]
            # Paste the value of P(z) in the ij-th bin into the (i, j, l)-th
            # element of self.Pijk


   def interpolate_LRbins(self, dM_new, dlogR_new):
      """
      Interpolate P(z) at finer M1500 and/or logR bins.
      Do not interpolate to finer redshift bins... this should be done at 
      the time the grid was constructed.
      """
      # copy the native values of self.dM and self.dlogR
      self.dM_old = self.dM
      self.dlogR_old = self.dlogR
      self.Pijk_old = self.Pijk.copy()
      self.M0_old_arr = self.M0_arr.copy()
      self.logR0_old_arr = self.logR0_arr.copy()
      self.Mc_old_arr = self.Mc_arr.copy()
      self.logRc_old_arr = self.logRc_arr.copy()
      self.nbins_old = self.nbins
      self.Mlolims_old = self.Mlolims.copy()
      self.Mhilims_old = self.Mhilims.copy()
      self.logRlolims_old = self.logRlolims.copy()
      self.logRhilims_old = self.logRhilims.copy()

      # assign new values of dM and dlogR
      self.dM = dM_new
      self.dlogR = dlogR_new
      self.M0_arr = N.arange(self.Mlims[0], self.Mlims[-1], dM_new)
      # the new array of lower-edges in M1500
      self.logR0_arr = N.arange(self.logRlims[0], self.logRlims[-1], 
                              dlogR_new)
      # now calculate the new values of the centers of each bin
      self.Mc_arr = self.M0_arr + self.dM / 2.
      self.logRc_arr = self.logR0_arr + self.dlogR / 2.
      self.Mlolims = self.M0_arr.copy()
      self.Mhilims = self.Mlolims + self.dM
      self.logRlolims = self.logR0_arr.copy()
      self.logRhilims = self.logRlolims + self.dlogR
      self.nbins = [len(self.Mlolims), len(self.logRlolims)]

      # the new array of lower-edges in logR
      self.Pijk_new = N.zeros([len(self.M0_arr), len(self.logR0_arr),
                            len(self.zarr)])
      # Now start interpolating at each redshift
      for l in range(len(self.zarr)):
         Pz = self.Pijk_old[:,:,l]
         Pz_interp = RectBivariateSpline(self.Mc_old_arr, self.logRc_old_arr,
                                         Pz)
         # a RectBivariateSpline class instance that will calculate the
         # interpolated P(z)
         self.Pijk_new[:,:,l] = Pz_interp(self.Mc_arr, self.logRc_arr)

      # Make sure self.Pijk_new is between 0 and 1
      self.Pijk_new = N.minimum(self.Pijk_new, 1.0)
      self.Pijk_new = N.maximum(self.Pijk_new, 0.0)
      # Now copy the P(z) in each new bin into self.zdist
      self.zdist = {}
      for i in range(self.nbins[0]):
         for j in range(self.nbins[1]):
            M0_bin = self.Mlolims[i]
            M1_bin = self.Mhilims[i]
            logR0_bin = self.logRlolims[j]
            logR1_bin = self.logRhilims[j]
            print "Interpolating [%.1f, %.1f]" % (M0_bin, logR0_bin)
            #bincrit = (self.x.M1500_in>=M0_bin) & (self.x.M1500_in<M1_bin)
            #bincrit = bincrit & (N.log10(self.x.re_in)>=logR0_bin) & \
            #                    (N.log10(self.x.re_in)<logR1_bin)
            #ninput_bin = sum(bincrit)
            #ndetect_bin = sum((self.x.detect==True) & (bincrit==True))
            #nselect_bin = sum((self.x.dropcrit==True) & (bincrit==True))
            ninput_bin = 0
            ndetect_bin = 0
            nselect_bin = 0
            self.zdist[(i,j)] = zdist(self.Mlolims[i], self.Mhilims[i],
                                      self.logRlolims[j], self.logRhilims[j],
                                      self.zarr, self.Pijk_new[i,j,:],
                                      ndetect_bin, nselect_bin, ninput_bin)


   def write(self, filename):
      delattr(self, 'x')
      f = open(filename, 'wb')
      cPickle.dump(self, f, 2)
      f.close()




