#!/usr/bin/env python


# Cosmology tools; use a cosmoclass object to store parameters

from numpy import *
import scipy
from scipy import integrate

C=2.99792458e10   # in cm/s
C_KM = C / 1.e5  # in km/s
SEC_PER_YEAR=3.15569e7
CM_PER_MPC=3.08568025e24
CM_PER_GPC=CM_PER_MPC*1.e3
CM3_PER_MPC3 = (CM_PER_MPC*CM_PER_MPC*CM_PER_MPC)
CM3_PER_GPC3 = CM_PER_GPC**3
CM_PER_KPC=3.08568025e21
SKYARCMIN       = 1.4850936e8   # arcmin in 4pi steradian
unitlist = ['cm', 'cm3', 'Mpc', 'Mpc3', 'Gpc3', 'kpc', 'yr', 'Gyr']  # understood distance units

class cosmoclass():
   def __init__(self, H0, omega_m, omega_l):
      # define the basic cosmology (cosmography) parameters
      self.H0 = H0
      self.H0_hertz = H0_hertz(H0)
      self.omega_tot = omega_m + omega_l
      self.omega_m = omega_m
      self.omega_l = omega_l
      self.omega_k = 1. - self.omega_tot
      self.D_H = C / self.H0_hertz  # Hubble distance in cm
      self.funcs = {'printpars()': self.printpars.__doc__,
                    'E(z)': self.E.__doc__,
                    'comoving_distance(z,unit)': self.comoving_distance.__doc__,
                    'pm_distance(z,unit)': self.pm_distance.__doc__,
                    'ang_diam_distance(z,unit)': self.ang_diam_distance.__doc__,
                    'D_A12(z1,z2,unit)': self.D_A12.__doc__,
                    'adsize(z,angsize,unit)': self.adsize.__doc__,
                    'lumdist(z,unit)': self.lumdist.__doc__,
                    'distmod(z)': self.distmod.__doc__,
                    'lookbacktime(z,unit)': self.lookbacktime.__doc__,
                    'age_from_z(z,zform,unit)': self.age_from_z.__doc__,
                    'comoving_volume(z,unit)': self.comoving_volume.__doc__,
                    'proper_volume(z,unit)': self.proper_volume.__doc__}
                    
   
   def __call__(self): # print out all of the available functions
      for k in self.funcs.keys():
         print '%s:' % k, self.funcs[k]


   def printpars(self):
      """
      Print the cosmology parameters used.
      printpars()
      """
      print "H_0 = %.1f Mpc/km/s" % self.H0
      print "Omega_tot = %.3f" % self.omega_tot
      print "Omega_m = %.3f" % self.omega_m
      print "Omega_l = %.3f" % self.omega_l
      print "Omega_k = %.3f" % self.omega_k

   def E(self,z):
      """
      The function relating H0 and H(z): H(z) = H0 * E(z)
      E(z)
      """
      x = sqrt(self.omega_m*(1.+z)**3 + self.omega_k*(1.+z)**2 + self.omega_l)
      return x

   def Hz(self,z):
      """
      H(z) = H0 * E(z)
      Hz(z)
      """
      return self.H0 * self.E(z)

   def comoving_distance(self, z, unit='cm'):
      """
      The line-of-sight integrated comoving distance from redshift 0 out to z
      Different from the Proper Motion Distance, which is the transverse comoving distance
      between 2 objects at the same redshift. Comoving distance in cm.
      **For a flat universe, line-of-sight comoving distance is the same as transverse
      comoving distance or proper motion distance.**
      comoving_distance(z, unit='cm')
      """
      checkunit(unit)
      integral = integrate.quad(lambda x: 1./self.E(x), 0., z)
      #print integral
      D_C = self.D_H * integral[0]
      if unit == 'cm':
         return D_C
      elif unit == 'Mpc':
         return D_C / CM_PER_MPC
   
   def pm_distance(self, z, unit='cm'):
      """
      The comoving distance between two events at the same redshift or distance but separated
      on the sky by some angle d_theta is D_M*d_theta and the transverse comoving distance
      D_M (so-denoted for a reason explained below) is simply related to the line-of-sight
      comving distance D_C:
      ...equations...
      ...............
      The comoving distance happens to be equivalent to the proper motion distance (hence
      the name D_M), defined as the ratio of the actual transverse velocity (in distance over
      time) of an object to its proper motion (in radians per unit time) (Weinberg, 1972, 
      pp 423-424). 
      pm_distance(z, unit='cm')
      """
      checkunit(unit)
      if z < 0.0001:
         dm = z * C / ((1.+z)*self.H0)  # in cm
         if unit == 'cm': return dm
         elif unit == 'Mpc': return dm / CM_PER_MPC
      D_C = self.comoving_distance(z)
      if self.omega_k == 0:
         dm = D_C
      elif self.omega_k > 0:
         ok = sqrt(abs(self.omega_k))
         dm = (self.D_H / ok) * (sinh(ok * D_C / self.D_H))
      else:  # omega_k < 0
         ok = sqrt(abs(self.omega_k))
         dm = (self.D_H / ok) * sin(ok * D_C / self.D_H)
      if unit == 'cm': return dm
      elif unit == 'Mpc': return dm / CM_PER_MPC 

   def ang_diam_distance(self, z, unit='cm'):
      """
      Angular diameter distance at z (per str)
      ang_diam_distance(z, unit='cm')
      """
      checkunit(unit)
      D_A = self.pm_distance(z, unit='cm') / (1.+z)
      if unit == 'cm': return D_A
      elif unit == 'Mpc': return D_A / CM_PER_MPC

   def D_A12(self, z1, z2, unit='Mpc'):
      checkunit(unit)
      if self.omega_k < 0: 
         raise ValueError, 'Omega_k < 0 case not yet implemented.'
      if unit not in unitlist:
         raise ValueError, 'Does not understand the unit.'
      D_M1 = self.pm_distance(z1, unit='cm')
      D_M2 = self.pm_distance(z2, unit='cm')
      term1 = 1./(1.+z2)
      term2 = D_M2 * sqrt(1. + self.omega_k * D_M1**2 / self.D_H**2)
      term3 = D_M1 * sqrt(1. + self.omega_k * D_M2**2 / self.D_H**2)
      da12 = term1 * (term2 - term3)
      if unit == 'cm': return da12
      elif unit == 'Mpc': return da12 / CM_PER_MPC
      

   def adsize(self, z, angsize, unit='cm'):
      """
      Calculates the physical size given the angular size in arcsec.
      adsize(z, angsize, unit='cm')
      """
      checkunit(unit)
      ang_rad = (angsize/3600.) * (pi/180.)
      D_A = self.ang_diam_distance(z)
      if unit == 'cm': return ang_rad * D_A
      elif unit == 'kpc': return ang_rad * D_A / CM_PER_KPC

   def lumdist(self, z, unit='cm'):
      """
      Luminosity distance out to z.
      lumdist(z, unit='cm')
      """
      checkunit(unit)
      D_M = self.pm_distance(z, unit='cm')
      if unit == 'cm': return (1.+z) * D_M
      elif unit == 'Mpc': return (1.+z) * D_M / CM_PER_MPC

   def distmod(self, z):
      """
      Returns the distance modulus.
      distmod(z)
      """
      D_L = self.lumdist(z)
      return cm_to_distmod(D_L)

   def lookbacktime(self, z, unit='Gyr'):
      """
      Calculates the lookback time in years or Gyrs.
      Be careful of numerical error: when z >= 1000. it might become unreliable!
      lookbacktime(z, unit='Gyr')
      """
      checkunit(unit)
      lb = integrate.quad(lambda x: 1./((1.+x)*self.E(x)), 0., z)
      t_l = (1./self.H0_hertz) * lb[0]   # in seconds
      if unit == 'yr':  return t_l / SEC_PER_YEAR
      elif unit == 'Gyr': return t_l / (1.e9 * SEC_PER_YEAR)

   def age_from_z(self, z, zform, unit='Gyr'):
      """
      Calculate the age of an object at z if it was formed at zform.
      age_from_z(z, zform, unit='Gyr')
      """
      checkunit(unit)
      t = self.lookbacktime(z, unit='yr')
      tf = self.lookbacktime(zform, unit='yr')
      dt = tf - t
      if unit == 'yr': return dt
      elif unit == 'Gyr': return dt/1.e9
         
   def comoving_volume(self, z, unit='Gpc3'):
      """
      Calculates the comoving volume out to redshift z over the entire sky.
      comoving_volume(z, unit='Gpc3')
      """
      checkunit(unit)
      D_M = self.pm_distance(z, unit='cm')
      if self.omega_k == 0:
         Vc = (4.*pi/3.) * D_M**3
      else:
         ok = sqrt(abs(self.omega_k))
         term1 = 4.*pi*self.D_H**3 / (2.*self.omega_k)
         term2 = (D_M/self.D_H) * sqrt(1.+self.omega_k*D_M**2/self.D_H**2)
         if self.omega_k > 0:
            term3 = (1./ok) * arcsinh(ok*D_M/self.D_H)
         else:
            term3 = (1./ok) * arcsin(ok*D_M/self.D_H)
         Vc = term1 * (term2 - term3)
      if unit == 'cm3': return Vc
      elif unit == 'Mpc3': return Vc/CM3_PER_MPC3
      elif unit == 'Gpc3': return Vc/CM3_PER_GPC3

   def comoving_volume2(self, z, unit='Gpc3'):
      checkunit(unit)
      ratio = 1.00
      WK = 1.0 - self.omega_m - self.omega_l
      DCMR = self.comoving_distance(z, unit='Mpc') * self.H0 / C_KM
      x = sqrt(abs(WK))*DCMR
      if x > 0.1:
         if WK > 0:
            ratio = (0.125*(exp(2.*x)-exp(-2.*x))-x/2.)/(x*x*x/3.)
         else:
            ratio = (x/2. - sin(2.*x)/4.)/(x*x*x/3.)
      else:
         y = x*x
         if WK < 0: 
            y = -y
         ratio = 1. + y/5. + (2./105.)*y*y
      VCM = ratio*DCMR*DCMR*DCMR/3.
      V_Gpc = 4.*pi*((0.001*C_KM/self.H0)**3)*VCM
      if unit == 'Gpc3':
         return V_Gpc
      elif unit == 'Mpc3':
         return V_Gpc * 1.e9
      else:
         print "unit %s is not yet supported." % unit
         return 0


   def proper_volume(self, z, unit='Gpc3'):
      """
      Calculates the proper volume out to redshift z over the entire sky.
      ** needs check!! **
      proper_volume(z, unit='Gpc3')
      """
      checkunit(unit)
      def dVp(zz):
         # differential proper volume element over 4 pi
         # it's equal to differential comovoing volume element / (1+z)**3
         D_A = self.ang_diam_distance(zz, unit='cm')
         dv = 4. * pi * self.D_H * D_A**2 / ((1.+zz) * self.E(zz))
         return dv

      Vp = integrate.quad(dVp, 0., z)[0]
      if unit == 'cm3': return Vp
      elif unit == 'Mpc3': return Vp / CM3_PER_MPC3
      elif unit == 'Gpc3': return Vp / CM3_PER_GPC3


def checkunit(unit):
   if unit not in unitlist:
      raise ValueError, "Does not understand the unit."


def H0_hertz(H0):
    return H0*3.240777649e-20


def cm_to_distmod(d):
    """distmod(d) - returns distance modulus from distance in cm"""
    return(5.0*log10(d/CM_PER_MPC)+25.0)


