#!/usr/bin/env python

from numpy import *
from mypytools import cosmoclass
import scipy
from scipy import integrate

"""
Calculate the effective volume of the survey.

The effective volume is calculated as

V_eff = \int_{z1}^{z2} C(z)* (dV(z)/dz) * dz

where C(z) is the selection completeness between z and z+dz, 
(dV/dz) is the differential volume element.

zarray and C(z) are the input parameters to the function. C(z) is determined from the
dropout selection simulation. zarray is the array of redshifts where the value of C(z) is
given. By default C(z) applies between z-dz/2 and z+dz/2, where dz is the step in the 
zarray.
"""

unitlist = ['Mpc3', 'Gpc3']
# cosmology parameters
H0 = 72.0   # km/s/Mpc
omega_m = 0.27 
omega_l = 0.73

# define the cosmoclass object
cc = cosmoclass.cosmoclass(H0, omega_m, omega_l)


def Veff(zarray, Cz, unit='Mpc3', vol='comoving'):
   # calculates the effective surveyed volume
   # can choose between comoving volume and proper volume
   dz = zarray[1] - zarray[0]
   V = 0.
   for i in range(len(zarray)):
      if i == 0:
         zi1 = zarray[0]
         zi2 = zarray[0] + dz/2.
      elif i == len(zarray) - 1:
         zi1 = zarray[-1] - dz/2.
         zi2 = zarray[-1]
      else:
         zi1 = zarray[i] - dz/2.
         zi2 = zarray[i] + dz/2.
      if vol == 'comoving':
         dV = cc.comoving_volume(zi2, unit=unit) - cc.comoving_volume(zi1, unit=unit)
      elif vol == 'proper':
         dV = cc.proper_volume(zi2, unit=unit) - cc.proper_volume(zi1, unit=unit)
      dV = Cz[i] * dV
      V += dV

   return V


def arcmin_to_rad(arcmin):
   deg = arcmin / 60.
   rad = deg * pi / 180.
   return rad

def arcmin2_to_str(arcmin2):
   """from http://www.uniteasy.com/en/unitsCon/solid_angle.htm
   """
   return 8.46159e-8 * arcmin2
