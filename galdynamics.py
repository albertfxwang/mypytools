#!/usr/bin/env python

""" 
Utility functions of galactic dynamics.
v1 - 02/23/12
"""
from numpy import *

C=2.99792458e10   # in cm/s
SEC_PER_YEAR=3.15569e7
CM_PER_MPC=3.08568025e24
CM_PER_GPC=CM_PER_MPC*1.e3
CM3_PER_MPC3 = (CM_PER_MPC*CM_PER_MPC*CM_PER_MPC)
CM3_PER_GPC3 = CM_PER_GPC**3
CM_PER_KPC=3.08568025e21
SKYARCMIN       = 1.4850936e8   # arcmin in 4pi steradian
MSOLAR_gm = 1.9891e33   # solar mass in gm
G = 6.674e-8  # cm***3 / gm / sec**2

def t_dyn(R, M, runit='kpc', munit='msolar', tunit='year'):
   # dynamical time: t_dyn = R_eq / v = sqrt(R**3 / 8GM) 
   # first convert to CGS units
   if runit == 'kpc':
      R = R * CM_PER_KPC
   if munit == 'msolar':
      M = M * MSOLAR_gm
   t = sqrt(R**3 / (8 * G * M))
   if tunit == 'year':
      t = t / SEC_PER_YEAR
   return t

