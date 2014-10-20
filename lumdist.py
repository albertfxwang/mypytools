#!/usr/bin/python

##### a tool to calculate luminosity distance #####
##### Use equation 10 in Szydlowski, Kurek, and Krawiec, Phy. Lett. B642, 171 (2006) #####

from numpy import *
from scipy.integrate import romberg

H0_def = 71.0     ### default current value of Hubble's constant in km/s/Mpc
c = 299792.458   ### default value of speed of light in km/s

def H_LamCDM(Omega_m0,z,H0=H0_def):             ### define Lambda-CDM cosmology model H(z)
    H_sq = (H0**2)*(Omega_m0*((1.0+z)**3)+(1.0-Omega_m0))
    return sqrt(H_sq)

def Omega_k0(k,H0=H0_def):
    a0 = 1.0
    q = -k/((H0*a0)**2)
    return q

def FF(x,k):
    if k > 0:
        return sin(x)
    elif k < 0:
        return sinh(x)
    #else:
        #return sinh(x)

def dL(z,k,Omega_m0,H0=H0_def,Hz=H_LamCDM):     ### calculate luminosity distance in Mpc
    def h(x):
        return 1.0/Hz(Omega_m0,x,H0)
    p = H0*sqrt(abs(Omega_k0(k,H0)))*romberg(h,0,z)
    if k:
        q = (1.0+z)*(c/H0)*(1.0/sqrt(abs(Omega_k0(k,H0))))*FF(p,k)
    else:
        q = (1.0+z)*(c/H0)*H0*romberg(h,0,z)
    return q

def Lookback(z1,z2,Omega_m0,Hz=H_LamCDM,H0=H0_def):       ### calculating the lookback time between z1 and z2
    def h(z):
        return 1/((1.0+z)*Hz(Omega_m0,z,H0))
    if z1 >= z2:
        q = romberg(h,z2,z1)*H0
    else:
        q = romberg(h,z1,z2)*H0
    time = q * (100/H0)*9.78                     ### age is in Gyr
    return time
