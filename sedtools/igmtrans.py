# renamed to igmtrans.py from meiksin.py -- K-H. Huang 02/13/12
# to include Madau IGM attenuation (the data file is in rest-frame wavelength)
#from pygoods import sextutils
from pygoods.sextutils import *
from numpy import *
import os
import pysynphot as S
from pysynphot import spectrum

pytools = os.getenv('mypytools')

def lyman(n):
   invlam = 1.0972e7*(1.-1/n**2)
   lam = 1.e10/invlam
   return(lam)


def init_meiksin():
    zs = arange(1.5,7.5,0.5)

    transmission = {}

    s = sextractor('%s/sedtools/meiksin_transmission.dat' % pytools)
    transmission[1.5] = s.z15
    transmission[2.0] = s.z20
    transmission[2.5] = s.z25
    transmission[3.0] = s.z30
    transmission[3.5] = s.z35
    transmission[4.0] = s.z40
    transmission[4.5] = s.z45
    transmission[5.0] = s.z50
    transmission[5.5] = s.z55
    transmission[6.0] = s.z60
    transmission[6.5] = s.z65
    transmission[7.0] = s.z70
    return s,zs,transmission

def init_madau():
    zs = arange(1.5,7.5,0.5)

    transmission = {}

    s = sextractor('%s/sedtools/madau_transmission.dat' % pytools)
    transmission[1.5] = s.z15
    transmission[2.0] = s.z20
    transmission[2.5] = s.z25
    transmission[3.0] = s.z30
    transmission[3.5] = s.z35
    transmission[4.0] = s.z40
    transmission[4.5] = s.z45
    transmission[5.0] = s.z50
    transmission[5.5] = s.z55
    transmission[6.0] = s.z60
    transmission[6.5] = s.z65
    transmission[7.0] = s.z70
    return s,zs,transmission
    
def restplot(s,zs,transmission):
    waves = array([900,920,1200.])
    dl = 1.
    
    t1 = 0.*zs
    t2 = 0.*zs
    t3 = 0.*zs
    tr = [t1,t2,t3]

    for i in range(len(waves)):
       j = 0
       for z in zs:
         w = waves[i]
         ifilt = (s.lam/(1+z) > w-dl) & (s.lam/(1+z) <= w+dl)
         trans = transmission[z].compress(ifilt).mean()
         # print z,waves[i],trans
         tr[i][j] = trans
         j += 1

    plot(zs,t1)
    plot(zs,t2)
    plot(zs,t3)

def interpolate_transmission(l,lam,trans):
    if l > 1216.: # return 1 at the longest wavelenths
       return 1.0
    if l < lam[0]: # extrapolate at the shortest wavelengths
       i0 = 0
       i1 = 1
       rt = trans[i0]+(l-lam[i0])*(trans[i1]-trans[i0])/(lam[i1]-lam[i0])
    else: # Take the next highest value
       i1 = searchsorted(lam,l)
       rt = trans[i1]
    return rt
       
def resttrans(s,zs,transmission):
    rtdict = {}
    restwave = array([250,300,350,400,450,500,550,600,650,700,
    720,740,760,780,800,820,840,860,880,900,
    909,910.,910.5,911.,911.5,912.,912.5,913.,914.,915.,916., 917., 
    918,919,920,921,922,923,924,925,926,927,928,929,930,931,936,937,938,940,945,948,
    949,950,952,960,970,971,972,973,975,980,990,1000,1010,1020,1022,1024,1025,1026,1027,
    1030,1050,1070,1090,1125,1150,1175,1200,1210,1211,1212,1213,1214,1215,1216,1217,
    1300,1500,2000,4000,6000,8000,10000,15000,20000,40000,60000,80000])
    for z in zs:
        rt = 0.*restwave
        # rt = map(lambda w: interpolate_transmission(w,s.lam/(1.+z),transmission[z]),
        #          restwave)
        for i in range(len(restwave)):
            rt[i] = interpolate_transmission(restwave[i],s.lam/(1+z),transmission[z])
        rtdict[z] = array(rt)
    return restwave,rtdict

def make_rtarray(zs,rtdict):
    zarray = arange(0.,10.,0.05)
    transmission_list = []
    for z in zarray:
        if z < 1.5:
            i1 = 1
            i0 = 0
        if z >= 7.0:
            i1 = 11
            i0 = 10
        if z > 1.5 and z < 7.0:
            i0 = int((z-1.5)*2.)
            i1 = i0+1
        # print z,i0,i1,zs[i0],zs[i1]
        rt0 = maximum(rtdict[zs[i0]],1.e-10)
        rt1 = maximum(rtdict[zs[i1]],1.e-10)

        slope = (log(rt1)-log(rt0))/(zs[i1]-zs[i0])
        logtransmission = log(rt0)+(z-zs[i0])*slope
        transmission = exp(logtransmission)
        transmission = maximum(transmission,0.)
        transmission = minimum(transmission,1.)
        transmission_list += [transmission]
    return zarray,transmission_list

def get_trans(z,zarray,transmission_list):
     i = int(z/0.05)
     if i >= len(zarray):
        i = len(zarray)-1
     return transmission_list[i]

class meiksin:
     def __init__(self):
         s,zs,transmission = init_meiksin()
         self.s = s
         self.zs = zs
         self.transmission = transmission
         self.restwave,self.rtdict = resttrans(s,zs,transmission)
         self.zarray, self.transmission_list = make_rtarray(zs,self.rtdict)
     def __call__(self,z):
         t = get_trans(z,self.zarray,self.transmission_list)
         m = S.ArrayBandpass(wave=self.restwave,throughput=t)
         #return self.restwave,t
         return m   

class madau:
     def __init__(self):
         s,zs,transmission = init_madau()
         self.s = s
         self.zs = zs
         self.transmission = transmission
         self.restwave,self.rtdict = resttrans(s,zs,transmission)
         self.zarray, self.transmission_list = make_rtarray(zs,self.rtdict)
     def __call__(self,z):
         t = get_trans(z,self.zarray,self.transmission_list)
         m = S.ArrayBandpass(wave=self.restwave,throughput=t)
         #return self.restwave,t
         return m