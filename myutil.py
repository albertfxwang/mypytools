#!/usr/bin/python

from numpy import *
#from pyraf import iraf
#from iraf import stsdas,hst_calib,ttools,synphot
#import my_pysynphot as S
#import arrayspec as A
import os

pc_cm = 3.0857e18

def closeint(x):
    if fmod(x,1)>=0.5:
        return int(ceil(x))
    else:
        return int(floor(x))

def closepar(a,b):                   ### function to search the closest element to b in list a
    a = sort(a)
    i = searchsorted(a,b)
    if i == len(a):                  ### b is larger than any element in a
        b = a[-1]
    elif i == 0:                     ### b is smaller than any element in a
        b = a[0]
    else:
        if fabs(b-a[i-1])<=fabs(b-a[i]):
            b = a[i-1]
        else: b = a[i]
    return b

def searchpar(Z,tau,rcy,age,ext):    ### function that returns the closest parameter value
    taumin = agemin = 70
    taumax = agemax = 98
    taustep = 2
    rcymax = 35
    rcystep = 5
    Zlist = arange(2,8)
    taulist = arange(taumin,taumax,taustep)
    rcylist = arange(0,rcymax,rcystep)
    extlist = arange(0,60,10)
    Z = closepar(Zlist,Z)
    rcy = closepar(rcylist,rcy)
    tau = closepar(taulist,tau)
    a = arange(agemin,tau-1,1)
    b = arange(tau,agemax,2)
    agelist = concatenate((a,b))  ### agelist be constructed after tau is selected!!
    age = closepar(agelist,age)
    ext = closepar(extlist,ext)  ### discrete extinction
    #if ext > 0.5: ext = 0.5       ### continuous extinction
    #elif ext < 0.0: ext = 0.0

    return (Z,tau,rcy,age,ext)
    
#def convertfit(outfile,cdfile,infile):
#    iraf.tcreate(outfile,cdfile,infile)
    
#def renorm_mag(spec,band,normval,magform):
    # input spec as pysynphot spectrum object
    # band as pysynphot bandpass object
    # normval is the value of abmag that one wants to normalize to
    # will return a renormalized spectrum object AND normalization factor
#    if magform == 'abmag':
#        M = A.ABmag
#    elif magform == 'vegamag':
#        M = A.Vegamag
#    else:
#        raise('unsupported magnitude system')
#    m2 = M(spec,band)
#    factor = 10**(-0.4*(normval - m2))
#    sp = spec*factor
#    return sp, factor

#def renorm_absmag(spec,band,normval,magform,parsec=pc_cm):
    # normalize the spectrum to a certain absolute magnitude
    # placed at 10 pc
#    if magform == 'abmag':
#        M = A.ABmag
#    elif magform == 'vegamag':
#        M = A.Vegamag
#    else:
#        raise('unsupported magnitude system')
#    m2 = M(spec,band)
#    absnormmag = 2.5*log10(4*pi*(10*parsec)**2)
#    mag = m2 + absnormmag   # mag of spec at 10pc
#    normmag = normval - mag
#    normfactor = 10**(-0.4*normmag)
#    sp = spec*normfactor
#    return normmag, normfactor

def trapezoidIntegration(x,y):
    """ Input x and y are numpy arrays"""
    npoints = x.size
    indices = arange(npoints)[:-1]
    deltas = x[indices+1] - x[indices]
    integrand = 0.5*(y[indices+1] + y[indices])*deltas

    return integrand.sum()

def chisquare(model_array, data_array, error_array):
    '''Calculates chi-square given arrays of model output, data points
    and uncertainties of data points. If multiple sources of error exist,
    must add them (in quadrature if neccessary) beforehand'''
    if len(data_array) != len(model_array):
        raise IndexError, 'data and error length do not match'
    temp = (model_array - data_array)/error_array
    temp = temp**2
    chisq = sum(temp)
    return chisq
    
def ffindall(strarray, dir=os.getcwd()):
    '''Finds and returns a list of file names in directory dir that
    contains ALL the strings in str'''
    flist = []
    for f in os.listdir(dir):
        n = 1
        for str in strarray:
            if str not in f:
                n = 0
                break
        if n == 1:
            flist.append(f)
    return flist

def uJy_to_ABmag(flux):
    ab = -2.5*(log10(flux)-29.0)-48.6
    return ab

def ABmag_to_uJy(abmag):
    f = 1.e29 * 10**(-(abmag + 48.6) / 2.5)
    return f

def writeprogress(pdir):
    '''
    Writes a template progress text file for a new directory.
    Input pdir is the directory name.
    '''
    if pdir[-1] == '/':
        pdir = pdir[:-1]
    f = open(pdir+'/progress.txt', 'w')
    pline = '#' * 60
    f.write('Last Updated:\n')
    f.write('Author: Kuang-Han Huang\n')
    f.write('\n')
    f.write('%s\n' % pline)
    f.write('  Description:\n')
    f.write('%s\n' % pline)
    f.write('\n\n')
    f.write('%s\n' % pline)
    f.write('  Works together with the following directory:\n')
    f.write('%s\n' % pline)
    f.write('\n\n')
    f.write('%s\n' % pline)
    f.write('  Completed Work\n')
    f.write('%s\n' % pline)
    f.write('\n\n')
    f.write('%s\n' % pline)
    f.write('  What to do next?\n')
    f.write('%s\n' % pline)
    f.flush()
    f.close()
    print 'Done writing progress in %s' % pdir


def deg2hms(deg):
   x = deg/15
   hh = int(floor(x))
   x = x % 1
   x = x * 60
   mm = int(floor(x))
   x = x % 1
   ss = x * 60
   return hh,mm,ss

def deg2dms(deg):
   if deg > 0.:
      dd = int(floor(deg))
      x = deg - floor(deg)
   else:
      dd = int(ceil(deg))
      x = deg - ceil(deg)
   x = abs(x) * 60
   mm = int(floor(x))
   x = x - floor(x)
   ss = x * 60.
   return dd,mm,ss

def hms2deg(hh,mm,ss):
   deg = hh*15.0 + mm*15.0/60.0 + ss*15.0/3600.0
   return deg

def dms2deg(dd,mm,ss):
   if dd>0:
      deg = dd + mm/60.0 + ss/3600.0
   else:
      deg = dd - mm/60.0 - ss/3600.0
   return deg

def surface_brightness(mag, A):
   """
   Calculates the surface brightness of an object given the total magnitude and
   area in arcsec^2
   The result is in the unit of mag/arcsec^2
   """
   s = m + 2.5 * log10(A)
   return s

def flux2mag(flux, zpt):
   if type(flux) == type(1):
      flux2 = float(flux)
   elif type(flux) == type(1.0):
      flux2 = flux
   elif type(flux) == type(arange(1)):
      flux2 = maximum(flux, 1.e-50)
   mag = zpt - 2.5 * log10(flux2)
   if type(flux) == type(arange(1)):
      putmask(mag, flux < 0, 99.0)
   return mag

def uJy_to_abmag(flux):
  if flux <= 0.:
    return 99.0
  else:
    return	23.9 - 2.5 * log10(flux)

def abmag_to_uJy(mag):
  return 10.**((23.9-mag)/2.5)
	
