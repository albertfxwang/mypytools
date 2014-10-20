#!/usr/bin/env python

# Compute magnitude vs. redshift using pysynphot
# My version of magredshift.py -- by Kuang

import pysynphot as S
from pysynphot import spectrum
from pysynphot.spectrum import TabularSourceSpectrum, ArraySourceSpectrum, ArraySpectralElement
from pysynphot import units
#import meiksin
import igmtrans
from arrayspec import ABmag,Vegamag
import sys
from numpy import *
#import parseconfig
import os
import myutil
import cosmo

pc_cm = 3.0857e18
tenpc = 4.*pi*(10.*pc_cm)**2  # 4*pi*R**2 with R = 10 parsec in centimeter


def readin(file):
    """ Read in an ASCII table of flux vs. redshift """
    s = S.FileSpectrum(file)
    return s

def redshift(spec,z,H0=70.,omega_tot=1.0,omega_m=0.27,omega_l=0.73):
    """Including Cosmology in redshifting spectrum
       Note that the input spectrum should have 'flux unit' in 
       luminosity/lam, since it will be divided by luminosity distance
       later
    """
    # This method does two things:
    # First, shift the wavelengths to higher values by *(1+z)
    # Second, divide the flux by (1+z) b/c wavelength interval is smaller
    # at higher z
    # Third, divide by luminosity distance if H0 > 0. In that case the "flux" unit
    # of spec is actually luminosity (e.g. erg/sec)
    # begin by building a copy of spec, but in Tabular format.
    copy = TabularSourceSpectrum()
    copy.fluxunits = spec.fluxunits    # original flux unit
    copy.waveunits = spec.waveunits    # original wave unit
    fluxunits = spec.fluxunits
    waveunits = spec.waveunits
    copy._wavetable = spec.GetWaveSet()
    copy._fluxtable = spec(copy._wavetable)   # in internal unit

    # flux expressed in photnu is invariant(??) wrt redshift; convert
    # the copy to photnu and extract its flux array.
    copy.convert('flam')  # IN FACT IT SHOULD BE LUMINOSITY/LAM
    # native flux unit is photlam so try to do the calculation in photlam
    #copy.convert('photlam')
    copy.convert('angstrom')   # note that internal flux(in photlam) is unchanged
    (lam, flam) = copy.getArrays()
    #(lam, photlam) = copy.getArrays()

    # convert wavelength array in the copy, to the redshifted frame.
    newlam = lam * (1.+z)    # all lambdas are shifted to a longer value
    copy._wavetable = newlam   # modifies internal wavelength
    #copy._wavetable = newlam
    redshifted_flux = flam / (1. + z)   # lambda interval is shorter at higher z
    #redshifted_flux = photlam / (1. + z)
    # If H0 > 0, then apply cosmology
    if H0 > 0:
       distance = cosmo.lumdist(z,H0,omega_tot,omega_m,omega_l) # in cm
       redshifted_flux = redshifted_flux/(4*pi*distance**2)
    # now comes the trick: put back the flux array in photnu units
    # into the copy, and convert it back to internal units.
    copy._fluxtable = redshifted_flux   # now fluxtable in native photlam       
    copy.fluxunits = units.Units('flam')
    #copy.fluxunits = units.Units('photlam')
    copy.waveunits = units.Units('angstrom')

    copy.ToInternal()
    # skip this step for now
    copy.convert(fluxunits.name)
    copy.convert(waveunits.name)

    return copy


def mag_redshift(spec,z,band,magform='abmag',normmag=None,normband=None,
    igmroutine=igmtrans.meiksin,
    H0=70.,omega_tot=1.0,omega_m=0.27,omega_l=0.73):
    # spec: input spectrum in luminosity/lambda in REST-FRAME
    # any dust extinction should be applied outside this function
    # z: redshift
    # Initialize the IGM routine
    if igmroutine:
        getigm = igmroutine()
    if magform == 'abmag':
        Mag = ABmag
    elif magform == 'vegamag':
        Mag = Vegamag
    else: raise('unsupoorted magnitude system')
    if normmag:
        # normalize to the ABSOLUTE magnitude in normband at normmag IN REST-FRAME
        mag0 = Mag(spec,normband) 
        mag1 = normmag-mag0
        nf = 10.**(mag1/-2.5)
        spec = spec * nf * tenpc # need to convert to luminosity for redshift function
    #w,trans = getigm(normredshift)
    #igm = ArraySpectralElement(w,trans) 
    igm = getigm(z)
    attenuated_spec = spec*igm   # attenuated by IGM **before** redshift

    # now redshift the spectrum using the function defined above 
    redshifted_spec = redshift(attenuated_spec,z,H0=H0,omega_tot=omega_tot,omega_m=omega_m,omega_l=omega_l)
    mag = Mag(redshifted_spec,band)
    return mag, redshifted_spec, attenuated_spec

    
def readbands(file):
    f = open(file,'r')
    lines = f.readlines()
    f.close()
    bandnames = []
    bands = []
    for l in lines:
        a = l.split()
        print a[0],a[1]
        if a[1].find('.fits') < 0:
             bands += [S.ObsBandpass(a[1])]
        else:
             bands += [spectrum.TabularSpectralElement(fileName=a[1])]
        bandnames += [a[0]]
    return bandnames,bands   

class specfile:
    def __init__(self,file):
        f = open(file,'r')
        self.lines = f.readlines()
        self.pointer = 0
        f.close
    def __iter__(self):
        return self
    def next(self):
        if self.pointer >= len(self.lines):
            raise StopIteration
        a = self.lines[self.pointer].split()
        s = S.FileSpectrum(a[0])
        s.otherstuff = " ".join(a[1:])
        self.pointer += 1
        return s

if __name__ == "__main__":
    # Get the parameter file name from the command line
    if len(sys.argv) != 2:
        print "usage: magredshift parfile.par"
        sys.exit()
    c = parseconfig.parseconfig(sys.argv[1])

    # Read the parameter file
    try:
        speclistfile = c['SPECTRA']
        spectra = specfile(speclistfile)
    except:
        spectra = [S.FileSpectrum(c['SPECTRUM'])]
        spectra[0].otherstuff = ""
    bandlistfile = c['BANDS']
    bandnames,bands = readbands(bandlistfile)
#   try:
#       bandlistfile = c['BANDS']
#       bandnames,bands = readbands(bandlistfile)
#   except:
#       band = c['BAND']
#       a = band.split(',')
#       bandnames = [a[-1]]
#       bands = [S.ObsBandpass(band)]
    try:
        redshifts = array(c['REDSHIFTS'])
    except:
        z = c['REDSHIFT_RANGE']
        redshifts = arange(z[0],z[1],z[2])
    if c.has_key('NORM_BAND'):
        normbandlist = c['NORM_BAND']
        normbandstring = ','.join(normbandlist)
        normband = S.ObsBandpass(normbandstring)
        normvalue = c['NORM_ABMAG']
        normredshift = c['NORM_REDSHIFT']
        norm = [normband,normvalue,normredshift]
#       print norm
    else:
        norm = []
    # Get cosmological parameters
    H0 = 70.
    omega_tot = 1.0
    omega_m = 0.27
    omega_l = 0.73
    if c.has_key('H0'): H0 = c['H0']
    if c.has_key('Omega_tot'): omega_tot = c['Omega_tot']
    if c.has_key('Omega_m'): omega_m = c['Omega_m']
    if c.has_key('Omega_lambda'): omega_l = c['Omega_lambda']

    print "%15.15s %15.5s %15.5s" % ("Filename","z","parameters"),
    for j in range(len(bands)):
        print " %8.8s" % (bandnames[j]),
    print ""
    for s in spectra:
        mags = mag_redshift(s,redshifts,bands,norm=norm,
          H0=H0,omega_tot=omega_tot,omega_m=omega_m,omega_l=omega_l)
        for i in range(len(redshifts)):
            print "%15.15s" % (os.path.basename(s.filename)),
            print "%15.5g " % (redshifts[i]),
            print "%-15.15s " % (s.otherstuff),
            for j in range(len(bands)):
                print "%8.3f " % (mags[i,j]),
            print ""

