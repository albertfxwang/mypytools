#!/usr/bin/env python

from numpy import *
import pysynphot as S

# solar mass in kg
M_sun_kg = 1.9891e30  
# AU in cm
AU_cm = 1.496e13
# solar luminosity in W
L_solar_W = 3.846e26
# solar luminosity in erg/s
L_solar_erg_s = L_solar_W * 1.e7
# V-band absolute magnitude
M_V = 4.81
# parsec in cm
parsec_cm = 3.0857e18
# standard bolometric M/L in kg/W
Y_kg_W_bol = 5133.0
# standard bolometric M/L in kg/erg/s
Y_kg_erg_s = 5.133e-4

# reference spectrum of the Sun
sun_spec = S.FileSpectrum('/Users/khuang/lib/pysyn_cdbs/calspec/sun_reference_stis_001.fits')

def mass2light_absmag(absmag, mass_solar, lam, w=100.):
   # calculate the mass-to-light ratio in solar unit from absolute magnitude 
   # around the central wavelength lambda
   # lam: in Angstrom
   # absmag: absolute magnitude
   # mass_solar: mass in solar units
   # width: width of the uniform filter to calculate L
   f_nu = 10.**(-0.4*(absmag+48.6))
   # f_lambda = f_nu * c / lambda**2, in Angstrom
   f_lambda = f_nu * 3.e18 / lam**2  
   # total luminosity around lam within w Angstroms, in erg/s
   L_lambda_erg_s = f_lambda * w * 4 * pi * (10.*parsec_cm)**2
   L_lambda_W = L_lambda_erg_s / 1.e7
   # mass-to-light ratio 
   Y_kg_W = mass_solar * M_sun_kg / L_lambda_W
   print "Y_kg_W", Y_kg_W
   # Now calculate the M/L of Sun at this wavelength, within the same width
   L_lambda_erg_s_sun = sun_spec.sample(lam) * w * (4 * pi * AU_cm**2)
   L_lambda_W_sun = L_lambda_erg_s_sun / 1.e7
   Y_kg_W_sun = M_sun_kg / L_lambda_W_sun
   print "Y_kg_W_sun", Y_kg_W_sun
   # return M/L in solar units
   return Y_kg_W / Y_kg_W_sun
