#!/usr/bin/env python

import numpy as np
import scipy
from scipy import integrate, special, optimize
from rebin import rebin
from types import *

"""
Utilities of calculations related to Sersic profiles.
Sersic profile is expressed as

I(r) = I(0) * exp( -bn * (r/Re)**(1/n) )

where bn ~ 2n - 0.324 for 1<=n<=15 (see Trujillo et al., 2001, MNRAS, 326, 869 for
the exact form).

The total flux (integrated to infinity) is

L = kn * Ie * Re**2 = I(0) * ( 2*pi*n / bn**2n ) * Gamma(2n) * Re**2

where Gamma() is the gamma function. One can then use the total magnitude to calculate
I(0).
"""

def sersic_bn(n):
	"""
	Return b_n in the definition of Sersic profile. If 1 <= n <= 15, use the fitting 
	relation bn = 2n - 0.324; otherwise, solve for the exact n by 
	Gamma(2n) = 2*gammainc(2n, bn)
	"""
	if (n>=1) & (n<=15):
		return 2.*n-0.324
	else:
		# need to solve for the exact value of bn
		def sigma(x, n):
			return np.abs(special.gamma(2.*n) - 2.*special.gammainc(2.*n, x))
		bmin = optimize.fmin(sigma, np.maximum(0.1,2.*n-0.324), args=(n,), 
			full_output=0, disp=False, xtol=1.e-6, ftol=1.e-6)
		if n < 0.3:
			print "Warning: error in b_n can be huge!"
		return bmin[0]
	
def sersic_func1d(r, mag, n, Reff, zp=0.):
	"""
	sersic_func(r, mag, n, Re, zp=0.)
	Re is in pixels.
	r -- an array of radius where one evaluates I(r)
	mag -- total magnitude
	n -- Sersic index
	Re -- effective radius
	zp -- photometric zeropoint
	"""
	bn = sersic_bn(n)
	flux_tot = 10.**(-0.4*(mag-zp))  # total flux
	I0 = flux_tot * bn**(2.*n) / (2.*np.pi*n*special.gamma(2.*n)*Reff**2)  
	# I0 is the central surface brightness, in units of flux_unit/pixel
	r = np.abs(r)  # use the absolute value of r
	f = I0 * np.exp( -1*bn*(r/Reff)**(1./n) )
	return f

class sersic1d():
	"""
	Define a class for Sersic functions.
	"""
	def __init__(self, n=1.0, mag=25.0, Reff=4.0, zp=24., pixscale=0.06, 
		scaleunit='arcsec'):
		self.n = n
		self.mag = mag
		self.Reff = Reff  # in units of pixels
		self.zp = zp
		self.func = np.vectorize(lambda r: sersic_func1d(r, self.mag, self.n, self.Reff,
			zp=self.zp))
		self.pixscale = pixscale
		self.scaleunit = scaleunit
		
	def __call__(self, xarr, origin=0., oversample=10):
		""" 
		Evaluate I(r). Only the 1-D case is implemented for now. (12/18/12)
		xarr -- pixel coordinates; starts with zero; each value represents the lower edge
		        of each pixel.
		origin -- the coordinate of the center of the profile, in pixels
		oversample -- the oversampling factor (the factor by which we subdivide the pixels
		              to attain higher accuracy).
		"""
		self.xarr = xarr  # xarr must be regularly spaced
		dx = xarr[1] - xarr[0]  # this should be 1...
		#assert xarr[0] == IntType  
		xarr_ctr = xarr + 0.5  # the pixel centers
		if oversample > 1:
			# calculate I(r) on a new grid
			dx_new = float(dx)/oversample
			xarr_new = np.arange(xarr[0], xarr[-1]+dx, dx_new)
			xarr_ctr_new = xarr_new + dx_new/2.
			# calculate I(r) at the centers of each pixel, using the origin
			r = xarr_ctr_new - origin
			f = self.func(r)
			# now rebin the array
			f = rebin(f, len(self.xarr))
		else:
			r = xarr_ctr - origin
			f = self.func(r)   # evaluate at the centers of each pixel
		return f
			
	def update(self, n=None, mag=None, Reff=None):
		"""
		Updates the values of n, mag, and Reff.
		"""
		if n != None:
			self.n = n
		if mag != None:
			self.mag = mag
		if Reff != None:
			self.Reff = Reff
		self.func = np.vectorize(lambda r: sersic_func(r, self.mag, self.n, self.Reff,
			zp=self.zp))
			
	def integrate(self, R1, R0=0.):
		"""
		Integrate total flux from R0 to R1. Return the integrated magnitude.
		"""
		F = integrate.romberg((lambda x: 2.*np.pi*x*self.__call__(x)), R0, R1)
		return self.zp - 2.5 * np.log10(F)