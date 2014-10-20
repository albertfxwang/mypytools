#!/usr/bin/env python

from numpy import *
from mypytools import cosmoclass
# calculate logR0 at M* given logR0 at M0

cc = cosmoclass.cosmoclass(H0=70.0, omega_m=0.3, omega_l=0.7)

def calc_logr0_mstar(logr0, mstar, beta, z, M0=-21.):
	# return logR0 and R0 in kpc at redshift z
	logr = logr0 - 0.4 * beta * (mstar-M0)
	r_pix = 10.**logr   # in pixels
	r_arcsec = r_pix * 0.03
	r_kpc = cc.adsize(z, r_arcsec, unit='kpc')
	return logr, r_kpc
	
	