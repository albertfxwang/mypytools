#!/usr/bin/env python

from numpy import *
import pysynphot as S

# calculate Madau's correction for IGM absorption from equation 12 and 13 in Madau 1995
# assume transmission = 0 at wavelengths shortward of Lyman limit
# return a pysynphot TabularSpectralElement class instance

rest_waveset = arange(900., 1300., 1.)   # a set of rest-frame wavelength 

def tau2(lam_obs, z_em):
	# the effective optical depth due to Ly-alpha (2 -> 1)
	# lam_obs: a numpy array of observed wavelength (in A)
	# z_em: the emission redshift
	A2 = 0.0036
	T2 = A2*(lam_obs/1216.)**3.46   # the effective optical depth
	T2 = where(lam_obs > (1216.*(1.+z_em)), 0., T2)
	return T2
	
def tau3(lam_obs, z_em):
	# the effective optical depth due to Ly-beta (3 -> 1)
	# lam_obs: a numpy array of observed wavelength (in A)
	# z_em: the emission redshift
	A3 = 1.7e-3
	T3 = A3*(lam_obs/1026.)**3.46
	T3 = where(lam_obs > (1026.*(1.+z_em)), 0., T3)
	return T3
	
def tau4(lam_obs, z_em):
	# the effective optical depth due to Ly-gamma (4 -> 1)
	# lam_obs: a numpy array of observed wavelength (in A)
	# z_em: the emission redshift
	A4 = 1.2e-3
	T4 = A4*(lam_obs/973.)**3.46
	T4 = where(lam_obs > (973.*(1.+z_em)), 0., T4)
	return T4
	
def tau5(lam_obs, z_em):
	# the effective optical depth due to Ly-delta (5 -> 1)
	# lam_obs: a numpy array of observed wavelength (in A)
	# z_em: the emission redshift
	A5 = 9.3e-4
	T5 = A5*(lam_obs/950.)**3.46
	T5 = where(lam_obs > (950.*(1.+z_em)), 0., T5)
	return T5
	
def tau_photoelectric(lam_obs, z_em):
	# calculate the effective optical depth due to photoelectric absorption
	# shortward of Lyman limit (912A restframe) using the fitting formula in footnote 3
	x_em = 1.+z_em
	x_c = lam_obs / 912.
	tau = 0.25*(x_c**3)*(x_em**0.46-x_c**0.46) + 9.4*(x_c**1.5)*(x_em**0.18-x_c**0.18)
	tau = tau - 0.7*(x_c**3)*(x_c**-1.32-x_em**-1.32) - 0.023*(x_em**1.68-x_c**1.68)
	tau = where(lam_obs<=912., 10., tau)
	tau = where(lam_obs>=(912.*(1.+z_em)), 0., tau)
	return tau


def make_bandpass(lam_obs, z_em):
	lam_ly = [1216*(1.+z_em)+0.01,1216*(1.+z_em)-0.01,1026*(1.+z_em)+0.01,1026*(1.+z_em)-0.01]
	lam_ly += [973*(1.+z_em)+0.01,973*(1.+z_em)-0.01,950*(1.+z_em)+0.01,950*(1.+z_em)-0.01]
	lam_obs = concatenate([lam_obs,lam_ly])
	lam_obs = unique1d(lam_obs)
	lam_obs = sort(lam_obs)
	tau_eff = tau2(lam_obs,z_em)+tau3(lam_obs,z_em)+tau4(lam_obs,z_em)+tau5(lam_obs,z_em)
	#trans = exp(-1.*tau_eff)	
	trans = where(lam_obs<(912.*(1.+z_em)), 0., exp(-1.*tau_eff))  
	# set transmitted flux to zero beyond Lyman limit
	bp = S.ArrayBandpass(wave=lam_obs, throughput=trans)
	return bp
	
def calc_transmission_table(filename='madau_transmission.dat'):
	# calculate the IGM transmission table as a function of OBSERVED wavelength using the 
	# Madau 1995 recipe; only include contributions up to Ly-delta because I can't find
	# the other coefficients for higher-order Lyman-series lines
	zarray = arange(1.5, 7.5, 0.5)  # the redshift array
	lam_obs = arange(1730., 9730., 1.)   # the array of observed wavelength
	transmission = {}
	for z in zarray:
		zstr = '%.1f' % z
		tau_eff = tau2(lam_obs,z)+tau3(lam_obs,z)+tau4(lam_obs,z)+tau5(lam_obs,z)+tau_photoelectric(lam_obs,z)
		trans = where(lam_obs<912., 0., exp(-1.*tau_eff))  # set trans=0 beyond Lyman limit
		transmission[zstr] = trans
		
	# Now write to the file
	headers = []
	headers += ['lam']
	for i in range(len(zarray)):
		headers += ['z%2d' % int(round(10.*zarray[i]))]
	f = open(filename, 'w')
	for j in range(len(headers)):
		f.write('# %d %s\n' % ((j+1), headers[j]))
	for k in range(len(lam_obs)):
		f.write('%f ' % lam_obs[k])
		for j in range(len(zarray)):
			zstr = '%.1f' % zarray[j]
			f.write('%f ' % transmission[zstr][k])
		f.write('\n')
		
	f.close()