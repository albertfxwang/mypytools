#!/usr/bin/env python

from numpy import *
from mypytools import magredshift, igmtrans
import pysynphot as S
from pygoods import *
import matplotlib.pyplot as plt
import lbg_colorcrit as lcc

bands = ['b','v','i','z']
bandpasses = {}
bandpasses['b'] = S.ObsBandpass('acs,wfc2,f435w')
bandpasses['v'] = S.ObsBandpass('acs,wfc2,f606w')
bandpasses['i'] = S.ObsBandpass('acs,wfc2,f775w')
bandpasses['z'] = S.ObsBandpass('acs,wfc2,f850lp')

g04_bdrops = lcc.bdrops_colorcrit()
g04_bdrops('g04', -1., 12.0)
g04_vdrops = lcc.vdrops_colorcrit()
g04_vdrops('g04', -1., 4.5)

def calc_mags_z(zarray=arange(3.0,6.1,0.1), bands=['b','v','i','z']):
	# calculate the B-V v.s. V-z color-color track for the fiducial LBG template SED
	# use E(B-V) = 0.0, 0.15, and 0.3 as the value for dust extinction
	sp_ext0 = S.FileSpectrum('mylbg_sfr10.sed')
	ext15 = S.Extinction(0.15, 'xgal')
	ext30 = S.Extinction(0.30, 'xgal')
	sp_ext15 = sp_ext0 * ext15
	sp_ext30 = sp_ext0 * ext30
	mag_meik = {}
	mag_mad = {}
	mag_meik['zarray'] = zarray
	mag_mad['zarray'] = zarray
	for b in bands:
		mag_meik[b] = {}
		mag_mad[b] = {}
	for b in bands:
		mag_meik[b]['ext0'] = zeros(len(zarray))
		mag_meik[b]['ext15'] = zeros(len(zarray))
		mag_meik[b]['ext30'] = zeros(len(zarray))
		mag_mad[b]['ext0'] = zeros(len(zarray))
		mag_mad[b]['ext15'] = zeros(len(zarray))
		mag_mad[b]['ext30'] = zeros(len(zarray))
	
	for i in range(len(zarray)):
		for b in bands:
			mag_meik[b]['ext0'][i] = magredshift.mag_redshift(sp_ext0, zarray[i], bandpasses[b], 
				igmroutine=igmtrans.meiksin)[0]
			mag_meik[b]['ext15'][i] = magredshift.mag_redshift(sp_ext15, zarray[i], bandpasses[b],
				igmroutine=igmtrans.meiksin)[0]
			mag_meik[b]['ext30'][i] = magredshift.mag_redshift(sp_ext30, zarray[i], bandpasses[b],
				igmroutine=igmtrans.meiksin)[0]
			mag_mad[b]['ext0'][i] = magredshift.mag_redshift(sp_ext0, zarray[i], bandpasses[b],
				igmroutine=igmtrans.madau)[0]
			mag_mad[b]['ext15'][i] = magredshift.mag_redshift(sp_ext15, zarray[i], bandpasses[b],
				igmroutine=igmtrans.madau)[0]
			mag_mad[b]['ext30'][i] = magredshift.mag_redshift(sp_ext30, zarray[i], bandpasses[b],
				igmroutine=igmtrans.madau)[0]
	return mag_meik, mag_mad

def plot_bdrops_colortrack(mag_meik, mag_mad):
	plt.figure()
	# E(B-V) = 0
	bmv_meik_ext0 = mag_meik['b']['ext0'] - mag_meik['v']['ext0']
	vmz_meik_ext0 = mag_meik['v']['ext0'] - mag_meik['z']['ext0']
	# E(B-V) = 0.15
	bmv_meik_ext15 = mag_meik['b']['ext15'] - mag_meik['v']['ext15']
	vmz_meik_ext15 = mag_meik['v']['ext15'] - mag_meik['z']['ext15']
	# E(B-V) = 0.30
	bmv_meik_ext30 = mag_meik['b']['ext30'] - mag_meik['v']['ext30']
	vmz_meik_ext30 = mag_meik['v']['ext30'] - mag_meik['z']['ext30']
	zarray = around(mag_meik['zarray'],1)
	plt.plot(vmz_meik_ext0, bmv_meik_ext0, marker='x', ls='-', c='blue', label='Meiksin; E(B-V)=0')
	plt.plot(vmz_meik_ext15, bmv_meik_ext15, marker='x', ls=':', c='blue', label='Meiksin; E(B-V)=0.15')
	plt.plot(vmz_meik_ext30, bmv_meik_ext30, marker='x', ls='-.', c='blue', label='Meiksin; E(B-V)=0.30')
	for z in [3.0,4.0,5.0,6.0]:
		plt.plot(vmz_meik_ext0[zarray==z], bmv_meik_ext0[zarray==z], marker='o', c='blue', ms=10, mfc='none')
		
	# E(B-V) = 0
	bmv_mad_ext0 = mag_mad['b']['ext0'] - mag_mad['v']['ext0']
	vmz_mad_ext0 = mag_mad['v']['ext0'] - mag_mad['z']['ext0']
	# E(B-V) = 0.15
	bmv_mad_ext15 = mag_mad['b']['ext15'] - mag_mad['v']['ext15']
	vmz_mad_ext15 = mag_mad['v']['ext15'] - mag_mad['z']['ext15']
	# E(B-V) = 0.30
	bmv_mad_ext30 = mag_mad['b']['ext30'] - mag_mad['v']['ext30']
	vmz_mad_ext30 = mag_mad['v']['ext30'] - mag_mad['z']['ext30']
	# Now plot!
	plt.plot(vmz_mad_ext0, bmv_mad_ext0, marker='x', ls='-', c='green', label='Madau; E(B-V)=0')
	plt.plot(vmz_mad_ext15, bmv_mad_ext15, marker='x', ls=':', c='green', label='Madau; E(B-V)=0.15')
	plt.plot(vmz_mad_ext30, bmv_mad_ext30, marker='x', ls='-.', c='green', label='Madau; E(B-V)=0.30')
	for z in [3.0,4.0,5.0,6.0]:
		plt.plot(vmz_mad_ext0[zarray==z], bmv_mad_ext0[zarray==z], marker='o', c='green', ms=10, mfc='none')
	g04_bdrops.plot(color='black')
	plt.xlabel('V - z')
	plt.ylabel('B - V')
	plt.legend(loc=4)
		
		
def plot_vdrops_colortrack(mag_meik, mag_mad):
	plt.figure()
	# E(B-V) = 0
	vmi_meik_ext0 = mag_meik['v']['ext0'] - mag_meik['i']['ext0']
	imz_meik_ext0 = mag_meik['i']['ext0'] - mag_meik['z']['ext0']
	# E(B-V) = 0.15
	vmi_meik_ext15 = mag_meik['v']['ext15'] - mag_meik['i']['ext15']
	imz_meik_ext15 = mag_meik['i']['ext15'] - mag_meik['z']['ext15']
	# E(B-V) = 0.30
	vmi_meik_ext30 = mag_meik['v']['ext30'] - mag_meik['i']['ext30']
	imz_meik_ext30 = mag_meik['i']['ext30'] - mag_meik['z']['ext30']
	zarray = around(mag_meik['zarray'],1)
	plt.plot(imz_meik_ext0, vmi_meik_ext0, marker='x', ls='-', c='blue', label='Meiksin; E(B-V)=0')
	plt.plot(imz_meik_ext15, vmi_meik_ext15, marker='x', ls=':', c='blue', label='Meiksin; E(B-V)=0.15')
	plt.plot(imz_meik_ext30, vmi_meik_ext30, marker='x', ls='-.', c='blue', label='Meiksin; E(B-V)=0.30')
	for z in [3.0,4.0,5.0,6.0]:
		plt.plot(imz_meik_ext0[zarray==z], vmi_meik_ext0[zarray==z], marker='o', c='blue', ms=10, mfc='none')
		
	# E(B-V) = 0
	vmi_mad_ext0 = mag_mad['v']['ext0'] - mag_mad['i']['ext0']
	imz_mad_ext0 = mag_mad['i']['ext0'] - mag_mad['z']['ext0']
	# E(B-V) = 0.15
	vmi_mad_ext15 = mag_mad['v']['ext15'] - mag_mad['i']['ext15']
	imz_mad_ext15 = mag_mad['i']['ext15'] - mag_mad['z']['ext15']
	# E(B-V) = 0.30
	vmi_mad_ext30 = mag_mad['v']['ext30'] - mag_mad['i']['ext30']
	imz_mad_ext30 = mag_mad['i']['ext30'] - mag_mad['z']['ext30']
	# Now plot!
	plt.plot(imz_mad_ext0, vmi_mad_ext0, marker='x', ls='-', c='green', label='Madau; E(B-V)=0')
	plt.plot(imz_mad_ext15, vmi_mad_ext15, marker='x', ls=':', c='green', label='Madau; E(B-V)=0.15')
	plt.plot(imz_mad_ext30, vmi_mad_ext30, marker='x', ls='-.', c='green', label='Madau; E(B-V)=0.30')
	for z in [3.0,4.0,5.0,6.0]:
		plt.plot(imz_mad_ext0[zarray==z], vmi_mad_ext0[zarray==z], marker='o', c='green', ms=10, mfc='none')
	g04_vdrops.plot(color='black')
	plt.xlabel('i - z')
	plt.ylabel('V - i')
	plt.legend(loc=4)